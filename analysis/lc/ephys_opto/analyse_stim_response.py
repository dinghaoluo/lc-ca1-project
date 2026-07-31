# -*- coding: utf-8 -*-
'''
Created on Thu 13 Feb 14:36:41 2025

analyse and plot stim responses for LC stim recordings

@author: Dinghao Luo
'''

#%% imports
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import fftconvolve

repo_root = Path(__file__).resolve().parents[3]
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))

from behaviour_functions import process_txt
from plotting_functions import plot_violin_with_scatter
from common_functions import gaussian_kernel_unity, mpl_formatting
from console_formatting import print_session
import project_paths as pp
mpl_formatting()

import rec_list
paths = rec_list.pathLCopt


#%% paths and parameters
mice_exp_stem    = pp.MICEEXP_ROOT
all_session_stem = pp.LC_EPHYS_STEM / 'all_sessions'
stim_raster_stem = pp.LC_EPHYS_FIGURES_STEM / 'stim_effects' / 'single_cell_stim_rasters'
stim_raster_stem.mkdir(parents=True, exist_ok=True)

BEF = 3  # s
AFT = 7  # s
SAMP_FREQ = 1250  # hz

MAX_LENGTH = int((BEF + AFT) * SAMP_FREQ)
XAXIS = np.arange(MAX_LENGTH) / SAMP_FREQ - BEF

AMP_WINDOW_LOW_S  = 0  # s
AMP_WINDOW_HIGH_S = 1
AMP_WINDOW_LOW    = int((AMP_WINDOW_LOW_S + BEF) * SAMP_FREQ)
AMP_WINDOW_HIGH   = int((AMP_WINDOW_HIGH_S + BEF) * SAMP_FREQ)

SIGMA_SPIKE = int(SAMP_FREQ * 0.10)  # 50 ms
GAUS_SPIKE = gaussian_kernel_unity(SIGMA_SPIKE, GPU_AVAILABLE=False)


#%% load data
cell_prop = pd.read_pickle(pp.LC_EPHYS_STEM / 'LC_all_cell_profiles.pkl')

tagged_keys, putative_keys = [], []
for clu in cell_prop.itertuples():
    if clu.identity == 'tagged':
        tagged_keys.append(clu.Index)
    if clu.identity == 'putative':
        putative_keys.append(clu.Index)


#%% main loop
all_ctrl_amps = []
all_stim_amps = []

for path in paths:
    recname = Path(path).name
    print_session(recname)

    txtpath = mice_exp_stem / f'ANMD{recname[1:4]}r' / recname[:-3] / recname / f'{recname}T.txt'
    beh = process_txt(txtpath)

    # stim condition per trial
    stim_cds = [trial[15] for trial in beh['trial_statements']][1:]
    stim_idx = [ti for ti, cond in enumerate(stim_cds) if cond != '0']
    ctrl_idx = [ti + 2 for ti in stim_idx]

    clu_iter = list(cell_prop[cell_prop['recname'] == recname].index)
    rasters_run_disk = np.load(
        all_session_stem / recname / f'{recname}_all_rasters_run.npy',
        allow_pickle=True
    ).item()

    for clu in clu_iter:
        if clu not in tagged_keys:
            continue

        # load run-aligned
        run_aligned = rasters_run_disk[clu][ctrl_idx]

        # load stim-aligned
        stim_aligned = rasters_run_disk[clu][stim_idx]

        ctrl_mean = np.nanmean(run_aligned, axis=0)
        stim_mean = np.nanmean(stim_aligned, axis=0)
        ctrl_prof = fftconvolve(ctrl_mean, GAUS_SPIKE, mode='same') * SAMP_FREQ
        stim_prof = fftconvolve(stim_mean, GAUS_SPIKE, mode='same') * SAMP_FREQ

        # collect for summary stats
        all_ctrl_amps.append(np.nanmean(ctrl_prof[AMP_WINDOW_LOW:AMP_WINDOW_HIGH]))
        all_stim_amps.append(np.nanmean(stim_prof[AMP_WINDOW_LOW:AMP_WINDOW_HIGH]))

        # PLOTTING
        # use stim range for BOTH twin axes
        ymax = max(np.nanmax(stim_prof), np.nanmax(ctrl_prof)) * 1.05

        fig, axs = plt.subplots(2, 1, figsize=(1.8, 2.25), sharex=True)

        stim_rows, stim_cols = np.where(stim_aligned > 0)

        axs[0].scatter(
            stim_cols / SAMP_FREQ - BEF,
            stim_rows + 1,
            s=0.8,
            color='lightsteelblue',
            ec='none'
            )

        ctrl_rows, ctrl_cols = np.where(run_aligned > 0)

        axs[1].scatter(
            ctrl_cols / SAMP_FREQ - BEF,
            ctrl_rows + 1,
            s=0.8,
            color='grey',
            ec='none'
            )

        axs[0].set(title=f'{clu}\nStim.', ylabel='Trial #')
        axs[1].set(title='Ctrl.', ylabel='Trial #', xlabel='Time from run onset (s)')

        # overlays (same y scale)
        axt0 = axs[0].twinx()
        axt0.plot(XAXIS, stim_prof, color='royalblue', lw=1)
        axt0.set_ylim(0, ymax)
        axt0.set(ylabel='Firing rate (Hz)')
        axt0.spines['top'].set_visible(False)

        axt1 = axs[1].twinx()
        axt1.plot(XAXIS, ctrl_prof, color='k', lw=1)
        axt1.set_ylim(0, ymax)
        axt1.set(ylabel='Firing rate (Hz)')
        axt1.spines[['top', 'left', 'bottom', 'right']].set_visible(False)

        for ax in axs:
            ax.set(xlim=(-1, 4), xticks=(0, 2, 4))

        for ax in axs:
            ax.spines['top'].set_visible(False)
        axt0.spines['top'].set_visible(False)
        axt1.spines['top'].set_visible(False)

        fig.tight_layout()

        tag = 'tagged' if clu in tagged_keys else 'putative'
        for ext in ['.png', '.pdf']:
            fig.savefig(
                stim_raster_stem / f'{clu}_{tag}{ext}',
                dpi=300,
                bbox_inches='tight'
            )
        plt.close(fig)


#%% summary stats
plot_violin_with_scatter(
    all_ctrl_amps,
    all_stim_amps,
    'grey',
    'royalblue',
    ylabel='Firing rate (Hz)',
    xticklabels=['Ctrl.', 'Stim.'],
    print_statistics=True,
    save=True,
    savepath=pp.LC_EPHYS_FIGURES_STEM / 'stim_effects' / 'tagged_putative_ctrl_stim_violin'
)
