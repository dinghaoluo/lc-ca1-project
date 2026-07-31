# -*- coding: utf-8 -*-
'''
Created on Fri 20 Dec 17:30:12 2024
Modified on 10 May Sat 2025

plot rasters of HPC cells in ctrl and stim trials
modified to label PyrUp and PyrDown cells

@author: Dinghao Luo
'''

#%% imports
import numpy as np
import pickle
import pandas as pd
import matplotlib.pyplot as plt
import sys
from pathlib import Path

repo_root = Path(__file__).resolve().parents[2]
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))
from common_functions import mpl_formatting
from console_formatting import print_session
import project_paths as pp
mpl_formatting()

import rec_list
pathHPCLC = rec_list.pathHPCLCopt
pathHPCLCterm = rec_list.pathHPCLCtermopt
paths = pathHPCLC + pathHPCLCterm


#%% load dataframe
print('loading dataframe...')
HPC_STEM = pp.HPC_EPHYS_STEM
HPCLC_BEH_STEM = pp.behaviour_experiment_stem('HPCLC')
HPCLCTERM_BEH_STEM = pp.behaviour_experiment_stem('HPCLCterm')

cell_profiles = pd.read_pickle(HPC_STEM / 'hpc_all_profiles.pkl')
df_pyr = cell_profiles[cell_profiles['cell_identity']=='pyr']  # pyramidal only


#%% parameters
time_bef = 1  # second
time_aft = 4
samp_freq = 1250  # Hz
PYRUP_CLASS = 'run-onset ON'
PYRDOWN_CLASS = 'run-onset OFF'


#%% main
for path in paths[14:]:
    recname = path[-17:]
    print_session(recname)

    curr_df_pyr = df_pyr[df_pyr['recname']==recname]
    rasters = np.load(
        HPC_STEM / 'all_sessions' / recname / f'{recname}_all_rasters.npy',
        allow_pickle=True
        ).item()

    beh_stem = HPCLC_BEH_STEM if path in pathHPCLC else HPCLCTERM_BEH_STEM
    with open(beh_stem / f'{recname}.pkl', 'rb') as f:
        beh = pickle.load(f)

    stim_conds = [t[15] for t in beh['trial_statements']][1:]
    stim_idx = [trial for trial, cond in enumerate(stim_conds)
                if cond!='0']
    ctrl_idx = [trial+2 for trial in stim_idx]
    baseline_idx = list(np.arange(stim_idx[0]))

    # extract stim times
    pulse_times = []
    for pulse in beh['pulse_times']:
        pulse_arr = np.asarray(pulse, dtype=float).ravel()
        pulse_times.extend(pulse_arr[np.isfinite(pulse_arr)])
    pulse_times = np.sort(np.asarray(pulse_times, dtype=float))

    run_onsets = np.asarray(beh['run_onsets'], dtype=float)
    pulse_times_aligned = []
    for trial in stim_idx:
        run_idx = trial + 1
        if run_idx >= len(run_onsets) or not np.isfinite(run_onsets[run_idx]):
            pulse_times_aligned.append(np.nan)
            continue

        run = run_onsets[run_idx]
        next_run = np.nan
        for next_run_candidate in run_onsets[run_idx + 1:]:
            if np.isfinite(next_run_candidate):
                next_run = next_run_candidate
                break
        trial_end = next_run if np.isfinite(next_run) else run + time_aft * 1000
        trial_pulses = pulse_times[
            (pulse_times >= run - time_bef * 1000)
            & (pulse_times < trial_end)
        ]

        if len(trial_pulses) == 0:
            pulse_times_aligned.append(np.nan)
        else:
            pulse_times_aligned.append(trial_pulses[0] - run)

    curr_df_pyrup = curr_df_pyr[curr_df_pyr['class_ctrl'] == PYRUP_CLASS]
    curr_df_pyrdown = curr_df_pyr[curr_df_pyr['class_ctrl'] == PYRDOWN_CLASS]
    pyr_list = curr_df_pyr.index.tolist()
    pyrup_list = curr_df_pyrup.index.tolist()
    pyrdown_list = curr_df_pyrdown.index.tolist()

    tot_time = 1250 + 5000  # 1 s before, 4 s after

    for cluname in pyr_list:
        if cluname in pyrup_list:
            clustr = f'{cluname} PyrUp'
        elif cluname in pyrdown_list:
            clustr = f'{cluname} PyrDown'
        else:
            clustr = f'{cluname} other'

        raster = rasters[cluname]

        ctrl_matrix = raster[ctrl_idx]
        stim_matrix = raster[stim_idx]

        # plotting
        fig, axs = plt.subplots(2, 1, figsize=(2.1,2.1))
        fig.suptitle(clustr, fontsize=10)

        for line in range(len(ctrl_idx)):
            axs[0].scatter(np.where(ctrl_matrix[line]==1)[0]/samp_freq-3,
                           [line+1]*int(sum(ctrl_matrix[line])),
                           c='grey', ec='none', s=1)
            axs[1].scatter(np.where(stim_matrix[line]==1)[0]/samp_freq-3,
                           [line+1]*int(sum(stim_matrix[line])),
                           c='royalblue', ec='none', s=1)

        for i in range(2):
            axs[i].set(xticks=[0,2,4], xlim=(-1, 4),
                       ylabel='trial #')
            for p in ['top', 'right']:
                axs[i].spines[p].set_visible(False)

        # only set xlabel for ax 1
        axs[1].set(xlabel='time from run onset (s)')

        fig.tight_layout()

        # save to rec folder
        for ext in ['.png', '.pdf']:
            filepath = pp.HPC_EPHYS_FIGURES_STEM / 'all_sessions' / recname / 'rasters_ctrl_stim_pyr' / f'{clustr}{ext}'
            filepath.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(filepath, dpi=300, bbox_inches='tight')

        # stim lines
        for line in range(len(ctrl_idx)):
            if line >= len(pulse_times_aligned) or not np.isfinite(pulse_times_aligned[line]):
                continue
            axs[1].plot([pulse_times_aligned[line]/samp_freq, pulse_times_aligned[line]/samp_freq],
                        [line, line+1], c='red', lw=1)

        # save figure to common folder
        for ext in ['.png', '.pdf']:
            if path in pathHPCLC:
                filepath = pp.HPC_EPHYS_FIGURES_STEM / 'single_cell_ctrl_stim_rasters' / 'HPC_LC_pyr' / f'{clustr}{ext}'
                filepath.parent.mkdir(parents=True, exist_ok=True)
                fig.savefig(filepath, dpi=300, bbox_inches='tight')
            if path in pathHPCLCterm:
                filepath = pp.HPC_EPHYS_FIGURES_STEM / 'single_cell_ctrl_stim_rasters' / 'HPC_LCterm_pyr' / f'{clustr}{ext}'
                filepath.parent.mkdir(parents=True, exist_ok=True)
                fig.savefig(filepath, dpi=300, bbox_inches='tight')

        plt.close(fig)
