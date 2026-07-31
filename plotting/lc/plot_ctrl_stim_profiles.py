# -*- coding: utf-8 -*-
'''
Created on Fri Feb 28 11:30:26 2025

plot ctrl vs stim profiles for LC cells

@author: Dinghao Luo
'''

#%% imports
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
from pathlib import Path

repo_root = Path(__file__).resolve().parents[2]
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))
from common_functions import mpl_formatting
import project_paths as pp
mpl_formatting()


#%% main
cell_profiles = pd.read_pickle(
    pp.LC_EPHYS_STEM / 'LC_all_cell_profiles.pkl'
    )
single_cell_stim_stem = pp.LC_EPHYS_FIGURES_STEM / 'single_cell_stim_response'
single_cell_stim_stem.mkdir(parents=True, exist_ok=True)

samp_freq = 1250  # Hz
xaxis = np.arange(samp_freq*5) / samp_freq - 1  # -1~4
for clu in cell_profiles.itertuples():
    cluname = clu.Index

    stim_mean = clu.stim_mean
    stim_sem = clu.stim_sem
    ctrl_mean = clu.ctrl_mean
    ctrl_sem = clu.ctrl_sem

    fig, ax = plt.subplots(figsize=(2.55,1.8))
    ctrlln, = ax.plot(xaxis, ctrl_mean[3750-1250:3750+1250*4],
                      color='grey')
    ax.fill_between(xaxis, (ctrl_mean+ctrl_sem)[3750-1250:3750+1250*4],
                           (ctrl_mean-ctrl_sem)[3750-1250:3750+1250*4],
                    color='grey', edgecolor='none', alpha=.35)
    stimln, = ax.plot(xaxis, stim_mean[3750-1250:3750+1250*4],
                      color='royalblue')
    ax.fill_between(xaxis, (stim_mean+stim_sem)[3750-1250:3750+1250*4],
                           (stim_mean-stim_sem)[3750-1250:3750+1250*4],
                    color='royalblue', edgecolor='none', alpha=.35)

    ax.set(title=cluname,
           xlabel='time from run-onset (s)',
           ylabel='spike rate (Hz)')

    ax.legend((ctrlln, stimln), ('ctrl.', 'stim.'),
              frameon=False)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)

    fig.savefig(
        single_cell_stim_stem / f'{cluname}.png',
        dpi=300,
        bbox_inches='tight'
        )
