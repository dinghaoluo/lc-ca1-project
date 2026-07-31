# -*- coding: utf-8 -*-
'''
Created on Mon 20 Nov 14:55:04 2023
Originally named ctrl_stim_lick_properties.py
Modified on 11 Sept 2025

quantify lick density etc
modified to plot also the lick distributions (ctrl vs stim)

@author: Dinghao Luo

'''

#%% imports
from pathlib import Path
import sys

import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio

repo_root = Path(__file__).resolve().parents[3]
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))

from plotting_functions import plot_violin_with_scatter
from common_functions import mpl_formatting
from console_formatting import print_session
import project_paths as pp
mpl_formatting()

import rec_list
paths = rec_list.pathLCopt


#%% path stems
exp_stem = pp.MICEEXP_ROOT


#%% parameters
SAMP_FREQ = 1250  # Hz


#%% main
stim_licks = []
ctrl_licks = []

stim_std_med = []; stim_std_mean = []
ctrl_std_med = []; ctrl_std_mean = []

for path in paths:
    recname = Path(path).name
    print_session(recname)

    run_file_path = exp_stem / f'ANMD{recname[1:5]}' / recname[:14] / recname / f'{recname}_DataStructure_mazeSection1_TrialType1_alignRun_msess1.mat'
    alignRun = sio.loadmat(run_file_path)

    licks = alignRun['trialsRun']['lickLfpInd'][0][0][0][1:]
    starts = alignRun['trialsRun']['startLfpInd'][0][0][0][1:]
    pumps = alignRun['trialsRun']['pumpLfpInd'][0][0][0][1:]
    tot_trial = licks.shape[0]
    for trial in range(tot_trial):
        if len(pumps[trial])>0:
            pumps[trial] = pumps[trial][0] - starts[trial]
        else:
            pumps[trial] = [np.nan]

    beh_file_path = exp_stem / f'ANMD{recname[1:5]}' / recname[:14] / recname / f'{recname}_DataStructure_mazeSection1_TrialType1_Info.mat'
    behInfo = sio.loadmat(beh_file_path)['beh']

    stim_idx = np.squeeze(np.where(behInfo['pulseMethod'][0][0][0]!=0))-1
    ctrl_idx = np.arange(stim_idx[0])

    tot_trial = len(behInfo['pulseMethod'][0][0][0])-1

    all_licks = []
    for trial in range(tot_trial):
        # only if the animal does not lick in the first second (carry-over licks) and only include pre-consumption licks
        curr_start = starts[trial]

        curr_licks = [(l[0] - curr_start) / SAMP_FREQ for l in licks[trial]
                      if (l[0] - curr_start) > SAMP_FREQ]
        all_licks.append(curr_licks)

    stim_licks.extend([lk for trial, lk in enumerate(all_licks)
                       if trial in stim_idx and lk])
    ctrl_licks.extend([lk for trial, lk in enumerate(all_licks)
                       if trial in ctrl_idx and lk])

    stim_std = [np.std(lk) for trial, lk in enumerate(all_licks) if trial in stim_idx]
    ctrl_std = [np.std(lk) for trial, lk in enumerate(all_licks) if trial in ctrl_idx]

    stim_std_med.append(np.nanmedian(stim_std))
    ctrl_std_med.append(np.nanmedian(ctrl_std))

    stim_std_mean.append(np.nanmean(stim_std))
    ctrl_std_mean.append(np.nanmean(ctrl_std))


#%% std comparison
savepath = pp.LC_OPTO_EPHYS_FIGURES_STEM / 'opto_licktime_020' / 'lick_std_summary_median'
savepath.parent.mkdir(parents=True, exist_ok=True)
plot_violin_with_scatter(ctrl_std_med, stim_std_med, 'grey', 'royalblue',
                         paired=True,
                         title='lick std. ctrl. v stim.',
                         xticklabels=['ctrl.', 'stim.'],
                         ylabel='lick std.',
                         save=True, savepath=savepath, dpi=300)


#%% lick profiles
plt.figure(figsize=(6,4))

ctrl_all = np.concatenate(ctrl_licks)
stim_all = np.concatenate(stim_licks)

# 50 ms bins over the first 10 s after run onset
bins = np.arange(0, 10, 0.05)

plt.hist(ctrl_all, bins=bins, density=True, histtype='step',
         color='grey', label='ctrl')
plt.hist(stim_all, bins=bins, density=True, histtype='step',
         color='royalblue', label='stim')

plt.xlabel('time from run onset (s)')
plt.ylabel('lick density (a.u.)')
plt.title('mean lick distribution (ctrl vs stim)')
plt.legend(frameon=False)
plt.tight_layout()
plt.show()
