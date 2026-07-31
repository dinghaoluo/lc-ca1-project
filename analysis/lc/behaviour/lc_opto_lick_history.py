# -*- coding: utf-8 -*-
'''
Created on Fri 18 Aug 17:41:33 2023
Originally named lick_history_dependency.py

analyse the history dependency of 1st lick timing in opto-sessions (all)

@author: Dinghao Luo

'''

#%% imports
import numpy as np
import matplotlib.pyplot as plt; plt.rcParams['font.family'] = 'Arial'
import scipy.io as sio
from scipy.stats import linregress  # median used
import sys
from pathlib import Path

repo_root = Path(__file__).resolve().parents[3]
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))

from console_formatting import print_session
import rec_list
import project_paths as pp
pathOpt = rec_list.pathLCopt

# 0-2-0
sess_list = [sess[-17:] for sess in pathOpt]

# n_bst = 1000  # hyperparameter for bootstrapping


#%% main
for sessname in sess_list:
    print_session(sessname)

    rec_stem = pp.MICEEXP_ROOT / f'ANMD{sessname[1:5]}' / sessname[:14] / sessname[:17]
    # infofilename = rec_stem / f'{sessname[:17]}_DataStructure_mazeSection1_TrialType1_Info.mat'

    # Info = sio.loadmat(infofilename)
    # pulseMethod = Info['beh'][0][0]['pulseMethod'][0]

    # # stim info
    # tot_stims = len([t for t in pulseMethod if t!=0])
    # stim_cond = pulseMethod[np.where(pulseMethod!=0)][0]  # check stim condition
    # stim = [i for i, e in enumerate(pulseMethod) if e==stim_cond]

    # licks
    lickfilename = rec_stem / f'{sessname[:17]}_DataStructure_mazeSection1_TrialType1_alignRun_msess1.mat'
    alignRun = sio.loadmat(lickfilename)

    # ignore all 1st trials since it is before counting starts and is an empty cell
    licks = alignRun['trialsRun']['lickLfpInd'][0][0][0][1:]
    starts = alignRun['trialsRun']['startLfpInd'][0][0][0][1:]

    tot_trial = licks.shape[0]
    first_licks = []

    for trial in range(tot_trial):
        lk = [l[0] for l in licks[trial] if l-starts[trial] > 1250]  # exclude licks in the 1st second, as they could be carry-over licks from the last trial
        if len(lk)!=0:  # append only if there is licks in this trial
            first_licks.append((lk[0]-starts[trial])/1250)
        else:
            first_licks.append(0)

    filt_curr_trials = []
    filt_prev_trials = []
    # filtering to remove outliers using 3 standard deviations
    mean_licks = np.mean(first_licks)
    std_licks = np.std(first_licks)
    upper = mean_licks + 2*std_licks
    lower = mean_licks - 2*std_licks
    if lower<=0:  # to make sure that trials with no licks (see above) do not get in the dataset
        lower = 0
    for i in range(1, len(first_licks)):
        if first_licks[i]<upper and first_licks[i]>lower and first_licks[i-1]<upper and first_licks[i-1]>lower:
            filt_curr_trials.append(first_licks[i])
            filt_prev_trials.append(first_licks[i-1])

    # linear regression comparing delta in lick times curr. v prev.
    filt_delta = [filt_curr_trials[i]-filt_prev_trials[i] for i in range(len(filt_curr_trials))]
    results = linregress(filt_prev_trials, filt_delta)
    pval = results[3]
    slope = results[0]; intercept = results[1]
    xmin = min(filt_prev_trials); xmax = max(filt_prev_trials)
    ymin = min(filt_delta); ymax = max(filt_delta)

    fig, ax = plt.subplots(figsize=(4,2.5))
    ax.set(xlabel='prev. trial t. 1st lick (s)',
           ylabel='change in t. 1st-lick\non curr. trial (s)',
           title='{}\n1st-lick time history dependency\npval={}, slope={}'.format(sessname, np.round(pval,4), np.round(slope,4)),
           xlim=(xmin-.5, xmax+.5), ylim=(ymin-.5, ymax+.5))
    ax.scatter(filt_prev_trials, filt_delta, color='grey', s=3)
    ax.plot([xmin-.1, xmax+.1], [intercept+(xmin-.1)*slope, intercept+(xmax+.1)*slope], color='k')
    plt.show()

    output_stem = pp.LC_OPTO_EPHYS_FIGURES_STEM / 'history_dependency'
    output_stem.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_stem / f'{sessname}.png',
                dpi=300,
                bbox_inches='tight')
