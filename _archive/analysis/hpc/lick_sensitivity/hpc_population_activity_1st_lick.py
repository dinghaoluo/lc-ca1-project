# -*- coding: utf-8 -*-
"""
Created on Thu Aug 22 12:22:43 2024

correlates run-onset population drifts (single-trial) with behavioural
    parameters

@author: Dinghao Luo

"""

#%% imports
import sys
from pathlib import Path

repo_root = Path(__file__).resolve().parent
while not (repo_root / 'utils').exists() and repo_root != repo_root.parent:
    repo_root = repo_root.parent
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))

import scipy.io as sio
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import linregress, zscore


#%% run HPC-LC or HPC-LCterm
HPC_LC = 1

# load paths to recordings
import rec_list
import project_paths as pp
if HPC_LC:
    pathHPC = rec_list.pathHPCLCopt
elif not HPC_LC:
    pathHPC = rec_list.pathHPCLCtermopt
OUTPUT_DIR = pp.HPC_ALL_FIGURES_STEM / 'population_v_1st_licks'


#%% pre-post ratio dataframe
df = pd.read_pickle(pp.HPC_ALL_STEM / 'HPC_LC_stim_stimcont_diff_profiles_pyr_only.pkl')


#%% main
all_pop_peak_z = {}
all_pop_peak_baseline_z = {}
all_first_licks_z = {}
all_first_licks_baseline_z = {}

for pathname in pathHPC:
    recname = pathname[-17:]

    if recname=='A063r-20230708-02' or recname=='A063r-20230708-01':    # lick detection problems
        continue
    print(recname)

    trains = list(np.load(pp.HPC_ALL_STEM / recname / f'HPC_all_info_{recname}.npy',
                          allow_pickle=True).item().values())
    tot_trial = len(trains[0])

    # determine if each cell is pyramidal or intern
    info = sio.loadmat('{}\{}_DataStructure_mazeSection1_TrialType1_Info.mat'.format(pathname, recname))
    rec_info = info['rec'][0][0]
    spike_rate = rec_info['firingRate'][0]
    intern_id = rec_info['isIntern'][0]
    pyr_id = [not(clu) for clu in intern_id]
    tot_clu = len(pyr_id)
    tot_pyr = sum(pyr_id)

    # classify PyrUp/PyrDown neurones
    pyrup = []; pyrdown = []
    for cluname, row in df.iterrows():
        if cluname.split(' ')[0]==recname:
            clu_ID = int(cluname.split(' ')[1][3:])
            if row['ctrl_pre_post']>=1.25:
                pyrdown.append(clu_ID-2)
            if row['ctrl_pre_post']<=.8:
                pyrup.append(clu_ID-2)  # HPCLC_all_train.py adds 2 to the ID, so we subtracts 2 here

    # behaviour parameters
    info = sio.loadmat('{}\{}_DataStructure_mazeSection1_TrialType1_Info.mat'.format(pathname, recname))
    beh_info = info['beh'][0][0]
    behPar = sio.loadmat('{}\{}_DataStructure_mazeSection1_TrialType1_behPar_msess1.mat'.format(pathname, recname))
    stimOn = behPar['behPar']['stimOn'][0][0][0][1:]
    stim_trials = np.where(stimOn!=0)[0]+1
    cont_trials = stim_trials+2
    tot_trial = len(stimOn)
    tot_baseline = stim_trials[0]

    fig, ax = plt.subplots(figsize=(3,2.5)); xaxis = np.arange(1250*5)/1250-1
    pyrup_peak = []
    pyrup_peak_baseline = []
    for trial in range(tot_trial):
        curr_pop_peak = []
        curr_pop_prof = []
        for pyr in pyrup:
            curr_train = trains[pyr][trial]
            if len(curr_train)>3750+1250*1:
                curr_peak = np.nanmean(trains[pyr][trial][3750:3750+1250*1])*1250 - spike_rate[pyr]
            else:
                curr_peak = np.nanmean(trains[pyr][trial][3750:])*1250 - spike_rate[pyr]
            if len(curr_train)>3750+4*1250:  # plotting
                curr_pop_prof.append(curr_train[2500:2500+5*1250])
            curr_pop_peak.append(curr_peak)
        if len(curr_train)>3750+4*1250:
            ax.plot(xaxis, np.mean(curr_pop_prof, axis=0), lw=1, alpha=.1)
        pyrup_peak.append(np.nanmean(curr_pop_peak))
        if trial<stim_trials[0]:
            pyrup_peak_baseline.append(np.nanmean(curr_pop_peak))
    fig.tight_layout()
    plt.show(fig)
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_DIR / f'{recname}_single_trial.png',
                dpi=150, bbox_inches='tight')

    # behaviour
    lickfilename = (
        pp.MICEEXP_ROOT / f'ANMD{recname[1:5]}' / recname[:14] / recname[:17]
        / f'{recname[:17]}_DataStructure_mazeSection1_TrialType1_alignRun_msess1.mat'
        )
    alignRun = sio.loadmat(lickfilename)

    # ignore all 1st trials since it is before counting starts and is an empty cell
    licks = alignRun['trialsRun']['lickLfpInd'][0][0][0][1:]
    starts = alignRun['trialsRun']['startLfpInd'][0][0][0][1:]
    pumps = alignRun['trialsRun']['pumpLfpInd'][0][0][0][1:]

    first_licks = []
    for trial in range(tot_trial):
        lk = [l[0] for l in licks[trial] if l-starts[trial] > 1250]  # exclude licks in the 1st second, as they could be carry-over licks from the last trial
        if len(lk)!=0:  # append only if there is licks in this trial
            first_licks.append((lk[0]-starts[trial])/1250)
        else:
            first_licks.append(0)
    first_licks_baseline = first_licks[:stim_trials[0]]

    # filtering
    delete = []
    for trial in range(tot_trial):
        if first_licks[trial]>10 or first_licks[trial]<=1.5:
            delete.append(trial)
    first_licks = [v for i, v in enumerate(first_licks) if i not in delete]
    first_licks_baseline = [v for i, v in enumerate(first_licks_baseline) if i not in delete]
    pyrup_peak = [v for i, v in enumerate(pyrup_peak) if i not in delete]
    pyrup_peak_baseline = [v for i, v in enumerate(pyrup_peak_baseline) if i not in delete]

    # statistics
    results = linregress(pyrup_peak, first_licks)
    pval = results[3]
    slope = results[0]; intercept = results[1]; rvalue = results[2]
    xmin = min(pyrup_peak); xmax = max(pyrup_peak)
    ymin = min(first_licks); ymax = max(first_licks)

    # plotting
    fig, ax = plt.subplots(figsize=(3,3))
    ax.scatter(pyrup_peak, first_licks, s=2, c='grey')
    ax.plot([xmin, xmax], [intercept+(xmin-.1)*slope, intercept+(xmax+.1)*slope],
            color='k', lw=1)
    ax.set(title='r={}\p={}'.format(round(slope, 5), round(pval, 5)),
           xlabel='PyrUp run-onset modulation',
           ylabel='delay to 1st-licks (s)')
    for s in ['top', 'right']: ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / f'{recname}_all_trias.png',
                dpi=300)
    plt.close()


    # statistics (baseline)
    results_baseline = linregress(pyrup_peak_baseline, first_licks_baseline)
    pval_baseline = results_baseline[3]
    slope_baseline = results_baseline[0]; intercept_baseline = results_baseline[1]; rvalue_baseline = results_baseline[2]
    xmin = min(pyrup_peak_baseline); xmax = max(pyrup_peak_baseline)
    ymin = min(first_licks_baseline); ymax = max(first_licks_baseline)

    # plotting (baseline)
    fig, ax = plt.subplots(figsize=(3,3))
    ax.scatter(pyrup_peak_baseline, first_licks_baseline, s=2, c='grey')
    ax.plot([xmin, xmax], [intercept_baseline+(xmin-.1)*slope_baseline, intercept_baseline+(xmax+.1)*slope_baseline],
            color='k', lw=1)
    ax.set(title='r={}\p={}'.format(round(slope_baseline, 5), round(pval_baseline, 5)),
           xlabel='PyrUp run-onset modulation',
           ylabel='delay to 1st-licks (s)')
    for s in ['top', 'right']: ax.spines[s].set_visible(False)
    fig.tight_layout()
    fig.savefig(OUTPUT_DIR / f'{recname}_baseline.png',
                dpi=300)
    plt.close()

    all_pop_peak_z[recname] = zscore(first_licks)
    all_pop_peak_baseline_z[recname] = zscore(first_licks_baseline)
    all_first_licks_z[recname] = zscore(first_licks)
    all_first_licks_baseline_z[recname] = zscore(first_licks_baseline)
