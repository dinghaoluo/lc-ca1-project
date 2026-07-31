# -*- coding: utf-8 -*-
"""
Created on Thu Sep 21 16:37:16 2023
Modified on Mon 10 Mar 17:45:12 2025:
    modified to work on HPC cells

loop over all cells for early v late trials

@author: Dinghao Luo

"""

#%% imports
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

import scipy.io as sio
import pandas as pd
import sys

repo_root = Path(__file__).resolve().parent
while not (repo_root / 'utils').exists() and repo_root != repo_root.parent:
    repo_root = repo_root.parent
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))
from plotting_functions import plot_violin_with_scatter
from decay_fitting_functions import detect_min_max, compute_tau, plot_fit_compare
import project_paths as pp

import rec_list
# paths = rec_list.pathHPCLCopt + rec_list.pathHPCLCtermopt
paths = rec_list.pathHPCLCopt


#%% parameters
SAMP_FREQ = 1250  # Hz
TIME = np.arange(-SAMP_FREQ, SAMP_FREQ*6)/SAMP_FREQ  # 7 seconds

early_colour = (.804, .267, .267)  # early trials
late_colour = (.545, 0, 0)  # late trials
FIRST_LICK_ANALYSIS_DIR = pp.HPC_EPHYS_STEM / 'first_lick_analysis'
FIRST_LICK_FIGURES_DIR = pp.HPC_EPHYS_FIGURES_STEM / 'first_lick_analysis'
DECAY_COMPARE_DIR = FIRST_LICK_FIGURES_DIR / 'single_cell_decay_constant_early_v_late'
PYRUP_CLASS = 'run-onset ON'
PYRDOWN_CLASS = 'run-onset OFF'
UNRESPONSIVE_CLASS = 'run-onset unresponsive'


def get_savepath(filepath):
    filepath.parent.mkdir(parents=True, exist_ok=True)
    return str(filepath)


#%% load data
print('loading dataframes...')

cell_profiles = pd.read_pickle(
    pp.HPC_EPHYS_STEM / 'HPC_all_profiles.pkl'
    )
df_pyr = cell_profiles[cell_profiles['cell_identity']=='pyr']  # pyramidal only

beh_df = pd.concat((
    pd.read_pickle(
        pp.behaviour_session_pickle('all_HPCLC_sessions.pkl')
        ),
    pd.read_pickle(
        pp.behaviour_session_pickle('all_HPCLCterm_sessions.pkl')
        )
    ))


#%% main loop
early_pyrup_profs = []
late_pyrup_profs = []
early_pyrdown_profs = []
late_pyrdown_profs = []
early_mid_pyrup_profs = []
late_mid_pyrup_profs = []
early_mid_pyrdown_profs = []
late_mid_pyrdown_profs = []

tau_values_early = []
fit_results_early = []
tau_values_late = []
fit_results_late = []

for path in paths:
    recname = path[-17:]

    # get lick times
    curr_beh_df = beh_df.loc[recname]  # subselect in read-only
    run_onsets = curr_beh_df['run_onsets'][1:]
    licks = [
        [(l-run_onset)/1000 for l in trial]  # convert from ms to s
        if len(trial)!=0 else np.nan
        for trial, run_onset in zip(
                curr_beh_df['lick_times'][1:],
                run_onsets
                )
        ]
    first_licks = np.asarray(
        [next((l for l in trial if l > 1), np.nan)  # >1 to prevent carry-over licks
        if isinstance(trial, list) else np.nan
        for trial in licks]
        )

    # get bad trials
    behPar = sio.loadmat(
        pp.MICEEXP_ROOT / f'ANMD{recname[1:5]}' / recname[:14] / recname
        / f'{recname}_DataStructure_mazeSection1_TrialType1_behPar_msess1.mat'
        )
    bad_beh_ind = np.where(behPar['behPar'][0]['indTrBadBeh'][0]==1)[1]-1

    # get early and late lick trials (that are not stim. trials)
    stim_trials = np.where(
        np.asarray([
            trial[15] for trial
            in curr_beh_df['trial_statements']
            ])!='0'
        )[0]
    valid_trials = [i for i in range(len(first_licks))
                    if i not in stim_trials
                    and i not in bad_beh_ind
                    and not np.isnan(first_licks[i])]

    if len(valid_trials) < 50:
        continue

    # load spike trains
    print(f'\n{recname}')
    trains = np.load(
        pp.HPC_EPHYS_STEM / 'all_sessions' / recname / f'{recname}_all_trains.npy',
        allow_pickle=True
        ).item()

    sorted_trials = sorted(valid_trials,
                           key=lambda i: first_licks[i])[10:-10]  # avoid extremities

    early_trials = sorted_trials[:10]
    early_mid_trials = sorted_trials[:int(len(sorted_trials)/2)]
    late_trials = sorted_trials[-10:]
    late_mid_trials = sorted_trials[int(len(sorted_trials)/2):]

    # get cell spiking data
    curr_df_pyr = df_pyr[df_pyr['recname']==recname]

    for cluname, session in curr_df_pyr.iterrows():
        train = trains[cluname]
        early_trains = [
            train[trial, :] for trial
            in early_trials
            ]
        late_trains = [
            train[trial, :] for trial
            in late_trials
            ]
        early_mid_trains = [
            train[trial, :] for trial
            in early_mid_trials
            ]
        late_mid_trains = [
            train[trial, :] for trial
            in late_mid_trials
            ]

        early_mean = np.mean(early_trains, axis=0)
        late_mean = np.mean(late_trains, axis=0)

        # PyrUp and PyrDown only
        if session['class'] == PYRUP_CLASS:
            early_pyrup_profs.append(early_mean)
            late_pyrup_profs.append(late_mean)
            early_mid_pyrup_profs.append(np.mean(early_mid_trains, axis=0))
            late_mid_pyrup_profs.append(np.mean(late_mid_trains, axis=0))
        elif session['class'] == PYRDOWN_CLASS:
            early_pyrdown_profs.append(early_mean)
            late_pyrdown_profs.append(late_mean)
            early_mid_pyrdown_profs.append(np.mean(early_mid_trains, axis=0))
            late_mid_pyrdown_profs.append(np.mean(late_mid_trains, axis=0))

        # for decay time calculation
        if session['class'] != UNRESPONSIVE_CLASS:
            mean_prof_early = early_mean[SAMP_FREQ*(3-1):SAMP_FREQ*(3+6)]
            mean_prof_late = late_mean[SAMP_FREQ*(3-1):SAMP_FREQ*(3+6)]

            peak_idx_early = detect_min_max(mean_prof_early,
                                            session['class'],
                                            run_onset_bin=1250)
            peak_idx_late = detect_min_max(mean_prof_late,
                                           session['class'],
                                           run_onset_bin=1250)

            tau_early, fit_params_early = compute_tau(
                TIME, mean_prof_early, peak_idx_early, session['class']
                )
            tau_late, fit_params_late = compute_tau(
                TIME, mean_prof_late, peak_idx_late, session['class']
                )

            # early and late plot
            plot_fit_compare(TIME,
                             mean_prof_early, peak_idx_early, fit_params_early,
                             mean_prof_late, peak_idx_late, fit_params_late,
                             cluname, session['class'],
                             filename=DECAY_COMPARE_DIR / f"{cluname} {session['class']}.png")

            tau_values_early.append(tau_early)
            fit_results_early.append(fit_params_early)
            tau_values_late.append(tau_late)
            fit_results_late.append(fit_params_late)


#%% save for further processing
FIRST_LICK_ANALYSIS_DIR.mkdir(parents=True, exist_ok=True)
np.save(FIRST_LICK_ANALYSIS_DIR / 'early_PyrUp_profs.npy',
        early_pyrup_profs,
        allow_pickle=True)
np.save(FIRST_LICK_ANALYSIS_DIR / 'early_PyrDown_profs.npy',
        early_pyrdown_profs,
        allow_pickle=True)
np.save(FIRST_LICK_ANALYSIS_DIR / 'late_PyrUp_profs.npy',
        late_pyrup_profs,
        allow_pickle=True)
np.save(FIRST_LICK_ANALYSIS_DIR / 'late_PyrDown_profs.npy',
        late_pyrdown_profs,
        allow_pickle=True)


#%% compute mean and sem
early_pyrup_profs_mean = np.mean(early_pyrup_profs, axis=0)
late_pyrup_profs_mean = np.mean(late_pyrup_profs, axis=0)
early_pyrdown_profs_mean = np.mean(early_pyrdown_profs, axis=0)
late_pyrdown_profs_mean = np.mean(late_pyrdown_profs, axis=0)
early_mid_pyrup_profs_mean = np.mean(early_mid_pyrup_profs, axis=0)
late_mid_pyrup_profs_mean = np.mean(late_mid_pyrup_profs, axis=0)
early_mid_pyrdown_profs_mean = np.mean(early_mid_pyrdown_profs, axis=0)
late_mid_pyrdown_profs_mean = np.mean(late_mid_pyrdown_profs, axis=0)

from scipy.stats import sem
early_pyrup_profs_sem = sem(early_pyrup_profs, axis=0)
late_pyrup_profs_sem = sem(late_pyrup_profs, axis=0)
early_pyrdown_profs_sem = sem(early_pyrdown_profs, axis=0)
late_pyrdown_profs_sem = sem(late_pyrdown_profs, axis=0)
early_mid_pyrup_profs_sem = sem(early_mid_pyrup_profs, axis=0)
late_mid_pyrup_profs_sem = sem(late_mid_pyrup_profs, axis=0)
early_mid_pyrdown_profs_sem = sem(early_mid_pyrdown_profs, axis=0)
late_mid_pyrdown_profs_sem = sem(late_mid_pyrdown_profs, axis=0)

early_pyrup_ratios = [np.nanmean(prof[3750-1250:3750])/np.nanmean(prof[3750:3750+1250])
                      if sum(prof[3750:3750+1250])>0 else 1
                      for prof in early_pyrup_profs]
late_pyrup_ratios = [np.nanmean(prof[3750-1250:3750])/np.nanmean(prof[3750:3750+1250])
                     if sum(prof[3750:3750+1250])>0 else 1
                     for prof in late_pyrup_profs]

outlier_mask = [i for i in range(len(early_pyrup_ratios))
                if early_pyrup_ratios[i] > 10 or late_pyrup_ratios[i] > 10]


early_pyrup_ratios = [v for i, v in enumerate(early_pyrup_ratios)
                      if i not in outlier_mask]
late_pyrup_ratios = [v for i, v in enumerate(late_pyrup_ratios)
                     if i not in outlier_mask]

early_mid_pyrup_ratios = [np.nanmean(prof[3750-1250:3750])/np.nanmean(prof[3750:3750+1250])
                          if sum(prof[3750:3750+1250])>0 else 1
                          for prof in early_mid_pyrup_profs]
late_mid_pyrup_ratios = [np.nanmean(prof[3750-1250:3750])/np.nanmean(prof[3750:3750+1250])
                         if sum(prof[3750:3750+1250])>0 else 1
                         for prof in late_mid_pyrup_profs]

outlier_mid_mask = [i for i in range(len(early_mid_pyrup_ratios))
                    if early_mid_pyrup_ratios[i] > 5 or late_mid_pyrup_ratios[i] > 5]

early_mid_pyrup_ratios = [v for i, v in enumerate(early_mid_pyrup_ratios)
                          if i not in outlier_mid_mask]
late_mid_pyrup_ratios = [v for i, v in enumerate(late_mid_pyrup_ratios)
                         if i not in outlier_mid_mask]

plot_violin_with_scatter(early_mid_pyrup_ratios, late_mid_pyrup_ratios, 'orange', 'darkred')


#%% decay time analysis
tau_values_early, tau_values_late = zip(
    *[(x, y) for x, y in zip(tau_values_early, tau_values_late)
      if x is not None and y is not None]
    )

tau_values_early_pyrup, tau_values_late_pyrup = zip(
    *[(x, y) for x, y in zip(tau_values_early, tau_values_late)
      if 0 < x < 5 and 0 < y < 5]
    )
tau_values_early_pyrdown, tau_values_late_pyrdown = zip(
    *[(x, y) for x, y in zip(tau_values_early, tau_values_late)
      if x < 0 and y < 0]
    )

plot_violin_with_scatter(tau_values_early_pyrup, tau_values_late_pyrup,
                         'lightcoral', 'firebrick',
                         xticklabels=['early\n$1^{st}$-lick', 'late\n$1^{st}$-lick'],
                         ylabel='τ (s)',
                         title='PyrUp',
                         showmeans=True,
                         showmedians=False,
                         showscatter=False,
                         ylim=(0,5),
                         save=True,
                         savepath=get_savepath(FIRST_LICK_FIGURES_DIR / 'PyrUp_decay_constant'))

plot_violin_with_scatter(tau_values_early_pyrdown, tau_values_late_pyrdown,
                         'thistle', 'purple',
                         xticklabels=['early\n$1^{st}$-lick', 'late\n$1^{st}$-lick'],
                         ylabel='τ (s)',
                         title='PyrDown',
                         showmeans=True,
                         showmedians=False,
                         showscatter=False,
                         save=True,
                         savepath=get_savepath(FIRST_LICK_FIGURES_DIR / 'PyrDown_decay_constant'))


#%% plotting
xaxis = np.arange(-1250, 1250*4)/1250

fig, ax = plt.subplots(figsize=(3,2))

ax.plot(xaxis,
        early_pyrup_profs_mean[3750-1250:3750+1250*4])
ax.plot(xaxis,
        late_pyrup_profs_mean[3750-1250:3750+1250*4],
        color='red')
ax.fill_between(
    xaxis,
    early_pyrup_profs_mean[3750-1250:3750+1250*4]+early_pyrup_profs_sem[3750-1250:3750+1250*4],
    early_pyrup_profs_mean[3750-1250:3750+1250*4]-early_pyrup_profs_sem[3750-1250:3750+1250*4],
    alpha=.35
    )
ax.fill_between(
    xaxis,
    late_pyrup_profs_mean[3750-1250:3750+1250*4]+late_pyrup_profs_sem[3750-1250:3750+1250*4],
    late_pyrup_profs_mean[3750-1250:3750+1250*4]-late_pyrup_profs_sem[3750-1250:3750+1250*4],
    alpha=.35, color='red', edgecolor='none'
    )
