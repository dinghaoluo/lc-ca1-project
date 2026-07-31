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
import scipy.io as sio
import pandas as pd
import sys
from pathlib import Path
import matplotlib.cm as cm

repo_root = Path(__file__).resolve().parent
while not (repo_root / 'utils').exists() and repo_root != repo_root.parent:
    repo_root = repo_root.parent
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))
from plotting_functions import plot_violin_with_scatter
import project_paths as pp

import rec_list
paths = rec_list.pathHPCLCopt + rec_list.pathHPCLCtermopt

# pre_post ratio thresholds
run_onset_activated_thres = 0.80
run_onset_inhibited_thres = 1.25
PYRUP_CLASS = 'run-onset ON'
PYRDOWN_CLASS = 'run-onset OFF'
UNRESPONSIVE_CLASS = 'run-onset unresponsive'


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


#%% functions
def classify_run_onset_activation_ratio(train,
                                        run_onset_activated_thres,
                                        run_onset_inhibited_thres):
    """
    classify run-onset activation ratio based on pre- and post-run periods.

    parameters:
    - train: array of firing rates over time.
    - run_onset_activated_thres: threshold for classifying activation.
    - run_onset_inhibited_thres: threshold for classifying inhibition.
    - samp_freq: sampling frequency in Hz, default is 1250.
    - run_onset_bin: bin marking the run onset, default is 3750.

    returns:
    - ratio: pre/post activation ratio.
    - response_class: string indicating the activation class (PyrUp, PyrDown, unresponsive).
    """
    samp_freq = 1250
    run_onset_bin = 3750

    pre = np.nanmean(train[int(run_onset_bin-samp_freq*1.5):int(run_onset_bin-samp_freq*.5)])
    post = np.nanmean(train[int(run_onset_bin+samp_freq*.5):int(run_onset_bin+samp_freq*1.5)])
    ratio = pre/post
    if ratio < run_onset_activated_thres:
        response_class = PYRUP_CLASS
    elif ratio > run_onset_inhibited_thres:
        response_class = PYRDOWN_CLASS
    else:
        response_class = UNRESPONSIVE_CLASS

    return ratio, response_class


#%% main loop
early_pyrup_props = []
late_pyrup_props = []
early_pyrdown_props = []
late_pyrdown_props = []

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
    late_trials = sorted_trials[-10:]

    # get cell spiking data
    curr_df_pyr = df_pyr[df_pyr['recname']==recname]

    early_pyrup_count = 0
    late_pyrup_count = 0
    early_pyrdown_count = 0
    late_pyrdown_count = 0

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
        early_prof = np.mean(early_trains, axis=0)
        late_prof = np.mean(late_trains, axis=0)

        early_ratio, early_response_class = classify_run_onset_activation_ratio(
            early_prof,
            run_onset_activated_thres,
            run_onset_inhibited_thres
            )
        late_ratio, late_response_class = classify_run_onset_activation_ratio(
            late_prof,
            run_onset_activated_thres,
            run_onset_inhibited_thres
            )

        response_class = session['class']
        if early_response_class == PYRUP_CLASS:
            early_pyrup_count+=1
        elif early_response_class == PYRDOWN_CLASS:
            early_pyrdown_count+=1
        if late_response_class == PYRUP_CLASS:
            late_pyrup_count+=1
        elif late_response_class == PYRDOWN_CLASS:
            late_pyrdown_count+=1

        # plot_colours = [cm.Reds(i / (len(sorted_trials)-1))
        #                 for i in range(len(sorted_trials))]
        # if response_class != UNRESPONSIVE_CLASS:
        #     fig, ax = plt.subplots(figsize=(3,2))
        #     for i, trial in enumerate(sorted_trials):
        #         ax.plot(np.arange(-1250, 1250*4)/1250,
        #                 train[trial, 3750-1250:3750+4*1250],
        #                 color=plot_colours[i],
        #                 linewidth=1)
        #     ax.set(xlabel='time from run-onset (s)',
        #            ylabel='spike rate (Hz)',
        #            title=f'{cluname}\n{response_class}')

        #     outdir = pp.HPC_EPHYS_STEM / 'first_lick_alignment' / 'single_cell_PyrUp_PyrDown'
        #     fig.savefig(outdir / f'{cluname} {response_class}.png',
        #                 dpi=300,
        #                 bbox_inches='tight')

    early_pyrup_props.append(early_pyrup_count/len(curr_df_pyr))
    early_pyrdown_props.append(early_pyrdown_count/len(curr_df_pyr))
    late_pyrup_props.append(late_pyrup_count/len(curr_df_pyr))
    late_pyrdown_props.append(late_pyrdown_count/len(curr_df_pyr))


#%% compute results
plot_violin_with_scatter(early_pyrup_props, late_pyrup_props, 'orange', 'darkred')
