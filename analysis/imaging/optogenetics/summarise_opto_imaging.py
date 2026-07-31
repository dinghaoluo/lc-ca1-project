# -*- coding: utf-8 -*-
'''
Created on Thu Sep  4 17:18:06 2025

summarise optogenetics data from imaging animals

@author: Dinghao Luo
'''

#%% imports
import sys
from pathlib import Path
import pickle

import numpy as np

repo_root = Path(__file__).resolve().parents[3]
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))

from common_functions import mpl_formatting
from console_formatting import print_session, print_status
mpl_formatting()

import plotting_functions as pf
import project_paths as pp


#%% recording list
import rec_list

paths = rec_list.pathdLightLCOpto

BEHAVIOUR_STEM = pp.behaviour_experiment_stem('HPCdLightLCOpto')
OUTPUT_STEM = pp.BEHAVIOUR_LC_OPTO_FIGURES_STEM / 'dLight'
OUTPUT_STEM.mkdir(parents=True, exist_ok=True)


#%% 020
ctrl_lick_times_020 = []
stim_lick_times_020 = []
ctrl_lick_distances_020 = []
stim_lick_distances_020 = []

ctrl_mean_speeds_020 = []
stim_mean_speeds_020 = []
ctrl_perc_rew_020 = []
stim_perc_rew_020 = []

for path in paths:
    recname = Path(path).name
    print_session(recname)
    with open(BEHAVIOUR_STEM / f'{recname}.pkl', 'rb') as f:
        curr_txt = pickle.load(f)

    optogenetic_protocol = [t[15] for t in curr_txt['trial_statements']]
    if '2' not in optogenetic_protocol:
        print_status('skipped', 'no optogenetic stim. in this session')
        continue

    # the block starts at the first non-zero protocol and ends before three
    # consecutive zero trials; otherwise it continues through the final trial
    start_idx = next(
        trial_idx
        for trial_idx, protocol in enumerate(optogenetic_protocol)
        if protocol != '0'
    )
    end_idx = len(optogenetic_protocol)
    for trial_idx in range(start_idx, len(optogenetic_protocol) - 2):
        if (
                optogenetic_protocol[trial_idx] == '0'
                and optogenetic_protocol[trial_idx + 1] == '0'
                and optogenetic_protocol[trial_idx + 2] == '0'
                ):
            end_idx = trial_idx
            break

    baseline_trials = list(range(start_idx))
    opto_stim_trials = list(range(start_idx, end_idx, 3))
    opto_stim_trials = [i for i in opto_stim_trials if i < len(optogenetic_protocol) - 1]
    opto_ctrl_trials = [i + 2 for i in opto_stim_trials]
    opto_ctrl_trials = [i for i in opto_ctrl_trials if i < len(optogenetic_protocol) - 1]

    if len(baseline_trials) < 2 or not opto_stim_trials or not opto_ctrl_trials:
        print_status('skipped', 'incomplete trial groups')
        continue

    ctrl_lick_times = []
    stim_lick_times = []
    ctrl_lick_distances = []
    stim_lick_distances = []
    for trial in baseline_trials:
        start_time = np.squeeze(curr_txt['run_onsets'][trial])
        licks = curr_txt['lick_times'][trial]
        ctrl_lick_times.append(
            [(lick[0] - start_time) / 1000 for lick in licks if lick[0] > start_time + 1000]
        )
    for trial in opto_stim_trials:
        start_time = np.squeeze(curr_txt['run_onsets'][trial])
        licks = curr_txt['lick_times'][trial]
        stim_lick_times.append(
            [(lick[0] - start_time) / 1000 for lick in licks if lick[0] > start_time + 1000]
        )
    for trial in baseline_trials:
        ctrl_lick_distances.append(
            [lick for lick in curr_txt['lick_distances_aligned'][trial] if lick > 30]
        )
    for trial in opto_stim_trials:
        stim_lick_distances.append(
            [lick for lick in curr_txt['lick_distances_aligned'][trial] if lick > 30]
        )

    ctrl_first_lick_times = [licks[0] for licks in ctrl_lick_times if licks]
    stim_first_lick_times = [licks[0] for licks in stim_lick_times if licks]
    ctrl_first_lick_distances = [licks[0] for licks in ctrl_lick_distances if licks]
    stim_first_lick_distances = [licks[0] for licks in stim_lick_distances if licks]

    if (
            not ctrl_first_lick_times
            or not stim_first_lick_times
            or not ctrl_first_lick_distances
            or not stim_first_lick_distances
            ):
        print_status('skipped', 'no valid first licks')
        continue

    ctrl_first_lick_time = np.median(ctrl_first_lick_times)
    stim_first_lick_time = np.median(stim_first_lick_times)
    ctrl_first_lick_distance = np.median(ctrl_first_lick_distances)
    stim_first_lick_distance = np.median(stim_first_lick_distances)

    # mean speed
    speed_times_aligned = curr_txt['speed_times_aligned']
    ctrl_mean_speeds = []
    for trial in baseline_trials:
        speed_times = speed_times_aligned[trial]
        if speed_times:
            speeds = [s[1] for s in speed_times]
            ctrl_mean_speeds.append(np.mean(speeds))

    stim_mean_speeds = []
    for trial in opto_stim_trials:
        speed_times = speed_times_aligned[trial]
        if speed_times:
            speeds = [s[1] for s in speed_times]
            stim_mean_speeds.append(np.mean(speeds))

    if not ctrl_mean_speeds or not stim_mean_speeds:
        print_status('skipped', 'no valid speed trials')
        continue

    if 2 < ctrl_first_lick_time < 10 and 2 < stim_first_lick_time < 10:
        ctrl_lick_times_020.append(ctrl_first_lick_time)
        stim_lick_times_020.append(stim_first_lick_time)

    if ctrl_first_lick_distance > 30 and stim_first_lick_distance > 30:
        ctrl_lick_distances_020.append(ctrl_first_lick_distance)
        stim_lick_distances_020.append(stim_first_lick_distance)

    ctrl_mean_speeds_020.append(np.mean(ctrl_mean_speeds))
    stim_mean_speeds_020.append(np.mean(stim_mean_speeds))

    # percent rewarded
    reward_times = curr_txt['reward_times']
    ctrl_rewarded = [not np.isnan(reward_times[trial]) for trial in baseline_trials[1:]]
    ctrl_reward_perc = sum(ctrl_rewarded)/len(baseline_trials[1:])
    stim_rewarded = [not np.isnan(reward_times[trial]) for trial in opto_stim_trials]
    stim_reward_perc = sum(stim_rewarded)/len(opto_stim_trials)

    ctrl_perc_rew_020.append(ctrl_reward_perc)
    stim_perc_rew_020.append(stim_reward_perc)

pf.plot_violin_with_scatter(
    ctrl_lick_times_020, stim_lick_times_020,
    'grey', 'royalblue',
    xticklabels=['ctrl.', 'stim.'],
    ylabel='first-lick time (s)',
    save=True,
    savepath=OUTPUT_STEM / '020_lick_times',
    dpi=300,
    show=False,
    close=True
    )

pf.plot_violin_with_scatter(
    ctrl_lick_distances_020, stim_lick_distances_020,
    'grey', 'royalblue',
    xticklabels=['ctrl.', 'stim.'],
    ylabel='first-lick distance (cm)',
    save=True,
    savepath=OUTPUT_STEM / '020_lick_distances',
    dpi=300,
    show=False,
    close=True
    )

pf.plot_violin_with_scatter(
    ctrl_mean_speeds_020, stim_mean_speeds_020,
    'grey', 'royalblue',
    xticklabels=['ctrl.', 'stim.'],
    ylabel='mean speed (cm/s)',
    save=True,
    savepath=OUTPUT_STEM / '020_mean_speed',
    dpi=300,
    show=False,
    close=True
    )

rmv = []
for i in range(len(ctrl_perc_rew_020)):
    if ctrl_perc_rew_020[i] < 0.5 or stim_perc_rew_020[i] < 0.5:
        rmv.append(i)
ctrl_perc_rew_020 = [trial for i, trial in enumerate(ctrl_perc_rew_020) if i not in rmv]
stim_perc_rew_020 = [trial for i, trial in enumerate(stim_perc_rew_020) if i not in rmv]

pf.plot_violin_with_scatter(
    ctrl_perc_rew_020, stim_perc_rew_020,
    'grey', 'royalblue',
    xticklabels=['ctrl.', 'stim.'],
    ylabel='reward perc.',
    ylim=(.5, 1.02),
    save=True,
    savepath=OUTPUT_STEM / '020_perc_rew',
    dpi=300,
    show=False,
    close=True
    )
