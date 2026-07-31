# -*- coding: utf-8 -*-
'''
Created on Mon Oct 13 12:03:43 2025
Modified on Friday to get reward-aligned firing profiles
Duplicated and modified on Thursday 8 Jan 2026 to
    correlate w first lick time (n+1 trial)

Correlation between time since last reward and first-lick time

@author: Dinghao Luo
'''

#%% imports
from pathlib import Path
import sys

import numpy as np
import pickle

repo_root = Path(__file__).resolve().parents[2]
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))

from common_functions import mpl_formatting
from console_formatting import print_files_saved, print_session, print_statistics_section
import glm_functions as gf
import project_paths as pp
mpl_formatting()

import rec_list
paths = rec_list.pathLC


#%% paths and parameters
LC_beh_stem = pp.behaviour_experiment_stem('LC')
output_path = (
    pp.BEHAVIOUR_STEM
    / 'first_lick_since_last_reward'
    / 'time_since_reward_vs_first_lick_payload.npy'
    )
SAMP_FREQ_BEH = 1000


#%% analysis
all_t_since = []
all_first_licks = []
animals = set()
n_sess = 0

for path in paths:
    recname = Path(path).name
    print_session(recname)
    n_sess += 1
    animals.add(recname.split('-')[0])

    with open(LC_beh_stem / f'{recname}.pkl', 'rb') as f:
        beh = pickle.load(f)

    lick_times = beh['lick_times_aligned'][1:]
    reward_times = beh['reward_times'][1:]
    run_onsets = beh['run_onsets'][1:]
    trial_statements = beh['trial_statements'][1:]

    opto_idx = [
        i for i, trial in enumerate(trial_statements)
        if trial[15] != '0'
    ]
    valid_trials = [
        trial for trial, run_onset in enumerate(run_onsets[:-1])
        if (
            trial not in opto_idx
            and trial - 1 not in opto_idx
            and not np.isnan(run_onset)
        )
    ]

    for trial_idx in valid_trials:
        onset_time = run_onsets[trial_idx] / SAMP_FREQ_BEH
        time_since_reward = gf.time_since_last_reward(
            reward_times,
            onset_time,
            trial_idx
        )
        if np.isnan(time_since_reward) or not 0 <= time_since_reward <= 8:
            continue

        curr_licks = lick_times[trial_idx]
        if not isinstance(curr_licks, list) or len(curr_licks) <= 1:
            continue

        first_lick = curr_licks[0] / SAMP_FREQ_BEH
        if not 0.5 <= first_lick <= 8:
            continue

        all_t_since.append(time_since_reward)
        all_first_licks.append(first_lick)

x = np.asarray(all_t_since, dtype=float)
y = np.asarray(all_first_licks, dtype=float)
print_statistics_section()
print(f'n = {len(animals)}')
print(f'n_sess = {n_sess}')

output_path.parent.mkdir(parents=True, exist_ok=True)
np.save(
    output_path,
    {
        'x': x,
        'y': y,
        'animals': sorted(animals),
        'n_sess': n_sess
        }
    )
print_files_saved([
    ('first-lick data', output_path),
])
