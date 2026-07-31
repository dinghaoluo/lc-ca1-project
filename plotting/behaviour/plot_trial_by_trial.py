# -*- coding: utf-8 -*-
'''
Created on Thu Mar  2 14:33:49 2023

plot run bouts for the selected example behaviour session

@author: Dinghao Luo
'''

#%% imports
from pathlib import Path
import pickle
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

repo_root = Path(__file__).resolve().parents[2]
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))

from common_functions import mpl_formatting
from console_formatting import print_session
import project_paths as pp
mpl_formatting()


#%% parameters
EXAMPLE_EXP_NAME = 'LCHPCGCaMP'
EXAMPLE_RECNAME = 'A106i-20250127-02'  # original rec_list.pathLCHPCGCaMP[43:44]

PRE_WINDOW_MS = 1000
POST_WINDOW_MS = 1000
TRIAL_STEP = 3


#%% load data
recname = EXAMPLE_RECNAME
print_session(recname)

beh_path = pp.behaviour_experiment_stem(EXAMPLE_EXP_NAME) / f'{recname}.pkl'
with open(beh_path, 'rb') as f:
    data = pickle.load(f)

run_bout_table = pd.read_csv(
    pp.RUN_BOUTS_STEM / f'{recname}_run_bouts_py.csv'
    )

output_dir = (
    pp.BEHAVIOUR_FIGURES_STEM
    / 'trial_profiles'
    / pp.behaviour_experiment_folder(EXAMPLE_EXP_NAME)
    / recname
    )
output_dir.mkdir(parents=True, exist_ok=True)

speed_times_ms = np.asarray(data['upsampled_timestamps_ms'], dtype=float)
speed_cm_s = np.asarray(data['upsampled_speed_cm_s'], dtype=float)

trial_end_times = np.asarray(
    [float(statement[1]) for statement in data['new_trial_statements']],
    dtype=float,
)
run_onsets = np.asarray(data['run_onsets'], dtype=float)
cue_times = np.asarray(
    [cue for trial_cues in data['start_cue_times'] for cue in trial_cues],
    dtype=float,
)
reward_times = np.asarray(data['reward_times'], dtype=float)
lick_times = np.asarray(
    [lick for trial_licks in data['lick_times'] for lick in trial_licks],
    dtype=float,
)
run_bout_starts = pd.to_numeric(
    run_bout_table['run_start_time'], errors='coerce'
    ).dropna().to_numpy(dtype=float)
# run-bout times are in seconds; behaviour event times are in milliseconds
if np.nanmax(run_bout_starts) < np.nanmax(speed_times_ms) / 10:
    run_bout_starts = run_bout_starts * 1000


#%% plotting
n_trials = len(run_onsets)
print('plotting...')
for t in np.arange(2, n_trials - 2, TRIAL_STEP):
    window_anchor_ms = run_onsets[t]
    window_end_ms = trial_end_times[t + 2]

    window_start_ms = window_anchor_ms - PRE_WINDOW_MS
    window_end_ms = window_end_ms + POST_WINDOW_MS
    speed_mask = (
        (speed_times_ms >= window_start_ms)
        & (speed_times_ms <= window_end_ms)
        )

    xaxis = (speed_times_ms[speed_mask] - window_start_ms) / 1000
    speed_window = speed_cm_s[speed_mask]
    speed_ylim = np.nanmax(speed_window) * 1.2

    fig_width = len(speed_window) / 5000
    fig, ax = plt.subplots(figsize=(fig_width, 1.2))
    ax.set(
        xlabel='time (s)',
        ylabel='speed (cm/s)',
        ylim=(0, speed_ylim),
        xlim=(0, xaxis[-1]),
        title=f'{recname} trials {t-1} to {t+1}',
        )
    ax.plot(xaxis, speed_window, color='royalblue', label='speed')

    cue_t = cue_times[
        (cue_times >= window_start_ms) & (cue_times <= window_end_ms)
        ]
    ax.vlines(
        (cue_t - window_start_ms) / 1000,
        0,
        speed_ylim,
        'darkgrey',
        )

    reward_t = reward_times[
        (reward_times >= window_start_ms) & (reward_times <= window_end_ms)
        ]
    ax.vlines(
        (reward_t - window_start_ms) / 1000,
        speed_ylim,
        speed_ylim * 0.95,
        'forestgreen',
        linewidth=1.5,
        zorder=10,
        )

    run_onset_t = run_onsets[
        (run_onsets >= window_start_ms) & (run_onsets <= window_end_ms)
        ]
    ax.vlines(
        (run_onset_t - window_start_ms) / 1000,
        0,
        speed_ylim,
        'red',
        linestyle='dashed',
        )

    run_bout_t = run_bout_starts[
        (run_bout_starts >= window_start_ms)
        & (run_bout_starts <= window_end_ms)
        ]
    ax.vlines(
        (run_bout_t - window_start_ms) / 1000,
        0,
        speed_ylim,
        'green',
        linestyle='dashed',
        zorder=10,
        )

    licks_t = lick_times[
        (lick_times >= window_start_ms) & (lick_times <= window_end_ms)
        ]
    ax.vlines(
        (licks_t - window_start_ms) / 1000,
        speed_ylim * 0.96,
        speed_ylim * 0.99,
        'magenta',
        linewidth=.8,
        )

    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)

    rb_tag = ' rb' if len(run_bout_t) > 0 else ''
    save_path = output_dir / f'trials_{t-1}_to_{t+1}{rb_tag}'
    for ext in ['.pdf', '.png']:
        fig.savefig(f'{save_path}{ext}', dpi=300, bbox_inches='tight')
    plt.close(fig)
