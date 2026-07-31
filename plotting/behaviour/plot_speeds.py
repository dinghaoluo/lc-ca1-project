# -*- coding: utf-8 -*-
'''
Created on Mon 28 Apr 16:37:41 2025

plot the speed profiles

@author: Dinghao Luo

'''

#%% imports
import pickle
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys
import numpy as np
from scipy.stats import sem
from pathlib import Path

repo_root = Path(__file__).resolve().parents[2]
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))
from common_functions import mpl_formatting
import project_paths as pp
mpl_formatting()

import rec_list


#%% load data
exp_name = 'LC'  # HPCLC, HPCLCterm, LC, HPCGRABNE, LCHPCGCaMP
figure_exp_name = exp_name.lower()
output_stem = pp.BEHAVIOUR_FIGURES_STEM / 'speed_profiles'
output_stem.mkdir(parents=True, exist_ok=True)

experiment_paths = {
    'HPCLC': rec_list.pathHPCLCopt,
    'HPCLCterm': rec_list.pathHPCLCtermopt,
    'LC': rec_list.pathLC,
    'HPCGRABNE': rec_list.pathHPCGRABNE,
    'LCHPCGCaMP': rec_list.pathLCHPCGCaMP,
    }

beh_stem = pp.behaviour_experiment_stem(exp_name)
behaviour_sessions = []
for path in experiment_paths[exp_name]:
    recname = Path(path).name
    beh_path = beh_stem / f'{recname}.pkl'
    with open(beh_path, 'rb') as f:
        behaviour_sessions.append(pickle.load(f))


#%% load behav data
trial_speeds = []
trial_speeds_distances = []

for beh in behaviour_sessions:
    run_onsets = beh['run_onsets']
    speed_times = beh['speed_times_aligned']
    speed_distances = beh['speed_distances_aligned']

    for i, run_onset in enumerate(run_onsets):
        if run_onset == -1:
            continue

        trial_speeds.append(
            np.array([s for t, s in speed_times[i]], dtype=float)
            )

        trial_speeds_distances.append(
            np.asarray(speed_distances[i], dtype=float)
            )

# determine max trial length
max_len = max(len(speeds) for speeds in trial_speeds)

# initialise with nan
all_speeds = np.full((len(trial_speeds), max_len), np.nan)

# fill in speed values
for i, speeds in enumerate(trial_speeds):
    all_speeds[i, :len(speeds)] = speeds

# compute mean ignoring nans
mean_speed = np.nanmean(all_speeds, axis=0)
sem_speed = sem(all_speeds, axis=0, nan_policy='omit')

max_distance_len = max(len(speeds) for speeds in trial_speeds_distances)
all_speeds_distances = np.full(
    (len(trial_speeds_distances), max_distance_len),
    np.nan
    )
for i, speeds in enumerate(trial_speeds_distances):
    all_speeds_distances[i, :len(speeds)] = speeds

mean_speed_distances = np.nanmean(all_speeds_distances, axis=0)
sem_speed_distances = sem(all_speeds_distances, axis=0, nan_policy='omit')

# time axis based on sampling rate
sampling_rate = 50  # Hz
dt = 1 / sampling_rate
time_axis = np.arange(max_len) * dt

distance_axis = np.arange(max_distance_len) / 10

#%% main
fig, ax = plt.subplots(figsize=(1.9,1.7))

ax.plot(time_axis, mean_speed, color='navy', lw=1)
ax.fill_between(time_axis, mean_speed+sem_speed,
                           mean_speed-sem_speed,
                color='navy', alpha=.2)

ax.set(xlabel='time from run-onset (s)', xlim=(0, 4), xticks=[0,2,4],
       ylabel='speed (cm/s)', ylim=(0, 70))
for s in ['top','right']:
    ax.spines[s].set_visible(False)

fig.tight_layout()
for ext in ['.png', '.pdf']:
    fig.savefig(
        output_stem / f'speed_profile_time_{figure_exp_name}{ext}',
        dpi=300,
        bbox_inches='tight')
plt.close(fig)

fig, ax = plt.subplots(figsize=(1.9,1.7))

ax.plot(distance_axis, mean_speed_distances, color='navy', lw=1)
ax.fill_between(distance_axis, mean_speed_distances+sem_speed_distances,
                               mean_speed_distances-sem_speed_distances,
                color='navy', alpha=.2)

ax.set(xlim=(0, 200), xlabel='dist. from run-onset (cm)',
       ylabel='speed (cm/s)')
for s in ['top','right']:
    ax.spines[s].set_visible(False)

fig.tight_layout()
for ext in ['.png', '.pdf']:
    fig.savefig(
        output_stem / f'speed_profile_distance_{figure_exp_name}{ext}',
        dpi=300,
        bbox_inches='tight')
plt.close(fig)
