# -*- coding: utf-8 -*-
'''
Created on Tue Sep 10 17:46:46 2024

plot the lick profile of GRABNE sessions

@author: Dinghao Luo
'''

#%% imports
from pathlib import Path
import pickle
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

repo_root = Path(__file__).resolve().parents[2]
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))

import project_paths as pp
import rec_list


#%% load behav data
behaviour = {}
behaviour_stem = pp.behaviour_experiment_stem('HPCGRABNE')
for path in rec_list.pathHPCGRABNE:
    recname = Path(path).name
    with open(behaviour_stem / f'{recname}.pkl', 'rb') as f:
        behaviour[recname] = pickle.load(f)

output_stem = pp.BEHAVIOUR_FIGURES_STEM
output_stem.mkdir(parents=True, exist_ok=True)

licks = []
for payload in behaviour.values():
    for lick_distances in payload['lick_distances_aligned']:
        lick_distances = np.asarray(lick_distances, dtype=float).ravel()
        licks.extend(lick_distances[np.isfinite(lick_distances)])


#%% main
fig, ax = plt.subplots(figsize=(2.4,2))

ax.hist(licks, bins=100, range=(0, 220), density=True, color='orchid')

ax.set(xlim=(30,220), xlabel='distance (cm)',
       ylim=(0, 0.05), yticks=[0, 0.04], ylabel='histogram of licks')
for s in ['top','right']: ax.spines[s].set_visible(False)

fig.tight_layout()
fig.savefig(output_stem / 'lick_profile_GRABNE.png',
            dpi=300,
            bbox_inches='tight')
fig.savefig(output_stem / 'lick_profile_GRABNE.pdf',
            bbox_inches='tight')
plt.close(fig)

#%% same to time
licks = []
for payload in behaviour.values():
    lick_series = payload['lick_times']
    starts = payload['run_onsets']
    tot_trials = min(len(lick_series), len(starts))
    for trial in range(tot_trials):
        start = starts[trial]
        if not np.isfinite(start):
            continue
        for t in lick_series[trial]:
            licks.append((t[0]-start)/10)


#%% main
fig, ax = plt.subplots(figsize=(2.4,2))

ax.hist(licks, bins=100, range=(0,1000), density=True, color='orchid')

ax.set(xlabel='time (s)',
       ylabel='histogram of licks')
for s in ['top','right']: ax.spines[s].set_visible(False)

fig.tight_layout()
fig.savefig(output_stem / 'lick_profile_GRABNE_time.png',
            dpi=300,
            bbox_inches='tight')
fig.savefig(output_stem / 'lick_profile_GRABNE_time.pdf',
            bbox_inches='tight')
plt.close(fig)
