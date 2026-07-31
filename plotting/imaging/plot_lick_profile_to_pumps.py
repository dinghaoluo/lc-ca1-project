# -*- coding: utf-8 -*-
'''
Created on Tue Sep 10 17:46:46 2024
Modified on Fri  Sept 20 15:15:12 2024 to plot lick to pump distribution

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

# updated to eliminate nested for loops, 3 Oct 2024 Dinghao
lick_to_pumps = []
for payload in behaviour.values():
    lick_series = payload['lick_times']
    pump_series = payload['reward_times']
    tot_trials = min(len(lick_series), len(pump_series))
    for trial in range(tot_trials):
        pump = pump_series[trial]
        if not np.isfinite(pump) or not lick_series[trial]:
            continue
        for lick in lick_series[trial]:
            lick_to_pumps.append((lick-pump)/1250)


#%% main
fig, ax = plt.subplots(figsize=(2,1.7))

# create histogram
ax.hist(lick_to_pumps, bins=60, range=(-5, 1), density=True, color='orchid')

ax.set(xlim=(-5, 1), xlabel='time to reward (s)',
       yticks=[0, 0.5], ylabel='histogram of licks')
for s in ['top','right']: ax.spines[s].set_visible(False)

fig.tight_layout()
for ext in ['png', 'pdf']:
    fig.savefig(output_stem / f'lick_to_pumps_GRABNE.{ext}',
                dpi=300,
                bbox_inches='tight')
plt.close(fig)
