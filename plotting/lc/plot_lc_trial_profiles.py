# -*- coding: utf-8 -*-
'''
Created on Tue Oct  1 17:45:35 2024
Originally named plot_trials_lc.py

plot trial-profiles for LC animals

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

# plotting parameters
plt.rcParams['font.family'] = 'Arial'
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

repo_root = Path(__file__).resolve().parents[2]
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))
from console_formatting import print_session
import project_paths as pp
import rec_list


#%% load data
behaviour = []
behaviour_stem = pp.behaviour_experiment_stem('LC')

for path in rec_list.pathLC:
    recname = Path(path).name
    payload_path = behaviour_stem / f'{recname}.pkl'
    with open(payload_path, 'rb') as f:
        behaviour.append((recname, pickle.load(f)))


#%% main
for recname, row in behaviour:
    print_session(recname)

    output_stem = pp.BEHAVIOUR_FIGURES_STEM / 'trial_profiles' / recname
    output_stem.mkdir(parents=True, exist_ok=True)

    speeds = row['speed_times_aligned']
    tot_trial = len(speeds)

    fig, ax = plt.subplots(figsize=(5,.75))
    speed_max = 0
    holdout_licks = []; holdout_rews = []
    for trial in range(tot_trial):
        if np.mod(trial, 10)==0 and trial!=0:  # plot ten trials at a time
            speed_ylim = max(speed_max+speed_max*.11, 1)
            for lick in holdout_licks:
                ax.vlines(x=lick, ymin=speed_max+speed_max*.01, ymax=speed_max+speed_max*.1, linewidth=.5, color='orchid')
            for rew in holdout_rews:
                ax.vlines(x=rew, ymin=speed_max+speed_max*.01, ymax=speed_max+speed_max*.1, linewidth=.5, color='darkgreen')
            for s in ['top','right']: ax.spines[s].set_visible(False)
            ax.set(xlabel='time (s)',
                   ylabel='speed (cm/s)', ylim=(0, speed_ylim),
                   title='{}_{}'.format(trial-9, trial))
            fig.savefig(output_stem / f'{trial-9}_{trial}.png',
                        dpi=200, bbox_inches='tight')
            fig.savefig(output_stem / f'{trial-9}_{trial}.pdf',
                        bbox_inches='tight')
            plt.close(fig)

            fig, ax = plt.subplots(figsize=(5,.75))
            speed_max = 0
            holdout_licks = []; holdout_rews = []

        curr_speeds = np.asarray(speeds[trial])
        speed_max = max([speed_max, max(curr_speeds[:,1])])
        ax.plot(curr_speeds[:,0]/1000, curr_speeds[:,1], linewidth=1, color='k')

        run_onset = row['run_onsets'][trial]
        if np.isfinite(run_onset):
            curr_start = run_onset/1000
            ax.axvline(x=curr_start, linewidth=1, color='red', alpha=.8)

        if not not(row['lick_times'][trial]):  # if the current trial has any lick
            curr_licks = np.asarray(row['lick_times'][trial])[:,0]/1000
            holdout_licks = np.concatenate((holdout_licks, curr_licks))

        curr_rew = row['reward_times'][trial]
        if not np.isnan(curr_rew):  # if the current trial has a reward
            holdout_rews.append(curr_rew/1000)
