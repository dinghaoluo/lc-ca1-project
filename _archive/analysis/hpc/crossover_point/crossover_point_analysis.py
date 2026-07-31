# -*- coding: utf-8 -*-
"""
Created on Thu Mar 27 17:51:40 2025

analyse the crossover point of PyrUp and PyrDown cells at single-trial level

@author: Dinghao Luo
"""

#%% imports
import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio
import pandas as pd
import sys
from pathlib import Path

repo_root = Path(__file__).resolve().parent
while not (repo_root / 'utils').exists() and repo_root != repo_root.parent:
    repo_root = repo_root.parent
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))

import rec_list
paths = rec_list.pathHPCLCopt + rec_list.pathHPCLCtermopt
# paths = rec_list.pathHPCLCopt

from plotting_functions import plot_violin_with_scatter

import hpc_ephys_support as support
import project_paths as pp


#%% parameters
SAMP_FREQ = 1250  # Hz
TIME = np.arange(-SAMP_FREQ*3, SAMP_FREQ*7)/SAMP_FREQ  # 8 seconds
PYRUP_CLASS = 'run-onset ON'
PYRDOWN_CLASS = 'run-onset OFF'


#%% load data
HPC_EPHYS_STEM = pp.HPC_EPHYS_STEM
HPC_ALL_SESSIONS_STEM = HPC_EPHYS_STEM / 'all_sessions'
BEHAVIOUR_STEM = pp.BEHAVIOUR_STEM


#%% main loop
def main():
    print('loading dataframes...')

    cell_profiles = pd.read_pickle(
        HPC_EPHYS_STEM / 'HPC_all_profiles.pkl'
        )
    df_pyr = cell_profiles[cell_profiles['cell_identity']=='pyr']  # pyramidal only

    beh_df = pd.concat((
        pd.read_pickle(BEHAVIOUR_STEM / 'all_HPCLC_sessions.pkl'),
        pd.read_pickle(BEHAVIOUR_STEM / 'all_HPCLCterm_sessions.pkl')
        ))

    for path in paths:
        recname = path[-17:]
        print(recname)

        # get cells
        curr_df = df_pyr[df_pyr['recname']==recname]
        curr_pyrup_idx = [int(s.split(' ')[-3][3:])-2 for s in curr_df.index
                          if curr_df.loc[s]['class'] == PYRUP_CLASS]
        curr_pyrdown_idx = [int(s.split(' ')[-3][3:])-2 for s in curr_df.index
                            if curr_df.loc[s]['class'] == PYRDOWN_CLASS]

        # filtering
        if len(curr_pyrup_idx) < 5 or len(curr_pyrdown_idx) < 5:
            continue

        # get lick times
        curr_beh_df = beh_df.loc[recname]  # subselect in read-only
        run_onsets = curr_beh_df['run_onsets'][1:]
        licks = [
            [(l-run_onset) for l in trial]
            if len(trial)!=0 else np.nan
            for trial, run_onset in zip(
                    curr_beh_df['lick_times'][1:],
                    run_onsets
                    )
            ]
        first_licks = np.asarray(
            [next((l for l in trial if l > 1), np.nan)  # >1 to eliminate carry-over licks
            if isinstance(trial, list) else np.nan
            for trial in licks]
            )

        stim_trials = np.where(
            np.asarray([
                trial[15] for trial
                in curr_beh_df['trial_statements']
                ])!='0'
            )[0]
        ctrl_trials = stim_trials+2

        # get current session spike trains
        clu_list, trains = support.load_train(
            str(HPC_ALL_SESSIONS_STEM / recname / f'{recname}_all_trains.npy')
            )

        # trials
        trials = np.arange(len(trains[0]))

        # main loop
        pyrup_mean_aligned = np.zeros((len(trials), 1250*6))
        pyrdown_mean_aligned = np.zeros((len(trials), 1250*6))
        for trial in trials:
            try:
                first_lick_bin = int(first_licks[trial]/1000*1250-3750)
                curr_pyrup_profiles = [trains[clu][trial] for clu in curr_pyrup_idx]
                curr_pyrup_mean = np.mean(curr_pyrup_profiles, axis=0)

                curr_pyrdown_profiles = [trains[clu][trial] for clu in curr_pyrdown_idx]
                curr_pyrdown_mean = np.mean(curr_pyrdown_profiles, axis=0)

                curr_pyrup_mean_aligned = curr_pyrup_mean[first_lick_bin-4*1250:first_lick_bin+2*1250]
                curr_pyrdown_mean_aligned = curr_pyrdown_mean[first_lick_bin-4*1250:first_lick_bin+2*1250]
                pyrup_mean_aligned[trial, :] = curr_pyrup_mean_aligned
                pyrdown_mean_aligned[trial, :] = curr_pyrdown_mean_aligned
            except ValueError:
                pyrup_mean_aligned[trial, :] = np.nan
                pyrdown_mean_aligned[trial, :] = np.nan

        pyrup_mean_aligned_mean = np.nanmean(pyrup_mean_aligned, axis=0)
        pyrdown_mean_aligned_mean = np.nanmean(pyrdown_mean_aligned, axis=0)

        fig, ax = plt.subplots(figsize=(3,2))
        ax.plot(np.arange(-1250*4, 1250*2)/1250,
                pyrup_mean_aligned_mean)
        ax.plot(np.arange(-1250*4, 1250*2)/1250,
                pyrdown_mean_aligned_mean)
        ax.set_title(recname)


if __name__ == '__main__':
    main()
