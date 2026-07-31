# -*- coding: utf-8 -*-
'''
Created on Mon Apr  1 17:34:28 2024
Modified on Tue 10 Dec 2024
Modified on 21 Jan 2026

generate profiles for all pyramidal cells in hippocampus recordings,
    segregated into baseline, ctrl and stim trials
A040 and earlier use Hz-scaled all_trains; A041 onwards retain raw
    SciPy-convolved arrays for the first-lick analyses

dependent on hpc_all_extract_raphi.py

@author: Dinghao Luo
'''


#%% to-do
# add depth to the dataframe?
# depth = sio.loadmat('{}\\{}_DataStructure_mazeSection1_TrialType1_Depth.mat'.format(pathname, recname))['depthNeu'][0]
# rel_depth = depth['relDepthNeu'][0][0]

#%% imports
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import os
from pathlib import Path
import sys

from time import time
from datetime import timedelta

import pickle
import pandas as pd
import numpy as np
import scipy.io as sio
from scipy.stats import sem
from tqdm import tqdm

repo_root = Path(__file__).resolve().parents[2]
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))

import hpc_ephys_support as support
from common_functions import mpl_formatting
from console_formatting import print_files_saved, print_session
import project_paths as pp
mpl_formatting()


#%% paths and parameters
HPC_STEM = pp.HPC_EPHYS_STEM
TRAIN_STEM = HPC_STEM / 'all_sessions_raphi'
BEH_STEM = pp.behaviour_experiment_stem('HPCRaphi')
DF_PATH = HPC_STEM / 'hpc_all_profiles_raphi.pkl'


#%% dataframe initialisation/loading
ALT_PRE_1S_WINDOW = (-1.0, 0.0)
ALT_PRE_HALF_S_WINDOW = (-0.5, 0.0)

ALT_PROFILE_COLUMNS = [
    'pre_post_1s',
    'pre_post_stim_1s',
    'pre_post_ctrl_1s',
    'class_1s',
    'class_stim_1s',
    'class_ctrl_1s',
    'pre_post_half_s',
    'pre_post_stim_half_s',
    'pre_post_ctrl_half_s',
    'class_half_s',
    'class_stim_half_s',
    'class_ctrl_half_s',
]

PROFILE_COLUMNS = [
    'recname',        # Axxxr-202xxxxx-0x
    'cell_identity',  # str, 'pyr' or 'int'
    'depth',          # depth relative to layer centre
    'firing_rate',    # in Hz
    'place_cell',
    'pre_post',       # post/pre ([.5:1.5]/[-1.5:-.5])
    'pre_post_stim',  # in stim trials
    'pre_post_ctrl',  # in stim-ctrl trials
    'class',          # run-onset activated/inhibited/unresponsive
    'class_stim',
    'class_ctrl',
    'var',            # trial-by-trial variability in firing
    'var_stim',
    'var_ctrl',
    'SI',             # spatial information
    'SI_stim',
    'SI_ctrl',
    'TI',             # temporal information
    'TI_stim',
    'TI_ctrl',
    'prof_mean',      # mean firing profile
    'prof_sem',
    'prof_stim_mean',
    'prof_stim_sem',
    'prof_ctrl_mean',
    'prof_ctrl_sem',
    *ALT_PROFILE_COLUMNS,
]

def save_profiles(df):
    DF_PATH.parent.mkdir(parents=True, exist_ok=True)
    df.to_pickle(DF_PATH)
    print_files_saved([
        ('dataframe', DF_PATH),
    ])


#%% load paths to recordings
import rec_list
paths = rec_list.pathHPC_Raphi


#%% parameters
# behaviour
track_length  = 200  # in cm
bin_size      = 0.1  # in cm

# ephys
SAMP_FREQ   = 1250  # in Hz
MAX_TIME    = 10  # collect (for each trial) a maximum of 10 s of spiking-profile
MAX_SAMPLES = SAMP_FREQ * MAX_TIME

# pre_post ratio thresholds
run_onset_activated_thres = 2/3
run_onset_inhibited_thres = 3/2


#%% main
def process_recording(path, show_cell_progress=False, verbose=False):
    recname = Path(path).name
    t0 = time()

    beh_path = BEH_STEM / f'{recname}.pkl'
    with open(beh_path, 'rb') as f:
        beh = pickle.load(f)

    speeds = support.load_speeds(beh)

    distance_bins = np.arange(0, track_length + bin_size, bin_size)
    occupancy = [
        support.calculate_occupancy(s, dt=.02, distance_bins=distance_bins)
        for s in speeds
    ]

    train_path = TRAIN_STEM / recname / f'{recname}_all_trains.npy'
    train_file = np.load(train_path, allow_pickle=True).item()
    clu_list = list(train_file.keys())
    trains = list(train_file.values())

    session_path = Path(path)
    train_dist_path = None
    train_run = 0
    for maze_sess_candidate in range(10):
        candidate_path = session_path / f'{recname}_DataStructure_mazeSection1_TrialType1_convSpikesDistAligned_msess{maze_sess_candidate}_Run0.mat'
        if candidate_path.exists():
            maze_sess = maze_sess_candidate
            train_dist_path = candidate_path
            break

    if train_dist_path is None:
        train_run = 1
        for maze_sess_candidate in range(10):
            candidate_path = session_path / f'{recname}_DataStructure_mazeSection1_TrialType1_convSpikesDistAligned_msess{maze_sess_candidate}_Run1.mat'
            if candidate_path.exists():
                maze_sess = maze_sess_candidate
                train_dist_path = candidate_path
                if verbose:
                    print(f'run1 distance-spike file used for {recname} due to absence of run0 distance-spike file')
                break

    if train_dist_path is None:
        raise FileNotFoundError(f'No Run0 or Run1 distance-spike file found for {recname}')

    place_cell_info_path = session_path / f'{recname}_DataStructure_mazeSection1_TrialType1_FieldSpCorrAligned_Run{maze_sess}_Run{train_run}.mat'
    if not place_cell_info_path.exists() and train_run == 0:
        run1_place_cell_info_path = session_path / f'{recname}_DataStructure_mazeSection1_TrialType1_FieldSpCorrAligned_Run{maze_sess}_Run1.mat'
        if run1_place_cell_info_path.exists():
            place_cell_info_path = run1_place_cell_info_path
            if verbose:
                print(f'run1 place-cell file used for {recname} due to absence of run0 place-cell file; distance spikes remain run0')

    if not place_cell_info_path.exists():
        raise FileNotFoundError(f'No matching place-cell file found for {recname}: {place_cell_info_path}')

    trains_dist = support.load_dist_spike_array(train_dist_path)

    cell_info_path = Path(path) / f'{recname}_DataStructure_mazeSection1_TrialType1_Info.mat'
    cell_identities, spike_rates = support.get_cell_info(cell_info_path)

    # MATLAB cell indices are 1-indexed; profile arrays are 0-indexed
    place_cell_idx = sio.loadmat(str(place_cell_info_path))[
        'fieldSpCorrSessNonStimGood'
    ][0][0]['indNeuron'][0] - 1

    depth_info_path = Path(path) / f'{recname}_DataStructure_mazeSection1_TrialType1_Depth.mat'
    depths = sio.loadmat(str(depth_info_path))['depthNeu'][0]['relDepthNeu'][0][0]

    beh_MATLAB_path = Path(path) / f'{recname}_DataStructure_mazeSection1_TrialType1_behPar_msess{maze_sess}.mat'
    (
        baseline_idx,
        stim_idx,
        ctrl_idx
    ) = support.get_trialtype_idx_MATLAB(beh_MATLAB_path)

    baseline_idx = baseline_idx[:-1]
    stim_idx = [t - 1 for t in stim_idx]
    ctrl_idx = [t - 1 for t in ctrl_idx]
    session_rows = {}

    for clu in tqdm(
        range(len(cell_identities)),
        desc=f'{recname} profiles',
        disable=not show_cell_progress,
    ):
        cell_identity = cell_identities[clu]
        if cell_identity == 'putative_pyr':
            if 0.15 < spike_rates[clu] < 7:  # modified 31 Mar 2025
                cell_identity = 'pyr'
            else:
                cell_identity = 'other'

        depth = depths[clu]

        baseline_matrix = support.get_trial_matrix(
            trains, baseline_idx, MAX_SAMPLES, clu)
        ctrl_matrix = support.get_trial_matrix(
            trains, ctrl_idx, MAX_SAMPLES, clu)
        stim_matrix = support.get_trial_matrix(
            trains, stim_idx, MAX_SAMPLES, clu)

        # some recordings have no ctrl or stim trials, so empty-slice warnings can be expected
        baseline_mean = np.nanmean(baseline_matrix, axis=0)
        ctrl_mean = np.nanmean(ctrl_matrix, axis=0)
        stim_mean = np.nanmean(stim_matrix, axis=0)

        baseline_sem = sem(baseline_matrix, axis=0)
        ctrl_sem = sem(ctrl_matrix, axis=0)
        stim_sem = sem(stim_matrix, axis=0)

        (
            baseline_run_onset_ratio,
            baseline_run_onset_ratiotype
        ) = support.classify_run_onset_activation_ratio(
            baseline_mean,
            run_onset_activated_thres,
            run_onset_inhibited_thres
        )
        if ctrl_idx:
            (
                ctrl_run_onset_ratio,
                ctrl_run_onset_ratiotype
            ) = support.classify_run_onset_activation_ratio(
                ctrl_mean,
                run_onset_activated_thres,
                run_onset_inhibited_thres
            )
            (
                stim_run_onset_ratio,
                stim_run_onset_ratiotype
            ) = support.classify_run_onset_activation_ratio(
                stim_mean,
                run_onset_activated_thres,
                run_onset_inhibited_thres
            )
        else:
            ctrl_run_onset_ratio = np.nan
            ctrl_run_onset_ratiotype = np.nan
            stim_run_onset_ratio = np.nan
            stim_run_onset_ratiotype = np.nan

        baseline_ratio_1s, baseline_class_1s = support.classify_run_onset_activation_ratio(
            baseline_mean,
            run_onset_activated_thres,
            run_onset_inhibited_thres,
            pre_window=ALT_PRE_1S_WINDOW
        )
        baseline_ratio_half_s, baseline_class_half_s = support.classify_run_onset_activation_ratio(
            baseline_mean,
            run_onset_activated_thres,
            run_onset_inhibited_thres,
            pre_window=ALT_PRE_HALF_S_WINDOW
        )
        if ctrl_idx:
            ctrl_ratio_1s, ctrl_class_1s = support.classify_run_onset_activation_ratio(
                ctrl_mean,
                run_onset_activated_thres,
                run_onset_inhibited_thres,
                pre_window=ALT_PRE_1S_WINDOW
            )
            stim_ratio_1s, stim_class_1s = support.classify_run_onset_activation_ratio(
                stim_mean,
                run_onset_activated_thres,
                run_onset_inhibited_thres,
                pre_window=ALT_PRE_1S_WINDOW
            )
            ctrl_ratio_half_s, ctrl_class_half_s = support.classify_run_onset_activation_ratio(
                ctrl_mean,
                run_onset_activated_thres,
                run_onset_inhibited_thres,
                pre_window=ALT_PRE_HALF_S_WINDOW
            )
            stim_ratio_half_s, stim_class_half_s = support.classify_run_onset_activation_ratio(
                stim_mean,
                run_onset_activated_thres,
                run_onset_inhibited_thres,
                pre_window=ALT_PRE_HALF_S_WINDOW
            )
        else:
            ctrl_ratio_1s = np.nan
            ctrl_class_1s = np.nan
            stim_ratio_1s = np.nan
            stim_class_1s = np.nan
            ctrl_ratio_half_s = np.nan
            ctrl_class_half_s = np.nan
            stim_ratio_half_s = np.nan
            stim_class_half_s = np.nan

        baseline_var = support.compute_trial_by_trial_variability(baseline_matrix)
        if ctrl_idx:
            ctrl_var = support.compute_trial_by_trial_variability(ctrl_matrix)
            stim_var = support.compute_trial_by_trial_variability(stim_matrix)
        else:
            ctrl_var = stim_var = np.nan

        baseline_SI = [support.compute_spatial_information(
            trains_dist[clu][trial], occupancy[trial])
            for trial in baseline_idx]
        ctrl_SI = [support.compute_spatial_information(
            trains_dist[clu][trial], occupancy[trial])
            for trial in ctrl_idx]
        stim_SI = [support.compute_spatial_information(
            trains_dist[clu][trial], occupancy[trial])
            for trial in stim_idx]

        baseline_TI = [support.compute_temporal_information(
            trains[clu][trial][SAMP_FREQ*3:],
            bin_size_steps=1
            ) for trial in baseline_idx
            if trains[clu][trial] is not None]
        ctrl_TI = [support.compute_temporal_information(
            trains[clu][trial][SAMP_FREQ*3:],
            bin_size_steps=1) for trial in ctrl_idx
            if trains[clu][trial] is not None]
        stim_TI = [support.compute_temporal_information(
            trains[clu][trial][SAMP_FREQ*3:],
            bin_size_steps=1) for trial in stim_idx
            if trains[clu][trial] is not None]

        cluname = clu_list[clu]
        session_rows[cluname] = [
            recname,
            cell_identity,
            depth,
            spike_rates[clu],
            clu in place_cell_idx,
            baseline_run_onset_ratio,
            stim_run_onset_ratio,
            ctrl_run_onset_ratio,
            baseline_run_onset_ratiotype,
            stim_run_onset_ratiotype,
            ctrl_run_onset_ratiotype,
            baseline_var,
            stim_var,
            ctrl_var,
            baseline_SI,
            stim_SI,
            ctrl_SI,
            baseline_TI,
            stim_TI,
            ctrl_TI,
            baseline_mean,
            baseline_sem,
            stim_mean,
            stim_sem,
            ctrl_mean,
            ctrl_sem,
            baseline_ratio_1s,
            stim_ratio_1s,
            ctrl_ratio_1s,
            baseline_class_1s,
            stim_class_1s,
            ctrl_class_1s,
            baseline_ratio_half_s,
            stim_ratio_half_s,
            ctrl_ratio_half_s,
            baseline_class_half_s,
            stim_class_half_s,
            ctrl_class_half_s,
        ]

    session_df = pd.DataFrame.from_dict(
        session_rows,
        orient='index',
        columns=PROFILE_COLUMNS,
    ).astype('object')
    return recname, session_df, time() - t0

def main(argv=None):
    parser = argparse.ArgumentParser(description='generate Raphi HPC profiles')
    parser.add_argument(
        '--n-workers',
        '--workers',
        dest='n_workers',
        type=int,
        default=max(1, min(6, (os.cpu_count() or 1) - 1)),
        help='number of recordings to process in parallel; use 1 for sequential processing',
    )
    parser.add_argument(
        '--verbose',
        action='store_true',
        help='print per-recording details',
    )
    args = parser.parse_args(argv)

    print('not using GPU acceleration; profiles are CPU/session parallel\n')

    df = pd.DataFrame({column: pd.Series(dtype='object') for column in PROFILE_COLUMNS})
    n_workers = min(max(1, args.n_workers), len(paths))
    if n_workers == 1:
        path_iter = paths if args.verbose else tqdm(paths, desc='recordings')
        for path in path_iter:
            if not args.verbose:
                path_iter.set_postfix_str(Path(path).name)
            else:
                print_session(Path(path).name)
            recname, session_df, elapsed = process_recording(
                path,
                show_cell_progress=args.verbose,
                verbose=args.verbose,
            )
            df = pd.concat([df, session_df], axis=0)
            save_profiles(df)
            print(f'{recname} done in {timedelta(seconds=int(elapsed))}\n')
    else:
        print(f'processing {len(paths)} recordings with {n_workers} workers')
        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            future_to_recname = {
                executor.submit(process_recording, path, False, args.verbose): Path(path).name
                for path in paths
            }
            for future in tqdm(
                    as_completed(future_to_recname),
                    total=len(future_to_recname),
                    desc='recordings',
                ):
                recname, session_df, elapsed = future.result()
                df = pd.concat([df, session_df], axis=0)
                save_profiles(df)
                print(f'{recname} done in {timedelta(seconds=int(elapsed))}\n')

    save_profiles(df)

if __name__ == '__main__':
    main()
