# -*- coding: utf-8 -*-
'''
Created on Mon Apr  1 17:34:28 2024
Modified on Tue 10 Dec 2024

generate profiles for all pyramidal cells in hippocampus recordings,
    segregated into baseline, ctrl and stim trials

dependent on hpc_all_extract.py

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

from tqdm import tqdm
from time import time
from datetime import timedelta
import sys
import pickle
import pandas as pd
import numpy as np
import scipy.io as sio
from scipy.stats import sem

repo_root = Path(__file__).resolve().parents[2]
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))
from common_functions import mpl_formatting
from console_formatting import print_files_saved, print_session
mpl_formatting()

import hpc_ephys_support as support
import project_paths as pp


#%% load paths to recordings
import rec_list


#%% parameters
OUTPUT_PATH = pp.HPC_EPHYS_STEM / 'hpc_all_profiles.pkl'
TRACK_LENGTH = 200  # in cm
BIN_SIZE = 0.1  # in cm
RUN_ONSET_BIN = 3750  # in bins
SAMP_FREQ = 1250  # in Hz
MAX_TIME = 10  # collect (for each trial) a maximum of 10 s of spiking-profile
MAX_SAMPLES = SAMP_FREQ * MAX_TIME
RUN_ONSET_ACTIVATED_THRES = 2 / 3
RUN_ONSET_INHIBITED_THRES = 3 / 2
DISTANCE_BINS = np.arange(0, TRACK_LENGTH + BIN_SIZE, BIN_SIZE)
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
    'rectype',
    'recname',
    'cell_identity',
    'depth',
    'spike_rate',
    'place_cell',
    'pre_post',
    'pre_post_stim',
    'pre_post_ctrl',
    'class',
    'class_stim',
    'class_ctrl',
    'pre_post_MATLAB',
    'pre_post_stim_MATLAB',
    'pre_post_ctrl_MATLAB',
    'class_MATLAB',
    'class_stim_MATLAB',
    'class_ctrl_MATLAB',
    'var',
    'var_stim',
    'var_ctrl',
    'SI',
    'SI_stim',
    'SI_ctrl',
    'TI',
    'TI_stim',
    'TI_ctrl',
    'var_MATLAB',
    'var_stim_MATLAB',
    'var_ctrl_MATLAB',
    'SI_MATLAB',
    'SI_stim_MATLAB',
    'SI_ctrl_MATLAB',
    'TI_MATLAB',
    'TI_stim_MATLAB',
    'TI_ctrl_MATLAB',
    'prof_mean',
    'prof_sem',
    'prof_stim_mean',
    'prof_stim_sem',
    'prof_ctrl_mean',
    'prof_ctrl_sem',
    'prof_mean_MATLAB',
    'prof_sem_MATLAB',
    'prof_stim_mean_MATLAB',
    'prof_stim_sem_MATLAB',
    'prof_ctrl_mean_MATLAB',
    'prof_ctrl_sem_MATLAB',
    'prof_good_mean',
    'prof_good_sem',
    'prof_bad_mean',
    'prof_bad_sem',
    'prof_good_mean_matlab',
    'prof_good_sem_matlab',
    'prof_bad_mean_matlab',
    'prof_bad_sem_matlab',
    'MI',
    'MI_early',
    'MI_late',
    'MI_shuf',
    'MI_shuf_early',
    'MI_shuf_late',
    *ALT_PROFILE_COLUMNS,
]


#%% profile extraction and saving
def save_profiles(df):
    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    df.to_pickle(OUTPUT_PATH)
    print_files_saved([
        ('dataframe', OUTPUT_PATH),
    ])

def compute_trial_profiles(trains, trial_indices, clu):
    matrix = support.get_trial_matrix(trains, trial_indices, MAX_SAMPLES, clu)
    mean_profile = np.nanmean(matrix, axis=0)
    sem_profile = sem(matrix, axis=0)
    return matrix, mean_profile, sem_profile

def load_session_inputs(pathname, recname, prefix):
    beh_path = pp.behaviour_experiment_stem(prefix) / f'{recname}.pkl'
    with open(beh_path, 'rb') as f:
        beh_df = pickle.load(f)

    speeds = support.load_speeds(beh_df)
    # trial 1 is empty and not included in the spike train
    bad_trial_map = beh_df['bad_trials']
    good_idx = [
        trial - 1 for trial, quality in enumerate(bad_trial_map)
        if not quality and trial > 0
    ]
    bad_idx = [
        trial - 1 for trial, quality in enumerate(bad_trial_map)
        if quality and trial > 0
    ]
    good_idx_matlab, bad_idx_matlab = support.get_good_bad_idx_MATLAB(pathname)

    occupancy = [
        support.calculate_occupancy(speed_row, dt=.02, distance_bins=DISTANCE_BINS)
        for speed_row in speeds
    ]

    train_path = pp.HPC_EPHYS_STEM / 'all_sessions' / recname / f'{recname}_all_trains.npy'
    train_file = np.load(str(train_path), allow_pickle=True).item()
    clu_list = list(train_file.keys())
    trains = list(train_file.values())
    trains_dist = support.load_dist_spike_array(
        r'{}\{}_DataStructure_mazeSection1_TrialType1_convSpikesDistAligned_msess1_Run0.mat'.format(
            pathname,
            recname
        )
    )
    cell_identities, spike_rates = support.get_cell_info(
        r'{}\{}_DataStructure_mazeSection1_TrialType1_Info.mat'.format(pathname, recname)
    )
    place_cell_info_path = r'{}\{}_DataStructure_mazeSection1_TrialType1_FieldSpCorrAligned_Run1_Run0.mat'.format(
        pathname,
        recname
    )
    # MATLAB cell indices are 1-indexed; profile arrays are 0-indexed
    place_cell_idx = sio.loadmat(place_cell_info_path)[
        'fieldSpCorrSessNonStimGood'
    ][0][0]['indNeuron'][0] - 1
    depth_info_path = r'{}\{}_DataStructure_mazeSection1_TrialType1_Depth.mat'.format(
        pathname, recname
    )
    depths = sio.loadmat(depth_info_path)['depthNeu'][0]['relDepthNeu'][0][0]

    (
        baseline_idx_matlab,
        stim_idx_matlab,
        ctrl_idx_matlab
    ) = support.get_trialtype_idx_MATLAB(
        r'{}\{}_DataStructure_mazeSection1_TrialType1_behPar_msess1.mat'.format(
            pathname,
            recname
        )
    )

    stim_conds = [trial[15] for trial in beh_df['trial_statements']][1:]
    stim_idx = [trial for trial, cond in enumerate(stim_conds) if cond != '0']
    if stim_idx:
        baseline_idx = list(range(stim_idx[0]))
        ctrl_idx = [idx + 2 for idx in stim_idx if idx + 2 < len(stim_conds)]
    else:
        baseline_idx = list(range(len(stim_conds)))
        ctrl_idx = []

    return {
        'prefix': prefix,
        'recname': recname,
        'beh_df': beh_df,
        'good_idx': good_idx,
        'bad_idx': bad_idx,
        'good_idx_matlab': good_idx_matlab,
        'bad_idx_matlab': bad_idx_matlab,
        'occupancy': occupancy,
        'clu_list': clu_list,
        'trains': trains,
        'trains_dist': trains_dist,
        'cell_identities': cell_identities,
        'spike_rates': spike_rates,
        'place_cell_idx': place_cell_idx,
        'depths': depths,
        'baseline_idx': baseline_idx,
        'stim_idx': stim_idx,
        'ctrl_idx': ctrl_idx,
        'baseline_idx_matlab': baseline_idx_matlab,
        'stim_idx_matlab': stim_idx_matlab,
        'ctrl_idx_matlab': ctrl_idx_matlab
    }

def build_cell_profile_row(clu, session_inputs):
    trains = session_inputs['trains']
    trains_cell = trains[clu]
    trains_dist_cell = session_inputs['trains_dist'][clu]

    cell_identity = session_inputs['cell_identities'][clu]
    if cell_identity == 'putative_pyr':
        if 0.15 < session_inputs['spike_rates'][clu] < 7:  # modified 31 Mar 2025
            cell_identity = 'pyr'
        else:
            cell_identity = 'other'
    depth = session_inputs['depths'][clu]

    baseline_matrix, baseline_mean, baseline_sem = compute_trial_profiles(
        trains,
        session_inputs['baseline_idx'],
        clu
    )
    ctrl_matrix, ctrl_mean, ctrl_sem = compute_trial_profiles(
        trains,
        session_inputs['ctrl_idx'],
        clu
    )
    stim_matrix, stim_mean, stim_sem = compute_trial_profiles(
        trains,
        session_inputs['stim_idx'],
        clu
    )

    baseline_matrix_matlab, baseline_mean_matlab, baseline_sem_matlab = compute_trial_profiles(
        trains,
        session_inputs['baseline_idx_matlab'],
        clu
    )
    ctrl_matrix_matlab, ctrl_mean_matlab, ctrl_sem_matlab = compute_trial_profiles(
        trains,
        session_inputs['ctrl_idx_matlab'],
        clu
    )
    stim_matrix_matlab, stim_mean_matlab, stim_sem_matlab = compute_trial_profiles(
        trains,
        session_inputs['stim_idx_matlab'],
        clu
    )

    baseline_run_onset_ratio, baseline_run_onset_ratiotype = support.classify_run_onset_activation_ratio(
        baseline_mean,
        RUN_ONSET_ACTIVATED_THRES,
        RUN_ONSET_INHIBITED_THRES
    )
    ctrl_run_onset_ratio, ctrl_run_onset_ratiotype = support.classify_run_onset_activation_ratio(
        ctrl_mean,
        RUN_ONSET_ACTIVATED_THRES,
        RUN_ONSET_INHIBITED_THRES
    )
    stim_run_onset_ratio, stim_run_onset_ratiotype = support.classify_run_onset_activation_ratio(
        stim_mean,
        RUN_ONSET_ACTIVATED_THRES,
        RUN_ONSET_INHIBITED_THRES
    )

    baseline_ratio_1s, baseline_class_1s = support.classify_run_onset_activation_ratio(
        baseline_mean,
        RUN_ONSET_ACTIVATED_THRES,
        RUN_ONSET_INHIBITED_THRES,
        pre_window=ALT_PRE_1S_WINDOW
    )
    ctrl_ratio_1s, ctrl_class_1s = support.classify_run_onset_activation_ratio(
        ctrl_mean,
        RUN_ONSET_ACTIVATED_THRES,
        RUN_ONSET_INHIBITED_THRES,
        pre_window=ALT_PRE_1S_WINDOW
    )
    stim_ratio_1s, stim_class_1s = support.classify_run_onset_activation_ratio(
        stim_mean,
        RUN_ONSET_ACTIVATED_THRES,
        RUN_ONSET_INHIBITED_THRES,
        pre_window=ALT_PRE_1S_WINDOW
    )

    baseline_ratio_half_s, baseline_class_half_s = support.classify_run_onset_activation_ratio(
        baseline_mean,
        RUN_ONSET_ACTIVATED_THRES,
        RUN_ONSET_INHIBITED_THRES,
        pre_window=ALT_PRE_HALF_S_WINDOW
    )
    ctrl_ratio_half_s, ctrl_class_half_s = support.classify_run_onset_activation_ratio(
        ctrl_mean,
        RUN_ONSET_ACTIVATED_THRES,
        RUN_ONSET_INHIBITED_THRES,
        pre_window=ALT_PRE_HALF_S_WINDOW
    )
    stim_ratio_half_s, stim_class_half_s = support.classify_run_onset_activation_ratio(
        stim_mean,
        RUN_ONSET_ACTIVATED_THRES,
        RUN_ONSET_INHIBITED_THRES,
        pre_window=ALT_PRE_HALF_S_WINDOW
    )

    baseline_run_onset_ratio_matlab, baseline_run_onset_ratiotype_matlab = support.classify_run_onset_activation_ratio(
        baseline_mean_matlab,
        RUN_ONSET_ACTIVATED_THRES,
        RUN_ONSET_INHIBITED_THRES
    )
    ctrl_run_onset_ratio_matlab, ctrl_run_onset_ratiotype_matlab = support.classify_run_onset_activation_ratio(
        ctrl_mean_matlab,
        RUN_ONSET_ACTIVATED_THRES,
        RUN_ONSET_INHIBITED_THRES
    )
    stim_run_onset_ratio_matlab, stim_run_onset_ratiotype_matlab = support.classify_run_onset_activation_ratio(
        stim_mean_matlab,
        RUN_ONSET_ACTIVATED_THRES,
        RUN_ONSET_INHIBITED_THRES
    )

    MI, MI_early, MI_late = support.compute_modulation_index(
        ctrl_mean,
        stim_mean,
        run_onset_bin=RUN_ONSET_BIN
    )
    MI_shuf, MI_early_shuf, MI_late_shuf = support.compute_modulation_index_shuf(
        ctrl_matrix,
        stim_matrix,
        run_onset_bin=RUN_ONSET_BIN
    )

    baseline_var = support.compute_trial_by_trial_variability(baseline_matrix)
    ctrl_var = support.compute_trial_by_trial_variability(ctrl_matrix)
    stim_var = support.compute_trial_by_trial_variability(stim_matrix)

    baseline_var_matlab = support.compute_trial_by_trial_variability(baseline_matrix_matlab)
    ctrl_var_matlab = support.compute_trial_by_trial_variability(ctrl_matrix_matlab)
    stim_var_matlab = support.compute_trial_by_trial_variability(stim_matrix_matlab)

    baseline_SI = [
        support.compute_spatial_information(trains_dist_cell[trial], session_inputs['occupancy'][trial])
        for trial in session_inputs['baseline_idx']
    ]
    ctrl_SI = [
        support.compute_spatial_information(trains_dist_cell[trial], session_inputs['occupancy'][trial])
        for trial in session_inputs['ctrl_idx']
    ]
    stim_SI = [
        support.compute_spatial_information(trains_dist_cell[trial], session_inputs['occupancy'][trial])
        for trial in session_inputs['stim_idx']
    ]

    baseline_SI_matlab = [
        support.compute_spatial_information(trains_dist_cell[trial], session_inputs['occupancy'][trial])
        for trial in session_inputs['baseline_idx_matlab']
    ]
    ctrl_SI_matlab = [
        support.compute_spatial_information(trains_dist_cell[trial], session_inputs['occupancy'][trial])
        for trial in session_inputs['ctrl_idx_matlab']
    ]
    stim_SI_matlab = [
        support.compute_spatial_information(trains_dist_cell[trial], session_inputs['occupancy'][trial])
        for trial in session_inputs['stim_idx_matlab']
    ]

    baseline_TI = [
        support.compute_temporal_information(trains_cell[trial][SAMP_FREQ * 3:], bin_size_steps=1)
        for trial in session_inputs['baseline_idx']
        if trains_cell[trial] is not None
    ]
    ctrl_TI = [
        support.compute_temporal_information(trains_cell[trial][SAMP_FREQ * 3:], bin_size_steps=1)
        for trial in session_inputs['ctrl_idx']
        if trains_cell[trial] is not None
    ]
    stim_TI = [
        support.compute_temporal_information(trains_cell[trial][SAMP_FREQ * 3:], bin_size_steps=1)
        for trial in session_inputs['stim_idx']
        if trains_cell[trial] is not None
    ]

    baseline_TI_matlab = [
        support.compute_temporal_information(trains_cell[trial][SAMP_FREQ * 3:], bin_size_steps=1)
        for trial in session_inputs['baseline_idx_matlab']
        if trains_cell[trial] is not None
    ]
    ctrl_TI_matlab = [
        support.compute_temporal_information(trains_cell[trial][SAMP_FREQ * 3:], bin_size_steps=1)
        for trial in session_inputs['ctrl_idx_matlab']
        if trains_cell[trial] is not None
    ]
    stim_TI_matlab = [
        support.compute_temporal_information(trains_cell[trial][SAMP_FREQ * 3:], bin_size_steps=1)
        for trial in session_inputs['stim_idx_matlab']
        if trains_cell[trial] is not None
    ]

    good_matrix = support.get_trial_matrix(trains, session_inputs['good_idx'], MAX_SAMPLES, clu)
    good_mean = np.nanmean(good_matrix, axis=0) if session_inputs['good_idx'] else np.array([])
    good_sem = sem(good_matrix, axis=0) if session_inputs['good_idx'] else np.array([])
    bad_matrix = support.get_trial_matrix(trains, session_inputs['bad_idx'], MAX_SAMPLES, clu)
    bad_mean = np.nanmean(bad_matrix, axis=0) if session_inputs['bad_idx'] else np.array([])
    bad_sem = sem(bad_matrix, axis=0) if session_inputs['bad_idx'] else np.array([])
    good_matrix_matlab = support.get_trial_matrix(trains, session_inputs['good_idx_matlab'], MAX_SAMPLES, clu)
    good_mean_matlab = np.nanmean(good_matrix_matlab, axis=0) if session_inputs['good_idx_matlab'] else np.array([])
    good_sem_matlab = sem(good_matrix_matlab, axis=0) if session_inputs['good_idx_matlab'] else np.array([])
    bad_matrix_matlab = support.get_trial_matrix(trains, session_inputs['bad_idx_matlab'], MAX_SAMPLES, clu)
    bad_mean_matlab = np.nanmean(bad_matrix_matlab, axis=0) if session_inputs['bad_idx_matlab'] else np.array([])
    bad_sem_matlab = sem(bad_matrix_matlab, axis=0) if session_inputs['bad_idx_matlab'] else np.array([])

    row = {
        'rectype': session_inputs['prefix'],
        'recname': session_inputs['recname'],
        'cell_identity': cell_identity,
        'depth': depth,
        'spike_rate': session_inputs['spike_rates'][clu],
        'place_cell': clu in session_inputs['place_cell_idx'],
        'pre_post': baseline_run_onset_ratio,
        'pre_post_stim': stim_run_onset_ratio,
        'pre_post_ctrl': ctrl_run_onset_ratio,
        'class': baseline_run_onset_ratiotype,
        'class_stim': stim_run_onset_ratiotype,
        'class_ctrl': ctrl_run_onset_ratiotype,
        'pre_post_1s': baseline_ratio_1s,
        'pre_post_stim_1s': stim_ratio_1s,
        'pre_post_ctrl_1s': ctrl_ratio_1s,
        'class_1s': baseline_class_1s,
        'class_stim_1s': stim_class_1s,
        'class_ctrl_1s': ctrl_class_1s,
        'pre_post_half_s': baseline_ratio_half_s,
        'pre_post_stim_half_s': stim_ratio_half_s,
        'pre_post_ctrl_half_s': ctrl_ratio_half_s,
        'class_half_s': baseline_class_half_s,
        'class_stim_half_s': stim_class_half_s,
        'class_ctrl_half_s': ctrl_class_half_s,
        'pre_post_MATLAB': baseline_run_onset_ratio_matlab,
        'pre_post_stim_MATLAB': stim_run_onset_ratio_matlab,
        'pre_post_ctrl_MATLAB': ctrl_run_onset_ratio_matlab,
        'class_MATLAB': baseline_run_onset_ratiotype_matlab,
        'class_stim_MATLAB': stim_run_onset_ratiotype_matlab,
        'class_ctrl_MATLAB': ctrl_run_onset_ratiotype_matlab,
        'var': baseline_var,
        'var_stim': stim_var,
        'var_ctrl': ctrl_var,
        'SI': baseline_SI,
        'SI_stim': stim_SI,
        'SI_ctrl': ctrl_SI,
        'TI': baseline_TI,
        'TI_stim': stim_TI,
        'TI_ctrl': ctrl_TI,
        'var_MATLAB': baseline_var_matlab,
        'var_stim_MATLAB': stim_var_matlab,
        'var_ctrl_MATLAB': ctrl_var_matlab,
        'SI_MATLAB': baseline_SI_matlab,
        'SI_stim_MATLAB': stim_SI_matlab,
        'SI_ctrl_MATLAB': ctrl_SI_matlab,
        'TI_MATLAB': baseline_TI_matlab,
        'TI_stim_MATLAB': stim_TI_matlab,
        'TI_ctrl_MATLAB': ctrl_TI_matlab,
        'prof_mean': baseline_mean,
        'prof_sem': baseline_sem,
        'prof_stim_mean': stim_mean,
        'prof_stim_sem': stim_sem,
        'prof_ctrl_mean': ctrl_mean,
        'prof_ctrl_sem': ctrl_sem,
        'prof_mean_MATLAB': baseline_mean_matlab,
        'prof_sem_MATLAB': baseline_sem_matlab,
        'prof_stim_mean_MATLAB': stim_mean_matlab,
        'prof_stim_sem_MATLAB': stim_sem_matlab,
        'prof_ctrl_mean_MATLAB': ctrl_mean_matlab,
        'prof_ctrl_sem_MATLAB': ctrl_sem_matlab,
        'prof_good_mean': good_mean,
        'prof_good_sem': good_sem,
        'prof_bad_mean': bad_mean,
        'prof_bad_sem': bad_sem,
        'prof_good_mean_matlab': good_mean_matlab,
        'prof_good_sem_matlab': good_sem_matlab,
        'prof_bad_mean_matlab': bad_mean_matlab,
        'prof_bad_sem_matlab': bad_sem_matlab,
        'MI': MI,
        'MI_early': MI_early,
        'MI_late': MI_late,
        'MI_shuf': MI_shuf,
        'MI_shuf_early': MI_early_shuf,
        'MI_shuf_late': MI_late_shuf
    }

    return session_inputs['clu_list'][clu], row

def process_recording(pathname, show_cell_progress=False):
    recname = Path(pathname).name
    t0 = time()
    if pathname in rec_list.pathHPCLCopt:
        prefix = 'HPCLC'
    else:
        prefix = 'HPCLCterm'
    session_inputs = load_session_inputs(pathname, recname, prefix)
    session_rows = {}

    for clu in tqdm(
        range(len(session_inputs['cell_identities'])),
        desc=f'{recname} profiles',
        disable=not show_cell_progress
    ):
        cluname, row = build_cell_profile_row(clu, session_inputs)
        session_rows[cluname] = [row[column] for column in PROFILE_COLUMNS]

    session_df = pd.DataFrame.from_dict(
        session_rows,
        orient='index',
        columns=PROFILE_COLUMNS,
    ).astype('object')
    return recname, session_df, time() - t0


#%% main

def main(argv=None):
    parser = argparse.ArgumentParser(description='generate HPC profiles')
    parser.add_argument(
        '--n-workers',
        '--workers',
        dest='n_workers',
        type=int,
        default=max(1, min(6, (os.cpu_count() or 1) - 1)),
        help='number of recordings to process in parallel; use 1 for sequential processing',
    )
    args = parser.parse_args(argv)

    print('not using GPU acceleration; profiles are CPU/session parallel\n')

    df = pd.DataFrame({column: pd.Series(dtype='object') for column in PROFILE_COLUMNS})
    recording_paths = rec_list.pathHPCLCopt + rec_list.pathHPCLCtermopt
    n_workers = min(max(1, args.n_workers), len(recording_paths))
    if n_workers == 1:
        for pathname in recording_paths:
            print_session(Path(pathname).name)
            recname, session_df, elapsed = process_recording(
                pathname,
                show_cell_progress=True,
            )
            df = pd.concat([df, session_df], axis=0)
            save_profiles(df)
            print(f'{recname} done in {timedelta(seconds=int(elapsed))}\n')
    else:
        print(f'processing {len(recording_paths)} recordings with {n_workers} workers')
        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            future_to_recname = {
                executor.submit(process_recording, pathname, False): Path(pathname).name
                for pathname in recording_paths
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
