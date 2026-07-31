# -*- coding: utf-8 -*-
'''
Created on Thu Mar 27 16:14:04 2025

cell metrics and trial data for the hippocampal ephys analyses

@author: Dinghao Luo
'''

#%% imports
from pathlib import Path

import numpy as np
import scipy.io as sio
import mat73

PYRUP_CLASS = 'run-onset ON'
PYRDOWN_CLASS = 'run-onset OFF'
UNRESPONSIVE_CLASS = 'run-onset unresponsive'


#%% cell metrics and trial data
def calculate_occupancy(speeds, dt, distance_bins):
    '''calculate spatial-bin occupancy from speed samples and a fixed time step'''
    speeds = np.asarray(speeds)
    cumulative_distance = np.cumsum(speeds * dt)
    bin_indices = np.digitize(cumulative_distance, distance_bins) - 1
    occupancy = np.bincount(bin_indices, minlength=len(distance_bins)) * dt

    return occupancy

def classify_run_onset_activation_ratio(train,
                                        run_onset_activated_thres,
                                        run_onset_inhibited_thres,
                                        pre_window=(-1.5, -0.5),
                                        run_onset_bin=3750,
                                        SAMP_FREQ=1250):
    '''
    classify run-onset modulation from the pre/post firing-rate ratio

    pre_window is in seconds relative to run onset.
    '''
    pre_start_s, pre_stop_s = pre_window
    pre_start = int(run_onset_bin + SAMP_FREQ * pre_start_s)
    pre_stop = int(run_onset_bin + SAMP_FREQ * pre_stop_s)
    pre = np.nanmean(train[pre_start:pre_stop])
    post = np.nanmean(train[int(run_onset_bin+SAMP_FREQ*.5):int(run_onset_bin+SAMP_FREQ*1.5)])
    ratio = pre/post
    if ratio < run_onset_activated_thres:
        ratiotype = PYRUP_CLASS
    elif ratio > run_onset_inhibited_thres:
        ratiotype = PYRDOWN_CLASS
    else:
        ratiotype = UNRESPONSIVE_CLASS

    return ratio, ratiotype

def compute_modulation_index(ctrl,
                             stim,
                             span=4,
                             run_onset_bin=3750,
                             SAMP_FREQ=1250):
    '''compute full, early, and late stimulation-control modulation indices'''
    ctrl_full = np.nanmean(ctrl[run_onset_bin:run_onset_bin+SAMP_FREQ*span])
    stim_full = np.nanmean(stim[run_onset_bin:run_onset_bin+SAMP_FREQ*span])
    MI_full = (stim_full-ctrl_full) / (stim_full+ctrl_full)

    demarc = int(run_onset_bin+SAMP_FREQ*span/2)

    ctrl_early = np.nanmean(ctrl[run_onset_bin:demarc])
    stim_early = np.nanmean(stim[run_onset_bin:demarc])
    MI_early = (stim_early-ctrl_early) / (stim_early+ctrl_early)
    ctrl_late = np.nanmean(ctrl[demarc:run_onset_bin+SAMP_FREQ*span])
    stim_late = np.nanmean(stim[demarc:run_onset_bin+SAMP_FREQ*span])
    MI_late = (stim_late-ctrl_late) / (stim_late+ctrl_late)

    return MI_full, MI_early, MI_late

def compute_modulation_index_shuf(ctrl_matrix,
                                  stim_matrix,
                                  span=4,
                                  bootstrap=100,
                                  run_onset_bin=3750,
                                  SAMP_FREQ=1250):
    '''compute modulation indices after shuffling control and stimulation trials'''
    ctrl_matrix = np.asarray(ctrl_matrix)
    stim_matrix = np.asarray(stim_matrix)
    tot_trials = ctrl_matrix.shape[0]

    pooled_matrix = np.vstack((ctrl_matrix, stim_matrix))
    pooled_idx = np.arange(ctrl_matrix.shape[0]+stim_matrix.shape[0])
    shuf_ctrl_idx = np.zeros((bootstrap, tot_trials), dtype=int)
    shuf_stim_idx = np.zeros((bootstrap, tot_trials), dtype=int)
    for i in range(bootstrap):  # shuffle n times
        shuf = np.random.permutation(pooled_idx)
        shuf_ctrl_idx[i, :] = shuf[:tot_trials]
        shuf_stim_idx[i, :] = shuf[tot_trials:]
    shuf_ctrl_mean = pooled_matrix[shuf_ctrl_idx, :].mean(axis=0).mean(axis=0)
    shuf_stim_mean = pooled_matrix[shuf_stim_idx, :].mean(axis=0).mean(axis=0)

    ctrl_full = np.nanmean(shuf_ctrl_mean[run_onset_bin:run_onset_bin+SAMP_FREQ*span])
    stim_full = np.nanmean(shuf_stim_mean[run_onset_bin:run_onset_bin+SAMP_FREQ*span])
    MI_full = (stim_full-ctrl_full) / (stim_full+ctrl_full)

    demarc = int(run_onset_bin+SAMP_FREQ*span/2)

    ctrl_early = np.nanmean(shuf_ctrl_mean[run_onset_bin:demarc])
    stim_early = np.nanmean(shuf_stim_mean[run_onset_bin:demarc])
    MI_early = (stim_early-ctrl_early) / (stim_early+ctrl_early)
    ctrl_late = np.nanmean(shuf_ctrl_mean[demarc:run_onset_bin+SAMP_FREQ*span])
    stim_late = np.nanmean(shuf_stim_mean[demarc:run_onset_bin+SAMP_FREQ*span])
    MI_late = (stim_late-ctrl_late) / (stim_late+ctrl_late)

    return MI_full, MI_early, MI_late  # this is shuffled

def compute_spatial_information(spike_counts, occupancy,
                                GPU_AVAILABLE=False):
    '''compute Skaggs spatial information in bits per spike'''
    # convert cupy arrays back before numpy indexing
    if GPU_AVAILABLE:
        occupancy = occupancy.get()

    valid_bins = occupancy > 0

    spike_counts = spike_counts[valid_bins]
    occupancy = occupancy[valid_bins]

    total_time = np.sum(occupancy)
    p_x = occupancy / total_time

    lambda_x = spike_counts / occupancy

    lambda_bar = np.sum(lambda_x * p_x)

    nonzero = lambda_x > 0
    spatial_info = np.nansum(
        p_x[nonzero]
        * (lambda_x[nonzero] / lambda_bar)
        * np.log2(lambda_x[nonzero] / lambda_bar)
    )

    return spatial_info

def compute_temporal_information(spike_times, bin_size_steps):
    '''compute temporal information for a neuron sampled at 1250 Hz'''
    spike_times = np.asarray(spike_times)

    total_steps = len(spike_times)
    num_bins = total_steps // bin_size_steps
    bin_edges = np.arange(0, total_steps + bin_size_steps, bin_size_steps)

    spike_counts, _ = np.histogram(spike_times, bins=bin_edges)

    # equal-width bins have uniform occupancy
    p_t = np.ones(num_bins) / num_bins

    lambda_t = spike_counts / (bin_size_steps / 1250)  # convert bin size to seconds

    lambda_bar = np.sum(lambda_t * p_t)

    nonzero = lambda_t > 0
    temporal_info = np.nansum(
        p_t[nonzero]
        * (lambda_t[nonzero] / lambda_bar)
        * np.log2(lambda_t[nonzero] / lambda_bar)
    )

    return temporal_info

def compute_trial_by_trial_variability(train):
    '''compute trial-wise variability as one minus median pairwise correlation'''
    if len(train) == 0:
        return np.nan

    max_length = max(len(trace) for trace in train)
    train = [np.asarray(trace[:max_length], dtype=float) for trace in train]

    num_trials = len(train)
    corr_matrix = np.full((num_trials, num_trials), np.nan)

    for i in range(num_trials):
        for j in range(i + 1, num_trials):
            if np.nanstd(train[i]) == 0 or np.nanstd(train[j]) == 0:
                corr_matrix[i, j] = np.nan
            else:
                corr_matrix[i, j] = np.corrcoef(train[i], train[j])[0, 1]

    corr_values = corr_matrix[np.triu_indices(num_trials, k=1)]
    if corr_values.size == 0 or np.all(np.isnan(corr_values)):
        return np.nan

    return 1 - np.nanmedian(corr_values)

def get_cell_info(info_filename):
    '''
    load cell type identities and spike rates from a MATLAB info file.

    parameters:
    - info_filename: path to the MATLAB .mat file containing cell information.

    returns:
    - cell_identities: list of strings ('pyr' or 'int') indicating cell types.
    - spike_rates: 1d array of spike rates for all cells.
    '''
    info = sio.loadmat(str(info_filename))
    # rec_info = info['rec'][0][0]
    autocorr = info['autoCorr'][0][0]

    # spike_rates = rec_info['firingRate'][0]

    # use the pipeline-calculated FR
    filestem = str(info_filename).split('_Info')[0]
    spike_rates = sio.loadmat(
        rf'{filestem}_FR_Run1.mat'
        )['mFRStruct']['mFR'][0][0].flatten()

    is_pyr = autocorr['isPyrneuron'][0]
    cell_identities = ['putative_pyr' if i else 'int' for i in is_pyr]

    return cell_identities, spike_rates

def resolve_behpar_path(pathname, maze_sess=None):
    '''
    find the MATLAB behPar file for a session or candidate behPar filename.

    Raphi's files can be named without the maze-session suffix, while the
    standard HPC files usually use `_behPar_msessN.mat`.
    '''
    pathname = Path(pathname)

    if pathname.suffix == '.mat':
        candidates = [pathname]
        name = pathname.name

        if '_behPar_msess' in name:
            stem = name.split('_behPar_msess')[0]
            candidates.append(pathname.with_name(f'{stem}_behPar.mat'))
        elif name.endswith('_behPar.mat'):
            stem = name[:-len('_behPar.mat')]
            candidates.extend([
                pathname.with_name(f'{stem}_behPar_msess{sess}.mat')
                for sess in range(10)
            ])
    else:
        recname = pathname.name
        if maze_sess is None:
            candidates = [
                pathname / f'{recname}_DataStructure_mazeSection1_TrialType1_behPar_msess{sess}.mat'
                for sess in range(10)
            ]
        else:
            candidates = [
                pathname / f'{recname}_DataStructure_mazeSection1_TrialType1_behPar_msess{maze_sess}.mat'
            ]

        candidates.append(
            pathname / f'{recname}_DataStructure_mazeSection1_TrialType1_behPar.mat'
        )

    for candidate in candidates:
        if candidate.exists():
            return candidate

    raise FileNotFoundError(
        'could not find MATLAB behPar file; checked:\n'
        + '\n'.join(str(candidate) for candidate in candidates)
    )

def get_good_bad_idx_MATLAB(pathname):
    '''
    extract indices of good and bad trials from a MATLAB behavioural parameter file.

    parameters:
    - pathname: full path to the session folder containing the MATLAB behaviour file.

    returns:
    - good_idx_matlab: list of indices (0-based) for trials marked as good.
    - bad_idx_matlab: list of indices (0-based) for trials marked as bad.
    '''
    beh_parameter_file_path = resolve_behpar_path(pathname)
    beh_parameter_file = sio.loadmat(beh_parameter_file_path)

    # same as the previous function
    bad_idx_matlab = [trial-1 for trial, quality
                      in enumerate(beh_parameter_file['behPar'][0]['indTrBadBeh'][0][0])
                      if quality and trial>0]
    good_idx_matlab = [trial-1 for trial, quality
                       in enumerate(beh_parameter_file['behPar'][0]['indTrBadBeh'][0][0])
                       if not quality and trial>0]

    return good_idx_matlab, bad_idx_matlab

def get_trial_matrix(trains, trialtype_idx, max_samples, clu):
    '''
    get the trial matrix for a given cluster and trial type indices.

    parameters:
    - trains: list of spike trains for all clusters.
    - trialtype_idx: list of trial indices to include.
    - max_samples: maximum number of samples per trial.
    - clu: cluster identifier.

    returns:
    - temp_matrix: matrix of spike trains for the specified trials and cluster.
    '''
    if len(trialtype_idx)==0:  # if there is no trial in the list
        return np.nan

    temp_matrix = np.zeros((len(trialtype_idx), max_samples))
    for idx, trial in enumerate(trialtype_idx):
        try:
            trial_length = len(trains[clu][trial])
        except TypeError:
            trial_length = 0
        if 0 < trial_length < max_samples:
            temp_matrix[idx, :trial_length] = np.asarray(trains[clu][trial][:])
        elif trial_length > 0:
            temp_matrix[idx, :] = np.asarray(trains[clu][trial][:max_samples])
    return temp_matrix

def get_trialtype_idx_MATLAB(beh_filename):
    '''
    get ctrl/stim trial indices from a MATLAB behPar file.

    parameters:
    - beh_filename: path to the MATLAB behaviour file.

    returns:
    - baseline_idx: indices for baseline trials.
    - stim_idx: indices for stimulation trials.
    - ctrl_idx: indices for control trials.
    '''
    beh_par = sio.loadmat(str(resolve_behpar_path(beh_filename)))
    stim_on = beh_par['behPar']['stimOn'][0][0][0]
    max_length = len(stim_on) - 1
    stim_idx = np.where(stim_on != 0)[0]

    if len(stim_idx) == 0:
        return np.arange(1, len(stim_on)), np.array([], dtype=int), np.array([], dtype=int)

    baseline_idx = np.arange(1, stim_idx[0])
    ctrl_idx = stim_idx + 2
    ctrl_idx = ctrl_idx[ctrl_idx < max_length]
    return baseline_idx, stim_idx, ctrl_idx

def load_speeds(beh_series):
    speed_times = beh_series['speed_times_aligned']
    new_speed_times = []
    for trial in range(1, len(speed_times)):  # trial 1 is empty and not included in spike trains
        curr_speed_times = speed_times[trial]
        curr_aligned = [s[1] for s in curr_speed_times]
        new_speed_times.append(curr_aligned)
    return new_speed_times

def load_dist_spike_array(dist_filename):
    dist_mat = mat73.loadmat(str(dist_filename))['filteredSpikeDistArrayRun']
    trains_dist = []  # use a list to mirror the structure of trains.npy
    for clu in range(len(dist_mat)):
        trains_dist.append(dist_mat[clu][1:])  # trial 1 is empty
    return trains_dist
