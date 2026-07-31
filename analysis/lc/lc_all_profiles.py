# -*- coding: utf-8 -*-
'''
Created on Thu Jun  8 09:26:49 2023

compiling cell properties into a dataframe

@author: Dinghao Luo

'''

#%% imports
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import timedelta
import os
from pathlib import Path
from time import time
import numpy as np
import pandas as pd
import scipy.io as sio
import sys
from scipy.stats import pearsonr, sem
from tqdm import tqdm

repo_root = Path(__file__).resolve().parents[2]
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))

import peak_detection_functions as pdf  # once i imported this as pd... stupid
from common_functions import gaussian_kernel_unity, get_trialtype_indices_from_stim_conds
from console_formatting import print_files_saved
import project_paths as pp

import rec_list

paths = rec_list.pathLC

OUTPUT_PATH = pp.LC_EPHYS_STEM / 'LC_all_cell_profiles.pkl'
PLOT_PAYLOAD_PATH = pp.LC_EPHYS_STEM / 'LC_all_profiles_plot_payload.npy'
PROFILE_COLUMNS = [
    'sessname',
    'identity',
    'spike_rate',
    'acg',
    'waveform',
    'run_onset_peak',
    'speed_rate_r',
    'speed_rate_p',
    'lick_sensitive',
    'lick_sensitive_type',
    'lick_sensitive_signif',
    'mean_profile',
    'sem_profile',
    'baseline_mean',
    'baseline_sem',
    'stim_mean',
    'stim_sem',
    'ctrl_mean',
    'ctrl_sem'
    ]


#%% GPU acceleration
# CPU initially looked faster because this loop was not parallelised enough;
# the present GPU path runs one cell in about 2 s instead of roughly 30 s
try:
    import cupy as cp
    GPU_AVAILABLE = cp.cuda.runtime.getDeviceCount() > 0
except ModuleNotFoundError:
    GPU_AVAILABLE = False
    GPU_MESSAGE = 'GPU-acceleration unavailable'

if GPU_AVAILABLE:
    xp = cp
    from common_functions import sem_gpu as sem
    name = cp.cuda.runtime.getDeviceProperties(0)['name'].decode('UTF-8')
    GPU_MESSAGE = f'GPU-acceleration with {str(name)} and cupy'
else:
    xp = np
    from scipy.stats import sem
    GPU_MESSAGE = 'GPU-acceleration unavailable'


#%% per-cell calculations


def compute_lick_info(beh_file: np.ndarray,
                      align_run_file: np.ndarray,
                      samp_freq=1250) -> tuple:
    '''
    extracts lick timing information and filters trials based on behaviour.

    parameters:
    - beh_file (np.ndarray): behavioural data file containing trial parameters
    - align_run_file (np.ndarray): file containing aligned running and licking data
    - samp_freq (int): sampling frequency in Hz, default 1250

    returns:
    - tuple:
        1. int: total number of good non-stimulation trials
        2. list: list of indices for good non-stimulation trials
        3. xp.ndarray: end indices of lick-aligned trials
        4. xp.ndarray: first lick times relative to trial start
    '''
    licks = align_run_file['trialsRun']['lickLfpInd'][0][0][0][1:]
    starts = align_run_file['trialsRun']['startLfpInd'][0][0][0][1:]
    tot_trial = licks.shape[0]

    stim_on = np.asarray(beh_file['behPar']['stimOn'][0][0][0][1:])
    stim_idx = np.where(stim_on != 0)[0] + 1

    bad = np.asarray(beh_file['behPar'][0]['indTrBadBeh'][0][0])
    bad_idx = list(np.where(bad == 1)[0] - 1) if bad.sum() != 0 else []

    first_licks = []
    for trial in range(tot_trial):
        lk = [l for l in licks[trial] if l - starts[trial] > samp_freq]
        if len(lk) == 0:
            first_licks.append(np.nan)
            bad_idx.append(trial)
        else:
            first_licks.extend(lk[0] - starts[trial])

    trial_list = [
        trial for trial in np.arange(tot_trial)
        if trial not in bad_idx and trial not in stim_idx
    ]
    tot_trial_good_nonstim = len(trial_list)

    first_licks = xp.asarray([first_licks[trial] for trial in trial_list])
    end_indices = first_licks + 3 * samp_freq + 3 * samp_freq

    return tot_trial_good_nonstim, trial_list, end_indices, first_licks

def compute_lick_sensitivity(cluname: str,
                             trains: xp.ndarray,
                             rasters: xp.ndarray,
                             identity: str,
                             tot_trial_good_nonstim: int,
                             trial_list: list,
                             end_indices: xp.ndarray,
                             first_licks: xp.ndarray,
                             plot_payloads,
                             samp_freq=1250,
                             around=6,
                             bootstrap=1000,
                             GPU_AVAILABLE=False) -> tuple:
    '''
    computes the lick sensitivity of a cell based on spike trains and rasters.

    parameters:
    - cluname (str): identifier for the cell
    - trains (xp.ndarray): spike train data (trials x timepoints)
    - rasters (xp.ndarray): raster data (trials x timepoints)
    - identity (str): identity of the cell ('tagged', 'putative', etc.)
    - tot_trial_good_nonstim (int): number of good non-stimulation trials
    - trial_list (list): indices of good non-stimulation trials
    - end_indices (xp.ndarray): end indices of lick-aligned trials
    - first_licks (xp.ndarray): first lick times relative to trial start
    - samp_freq (int): sampling frequency in Hz, default 1250
    - around (int): time window in seconds, default 6
    - bootstrap (int): number of bootstrap iterations, default 1000
    - GPU_AVAILABLE (bool): if true, uses GPU for calculations, default false

    returns:
    - tuple:
        1. bool: whether the cell is lick-sensitive
        2. str: type of sensitivity ('ON', 'OFF', or unresponsive)
        3. str: significance of sensitivity ('***', '**', '*', or 'n.s.')
    '''
    aligned_prof = xp.zeros((tot_trial_good_nonstim, samp_freq * around))
    aligned_rasters = xp.zeros((tot_trial_good_nonstim, samp_freq * around))
    for i, trial in enumerate(trial_list):
        train = xp.asarray(trains[trial])
        raster = xp.asarray(rasters[trial])
        if end_indices[i] <= len(train):
            aligned_prof[i, :] = train[first_licks[i]:first_licks[i] + 6 * samp_freq]
        else:
            aligned_prof[i, :] = xp.pad(
                train[first_licks[i]:],
                (0, 6 * samp_freq - len(train[first_licks[i]:])),
                mode='constant'
            )
        if end_indices[i] <= len(raster):
            aligned_rasters[i, :] = raster[first_licks[i]:first_licks[i] + 6 * samp_freq]
        else:
            aligned_rasters[i, :] = xp.pad(
                raster[first_licks[i]:],
                (0, 6 * samp_freq - len(raster[first_licks[i]:])),
                mode='constant'
            )

    # bootstrap post/pre ratios after independent circular trial shifts
    bootstrap_xp = cp if GPU_AVAILABLE else np
    length = 6 * 1250
    shuf_ratio = bootstrap_xp.zeros(bootstrap)
    indices = bootstrap_xp.arange(length)
    for shuf in range(bootstrap):
        rand_shifts = bootstrap_xp.random.randint(1, length, tot_trial_good_nonstim)
        shifted_indices = (indices[None, :] - rand_shifts[:, None]) % length
        shuf_arr = aligned_prof[
            bootstrap_xp.arange(tot_trial_good_nonstim)[:, None],
            shifted_indices,
        ]
        shuf_result = bootstrap_xp.mean(shuf_arr, axis=0)
        shuf_ratio[shuf] = (
            bootstrap_xp.sum(shuf_result[length // 2:length // 2 + 1250]) /
            bootstrap_xp.sum(shuf_result[length // 2 - 1250:length // 2])
        )
    # percentile thresholds: [99.9, 99, 95, 50, 5, 1, .1]
    shuf_ratios = bootstrap_xp.percentile(
        shuf_ratio, [99.9, 99, 95, 50, 5, 1, .1], axis=0
    )

    aligned_prof_mean = xp.nanmean(aligned_prof, axis=0)
    aligned_prof_sem = sem(aligned_prof, axis=0)
    true_ratio = xp.sum(
        aligned_prof_mean[3 * samp_freq: 3 * samp_freq + 1 * samp_freq]
    ) / xp.sum(
        aligned_prof_mean[3 * samp_freq - 1 * samp_freq: 3 * samp_freq]
    )

    if true_ratio >= shuf_ratios[2]:
        lick_sensitive = True
        lick_sensitive_type = 'ON'
        suffix = 'lick-ON'
        for i, ratio in enumerate(shuf_ratios[:3]):
            if true_ratio > ratio:
                break
        if i == 0:
            lick_sensitive_signif = '***'
        if i == 1:
            lick_sensitive_signif = '**'
        if i == 2:
            lick_sensitive_signif = '*'
    elif true_ratio <= shuf_ratios[-3]:
        lick_sensitive = True
        lick_sensitive_type = 'OFF'
        suffix = 'lick-OFF'
        for i, ratio in enumerate(reversed(shuf_ratios[4:])):
            if true_ratio < ratio:
                break
        if i == 0:
            lick_sensitive_signif = '***'
        if i == 1:
            lick_sensitive_signif = '**'
        if i == 2:
            lick_sensitive_signif = '*'
    else:
        lick_sensitive = False
        lick_sensitive_type = np.nan
        suffix = 'unresponsive'
        lick_sensitive_signif = 'n.s.'

    raster_arr = [
        (raster.nonzero()[0] - 3 * samp_freq).get() / samp_freq if GPU_AVAILABLE else (raster.nonzero()[0] - 3 * samp_freq) / samp_freq
        for raster in aligned_rasters
    ]

    plot_payloads['lick_sensitivity'][cluname] = {
        'cluname': cluname,
        'suffix': suffix,
        'identity': identity,
        'lick_sensitive_signif': lick_sensitive_signif,
        'aligned_rasters': raster_arr,
        'aligned_prof_mean': aligned_prof_mean.get() if GPU_AVAILABLE else np.asarray(aligned_prof_mean),
        'aligned_prof_sem': aligned_prof_sem.get() if GPU_AVAILABLE else np.asarray(aligned_prof_sem),
        'shuf_ratios': shuf_ratios.get() if GPU_AVAILABLE else np.asarray(shuf_ratios),
        'true_ratio': float(true_ratio.get() if GPU_AVAILABLE else true_ratio)
        }

    return lick_sensitive, lick_sensitive_type, lick_sensitive_signif

def compute_speeds(align_run_file: xp.ndarray,
                   samp_freq=1250) -> xp.ndarray:
    '''
    compute the speed of all trials and apply gaussian smoothing.

    parameters:
    - align_run_file: the file containing aligned running data.
    - samp_freq: the sampling frequency in Hz, default is 1250.

    returns:
    - an array of truncated and smoothed speed values for each trial.
    '''
    speed_time_bef = align_run_file['trialsRun'][0]['speed_MMsecBef'][0][0][1:]
    speed_time = align_run_file['trialsRun'][0]['speed_MMsec'][0][0][1:]

    speed_time_all = np.empty(shape=speed_time.shape[0], dtype='object')
    for i in range(speed_time.shape[0]):
        bef = speed_time_bef[i]
        aft = speed_time[i]
        speed_time_all[i] = np.concatenate([bef, aft])
        speed_time_all[i][speed_time_all[i] < 0] = 0
    gaus_speed = gaussian_kernel_unity(samp_freq / 100)
    speed_time_conv = [np.convolve(np.squeeze(trial), gaus_speed, mode='same')
                       for trial in speed_time_all]

    speed_trunc = np.zeros((len(speed_time_conv), 5 * samp_freq))
    for trial, speeds in enumerate(speed_time_conv):
        if len(speeds) > (3 + 4) * samp_freq:
            speed_trunc[trial, :] = speeds[2 * samp_freq: 7 * samp_freq]
        else:
            speed_trunc[trial, :len(speeds) - 2 * samp_freq] = speeds[2 * samp_freq:]

    return speed_trunc

def compute_speed_rate_corr(trains: xp.ndarray,
                            speed_trunc: xp.ndarray,
                            samp_freq=1250) -> tuple:
    '''
    calculate the correlation between running speed and spike rate for a single cell.

    parameters:
    - cluname: identifier for the cell.
    - trains: spike train data (trials x timepoints).
    - speed_trunc: truncated and smoothed speed data (trials x timepoints).
    - samp_freq: sampling frequency in Hz, default is 1250.

    returns:
    - tuple: (mean correlation coefficient, mean p-value).
    '''
    train_trunc = np.zeros((len(trains), 5 * samp_freq))
    for trial, train in enumerate(trains):
        if len(train) > (3 + 4) * samp_freq:
            train_trunc[trial, :] = train[2 * samp_freq: 7 * samp_freq]
        else:
            train_trunc[trial, :len(train) - 2 * samp_freq] = train[2 * samp_freq:]

    rate_speed_corr = np.zeros((len(trains), 2))
    for trial in range(len(trains)):
        rate_speed_corr[trial, :] = pearsonr(speed_trunc[trial], train_trunc[trial])

    return np.mean(rate_speed_corr[:, 0]), np.mean(rate_speed_corr[:, 1])

def get_trial_matrix(trains, trialtype_idx, max_samples=1250 * 8):
    '''
    get the trial matrix for a given trial type indices.

    parameters:
    - trains: list of spike trains for one cell.
    - trialtype_idx: list of trial indices to include.
    - max_samples: maximum number of samples per trial.

    returns:
    - temp_matrix: matrix of spike trains for the specified trials.
    '''
    if len(trialtype_idx) == 0:
        return xp.nan

    temp_matrix = xp.zeros((len(trialtype_idx), max_samples))

    trial_lengths = [
        len(trains[t]) if isinstance(trains[t], (list, np.ndarray)) else 0
        for t in trialtype_idx
    ]

    for idx, trial in enumerate(trialtype_idx):
        if 0 < trial_lengths[idx] < max_samples:
            temp_matrix[idx, :trial_lengths[idx]] = xp.array(trains[trial][:])
        elif trial_lengths[idx] > 0:
            temp_matrix[idx, :] = xp.array(trains[trial][:max_samples])
    return temp_matrix

def process_session(path: str,
                    behaviour: pd.DataFrame,
                    kmeans: dict,
                    max_sample=1250 * 8,
                    verbose=False):
    recname = path[-17:]
    t0 = time()
    log_lines = [recname]
    plot_payloads = {
        'peak_detection': {},
        'lick_sensitivity': {}
        }

    beh = behaviour.loc[recname]
    stim_conds = [trial[15] for trial in beh['trial_statements']][1:]
    baseline_idx, stim_idx, ctrl_idx = get_trialtype_indices_from_stim_conds(stim_conds)

    session_stem = pp.LC_EPHYS_STEM / 'all_sessions' / recname
    all_trains = np.load(
        session_stem / f'{recname}_all_trains.npy',
        allow_pickle=True,
    ).item()
    all_rasters = np.load(
        session_stem / f'{recname}_all_rasters.npy',
        allow_pickle=True,
    ).item()
    all_identities = np.load(
        rf'{session_stem}\{recname}_all_identities.npy',
        allow_pickle=True
    ).item()
    all_acgs = np.load(
        rf'{session_stem}\{recname}_all_ACGs.npy',
        allow_pickle=True
    ).item()
    all_waveforms = np.load(
        rf'{session_stem}\{recname}_all_waveforms.npy',
        allow_pickle=True
    ).item()

    raw_session_stem = pp.MICEEXP_ROOT / f'ANMD{recname[1:5]}' / recname[:14] / recname
    spike_rate_file = sio.loadmat(
        raw_session_stem / f'{recname}_DataStructure_mazeSection1_TrialType1_FR_Ctrl_Run0_mazeSess1.mat'
    )
    beh_file = sio.loadmat(
        raw_session_stem / f'{recname}_DataStructure_mazeSection1_TrialType1_behPar_msess1.mat'
    )
    align_run_file = sio.loadmat(
        raw_session_stem / f'{recname}_DataStructure_mazeSection1_TrialType1_alignRun_msess1.mat'
    )
    speeds = compute_speeds(align_run_file)

    tagged_keys = [clu for clu in all_identities if all_identities[clu] == 1]
    log_lines.append(f'{recname} ({len(all_trains)} cells)')
    session_rows = {}

    for cluname in all_trains:
        if verbose:
            log_lines.append(cluname)
        trains = all_trains[cluname]
        rasters = all_rasters[cluname]
        waveform = all_waveforms[cluname]
        acg = all_acgs[cluname]

        mean_profile = xp.mean(xp.array([train[:max_sample] for train in trains]), axis=0)
        sem_profile = sem(xp.array([train[:max_sample] for train in trains]), axis=0)
        baseline_mean = xp.mean(get_trial_matrix(trains, baseline_idx), axis=0)
        baseline_sem = sem(get_trial_matrix(trains, baseline_idx), axis=0)
        if GPU_AVAILABLE:
            mean_profile = mean_profile.get()
            sem_profile = sem_profile.get()
            baseline_mean = baseline_mean.get()
            baseline_sem = baseline_sem.get()

        if stim_idx:
            stim_mean = xp.mean(get_trial_matrix(trains, stim_idx), axis=0)
            stim_sem = sem(get_trial_matrix(trains, stim_idx), axis=0)
            ctrl_mean = xp.mean(get_trial_matrix(trains, ctrl_idx), axis=0)
            ctrl_sem = sem(get_trial_matrix(trains, ctrl_idx), axis=0)
            if GPU_AVAILABLE:
                stim_mean = stim_mean.get()
                stim_sem = stim_sem.get()
                ctrl_mean = ctrl_mean.get()
                ctrl_sem = ctrl_sem.get()
        else:
            stim_mean, ctrl_mean, stim_sem, ctrl_sem = [], [], [], []

        spike_rate_ctrl = spike_rate_file['mFRStructSessCtrl']['mFR'][0][0][0]
        clu_id = int(cluname.split('clu')[1])
        spike_rate = spike_rate_ctrl[clu_id - 2]
        # tagged first, then k-means; putative cells above 10 Hz are classed as other
        if cluname in tagged_keys:
            identity = 'tagged'
        elif kmeans[cluname] == 0:
            identity = 'other'
        elif kmeans[cluname] == 1:
            identity = 'putative' if spike_rate <= 10 else 'other'
        else:
            identity = None

        trial_types = beh_file['behPar']['stimOn'][0][0][0]
        first_stim = next((trial for trial, stim_on in enumerate(trial_types) if stim_on), -1)
        peak, mean_prof, shuf_prof = pdf.peak_detection(
            trains,
            first_stim=first_stim,
            around=4,
            peak_width=2,
            bootstrap=5000,
            GPU_AVAILABLE=GPU_AVAILABLE,
            VERBOSE=False
        )
        plot_payloads['peak_detection'][cluname] = {
            'cluname': cluname,
            'mean_prof': mean_prof,
            'shuf_prof': shuf_prof,
            'peak': peak,
            'identity': identity
            }

        speed_rate_r, speed_rate_p = compute_speed_rate_corr(trains, speeds)
        tot_trial_good_nonstim, trial_list, end_indices, first_licks = compute_lick_info(
            beh_file,
            align_run_file
        )
        lick_sensitive, lick_sensitivity_type, lick_sensitivity_signif = compute_lick_sensitivity(
            cluname,
            trains,
            rasters,
            identity,
            tot_trial_good_nonstim,
            trial_list,
            end_indices,
            first_licks,
            bootstrap=5000,
            GPU_AVAILABLE=GPU_AVAILABLE,
            plot_payloads=plot_payloads
        )

        session_rows[cluname] = [
            recname,
            identity,
            spike_rate,
            acg,
            waveform,
            peak,
            speed_rate_r,
            speed_rate_p,
            lick_sensitive,
            lick_sensitivity_type,
            lick_sensitivity_signif,
            mean_profile,
            sem_profile,
            baseline_mean,
            baseline_sem,
            stim_mean,
            stim_sem,
            ctrl_mean,
            ctrl_sem
            ]

    session_df = pd.DataFrame.from_dict(
        session_rows,
        orient='index',
        columns=PROFILE_COLUMNS,
        ).astype('object')
    return recname, session_df, plot_payloads, log_lines, time() - t0


#%% main loop


def main(argv=None):
    parser = argparse.ArgumentParser(description='compile LC cell profiles')
    parser.add_argument(
        '--n-workers',
        '--workers',
        dest='n_workers',
        type=int,
        default=min(6, os.cpu_count() - 1),
        help='number of sessions to process in parallel; use 1 for sequential processing',
        )
    parser.add_argument(
        '-v',
        '--verbose',
        '-verbose',
        action='store_true',
        help='print each session and cell',
        )
    args = parser.parse_args(argv)
    n_workers = args.n_workers
    fname = OUTPUT_PATH
    max_sample = 1250 * 8

    if args.verbose:
        print(GPU_MESSAGE)

    df = pd.DataFrame({column: pd.Series(dtype='object') for column in PROFILE_COLUMNS})
    behaviour_stem = pp.behaviour_experiment_stem('LC')
    sessions = {
        session_path.stem: pd.read_pickle(session_path)
        for session_path in sorted(behaviour_stem.glob('*.pkl'))
    }
    behaviour = pd.DataFrame.from_dict(sessions, orient='index')
    behaviour.index.name = 'recname'
    kmeans = np.load(
        pp.LC_EPHYS_STEM / 'UMAP' / 'LC_all_UMAP_kmeans.npy',
        allow_pickle=True
    ).item()
    plot_payloads = {
        'peak_detection': {},
        'lick_sensitivity': {}
        }

    n_workers = min(n_workers, len(paths))
    if n_workers == 1:
        for path in paths:
            recname, session_df, session_payloads, log_lines, elapsed = process_session(
                path,
                behaviour,
                kmeans,
                max_sample=max_sample,
                verbose=args.verbose
            )
            df = pd.concat([df, session_df], axis=0)
            for payload_name, payload in session_payloads.items():
                plot_payloads[payload_name].update(payload)
            if args.verbose:
                for line in log_lines:
                    print(line)

            df.to_pickle(fname)
            print_files_saved([('dataframe', fname)])
            np.save(PLOT_PAYLOAD_PATH, plot_payloads)
            print(f'{recname} done in {str(timedelta(seconds=int(elapsed)))}\n')
    else:
        print(f'processing {len(paths)} sessions with {n_workers} workers')
        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            futures = [
                executor.submit(
                    process_session,
                    path,
                    behaviour,
                    kmeans,
                    max_sample,
                    args.verbose,
                )
                for path in paths
            ]
            for future in tqdm(
                as_completed(futures),
                total=len(futures),
                desc='sessions',
            ):
                recname, session_df, session_payloads, log_lines, elapsed = future.result()
                df = pd.concat([df, session_df], axis=0)
                for payload_name, payload in session_payloads.items():
                    plot_payloads[payload_name].update(payload)
                if args.verbose:
                    for line in log_lines:
                        print(line)

                df.to_pickle(fname)
                print_files_saved([('dataframe', fname)])
                np.save(PLOT_PAYLOAD_PATH, plot_payloads)
                print(f'{recname} done in {str(timedelta(seconds=int(elapsed)))}\n')

if __name__ == '__main__':
    main(sys.argv[1:])
