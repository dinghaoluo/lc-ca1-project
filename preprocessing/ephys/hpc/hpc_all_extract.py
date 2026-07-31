# -*- coding: utf-8 -*-
'''
Created on Mon 10 July 10:02:32 2023
Modified on 24 Nov 2025 to mirror the LC processing pipeline

pool all cells from all recording sessions
modified 11 Dec 2024 to process with all trials (not skipping trial 0) and
    added GPU support
    - memory leakage problems on GPU, 20 Dec 2024
    - process each recording inside a function so its GPU references are
        destroyed afterwards, 26 Dec 2024

@author: Dinghao Luo

'''

#%% imports
import argparse
from pathlib import Path
import sys

import h5py
import scipy.io as sio
import mat73
from scipy.signal import fftconvolve
from tqdm import tqdm
from time import time
from datetime import timedelta
import gc

repo_root = Path(__file__).resolve().parents[3]
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))

from common_functions import gaussian_kernel_unity, get_GPU_availability
from console_formatting import print_files_saved, print_session
import project_paths as pp

import rec_list


#%% GPU acceleration
cp, GPU_AVAILABLE, _ = get_GPU_availability()

if GPU_AVAILABLE:
    import cupyx.scipy.signal as cpss  # for GPU
    import numpy as np
    mempool = cp.get_default_memory_pool()
    pinned_mempool = cp.get_default_pinned_memory_pool()
else:
    import numpy as np


#%% parameters
SAMP_FREQ   = 1250  # Hz
SIGMA_SPIKE = int(SAMP_FREQ * 0.05)  # 50 ms
GAUS_SPIKE  = gaussian_kernel_unity(SIGMA_SPIKE, GPU_AVAILABLE)
SIGMA_SPIKE_ALL_TRAINS = int(SAMP_FREQ * 0.03)  # 30 ms
GAUS_SPIKE_ALL_TRAINS = gaussian_kernel_unity(
    SIGMA_SPIKE_ALL_TRAINS,
    GPU_AVAILABLE
    )
GAUS_SPIKE_CPU = gaussian_kernel_unity(SIGMA_SPIKE, False)

MAX_LENGTH = 12500  # samples

BEF = 3  # seconds before
AFT = 7  # seconds after

MICEEXP_ROOT = pp.MICEEXP_ROOT
HPC_ALL_SESSIONS_STEM = pp.HPC_EPHYS_STEM / 'all_sessions'
ARCHIVE_EXTRA_RECNAMES = [
    'A063r-20230712-01',
    'A069r-20230914-02',
    ]
PROFILE_INPUT_SKIP_RECNAMES = {
    'A063r-20230712-01',
    'A069r-20230914-02',
    'A070r-20231108-01',
    'A076r-20231213-01',
    'A076r-20231214-01',
    'A083r-20240311-01',
    }

paths = (
    rec_list.pathHPCLCopt
    + rec_list.pathHPCLCtermopt
    + [
        str(MICEEXP_ROOT / f'ANMD{recname[1:5]}' / recname[:14] / recname)
        for recname in ARCHIVE_EXTRA_RECNAMES
    ]
    )


#%% main
def release_gpu_memory():
    cp.cuda.Stream.null.synchronize()
    mempool.free_all_blocks()
    pinned_mempool.free_all_blocks()

def convolve_raster(raster, gpu_kernel):
    raster_gpu = cp.asarray(raster)
    result = (
        cpss.fftconvolve(
            raster_gpu,
            gpu_kernel[None, None, :],
            mode='same'
        ) * SAMP_FREQ
    ).get()
    del raster_gpu
    release_gpu_memory()
    return result

def process_session(path, verbose=False):
    recname = Path(path).name

    rec_stem = MICEEXP_ROOT / f'ANMD{recname[1:5]}' / recname[:14] / recname

    # aligned behavioural landmarks; use the raw recording folder first, then the in-repo copy
    aligned_run_path = pp.resolve_matlab_pipeline_file(
        rec_stem / f'{recname}_DataStructure_mazeSection1_TrialType1_alignRun_msess1.mat',
        recname,
    )
    aligned_cue_path = pp.resolve_matlab_pipeline_file(
        rec_stem / f'{recname}_DataStructure_mazeSection1_TrialType1_alignCue_msess1.mat',
        recname,
    )
    aligned_rew_path = pp.resolve_matlab_pipeline_file(
        rec_stem / f'{recname}_DataStructure_mazeSection1_TrialType1_alignRew_msess1.mat',
        recname,
    )
    # spike file paths
    clu_paths = [rec_stem / f'{recname}.clu.{probe}' for probe in range(1,7)]
    res_paths = [rec_stem / f'{recname}.res.{probe}' for probe in range(1,7)]

    # get cluname
    filename = pp.resolve_matlab_pipeline_file(
        rec_stem / f'{recname}_BehavElectrDataLFP.mat',
        recname,
    )
    spikefile_name = pp.resolve_matlab_pipeline_file(
        rec_stem / f'{recname}_DataStructure_mazeSection1_TrialType1_alignedSpikesPerNPerT_msess1_Run0.mat',
        recname,
    )
    sess_stem = HPC_ALL_SESSIONS_STEM / recname
    if verbose:
        print_session(recname)

    BehavLFP = mat73.loadmat(filename)
    Clu = BehavLFP['Clu']
    shank = Clu['shank']
    localClu = Clu['localClu']

    # alignment timepoints
    aligned_run = sio.loadmat(aligned_run_path)['trialsRun'][0][0]
    run_onsets = aligned_run['startLfpInd'][0][1:]  # discard the first trial which is empty

    aligned_cue = sio.loadmat(aligned_cue_path)['trialsCue'][0][0]
    cue_onsets = aligned_cue['startLfpInd'][0][1:]  # similar to above

    aligned_rew = sio.loadmat(aligned_rew_path)['trialsRew'][0][0]
    rew_onsets = aligned_rew['startLfpInd'][0][1:]  # this marks the last trial's reward

    tot_trials = len(run_onsets)
    if len(rew_onsets) != tot_trials or len(cue_onsets) != tot_trials:
        if verbose:
            print('warning: onsets of different lengths')

    ## ---- spike reading ---- ##
    clusters    = np.loadtxt(clu_paths[0], dtype=int, skiprows=1)  # initiate
    spike_times = np.loadtxt(res_paths[0], dtype=int) / (20_000 / SAMP_FREQ)
    valid = (clusters != 0) & (clusters != 1)  # filter out the MUA and noise
    clusters    = clusters[valid]
    clusters    = np.array([clu - 1 for clu in clusters])
    spike_times = spike_times[valid]
    spike_times = spike_times.astype(int)

    for probe in range(1, 6):
        last_cluster = clusters.max()

        if not clu_paths[probe].exists() or not res_paths[probe].exists():
            continue
        new_clusters    = np.loadtxt(clu_paths[probe], dtype=int, skiprows=1)
        new_spike_times = np.loadtxt(res_paths[probe], dtype=int) / (20_000 / SAMP_FREQ)
        valid = (new_clusters != 0) & (new_clusters != 1)
        new_clusters    = new_clusters[valid]
        new_clusters    = [clu - 1 for clu in new_clusters]
        new_clusters    = np.array([clu + last_cluster for clu in new_clusters])
        new_spike_times = new_spike_times[valid]
        new_spike_times = new_spike_times.astype(int)

        clusters    = np.concatenate((clusters, new_clusters))
        spike_times = np.concatenate((spike_times, new_spike_times))
    ## ---- spike reading ends ---- ##
    unique_clus = [clu for clu in np.unique(clusters)]
    clu_to_row = {clu: i for i, clu in enumerate(unique_clus)}  # map cluster ID to row index

    max_time = spike_times.max() + 1  # +1 to make sure last time index is included
    spike_map = np.zeros((len(unique_clus), max_time), dtype=int)

    for t, clu in zip(spike_times, clusters):
        row = clu_to_row[clu]
        spike_map[row, t] = 1  # set spike bin to 1

    # spike times are now in milliseconds

    all_trains_run = {}
    all_trains_cue = {}
    all_trains_rew = {}
    all_trains = {}
    all_rasters_run = {}
    all_rasters_cue = {}
    all_rasters_rew = {}
    all_rasters = {}

    spike_file_handle = h5py.File(spikefile_name, 'r')
    spike_time_file = spike_file_handle['trialsRunSpikes']

    time_bef = spike_time_file['TimeBef']
    time_aft = spike_time_file['Time']
    tot_clu = time_aft.shape[1]
    tot_trial = time_aft.shape[0]-1  # trial 1 is empty

    rasters = np.zeros((tot_clu, tot_trial, MAX_LENGTH), dtype=np.uint16)
    rasters_run = np.zeros((tot_clu, tot_trial, MAX_LENGTH), dtype=np.uint16)
    rasters_rew = np.zeros_like(rasters_run)
    rasters_cue = np.zeros_like(rasters_run)

    for clu in tqdm(
        range(tot_clu),
        desc='generating spike array',
        disable=not verbose,
    ):
        for trial in range(1, tot_trial+1):
            before = (
                spike_time_file[time_bef[trial, clu]][0]
                if not isinstance(spike_time_file[time_bef[trial, clu]][0], np.uint64)
                else []
                )
            after = (
                spike_time_file[time_aft[trial, clu]][0]
                if not isinstance(spike_time_file[time_aft[trial, clu]][0], np.uint64)
                else []
                )
            spikes = np.asarray(np.concatenate((before, after)), dtype=np.int32)
            spikes = spikes + BEF*SAMP_FREQ
            spikes = spikes[spikes < MAX_LENGTH]
            rasters[clu, trial-1, spikes] = 1

        for trial in range(tot_trials):
            run_x0 = max(run_onsets[trial] - BEF*SAMP_FREQ, 0)
            run_x1 = min(run_onsets[trial] + AFT*SAMP_FREQ, max_time)
            rasters_run[clu, trial, : run_x1-run_x0] = spike_map[clu, run_x0 : run_x1]

            rew_x0 = max(rew_onsets[trial] - BEF*SAMP_FREQ, 0)
            rew_x1 = min(rew_onsets[trial] + AFT*SAMP_FREQ, max_time)
            rasters_rew[clu, trial, : rew_x1-rew_x0] = spike_map[clu, rew_x0 : rew_x1]

            cue_x0 = max(cue_onsets[trial] - BEF*SAMP_FREQ, 0)
            cue_x1 = min(cue_onsets[trial] + AFT*SAMP_FREQ, max_time)
            rasters_cue[clu, trial, : cue_x1-cue_x0] = spike_map[clu, cue_x0 : cue_x1]

    spike_file_handle.close()
    del spike_time_file, time_bef, time_aft

    t0 = time()
    if GPU_AVAILABLE:
        # Upload one raster stack at a time to keep GPU peak memory bounded.
        trains = convolve_raster(
            rasters,
            GAUS_SPIKE_ALL_TRAINS,
        )
        trains_run = convolve_raster(
            rasters_run,
            GAUS_SPIKE,
        )
        trains_rew = convolve_raster(
            rasters_rew,
            GAUS_SPIKE,
        )
        trains_cue = convolve_raster(
            rasters_cue,
            GAUS_SPIKE,
        )
        if verbose:
            print(
                'convolution done in '
                f'{str(timedelta(seconds=int(time() - t0)))} s')
    else:
        # CPU convolution using SciPy's FFT-based convolution for better performance
        trains = fftconvolve(
            rasters, GAUS_SPIKE_ALL_TRAINS[None, None, :],
            mode='same'
            ) * SAMP_FREQ
        trains_run = fftconvolve(
            rasters_run, GAUS_SPIKE[None, None, :],
            mode='same'
            ) * SAMP_FREQ
        trains_rew = fftconvolve(
            rasters_rew, GAUS_SPIKE[None, None, :],
            mode='same'
            ) * SAMP_FREQ
        trains_cue = fftconvolve(
            rasters_cue, GAUS_SPIKE[None, None, :],
            mode='same'
            ) * SAMP_FREQ

        if verbose:
            print(
                'convolution on CPU done in '
                f'{str(timedelta(seconds=int(time() - t0)))} s')

    for clu in range(tot_clu):
        cluname = f'{recname} clu{clu+2} {int(shank[clu])} {int(localClu[clu])}'

        all_trains[cluname] = trains[clu]
        all_rasters[cluname] = rasters[clu]

        all_trains_run[cluname] = trains_run[clu]
        all_rasters_run[cluname] = rasters_run[clu]

        all_trains_rew[cluname] = trains_rew[clu]
        all_rasters_rew[cluname] = rasters_rew[clu]

        all_trains_cue[cluname] = trains_cue[clu]
        all_rasters_cue[cluname] = rasters_cue[clu]

    # smooth entire spike map for all clusters
    if verbose:
        print('smoothing full spike maps...')
    if GPU_AVAILABLE:
        spike_map_gpu = cp.asarray(spike_map)
        smoothed_spike_map = (
            cpss.fftconvolve(
                spike_map_gpu,
                GAUS_SPIKE[None, :],
                mode='same'
            ) * SAMP_FREQ
        ).get()
        del spike_map_gpu
        release_gpu_memory()
    else:
        smoothed_spike_map = fftconvolve(
            spike_map,
            GAUS_SPIKE_CPU[None, :],
            mode='same'
        ) * SAMP_FREQ

    if verbose:
        print('done; saving...')
    sess_stem.mkdir(parents=True, exist_ok=True)
    if recname not in PROFILE_INPUT_SKIP_RECNAMES:
        np.save(
            sess_stem / f'{recname}_all_trains.npy',
            all_trains
            )
        np.save(
            sess_stem / f'{recname}_all_rasters.npy',
            all_rasters
            )
    np.save(
        sess_stem / f'{recname}_all_trains_run.npy',
        all_trains_run
        )
    np.save(
        sess_stem / f'{recname}_all_rasters_run.npy',
        all_rasters_run
        )
    np.save(
        sess_stem / f'{recname}_all_trains_rew.npy',
        all_trains_rew
        )
    np.save(
        sess_stem / f'{recname}_all_rasters_rew.npy',
        all_rasters_rew
        )
    np.save(
        sess_stem / f'{recname}_all_trains_cue.npy',
        all_trains_cue
        )
    np.save(
        sess_stem / f'{recname}_all_rasters_cue.npy',
        all_rasters_cue
        )
    np.save(
        sess_stem / f'{recname}_smoothed_spike_map.npy',
        smoothed_spike_map
        )
    if verbose:
        print_files_saved([
            ('session folder', sess_stem),
        ])
        print(f'elapsed = {str(timedelta(seconds=int(time() - t0)))}')
    if GPU_AVAILABLE:
        release_gpu_memory()
    gc.collect()
    return sess_stem

def main(argv=None):
    parser = argparse.ArgumentParser(
        description='extract HPC spike trains aligned to behavioural landmarks.'
    )
    parser.add_argument(
        '--recording',
        dest='recording_filter',
        help='optional substring filter for recname or recording path',
    )
    parser.add_argument(
        '--start-at',
        dest='start_at',
        help='optional recname/path substring; process this recording and all later recordings in rec_list order',
    )
    parser.add_argument(
        '--verbose',
        action='store_true',
        help='print per-recording details',
    )
    args = parser.parse_args(argv)

    session_paths = paths
    if args.start_at:
        start_idx = next(
            (
                idx for idx, rec_path in enumerate(session_paths)
                if args.start_at in rec_path
                or args.start_at in Path(rec_path).name
            ),
            None,
        )
        if start_idx is None:
            raise ValueError(f'no recordings matched start-at filter: {args.start_at}')
        session_paths = session_paths[start_idx:]

    if args.recording_filter:
        session_paths = [
            rec_path for rec_path in session_paths
            if args.recording_filter in rec_path
            or args.recording_filter in Path(rec_path).name
        ]
        if not session_paths:
            raise ValueError(f'no recordings matched filter: {args.recording_filter}')

    if args.verbose:
        for path in session_paths:
            process_session(path, verbose=True)
    else:
        session_iter = tqdm(session_paths, desc='recordings')
        for path in session_iter:
            session_iter.set_postfix_str(Path(path).name)
            process_session(path)

if __name__ == '__main__':
    main()
