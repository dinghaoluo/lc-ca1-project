# -*- coding: utf-8 -*-
'''
Created on Mon 10 July 10:02:32 2023
Modified on 29 Aug 2025

pool all cells from all recording sessions
modified 11 Dec 2024 to process with all trials (not skipping trial 0) and
    added GPU support
    - memory leakage problems on GPU, 20 Dec 2024
    - process each recording inside a function so its GPU references are
        destroyed afterwards, 26 Dec 2024

modified to be used on Raphi's recordings

@author: Dinghao Luo

'''

#%% imports
import argparse
from pathlib import Path
import sys

import h5py
import mat73
import scipy.io as sio
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

import rec_list
paths = rec_list.pathHPC_Raphi

from common_functions import gaussian_kernel_unity, get_GPU_availability
from console_formatting import print_files_saved, print_session
import project_paths as pp


#%% GPU acceleration
cp, GPU_AVAILABLE, _ = get_GPU_availability()

if GPU_AVAILABLE:
    xp = cp
    import cupyx.scipy.signal as cpss  # for GPU
    import numpy as np
    mempool = cp.get_default_memory_pool()
    pinned_mempool = cp.get_default_pinned_memory_pool()
else:
    import numpy as np
    xp = np


#%% parameters
SAMP_FREQ   = 1250  # Hz
SIGMA_SPIKE = int(SAMP_FREQ * 0.05)  # 50 ms
GAUS_SPIKE  = gaussian_kernel_unity(SIGMA_SPIKE, GPU_AVAILABLE)
SIGMA_SPIKE_ALL_TRAINS = int(SAMP_FREQ * 0.03)  # 30 ms
GAUS_SPIKE_ALL_TRAINS = gaussian_kernel_unity(
    SIGMA_SPIKE_ALL_TRAINS,
    GPU_AVAILABLE
    )
GAUS_SPIKE_ALL_TRAINS_CPU = gaussian_kernel_unity(
    SIGMA_SPIKE_ALL_TRAINS,
    False
    )

MAX_LENGTH = 12500  # samples

BEF = 3  # seconds before
AFT = 7  # seconds after

HPC_ALL_SESSIONS_RAPHI_STEM = pp.HPC_EPHYS_STEM / 'all_sessions_raphi'


#%% main
def clear_gpu_memory():
    gc.collect()
    cp.cuda.runtime.deviceSynchronize()
    cp.fft.config.get_plan_cache().clear()
    mempool.free_all_blocks()
    pinned_mempool.free_all_blocks()
    gc.collect()

def process_session(path, verbose=False, all_trains_only=False):
    recname = Path(path).name

    rec_stem = Path(path)

    # Raphi's data require maze_sess to access
    # aligned behavioural landmarks
    aligned_run_path = None
    for maze_sess_idx in range(10):
        candidate_run_path = rec_stem / f'{recname}_DataStructure_mazeSection1_TrialType1_alignRun_msess{maze_sess_idx}.mat'
        if candidate_run_path.exists():
            aligned_run_path = candidate_run_path
            break
    if aligned_run_path is None:
        raise FileNotFoundError(f'no aligned run file found for {recname}')
    aligned_cue_path = rec_stem / f'{recname}_DataStructure_mazeSection1_TrialType1_alignCue_msess{maze_sess_idx}.mat'
    aligned_rew_path = rec_stem / f'{recname}_DataStructure_mazeSection1_TrialType1_alignRew_msess{maze_sess_idx}.mat'

    # spike file paths
    clu_paths = [rec_stem / f'{recname}.clu.{probe}' for probe in range(1,7)]
    res_paths = [rec_stem / f'{recname}.res.{probe}' for probe in range(1,7)]

    # get cluname
    filename = Path(path) / f'{recname}_BehavElectrDataLFP.mat'
    spikefile_name = rec_stem / f'{recname}_DataStructure_mazeSection1_TrialType1_alignedSpikesPerNPerT_msess{maze_sess_idx}_Run0.mat'

    sess_stem = HPC_ALL_SESSIONS_RAPHI_STEM / recname
    sess_stem.mkdir(parents=True, exist_ok=True)
    all_trains_path = sess_stem / f'{recname}_all_trains.npy'
    all_rasters_path = sess_stem / f'{recname}_all_rasters.npy'
    if verbose:
        print_session(recname)

    BehavLFP = mat73.loadmat(filename)
    Clu = BehavLFP['Clu']
    shank = Clu['shank']
    localClu = Clu['localClu']

    with h5py.File(spikefile_name, 'r') as spike_file_handle:
        spike_time_file = spike_file_handle['trialsRunSpikes']
        time_bef = spike_time_file['TimeBef']
        time_aft = spike_time_file['Time']
        tot_clu = time_aft.shape[1]
        tot_trial = time_aft.shape[0]-1  # trial 1 is empty
        all_trains_scale = SAMP_FREQ if int(recname[1:4]) <= 40 else 1

        if GPU_AVAILABLE and int(recname[1:4]) <= 40:
            rasters_gpu = xp.zeros(
                (tot_clu, tot_trial, MAX_LENGTH),
                dtype=xp.uint16
                )
            for clu in tqdm(
                range(tot_clu),
                desc='generating all-trial spike array (GPU)',
                disable=not verbose,
            ):
                for trial in range(1, tot_trial+1):
                    before = spike_time_file[time_bef[trial, clu]][0]
                    after = spike_time_file[time_aft[trial, clu]][0]
                    if isinstance(before, np.uint64):
                        before = []
                    if isinstance(after, np.uint64):
                        after = []
                    spikes = (
                        cp.asarray(
                            np.asarray(
                                np.concatenate((before, after)),
                                dtype=np.int32,
                                ),
                            dtype=cp.int32
                            )
                        + BEF * SAMP_FREQ
                        )
                    spikes = spikes[spikes < MAX_LENGTH]
                    rasters_gpu[clu, trial-1, spikes] = 1

            trains = (
                cpss.fftconvolve(
                    rasters_gpu,
                    GAUS_SPIKE_ALL_TRAINS[None, None, :],
                    mode='same'
                    )
                * all_trains_scale
                ).get()
            rasters = rasters_gpu.get()
            del rasters_gpu
            clear_gpu_memory()
        else:
            rasters = np.zeros(
                (tot_clu, tot_trial, MAX_LENGTH),
                dtype=np.uint16
                )
            for clu in tqdm(
                range(tot_clu),
                desc='generating all-trial spike array (CPU)',
                disable=not verbose,
            ):
                for trial in range(1, tot_trial+1):
                    before = spike_time_file[time_bef[trial, clu]][0]
                    after = spike_time_file[time_aft[trial, clu]][0]
                    if isinstance(before, np.uint64):
                        before = []
                    if isinstance(after, np.uint64):
                        after = []
                    spikes = (
                        np.asarray(
                            np.concatenate((before, after)),
                            dtype=np.int32,
                            )
                        + BEF * SAMP_FREQ
                        )
                    spikes = spikes[spikes < MAX_LENGTH]
                    rasters[clu, trial-1, spikes] = 1

            trains = (
                fftconvolve(
                    rasters,
                    GAUS_SPIKE_ALL_TRAINS_CPU[None, None, :],
                    mode='same'
                    )
                * all_trains_scale
                )

    all_trains = {}
    all_rasters = {}
    for clu in range(tot_clu):
        cluname = f'{recname} clu{clu+2} {int(shank[clu])} {int(localClu[clu])}'
        all_trains[cluname] = trains[clu]
        all_rasters[cluname] = rasters[clu]

    np.save(
        all_trains_path,
        all_trains
        )
    np.save(
        all_rasters_path,
        all_rasters
        )
    del trains, rasters, all_trains, all_rasters
    gc.collect()

    if all_trains_only:
        if verbose:
            print_files_saved([
                ('session folder', sess_stem),
            ])
        return sess_stem

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
    all_rasters_run = {}
    all_rasters_cue = {}
    all_rasters_rew = {}

    if GPU_AVAILABLE:
        rasters_run_gpu = xp.zeros((tot_clu, tot_trial, MAX_LENGTH), dtype=xp.uint16)
        rasters_rew_gpu = xp.zeros_like(rasters_run_gpu)
        rasters_cue_gpu = xp.zeros_like(rasters_run_gpu)

        for clu in tqdm(
            range(tot_clu),
            desc='generating spike array (GPU)',
            disable=not verbose,
        ):
            for trial in range(tot_trials):
                run_x0 = max(run_onsets[trial] - BEF*SAMP_FREQ, 0)
                run_x1 = min(run_onsets[trial] + AFT*SAMP_FREQ, max_time)
                rasters_run_gpu[clu, trial, : run_x1-run_x0] = cp.asarray(
                    spike_map[clu, run_x0 : run_x1])

                rew_x0 = max(rew_onsets[trial] - BEF*SAMP_FREQ, 0)
                rew_x1 = min(rew_onsets[trial] + AFT*SAMP_FREQ, max_time)
                rasters_rew_gpu[clu, trial, : rew_x1-rew_x0] = cp.asarray(
                    spike_map[clu, rew_x0 : rew_x1])

                cue_x0 = max(cue_onsets[trial] - BEF*SAMP_FREQ, 0)
                cue_x1 = min(cue_onsets[trial] + AFT*SAMP_FREQ, max_time)
                rasters_cue_gpu[clu, trial, : cue_x1-cue_x0] = cp.asarray(
                    spike_map[clu, cue_x0 : cue_x1])

    else:
        rasters_run = xp.zeros((tot_clu, tot_trial, MAX_LENGTH), dtype=xp.uint16)
        rasters_rew = xp.zeros_like(rasters_run)
        rasters_cue = xp.zeros_like(rasters_run)

        for clu in tqdm(
            range(tot_clu),
            desc='generating spike array (CPU)',
            disable=not verbose,
        ):
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

    t0 = time()
    if GPU_AVAILABLE:
        # GPU-accelerated convolution using CuPy
        trains_run = (cpss.fftconvolve(
            rasters_run_gpu, GAUS_SPIKE[None, None, :],
            mode='same'
            ) * SAMP_FREQ).get()
        trains_rew = (cpss.fftconvolve(
            rasters_rew_gpu, GAUS_SPIKE[None, None, :],
            mode='same'
            ) * SAMP_FREQ).get()
        trains_cue = (cpss.fftconvolve(
            rasters_cue_gpu, GAUS_SPIKE[None, None, :],
            mode='same'
            ) * SAMP_FREQ).get()

        rasters_run = rasters_run_gpu.get()
        rasters_rew = rasters_rew_gpu.get()
        rasters_cue = rasters_cue_gpu.get()
        del rasters_run_gpu, rasters_rew_gpu, rasters_cue_gpu
        clear_gpu_memory()

        if verbose:
            print(
                'convolution on GPU done in '
                f'{str(timedelta(seconds=int(time() - t0)))} s')
    else:
        # CPU convolution using SciPy's FFT-based convolution for better performance
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

        all_trains_run[cluname] = trains_run[clu]
        all_rasters_run[cluname] = rasters_run[clu]

        all_trains_rew[cluname] = trains_rew[clu]
        all_rasters_rew[cluname] = rasters_rew[clu]

        all_trains_cue[cluname] = trains_cue[clu]
        all_rasters_cue[cluname] = rasters_cue[clu]

    if verbose:
        print('done; saving...')
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
    if verbose:
        print_files_saved([
            ('session folder', sess_stem),
        ])
        print(f'elapsed = {str(timedelta(seconds=int(time() - t0)))}')
    return sess_stem

def main(argv=None):
    parser = argparse.ArgumentParser(
        description='extract Raphi/passive HPC spike trains aligned to behavioural landmarks.'
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
        '--all-trains-only',
        action='store_true',
        help='only write all_trains/all_rasters outputs',
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
        session_iter = session_paths
    else:
        session_iter = tqdm(session_paths, desc='recordings')

    for path in session_iter:
        if not args.verbose:
            session_iter.set_postfix_str(Path(path).name)
        if GPU_AVAILABLE:
            clear_gpu_memory()
        process_session(
            path,
            verbose=args.verbose,
            all_trains_only=args.all_trains_only,
        )
        if GPU_AVAILABLE:
            clear_gpu_memory()

if __name__ == '__main__':
    main()
