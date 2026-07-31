# -*- coding: utf-8 -*-
'''
Created on Thu Apr 24 15:08:04 2025
Originally named plot_raw traces_axon_gcamp.py

plot the raw trace aligned with speeds for LC-HPC axon-GCaMP data,
with distance resetting at trial start (180 cm per trial max).

@author: Dinghao Luo
'''

#%% imports
import sys
import gc
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pickle

repo_root = Path(__file__).resolve().parents[2]
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))

from common_functions import mpl_formatting, smooth_convolve
from console_formatting import print_session
import project_paths as pp
mpl_formatting()

import rec_list
paths = rec_list.pathLCHPCGCaMP


#%% parameters
ROI_size_threshold = 500  # pixels
smoothing_sigma = 100  # in ms


#%% main
for path in paths[48:]:
    recname = Path(path).name
    print_session(recname)

    # load fluorescence data
    plane_stem = pp.resolve_suite2p_session_stem(path) / 'plane0'
    F = np.load(plane_stem / 'F.npy')

    valid_ROIs_coord_dict = np.load(
        pp.LCHPC_AXON_STEM
        / 'all_sessions'
        / recname
        / 'processed_data'
        / 'valid_ROIs_coord_dict.npy',
        allow_pickle=True
    ).item()
    ROIs = [key for key, value in valid_ROIs_coord_dict.items()
            if len(value[0]) > ROI_size_threshold]
    ROI_idx = [int(roi.split(' ')[-1]) for roi in ROIs]

    # load behaviour data
    with open(pp.behaviour_experiment_stem('LCHPCGCaMP') / f'{recname}.pkl', 'rb') as f:
        beh = pickle.load(f)

    frame_times = beh['frame_times']
    first_frame_time = frame_times[0]
    last_frame_time = frame_times[-1]

    times, speeds = zip(
        *[(t, s) for t, s
          in zip(beh['upsampled_timestamps_ms'], beh['upsampled_speed_cm_s'])
          if first_frame_time < t < last_frame_time]
        )

    # distance interpolation with resetting at trial starts
    run_onsets = beh['run_onsets']
    run_onsets = [onset for onset in run_onsets
                  if first_frame_time < onset < last_frame_time]

    # generate distance array following run_onsets
    distances = np.zeros_like(speeds)
    run_onsets_idx = [np.searchsorted(times, onset) for onset in run_onsets]

    # start generating accumulated distances
    trial_idx = 0
    accumulating = False
    distance = 0
    for idx in range(len(times)):
        if trial_idx < len(run_onsets_idx) and idx == run_onsets_idx[trial_idx]:
            distance = 0
            accumulating = True
            trial_idx += 1

        if accumulating:
            distance += speeds[idx] / 1000  # speeds in cm/s but sampled at 1000 Hz
            if distance >= 180:  # 180 cm
                distance = 0
                accumulating = False

        distances[idx] = distance

    speeds_smoothed = smooth_convolve(speeds, 300)

    # output folder
    save_dir = pp.LCHPC_AXON_FIGURES_STEM / 'all_sessions' / recname / 'ROI_distance_dFF_speed'
    save_dir.mkdir(parents=True, exist_ok=True)

    ## plotting
    window_size = 100 * 1000  # 100 s
    n_windows = (len(times) + window_size - 1) // window_size
    for roi in ROI_idx:
        F_ROI = F[roi, :]
        n_frames = min(len(frame_times), len(F_ROI))
        F_ROI = F_ROI[:n_frames]
        F_ROI_interp = smooth_convolve(
            np.interp(times, frame_times[:n_frames], F_ROI),
            sigma=smoothing_sigma
        )

        for fig_idx in range(0, n_windows, 4):  # 4 windows per figure
            # four windows, each with three panels and one blank row
            fig, axs = plt.subplots(16, 1, figsize=(12, 28), sharex=False)
            fig.subplots_adjust(hspace=0.4)  # smaller hspace globally

            for i in range(4):  # 4 windows per figure
                window_idx = fig_idx + i
                if window_idx >= n_windows:
                    break  # no more windows

                start_idx = window_idx * window_size
                end_idx = min((window_idx + 1) * window_size, len(times))

                times_plot = np.array(times[start_idx:end_idx]) / 1000  # ms to s

                base = i * 4  # shift every window block by 4 rows

                # distance
                axs[base + 0].plot(times_plot - times_plot[0],
                                   distances[start_idx:end_idx] / 100,
                                   color='black')
                axs[base + 0].set_ylabel('distance (m)')
                axs[base + 0].set_title(f'ROI {roi} | Window {window_idx}')

                # fluorescence
                axs[base + 1].plot(times_plot - times_plot[0],
                                   F_ROI_interp[start_idx:end_idx],
                                   color='green')
                axs[base + 1].set_ylabel('dF/F')

                # speed
                axs[base + 2].plot(times_plot - times_plot[0],
                                   speeds_smoothed[start_idx:end_idx],
                                   color=(0.4, 0.5, 0.0))
                axs[base + 2].set_ylabel('speed (cm/s)')
                axs[base + 2].set_xlabel('time (s)')

                # blank panel (axs[base + 3])
                axs[base + 3].axis('off')  # turn off the blank panel

            # save figure
            fig.savefig(
                save_dir / f'roi_{roi}_{int(fig_idx/4)+1}.png',
                dpi=300,
                bbox_inches='tight'
                )
            plt.close(fig)

        # drop large per-ROI buffers
        del F_ROI, F_ROI_interp
        gc.collect()

    # beh dict seems to persist for whatever reason causing memory issues
    del beh, frame_times, times, speeds, distances, speeds_smoothed, run_onsets, run_onsets_idx
    gc.collect()
