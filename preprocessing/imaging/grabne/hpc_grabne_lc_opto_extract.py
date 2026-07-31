# -*- coding: utf-8 -*-
'''
Created on Mon Mar 31 13:02:49 2025
Modified on Tue 4 Nov 2025 to be used on GRABNE data

extract LC-opto GRABNE imaging data

@author: Dinghao Luo
'''

#%% imports
import argparse
import gc
from pathlib import Path
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import colormaps
from matplotlib.colors import TwoSlopeNorm
import tifffile
from tqdm import tqdm

repo_root = Path(__file__).resolve().parents[3]
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))

import behaviour_functions as bf
import imaging_pipeline_functions as ipf
from plotting_functions import plot_violin_with_scatter, add_scale_bar
from common_functions import mpl_formatting, get_GPU_availability
from console_formatting import print_session
import project_paths as pp
mpl_formatting()

import rec_list
paths = rec_list.pathGRABNELCOpto + rec_list.pathGRABNELCOptoDbhBlock
default_paths = paths

# GPU acceleration
cp, GPU_AVAILABLE, _ = get_GPU_availability()


#%% parameters
SAMP_FREQ = 30
SPATIAL_FILTER_CHUNK_SIZE = 128

# post-stim dispersion calculation
BIN_WIDTH = 0.1
SENSOR_LABEL = 'GRABNE'

# path stems
mice_exp_stem = pp.MICEEXP_ROOT
all_sess_stem = pp.HPC_GRABNE_LC_OPTO_STEM / 'all_sessions'
all_sess_fig_stem = pp.HPC_GRABNE_LC_OPTO_FIGURES_STEM / 'all_sessions'

def clear_memory(label=None):
    if label:
        print(f'clearing CPU/GPU memory after {label}...', flush=True)
    plt.close('all')
    gc.collect()
    if GPU_AVAILABLE:
        cp.cuda.Stream.null.synchronize()
        cp.get_default_memory_pool().free_all_blocks()
        cp.get_default_pinned_memory_pool().free_all_blocks()

def prepare_movie_channel(
        binpath, shape, recname, savepath, figurepath, make_plots,
        channel
        ):
    print(f'ch{channel}: loading movie into float32...', flush=True)
    mov_float = np.memmap(
        binpath, mode='r', dtype='int16', shape=shape
        ).astype(np.float32)

    print(f'ch{channel}: saving reference array...', flush=True)
    ref = ipf.plot_reference(
        mov_float, recname=recname, outpath=savepath, channel=channel,
        save_figure=make_plots, figure_outpath=figurepath
        )

    print(f'ch{channel}: computing whole-field raw trace...', flush=True)
    # raw traces replaced dF/F here because the dF/F baseline can cover the
    # PMT-off stimulation period, 5 Aug 2025
    raw_trace = np.sum(mov_float, axis=(1,2))

    del mov_float
    clear_memory(f'ch{channel} movie loading')

    mov = np.memmap(binpath, mode='r', dtype='int16', shape=shape)
    return mov, ref, raw_trace

def save_release_map_figure(
        savepath, recname, ref, ref2, release_map, suffix, title,
        positive_floor=-.001):
    vmin = np.nanpercentile(release_map, 1)
    vmax = np.nanpercentile(release_map, 99)

    if vmin >= 0:
        vmin = positive_floor
    if vmax <= 0:
        vmax = abs(positive_floor)

    norm = TwoSlopeNorm(vcenter=0, vmin=vmin, vmax=vmax)

    fig, axs = plt.subplots(1, 3, figsize=(12, 4))

    axs[0].imshow(ref, cmap='gray', interpolation='none')
    axs[0].set_title('channel 1', fontsize=10)
    axs[0].axis('off')

    axs[1].imshow(ref2, cmap='gray', interpolation='none')
    axs[1].set_title('channel 2', fontsize=10)
    axs[1].axis('off')

    im = axs[2].imshow(release_map, cmap='RdBu_r', norm=norm, interpolation='none')
    axs[2].set_title(title, fontsize=10)
    axs[2].axis('off')

    cbar = fig.colorbar(im, ax=axs[2], shrink=0.8, fraction=0.046, pad=0.04)
    cbar.set_label('ΔF/F ratio', fontsize=10)
    cbar.set_ticks([vmin, 0, vmax])

    fig.tight_layout()

    for ext in ['.png', '.pdf']:
        fig.savefig(
            savepath / f'{recname}_{suffix}{ext}',
            dpi=300,
            bbox_inches='tight'
        )
    plt.close(fig)

def save_release_map_tiff(savepath, recname, release_map, suffix):
    vmin = np.nanpercentile(release_map, 1)
    vmax = np.nanpercentile(release_map, 99)

    if vmin >= 0:
        vmin = -.001
    if vmax <= 0:
        vmax = .001

    norm = TwoSlopeNorm(vcenter=0, vmin=vmin, vmax=vmax)
    cmap = colormaps['RdBu_r']
    release_map_rgba = cmap(norm(release_map))
    release_map_rgb = (release_map_rgba[..., :3] * 255).astype(np.uint8)

    tifffile.imwrite(savepath / f'{recname}_{suffix}.tiff',
                     release_map_rgb)


#%% main
def process_session(path, make_plots=False):
    recname = Path(path).name
    print_session(recname)

    plane_stem = pp.resolve_suite2p_session_stem(path) / 'plane0'
    sessname = recname.replace('i', '')

    binpath = plane_stem / 'data.bin'
    bin2path = plane_stem / 'data_chan2.bin'
    opspath = plane_stem / 'ops.npy'
    txtpath = mice_exp_stem / f'ANMD{recname[1:4]}' / f'{sessname}T.txt'

    savepath = all_sess_stem / f'{recname}'
    savepath.mkdir(parents=True, exist_ok=True)
    processed_datapath = savepath / 'processed_data'
    processed_datapath.mkdir(parents=True, exist_ok=True)
    figurepath = all_sess_fig_stem / f'{recname}'
    if make_plots:
        figurepath.mkdir(parents=True, exist_ok=True)

    txt = ipf.process_txt_nobeh(txtpath)

    # load data
    ops = np.load(opspath, allow_pickle=True).item()
    tot_frames = ops['nframes']
    shape = tot_frames, ops['Ly'], ops['Lx']

    if make_plots:
        print('loading movies, saving references, and generating reference figures...')
    else:
        print('loading movies and saving reference arrays...')
    mov, ref, raw_trace = prepare_movie_channel(
        binpath, shape, recname, savepath, figurepath, make_plots,
        channel=1
        )
    mov2, ref2, raw_trace2 = prepare_movie_channel(
        bin2path, shape, recname, savepath, figurepath, make_plots,
        channel=2
        )

    tot_frames = mov.shape[0]  # once loaded, update tot_frames to be the max frame number, 16 June 2025

    # dFF traces are ONLY used for plotting figure 1 now
    print('computing dFF traces...')
    trace_dFF = ipf.calculate_dFF(raw_trace, sigma=300, t_axis=0,
                                  GPU_AVAILABLE=GPU_AVAILABLE)
    trace2_dFF = ipf.calculate_dFF(raw_trace2, sigma=300, t_axis=0,
                                   GPU_AVAILABLE=GPU_AVAILABLE)

    # behaviour
    print('processing .txt file...')
    if txt is None:
        txt = ipf.process_txt_nobeh(txtpath)
    frame_times = txt['frame_times']
    pulse_info = ipf.extract_opto_pulse_metadata(txt, SAMP_FREQ)
    pulse_trains = pulse_info['pulse_trains']
    pulse_width_ON = pulse_info['pulse_width_on']
    pulse_width = pulse_info['pulse_width']
    pulse_number = pulse_info['pulse_number']
    taper_enabled = pulse_info['taper_enabled']
    taper_duration = pulse_info['taper_duration']
    duty_cycle = pulse_info['duty_cycle']
    tot_pulses = pulse_info['total_pulses']
    last_time = pulse_info['total_train_duration_ms']

    print(f'\npulse ON time: {pulse_width_ON}')
    print(f'pulse width: {pulse_width}')
    print(f'pulse number: {pulse_number}')

    if taper_enabled:
        print(f'taper duration: {taper_duration}')

    print(f'pulse count: {len(pulse_trains[0])}')
    print(f'total train duration: {last_time} ms')
    print(f'predicted PMT shut-off at {last_time + 100} ms')

    ## -- PARAMETER DEFINITIONS
    # now defined within the function scope, since we reassign them later in an if statement
    BEF = 2
    AFT = 10
    TAXIS = np.arange(-BEF*SAMP_FREQ, AFT*SAMP_FREQ) / SAMP_FREQ
    BASELINE_IDX = (TAXIS >= -1.0) & (TAXIS <= -0.15)
    ## -- END PARAMETER DEFINITIONS

    # determine time bin mask
    # NOTE: stim_start means the start of the STIM_IDX period, which comes AFTER the stim
    PMT_BUFFER_FRAMES = 8  # frames
    PMT_BUFFER = PMT_BUFFER_FRAMES / SAMP_FREQ
    last_time_s = last_time / 1_000  # convert to seconds
    if last_time_s >= 5:  # for long stimulation attempts to elicit NE release, 30 Oct 2025
        stim_start = last_time_s + PMT_BUFFER
        stim_end   = stim_start + 10  # extract 10 s of activity for NE mean
        # reset BEF and AFT, along with TAXIS
        BEF = 10  # s
        AFT = 30
        TAXIS = np.arange(-BEF*SAMP_FREQ, AFT*SAMP_FREQ) / SAMP_FREQ
        BASELINE_IDX = (TAXIS >= -10.0) & (TAXIS <= -0.15)
    else:
        stim_start = last_time_s + PMT_BUFFER
        stim_end   = stim_start + 1

    STIM_IDX = (TAXIS >= stim_start) & (TAXIS < stim_end)  # note that this is the mask for extracting stim_mean

    pulse_frames = [
        [ipf.find_nearest(p, frame_times) for p in train]
        for train in pulse_trains
    ]
    pulse_start_frames = [p[0] for p in pulse_frames]
    valid_pulse_start_frames = [
        p for p in pulse_start_frames
        if (p - BEF * SAMP_FREQ >= 0) and (p + AFT * SAMP_FREQ <= tot_frames)
        ]
    tot_valid_pulses = len(valid_pulse_start_frames)

    # post-stim dispersion calculation, 10 Sept 2025
    BIN_START = stim_start
    BIN_END   = stim_start + 4
    bin_edges = np.arange(BIN_START, BIN_END + BIN_WIDTH, BIN_WIDTH)
    n_bins = len(bin_edges) - 1

    # pulse processing
    print('extracting data...')

    # checks $FM against tot_frame
    if tot_frames < len(frame_times) - 3 or tot_frames > len(frame_times):
        raise ValueError(f'check $FM; frame count mismatch for {recname}')

    # filter for opto artefact periods
    pulse_period_frames = np.concatenate(
        [np.arange(
            max(0, pulse_train[0]-3),
            min(pulse_train[-1]+int(pulse_width)*SAMP_FREQ+PMT_BUFFER_FRAMES, tot_frames)  # +15 as a buffer, half a second
            )
        for pulse_train in pulse_frames]
        )

    # filtering
    raw_trace[pulse_period_frames]  = np.nan
    raw_trace2[pulse_period_frames] = np.nan

    trace_dFF[pulse_period_frames]  = np.nan
    trace2_dFF[pulse_period_frames] = np.nan

    # raw traces aligned
    raw_aligned  = np.zeros((tot_valid_pulses, (BEF+AFT)*SAMP_FREQ), dtype=np.float32)
    raw2_aligned = np.zeros((tot_valid_pulses, (BEF+AFT)*SAMP_FREQ), dtype=np.float32)
    for i, p in enumerate(valid_pulse_start_frames):
        start = p - BEF * SAMP_FREQ
        end   = p + AFT * SAMP_FREQ
        raw_aligned[i, :]  = raw_trace[start:end]
        raw2_aligned[i, :] = raw_trace2[start:end]

    # dFF traces aligned
    trace_dFF_aligned = np.zeros((tot_valid_pulses, (BEF+AFT)*SAMP_FREQ), dtype=np.float32)
    trace2_dFF_aligned = np.zeros((tot_valid_pulses, (BEF+AFT)*SAMP_FREQ), dtype=np.float32)
    for i, p in enumerate(valid_pulse_start_frames):
        start = p - BEF * SAMP_FREQ
        end   = p + AFT * SAMP_FREQ
        trace_dFF_aligned[i, :] = trace_dFF[start:end]
        trace2_dFF_aligned[i, :] = trace2_dFF[start:end]
    trace_dFF_aligned_mean = np.nanmean(trace_dFF_aligned, axis=0)
    trace2_dFF_aligned_mean = np.nanmean(trace2_dFF_aligned, axis=0)
    trace_dFF_aligned_mean[60:83] = np.nan
    trace2_dFF_aligned_mean[60:83] = np.nan

    # calculate ratios
    # per‐trial raw means
    baseline_raw = np.nanmean(raw_aligned[:,  BASELINE_IDX], axis=1)
    stim_raw = np.nanmean(raw_aligned[:,  STIM_IDX], axis=1)
    baseline2_raw = np.nanmean(raw2_aligned[:, BASELINE_IDX], axis=1)
    stim2_raw = np.nanmean(raw2_aligned[:, STIM_IDX], axis=1)

    # per‐trial ΔF/F exactly like pixel dFF (stim − base) / |base|
    dFF = (stim_raw - baseline_raw) / np.abs(baseline_raw)
    dFF2 = (stim2_raw - baseline2_raw) / np.abs(baseline2_raw)

    # dFF comp
    baseline_dFF = np.nanmean(trace_dFF_aligned[:,  BASELINE_IDX], axis=1)
    stim_dFF = np.nanmean(trace_dFF_aligned[:,  STIM_IDX], axis=1)
    baseline2_dFF = np.nanmean(trace2_dFF_aligned[:, BASELINE_IDX], axis=1)
    stim2_dFF = np.nanmean(trace2_dFF_aligned[:, STIM_IDX], axis=1)

    plot_payload = {
        'sensor_label': SENSOR_LABEL,
        'samp_freq': SAMP_FREQ,
        'before': BEF,
        'after': AFT,
        'taxis': TAXIS,
        'trace_dFF_aligned': trace_dFF_aligned,
        'trace2_dFF_aligned': trace2_dFF_aligned,
        'trace_dFF_aligned_mean': trace_dFF_aligned_mean,
        'trace2_dFF_aligned_mean': trace2_dFF_aligned_mean,
        'dFF': dFF,
        'dFF2': dFF2,
        'baseline_dFF': baseline_dFF,
        'stim_dFF': stim_dFF,
        'baseline2_dFF': baseline2_dFF,
        'stim2_dFF': stim2_dFF,
        'pulse_number': pulse_number,
        'pulse_width': pulse_width,
        'pulse_width_on': pulse_width_ON,
        'taper_enabled': taper_enabled,
        'taper_duration': taper_duration,
        'duty_cycle': duty_cycle,
        'total_pulse_trains': tot_pulses,
        }
    np.save(processed_datapath / f'{recname}_opto_plot_payload.npy',
            plot_payload)

    if make_plots:
        ymin = np.nanmin(trace_dFF_aligned.T)
        ymin2 = np.nanmin(trace2_dFF_aligned.T)

        fig, axs = plt.subplots(3,1,figsize=(3.5,5),
                                sharex=True)

        axs[0].plot(TAXIS, trace_dFF_aligned_mean, color='green', linewidth=2)
        axs[0].plot(TAXIS, trace_dFF_aligned.T, color='green', alpha=.05)
        axs[0].text(-2, ymin-.01, 'GRABNE', ha='left', va='center', color='green', fontsize=10)

        axs[1].plot(TAXIS, trace2_dFF_aligned_mean, color='darkred', linewidth=2)
        axs[1].plot(TAXIS, trace2_dFF_aligned.T, color='darkred', alpha=.05)
        axs[1].text(-2, ymin2-.01, 'red ctrl.', ha='left', va='center', color='darkred', fontsize=10)

        add_scale_bar(axs[0], x_start=8.5, y_start=ymin-.015, x_len=1, y_len=.02)
        add_scale_bar(axs[1], x_start=8.5, y_start=ymin2-.015, x_len=1, y_len=.02)
        axs[1].text(8.93, ymin2-.02, '1 s', ha='center', va='top', fontsize=8)
        axs[1].text(8.45, ymin2-.005, '2% ΔF/F', ha='right', va='center', rotation='vertical', fontsize=8)

        stim_trace = np.zeros_like(TAXIS)
        for p in range(pulse_number):
            stim_onset = BEF + p * (pulse_width)
            stim_offset = stim_onset + pulse_width_ON
            stim_onset_idx = int(stim_onset * SAMP_FREQ)
            stim_offset_idx = int(stim_offset * SAMP_FREQ)
            stim_trace[stim_onset_idx:stim_offset_idx] = 1

        if taper_enabled and taper_duration > 0:
            taper_start_idx = stim_offset_idx
            taper_len = int((taper_duration / 1_000_000) * SAMP_FREQ)
            taper = np.linspace(1,0, taper_len, endpoint=True)
            stim_trace[taper_start_idx : taper_start_idx+taper_len] = taper[:len(stim_trace) - taper_start_idx]

        axs[2].plot(TAXIS, stim_trace, color='k', linewidth=1)

        axs[2].set_ylim(0, 2)
        axs[2].set_xticks([])
        axs[2].set_yticks([])

        for ax in axs:
            ax.set_xticks([])
            ax.set_yticks([])
            for s in ['top', 'right', 'left', 'bottom']:
                ax.spines[s].set_visible(False)

        if taper_enabled:
            title_str = (
                f'{recname}\n'
                f'duty cycle = {duty_cycle}\n'
                f'pulse width = {pulse_width} s\n'
                f'pulse(s) per pulse train = {pulse_number}\n'
                f'total pulse trains = {tot_pulses}\n'
                f'taper duration = {taper_duration}'
                )
        else:
            title_str = (
                f'{recname}\n'
                f'duty cycle = {duty_cycle}\n'
                f'pulse width = {pulse_width} s\n'
                f'pulse(s) per pulse train = {pulse_number}\n'
                f'total pulse trains = {tot_pulses}\n'
                )

        fig.suptitle(title_str)
        fig.tight_layout()

        for ext in ['.png', '.pdf']:
            fig.savefig(
                figurepath / f'{recname}_aligned_stim{ext}',
                dpi=300,
                bbox_inches='tight'
                )
        plt.close(fig)

        valid_mask = (~np.isnan(dFF)) & (~np.isnan(dFF2))
        dFF_valid = dFF[valid_mask]
        dFF2_valid = dFF2[valid_mask]
        plot_violin_with_scatter(
            dFF2_valid, dFF_valid,
            'darkred', 'green',
            xticklabels=['ref.', 'GRABNE'],
            ylabel='ΔF/F',
            title=recname,
            save=True,
            savepath=figurepath / f'{recname}_dFF_dFF2_violinplot',
            show=False,
            close=True,
            )

        baseline_dFF, stim_dFF = map(list, zip(*[
            (b, s) for b, s in zip(baseline_dFF, stim_dFF)
            if not np.isnan(b) and not np.isnan(s)
        ]))
        baseline2_dFF, stim2_dFF = map(list, zip(*[
            (b, s) for b, s in zip(baseline2_dFF, stim2_dFF)
            if not np.isnan(b) and not np.isnan(s)
        ]))
        plot_violin_with_scatter(
            baseline_dFF, stim_dFF,
            '#8CA082', 'green',
            xticklabels=['baseline', 'stim.'],
            ylabel='ΔF/F',
            title=recname,
            save=True,
            savepath=figurepath / f'{recname}_baseline_stim_violinplot',
            show=False,
            close=True,
            )
        plot_violin_with_scatter(
            baseline2_dFF, stim2_dFF,
            '#8C6464', 'darkred',
            xticklabels=['baseline\nch2', 'stim.\nch2'],
            ylabel='ΔF/F',
            title=recname,
            save=True,
            savepath=figurepath / f'{recname}_baseline_stim_ch2_violinplot',
            show=False,
            close=True,
            )

    ## pixel-wise extraction
    # spatial smoothing
    clear_memory(f'{recname} trace extraction')
    print('performing trial-window spatial filtering and pixel extraction...')

    # compute dF/F per pixel (stim. / baseline for raw F), 24 June 2025
    pixel_dFF = np.zeros((shape[1], shape[2], len(valid_pulse_start_frames)),
                         dtype=np.float32)
    pixel_dFF2 = np.zeros_like(pixel_dFF)

    # we still want F aligned
    pixel_F_aligned = np.zeros((len(valid_pulse_start_frames),
                                ((BEF+AFT) * SAMP_FREQ),
                                shape[1],
                                shape[2]),
                               dtype=np.float32)
    pixel_F2_aligned = np.zeros_like(pixel_F_aligned)
    all_bins_ch1 = []
    all_bins_ch2 = []

    pulse_iter = tqdm(valid_pulse_start_frames,
                      desc='pixel-wise stim extraction')
    for i, p in enumerate(pulse_iter):
        start = p - BEF * SAMP_FREQ
        end = p + AFT * SAMP_FREQ
        temp_F, temp_F2 = [
            ipf.spatial_gaussian_filter(
                np.array(movie[start:end], dtype=np.float32, copy=True),
                sigma_spatial=1,
                GPU_AVAILABLE=GPU_AVAILABLE,
                CHUNK=True,
                chunk_size=SPATIAL_FILTER_CHUNK_SIZE,
                inplace=True,
                show_progress=False,
                )
            for movie in (mov, mov2)
            ]

        # save F aligned too
        pixel_F_aligned[i, :, :, :] = temp_F
        pixel_F2_aligned[i, :, :, :] = temp_F2

        stim_mean = np.mean(temp_F[STIM_IDX, :, :], axis=0)
        baseline_mean = np.mean(temp_F[BASELINE_IDX, :, :], axis=0)
        dFF = (stim_mean - baseline_mean) / np.abs(baseline_mean)
        dFF[np.abs(dFF) > 10] = np.nan  # hard cap
        pixel_dFF[:, :, i] = dFF

        stim_mean2 = np.mean(temp_F2[STIM_IDX, :, :], axis=0)
        baseline_mean2 = np.mean(temp_F2[BASELINE_IDX, :, :], axis=0)
        dFF2 = (stim_mean2 - baseline_mean2) / np.abs(baseline_mean2)
        dFF2[np.abs(dFF2) > 10] = np.nan
        pixel_dFF2[:, :, i] = dFF2

        trial_bins_ch1 = np.zeros((shape[1], shape[2], n_bins),
                                  dtype=np.float32)
        trial_bins_ch2 = np.zeros_like(trial_bins_ch1)

        for b in range(n_bins):
            bin_mask = (TAXIS >= bin_edges[b]) & (TAXIS < bin_edges[b+1])

            stim_mean = np.mean(temp_F[bin_mask, :, :], axis=0)
            baseline_mean = np.mean(temp_F[BASELINE_IDX, :, :], axis=0)
            dFF_bin = (stim_mean - baseline_mean) / np.abs(baseline_mean)
            dFF_bin[np.abs(dFF_bin) > 10] = np.nan
            trial_bins_ch1[:, :, b] = dFF_bin

            stim_mean2 = np.mean(temp_F2[bin_mask, :, :], axis=0)
            baseline_mean2 = np.mean(temp_F2[BASELINE_IDX, :, :], axis=0)
            dFF_bin2 = (stim_mean2 - baseline_mean2) / np.abs(baseline_mean2)
            dFF_bin2[np.abs(dFF_bin2) > 10] = np.nan
            trial_bins_ch2[:, :, b] = dFF_bin2

        all_bins_ch1.append(trial_bins_ch1)
        all_bins_ch2.append(trial_bins_ch2)

        del temp_F, temp_F2, trial_bins_ch1, trial_bins_ch2
        clear_memory()

    np.save(processed_datapath / f'{recname}_pixel_dFF_stim.npy',
            pixel_dFF)
    np.save(processed_datapath / f'{recname}_pixel_dFF_ch2_stim.npy',
            pixel_dFF2)

    # save F aligned
    np.save(processed_datapath / f'{recname}_pixel_F_aligned.npy',
            pixel_F_aligned)
    np.save(processed_datapath / f'{recname}_pixel_F2_aligned.npy',
            pixel_F2_aligned)

    # generate mean dFF release map as proxy for t-map, 24 June 2025
    print('computing mean release map...')
    release_map = np.nanmean(pixel_dFF, axis=2)  # shape: (y, x)
    release_map2 = np.nanmean(pixel_dFF2, axis=2)

    # save map matrices
    np.save(processed_datapath / f'{recname}_release_map.npy', release_map)
    np.save(processed_datapath / f'{recname}_release_map_ch2.npy', release_map2)

    if make_plots:
        save_release_map_figure(
            figurepath, recname, ref, ref2, release_map,
            'release_map', 'CH1: stim / baseline (mean)'
            )
        save_release_map_tiff(figurepath, recname, release_map, 'release_map')
        save_release_map_figure(
            figurepath, recname, ref, ref2, release_map2,
            'release_map_ch2', 'CH2: stim / baseline (mean)'
            )
        save_release_map_tiff(figurepath, recname, release_map2, 'release_map_ch2')

    # dispersion rate calculation, 10 Sept 2025
    print('calculating binned ratios (for dispersion rate analysis)...')

    # average across trials → final shape (y, x, n_bins)
    pixel_dFF_bins  = np.nanmean(np.stack(all_bins_ch1, axis=-1), axis=-1)
    pixel_dFF2_bins = np.nanmean(np.stack(all_bins_ch2, axis=-1), axis=-1)

    # save arrays
    np.save(processed_datapath / f'{recname}_pixel_dFF_bins.npy', pixel_dFF_bins)
    np.save(processed_datapath / f'{recname}_pixel_dFF_ch2_bins.npy', pixel_dFF2_bins)

    ## compute dF/F per pixel IF BEHAVIOUR (run-onset / baseline), 24 June 2025
    if txt['behaviour']:
        print('behaviour session; compiling run-onset dFF dict...')
        txt = bf.process_behavioural_data_imaging(txtpath)
        run_onsets = txt['run_onset_frames']
        stim_conds = [t[15] for t in txt['trial_statements']]
        stim_idx = [trial for trial, cond in enumerate(stim_conds)
                    if cond!='0']
        stim_idx_et = [trial + 1 for trial in stim_idx if trial + 1 < len(run_onsets)]

        # new axes for run
        BASELINE_IDX_RUN = (TAXIS >= -1.0) & (TAXIS <= -0.15)
        RUN_IDX = (TAXIS >= 0.15) & (TAXIS <= 1.0)

        run_onsets = [f for trial, f in enumerate(run_onsets)
                      if not np.isnan(f)
                      and trial not in stim_idx and trial not in stim_idx_et
                      and f > BEF*SAMP_FREQ
                      and f < tot_frames - AFT*SAMP_FREQ]

        pixel_dFF_run = np.zeros((shape[1], shape[2], len(run_onsets)),
                                 dtype=np.float32)
        pixel_dFF_run2 = np.zeros_like(pixel_dFF_run)
        run_iter = tqdm(run_onsets, desc='pixel-wise run extraction')
        for i, f in enumerate(run_iter):
            start = f - BEF * SAMP_FREQ
            end = f + AFT * SAMP_FREQ
            temp_F, temp_F2 = [
                ipf.spatial_gaussian_filter(
                    np.array(movie[start:end], dtype=np.float32, copy=True),
                    sigma_spatial=1,
                    GPU_AVAILABLE=GPU_AVAILABLE,
                    CHUNK=True,
                    chunk_size=SPATIAL_FILTER_CHUNK_SIZE,
                    inplace=True,
                    show_progress=False,
                    )
                for movie in (mov, mov2)
                ]

            run_mean = np.mean(temp_F[RUN_IDX, :, :], axis=0)
            prerun_mean = np.mean(temp_F[BASELINE_IDX_RUN, :, :], axis=0)
            dFF_run = (run_mean - prerun_mean) / np.abs(prerun_mean)
            dFF_run[np.abs(dFF_run) > 10] = np.nan
            pixel_dFF_run[:, :, i] = dFF_run

            run_mean2 = np.mean(temp_F2[RUN_IDX, :, :], axis=0)
            prerun_mean2 = np.mean(temp_F2[BASELINE_IDX_RUN, :, :], axis=0)
            dFF_run2 = (run_mean2 - prerun_mean2) / np.abs(prerun_mean2)
            dFF_run2[np.abs(dFF_run2) > 10] = np.nan
            pixel_dFF_run2[:, :, i] = dFF_run2

            del temp_F, temp_F2
            clear_memory()

        np.save(processed_datapath / f'{recname}_pixel_dFF_run.npy',
                pixel_dFF_run)
        np.save(processed_datapath / f'{recname}_pixel_dFF_ch2_run.npy',
                pixel_dFF_run2)

        # compute mean run-onset aligned release maps
        print('computing mean run-onset release map...')
        release_map_run = np.nanmean(pixel_dFF_run, axis=2)
        release_map_run2 = np.nanmean(pixel_dFF_run2, axis=2)

        # save map matrices
        np.save(processed_datapath / f'{recname}_release_map_run.npy', release_map_run)
        np.save(processed_datapath / f'{recname}_release_map_run_ch2.npy', release_map_run2)

        if make_plots:
            save_release_map_figure(
                figurepath, recname, ref, ref2, release_map_run,
                'release_map_run', 'run-onset / baseline (mean)',
                positive_floor=-.01
                )
            save_release_map_figure(
                figurepath, recname, ref, ref2, release_map_run2,
                'release_map_run_ch2', 'stim / baseline (mean)'
                )

    else:
        print('session with no behaviour; finishing...')

    # axon-only imaging session check and processing
    axon_only_folder = path + '_1100'
    if Path(axon_only_folder).exists():  # if we have a 1100-nm wavelength session
        print(f'axon-only recording found: {axon_only_folder}')

        axon_only_plane_stem = pp.resolve_suite2p_session_stem(axon_only_folder) / 'plane0'
        axon_only_bin_path = axon_only_plane_stem / 'data_chan2.bin'
        axon_only_ops_path = axon_only_plane_stem / 'ops.npy'

        if (axon_only_bin_path.exists() and
            axon_only_ops_path.exists()):
            ops_axon_only = np.load(
                axon_only_ops_path, allow_pickle=True
                ).item()
            tot_frames_axon_only = ops_axon_only['nframes']
            shape_axon_only = (tot_frames_axon_only,
                               ops_axon_only['Ly'], ops_axon_only['Lx'])

            print('loading axon-only movie...')
            mov_axon_only = np.memmap(
                axon_only_bin_path,
                mode='r', dtype='int16', shape=shape_axon_only
                ).astype(np.float32)

            print('computing and saving axon-only reference..')
            reference_axon_only = np.mean(mov_axon_only, axis=0)
            reference_axon_only = ipf.post_processing_suite2p_gui(
                reference_axon_only
                )

            if make_plots:
                fig, ax = plt.subplots(figsize=(4,4))
                ax.imshow(reference_axon_only,
                          aspect='auto', cmap='gist_gray', interpolation='none',
                          extent=[0, 512, 512, 0])

                ax.set(xlim=(0,512), ylim=(0,512))

                fig.suptitle('ref 1100 nm')
                fig.tight_layout()
                fig.savefig(figurepath / f'{recname}_ref_1100nm.png',
                            dpi=300,
                            bbox_inches='tight')
                plt.close(fig)

            np.save(processed_datapath / f'{recname}_ref_mat_1100nm.npy',
                    reference_axon_only)


#%% main
def main(argv=None):
    parser = argparse.ArgumentParser(
        description='extract GRABNE LC-opto imaging data.'
    )
    parser.add_argument(
        '--recording',
        dest='recording_filter',
        help='optional substring filter for recname or recording path',
    )
    parser.add_argument(
        '--plots',
        dest='make_plots',
        action='store_true',
        help='also regenerate figure outputs during extraction; default is data only',
    )
    parser.add_argument(
        '--no-plots',
        dest='make_plots',
        action='store_false',
        help=argparse.SUPPRESS,
    )
    parser.set_defaults(make_plots=False)
    args = parser.parse_args(argv)

    session_paths = default_paths
    if args.recording_filter:
        session_paths = [
            path for path in default_paths
            if args.recording_filter in path
            or args.recording_filter in Path(path).name
        ]
        if not session_paths:
            raise ValueError(f'no recordings matched filter: {args.recording_filter}')

    for path in session_paths:
        process_session(path, make_plots=args.make_plots)
        clear_memory(Path(path).name)

if __name__ == '__main__':
    main()
