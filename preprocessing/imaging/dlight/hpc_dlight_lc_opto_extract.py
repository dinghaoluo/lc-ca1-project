# -*- coding: utf-8 -*-
'''
Created on Mon Mar 31 13:02:49 2025
Modified on Tue 24 June 16:24:15 2025
Modified on Tue 25 Nov 2025

extract opto-LC stimulation + dLight imaging data
modification notes:
    - 24 June 2025: removed pixel-wise dFF calculation and replaced it with a
        simplified dFF (stim. / baseline for raw F) that is easy to compute
        and produces the exact same desideratum (the release map)
    - 25 Nov 2025: changed stim.-alignment method to be detection based on
        channel 2 (more stable than before)

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

import imaging_pipeline_functions as ipf
from plotting_functions import plot_violin_with_scatter, add_scale_bar
from common_functions import mpl_formatting, get_GPU_availability
from console_formatting import print_session
import project_paths as pp
mpl_formatting()

import rec_list
paths = rec_list.pathdLightLCOpto + \
        rec_list.pathdLightLCOptoCtrl + \
        rec_list.pathdLightLCOptoInh + \
        rec_list.pathdLightLCOptoDbhBlock
default_paths = paths[31:]

# GPU acceleration
cp, GPU_AVAILABLE, device_name = get_GPU_availability()


#%% parameters
pixel_wise_processing = True
SPATIAL_FILTER_CHUNK_SIZE = 128

SAMP_FREQ = 30

# post-stim dispersion calculation
BIN_WIDTH = 0.1
SENSOR_LABEL = 'dLight'

# path stems
mice_exp_stem = pp.MICEEXP_ROOT
all_sess_stem = pp.HPC_DLIGHT_LC_OPTO_STEM / 'all_sessions'
all_sess_fig_stem = pp.HPC_DLIGHT_LC_OPTO_ALL_SESSIONS_FIGURES_STEM

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
        channel, reference_frames
        ):
    print(f'ch{channel}: loading movie into float32...', flush=True)
    mov_float = np.memmap(
        binpath, mode='r', dtype='int16', shape=shape
        ).astype(np.float32)

    print(f'ch{channel}: saving reference array...', flush=True)
    ref = ipf.plot_reference(
        mov_float, recname=recname, frames=reference_frames,
        outpath=savepath, channel=channel,
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

def save_release_map_figures(savepath, recname, ref, ref2, release_map, release_map2):
    for suffix, release, title in [
            ('', release_map, 'CH1: stim / baseline (mean)'),
            ('_ch2', release_map2, 'CH2: stim / baseline (mean)'),
            ]:
        vmin = np.nanpercentile(release, 1)
        vmax = np.nanpercentile(release, 99)

        if vmin >= 0:
            vmin = -.001
        if vmax <= 0:
            vmax = .001

        norm = TwoSlopeNorm(vcenter=0, vmin=vmin, vmax=vmax)

        fig, axs = plt.subplots(1, 3, figsize=(12, 4))

        axs[0].imshow(ref, cmap='gray', interpolation='none')
        axs[0].set_title('channel 1', fontsize=10)
        axs[0].axis('off')

        axs[1].imshow(ref2, cmap='gray', interpolation='none')
        axs[1].set_title('channel 2', fontsize=10)
        axs[1].axis('off')

        im = axs[2].imshow(release, cmap='RdBu_r', norm=norm, interpolation='none')
        axs[2].set_title(title, fontsize=10)
        axs[2].axis('off')

        cbar = fig.colorbar(im, ax=axs[2], shrink=0.8, fraction=0.046, pad=0.04)
        cbar.set_label('ΔF/F ratio', fontsize=10)
        cbar.set_ticks([vmin, 0, vmax])

        fig.tight_layout()

        for ext in ['.png', '.pdf']:
            fig.savefig(
                savepath / f'{recname}_release_map{suffix}{ext}',
                dpi=300,
                bbox_inches='tight'
            )

        cmap = colormaps['RdBu_r']
        release_rgba = cmap(norm(release))
        release_rgb = (release_rgba[..., :3] * 255).astype(np.uint8)
        tifffile.imwrite(savepath / f'{recname}_release_map{suffix}.tiff',
                         release_rgb)
        plt.close(fig)


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

    whether_ctrl = '_ctrl' if path in rec_list.pathdLightLCOptoCtrl else ''
    savepath = all_sess_stem / f'{recname}{whether_ctrl}'
    savepath.mkdir(parents=True, exist_ok=True)
    processed_datapath = savepath / 'processed_data'
    processed_datapath.mkdir(parents=True, exist_ok=True)
    figurepath = all_sess_fig_stem / f'{recname}{whether_ctrl}'
    if make_plots:
        figurepath.mkdir(parents=True, exist_ok=True)

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
        channel=1, reference_frames=tot_frames
        )
    mov2, ref2, raw_trace2 = prepare_movie_channel(
        bin2path, shape, recname, savepath, figurepath, make_plots,
        channel=2, reference_frames=tot_frames
        )

    tot_frames = mov.shape[0]  # once loaded, update tot_frames to be the max frame number, 16 June 2025

    # dFF traces are ONLY used for plotting figure 1 now
    print('computing dFF traces...')
    trace_dFF  = ipf.calculate_dFF_percentile(raw_trace,
                                              t_axis=0,
                                              window_size=9000,
                                              pct=15,  # changed from default (10), 5 Feb 2026
                                              GPU_AVAILABLE=GPU_AVAILABLE,
                                              device_name=device_name,
                                              progress_desc='Ch1 whole-field dF/F baseline')
    trace2_dFF = ipf.calculate_dFF_percentile(raw_trace2,
                                              t_axis=0,
                                              window_size=9000,
                                              pct=15,  # changed from default (10), 5 Feb 2026
                                              GPU_AVAILABLE=GPU_AVAILABLE,
                                              device_name=device_name,
                                              progress_desc='Ch2 whole-field dF/F baseline')

    # ---------------------------
    # determine stim. timestamps
    # ---------------------------
    # we now use as the primary method derivatives of channel 2 signals to
    #   determine where step changes happened (both up and down), and only fall
    #   back to the text files for timestamps if this fails

    # we first read the .txt file to figure out how long the pulse trains are
    #   in this recording and to plot the example pulse trace
    print('retrieving stim. parameters...')
    txt = ipf.process_txt_nobeh(txtpath)

    # load frame times and check $FM against tot_frame
    frame_times = txt['frame_times']
    if tot_frames < len(frame_times) - 3 or tot_frames > len(frame_times):
        raise ValueError(f'check $FM; frame count mismatch for {recname}')
    pulse_info = ipf.extract_opto_pulse_metadata(txt, SAMP_FREQ)
    pulse_trains = pulse_info['pulse_trains']
    pulse_width_ON = pulse_info['pulse_width_on']
    pulse_width = pulse_info['pulse_width']
    pulse_number = pulse_info['pulse_number']
    taper_enabled = pulse_info['taper_enabled']
    taper_duration = pulse_info['taper_duration']
    duty_cycle = pulse_info['duty_cycle']
    min_interval_frames = pulse_info['min_interval_frames']

    ## now we detect stim onsets by step changes in trace2_dFF
    detected_onsets_raw, detected_offsets_raw = ipf.detect_step_pairs(
        trace2_dFF,
        zthr=10,
        min_interval_frames=min_interval_frames
        )
    detected_onsets  = np.array([f - 3 for f in detected_onsets_raw])  # 3-frame buffer for envelope
    detected_offsets = np.array([f + 3 for f in detected_offsets_raw])  # same as above
    valid_edges      = [True if off - on < min_interval_frames else False
                        for on, off in zip(detected_onsets, detected_offsets)]
    detected_onsets  = detected_onsets[valid_edges]
    detected_offsets = detected_offsets[valid_edges]
    detected_stim_durations = [off-on for off, on in zip(detected_offsets, detected_onsets)]

    # if detected and not over-detected
    if abs(len(detected_onsets) - len(pulse_trains)) > 10:  # if mismatched
        print('mismatch between detection and logged pulses; rerunning detection with a higher threshold...')
        detected_onsets_raw, detected_offsets_raw = ipf.detect_step_pairs(
            trace2_dFF,
            zthr=20,
            min_interval_frames=min_interval_frames
            )
        detected_onsets  = np.array([f - 3 for f in detected_onsets_raw])  # 1-frame buffer for envelope
        detected_offsets = np.array([f + 3 for f in detected_offsets_raw])  # same as above
        valid_edges      = [True if off - on < min_interval_frames else False
                            for on, off in zip(detected_onsets, detected_offsets)]
        detected_onsets  = detected_onsets[valid_edges]
        detected_offsets = detected_offsets[valid_edges]
        detected_stim_durations = [off-on for off, on in zip(detected_offsets, detected_onsets)]
        if abs(len(detected_onsets) - len(pulse_trains)) > 10:  # if still not matching
            print('mismatch persists; fall back to logged pulses.')
            pulse_frames = [
                [ipf.find_nearest(p, frame_times) for p in train]
                for train in pulse_trains
            ]
            detected_onsets  = np.array([pf[0] - 3 for pf in pulse_frames])
            detected_offsets = np.array([pf[-1] + 3 for pf in pulse_frames])
            valid_edges      = [True if off - on < min_interval_frames else False
                                for on, off in zip(detected_onsets, detected_offsets)]
            detected_onsets  = detected_onsets[valid_edges]
            detected_offsets = detected_offsets[valid_edges]
            detected_stim_durations = [off-on for off, on in zip(detected_offsets, detected_onsets)]
            max_stim_duration_s = max(detected_stim_durations) / SAMP_FREQ
        else:
            max_stim_duration_s = max(detected_stim_durations) / SAMP_FREQ
    elif 0 < len(detected_onsets):
        max_stim_duration_s = max(detected_stim_durations) / SAMP_FREQ
        detection_printout = (f'''Detected based on channel 2:
            {len(detected_onsets)} candidate stim. onset-offset pairs
            max stim. duration: {max(detected_stim_durations)} frames ({max_stim_duration_s} s)''')
        print(detection_printout)
    else:
        print('no detected onset-offset pairs; falling back to .txt pulse_times')
        pulse_frames = [
            [ipf.find_nearest(p, frame_times) for p in train]
            for train in pulse_trains
        ]
        detected_onsets  = np.array([pf[0] - 3 for pf in pulse_frames])
        detected_offsets = np.array([pf[-1] + 3 for pf in pulse_frames])
        valid_edges      = [True if off - on < min_interval_frames else False
                            for on, off in zip(detected_onsets, detected_offsets)]
        detected_onsets  = detected_onsets[valid_edges]
        detected_offsets = detected_offsets[valid_edges]
        detected_stim_durations = [off-on for off, on in zip(detected_offsets, detected_onsets)]
        max_stim_duration   = max(detected_stim_durations)
        max_stim_duration_s = max_stim_duration / SAMP_FREQ

    ## -- PARAMETER DEFINITIONS
    # now defined within the function scope, since we reassign them later in an if statement
    BEF = 2
    AFT = 10
    TAXIS = np.arange(-BEF*SAMP_FREQ, AFT*SAMP_FREQ) / SAMP_FREQ
    BASELINE_IDX = (TAXIS >= -1.0) & (TAXIS < 0)
    STIM_IDX = (TAXIS > max_stim_duration_s) & (TAXIS <= max_stim_duration_s + 1.0)  # note that this is the mask for extracting stim_mean
    ## -- END PARAMETER DEFINITIONS

    ## finally we filter out valid start frames
    valid_pulse_start_frames = [
        on for on, off in zip(detected_onsets, detected_offsets)
        if on - BEF*SAMP_FREQ >= 0 and off + AFT*SAMP_FREQ <= tot_frames
    ]
    tot_valid_pulses = len(valid_pulse_start_frames)
    print(f'valid stim. events for alignment: {tot_valid_pulses}')

    # post-stim dispersion calculation, 10 Sept 2025
    BIN_START = max_stim_duration_s
    BIN_END   = max_stim_duration_s + 4
    bin_edges = np.arange(BIN_START, BIN_END + BIN_WIDTH, BIN_WIDTH)
    n_bins = len(bin_edges) - 1
    # ---------------------------
    # end
    # ---------------------------

    # pulse processing
    print('extracting data...')

    # filter for opto artefact periods
    pulse_period_frames = np.concatenate([
        np.arange(on, off+1) for on, off in zip(detected_onsets, detected_offsets)
        if on - BEF*SAMP_FREQ >= 0 and off + AFT*SAMP_FREQ <= tot_frames
        ])

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
        raw_aligned[i, :]   = raw_trace[start:end]
        raw2_aligned[i, :] = raw_trace2[start:end]

    # dFF traces aligned
    trace_dFF_aligned  = np.zeros((tot_valid_pulses, (BEF+AFT)*SAMP_FREQ), dtype=np.float32)
    trace2_dFF_aligned = np.zeros((tot_valid_pulses, (BEF+AFT)*SAMP_FREQ), dtype=np.float32)
    for i, p in enumerate(valid_pulse_start_frames):
        start = p - BEF * SAMP_FREQ
        end   = p + AFT * SAMP_FREQ
        trace_dFF_aligned[i, :]  = trace_dFF[start:end]
        trace2_dFF_aligned[i, :] = trace2_dFF[start:end]

    # calculate the mean traces for plotting
    trace_dFF_aligned_mean  = np.nanmean(trace_dFF_aligned, axis=0)
    trace2_dFF_aligned_mean = np.nanmean(trace2_dFF_aligned, axis=0)

    # block out the stim period
    blocked_on  = BEF * SAMP_FREQ - 1  # -1 frame as a buffer
    blocked_off = BEF * SAMP_FREQ + int(max_stim_duration_s * SAMP_FREQ)
    trace_dFF_aligned_mean[blocked_on : blocked_off]  = np.nan
    trace2_dFF_aligned_mean[blocked_on : blocked_off] = np.nan

    # calculate ratios
    # per‐trial raw means
    baseline_raw  = np.nanmean(raw_aligned[:,  BASELINE_IDX], axis=1)
    stim_raw      = np.nanmean(raw_aligned[:,  STIM_IDX], axis=1)
    baseline2_raw = np.nanmean(raw2_aligned[:, BASELINE_IDX], axis=1)
    stim2_raw     = np.nanmean(raw2_aligned[:, STIM_IDX], axis=1)

    # per‐trial ΔF/F exactly like pixel dFF (stim − base) / |base|
    dFF  = (stim_raw - baseline_raw) / np.abs(baseline_raw)
    dFF2 = (stim2_raw - baseline2_raw) / np.abs(baseline2_raw)

    # dFF comp
    baseline_dFF  = np.nanmean(trace_dFF_aligned[:,  BASELINE_IDX], axis=1)
    stim_dFF      = np.nanmean(trace_dFF_aligned[:,  STIM_IDX], axis=1)
    baseline2_dFF = np.nanmean(trace2_dFF_aligned[:, BASELINE_IDX], axis=1)
    stim2_dFF     = np.nanmean(trace2_dFF_aligned[:, STIM_IDX], axis=1)

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
        'total_pulse_trains': len(detected_onsets),
        }
    np.save(processed_datapath / f'{recname}_opto_plot_payload.npy',
            plot_payload)

    # saving
    np.save(processed_datapath / f'{recname}_wholefield_dFF_stim.npy',
            trace_dFF_aligned_mean)
    np.save(processed_datapath / f'{recname}_wholefield_dFF2_stim.npy',
            trace2_dFF_aligned_mean)

    if make_plots:
        ymin  = np.nanmin(trace_dFF_aligned.T)
        ymin2 = np.nanmin(trace2_dFF_aligned.T)

        fig, axs = plt.subplots(3,1,figsize=(3.5,5),
                                sharex=True)

        axs[0].plot(TAXIS, trace_dFF_aligned_mean, color='green', linewidth=2)
        axs[0].plot(TAXIS, trace_dFF_aligned.T, color='green', alpha=.05)
        axs[0].text(-2, ymin-.01, 'dLight', ha='left', va='center', color='green', fontsize=10)

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
                f'total pulse trains = {len(detected_onsets)}\n'
                f'taper duration = {taper_duration}'
                )
        else:
            title_str = (
                f'{recname}\n'
                f'duty cycle = {duty_cycle}\n'
                f'pulse width = {pulse_width} s\n'
                f'pulse(s) per pulse train = {pulse_number}\n'
                f'total pulse trains = {len(detected_onsets)}\n'
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
            xticklabels=['ref.', 'dLight'],
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

    # ----------------------
    # pixel-wise extraction
    # ----------------------
    # truncate movies first to reduce load
    # spatial smoothing
    if not pixel_wise_processing:
        clear_memory(f'{recname} non-pixel extraction')
        return

    clear_memory(f'{recname} trace extraction')
    print('performing trial-window spatial filtering and pixel extraction...')

    # compute dF/F per pixel (stim. / baseline for raw F), 24 June 2025
    pixel_RI  = np.zeros((shape[1], shape[2], len(valid_pulse_start_frames)),
                         dtype=np.float32)
    pixel_RI2 = np.zeros_like(pixel_RI)

    # we want F and dF/F aligned
    pixel_F_aligned    = np.zeros((len(valid_pulse_start_frames),
                                   ((BEF+AFT) * SAMP_FREQ),
                                   shape[1],
                                   shape[2]),
                                  dtype=np.float32)
    pixel_F2_aligned   = np.zeros_like(pixel_F_aligned)

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

        # mask pulse frames in each aligned window so F_aligned blocks them out, 28 Nov 2025
        trial_pulse_frames = pulse_period_frames[
            (pulse_period_frames >= start) & (pulse_period_frames < end)
            ] - start
        temp_F[trial_pulse_frames, :, :] = np.nan
        temp_F2[trial_pulse_frames, :, :] = np.nan

        # save F aligned too
        pixel_F_aligned[i, :, :, :]    = temp_F
        pixel_F2_aligned[i, :, :, :]   = temp_F2

        stim_mean     = np.mean(temp_F[STIM_IDX, :, :], axis=0)
        baseline_mean = np.mean(temp_F[BASELINE_IDX, :, :], axis=0)
        RI = (stim_mean - baseline_mean) / np.abs(baseline_mean)
        RI[np.abs(RI) > 10] = np.nan  # hard cap
        pixel_RI[:, :, i] = RI

        stim_mean2 = np.mean(temp_F2[STIM_IDX, :, :], axis=0)
        baseline_mean2 = np.mean(temp_F2[BASELINE_IDX, :, :], axis=0)
        RI2 = (stim_mean2 - baseline_mean2) / np.abs(baseline_mean2)
        RI2[np.abs(RI2) > 10] = np.nan
        pixel_RI2[:, :, i] = RI2

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

    np.save(processed_datapath / f'{recname}_pixel_RI_stim.npy',
            pixel_RI)
    np.save(processed_datapath / f'{recname}_pixel_RI_ch2_stim.npy',
            pixel_RI2)

    # save aligned
    np.save(processed_datapath / f'{recname}_pixel_F_aligned.npy',
            pixel_F_aligned)
    np.save(processed_datapath / f'{recname}_pixel_F2_aligned.npy',
            pixel_F2_aligned)
    # np.save(processed_datapath / f'{recname}_pixel_dFF_aligned.npy',
    #         pixel_dFF_aligned)
    # np.save(processed_datapath / f'{recname}_pixel_dFF2_aligned.npy',
    #         pixel_dFF2_aligned)

    # generate mean dFF release map as proxy for t-map, 24 June 2025
    print('computing mean release map...')
    release_map = np.nanmean(pixel_RI, axis=2)  # shape: (y, x)
    release_map2 = np.nanmean(pixel_RI2, axis=2)

    # save map matrices
    np.save(processed_datapath / f'{recname}_release_map.npy', release_map)
    np.save(processed_datapath / f'{recname}_release_map_ch2.npy', release_map2)

    if make_plots:
        save_release_map_figures(
            figurepath, recname, ref, ref2, release_map, release_map2
            )

    # dispersion rate calculation, 10 Sept 2025
    print('calculating binned ratios (for dispersion rate analysis)...')
    # average across trials → final shape (y, x, n_bins)
    pixel_RI_bins  = np.nanmean(np.stack(all_bins_ch1, axis=-1), axis=-1)
    pixel_RI2_bins = np.nanmean(np.stack(all_bins_ch2, axis=-1), axis=-1)

    # save arrays
    np.save(processed_datapath / f'{recname}_pixel_RI_bins.npy', pixel_RI_bins)
    np.save(processed_datapath / f'{recname}_pixel_RI2_bins.npy', pixel_RI2_bins)

    ## ---- axon-only imaging session check and processing
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
    ## ---- axon only reference map ends


#%% main
def main(argv=None):
    parser = argparse.ArgumentParser(
        description='extract dLight LC-opto imaging data.'
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
