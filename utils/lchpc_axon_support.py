# -*- coding: utf-8 -*-
'''
Created on Fri Mar 21 14:21:07 2025

ROI geometry and reference-image plots for the LC-CA1 axon-GCaMP analyses

@author: Dinghao Luo
'''

#%% imports
import numpy as np
import os
from tqdm import tqdm
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.gridspec import GridSpec
from skimage.measure import find_contours

import imaging_pipeline_functions as ipf

from common_functions import (
    mpl_formatting,
    normalise,
    smooth_convolve,
    )
mpl_formatting()


#%% functions for lchpc_axon_all_extract.py

def calculate_and_plot_overlap_indices(
        ref_im,
        ref_ch2_im,
        stat,
        valid_rois,
        recname,
        figure_path,
        border=10
        ):
    '''calculate channel-overlap indices and save per-ROI validation plots'''
    overlap_indices = {}
    output_dir = os.path.join(figure_path, 'roi_ch2_validation')
    os.makedirs(output_dir, exist_ok=True)

    for roi in valid_rois:
        xpix = stat[roi]['xpix']
        ypix = stat[roi]['ypix']

        ch2_values = ref_ch2_im[ypix, xpix]
        ch1_values = ref_im[ypix, xpix]
        overlap_index = ch2_values.mean() / ch1_values.mean()
        overlap_indices[roi] = overlap_index

        x_min, x_max = xpix.min(), xpix.max()
        y_min, y_max = ypix.min(), ypix.max()
        x_center = (x_min + x_max) // 2
        y_center = (y_min + y_max) // 2
        half_span = max(x_max - x_min, y_max - y_min) // 2 + border
        x_min_sq = max(0, x_center - half_span)
        x_max_sq = min(ref_im.shape[1], x_center + half_span)
        y_min_sq = max(0, y_center - half_span)
        y_max_sq = min(ref_im.shape[0], y_center + half_span)

        ch1_sub = ref_im[y_min_sq:y_max_sq, x_min_sq:x_max_sq]
        ch2_sub = ref_ch2_im[y_min_sq:y_max_sq, x_min_sq:x_max_sq]

        ch1_proc = ipf.post_processing_suite2p_gui(ch1_sub)
        ch2_proc = ipf.post_processing_suite2p_gui(ch2_sub)

        roi_mask = np.zeros_like(ch1_sub, dtype=bool)
        for x, y in zip(xpix, ypix):
            x_rel = x - x_min_sq
            y_rel = y - y_min_sq
            if 0 <= x_rel < roi_mask.shape[1] and 0 <= y_rel < roi_mask.shape[0]:
                roi_mask[y_rel, x_rel] = True

        # find_contours returns row, column coordinates
        contours = find_contours(roi_mask.astype(float), level=0.5)

        fig, axes = plt.subplots(1, 3, figsize=(6, 2.2))
        for ax in axes:
            ax.axis('off')

        axes[0].imshow(ch1_proc, cmap='gray')  # raw channel 1
        axes[0].set_title('ch1 raw')

        axes[1].imshow(ch1_proc, cmap='gray')  # channel 1 with ROI overlay
        # axes[1].scatter(xpix - x_min_sq, ypix - y_min_sq,
        #                 color='limegreen', s=1, edgecolor='none', alpha=0.5)
        for contour in contours:
            axes[1].plot(contour[:, 1], contour[:, 0],
                         linewidth=1, color='limegreen')
        axes[1].set_title(f'ch1 + ROI ({round(overlap_index, 3)})')

        axes[2].imshow(ch2_proc, cmap='gray')  # channel 2
        axes[2].set_title('ch2')

        fig.suptitle(f'ROI {roi} overlap')
        fig.tight_layout()

        for ext in ['.png', '.pdf']:
            fig.savefig(os.path.join(output_dir, f'roi_{roi}{ext}'),
                        dpi=300,
                        bbox_inches='tight')
        plt.close(fig)

    return overlap_indices

def filter_valid_rois(stat):
    '''
    retain the longest merged ROI for overlapping merge histories

    Serial merges can produce subset duplicates from the same constituents; keep
    the longest merge in each overlap chain (13 Nov 2024, Dinghao).
    '''
    valid_rois_dict = {}

    sorted_rois = sorted(range(len(stat)), key=lambda roi: len(stat[roi]['imerge']), reverse=True)
    covered_constituents = set()

    for roi in sorted_rois:
        imerge_set = set(stat[roi]['imerge'])

        if not imerge_set.issubset(covered_constituents):
            valid_rois_dict[roi] = list(imerge_set)
            covered_constituents.update(imerge_set)

    return valid_rois_dict

def extract_spatial_median_roi_traces(
        mov,
        roi_coords_dict,
        size=5,
        GPU_AVAILABLE=False,
        chunk_size=500,
        output_dtype=np.float32):
    '''
    Apply spatial median filtering in chunks and keep only ROI pixel traces.

    roi_coords_dict values are expected to be [xpix, ypix], matching Suite2p's
    stat fields. Returned arrays have shape (time, pixels).
    '''
    T = mov.shape[0]
    roi_specs = []
    all_x = []
    all_y = []
    offset = 0

    for roi, coords in roi_coords_dict.items():
        xpix, ypix = coords
        xpix = np.asarray(xpix, dtype=np.intp)
        ypix = np.asarray(ypix, dtype=np.intp)

        n_pix = xpix.size
        traces = np.empty((T, n_pix), dtype=output_dtype)
        roi_specs.append((roi, offset, n_pix, traces))
        all_x.append(xpix)
        all_y.append(ypix)
        offset += n_pix

    all_x = np.concatenate(all_x)
    all_y = np.concatenate(all_y)

    if GPU_AVAILABLE:
        import cupy as cp
        import cupyx.scipy.ndimage

        all_x_gpu = cp.asarray(all_x)
        all_y_gpu = cp.asarray(all_y)
        iterator = tqdm(
            range(0, T, chunk_size),
            desc='chunk median-filtering ROI pixels on GPU...'
            )

        for start in iterator:
            end = min(start + chunk_size, T)
            chunk = cp.asarray(mov[start:end])
            filtered_chunk = cupyx.scipy.ndimage.median_filter(
                chunk, size=(1, size, size)
                )
            filtered_pixels = filtered_chunk[:, all_y_gpu, all_x_gpu].get()
            filtered_pixels = filtered_pixels.astype(output_dtype, copy=False)

            for _, roi_offset, n_pix, traces in roi_specs:
                traces[start:end] = filtered_pixels[
                    :, roi_offset:roi_offset + n_pix
                    ]

            del chunk, filtered_chunk, filtered_pixels
            cp.get_default_memory_pool().free_all_blocks()
    else:
        from scipy.ndimage import median_filter as cpu_median_filter
        iterator = tqdm(
            range(0, T, chunk_size),
            desc='chunk median-filtering ROI pixels on CPU...'
            )

        for start in iterator:
            end = min(start + chunk_size, T)
            filtered_chunk = cpu_median_filter(
                mov[start:end],
                size=(1, size, size)
                )
            filtered_pixels = filtered_chunk[:, all_y, all_x].astype(
                output_dtype,
                copy=False
                )

            for _, roi_offset, n_pix, traces in roi_specs:
                traces[start:end] = filtered_pixels[
                    :, roi_offset:roi_offset + n_pix
                    ]

    return {roi: traces for roi, _, _, traces in roi_specs}

def get_roi_coord_dict(
        ref_im, ref_ch2_im,
        stat, rois, recname, figure_path,
        plot=True):
    '''
    generate a dictionary of ROI pixel coordinates and optionally save a 3-panel reference plot.

    parameters:
    - ref_im: np.ndarray
        reference image for channel 1
    - ref_ch2_im: np.ndarray
        reference image for channel 2
    - stat: list of dict
        list of ROI dictionaries, each containing 'xpix' and 'ypix' with ROI pixel coordinates
    - rois: iterable
        list or set of ROI indices to include
    - recname: str
        name of the recording session
    - figure_path: str
        path to save the plot if plotting is enabled
    - plot: bool, optional
        whether to generate and save a 3-panel plot showing merged ROIs, channel 1, and channel 2 (default: True)

    returns:
    - roi_coord_dict: dict
        dictionary mapping ROI names (e.g., 'ROI 23') to their [xpix, ypix] coordinates
    '''
    if plot:
        fig, axs = plt.subplots(1, 3, figsize=(6, 2))
        fig.subplots_adjust(wspace=0.35, top=0.75)

        for ax in axs:
            ax.set(xlim=(0, 512), ylim=(0, 512))
            ax.set_aspect('equal')
            ax.set_xticks([])  # remove x ticks
            ax.set_yticks([])  # remove y ticks

        # custom colour maps for channels 1 and 2
        colors_ch1 = plt.cm.Greens(np.linspace(0, 0.8, 256))
        colors_ch2 = plt.cm.Reds(np.linspace(0, 0.8, 256))
        custom_cmap_ch1 = LinearSegmentedColormap.from_list('mycmap_ch1', colors_ch1)
        custom_cmap_ch2 = LinearSegmentedColormap.from_list('mycmap_ch2', colors_ch2)

        # display reference images in channels 1 and 2
        axs[0].set(title='merged ROIs')
        axs[1].imshow(ref_im, cmap=custom_cmap_ch1)
        axs[2].imshow(ref_ch2_im, cmap=custom_cmap_ch2)
        axs[1].set(title='axon-GCaMP')
        axs[2].set(title='Dbh:Ai14')

        for roi in rois:
            axs[0].scatter(stat[roi]['xpix'], stat[roi]['ypix'],
                           edgecolor='none', s=0.1, alpha=0.2)

        fig.suptitle(recname)
        for ext in ['.png', '.pdf']:
            fig.savefig(os.path.join(figure_path, f'rois_v_ref{ext}'), dpi=200)
        plt.close(fig)

    # store ROI coords in roi_dict
    roi_coord_dict = {}
    for roi in rois:
        roi_coord_dict[f'ROI {roi}'] = [stat[roi]['xpix'], stat[roi]['ypix']]

    return roi_coord_dict

def prepare_merged_roi_metadata(
        ref_im, ref_ch2_im, stat, recname, proc_path,
        plot=True, figure_path=None):
    '''
    prepare final/constituent ROI sets, ROI coordinate dictionaries, and
    validation plots for merged axon-GCaMP Suite2p outputs.
    '''
    valid_rois_dict = filter_valid_rois(stat)
    valid_rois = set(valid_rois_dict)
    figure_path = proc_path if figure_path is None else figure_path

    constituent_rois = {
        roi
        for sublist in valid_rois_dict.values()
        for roi in sublist
        }
    all_rois = valid_rois | constituent_rois

    if plot:
        calculate_and_plot_overlap_indices(
            ref_im, ref_ch2_im,
            stat, valid_rois, recname, figure_path
            )

    constituent_rois_coord_dict = get_roi_coord_dict(
        ref_im, ref_ch2_im,
        stat, constituent_rois, recname, figure_path,
        plot=plot
        )
    valid_rois_coord_dict = get_roi_coord_dict(
        ref_im, ref_ch2_im,
        stat, valid_rois, recname, figure_path,
        plot=plot
        )

    return (
        valid_rois_dict,
        valid_rois,
        constituent_rois,
        all_rois,
        constituent_rois_coord_dict,
        valid_rois_coord_dict,
        )

def save_roi_alignment_plots(
        *,
        ref_im,
        stat,
        roi,
        aligned_im,
        aligned_mean,
        aligned_sem,
        aligned_mean_ch2,
        aligned_sem_ch2,
        xaxis,
        before,
        after,
        total_events,
        out_dir,
        event_title,
        xlabel,
        overview_xlabel=None,
        overlay_xlabel=None,
        dual_axis_ch2=False
        ):
    '''
    save the standard single-ROI alignment figures.

    This is plotting-only. It receives precomputed aligned traces and summary
    statistics so callers can keep processing and figure generation separate.
    '''
    os.makedirs(out_dir, exist_ok=True)
    overview_xlabel = overview_xlabel or xlabel
    overlay_xlabel = overlay_xlabel or xlabel

    fig = plt.figure(figsize=(5, 2.5))
    gs = GridSpec(2, 2, width_ratios=[1, 1], height_ratios=[1, 1])
    ax1 = fig.add_subplot(gs[:, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, 1])
    fig.subplots_adjust(wspace=.4)

    ax1.imshow(ref_im, cmap='gist_gray')
    ax1.scatter(stat[roi]['xpix'], stat[roi]['ypix'],
                color='limegreen', edgecolor='none', s=.1, alpha=1)
    ax1.plot(stat[roi]['xcirc'], stat[roi]['ycirc'],
             linewidth=.5, color='white')
    ax1.set(xlim=(0, 512), ylim=(0, 512))

    ax2.imshow(aligned_im,
               cmap='Greys', extent=[-before, after, 0, total_events],
               aspect='auto')
    ax2.set(xticklabels=[],
            ylabel='trial #')

    ax3.plot(xaxis, aligned_mean, c='darkgreen', linewidth=1)
    ax3.fill_between(xaxis, aligned_mean + aligned_sem,
                     aligned_mean - aligned_sem,
                     color='darkgreen', alpha=.2, edgecolor='none')
    ax3.set(xlabel=overview_xlabel,
            ylabel='dF/F')

    fig.suptitle(f'ROI {roi} {event_title}')
    for ext in ('.png', '.pdf'):
        fig.savefig(os.path.join(out_dir, f'roi_{roi}{ext}'),
                    dpi=300,
                    bbox_inches='tight')
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(2.5, 1.7))
    if dual_axis_ch2:
        axt = ax.twinx()
        axt.set_zorder(0)
        ax.set_zorder(1)
        ax.patch.set_visible(False)

        ref_ca_ln, = axt.plot(xaxis, aligned_mean_ch2,
                              linewidth=1,
                              c='indianred', label='tdT', alpha=.5)
        axt.fill_between(xaxis, aligned_mean_ch2 + aligned_sem_ch2,
                         aligned_mean_ch2 - aligned_sem_ch2,
                         color='indianred', alpha=.1, edgecolor='none')

        ca_ln, = ax.plot(xaxis, aligned_mean,
                         linewidth=1,
                         c='darkgreen', label='GCaMP')
        ax.fill_between(xaxis, aligned_mean + aligned_sem,
                        aligned_mean - aligned_sem,
                        color='darkgreen', alpha=.2, edgecolor='none')

        handles = [ca_ln, ref_ca_ln]
        labels = [h.get_label() for h in handles]
        ax.legend(handles, labels, fontsize=6, frameon=False)
        axt.set(ylabel='dF/F')
        axt.spines['top'].set_visible(False)
    else:
        ref_ca_ln, = ax.plot(xaxis, aligned_mean_ch2,
                             linewidth=1,
                             c='indianred', label='tdT', alpha=.5)
        ax.fill_between(xaxis, aligned_mean_ch2 + aligned_sem_ch2,
                        aligned_mean_ch2 - aligned_sem_ch2,
                        color='indianred', alpha=.1, edgecolor='none')

        ca_ln, = ax.plot(xaxis, aligned_mean,
                         linewidth=1,
                         c='darkgreen', label='GCaMP',
                         zorder=10)
        ax.fill_between(xaxis, aligned_mean + aligned_sem,
                        aligned_mean - aligned_sem,
                        color='darkgreen', alpha=.2, edgecolor='none',
                        zorder=10)

        ax.legend(fontsize=6, frameon=False)
        ax.spines['right'].set_visible(False)

    ax.set(xlabel=overlay_xlabel,
           ylabel='dF/F')
    ax.spines['top'].set_visible(False)

    fig.suptitle(f'ROI {roi} {event_title}')
    for ext in ('.png', '.pdf'):
        fig.savefig(os.path.join(out_dir, f'roi_{roi}_w_ch2{ext}'),
                    dpi=300,
                    bbox_inches='tight')
    plt.close(fig)

def align_dual_channel_traces(
        ca,
        ref_ca,
        event_frames,
        before,
        after,
        sample_freq,
        total_frames
        ):
    '''
    align green and red-channel traces to monotonically increasing events.

    Returns raw aligned traces, display-normalised aligned traces, and the
    pre-truncation event count used by the archive plotting extent.
    '''
    filtered_frames = []
    last_frame = float('-inf')
    for frame in event_frames:
        if frame > last_frame:
            filtered_frames.append(frame)
            last_frame = frame

    total_events = len(filtered_frames)
    head = 0
    tail = len(filtered_frames)

    for frame in filtered_frames:
        if frame - before * sample_freq < 0:
            head += 1
        else:
            break

    for frame in reversed(filtered_frames):
        if frame + after * sample_freq > total_frames:
            tail -= 1
        else:
            break

    total_truncated = head + (len(filtered_frames) - tail)
    n_samples = (before + after) * sample_freq
    n_aligned = total_events - total_truncated

    aligned = np.zeros((n_aligned, n_samples))
    aligned_im = np.zeros((n_aligned, n_samples))
    aligned_ch2 = np.zeros((n_aligned, n_samples))
    aligned_im_ch2 = np.zeros((n_aligned, n_samples))

    for i, frame in enumerate(filtered_frames[head:tail]):
        segment = ca[frame - before * sample_freq:frame + after * sample_freq]
        segment_ch2 = ref_ca[frame - before * sample_freq:frame + after * sample_freq]
        aligned[i, :] = segment
        aligned_im[i, :] = normalise(smooth_convolve(segment))
        aligned_ch2[i, :] = segment_ch2
        aligned_im_ch2[i, :] = normalise(smooth_convolve(segment_ch2))

    return aligned, aligned_im, aligned_ch2, aligned_im_ch2, total_events
