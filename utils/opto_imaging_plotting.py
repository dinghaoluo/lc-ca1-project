# -*- coding: utf-8 -*-
'''
Created on 27 May 2026

plot saved LC-opto imaging outputs without reopening raw movies or rerunning
the extraction

@author: Dinghao Luo
'''

#%% imports
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib import colormaps
from matplotlib.colors import TwoSlopeNorm
import numpy as np
import tifffile

from plotting_functions import plot_violin_with_scatter, add_scale_bar


#%% constants
PLOT_PAYLOAD_SUFFIX = 'opto_plot_payload.npy'


#%% aligned-stim summary
def save_aligned_stim_summary(savepath, recname, payload):
    taxis = payload['taxis']
    trace_dFF_aligned = payload['trace_dFF_aligned']
    trace2_dFF_aligned = payload['trace2_dFF_aligned']
    trace_dFF_aligned_mean = payload['trace_dFF_aligned_mean']
    trace2_dFF_aligned_mean = payload['trace2_dFF_aligned_mean']
    sensor_label = payload['sensor_label']

    ymin = np.nanmin(trace_dFF_aligned.T)
    ymin2 = np.nanmin(trace2_dFF_aligned.T)

    fig, axs = plt.subplots(3, 1, figsize=(3.5, 5), sharex=True)

    axs[0].plot(taxis, trace_dFF_aligned_mean, color='green', linewidth=2)
    axs[0].plot(taxis, trace_dFF_aligned.T, color='green', alpha=.05)
    axs[0].text(
        -2, ymin - .01, sensor_label,
        ha='left', va='center', color='green', fontsize=10
        )

    axs[1].plot(taxis, trace2_dFF_aligned_mean, color='darkred', linewidth=2)
    axs[1].plot(taxis, trace2_dFF_aligned.T, color='darkred', alpha=.05)
    axs[1].text(
        -2, ymin2 - .01, 'red ctrl.',
        ha='left', va='center', color='darkred', fontsize=10
        )

    add_scale_bar(axs[0], x_start=8.5, y_start=ymin - .015, x_len=1, y_len=.02)
    add_scale_bar(axs[1], x_start=8.5, y_start=ymin2 - .015, x_len=1, y_len=.02)
    axs[1].text(8.93, ymin2 - .02, '1 s', ha='center', va='top', fontsize=8)
    axs[1].text(
        8.45, ymin2 - .005, '2% ΔF/F',
        ha='right', va='center', rotation='vertical', fontsize=8
        )

    before = payload['before']
    samp_freq = payload['samp_freq']
    pulse_number = payload['pulse_number']
    pulse_width = payload['pulse_width']
    pulse_width_on = payload['pulse_width_on']
    taper_enabled = payload['taper_enabled']
    taper_duration = payload['taper_duration']

    stim_trace = np.zeros_like(taxis)
    for p in range(pulse_number):
        stim_onset = before + p * pulse_width
        stim_offset = stim_onset + pulse_width_on
        stim_onset_idx = int(stim_onset * samp_freq)
        stim_offset_idx = int(stim_offset * samp_freq)
        stim_trace[stim_onset_idx:stim_offset_idx] = 1

    if taper_enabled and taper_duration > 0:
        taper_start_idx = stim_offset_idx
        taper_len = int((taper_duration / 1_000_000) * samp_freq)
        taper = np.linspace(1, 0, taper_len, endpoint=True)
        stim_trace[taper_start_idx : taper_start_idx + taper_len] = (
            taper[:len(stim_trace) - taper_start_idx]
            )

    axs[2].plot(taxis, stim_trace, color='k', linewidth=1)
    axs[2].set_ylim(0, 2)
    axs[2].set_xticks([])
    axs[2].set_yticks([])

    for ax in axs:
        ax.set_xticks([])
        ax.set_yticks([])
        for spine in ['top', 'right', 'left', 'bottom']:
            ax.spines[spine].set_visible(False)

    if payload['taper_enabled']:
        title_str = (
            f'{recname}\n'
            f"duty cycle = {payload['duty_cycle']}\n"
            f"pulse width = {payload['pulse_width']} s\n"
            f"pulse(s) per pulse train = {payload['pulse_number']}\n"
            f"total pulse trains = {payload['total_pulse_trains']}\n"
            f"taper duration = {payload['taper_duration']}"
            )
    else:
        title_str = (
            f'{recname}\n'
            f"duty cycle = {payload['duty_cycle']}\n"
            f"pulse width = {payload['pulse_width']} s\n"
            f"pulse(s) per pulse train = {payload['pulse_number']}\n"
            f"total pulse trains = {payload['total_pulse_trains']}\n"
            )
    fig.suptitle(title_str)
    fig.tight_layout()

    for ext in ['.png', '.pdf']:
        fig.savefig(
            savepath / f'{recname}_aligned_stim{ext}',
            dpi=300,
            bbox_inches='tight'
            )
    plt.close(fig)

    dFF = payload['dFF']
    dFF2 = payload['dFF2']
    valid_mask = (~np.isnan(dFF)) & (~np.isnan(dFF2))
    dFF_valid = dFF[valid_mask]
    dFF2_valid = dFF2[valid_mask]
    plot_violin_with_scatter(
        dFF2_valid, dFF_valid,
        'darkred', 'green',
        xticklabels=['ref.', sensor_label],
        ylabel='ΔF/F',
        title=recname,
        save=True,
        savepath=savepath / f'{recname}_dFF_dFF2_violinplot',
        show=False,
        close=True,
        )

    baseline_dFF = payload['baseline_dFF']
    stim_dFF = payload['stim_dFF']
    baseline2_dFF = payload['baseline2_dFF']
    stim2_dFF = payload['stim2_dFF']

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
        savepath=savepath / f'{recname}_baseline_stim_violinplot',
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
        savepath=savepath / f'{recname}_baseline_stim_ch2_violinplot',
        show=False,
        close=True,
        )


#%% release maps
def save_paired_release_map_figures(
        savepath, recname, ref, ref2, release_map, release_map2):
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
        tifffile.imwrite(
            savepath / f'{recname}_release_map{suffix}.tiff',
            release_rgb
            )
        plt.close(fig)

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

def plot_grabne_release_maps(savepath, recname, proc_data_path):
    release_path = proc_data_path / f'{recname}_release_map.npy'
    release2_path = proc_data_path / f'{recname}_release_map_ch2.npy'
    release_run_path = proc_data_path / f'{recname}_release_map_run.npy'
    release_run2_path = proc_data_path / f'{recname}_release_map_run_ch2.npy'

    ref = np.load(proc_data_path / 'ref_mat_ch1.npy', allow_pickle=True)
    ref2 = np.load(proc_data_path / 'ref_mat_ch2.npy', allow_pickle=True)

    release_map = np.load(release_path, allow_pickle=True)
    release_map2 = np.load(release2_path, allow_pickle=True)
    save_release_map_figure(
        savepath, recname, ref, ref2, release_map,
        'release_map', 'CH1: stim / baseline (mean)'
        )
    save_release_map_tiff(savepath, recname, release_map, 'release_map')
    save_release_map_figure(
        savepath, recname, ref, ref2, release_map2,
        'release_map_ch2', 'CH2: stim / baseline (mean)'
        )
    save_release_map_tiff(savepath, recname, release_map2, 'release_map_ch2')

    if release_run_path.exists() and release_run2_path.exists():
        release_map_run = np.load(release_run_path, allow_pickle=True)
        release_map_run2 = np.load(release_run2_path, allow_pickle=True)
        save_release_map_figure(
            savepath, recname, ref, ref2, release_map_run,
            'release_map_run', 'run-onset / baseline (mean)',
            positive_floor=-.01
            )
        save_release_map_figure(
            savepath, recname, ref, ref2, release_map_run2,
            'release_map_run_ch2', 'stim / baseline (mean)'
            )


#%% session entry
def plot_session_outputs(
        savepath, recname, release_mode,
        plot_stim_summary=True,
        plot_release_maps=True,
        plot_axon_reference=True,
        output_path=None):
    savepath = Path(savepath)
    output_path = savepath if output_path is None else Path(output_path)
    output_path.mkdir(parents=True, exist_ok=True)
    proc_data_path = savepath / 'processed_data'

    if plot_stim_summary:
        payload = np.load(
            proc_data_path / f'{recname}_{PLOT_PAYLOAD_SUFFIX}',
            allow_pickle=True,
            ).item()
        save_aligned_stim_summary(output_path, recname, payload)

    if plot_release_maps:
        if release_mode not in ('paired', 'grabne'):
            raise ValueError(f'unknown release map mode: {release_mode}')
        if release_mode == 'paired':
            release_path = proc_data_path / f'{recname}_release_map.npy'
            release2_path = proc_data_path / f'{recname}_release_map_ch2.npy'
            ref = np.load(proc_data_path / 'ref_mat_ch1.npy', allow_pickle=True)
            ref2 = np.load(proc_data_path / 'ref_mat_ch2.npy', allow_pickle=True)
            release_map = np.load(release_path, allow_pickle=True)
            release_map2 = np.load(release2_path, allow_pickle=True)
            save_paired_release_map_figures(
                output_path, recname, ref, ref2, release_map, release_map2
                )
        else:
            plot_grabne_release_maps(output_path, recname, proc_data_path)

    axon_ref_path = proc_data_path / f'{recname}_ref_mat_1100nm.npy'
    if plot_axon_reference and axon_ref_path.exists():
        reference_axon_only = np.load(axon_ref_path, allow_pickle=True)
        fig, ax = plt.subplots(figsize=(4, 4))
        ax.imshow(
            reference_axon_only,
            aspect='auto', cmap='gist_gray', interpolation='none',
            extent=[0, 512, 512, 0]
            )
        ax.set(xlim=(0, 512), ylim=(0, 512))
        fig.suptitle('ref 1100 nm')
        fig.tight_layout()
        fig.savefig(
            output_path / f'{recname}_ref_1100nm.png',
            dpi=300,
            bbox_inches='tight'
            )
        plt.close(fig)
