# -*- coding: utf-8 -*-
'''
Created on Thu Nov 12 10:20:36 2026

single-trial dLight amplitude maps and aligned trace examples

@author: Dinghao Luo
'''

#%% imports
import sys
from pathlib import Path

repo_root = Path(__file__).resolve().parents[3]
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))

import project_paths as pp
import rec_list

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from common_functions import mpl_formatting
mpl_formatting()
from matplotlib.patches import Rectangle
import numpy as np
import pandas as pd

from analyse_acr_consistency_aligned import (RAW_ROOT, recnames,
                                              load_recording, event_frames)
from grid_response_test import align_trials, keep_aligned_trials, trial_amplitudes

repo_root   = Path(__file__).resolve().parents[3]
FIGURE_ROOT = pp.DLIGHT_HOTSPOTS_FIGURE_ROOT


#%% selection
def choose_recordings():
    rows = []
    for recname in list(recnames):
        path = (RAW_ROOT / 'processed_dataframe_grid_free_dilation'
                / f'{recname}_profile_combined_dilation=0_pre(-1, 0)_post(0, 1)'
                '_ES=0.05_shuff95_0621.parquet')
        profile = pd.read_parquet(path)
        valid = (profile['dlight_valid'].astype(bool)
                 & profile['red_valid'].astype(bool)
                 & ~profile['edge'].astype(bool))
        profile = profile.loc[valid]
        if len(profile) >= 50:
            rows.append(dict(rec_id=recname, n_valid=len(profile),
                             pct_up=100 * profile['Up'].astype(bool).mean()))
    eligible = pd.DataFrame(rows)
    median = eligible.pct_up.median()
    chosen = [(label, eligible.loc[index]) for label, index in (
        ('highest', eligible.pct_up.idxmax()),
        ('closest to median', (eligible.pct_up - median).abs().idxmin()),
        ('lowest', eligible.pct_up.idxmin()))]
    return eligible, median, chosen


#%% figures
def plot_recording(recname):
    # the loader applies postprocess_dff to the full grid before subsetting
    behaviour, profile, traces = load_recording(recname)
    frames = event_frames(behaviour)
    run = np.asarray(behaviour['run_onset_frames'], dtype=float)
    reward = np.asarray(behaviour['reward_times'], dtype=float)
    non_stop = np.asarray(behaviour['non_stop_trials'])
    non_fullstop = np.asarray(behaviour['non_fullstop_trials'])
    n = min(len(run), len(reward), len(non_stop), len(non_fullstop))
    valid = (np.isfinite(run[:n]) & (run[:n] != -1)
             & np.isfinite(reward[:n]) & (non_stop[:n] == 0)
             & (non_fullstop[:n] == 0))
    trial_numbers = np.flatnonzero(valid) + 1
    assert np.array_equal(run[:n][valid].astype(int), frames)
    if len(frames) < 5:
        raise ValueError(f'{recname}: fewer than five valid trials')
    first = (len(frames) - 5) // 2
    examples = np.arange(first, first + 5)
    aligned, keep = keep_aligned_trials(align_trials(traces, frames))
    amplitude, baseline, _ = trial_amplitudes(aligned)
    average_amplitude = np.nanmean(amplitude, axis=1)
    centred = aligned - np.nanmean(baseline, axis=2)[:, :, None]
    average_trace = np.nanmean(centred, axis=1)
    order = np.argsort(-average_amplitude, kind='stable')
    up = profile['Up'].to_numpy(dtype=bool)
    coords = np.asarray(profile['roi_id'].tolist(), dtype=int)
    maps = [average_amplitude] + [amplitude[:, trial] for trial in examples]
    heatmaps = [average_trace[order]] + [centred[order, trial, :] for trial in examples]
    values = np.concatenate(maps)
    map_limit = np.percentile(np.abs(values[np.isfinite(values)]), 99)
    values = np.concatenate([heat.ravel() for heat in heatmaps])
    heat_limit = np.percentile(np.abs(values[np.isfinite(values)]), 99)
    cmap = plt.get_cmap('RdBu_r').copy()
    cmap.set_bad('0.7')
    titles = ['Trial average'] + [f'Trial {trial_numbers[trial]}' for trial in examples]
    fig = plt.figure(figsize=(17, 7.2), layout='constrained')
    grid = fig.add_gridspec(2, 7, width_ratios=[1, 1, 1, 1, 1, 1, 0.06],
                           height_ratios=[1, 1.25])
    for column, (values, heat, title) in enumerate(zip(maps, heatmaps, titles)):
        ax = fig.add_subplot(grid[0, column])
        spatial = np.full((32, 32), np.nan)
        spatial[coords[:, 0], coords[:, 1]] = values
        image_map = ax.imshow(spatial, cmap=cmap, vmin=-map_limit, vmax=map_limit,
                              origin='upper', interpolation='nearest')
        for y, x in coords[up]:
            ax.add_patch(Rectangle((x - 0.5, y - 0.5), 1, 1,
                                   facecolor='none', edgecolor='black', linewidth=0.65))
        ax.set_title(title, fontsize=10)
        ax.set_xlabel('Grid column')
        if column == 0:
            ax.set_ylabel('Grid row')
        ax.set_xticks([0, 15, 31])
        ax.set_yticks([0, 15, 31])
        ax = fig.add_subplot(grid[1, column])
        image_heat = ax.imshow(heat, cmap=cmap, vmin=-heat_limit, vmax=heat_limit,
                               origin='upper', aspect='auto', interpolation='nearest',
                               extent=(-2, 4, len(profile) - 0.5, -0.5))
        ax.axvline(0, color='black', linewidth=0.8)
        rows = np.flatnonzero(up[order])
        ax.scatter(np.full(len(rows), -0.025), rows, transform=ax.get_yaxis_transform(),
                   marker='>', s=8, color='black', linewidths=0, clip_on=False)
        ax.set_title(title, fontsize=10)
        ax.set_xlabel('Time from run onset (s)')
        if column == 0:
            ax.set_ylabel('ACR row (sorted by mean amplitude)')
        ax.set_xticks([-2, 0, 2, 4])
        ticks = np.unique([0, len(profile) // 2, len(profile) - 1])
        ax.set_yticks(ticks)
        ax.set_yticklabels(ticks + 1)
    fig.colorbar(image_map, cax=fig.add_subplot(grid[0, 6]),
                 label='Response amplitude (dF/F)')
    fig.colorbar(image_heat, cax=fig.add_subplot(grid[1, 6]),
                 label='Baseline-subtracted dF/F')
    path = FIGURE_ROOT / f'single_trial_examples_{recname}.png'
    if path.exists():
        raise FileExistsError(f'Existing figure left unchanged: {path}')
    fig.savefig(path, dpi=200, facecolor='white')
    plt.close(fig)
    trials = ', '.join(str(number) for number in trial_numbers[examples])
    positions = ', '.join(str(trial + 1) for trial in examples)
    missing = ', '.join(str(int(np.sum(~keep[:, trial]))) for trial in examples)
    return (f'{recname}: {len(profile)} valid ACRs, {100 * up.mean():.9g}% saved Up; '
            f'{len(frames)} valid trials; example behaviour trial numbers (one-based) '
            f'{trials}, at valid-list positions {positions}. Missing ACR windows in these '
            f'five trials: {missing}; ACRs with no kept trial for the average: '
            f'{int(np.sum(keep.sum(axis=1) == 0))}. Symmetric colour limits: amplitude '
            f'+/-{map_limit:.9g} dF/F, heatmap +/-{heat_limit:.9g} dF/F. '
            f'Figure: figures/dlight_hotspots/{path.name}.')


#%% run
def main():
    eligible, median, chosen = choose_recordings()
    print(f'Eligible recordings: {len(eligible)}; median saved Up={median:.9g}%', flush=True)
    for label, row in chosen:
        result = plot_recording(row.rec_id)
        print(f'{label}: {result}', flush=True)


if __name__ == '__main__':
    main()
