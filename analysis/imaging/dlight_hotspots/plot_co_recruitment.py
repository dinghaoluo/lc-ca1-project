# -*- coding: utf-8 -*-
'''
Created on Fri Nov 13 18:21:14 2026

plot trialwise Up versus other co-recruitment and recording correlations

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
original_plot_style = plt.rcParams.copy()
mpl_formatting()
import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

from check_smoothed_template_and_co_recruitment import (
    load_recording, event_frames, align_trials, keep_aligned_trials,
    trial_amplitudes, group_mean, correlation, OUTPUT_ROOT)
from plot_single_trial_examples import choose_recordings

repo_root = Path(__file__).resolve().parents[3]
FIGURE_PATH = repo_root / 'figures' / 'dlight_hotspots' / 'co_recruitment.png'


#%% trial means
def trial_means(recname):
    behaviour, profile, traces = load_recording(recname)
    frames = event_frames(behaviour)
    aligned, _ = keep_aligned_trials(align_trials(traces, frames))
    response, _, _ = trial_amplitudes(aligned)
    baseline = aligned[:, :, 30:60].mean(axis=2) - aligned[:, :, 0:30].mean(axis=2)
    up = profile['Up'].to_numpy(dtype=bool)
    return [(group_mean(values[up]), group_mean(values[~up]))
            for values in (response, baseline)]


#%% figure
def main():
    saved = pd.read_csv(OUTPUT_ROOT / 'co_recruitment.csv')
    assert len(saved) == 57 and saved.rec_id.nunique() == 57
    _, _, chosen = choose_recordings()
    fig = plt.figure(figsize=(14, 9.5), layout='constrained')
    grid = fig.add_gridspec(3, 3, width_ratios=[1, 1, 1.05])
    for row, (_, recording) in enumerate(chosen):
        recname = recording.rec_id
        pairs = trial_means(recname)
        all_values = np.concatenate([values for pair in pairs for values in pair])
        all_values = all_values[np.isfinite(all_values)]
        low, high = all_values.min(), all_values.max()
        margin = 0.06 * (high - low) if high > low else 0.1
        limits = (low - margin, high + margin)
        for column, (label, pair) in enumerate(zip(('response', 'baseline'), pairs)):
            up_mean, other_mean = pair
            finite = np.isfinite(up_mean) & np.isfinite(other_mean)
            x, y = up_mean[finite], other_mean[finite]
            r, n, _ = correlation(up_mean, other_mean)
            previous = saved.loc[saved.rec_id.eq(recname)].iloc[0]
            key = 'resp' if label == 'response' else 'base'
            assert np.isclose(r, previous[f'r_{key}'])
            assert n == previous[f'n_{key}_trials']
            ax = fig.add_subplot(grid[row, column])
            ax.scatter(x, y, s=12, color='0.35', alpha=0.65, linewidths=0,
                       label='Trial')
            if len(x) >= 2 and np.ptp(x) > 0:
                slope, intercept = np.linalg.lstsq(
                    np.column_stack((x, np.ones(len(x)))), y, rcond=None)[0]
                line = np.array([x.min(), x.max()])
                ax.plot(line, slope * line + intercept, color='#228b22',
                        linewidth=1.5, label='OLS fit')
            ax.set(xlim=limits, ylim=limits,
                   xlabel='Up ACRs, mean dF/F change',
                   ylabel='Other ACRs, mean dF/F change')
            ax.set_title(f'{recname}, {label} window\nr = {r:.3f}', fontsize=10)
            ax.spines[['top', 'right']].set_visible(False)
            if row == 0 and column == 0:
                ax.legend(frameon=False, fontsize=8)
        print(f'{recname}: response and baseline scatter panels, '
              f'n={sum(np.isfinite(pairs[0][0]) & np.isfinite(pairs[0][1]))} trials',
              flush=True)
    ax = fig.add_subplot(grid[:, 2])
    finite = np.isfinite(saved.r_base) & np.isfinite(saved.r_resp)
    assert finite.sum() == 57
    ax.scatter(saved.loc[finite, 'r_base'], saved.loc[finite, 'r_resp'],
               s=24, color='0.35', alpha=0.8, linewidths=0, label='Recording')
    ax.plot([-1, 1], [-1, 1], color='0.5', linestyle='--', linewidth=1,
            label='Identity')
    ax.set(xlim=(-1, 1), ylim=(-1, 1),
           xlabel='r, baseline window', ylabel='r, response window',
           xticks=[-1, -0.5, 0, 0.5, 1], yticks=[-1, -0.5, 0, 0.5, 1])
    ax.set_aspect('equal', adjustable='box')
    ax.spines[['top', 'right']].set_visible(False)
    ax.legend(frameon=False, fontsize=8)
    fig.savefig(FIGURE_PATH, dpi=200, facecolor='white')
    plt.close(fig)
    print(f'Saved {FIGURE_PATH}; final panel n={int(finite.sum())} recordings.', flush=True)


#%% run
if __name__ == '__main__':
    with threadpool_limits(limits=1):
        with plt.rc_context(original_plot_style):
            main()
