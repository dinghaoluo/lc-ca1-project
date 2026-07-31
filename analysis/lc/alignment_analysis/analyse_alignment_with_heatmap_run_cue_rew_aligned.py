# -*- coding: utf-8 -*-
'''
Created on Sun Jul 20 17:53:46 2025
Modified on 5 Dec 2025

plot population heatmap, but aligned to cue and rew
modified to include statistics

@author: Dinghao Luo
'''

#%% imports
from pathlib import Path
import sys

repo_root = Path(__file__).resolve().parents[3]

if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
utils_path = repo_root / 'utils'
if str(utils_path) not in sys.path:
    sys.path.insert(0, str(utils_path))

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pickle
from tqdm import tqdm
from scipy.stats import wilcoxon, ttest_rel, sem

import rec_list
paths = rec_list.pathLC

from common_functions import normalise, mpl_formatting, colour_putative, colour_tagged
from console_formatting import print_session, print_statistics_section
import project_paths as pp
mpl_formatting()


#%% paths
LC_stem = pp.LC_EPHYS_STEM
all_sess_stem = LC_stem / 'all_sessions'
beh_stem = pp.behaviour_experiment_stem('LC')
population_map_stem = LC_stem / 'population_maps'
population_map_stem.mkdir(parents=True, exist_ok=True)

CENTER = 1250
WINDOW = 313
BUMP = 0.02


#%% statistics and heatmaps

def _print_paired_stats(x, y, label_x, label_y):
    '''
    print mean +/- sem, median [iqr], wilcoxon and paired t-test.
    '''
    x = np.asarray(x, float)
    y = np.asarray(y, float)

    mask = np.isfinite(x) & np.isfinite(y)
    x = x[mask]
    y = y[mask]

    print(f'\n{label_x} vs {label_y}')
    print(f'n = {len(x)}')
    print(f'{label_x}: mean +/- sem = {np.mean(x):.3f} +/- {sem(x):.3f}')
    print(f'{label_y}: mean +/- sem = {np.mean(y):.3f} +/- {sem(y):.3f}')

    mx, q1x, q3x = np.median(x), *np.percentile(x, [25, 75])
    my, q1y, q3y = np.median(y), *np.percentile(y, [25, 75])
    print(f'{label_x}: median [iqr] = {mx:.3f} [{q1x:.3f}, {q3x:.3f}]')
    print(f'{label_y}: median [iqr] = {my:.3f} [{q1y:.3f}, {q3y:.3f}]')

    wstat, wp = wilcoxon(x, y)
    tstat, tp = ttest_rel(x, y)

    print(f'wilcoxon: W = {wstat:.3f}, p = {wp:.3e}')
    print(f'paired t: t = {tstat:.3f}, p = {tp:.3e}')

def _plot_heatmap(data, title, xlabel, save_stub, cmap, ytick_values):
    fig, ax = plt.subplots(figsize=(2.6, 2.1))
    ax.set(xlabel=xlabel, ylabel='Cell #')
    ax.set_aspect('equal')
    fig.suptitle(title)

    image = ax.imshow(
        data,
        aspect='auto',
        cmap=cmap,
        interpolation='none',
        extent=[-1, 4, 1, len(data)],
    )
    plt.colorbar(image, shrink=.5, ticks=[0, 1], label='Norm. spike rate')
    ax.set(yticks=ytick_values)

    for ext in ['.png', '.pdf']:
        fig.savefig(population_map_stem / f'{save_stub}{ext}', dpi=300, bbox_inches='tight')

#%% main

cell_profiles = pd.read_pickle(LC_stem / 'LC_all_cell_profiles.pkl')
tag_list = [clu for clu in cell_profiles.index if cell_profiles['identity'][clu] == 'tagged']
put_list = [clu for clu in cell_profiles.index if cell_profiles['identity'][clu] == 'putative']
peak_list = [clu for clu in cell_profiles.index if cell_profiles['run_onset_peak'][clu]]

all_tagged_run = []
all_tagged_cue = []
all_tagged_rew = []
all_putative_run = []
all_putative_cue = []
all_putative_rew = []
all_pooled_run = []

sess_tagged_peaks_run = {}
sess_tagged_peaks_cue = {}
sess_tagged_peaks_rew = {}
sess_put_peaks_run = {}
sess_put_peaks_cue = {}
sess_put_peaks_rew = {}
peak_run_peak_time = []

for path in paths:
    recname = Path(path).name
    print_session(recname)

    with open(beh_stem / f'{recname}.pkl', 'rb') as f:
        beh = pickle.load(f)
    stim_conds = [trial[15] for trial in beh['trial_statements']][1:]
    try:
        stim_start = stim_conds.index('2')
    except ValueError:
        stim_start = len(stim_conds)

    sess_path = all_sess_stem / recname
    trains_run = np.load(sess_path / f'{recname}_all_trains_run.npy', allow_pickle=True).item()
    trains_cue = np.load(sess_path / f'{recname}_all_trains_cue.npy', allow_pickle=True).item()
    trains_rew = np.load(sess_path / f'{recname}_all_trains_rew.npy', allow_pickle=True).item()

    for clu in tqdm(trains_cue, total=len(trains_cue)):
        if clu not in tag_list and clu not in put_list:
            continue

        mean_run = np.mean(trains_run[clu][:stim_start, 3750 - 1250:3750 + 1250 * 4], axis=0)
        mean_cue = np.mean(trains_cue[clu][:stim_start, 3750 - 1250:3750 + 1250 * 4], axis=0)
        mean_rew = np.mean(trains_rew[clu][:stim_start, 3750 - 1250:3750 + 1250 * 4], axis=0)

        all_pooled_run.append(mean_run)

        if clu in tag_list:
            sess_tagged_peaks_run.setdefault(recname, []).append(np.argmax(mean_run))
            sess_tagged_peaks_cue.setdefault(recname, []).append(np.argmax(mean_cue))
            sess_tagged_peaks_rew.setdefault(recname, []).append(np.argmax(mean_rew))
            all_tagged_run.append(mean_run)
            all_tagged_cue.append(mean_cue)
            all_tagged_rew.append(mean_rew)

        if clu in put_list:
            sess_put_peaks_run.setdefault(recname, []).append(np.argmax(mean_run))
            sess_put_peaks_cue.setdefault(recname, []).append(np.argmax(mean_cue))
            sess_put_peaks_rew.setdefault(recname, []).append(np.argmax(mean_rew))
            all_putative_run.append(mean_run)
            all_putative_cue.append(mean_cue)
            all_putative_rew.append(mean_rew)

        if clu in peak_list:
            peak_run_peak_time.append(np.argmax(mean_run[:2500]))

pooled_run_argmax = [np.argmax(clu) for clu in all_pooled_run]
pooled_run_sorted_idx = np.argsort(pooled_run_argmax)
all_pooled_run_sorted = np.array([normalise(all_pooled_run[clu]) for clu in pooled_run_sorted_idx])

tagged_run_argmax = [np.argmax(clu) for clu in all_tagged_run]
tagged_run_sorted_idx = np.argsort(tagged_run_argmax)
all_tagged_run_sorted = np.array([normalise(all_tagged_run[clu]) for clu in tagged_run_sorted_idx])

putative_run_argmax = [np.argmax(clu) for clu in all_putative_run]
putative_run_sorted_idx = np.argsort(putative_run_argmax)
all_putative_run_sorted = np.array([normalise(all_putative_run[clu]) for clu in putative_run_sorted_idx])

tagged_cue_argmax = [np.argmax(clu) for clu in all_tagged_cue]
tagged_cue_sorted_idx = np.argsort(tagged_cue_argmax)
all_tagged_cue_sorted = np.array([normalise(all_tagged_cue[clu]) for clu in tagged_cue_sorted_idx])

putative_cue_argmax = [np.argmax(clu) for clu in all_putative_cue]
putative_cue_sorted_idx = np.argsort(putative_cue_argmax)
all_putative_cue_sorted = np.array([normalise(all_putative_cue[clu]) for clu in putative_cue_sorted_idx])

tagged_rew_argmax = [np.argmax(clu) for clu in all_tagged_rew]
tagged_rew_sorted_idx = np.argsort(tagged_rew_argmax)
all_tagged_rew_sorted = np.array([normalise(all_tagged_rew[clu]) for clu in tagged_rew_sorted_idx])

putative_rew_argmax = [np.argmax(clu) for clu in all_putative_rew]
putative_rew_sorted_idx = np.argsort(putative_rew_argmax)
all_putative_rew_sorted = np.array([normalise(all_putative_rew[clu]) for clu in putative_rew_sorted_idx])

for cmap, suffix in [('Greys', ''), ('viridis', '_viridis')]:
    _plot_heatmap(all_pooled_run_sorted, 'Pooled Dbh+ cells', 'Time from run onset (s)', f'pooled_run_aligned{suffix}', cmap, [1, 200])
    _plot_heatmap(all_tagged_run_sorted, 'Tagged Dbh+ cells', 'Time from run onset (s)', f'tagged_run_aligned{suffix}', cmap, [1, 80])
    _plot_heatmap(all_putative_run_sorted, 'Putative Dbh+ cells', 'Time from run onset (s)', f'putative_run_aligned{suffix}', cmap, [1, 100, 200])
    _plot_heatmap(all_tagged_cue_sorted, 'Tagged Dbh+ cells', 'Time from cue (s)', f'tagged_cue_aligned{suffix}', cmap, [1, 80])
    _plot_heatmap(all_putative_cue_sorted, 'Putative Dbh+ cells', 'Time from cue (s)', f'putative_cue_aligned{suffix}', cmap, [1, 100, 200])
    _plot_heatmap(all_tagged_rew_sorted, 'Tagged Dbh+ cells', 'Time from reward (s)', f'tagged_rew_aligned{suffix}', cmap, [1, 80])
    _plot_heatmap(all_putative_rew_sorted, 'Putative Dbh+ cells', 'Time from reward (s)', f'putative_rew_aligned{suffix}', cmap, [1, 100, 200])

peaks_tagged_run = np.argmax(all_tagged_run_sorted, axis=1)
peaks_tagged_cue = np.argmax(all_tagged_cue_sorted, axis=1)
peaks_tagged_rew = np.argmax(all_tagged_rew_sorted, axis=1)
peaks_put_run = np.argmax(all_putative_run_sorted, axis=1)
peaks_put_cue = np.argmax(all_putative_cue_sorted, axis=1)
peaks_put_rew = np.argmax(all_putative_rew_sorted, axis=1)
proportions = [
    [np.mean((peaks_tagged_run >= CENTER - WINDOW) & (peaks_tagged_run <= CENTER + WINDOW)),
     np.mean((peaks_put_run >= CENTER - WINDOW) & (peaks_put_run <= CENTER + WINDOW))],
    [np.mean((peaks_tagged_cue >= CENTER - WINDOW) & (peaks_tagged_cue <= CENTER + WINDOW)),
     np.mean((peaks_put_cue >= CENTER - WINDOW) & (peaks_put_cue <= CENTER + WINDOW))],
    [np.mean((peaks_tagged_rew >= CENTER - WINDOW) & (peaks_tagged_rew <= CENTER + WINDOW)),
     np.mean((peaks_put_rew >= CENTER - WINDOW) & (peaks_put_rew <= CENTER + WINDOW))],
]

sess_p_tagged_run = []
sess_p_tagged_cue = []
sess_p_tagged_rew = []
sess_p_put_run = []
sess_p_put_cue = []
sess_p_put_rew = []
for path in paths:
    recname = Path(path).name
    peaks = np.array(sess_tagged_peaks_run.get(recname, []))
    if peaks.size > 0:
        sess_p_tagged_run.append(np.mean((peaks >= CENTER - WINDOW) & (peaks <= CENTER + WINDOW)))
    peaks = np.array(sess_tagged_peaks_cue.get(recname, []))
    if peaks.size > 0:
        sess_p_tagged_cue.append(np.mean((peaks >= CENTER - WINDOW) & (peaks <= CENTER + WINDOW)))
    peaks = np.array(sess_tagged_peaks_rew.get(recname, []))
    if peaks.size > 0:
        sess_p_tagged_rew.append(np.mean((peaks >= CENTER - WINDOW) & (peaks <= CENTER + WINDOW)))
    peaks = np.array(sess_put_peaks_run.get(recname, []))
    if peaks.size > 0:
        sess_p_put_run.append(np.mean((peaks >= CENTER - WINDOW) & (peaks <= CENTER + WINDOW)))
    peaks = np.array(sess_put_peaks_cue.get(recname, []))
    if peaks.size > 0:
        sess_p_put_cue.append(np.mean((peaks >= CENTER - WINDOW) & (peaks <= CENTER + WINDOW)))
    peaks = np.array(sess_put_peaks_rew.get(recname, []))
    if peaks.size > 0:
        sess_p_put_rew.append(np.mean((peaks >= CENTER - WINDOW) & (peaks <= CENTER + WINDOW)))

labels = ['run', 'cue', 'rew']
x = np.arange(len(labels))
height = .35
fig, ax = plt.subplots(figsize=(2, 3))
for i in range(len(labels)):
    ax.barh(i - height / 2, proportions[i][0], height,
            label='tagged' if i == 0 else '', color=colour_tagged)
    ax.barh(i + height / 2, proportions[i][1], height,
            label='putative' if i == 0 else '', color=colour_putative)

ytag_run, yput_run = -height / 2, height / 2
ytag_cue, yput_cue = 1 - height / 2, 1 + height / 2
ytag_rew, yput_rew = 2 - height / 2, 2 + height / 2
if len(sess_p_tagged_run) > 0:
    vals = np.array(sess_p_tagged_run, float)
    vals[vals == 0] = BUMP
    ax.scatter(vals, np.full(len(vals), ytag_run), s=8, color=colour_tagged, edgecolors='k')
if len(sess_p_put_run) > 0:
    vals = np.array(sess_p_put_run, float)
    vals[vals == 0] = BUMP
    ax.scatter(vals, np.full(len(vals), yput_run), s=8, color=colour_putative, edgecolors='k')
if len(sess_p_tagged_cue) > 0:
    vals = np.array(sess_p_tagged_cue, float)
    vals[vals == 0] = BUMP
    ax.scatter(vals, np.full(len(vals), ytag_cue), s=8, color=colour_tagged, edgecolors='k')
if len(sess_p_put_cue) > 0:
    vals = np.array(sess_p_put_cue, float)
    vals[vals == 0] = BUMP
    ax.scatter(vals, np.full(len(vals), yput_cue), s=8, color=colour_putative, edgecolors='k')
if len(sess_p_tagged_rew) > 0:
    vals = np.array(sess_p_tagged_rew, float)
    vals[vals == 0] = BUMP
    ax.scatter(vals, np.full(len(vals), ytag_rew), s=8, color=colour_tagged, edgecolors='k')
if len(sess_p_put_rew) > 0:
    vals = np.array(sess_p_put_rew, float)
    vals[vals == 0] = BUMP
    ax.scatter(vals, np.full(len(vals), yput_rew), s=8, color=colour_putative, edgecolors='k')

ax.set_yticks(x)
ax.set_yticklabels(labels)
ax.set_xlabel('prop. with peak +/-0.25 s')
ax.set_title('Peak proximity to alignment')
ax.legend(frameon=False, loc='upper right')
for spine in ['top', 'right', 'bottom']:
    ax.spines[spine].set_visible(False)
fig.tight_layout()
for ext in ['.png', '.pdf']:
    fig.savefig(population_map_stem / f'peak_proximity_bar{ext}', dpi=300, bbox_inches='tight')

print_statistics_section()
print('tagged LC (session-level):')
_print_paired_stats(sess_p_tagged_run, sess_p_tagged_cue, 'run', 'cue')
_print_paired_stats(sess_p_tagged_run, sess_p_tagged_rew, 'run', 'reward')
_print_paired_stats(sess_p_tagged_cue, sess_p_tagged_rew, 'cue', 'reward')
print('\nputative LC (session-level):')
_print_paired_stats(sess_p_put_run, sess_p_put_cue, 'run', 'cue')
_print_paired_stats(sess_p_put_run, sess_p_put_rew, 'run', 'reward')
_print_paired_stats(sess_p_put_cue, sess_p_put_rew, 'cue', 'reward')

_ = peak_run_peak_time
