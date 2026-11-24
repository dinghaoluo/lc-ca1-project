'''
Created on Mon Nov 24 2026

check split-trial Up and baseline response medians

@author: Dinghao Luo
'''

#%% imports
import sys
from collections import defaultdict
from pathlib import Path

repo_root = Path(__file__).resolve().parents[3]
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))

import project_paths as pp
import rec_list

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon
from threadpoolctl import threadpool_limits

from analyse_acr_consistency_aligned import load_recording, event_frames
from grid_response_test import (align_trials, classify_up, keep_aligned_trials,
                                quantify_grid_response)

OUTPUT_ROOT = pp.DLIGHT_HOTSPOTS_OUTPUT_ROOT
N_SHUFFLE = 1000
SEED = 42


#%% summaries
def median_iqr(values):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    return tuple(np.percentile(values, [25, 50, 75])) if len(values) else (np.nan, np.nan, np.nan)


def group_summary(values):
    values = {
        acr: np.mean(acr_values)
        for acr, acr_values in values.items()
        if len(acr_values)
    }
    finite = np.asarray(list(values.values()), dtype=float)
    low, median, high = median_iqr(finite)
    return median, low, high, int(len(finite)), values


#%% one recording
def recording_summary(recname):
    behaviour, profile, traces = load_recording(recname)
    frames = event_frames(behaviour)
    aligned, keep = keep_aligned_trials(align_trials(traces, frames))
    kept_trials = np.flatnonzero(np.any(keep, axis=0))
    if len(kept_trials) < 4:
        raise ValueError(f'{recname}: fewer than four kept trials')

    response = (
        np.nanmean(aligned[:, :, 60:90], axis=2)
        - np.nanmean(aligned[:, :, 30:60], axis=2)
    )
    baseline = (
        np.nanmean(aligned[:, :, 30:60], axis=2)
        - np.nanmean(aligned[:, :, 0:30], axis=2)
    )
    halves = np.array_split(np.random.RandomState(SEED).permutation(kept_trials), 2)
    group_values = {
        'up': {'response': defaultdict(list), 'baseline': defaultdict(list)},
        'other': {'response': defaultdict(list), 'baseline': defaultdict(list)},
    }
    labels = []
    for train_trials, held_trials in (halves, halves[::-1]):
        tested = quantify_grid_response(
            traces, frames[train_trials], n_shuffle=N_SHUFFLE, seed=SEED)
        up, _ = classify_up(
            tested['amplitude'], tested['effect_size'],
            tested['shuffle_amplitude'],
            np.zeros(len(profile), dtype=bool),
        )
        labels.append(up)
        for acr in range(len(profile)):
            valid = keep[acr, held_trials]
            response_values = response[acr, held_trials][valid]
            baseline_values = baseline[acr, held_trials][valid]
            label = 'up' if up[acr] else 'other'
            if len(response_values):
                group_values[label]['response'][acr].append(float(np.median(response_values)))
            if len(baseline_values):
                group_values[label]['baseline'][acr].append(float(np.median(baseline_values)))

    row = dict(
        rec_id=recname,
        animal=recname.split('-')[0],
        n_trials=len(frames),
        n_kept_trials=len(kept_trials),
        n_valid_acr=len(profile),
        n_discordant_up_calls=int(np.count_nonzero(labels[0] != labels[1])),
    )
    for label in ('up', 'other'):
        row[f'n_{label}_acr'] = len(set(group_values[label]['baseline']))
        for measure in ('response', 'baseline'):
            median, low, high, n, _ = group_summary(group_values[label][measure])
            row.update(
                **{
                    f'{measure}_{label}_median': median,
                    f'{measure}_{label}_q25': low,
                    f'{measure}_{label}_q75': high,
                    f'{measure}_{label}_n': n,
                }
            )
    return row


#%% report
def append_log(table):
    lines = [
        '',
        '## 24 November 2026: split-trial Up and baseline medians',
        '',
        f'I used all {table.n_trials.sum()} valid trials from '
        f'{len(table)} recordings in {table.animal.nunique()} animals, with '
        f'Jingyu\'s per-ACR full-window NaN rule. Each recording was split with '
        f'RandomState({SEED}), Up was called on one half, and response and baseline '
        f'medians were computed on the held-out half. The halves were swapped, and '
        f'per-ACR estimates were averaged within the group assigned by each direction.',
        '',
    ]
    for measure, label in (('response', 'response'), ('baseline', 'baseline')):
        up = table[f'{measure}_up_median'].to_numpy()
        other = table[f'{measure}_other_median'].to_numpy()
        low, median, high = median_iqr(up)
        other_low, other_median, other_high = median_iqr(other)
        lines.append(
            f'{label}, Up: median={median:.9g}, IQR [{low:.9g}, {high:.9g}]; '
            f'Other: median={other_median:.9g}, IQR '
            f'[{other_low:.9g}, {other_high:.9g}]; n={len(table)} recordings, '
            f'n={table.animal.nunique()} animals.'
        )
    paired = table[['animal', 'baseline_up_median', 'baseline_other_median']].dropna()
    up = paired.baseline_up_median.to_numpy()
    other = paired.baseline_other_median.to_numpy()
    test = wilcoxon(up, alternative='two-sided')
    lines.append(
        f'Baseline Up against 0: median={np.median(up):.9g}, IQR '
        f'[{np.percentile(up, 25):.9g}, {np.percentile(up, 75):.9g}]; '
        f'two-sided Wilcoxon W={test.statistic:.9g}, p={test.pvalue:.12g}; '
        f'n={len(up)} recordings, n={paired.animal.nunique()} animals.'
    )
    difference = up - other
    test = wilcoxon(difference, alternative='two-sided')
    lines.append(
        f'Baseline Up against Other: median paired difference={np.median(difference):.9g}, '
        f'IQR [{np.percentile(difference, 25):.9g}, '
        f'{np.percentile(difference, 75):.9g}]; two-sided paired Wilcoxon '
        f'W={test.statistic:.9g}, p={test.pvalue:.12g}; n={len(paired)} recordings, '
        f'n={paired.animal.nunique()} animals.'
    )
    lines += ['', 'Per-recording results: outputs/da_hotspots/up_baseline_split.csv.', '']
    with (OUTPUT_ROOT / 'PROGRESS.md').open('a', encoding='utf-8') as stream:
        stream.write('\n'.join(lines))


#%% run
def main():
    output = OUTPUT_ROOT / 'up_baseline_split.csv'
    if output.exists():
        raise FileExistsError(output)
    rows = []
    for recname in rec_list.recDLightDbh:
        row = recording_summary(recname)
        rows.append(row)
        print(f'{recname}: Up baseline={row["baseline_up_median"]:.6g}, '
              f'Other baseline={row["baseline_other_median"]:.6g}', flush=True)
    table = pd.DataFrame(rows)
    table.to_csv(output, index=False)
    append_log(table)


if __name__ == '__main__':
    with threadpool_limits(limits=1):
        main()
