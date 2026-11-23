'''
Created on Mon Nov 23 2026

check Up-ACR response variance across split trials

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
SPARSE_PROBABILITIES = (0.5, 0.25, 0.1)


#%% summaries
def median_iqr(values):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if not len(values):
        return np.nan, np.nan, np.nan
    return tuple(np.percentile(values, [25, 50, 75]))


def ratio_summary(values):
    low, median, high = median_iqr(values)
    return median, low, high, int(np.sum(np.isfinite(values)))


def held_ratios(response_change, baseline_change, keep, trials):
    ratios = np.full(response_change.shape[0], np.nan)
    for acr in range(response_change.shape[0]):
        valid = keep[acr, trials]
        response = response_change[acr, trials][valid]
        baseline = baseline_change[acr, trials][valid]
        if len(response) >= 2 and len(baseline) >= 2:
            baseline_sd = np.std(baseline, ddof=1)
            if np.isfinite(baseline_sd) and baseline_sd != 0:
                response_sd = np.std(response, ddof=1)
                ratios[acr] = response_sd / baseline_sd
    return ratios


def sparse_ratio(response, baseline, probability, rng):
    if len(response) < 2 or len(baseline) < 2:
        return np.nan
    baseline_sd = np.std(baseline, ddof=1)
    if not np.isfinite(baseline_sd) or baseline_sd == 0:
        return np.nan
    signal = (rng.random_sample((len(baseline), N_SHUFFLE)) < probability)
    signal = signal * (np.mean(response) / probability)
    predicted = baseline[:, None] + signal
    predicted_sd = np.std(predicted, axis=0, ddof=1)
    return float(np.nanmedian(predicted_sd / baseline_sd))


def direction_summary(response_change, baseline_change, keep, frames, traces,
                      train_trials, held_trials, profile):
    tested = quantify_grid_response(
        traces, frames[train_trials], n_shuffle=N_SHUFFLE, seed=SEED)
    up, _ = classify_up(
        tested['amplitude'],
        tested['effect_size'],
        tested['shuffle_amplitude'],
        np.zeros(len(profile), dtype=bool),
    )
    ratios = held_ratios(response_change, baseline_change, keep, held_trials)
    sparse = {probability: np.full(len(profile), np.nan)
              for probability in SPARSE_PROBABILITIES}
    rng = np.random.RandomState(SEED)
    for acr in np.flatnonzero(up):
        valid = keep[acr, held_trials]
        response = response_change[acr, held_trials][valid]
        baseline = baseline_change[acr, held_trials][valid]
        for probability in SPARSE_PROBABILITIES:
            sparse[probability][acr] = sparse_ratio(
                response, baseline, probability, rng)
    return up, ratios, sparse


#%% one recording
def recording_summary(recname):
    behaviour, profile, traces = load_recording(recname)
    frames = event_frames(behaviour)
    aligned, keep = keep_aligned_trials(align_trials(traces, frames))
    kept_trials = np.flatnonzero(np.any(keep, axis=0))
    if len(kept_trials) < 4:
        raise ValueError(f'{recname}: fewer than four kept trials')

    response_change = (
        np.nanmean(aligned[:, :, 60:90], axis=2)
        - np.nanmean(aligned[:, :, 30:60], axis=2)
    )
    baseline_change = (
        np.nanmean(aligned[:, :, 30:60], axis=2)
        - np.nanmean(aligned[:, :, 0:30], axis=2)
    )
    rng = np.random.RandomState(SEED)
    halves = np.array_split(rng.permutation(kept_trials), 2)

    group_ratios = {'up': defaultdict(list), 'other': defaultdict(list)}
    model_ratios = {probability: defaultdict(list)
                    for probability in SPARSE_PROBABILITIES}
    labels = []
    for train_trials, held_trials in (halves, halves[::-1]):
        up, ratios, sparse = direction_summary(
            response_change, baseline_change, keep, frames, traces,
            train_trials, held_trials, profile)
        labels.append(up)
        for acr, ratio in enumerate(ratios):
            if np.isfinite(ratio):
                label = 'up' if up[acr] else 'other'
                group_ratios[label][acr].append(ratio)
        for probability in SPARSE_PROBABILITIES:
            for acr, ratio in enumerate(sparse[probability]):
                if up[acr] and np.isfinite(ratio):
                    model_ratios[probability][acr].append(ratio)

    grouped = {
        label: np.asarray([np.mean(values) for values in ratios.values()])
        for label, ratios in group_ratios.items()
    }
    models = {
        probability: np.asarray([np.mean(values) for values in ratios.values()])
        for probability, ratios in model_ratios.items()
    }
    discordant = int(np.count_nonzero(labels[0] != labels[1]))
    row = dict(
        rec_id=recname,
        animal=recname.split('-')[0],
        n_trials=len(frames),
        n_kept_trials=len(kept_trials),
        n_valid_acr=len(profile),
        n_discordant_up_calls=discordant,
        n_up_acr=len(grouped['up']),
        n_other_acr=len(grouped['other']),
        n_up_estimates=sum(len(values) for values in group_ratios['up'].values()),
        n_other_estimates=sum(len(values) for values in group_ratios['other'].values()),
    )
    up_median, up_q25, up_q75, up_n = ratio_summary(grouped['up'])
    other_median, other_q25, other_q75, other_n = ratio_summary(grouped['other'])
    row.update(
        ratio_up_median=up_median,
        ratio_up_q25=up_q25,
        ratio_up_q75=up_q75,
        ratio_up_n=up_n,
        ratio_other_median=other_median,
        ratio_other_q25=other_q25,
        ratio_other_q75=other_q75,
        ratio_other_n=other_n,
    )
    for probability, values in models.items():
        label = f'sparse_{str(probability).replace(".", "p")}'
        median, q25, q75, n = ratio_summary(values)
        row.update(
            **{
                f'{label}_up_median': median,
                f'{label}_up_q25': q25,
                f'{label}_up_q75': q75,
                f'{label}_up_n': n,
            }
        )
    return row


#%% report
def append_log(table):
    lines = [
        '',
        '## 23 November 2026: split-trial Up variance',
        '',
        f'I used all {table.n_trials.sum()} valid trials from {len(table)} recordings '
        f'in {table.animal.nunique()} animals. Each recording was split with '
        f'RandomState({SEED}) after the per-ACR full-window NaN rule. Up was called '
        f'on one half with Jingyu\'s test, and the held-out half supplied the '
        f'response SD, baseline SD and their ratio. The halves were swapped, and '
        f'ratios were averaged per ACR within the group assigned by each direction. '
        f'Discordant half-split calls therefore contribute one held-out estimate to '
        f'each corresponding group.',
        '',
    ]
    for label in ('up', 'other'):
        values = table[f'ratio_{label}_median'].to_numpy()
        low, median, high = median_iqr(values)
        lines.append(
            f'{label.capitalize()} ratio: median={median:.9g}, IQR '
            f'[{low:.9g}, {high:.9g}]; n={len(values)} recordings, '
            f'n={table.animal.nunique()} animals.'
        )
    paired = table[['animal', 'ratio_up_median', 'ratio_other_median']].dropna()
    difference = paired.ratio_up_median.to_numpy() - paired.ratio_other_median.to_numpy()
    test = wilcoxon(difference, alternative='two-sided')
    lines.append(
        f'Up against Other: median paired difference={np.median(difference):.9g}, '
        f'IQR [{np.percentile(difference, 25):.9g}, '
        f'{np.percentile(difference, 75):.9g}]; two-sided paired Wilcoxon '
        f'W={test.statistic:.9g}, p={test.pvalue:.12g}; n={len(paired)} '
        f'recordings, n={paired.animal.nunique()} animals.'
    )
    up_against_one = paired.ratio_up_median.to_numpy() - 1
    test = wilcoxon(up_against_one, alternative='two-sided')
    lines.append(
        f'Up against 1: median difference={np.median(up_against_one):.9g}, '
        f'IQR [{np.percentile(up_against_one, 25):.9g}, '
        f'{np.percentile(up_against_one, 75):.9g}]; two-sided Wilcoxon '
        f'W={test.statistic:.9g}, p={test.pvalue:.12g}; n={len(paired)} '
        f'recordings, n={paired.animal.nunique()} animals.'
    )
    lines.append('')
    lines.append('Sparse-response models for Up ACRs, pooled over recording medians:')
    for probability in SPARSE_PROBABILITIES:
        label = f'sparse_{str(probability).replace(".", "p")}_up_median'
        values = table[label].to_numpy()
        low, median, high = median_iqr(values)
        lines.append(
            f'p={probability:g}: expected ratio median={median:.9g}, '
            f'IQR [{low:.9g}, {high:.9g}]; n={len(values)} recordings, '
            f'n={table.animal.nunique()} animals.'
        )
    lines += [
        '',
        'The sparse models add mean(response) / p to a random fraction p of '
        'held-out baseline-window values and add zero to the rest, with 1,000 '
        'random draws per Up ACR and the same held-out mean used in the '
        'corresponding observed ratio.',
        '',
        'Per-recording results: outputs/da_hotspots/up_variance_split.csv.',
        '',
    ]
    with (OUTPUT_ROOT / 'PROGRESS.md').open('a', encoding='utf-8') as stream:
        stream.write('\n'.join(lines))


#%% run
def main():
    output = OUTPUT_ROOT / 'up_variance_split.csv'
    if output.exists():
        raise FileExistsError(output)
    rows = []
    for recname in rec_list.recDLightDbh:
        row = recording_summary(recname)
        rows.append(row)
        print(f'{recname}: Up ratio={row["ratio_up_median"]:.6g}, '
              f'Other ratio={row["ratio_other_median"]:.6g}', flush=True)
    table = pd.DataFrame(rows)
    table.to_csv(output, index=False)
    append_log(table)


if __name__ == '__main__':
    with threadpool_limits(limits=1):
        main()
