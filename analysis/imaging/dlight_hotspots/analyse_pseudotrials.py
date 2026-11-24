'''
Created on Mon Nov 24 2026

analyse aligned response pseudotrials across trial-count scales

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

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, wilcoxon
from threadpoolctl import threadpool_limits

from analyse_acr_consistency_aligned import jaccard_summary, load_recording, event_frames
from grid_response_test import pooled_std, quantify_grid_response

OUTPUT_ROOT = pp.DLIGHT_HOTSPOTS_OUTPUT_ROOT
N_SHUFFLE = 1000
N_PERMUTATION = 100
SEED = 42
K_VALUES = (1, 2, 4, 8)


#%% helpers
def median_iqr(values):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    return tuple(np.percentile(values, [25, 50, 75])) if len(values) else (np.nan, np.nan, np.nan)


def safe_spearman(left, right):
    finite = np.isfinite(left) & np.isfinite(right)
    if finite.sum() < 3 or np.ptp(left[finite]) == 0 or np.ptp(right[finite]) == 0:
        return np.nan
    return float(spearmanr(left[finite], right[finite]).statistic)


def permute_trial_labels(response, keep, seed=SEED):
    rng = np.random.RandomState(seed)
    permuted = np.zeros((*response.shape, N_PERMUTATION), dtype=bool)
    for pseudotrial in range(keep.shape[1]):
        available = np.flatnonzero(keep[:, pseudotrial])
        n_response = int(response[available, pseudotrial].sum())
        if n_response == 0:
            continue
        if n_response == len(available):
            permuted[available, pseudotrial, :] = True
            continue
        scores = rng.random_sample((N_PERMUTATION, len(available)))
        chosen = np.argpartition(scores, n_response - 1, axis=1)[:, :n_response]
        permuted[available[chosen], pseudotrial,
                 np.arange(N_PERMUTATION)[:, None]] = True
    return permuted


def grouped_mean(values, groups, valid, k):
    result = np.zeros((values.shape[0], len(groups), values.shape[-1]), dtype=float)
    for member in range(k):
        result += np.nan_to_num(values[:, groups[:, member]], nan=0.0)
    result /= k
    if valid.ndim == 3:
        valid = valid[:, :, 0]
    result[~valid] = np.nan
    return result


#%% pseudotrials
def pseudotrial_result(aligned, keep, shuffled, groups, k):
    valid = np.all(keep[:, groups], axis=2)
    pseudo_aligned = grouped_mean(aligned, groups, valid, k)
    pseudo_shuffled = grouped_mean(shuffled, groups, valid, k)
    baseline = pseudo_aligned[:, :, 30:60]
    response = pseudo_aligned[:, :, 60:90]
    response_change = np.nanmean(response, axis=2) - np.nanmean(baseline, axis=2)
    baseline_change = np.nanmean(baseline, axis=2) - np.nanmean(
        pseudo_aligned[:, :, 0:30], axis=2)
    std = pooled_std(baseline, response)
    threshold = np.nanpercentile(
        pseudo_shuffled.reshape(pseudo_shuffled.shape[0], -1), 95, axis=1)
    effect = np.divide(response_change, std[:, None],
                        out=np.full(response_change.shape, np.nan),
                        where=np.isfinite(std[:, None]) & (std[:, None] != 0))
    baseline_effect = np.divide(baseline_change, std[:, None],
                                out=np.full(baseline_change.shape, np.nan),
                                where=np.isfinite(std[:, None]) & (std[:, None] != 0))
    responding = (response_change > threshold[:, None]) & (effect > 0.05) & valid
    baseline_responding = (
        (baseline_change > threshold[:, None])
        & (baseline_effect > 0.05)
        & valid
    )
    del pseudo_aligned, pseudo_shuffled
    return response_change, baseline_change, responding, valid, baseline_responding


def partition_metrics(response_change, baseline_change, responding,
                      valid, baseline_responding, up):
    observed_overlap, _, n_pairs = jaccard_summary(
        responding, responding[:, :, None], valid)
    permutation = permute_trial_labels(responding, valid)
    _, null_overlap, _ = jaccard_summary(responding, permutation, valid)
    del permutation

    n_up_valid = valid[up].sum(axis=0)
    up_mean_response = np.divide(np.nansum(response_change[up], axis=0),
                                 n_up_valid, out=np.full(valid.shape[1], np.nan),
                                 where=n_up_valid > 0)
    up_mean_baseline = np.divide(np.nansum(baseline_change[up], axis=0),
                                 n_up_valid, out=np.full(valid.shape[1], np.nan),
                                 where=n_up_valid > 0)
    other_response = responding[~up].sum(axis=0)
    other_baseline = baseline_responding[~up].sum(axis=0)
    responding_count = responding.sum(axis=0)
    up_fraction = np.divide(responding[up].sum(axis=0), responding_count,
                            out=np.full(valid.shape[1], np.nan),
                            where=responding_count > 0)
    return dict(
        observed_overlap=observed_overlap,
        null_overlap=float(np.nanmedian(null_overlap)),
        responding_up_fraction=float(np.nanmean(up_fraction)),
        rho_response=safe_spearman(up_mean_response, other_response),
        rho_baseline=safe_spearman(up_mean_baseline, other_baseline),
        n_pairs=n_pairs,
    )


def recording_summary(recname):
    behaviour, profile, traces = load_recording(recname)
    frames = event_frames(behaviour)
    result = quantify_grid_response(traces, frames, n_shuffle=N_SHUFFLE,
                                    seed=SEED, return_single_trial=True)
    aligned = result['aligned']
    keep = result['keep']
    shuffled = result['shuffle_single_trial_amplitude']
    n_trials = len(frames)
    n_kept_trials = int(keep.any(axis=0).sum())
    up = profile['Up'].to_numpy(dtype=bool)
    rows = []
    for k in K_VALUES:
        n_pseudotrial = n_trials // k
        partition_rows = []
        for seed in range(10):
            groups = np.random.RandomState(seed).permutation(n_trials)
            groups = groups[:n_pseudotrial * k].reshape(n_pseudotrial, k)
            values = pseudotrial_result(aligned, keep, shuffled, groups, k)
            partition_rows.append(partition_metrics(*values, up))
        row = dict(
            rec_id=recname, animal=recname.split('-')[0], k=k,
            n_trials=n_trials, n_kept_trials=n_kept_trials,
            n_valid_acr=len(profile), n_up_acr=int(up.sum()),
            n_partitions=10, n_pseudotrials=n_pseudotrial,
            n_discarded_remainder=n_trials % k,
            n_overlap_pairs=int(np.nanmean([x['n_pairs'] for x in partition_rows])),
        )
        for key in ('observed_overlap', 'null_overlap',
                    'responding_up_fraction', 'rho_response', 'rho_baseline'):
            values = np.asarray([x[key] for x in partition_rows], dtype=float)
            row[key] = float(np.nanmean(values))
            row[f'{key}_n_partitions'] = int(np.isfinite(values).sum())
        rows.append(row)
        print(f'{recname}, k={k}: observed={row["observed_overlap"]:.6g}, '
              f'null={row["null_overlap"]:.6g}, rho={row["rho_response"]:.6g}',
              flush=True)
    return rows


#%% log and run
def append_log(table):
    total_trials = table.groupby('rec_id').n_trials.first().sum()
    lines = [
        '',
        '## 24 November 2026: aligned response pseudotrials',
        '',
        f'I used all {total_trials} valid trials from {table.rec_id.nunique()} '
        f'recordings in {table.animal.nunique()} animals, with Jingyu\'s per-ACR '
        f'full-window NaN rule. Each k partition used disjoint random groups with '
        f'seeds 0 to 9, dropped the remainder, averaged all k aligned trials, and '
        f'required all k trials to be kept for that ACR. Pseudotrial thresholds used '
        f'1,000 independently shifted-trial averages with seed 42 and the 95th-'
        f'percentile and effect-size rules. The count-preserving overlap null used '
        f'{N_PERMUTATION} permutations per partition.',
        '',
    ]
    for k in K_VALUES:
        subset = table.loc[table.k.eq(k)]
        observed = subset.observed_overlap.to_numpy()
        null = subset.null_overlap.to_numpy()
        overlap_test = wilcoxon(observed - null, alternative='two-sided')
        z_difference = np.array([
            np.arctanh(response) - np.arctanh(baseline)
            for response, baseline in zip(subset.rho_response, subset.rho_baseline)
            if np.isfinite(response) and np.isfinite(baseline)
            and abs(response) < 1 and abs(baseline) < 1
        ])
        z_test = wilcoxon(z_difference, alternative='two-sided')
        low, median, high = median_iqr(observed)
        null_low, null_median, null_high = median_iqr(null)
        frac_low, frac_median, frac_high = median_iqr(subset.responding_up_fraction)
        lines += [
            f'k={k}, observed overlap: median={median:.9g}, IQR '
            f'[{low:.9g}, {high:.9g}]; null overlap median={null_median:.9g}, '
            f'IQR [{null_low:.9g}, {null_high:.9g}]; observed against null '
            f'W={overlap_test.statistic:.9g}, p={overlap_test.pvalue:.12g}; '
            f'n={len(subset)} recordings, n={subset.animal.nunique()} animals.',
            f'k={k}, Fisher-z(response) minus Fisher-z(baseline): median='
            f'{np.median(z_difference):.9g}, IQR [{np.percentile(z_difference,25):.9g}, '
            f'{np.percentile(z_difference,75):.9g}]; W={z_test.statistic:.9g}, '
            f'p={z_test.pvalue:.12g}; n={len(z_difference)} recordings, '
            f'n={subset.animal.nunique()} animals.',
            f'k={k}, responding Up fraction: median={frac_median:.9g}, IQR '
            f'[{frac_low:.9g}, {frac_high:.9g}]; n={len(subset)} recordings, '
            f'n={subset.animal.nunique()} animals.',
        ]
    lines += ['', 'Per-recording results: outputs/da_hotspots/pseudotrials.csv.', '']
    with (OUTPUT_ROOT / 'PROGRESS.md').open('a', encoding='utf-8') as stream:
        stream.write('\n'.join(lines))


def main():
    output = OUTPUT_ROOT / 'pseudotrials.csv'
    if output.exists():
        raise FileExistsError(output)
    rows = []
    for recname in rec_list.recDLightDbh:
        rows.extend(recording_summary(recname))
    table = pd.DataFrame(rows)
    table.to_csv(output, index=False)
    append_log(table)


if __name__ == '__main__':
    with threadpool_limits(limits=1):
        main()
