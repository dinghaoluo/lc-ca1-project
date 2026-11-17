# -*- coding: utf-8 -*-
'''
Created on Thu Nov 12 11:11:55 2026

trial-label permutation reference for early and late amplitude KS distances

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
from scipy.stats import ks_2samp, wilcoxon
from threadpoolctl import threadpool_limits

from analyse_early_late_ecdf_and_prediction import condition
from grid_response_test import align_trials, keep_aligned_trials, trial_amplitudes

repo_root   = Path(__file__).resolve().parents[3]
OUTPUT_ROOT = pp.DLIGHT_HOTSPOTS_OUTPUT_ROOT
N_PERMUTE   = 1000
SEED        = 42


#%% amplitude distributions
def group_amplitude(amplitude, trials):
    values = amplitude[:, trials]
    count = np.isfinite(values).sum(axis=1)
    return np.divide(np.nansum(values, axis=1), count,
                     out=np.full(len(values), np.nan), where=count > 0)


def ks_distance(first, second):
    first = first[np.isfinite(first)]
    second = second[np.isfinite(second)]
    if not len(first) or not len(second):
        return np.nan
    return float(ks_2samp(first, second).statistic)


#%% trial-label permutations
def main():
    output = OUTPUT_ROOT / 'early_late_ks_permutation.csv'
    if output.exists():
        raise FileExistsError(f'Existing output left unchanged: {output}')
    membership = pd.read_parquet(OUTPUT_ROOT / 'aligned_early_late_per_recording.parquet')
    recordings = membership.loc[membership.unmatched_included, 'rec_id'].tolist()
    selections = pd.read_csv(OUTPUT_ROOT / 'followup_early_late_trials.csv')
    selections = selections.loc[selections.variant == 'unmatched']
    saved = pd.read_csv(OUTPUT_ROOT / 'followup_early_late_amplitudes.csv',
                        float_precision='round_trip')
    saved = saved.loc[saved.variant == 'unmatched']
    rows = []
    for recname in recordings:
        behaviour, profile, traces = condition.load_recording(recname)
        selected = selections.loc[selections.rec_id == recname]
        early = selected.loc[selected.group == 'early', 'trial_index'].to_numpy(dtype=int)
        late = selected.loc[selected.group == 'late', 'trial_index'].to_numpy(dtype=int)
        pooled = np.r_[early, late]
        assert len(early) == len(late) and len(early) > 0
        assert len(np.unique(pooled)) == len(pooled)
        _, eligible_early, eligible_late, _ = condition.valid_trials(behaviour)
        assert set(early).issubset(eligible_early) and set(late).issubset(eligible_late)
        run = np.asarray(behaviour['run_onset_frames'], dtype=int)
        for label, trials in (('early', early), ('late', late)):
            stored = selected.loc[selected.group == label, 'run_onset_frame'].to_numpy(dtype=int)
            assert np.array_equal(run[trials], stored)
        cleaned, keep = keep_aligned_trials(align_trials(traces, run[pooled]))
        amplitude, _, _ = trial_amplitudes(cleaned)
        observed_early = group_amplitude(amplitude, np.arange(len(early)))
        observed_late = group_amplitude(amplitude, np.arange(len(early), len(pooled)))
        observed_d = ks_distance(observed_early, observed_late)
        reference = saved.loc[saved.rec_id == recname]
        coords = pd.MultiIndex.from_tuples(profile['roi_id'].tolist(), names=['roi_y', 'roi_x'])
        saved_values = []
        differences = []
        missing_agree = True
        for label, values in (('early', observed_early), ('late', observed_late)):
            old = reference.loc[reference.group == label].set_index(['roi_y', 'roi_x'])
            previous = old.reindex(coords).amplitude.to_numpy()
            saved_values.append(previous)
            missing_agree &= bool(np.array_equal(np.isfinite(values), np.isfinite(previous)))
            finite = np.isfinite(values) & np.isfinite(previous)
            differences.extend(np.abs(values[finite] - previous[finite]))
        saved_d = ks_distance(*saved_values)
        rng = np.random.RandomState(SEED)
        null_d = np.full(N_PERMUTE, np.nan)
        for draw in range(N_PERMUTE):
            reassigned = rng.permutation(len(pooled))
            first = group_amplitude(amplitude, reassigned[:len(early)])
            second = group_amplitude(amplitude, reassigned[len(early):])
            null_d[draw] = ks_distance(first, second)
        finite_null = null_d[np.isfinite(null_d)]
        exceed = int(np.sum(finite_null >= observed_d))
        p_value = (exceed + 1) / (len(finite_null) + 1) if np.isfinite(observed_d) and len(finite_null) else np.nan
        rows.append(dict(rec_id=recname, animal=recname.split('-')[0], n_early_trials=len(early),
                         n_late_trials=len(late), n_valid_acr=len(profile),
                         n_early_finite_acr=int(np.isfinite(observed_early).sum()),
                         n_late_finite_acr=int(np.isfinite(observed_late).sum()),
                         observed_D=observed_d, null_median=float(np.median(finite_null)),
                         null_95th=float(np.percentile(finite_null, 95)),
                         permutation_p=p_value, n_permutations=N_PERMUTE,
                         n_finite_permutations=len(finite_null), n_null_ge_observed=exceed,
                         saved_observed_D=saved_d, max_saved_amplitude_difference=max(differences),
                         saved_finite_mask_agrees=missing_agree))
        print(f'{recname}: {len(early)} trials per group, D={observed_d:.6g}, null median={rows[-1]["null_median"]:.6g}, p={p_value:.6g}', flush=True)
    table = pd.DataFrame(rows)
    table.to_csv(output, index=False)
    valid = (np.isfinite(table.observed_D) & np.isfinite(table.null_median)
             & np.isfinite(table.permutation_p))
    tested = table.loc[valid]
    difference = tested.observed_D.to_numpy() - tested.null_median.to_numpy()
    if np.all(difference == 0):
        statistic, p_value = 0.0, 1.0
    else:
        test = wilcoxon(difference, alternative='two-sided')
        statistic, p_value = float(test.statistic), float(test.pvalue)
    significant = int((tested.permutation_p < 0.05).sum())
    low, high = np.percentile(difference, [25, 75])
    lines = ['', '## 12 November 2026: early/late KS trial-label permutation reference', '',
             'I used the saved unmatched equal-number subsamples in followup_early_late_trials.csv for every included recording in analyse_early_late_ecdf_and_prediction.py, preserving the saved trial order and run-onset frames. All 57 recordings were retained, including recordings with small groups; the late-trial cutoff used for the separate split-half figure was not applied here.', '',
             'The imported loader prepares the full-grid traces and selects the saved valid ACRs. The imported alignment, complete-window NaN rule and trial_amplitudes define response [0, 1) s minus baseline [-1, 0) s around run onset. Each group averages its kept trials separately for each ACR, omitting non-finite ACR means separately from the two distributions. Observed D is the two-sample KS statistic on these finite distributions, over their full range, exactly as in plot_split_half_ecdf.py.', '',
             'For each recording I pooled only its two selected groups and used a fresh RandomState(42) to permute the trial indices 1,000 times, assigning the first original group size to one group and the remainder to the other. The same reassignment applies to every ACR, and the original missing windows travel with their trials. Recomputing means from the fixed per-trial amplitudes is equivalent to realigning each reassigned group because neither alignment nor the NaN rule depends on the group labels.', '',
             'The null summary is its median and 95th percentile. The upper-tail permutation p-value is (1 + number of null D >= observed D) / (1 + number of finite permutations), including ties; this Monte Carlo correction prevents zero p-values. All 1,000 permutations are expected in each denominator; any non-finite counts are retained in the CSV. Only the KS distance is used from the two-sample routine, not its parametric p-value.', '',
             f'Across recordings, {significant}/{len(tested)} ({100 * significant / len(tested):.6g}%) have permutation p < 0.05, against the requested nominal 5% reference. The finite permutation reference and ties can make the attainable rejection rate below 5%; no separate test of this fraction was requested.', '',
             f'Observed recording median D={tested.observed_D.median():.9g}; null-median recording median={tested.null_median.median():.9g}. Median paired observed D minus null median={np.median(difference):.9g}, recording IQR [{low:.9g}, {high:.9g}]. Paired two-sided Wilcoxon signed-rank test: W={statistic:g}, p={p_value:.9g}, n={len(tested)} recordings, n={tested.animal.nunique()} animals; non-finite recordings omitted={int((~valid).sum())}.', '',
             f'Agreement with saved amplitudes: maximum absolute difference={table.max_saved_amplitude_difference.max():.9g} dF/F; all finite masks agree={bool(table.saved_finite_mask_agrees.all())}; maximum absolute difference from saved observed KS D={np.max(np.abs(table.observed_D - table.saved_observed_D)):.9g}. Finite null counts range from {table.n_finite_permutations.min()} to {table.n_finite_permutations.max()}.', '',
             'Per-recording results:', '',
             '| recording | trials per group | observed D | null median | null 95th percentile | permutation p |',
             '| --- | ---: | ---: | ---: | ---: | ---: |']
    for row in table.itertuples():
        lines.append(f'| {row.rec_id} | {row.n_early_trials} | {row.observed_D:.9g} | {row.null_median:.9g} | {row.null_95th:.9g} | {row.permutation_p:.9g} |')
    lines += ['', 'Table: outputs/da_hotspots/early_late_ks_permutation.csv. Trial-label exchangeability is the assumption behind this reference: pooling preserves the recorded trials and ACR dependence while removing their early/late labels. No figures or other analyses were produced.']
    with (OUTPUT_ROOT / 'PROGRESS.md').open('a', encoding='utf-8') as stream:
        stream.write(chr(10).join(lines) + chr(10))
    print(f'Complete: {significant}/{len(tested)} p < 0.05; W={statistic:g}, p={p_value:.9g}, n={len(tested)}, animals={tested.animal.nunique()}', flush=True)


#%% run
if __name__ == '__main__':
    with threadpool_limits(limits=1):
        main()
