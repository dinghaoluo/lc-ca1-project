# -*- coding: utf-8 -*-
'''
Created on Wed Nov 11 15:20:24 2026

early and late trial-count sensitivity and single-trial amplitude variance

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
from scipy.stats import wilcoxon
from threadpoolctl import threadpool_limits

from analyse_acr_consistency_aligned import (load_recording, event_frames,
                                              recnames)
from grid_response_test import quantify_grid_response

repo_root   = Path(__file__).resolve().parents[3]
OUTPUT_ROOT = pp.DLIGHT_HOTSPOTS_OUTPUT_ROOT
N_SHUFFLE   = 1000
SEED        = 42


#%% early and late sensitivity
def early_late_sensitivity():
    saved = pd.read_parquet(OUTPUT_ROOT / 'aligned_early_late_per_recording.parquet')
    smaller = saved[['n_early', 'n_late']].min(axis=1)
    excluded = saved.loc[smaller < 10]
    retained = saved.loc[saved.unmatched_included & (smaller >= 10)]
    lines = ['', '## 11 November 2026: unmatched early/late sensitivity to trial count', '',
             f'I repeated the paired two-sided Wilcoxon tests of late minus early from aligned_early_late_per_recording.parquet, retaining the saved unmatched estimates and excluding recordings with fewer than 10 trials in the smaller behavioural group. This leaves {len(retained)} of {len(saved)} recordings, from {retained.animal.nunique()} animals; the cutoff precedes the per-ACR NaN rule, and no trial selections or response estimates were recalculated.', '',
             'Excluded recordings (early, late trial counts):']
    for row in excluded.itertuples():
        lines.append(f'- {row.rec_id}: early={row.n_early}, late={row.n_late}.')
    for measure, label, units in (('pct_up', 'percentage of Up ACRs', 'percentage points'),
                                  ('median_amplitude', 'median response amplitude', 'dF/F')):
        early = retained[f'unmatched_early_{measure}'].to_numpy(dtype=float)
        late = retained[f'unmatched_late_{measure}'].to_numpy(dtype=float)
        finite = np.isfinite(early) & np.isfinite(late)
        difference = late[finite] - early[finite]
        if np.all(difference == 0):
            statistic, p_value = 0.0, 1.0
        else:
            test = wilcoxon(difference, alternative='two-sided')
            statistic, p_value = float(test.statistic), float(test.pvalue)
        low, high = np.percentile(difference, [25, 75])
        n_animals = retained.loc[finite, 'animal'].nunique()
        lines.append(f'{label}: early recording median={np.median(early[finite]):.9g}, late recording median={np.median(late[finite]):.9g}; median paired late-minus-early difference={np.median(difference):.9g} {units}, recording IQR [{low:.9g}, {high:.9g}]; W={statistic:g}, p={p_value:.9g}, n={finite.sum()} recordings, n={n_animals} animals; non-finite pairs omitted={int((~finite).sum())}.')
    with (OUTPUT_ROOT / 'PROGRESS.md').open('a', encoding='utf-8') as stream:
        stream.write(chr(10).join(lines) + chr(10))
    print(chr(10).join(lines), flush=True)


#%% variance across kept trials
def amplitude_variance():
    rows = []
    for recname in list(recnames):
        behaviour, profile, traces = load_recording(recname)
        frames = event_frames(behaviour)
        result = quantify_grid_response(traces, frames, n_shuffle=N_SHUFFLE,
                                        seed=SEED, return_single_trial=True)
        keep = result['keep']
        kept = keep.sum(axis=1)
        up = profile['Up'].to_numpy(dtype=bool)
        observed_variance = np.full(len(profile), np.nan)
        shuffle_variance = np.full(len(profile), np.nan)
        for acr in range(len(profile)):
            if kept[acr] < 2:
                continue
            observed = result['trial_amplitude'][acr, keep[acr]]
            shuffled = result['shuffle_single_trial_amplitude'][acr, keep[acr], :]
            observed_variance[acr] = np.var(observed, ddof=1)
            shuffle_variance[acr] = np.median(np.var(shuffled, axis=0, ddof=1))
        ratio = np.divide(observed_variance, shuffle_variance,
                          out=np.full(len(profile), np.nan),
                          where=np.isfinite(shuffle_variance) & (shuffle_variance > 0))
        for label, group in (('Up', up), ('other', ~up)):
            finite = group & np.isfinite(ratio)
            median_ratio = float(np.median(ratio[finite])) if finite.any() else np.nan
            rows.append(dict(rec_id=recname, animal=recname.split('-')[0], group=label,
                             n_trials=len(frames), n_acr=int(group.sum()), n_ratio=int(finite.sum()),
                             n_fewer_two=int(np.sum(group & (kept < 2))),
                             n_zero_observed=int(np.sum(group & (observed_variance == 0))),
                             n_invalid_shuffle=int(np.sum(group & (kept >= 2)
                                                          & (~np.isfinite(shuffle_variance) | (shuffle_variance <= 0)))),
                             median_ratio=median_ratio))
        first, second = rows[-2], rows[-1]
        print(f'{recname}: {len(frames)} valid trials; Up median variance ratio={first["median_ratio"]:.6g} ({first["n_ratio"]} ACRs), other={second["median_ratio"]:.6g} ({second["n_ratio"]} ACRs)', flush=True)
        del result
    table = pd.DataFrame(rows)
    lines = ['', '## 11 November 2026: observed versus circular-shuffle amplitude variance', '',
             f'I used the validated CPU test with return_single_trial=True on all valid trials for {table.rec_id.nunique()} recordings from {table.animal.nunique()} animals, with 1,000 circular shuffles and seed 42 reset per recording. Valid trials have finite reward and run onset, run onset other than -1, and both stop flags zero, with the complete finite [-2, 4] s window retained separately for each saved valid ACR. Single-trial amplitude is response [0, 1] s minus baseline [-1, 0] s around run onset, in dF/F units.', '',
             'For each ACR I computed sample variance across kept trials (ddof=1), then the same variance for each shuffle and its median over 1,000 shuffles; the ACR ratio is observed variance divided by that median shuffle variance. Both variances use the same kept trials, so their ddof factors cancel in the ratio. ACRs with fewer than two kept trials or a non-positive or non-finite shuffle denominator cannot supply a ratio. Finite zero observed ratios are retained in the recording median; only recordings with a finite positive median ratio can enter the log-ratio test. Up groups use the saved classification and the other group contains all remaining valid ACRs.', '',
             'Each recording supplies the median of its ACR ratios separately for Up and other ACRs. The two-sided Wilcoxon test is applied to the natural log of these recording medians against zero, with recordings as the paired units; the IQR describes variation between recordings and is not a confidence interval.']
    for label in ('Up', 'other'):
        group = table.loc[table.group == label]
        valid = np.isfinite(group.median_ratio) & (group.median_ratio > 0)
        tested = group.loc[valid]
        log_ratio = np.log(tested.median_ratio.to_numpy())
        if np.all(log_ratio == 0):
            statistic, p_value = 0.0, 1.0
        else:
            test = wilcoxon(log_ratio, alternative='two-sided')
            statistic, p_value = float(test.statistic), float(test.pvalue)
        low, high = np.percentile(tested.median_ratio, [25, 75])
        lines.append(f'{label}: median of recording variance ratios={tested.median_ratio.median():.9g}, recording IQR [{low:.9g}, {high:.9g}]; median log ratio={np.median(log_ratio):.9g}; W={statistic:g}, p={p_value:.9g}, n={len(tested)} recordings, n={tested.animal.nunique()} animals, n={int(tested.n_ratio.sum())} ACR ratios across included recordings; recordings omitted from log test={int((~valid).sum())}.')
        lines.append(f'{label} exclusions and boundaries: {int(group.n_fewer_two.sum())} ACRs with fewer than two kept trials, {int(group.n_invalid_shuffle.sum())} with a non-positive or non-finite shuffle variance despite at least two kept trials, {int(group.n_zero_observed.sum())} zero observed variances. Non-finite ACR ratios omitted={int((group.n_acr - group.n_ratio).sum())}.')
    lines += ['', '| recording | valid trials | Up ACRs with ratio / valid | Up median ratio | other ACRs with ratio / valid | other median ratio |',
              '| --- | ---: | ---: | ---: | ---: | ---: |']
    for recname in table.rec_id.unique():
        group = table.loc[table.rec_id == recname].set_index('group')
        first, second = group.loc['Up'], group.loc['other']
        lines.append(f'| {recname} | {int(first.n_trials)} | {int(first.n_ratio)} / {int(first.n_acr)} | {first.median_ratio:.9g} | {int(second.n_ratio)} / {int(second.n_acr)} | {second.median_ratio:.9g} |')
    lines += ['', 'The saved Up grouping is descriptive because it was selected from trial-averaged responses in these recordings; the variance comparison uses those same trials. Circular shifts move the windows within the aligned trace, so this comparison measures variance at run onset relative to variance under that temporal reference. No responding threshold, injected amplitude, new figure or additional analysis was used.']
    with (OUTPUT_ROOT / 'PROGRESS.md').open('a', encoding='utf-8') as stream:
        stream.write(chr(10).join(lines) + chr(10))
    print(chr(10).join(lines[:12]), flush=True)


#%% run
if __name__ == '__main__':
    with threadpool_limits(limits=1):
        early_late_sensitivity()
        amplitude_variance()
