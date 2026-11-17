# -*- coding: utf-8 -*-
'''
Created on Wed Nov 11 09:12:02 2026

check trial coverage and null expectations in aligned dLight analyses

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

import traceback


import numpy as np
import pandas as pd
from scipy.stats import wilcoxon
from threadpoolctl import threadpool_limits

import compare_early_late_aligned as condition
import analyse_acr_consistency_aligned as consistency
from grid_response_test import (align_trials, keep_aligned_trials,
                                quantify_grid_response, classify_up)

ROOT = Path(__file__).resolve().parents[3]
OUT = ROOT / 'outputs' / 'da_hotspots'
N_SHUFFLE = 1000
N_PERMUTATION = 1000
SEED = 42
THRESHOLDS = (90, 95, 99)


#%% nulls
def leave_one_out_percentages(null_k, usable):
    p_values = np.empty(null_k.shape, dtype=float)
    for acr, counts in enumerate(null_k):
        at_least = len(counts) - np.searchsorted(np.sort(counts), counts, side='left')
        p_values[acr] = (at_least - 1) / (len(counts) - 1)
    return 100 * (p_values[usable] < 0.05).mean(axis=0)


def permute_trial_labels(response, keep, n_permutation=N_PERMUTATION, seed=SEED):
    rng = np.random.RandomState(seed)
    permuted = np.zeros((*response.shape, n_permutation), dtype=bool)
    for trial in range(keep.shape[1]):
        available = np.flatnonzero(keep[:, trial])
        n_response = int(response[available, trial].sum())
        if n_response == 0:
            continue
        if n_response == len(available):
            permuted[available, trial, :] = True
            continue
        scores = rng.random_sample((n_permutation, len(available)))
        chosen = np.argpartition(scores, n_response - 1, axis=1)[:, :n_response]
        permuted[available[chosen], trial, np.arange(n_permutation)[:, None]] = True
    return permuted


def interval(values):
    values = np.asarray(values, dtype=float)
    if not np.any(np.isfinite(values)):
        return np.nan, np.nan, np.nan
    low, median, high = np.nanpercentile(values, [2.5, 50, 97.5])
    return float(median), float(low), float(high)


def paired_test(rows, observed, null, name, threshold):
    usable = rows.loc[np.isfinite(rows[observed]) & np.isfinite(rows[null])]
    difference = usable[observed].to_numpy() - usable[null].to_numpy()
    if not len(difference):
        w, p = np.nan, np.nan
    elif np.all(difference == 0):
        w, p = 0.0, 1.0
    else:
        test = wilcoxon(difference, alternative='two-sided')
        w, p = float(test.statistic), float(test.pvalue)
    return {'test': name, 'threshold': threshold, 'n_recordings': len(usable),
            'n_animals': usable.animal.nunique(), 'W': w, 'p': p,
            'observed_recording_median': usable[observed].median(),
            'null_recording_median': usable[null].median(),
            'median_paired_difference': np.median(difference) if len(difference) else np.nan}


#%% counts and preservation checks
def first_draw_coverage(traces, behaviour, early, late, rng):
    n = min(len(early), len(late))
    first = {}
    for draw in range(condition.N_DRAWS):
        if len(early) <= len(late):
            groups = (early, list(rng.choice(late, size=n, replace=False)))
        else:
            groups = (list(rng.choice(early, size=n, replace=False)), late)
        if draw == 0:
            for label, trials in zip(('early', 'late'), groups):
                frames = np.asarray(behaviour['run_onset_frames'], dtype=int)[trials]
                aligned = align_trials(traces, frames)
                _, keep = keep_aligned_trials(aligned)
                first[f'{label}_first_draw_n_selected'] = n
                first[f'{label}_first_draw_mean_kept_per_acr'] = keep.sum(axis=1).mean()
    return first






def coverage_row(recname, animal, keep):
    counts = keep.sum(axis=1)
    return {'rec_id': recname, 'animal': animal, 'n_recordings': 1,
            'n_valid_acr': keep.shape[0], 'n_consistency_trials': keep.shape[1],
            'n_possible_acr_trials': keep.size, 'n_kept_acr_trials': int(keep.sum()),
            'fraction_kept': float(keep.mean()),
            'kept_trials_per_acr_min': int(counts.min()),
            'kept_trials_per_acr_median': float(np.median(counts)),
            'kept_trials_per_acr_max': int(counts.max()),
            'n_acr_fewer_than_10_kept': int(np.sum(counts < 10))}


#%% summary
def write_summary(tables, issues):
    counts, coverage, rates, up, response_count_null, overlap, variance, tests = tables
    pooled_coverage = coverage.loc[coverage.rec_id.eq('pooled')].iloc[0]
    pooled_up = up.loc[up.rec_id.eq('pooled')].iloc[0]
    lines = ['', '## 11 November 2026: check of aligned analyses', '',
             f'The check completed {len(counts)} recordings from {counts.animal.nunique()} animals, recomputing the kept windows and all 1,000 circular shuffles from the unchanged response module. No figures were produced.', '',
             'The first-draw coverage replays the global seed-42 subsampling sequence, in recording order with the speed-matched variant before the unmatched variant, consuming all 20 draws for every included variant. Excluded speed variants have no first draw and are recorded as missing rather than assigned a hypothetical draw.', '',
             'Coverage uses the consistency selector, including finite run onset other than -1, while behavioural valid-trial counts use reward and stop flags only. Each ACR-trial denominator includes every saved valid ACR and selected trial before the complete-window NaN rule. Pooled kept-trial quantiles concatenate the ACR counts across recordings.', '',
             f'Behaviour counts: {counts.n_behaviour_run_onset.sum()} trials, {counts.n_run_onset_minus1.sum()} with run onset -1, {counts.n_valid_trials.sum()} behaviourally valid trials; speed-matched inclusion={counts.speed_included.sum()}, unmatched inclusion={counts.unmatched_included.sum()}.']
    for reason, n in counts.loc[~counts.speed_included, 'speed_exclusion'].value_counts().items():
        lines.append(f'- Speed exclusion, {reason}: {n} recordings.')
    lines += ['', f'Pooled coverage: {pooled_coverage.n_kept_acr_trials:g}/{pooled_coverage.n_possible_acr_trials:g} ACR-trial pairs kept ({100 * pooled_coverage.fraction_kept:.3f}%); kept-trial count per ACR min/median/max={pooled_coverage.kept_trials_per_acr_min:g}/{pooled_coverage.kept_trials_per_acr_median:g}/{pooled_coverage.kept_trials_per_acr_max:g}; {pooled_coverage.n_acr_fewer_than_10_kept:g} ACRs have fewer than ten kept trials.', '',
              'Overall responding rates are ratios of total responding to total kept ACR-trials, with the shuffled rate averaged over the 1,000 shuffle columns; they are not medians of per-trial fractions.']
    for _, row in rates.loc[rates.rec_id.eq('pooled')].iterrows():
        lines.append(f'- {row.threshold:g}th threshold: observed={100 * row.observed_rate:.4f}%, shuffle mean={100 * row.shuffle_mean_rate:.4f}%, denominator={row.n_kept_acr_trials:g} kept ACR-trials.')
    lines += ['', f'Pooled Up classification: saved={pooled_up.saved_up_pct:.4f}%, CPU={pooled_up.cpu_up_pct:.4f}% across {pooled_up.n_valid_acr:g} ACRs; mean per-recording CPU minus saved difference={pooled_up.mean_recording_difference_pp:.4f} percentage points.', '',
              'For the corrected response-count null, each shuffle is treated as the observation and compared with the other 999 counts, using the fraction greater than or equal to it, including ties, with no added pseudocount. The observed p-values retain their original denominator of 1,000. Zero-kept-trial ACRs would be excluded from both percentages. The threshold and effect-size rule are held fixed.', '',
              'For the added Jaccard null, each permutation independently reallocates the observed responding labels among ACRs kept on each trial, preserving the responding count of that trial. Seed 42 is reset for each recording and threshold; the 1,000 independent columns therefore have reproducible coupling across thresholds. Both observed and null Jaccard calculations call the existing function, restricting each trial pair to ACRs kept on both trials and omitting zero unions separately in each matrix. Bands are the 2.5th and 97.5th percentiles.', '',
              'Paired two-sided Wilcoxon tests compare observed values with the corresponding recording-specific null median:']
    for _, row in tests.iterrows():
        lines.append(f'- {row.test}, {row.threshold:g}th threshold: observed recording median={row.observed_recording_median:.6g}, null recording median={row.null_recording_median:.6g}, W={row.W:g}, p={row.p:.8g}, n={row.n_recordings:g} recordings, n={row.n_animals:g} animals.')
    lines += ['', 'At the 95th threshold, responding-count variance uses population variance (ddof=0) across the complete consistency trial axis, including selected trials with no kept ACRs as zero responding; this follows the existing responding matrix rather than dropping those columns. The null median is taken over the 1,000 shuffle variances. No additional test of variance or responding rate was requested or performed.', '',
              f'Variance across recordings: observed median={variance.observed_variance.median():.6g}, shuffle-median variance median={variance.shuffle_variance_median.median():.6g}.', '',
              'CSV outputs under outputs/da_hotspots/: aligned_check_counts.csv, aligned_check_coverage.csv, aligned_check_responding_rates.csv, aligned_check_up_classification.csv, aligned_check_response_count_null.csv, aligned_check_overlap.csv, aligned_check_response_count_variance_95.csv, aligned_check_tests.csv, aligned_check_null_draws.csv and aligned_check_preserved_files.csv. The null-draw table retains the leave-one-out percentage for each recording, label-permutation overlap, circular-shuffle overlap and responding count for every threshold and replicate.', '',
              'Ambiguities and missing inputs:']
    lines += ['- ' + issue for issue in issues] if issues else ['- No missing input or mismatch with the saved aligned tables was found.']
    lines += ['', 'The existing measure-2 comparison against a fixed 5% remains unchanged in its original output; aligned_check_response_count_null.csv and aligned_check_tests.csv contain the discrete-null comparison requested here.']
    with (OUT / 'PROGRESS.md').open('a', encoding='utf-8') as handle:
        handle.write('\n'.join(lines) + '\n')


#%% run
def main():
    recordings = list(consistency.recnames)
    assert len(recordings) == 57
    saved_acr = pd.read_parquet(OUT / 'aligned_acr_consistency_per_acr.parquet')
    saved_rec = pd.read_parquet(OUT / 'aligned_acr_consistency_per_recording.parquet').set_index('rec_id')
    saved_conditions = pd.read_parquet(OUT / 'aligned_early_late_per_recording.parquet').set_index('rec_id')
    saved_validation = pd.read_parquet(OUT / 'aligned_validation_per_acr.parquet')
    counts_rows, coverage_rows, rate_rows, up_rows = [], [], [], []
    response_count_rows, overlap_rows, variance_rows, null_rows = [], [], [], []
    issues, all_kept_counts = [], []
    rng = np.random.RandomState(SEED)
    for recname in recordings:
        try:
            behaviour, profile, traces = consistency.load_recording(recname)
            animal = recname.split('-')[0]
            valid, early, late, _ = condition.valid_trials(behaviour)
            vectors = condition.speed_vectors(behaviour, len(valid))
            early_speed, late_speed, reason = condition.speed_match(early, late, vectors)
            if not early or not late:
                reason = 'no early or late trials'
            row = {'rec_id': recname, 'animal': animal}
            for key in ('run_onset_frames', 'reward_times', 'non_stop_trials',
                        'non_fullstop_trials', 'lick_times_aligned', 'speed_times_aligned'):
                row['n_behaviour_' + key.replace('_frames', '')] = len(behaviour[key])
            lengths = [len(behaviour[key]) for key in ('run_onset_frames', 'reward_times',
                       'non_stop_trials', 'non_fullstop_trials', 'lick_times_aligned', 'speed_times_aligned')]
            if len(set(lengths)) != 1:
                issues.append(f'{recname}: behaviour array lengths differ: {lengths}; existing selectors retained.')
            run = np.asarray(behaviour['run_onset_frames'], dtype=float)
            row.update(n_run_onset_minus1=int(np.sum(run == -1)),
                       n_valid_trials=int(valid.sum()),
                       n_valid_run_onset_minus1=int(np.sum(valid & (run[:len(valid)] == -1))),
                       n_early=len(early), n_late=len(late),
                       n_early_speed=len(early_speed), n_late_speed=len(late_speed),
                       speed_included=not bool(reason), speed_exclusion=reason,
                       unmatched_included=bool(early and late),
                       unmatched_exclusion='' if early and late else 'no early or late trials',
                       n_valid_acr=len(profile))
            for variant, eg, lg, included in (('speed', early_speed, late_speed, not bool(reason)),
                                               ('unmatched', early, late, bool(early and late))):
                values = first_draw_coverage(traces, behaviour, eg, lg, rng) if included else {
                    label + suffix: np.nan for label in ('early', 'late')
                    for suffix in ('_first_draw_n_selected', '_first_draw_mean_kept_per_acr')}
                row.update({variant + '_' + key: value for key, value in values.items()})
            for key in ('n_early', 'n_late', 'n_early_speed', 'n_late_speed', 'speed_included', 'unmatched_included'):
                if row[key] != saved_conditions.loc[recname, key]:
                    issues.append(f'{recname}: {key} differs from saved condition table.')
            frames = consistency.event_frames(behaviour)
            result = quantify_grid_response(traces, frames, n_shuffle=N_SHUFFLE, seed=SEED,
                                            return_single_trial=True)
            keep = result['keep']
            n_kept = keep.sum(axis=1)
            usable = n_kept > 0
            all_kept_counts.extend(n_kept.tolist())
            coverage_rows.append(coverage_row(recname, animal, keep))
            row['n_consistency_trials'] = len(frames)
            counts_rows.append(row)
            up_cpu, _ = classify_up(result['amplitude'], result['effect_size'],
                                    result['shuffle_amplitude'], profile['edge'].to_numpy())
            up_rows.append({'rec_id': recname, 'animal': animal, 'n_valid_acr': len(profile),
                            'n_saved_up': int(profile.Up.sum()), 'n_cpu_up': int(up_cpu.sum()),
                            'saved_up_pct': 100 * profile.Up.mean(), 'cpu_up_pct': 100 * up_cpu.mean(),
                            'cpu_minus_saved_pp': 100 * (up_cpu.mean() - profile.Up.mean())})
            old = saved_acr.loc[saved_acr.rec_id.eq(recname)].set_index(['roi_y', 'roi_x'])
            keys = pd.MultiIndex.from_tuples([tuple(x) for x in profile.roi_id], names=['roi_y', 'roi_x'])
            old = old.reindex(keys)
            if not np.array_equal(n_kept, old.n.to_numpy()):
                issues.append(f'{recname}: kept-trial counts differ from the saved consistency table.')
            validation = saved_validation.loc[saved_validation.rec_id.eq(recname)].set_index(['roi_y', 'roi_x']).reindex(keys)
            if not np.array_equal(up_cpu, validation.cpu_up.to_numpy()):
                issues.append(f'{recname}: CPU Up flags differ from the saved validation table.')
            for q in THRESHOLDS:
                _, real, shuffle = consistency.response_matrices(result, q)
                observed_k = real.sum(axis=1)
                null_k = shuffle.sum(axis=1)
                observed_p = (null_k >= observed_k[:, None]).mean(axis=1)
                observed_pct = 100 * (observed_p[usable] < 0.05).mean()
                null_pct = leave_one_out_percentages(null_k, usable)
                median, low, high = interval(null_pct)
                response_count_rows.append({'rec_id': recname, 'animal': animal, 'threshold': q,
                                     'n_acr': int(usable.sum()), 'observed_pct': observed_pct,
                                     'null_median_pct': median, 'null_low_pct': low, 'null_high_pct': high,
                                     'null_mean_pct': float(null_pct.mean())})
                null_responding = shuffle.sum(axis=(0, 1))
                rate_rows.append({'rec_id': recname, 'animal': animal, 'threshold': q,
                                  'n_kept_acr_trials': int(keep.sum()),
                                  'n_observed_responding': int(real.sum()),
                                  'shuffle_mean_n_responding': float(null_responding.mean()),
                                  'observed_rate': real.sum() / keep.sum(),
                                  'shuffle_mean_rate': null_responding.mean() / keep.sum()})
                observed_overlap, shuffle_overlap, n_pairs = consistency.jaccard_summary(real, shuffle, keep)
                permuted = permute_trial_labels(real, keep)
                _, permutation_overlap, _ = consistency.jaccard_summary(real, permuted, keep)
                pm, pl, ph = interval(permutation_overlap)
                sm, sl, sh = interval(shuffle_overlap)
                overlap_rows.append({'rec_id': recname, 'animal': animal, 'threshold': q,
                                     'n_observed_pairs': n_pairs, 'observed_overlap': observed_overlap,
                                     'permutation_null_median': pm, 'permutation_null_low': pl,
                                     'permutation_null_high': ph, 'shuffle_null_median': sm,
                                     'shuffle_null_low': sl, 'shuffle_null_high': sh,
                                     'saved_shuffle_null_median': saved_rec.loc[recname, f'overlap_null_median_{q}'],
                                     'n_finite_permutation': int(np.isfinite(permutation_overlap).sum()),
                                     'n_finite_shuffle': int(np.isfinite(shuffle_overlap).sum())})
                if not (np.array_equal(observed_k, old[f'k_{q}'].to_numpy())
                        and np.allclose(observed_p, old[f'empirical_p_{q}'].to_numpy(), equal_nan=True)
                        and np.isclose(observed_overlap, saved_rec.loc[recname, f'overlap_observed_{q}'])
                        and np.isclose(sm, saved_rec.loc[recname, f'overlap_null_median_{q}'])):
                    issues.append(f'{recname}, threshold {q}: recomputed response counts, p-values or overlap differ from saved outputs.')
                if not np.array_equal(permuted.sum(axis=0), np.repeat(real.sum(axis=0)[:, None], N_PERMUTATION, axis=1)):
                    raise RuntimeError('Permutation changed a trial responding count')
                if np.any(permuted & ~keep[:, :, None]):
                    raise RuntimeError('Permutation assigned a response to an unavailable ACR')
                if q == 95:
                    null_variance = np.var(shuffle.sum(axis=0), axis=0, ddof=0)
                    variance_rows.append({'rec_id': recname, 'animal': animal, 'threshold': q,
                                          'n_trials': len(frames), 'n_trials_with_kept_acr': int(keep.any(axis=0).sum()),
                                          'observed_variance': float(np.var(real.sum(axis=0), ddof=0)),
                                          'shuffle_variance_median': float(np.median(null_variance))})
                else:
                    null_variance = np.full(N_SHUFFLE, np.nan)
                null_rows.extend({'rec_id': recname, 'animal': animal, 'threshold': q, 'replicate': j,
                                  'leave_one_out_pct': null_pct[j],
                                  'permutation_overlap': permutation_overlap[j],
                                  'shuffle_overlap': shuffle_overlap[j],
                                  'shuffle_n_responding': int(null_responding[j]),
                                  'shuffle_response_count_variance_95': null_variance[j]}
                                 for j in range(N_SHUFFLE))
                del permuted, shuffle
            del result, traces
            print(f'{recname}: counts, coverage and three null comparisons complete', flush=True)
        except Exception as error:
            message = f'{recname}: {type(error).__name__}: {error}; remaining recordings continue.'
            issues.append(message)
            with (OUT / 'PROGRESS.md').open('a', encoding='utf-8') as handle:
                handle.write('\nMissing or failed input: ' + message + '\n')
            traceback.print_exc()
        pd.DataFrame(counts_rows).to_csv(OUT / 'aligned_check_counts.csv', index=False)
        pd.DataFrame(coverage_rows).to_csv(OUT / 'aligned_check_coverage.csv', index=False)
        pd.DataFrame(rate_rows).to_csv(OUT / 'aligned_check_responding_rates.csv', index=False)
        pd.DataFrame(response_count_rows).to_csv(OUT / 'aligned_check_response_count_null.csv', index=False)
        pd.DataFrame(overlap_rows).to_csv(OUT / 'aligned_check_overlap.csv', index=False)
        pd.DataFrame(variance_rows).to_csv(OUT / 'aligned_check_response_count_variance_95.csv', index=False)
    counts = pd.DataFrame(counts_rows)
    coverage = pd.DataFrame(coverage_rows)
    pooled = {'rec_id': 'pooled', 'animal': 'all', 'n_recordings': len(coverage),
              'n_consistency_trials': int(coverage.n_consistency_trials.sum())}
    for key in ('n_valid_acr', 'n_possible_acr_trials', 'n_kept_acr_trials', 'n_acr_fewer_than_10_kept'):
        pooled[key] = int(coverage[key].sum())
    pooled.update(fraction_kept=pooled['n_kept_acr_trials'] / pooled['n_possible_acr_trials'],
                  kept_trials_per_acr_min=np.min(all_kept_counts),
                  kept_trials_per_acr_median=np.median(all_kept_counts),
                  kept_trials_per_acr_max=np.max(all_kept_counts))
    coverage = pd.concat([coverage, pd.DataFrame([pooled])], ignore_index=True)
    rates = pd.DataFrame(rate_rows)
    for q in THRESHOLDS:
        rows = rates.loc[rates.threshold.eq(q)]
        denominator = rows.n_kept_acr_trials.sum()
        rates = pd.concat([rates, pd.DataFrame([{'rec_id': 'pooled', 'animal': 'all', 'threshold': q,
            'n_kept_acr_trials': denominator, 'n_observed_responding': rows.n_observed_responding.sum(),
            'shuffle_mean_n_responding': rows.shuffle_mean_n_responding.sum(),
            'observed_rate': rows.n_observed_responding.sum() / denominator,
            'shuffle_mean_rate': rows.shuffle_mean_n_responding.sum() / denominator}])], ignore_index=True)
    up = pd.DataFrame(up_rows)
    total_acr = up.n_valid_acr.sum()
    up = pd.concat([up, pd.DataFrame([{'rec_id': 'pooled', 'animal': 'all', 'n_valid_acr': total_acr,
        'n_saved_up': up.n_saved_up.sum(), 'n_cpu_up': up.n_cpu_up.sum(),
        'saved_up_pct': 100 * up.n_saved_up.sum() / total_acr,
        'cpu_up_pct': 100 * up.n_cpu_up.sum() / total_acr,
        'cpu_minus_saved_pp': 100 * (up.n_cpu_up.sum() - up.n_saved_up.sum()) / total_acr,
        'mean_recording_difference_pp': up.cpu_minus_saved_pp.mean()}])], ignore_index=True)
    response_count_null, overlap = pd.DataFrame(response_count_rows), pd.DataFrame(overlap_rows)
    variance = pd.DataFrame(variance_rows)
    test_rows = []
    for q in THRESHOLDS:
        test_rows.append(paired_test(response_count_null.loc[response_count_null.threshold.eq(q)], 'observed_pct',
                                     'null_median_pct', 'measure 2, leave-one-out count null', q))
        for label in ('permutation', 'shuffle'):
            test_rows.append(paired_test(overlap.loc[overlap.threshold.eq(q)], 'observed_overlap',
                                         label + '_null_median', 'measure 4, ' + label + ' null', q))
    tests = pd.DataFrame(test_rows)
    tables = counts, coverage, rates, up, response_count_null, overlap, variance, tests
    names = ('counts', 'coverage', 'responding_rates', 'up_classification', 'response_count_null',
             'overlap', 'response_count_variance_95', 'tests')
    for table, name in zip(tables, names):
        table.to_csv(OUT / ('aligned_check_' + name + '.csv'), index=False)
    pd.DataFrame(null_rows).to_csv(OUT / 'aligned_check_null_draws.csv', index=False)
    write_summary(tables, issues)
    print(tests.to_string(index=False), flush=True)


if __name__ == '__main__':
    # these small matrix products do not benefit from a large BLAS thread pool
    with threadpool_limits(limits=1):
        main()
