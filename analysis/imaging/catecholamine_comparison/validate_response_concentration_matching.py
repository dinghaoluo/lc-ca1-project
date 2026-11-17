# -*- coding: utf-8 -*-
'''
Created on Sat Aug 29 12:01:56 2026

validate the trial-, ROI- and reliability-matched response-concentration analysis

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

from _sha256 import sha256

import numpy as np
import pandas as pd
import scipy
from scipy.stats import mannwhitneyu, rankdata, wilcoxon


#%% paths and parameters
repo_root   = Path(__file__).resolve().parents[3]
INPUT_ROOT  = repo_root / 'data' / 'catecholamine_comparison'
OUTPUT_ROOT = repo_root / 'data' / 'catecholamine_comparison' / 'evidence'

SESSION_TABLE  = INPUT_ROOT / 'roi_release_metrics.csv'
SESSION_OUTPUT = OUTPUT_ROOT / 'full_response_concentration_matching_sessions.csv'
ANIMAL_OUTPUT  = OUTPUT_ROOT / 'full_response_concentration_matching_animals.csv'
RESULT_OUTPUT  = OUTPUT_ROOT / 'full_response_concentration_matching_comparisons.csv'
METHOD_OUTPUT  = OUTPUT_ROOT / 'full_response_concentration_matching_methods.csv'

PROCESSED_ROOTS = {
    'dLight': repo_root / 'data' / 'hpc_dlight' / 'lc_opto' / 'all_sessions',
    'nLight': repo_root / 'data' / 'hpc_nlight' / 'lc_opto' / 'all_sessions',
    }

N_TRIALS_MATCHED = 10
N_ROI_MATCHED    = 10
N_RESAMPLES      = 1000
N_PERMUTATIONS   = 50
N_BOOTSTRAPS     = 5000
RESAMPLE_SEED    = 20260829
BOOTSTRAP_SEED   = 20260830

METRICS = {
    'top10_excess_median': 'primary concentration endpoint',
    'split_half_excess_median': 'matched reliability diagnostic',
    }


#%% NumPy pickle compatibility
np_core = np.core
sys.modules.setdefault('numpy._core', np_core)
sys.modules.setdefault('numpy._core.multiarray', np_core.multiarray)


#%% metric helpers
def top_fraction_contribution(values, fraction=0.1):
    values = np.asarray(values, dtype=float)
    values = np.clip(values[np.isfinite(values)], 0, None)
    if values.size == 0 or np.sum(values) <= 0:
        return np.nan
    n_top = max(1, int(np.ceil(values.size * fraction)))
    return float(np.sum(np.sort(values)[-n_top:]) / np.sum(values))


def spearman_r(first, second):
    valid = np.isfinite(first) & np.isfinite(second)
    if np.sum(valid) < 3:
        return np.nan
    first_rank = rankdata(np.asarray(first)[valid], method='average')
    second_rank = rankdata(np.asarray(second)[valid], method='average')
    if np.ptp(first_rank) == 0 or np.ptp(second_rank) == 0:
        return np.nan
    return float(np.corrcoef(first_rank, second_rank)[0, 1])


def shuffled_spearman_median(first, second, rng):
    first_rank = rankdata(first, method='average')
    second_rank = rankdata(second, method='average')
    first_rank -= np.mean(first_rank)
    second_rank -= np.mean(second_rank)
    denominator = np.sqrt(
        np.sum(first_rank ** 2) * np.sum(second_rank ** 2)
        )
    if denominator == 0:
        return np.nan
    orders = np.argsort(
        rng.random((N_PERMUTATIONS, second.size)), axis=1
        )
    correlations = second_rank[orders] @ first_rank / denominator
    return float(np.median(correlations))


def session_seed(sensor, recname):
    token = f'{RESAMPLE_SEED}:{sensor}:{recname}'.encode('utf-8')
    return int.from_bytes(sha256(token).digest()[:4], 'little')


def load_roi_trial_matrix(sensor, recname):
    processed = PROCESSED_ROOTS[sensor] / recname / 'processed_data'
    response_path = processed / f'{recname}_pixel_RI_stim.npy'
    roi_path = processed / f'{recname}_ROI_dict.npy'
    response = np.load(response_path, mmap_mode='r')
    roi_dict = np.load(roi_path, allow_pickle=True).item()
    roi_trial = np.full((len(roi_dict), response.shape[2]), np.nan)
    for row, roi in enumerate(roi_dict.values()):
        roi_trial[row, :] = np.nanmean(
            response[roi['ypix'], roi['xpix'], :], axis=0
            )
    return roi_trial, response_path, roi_path


def resample_session(roi_trial, rng):
    results = {
        'top10_observed': np.full(N_RESAMPLES, np.nan),
        'top10_shuffle': np.full(N_RESAMPLES, np.nan),
        'top10_excess': np.full(N_RESAMPLES, np.nan),
        'split_half_observed': np.full(N_RESAMPLES, np.nan),
        'split_half_shuffle': np.full(N_RESAMPLES, np.nan),
        'split_half_excess': np.full(N_RESAMPLES, np.nan),
        }

    for index in range(N_RESAMPLES):
        roi_indices = rng.choice(
            roi_trial.shape[0], size=N_ROI_MATCHED, replace=False
            )
        trial_indices = rng.choice(
            roi_trial.shape[1], size=N_TRIALS_MATCHED, replace=False
            )
        sample = roi_trial[np.ix_(roi_indices, trial_indices)]

        observed = top_fraction_contribution(np.mean(sample, axis=1))
        trial_roi = sample.T
        orders = np.argsort(
            rng.random((N_PERMUTATIONS, N_TRIALS_MATCHED, N_ROI_MATCHED)),
            axis=2
            )
        spatial_shuffles = np.take_along_axis(
            np.broadcast_to(trial_roi, orders.shape), orders, axis=2
            )
        shuffled_means = np.mean(spatial_shuffles, axis=1)
        shuffled_positive = np.clip(shuffled_means, 0, None)
        shuffled_totals = np.sum(shuffled_positive, axis=1)
        n_top = max(1, int(np.ceil(N_ROI_MATCHED * .1)))
        shuffled_contributions = np.full(N_PERMUTATIONS, np.nan)
        positive_total = shuffled_totals > 0
        shuffled_contributions[positive_total] = (
            np.sum(
                np.sort(shuffled_positive[positive_total], axis=1)[:, -n_top:],
                axis=1
                )
            / shuffled_totals[positive_total]
            )
        shuffled = float(np.nanmedian(shuffled_contributions))
        results['top10_observed'][index] = observed
        results['top10_shuffle'][index] = shuffled
        if np.isfinite(observed) and np.isfinite(shuffled):
            results['top10_excess'][index] = observed - shuffled

        trial_order = rng.permutation(N_TRIALS_MATCHED)
        first_half = np.mean(sample[:, trial_order[:5]], axis=1)
        second_half = np.mean(sample[:, trial_order[5:]], axis=1)
        observed = spearman_r(first_half, second_half)
        shuffled = shuffled_spearman_median(first_half, second_half, rng)
        results['split_half_observed'][index] = observed
        results['split_half_shuffle'][index] = shuffled
        if np.isfinite(observed) and np.isfinite(shuffled):
            results['split_half_excess'][index] = observed - shuffled

    return results


def finite_summary(values):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0:
        return {'n': 0, 'median': np.nan, 'q1': np.nan, 'q3': np.nan}
    return {
        'n': int(values.size),
        'median': float(np.median(values)),
        'q1': float(np.quantile(values, 0.25)),
        'q3': float(np.quantile(values, 0.75)),
        }


def summarise_resamples(resamples):
    row = {}
    for metric, values in resamples.items():
        summary = finite_summary(values)
        for field, value in summary.items():
            row[f'{metric}_{field}'] = value
    return row


#%% session estimates
def build_session_table():
    source = pd.read_csv(SESSION_TABLE)
    if source.duplicated(['sensor', 'recname']).any():
        raise ValueError('roi_release_metrics.csv contains duplicated sessions')

    rows = []
    for source_row in source.itertuples(index=False):
        roi_trial, response_path, roi_path = load_roi_trial_matrix(
            source_row.sensor, source_row.recname
            )
        loaded_n_roi, loaded_n_trials = roi_trial.shape
        if loaded_n_roi != int(source_row.n_roi):
            raise ValueError(f'ROI count mismatch for {source_row.recname}')
        if loaded_n_trials != int(source_row.n_trials):
            raise ValueError(f'trial count mismatch for {source_row.recname}')

        reasons = []
        if loaded_n_roi < N_ROI_MATCHED:
            reasons.append(f'fewer than {N_ROI_MATCHED} ROIs')
        if loaded_n_trials < N_TRIALS_MATCHED:
            reasons.append(f'fewer than {N_TRIALS_MATCHED} trials')
        if not np.all(np.isfinite(roi_trial)):
            reasons.append('non-finite ROI-by-trial values')

        seed = session_seed(source_row.sensor, source_row.recname)
        exclusion_reason = '; '.join(reasons)
        output_row = {
            'transmitter': source_row.transmitter,
            'sensor': source_row.sensor,
            'recname': source_row.recname,
            'animal': source_row.animal,
            'source_n_roi': loaded_n_roi,
            'source_n_trials': loaded_n_trials,
            'matched_n_roi': N_ROI_MATCHED,
            'matched_n_trials': N_TRIALS_MATCHED,
            'n_resamples_requested': N_RESAMPLES,
            'session_seed': seed,
            'included_matched': len(reasons) == 0,
            'exclusion_reason': exclusion_reason,
            'full_session_top10_excess': source_row.top10_minus_shuffle,
            'full_session_split_half_excess': source_row.split_half_minus_shuffle,
            'response_source': str(response_path.relative_to(repo_root)),
            'roi_source': str(roi_path.relative_to(repo_root)),
            }
        if not reasons:
            output_row.update(summarise_resamples(
                resample_session(roi_trial, np.random.default_rng(seed))
                ))
        rows.append(output_row)
        status = 'included' if not reasons else f'excluded: {exclusion_reason}'
        print(f'{source_row.sensor} {source_row.recname}: {status}')

    sessions = pd.DataFrame(rows)
    eligible = sessions[sessions['included_matched']].copy()
    sensor_ranges = eligible.groupby('sensor')['full_session_split_half_excess'].agg(
        ['min', 'max']
        )
    overlap_low = float(sensor_ranges['min'].max())
    overlap_high = float(sensor_ranges['max'].min())
    if overlap_low > overlap_high:
        raise ValueError('the sensors have no overlapping full-session reliability')

    sessions['reliability_overlap_low'] = overlap_low
    sessions['reliability_overlap_high'] = overlap_high
    sessions['included_reliability_overlap'] = (
        sessions['included_matched']
        & sessions['full_session_split_half_excess'].between(
            overlap_low, overlap_high, inclusive='both'
            )
        )
    return sessions, overlap_low, overlap_high


#%% sensor comparisons
def values_by_animal(frame, metric):
    return [group[metric].dropna().to_numpy(dtype=float)
            for _, group in frame.groupby('animal')]


def resample_sensor_median(groups, rng):
    sampled_animals = rng.integers(0, len(groups), size=len(groups))
    animal_medians = [
        np.median(rng.choice(groups[index], size=len(groups[index]), replace=True))
        for index in sampled_animals
        ]
    return float(np.median(animal_medians))


def hierarchical_difference(frame, metric, seed):
    rng = np.random.default_rng(seed)
    dlight = values_by_animal(frame[frame['sensor'] == 'dLight'], metric)
    nlight = values_by_animal(frame[frame['sensor'] == 'nLight'], metric)
    differences = np.empty(N_BOOTSTRAPS)
    for index in range(N_BOOTSTRAPS):
        differences[index] = (
            resample_sensor_median(nlight, rng)
            - resample_sensor_median(dlight, rng)
            )
    return {
        'median': float(np.median(differences)),
        'ci_low': float(np.quantile(differences, 0.025)),
        'ci_high': float(np.quantile(differences, 0.975)),
        }


def within_sensor_test(values):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size == 0 or np.all(values == 0):
        return np.nan, np.nan
    result = wilcoxon(values, alternative='greater', method='approx')
    return float(result.statistic), float(result.pvalue)


def compare_sensors(frame, subset, metric, role, seed):
    current = frame[['sensor', 'animal', 'recname', metric]].dropna().copy()
    nlight = current.loc[current['sensor'] == 'nLight', metric].to_numpy(float)
    dlight = current.loc[current['sensor'] == 'dLight', metric].to_numpy(float)
    n_summary = finite_summary(nlight)
    d_summary = finite_summary(dlight)

    session_test = mannwhitneyu(
        nlight, dlight, alternative='two-sided', method='auto'
        )
    session_effect = (
        2 * session_test.statistic / (len(nlight) * len(dlight)) - 1
        )
    n_w, n_w_p = within_sensor_test(nlight)
    d_w, d_w_p = within_sensor_test(dlight)

    animal_values = current.groupby(['sensor', 'animal'], as_index=False)[
        metric
        ].median()
    n_animals = animal_values.loc[
        animal_values['sensor'] == 'nLight', metric
        ].to_numpy(float)
    d_animals = animal_values.loc[
        animal_values['sensor'] == 'dLight', metric
        ].to_numpy(float)
    animal_test = mannwhitneyu(
        n_animals, d_animals, alternative='two-sided', method='exact'
        )
    bootstrap = hierarchical_difference(current, metric, seed)

    return {
        'subset': subset,
        'metric': metric,
        'analysis_role': role,
        'nlight_n_sessions': n_summary['n'],
        'nlight_n_animals': len(n_animals),
        'nlight_median': n_summary['median'],
        'nlight_q1': n_summary['q1'],
        'nlight_q3': n_summary['q3'],
        'nlight_within_sensor_w': n_w,
        'nlight_within_sensor_one_sided_p': n_w_p,
        'dlight_n_sessions': d_summary['n'],
        'dlight_n_animals': len(d_animals),
        'dlight_median': d_summary['median'],
        'dlight_q1': d_summary['q1'],
        'dlight_q3': d_summary['q3'],
        'dlight_within_sensor_w': d_w,
        'dlight_within_sensor_one_sided_p': d_w_p,
        'session_median_difference_nlight_minus_dlight': (
            n_summary['median'] - d_summary['median']
            ),
        'session_mann_whitney_u_nlight': float(session_test.statistic),
        'session_two_sided_p': float(session_test.pvalue),
        'session_rank_biserial_nlight_minus_dlight': float(session_effect),
        'nlight_animal_median': float(np.median(n_animals)),
        'dlight_animal_median': float(np.median(d_animals)),
        'animal_median_difference_nlight_minus_dlight': float(
            np.median(n_animals) - np.median(d_animals)
            ),
        'animal_mann_whitney_u_nlight': float(animal_test.statistic),
        'animal_exact_two_sided_p': float(animal_test.pvalue),
        'hierarchical_bootstrap_difference_median': bootstrap['median'],
        'hierarchical_bootstrap_ci_low': bootstrap['ci_low'],
        'hierarchical_bootstrap_ci_high': bootstrap['ci_high'],
        'bootstrap_seed': seed,
        'n_bootstraps': N_BOOTSTRAPS,
        }


def build_comparison_tables(sessions):
    subsets = {
        'all_dimension_matched_sessions': sessions['included_matched'],
        'overlapping_full_session_reliability': sessions[
            'included_reliability_overlap'
            ],
        }
    comparison_rows = []
    animal_rows = []
    row_index = 0
    for subset, mask in subsets.items():
        current = sessions[mask].copy()
        for metric, role in METRICS.items():
            comparison_rows.append(compare_sensors(
                current, subset, metric, role, BOOTSTRAP_SEED + row_index
                ))
            row_index += 1

        for (sensor, animal), group in current.groupby(['sensor', 'animal']):
            animal_rows.append({
                'subset': subset,
                'sensor': sensor,
                'animal': animal,
                'n_sessions': len(group),
                'top10_excess_animal_median': group[
                    'top10_excess_median'
                    ].median(),
                'split_half_excess_animal_median': group[
                    'split_half_excess_median'
                    ].median(),
                'full_session_split_half_excess_animal_median': group[
                    'full_session_split_half_excess'
                    ].median(),
                })
    return pd.DataFrame(animal_rows), pd.DataFrame(comparison_rows)


#%% methods and output
def build_method_table(sessions, overlap_low, overlap_high):
    included = sessions[sessions['included_matched']]
    overlap = sessions[sessions['included_reliability_overlap']]
    methods = {
        'session_input': str(SESSION_TABLE.relative_to(repo_root)),
        'response_input_pattern': (
            'data/hpc_<sensor>/lc_opto/all_sessions/<session>/processed_data/'
            '<session>_pixel_RI_stim.npy'
            ),
        'roi_input_pattern': (
            'data/hpc_<sensor>/lc_opto/all_sessions/<session>/processed_data/'
            '<session>_ROI_dict.npy'
            ),
        'matched_trials': N_TRIALS_MATCHED,
        'matched_rois': N_ROI_MATCHED,
        'resamples_per_session': N_RESAMPLES,
        'null_permutations_per_resample': N_PERMUTATIONS,
        'resample_seed': RESAMPLE_SEED,
        'top10_definition': (
            'largest positive ROI mean divided by all positive ROI means '
            'after matching 10 ROIs'
            ),
        'top10_null': (
            'independent ROI-identity permutation within each sampled trial'
            ),
        'top10_session_endpoint': (
            'median paired observed-minus-permuted top-10% contribution'
            ),
        'split_half_definition': (
            'Spearman correlation between two random five-trial mean maps'
            ),
        'split_half_null': 'permuted ROI identity in the second half map',
        'split_half_session_endpoint': (
            'median paired observed-minus-permuted split-half correlation'
            ),
        'full_session_reliability_source': (
            'split_half_minus_shuffle in roi_release_metrics.csv'
            ),
        'reliability_overlap_rule': (
            'intersection of sensor ranges after trial and ROI eligibility'
            ),
        'reliability_overlap_low': overlap_low,
        'reliability_overlap_high': overlap_high,
        'bootstrap_resamples': N_BOOTSTRAPS,
        'bootstrap_seed_start': BOOTSTRAP_SEED,
        'within_sensor_test': 'one-sided Wilcoxon normal approximation',
        'dlight_dimension_matched_sessions': int(
            np.sum(included['sensor'] == 'dLight')
            ),
        'nlight_dimension_matched_sessions': int(
            np.sum(included['sensor'] == 'nLight')
            ),
        'dlight_reliability_overlap_sessions': int(
            np.sum(overlap['sensor'] == 'dLight')
            ),
        'nlight_reliability_overlap_sessions': int(
            np.sum(overlap['sensor'] == 'nLight')
            ),
        'numpy_version': np.__version__,
        'pandas_version': pd.__version__,
        'scipy_version': scipy.__version__,
        }
    return pd.DataFrame([
        {'parameter': key, 'value': value} for key, value in methods.items()
        ])


def save_new(frame, path):
    frame.to_csv(path, index=False)


#%% run
def main():
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    sessions, overlap_low, overlap_high = build_session_table()
    animals, comparisons = build_comparison_tables(sessions)
    methods = build_method_table(sessions, overlap_low, overlap_high)

    save_new(sessions, SESSION_OUTPUT)
    save_new(animals, ANIMAL_OUTPUT)
    save_new(comparisons, RESULT_OUTPUT)
    save_new(methods, METHOD_OUTPUT)

    print('\nfull-session reliability overlap: '
          f'{overlap_low:.6f} to {overlap_high:.6f}')
    print(comparisons.to_string(index=False))
    for path in [SESSION_OUTPUT, ANIMAL_OUTPUT, RESULT_OUTPUT, METHOD_OUTPUT]:
        print(f'saved {path.relative_to(repo_root)}')


if __name__ == '__main__':
    main()
