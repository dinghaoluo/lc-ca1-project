# -*- coding: utf-8 -*-
'''
Created on Sat Aug 29 11:55:22 2026

recompute optogenetic response concentration after matching trials, ROIs and reliability

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

import zlib

import numpy as np

# The processed arrays were written with a newer NumPy module path.
import numpy.core as _np_core
sys.modules.setdefault('numpy._core', _np_core)
sys.modules.setdefault('numpy._core.multiarray', _np_core.multiarray)
sys.modules.setdefault('numpy._core.numeric', _np_core.numeric)
try:
    sys.modules.setdefault('numpy._core._multiarray_umath',
                           _np_core._multiarray_umath)
except:
    pass

import pandas as pd
from scipy.stats import mannwhitneyu, rankdata, spearmanr, wilcoxon

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from common_functions import mpl_formatting
mpl_formatting()


#%% paths and parameters
repo_root = Path(__file__).resolve().parents[3]
DATA_ROOT = repo_root / 'data'
COMPARISON_ROOT = DATA_ROOT / 'catecholamine_comparison'
INPUT_TABLE = COMPARISON_ROOT / 'roi_release_metrics.csv'
OUTPUT_ROOT = COMPARISON_ROOT / 'evidence'
FIGURE_ROOT = repo_root / 'figures' / 'catecholamine_comparison' / 'evidence'

SENSOR_STEMS = {
    'dLight': DATA_ROOT / 'hpc_dlight' / 'lc_opto' / 'all_sessions',
    'nLight': DATA_ROOT / 'hpc_nlight' / 'lc_opto' / 'all_sessions',
    }

N_MATCHED_TRIALS = 10
N_MATCHED_ROIS = 10
N_RESAMPLES = 1000
N_SHUFFLES = 50
N_BOOTSTRAP = 5000
RNG_SEED = 20260829


#%% metric functions
def top_fraction_contribution(values, fraction=0.1):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    values = np.clip(values, 0, None)
    if values.size == 0 or np.nansum(values) <= 0:
        return np.nan
    n_top = max(1, int(np.ceil(values.size * fraction)))
    return float(np.sum(np.sort(values)[-n_top:]) / np.sum(values))


def spearman_value(first, second):
    valid = np.isfinite(first) & np.isfinite(second)
    if np.sum(valid) < 3:
        return np.nan
    if np.nanstd(first[valid]) == 0 or np.nanstd(second[valid]) == 0:
        return np.nan
    return float(spearmanr(first[valid], second[valid]).statistic)


def spatial_shuffle_excess(matrix, rng):
    observed_mean = np.nanmean(matrix, axis=1)
    observed_top = top_fraction_contribution(observed_mean)

    trial_roi = matrix.T
    orders = np.argsort(
        rng.random((N_SHUFFLES, matrix.shape[1], matrix.shape[0])), axis=2)
    permuted = np.take_along_axis(
        np.broadcast_to(trial_roi, orders.shape), orders, axis=2)
    shuffled_means = np.nanmean(permuted, axis=1)
    shuffled_positive = np.clip(shuffled_means, 0, None)
    n_top = max(1, int(np.ceil(matrix.shape[0] * .1)))
    shuffled_totals = np.nansum(shuffled_positive, axis=1)
    shuffled_numerator = np.nansum(
        np.sort(shuffled_positive, axis=1)[:, -n_top:], axis=1)
    shuffled_top = np.full(N_SHUFFLES, np.nan)
    np.divide(shuffled_numerator, shuffled_totals, out=shuffled_top,
              where=shuffled_totals > 0)

    trial_order = rng.permutation(matrix.shape[1])
    half = matrix.shape[1] // 2
    odd = np.nanmean(matrix[:, trial_order[:half]], axis=1)
    even = np.nanmean(matrix[:, trial_order[half:]], axis=1)
    observed_split = spearman_value(odd, even)
    split_orders = np.argsort(
        rng.random((N_SHUFFLES, even.size)), axis=1)
    if np.all(np.isfinite(odd)) and np.all(np.isfinite(even)):
        odd_rank = rankdata(odd)
        even_rank = rankdata(even)
        odd_rank -= np.mean(odd_rank)
        even_rank -= np.mean(even_rank)
        denominator = np.sqrt(
            np.sum(odd_rank ** 2) * np.sum(even_rank ** 2))
        if denominator > 0:
            shuffled_split = (
                even_rank[split_orders] @ odd_rank / denominator)
        else:
            shuffled_split = np.full(N_SHUFFLES, np.nan)
    else:
        shuffled_split = np.asarray([
            spearman_value(odd, even[order]) for order in split_orders])

    return {
        'top10_excess': observed_top - np.nanmedian(shuffled_top),
        'top10_observed': observed_top,
        'top10_shuffle_median': float(np.nanmedian(shuffled_top)),
        'split_half_excess': observed_split - np.nanmedian(shuffled_split),
        'split_half_observed': observed_split,
        'split_half_shuffle_median': float(np.nanmedian(shuffled_split)),
        }


def load_trial_matrix(sensor, recname):
    stem = SENSOR_STEMS[sensor] / recname / 'processed_data'
    pixel_path = stem / f'{recname}_pixel_RI_stim.npy'
    roi_path = stem / f'{recname}_ROI_dict.npy'
    pixel_ri = np.load(pixel_path, allow_pickle=True)
    roi_dict = np.load(roi_path, allow_pickle=True).item()

    rows = []
    for roi_id in roi_dict:
        roi = roi_dict[roi_id]
        xpix = np.asarray(roi['xpix'], dtype=int)
        ypix = np.asarray(roi['ypix'], dtype=int)
        rows.append(np.nanmean(pixel_ri[ypix, xpix, :], axis=0))

    matrix = np.asarray(rows, dtype=float)
    valid = np.sum(np.isfinite(matrix), axis=1) >= 5
    matrix = matrix[valid]
    if matrix.shape[0] < N_MATCHED_ROIS or matrix.shape[1] < N_MATCHED_TRIALS:
        raise ValueError(f'{recname}: {matrix.shape} cannot meet matching targets')
    return matrix


def session_resamples(row):
    matrix = load_trial_matrix(row['sensor'], row['recname'])
    seed = RNG_SEED + zlib.crc32(row['recname'].encode('utf-8'))
    rng = np.random.default_rng(seed)
    records = []
    for resample in range(N_RESAMPLES):
        roi_idx = rng.choice(matrix.shape[0], N_MATCHED_ROIS, replace=False)
        trial_idx = rng.choice(matrix.shape[1], N_MATCHED_TRIALS, replace=False)
        matched = matrix[np.ix_(roi_idx, trial_idx)]
        metrics = spatial_shuffle_excess(matched, rng)
        records.append({
            'sensor': row['sensor'],
            'transmitter': row['transmitter'],
            'recname': row['recname'],
            'animal': row['animal'],
            'resample': resample,
            'n_original_rois': matrix.shape[0],
            'n_original_trials': matrix.shape[1],
            **metrics,
            })
    return records


#%% clustered summaries
def session_summary(resamples, source):
    grouped = resamples.groupby(
        ['sensor', 'transmitter', 'recname', 'animal'], as_index=False)
    summary = grouped.agg(
        matched_top10_excess=('top10_excess', 'median'),
        matched_split_half_excess=('split_half_excess', 'median'),
        matched_top10_observed=('top10_observed', 'median'),
        matched_top10_shuffle=('top10_shuffle_median', 'median'),
        matched_split_half_observed=('split_half_observed', 'median'),
        matched_split_half_shuffle=('split_half_shuffle_median', 'median'),
        n_resamples=('resample', 'size'),
        n_original_rois=('n_original_rois', 'first'),
        n_original_trials=('n_original_trials', 'first'),
        )
    source_cols = ['sensor', 'recname', 'split_half_minus_shuffle',
                   'top10_minus_shuffle', 'n_roi', 'n_trials']
    return summary.merge(source[source_cols], on=['sensor', 'recname'],
                         how='left', validate='one_to_one')


def clustered_bootstrap(frame, metric, seed):
    rng = np.random.default_rng(seed)
    values = {}
    for sensor in ['nLight', 'dLight']:
        curr = frame[frame['sensor'] == sensor]
        values[sensor] = {
            animal: group[metric].to_numpy(dtype=float)
            for animal, group in curr.groupby('animal')
            }
    if not values['nLight'] or not values['dLight']:
        return np.array([], dtype=float)

    output = np.empty(N_BOOTSTRAP, dtype=float)
    for index in range(N_BOOTSTRAP):
        medians = {}
        for sensor in ['nLight', 'dLight']:
            animals = list(values[sensor])
            sampled_animals = rng.choice(animals, len(animals), replace=True)
            sampled_medians = []
            for animal in sampled_animals:
                sessions = values[sensor][animal]
                sampled = rng.choice(sessions, len(sessions), replace=True)
                sampled_medians.append(np.nanmedian(sampled))
            medians[sensor] = sampled_medians
        output[index] = (np.nanmedian(medians['nLight'])
                         - np.nanmedian(medians['dLight']))
    return output


def compare_frame(frame, label):
    metric = 'matched_top10_excess'
    rows = []
    for sensor in ['nLight', 'dLight']:
        vals = frame.loc[frame['sensor'] == sensor, metric].dropna()
        within_test = wilcoxon(vals, alternative='greater')
        rows.append({
            'analysis': label,
            'sensor': sensor,
            'n_sessions': len(vals),
            'n_animals': frame.loc[frame['sensor'] == sensor, 'animal'].nunique(),
            'median': vals.median(),
            'q1': vals.quantile(.25),
            'q3': vals.quantile(.75),
            'within_sensor_w': float(within_test.statistic),
            'within_sensor_p_value': float(within_test.pvalue),
            })

    nlight = frame.loc[frame['sensor'] == 'nLight', metric].dropna()
    dlight = frame.loc[frame['sensor'] == 'dLight', metric].dropna()
    nlight_reliability = frame.loc[
        frame['sensor'] == 'nLight', 'matched_split_half_excess'].dropna()
    dlight_reliability = frame.loc[
        frame['sensor'] == 'dLight', 'matched_split_half_excess'].dropna()
    animal = (frame.groupby(['sensor', 'animal'], as_index=False)[metric]
              .median())
    nlight_animal = animal.loc[animal['sensor'] == 'nLight', metric]
    dlight_animal = animal.loc[animal['sensor'] == 'dLight', metric]
    if len(nlight) and len(dlight):
        session_test = mannwhitneyu(
            nlight, dlight, alternative='two-sided', method='auto')
        session_effect = (
            2 * session_test.statistic / (len(nlight) * len(dlight)) - 1)
        animal_test = mannwhitneyu(
            nlight_animal, dlight_animal, alternative='two-sided',
            method='exact')
        animal_effect = (
            2 * animal_test.statistic
            / (len(nlight_animal) * len(dlight_animal)) - 1)
        reliability_test = mannwhitneyu(
            nlight_reliability, dlight_reliability,
            alternative='two-sided', method='auto')
        boot = clustered_bootstrap(frame, metric, RNG_SEED + 11)
        rows.append({
            'analysis': label,
            'sensor': 'contrast_nLight_minus_dLight',
            'n_sessions': len(nlight) + len(dlight),
            'n_animals': frame['animal'].nunique(),
            'median': float(nlight.median() - dlight.median()),
            'q1': np.nan,
            'q3': np.nan,
            'mann_whitney_u': float(session_test.statistic),
            'p_value': float(session_test.pvalue),
            'rank_biserial': float(session_effect),
            'animal_mann_whitney_u': float(animal_test.statistic),
            'animal_p_value': float(animal_test.pvalue),
            'animal_rank_biserial': float(animal_effect),
            'nlight_reliability_median': float(nlight_reliability.median()),
            'dlight_reliability_median': float(dlight_reliability.median()),
            'reliability_mann_whitney_u': float(reliability_test.statistic),
            'reliability_p_value': float(reliability_test.pvalue),
            'bootstrap_median': float(np.nanmedian(boot)),
            'bootstrap_ci_low': float(np.nanquantile(boot, .025)),
            'bootstrap_ci_high': float(np.nanquantile(boot, .975)),
            })
    return rows


#%% figure
def plot_summary(summary, overlap_low, overlap_high):
    colours = {'dLight': '#4f8f62', 'nLight': '#8b5ea7'}
    fig, axes = plt.subplots(1, 2, figsize=(8.0, 3.4))

    ax = axes[0]
    for xpos, sensor in enumerate(['dLight', 'nLight'], start=1):
        current = summary[summary['sensor'] == sensor]
        values = current['matched_top10_excess'].dropna().to_numpy()
        rng = np.random.default_rng(RNG_SEED + xpos)
        jitter = rng.normal(0, .045, values.size)
        ax.scatter(np.full(values.size, xpos) + jitter, values,
                   s=18, alpha=.55, color=colours[sensor], edgecolor='none')
        animal_values = current.groupby('animal')[
            'matched_top10_excess'].median().to_numpy()
        ax.scatter(np.full(animal_values.size, xpos), animal_values,
                   s=38, facecolor='white', edgecolor=colours[sensor],
                   marker='o', linewidth=1.2, zorder=4)
        ax.plot([xpos - .18, xpos + .18], [np.nanmedian(values)] * 2,
                color='black', linewidth=2, zorder=5)
    ax.set_xticks([1, 2])
    ax.set_xticklabels(['dLight', 'nLight'])
    ax.set_ylabel('Matched top-10% excess contribution')
    ax.set_title('10 trials, 10 ROIs')
    ax.spines[['top', 'right']].set_visible(False)

    ax = axes[1]
    for sensor in ['dLight', 'nLight']:
        current = summary[summary['sensor'] == sensor]
        ax.scatter(current['matched_split_half_excess'],
                   current['matched_top10_excess'], s=20,
                   color=colours[sensor], alpha=.65, label=sensor,
                   edgecolor='none')
    ax.axvspan(overlap_low, overlap_high, color='grey', alpha=.15,
               label='Reliability overlap')
    ax.set_xlabel('Matched split-half excess')
    ax.set_ylabel('Matched top-10% excess')
    ax.legend(frameon=False, fontsize=7)
    ax.spines[['top', 'right']].set_visible(False)

    fig.tight_layout()
    FIGURE_ROOT.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGURE_ROOT / 'response_concentration_matching.png', dpi=300,
                bbox_inches='tight')
    fig.savefig(FIGURE_ROOT / 'response_concentration_matching.svg',
                bbox_inches='tight')
    plt.close(fig)


#%% main
def main():
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    source = pd.read_csv(INPUT_TABLE)
    source = source[source['sensor'].isin(SENSOR_STEMS)].copy()
    source['matched_analysis_status'] = 'included'
    source.loc[source['n_roi'] < N_MATCHED_ROIS,
               'matched_analysis_status'] = 'fewer than 10 ROIs'
    source.loc[source['n_trials'] < N_MATCHED_TRIALS,
               'matched_analysis_status'] = 'fewer than 10 trials'
    both_low = ((source['n_roi'] < N_MATCHED_ROIS)
                & (source['n_trials'] < N_MATCHED_TRIALS))
    source.loc[both_low, 'matched_analysis_status'] = (
        'fewer than 10 ROIs and 10 trials')
    eligible = source[source['matched_analysis_status'] == 'included'].copy()

    all_records = []
    for _, row in eligible.iterrows():
        print(row['recname'], flush=True)
        all_records.extend(session_resamples(row))
    resamples = pd.DataFrame(all_records)
    summary = session_summary(resamples, eligible)

    overlap_low = max(
        summary.loc[summary['sensor'] == 'dLight',
                    'matched_split_half_excess'].quantile(.05),
        summary.loc[summary['sensor'] == 'nLight',
                    'matched_split_half_excess'].quantile(.05))
    overlap_high = min(
        summary.loc[summary['sensor'] == 'dLight',
                    'matched_split_half_excess'].quantile(.95),
        summary.loc[summary['sensor'] == 'nLight',
                    'matched_split_half_excess'].quantile(.95))
    summary['matched_reliability_overlap'] = summary[
        'matched_split_half_excess'].between(overlap_low, overlap_high)

    overlap = summary[summary['matched_reliability_overlap']].copy()
    result_rows = compare_frame(summary, 'trial- and ROI-matched')
    result_rows.extend(compare_frame(
        overlap, 'trial-, ROI- and reliability-overlap-matched'))
    results = pd.DataFrame(result_rows)
    results['matched_trials'] = N_MATCHED_TRIALS
    results['matched_rois'] = N_MATCHED_ROIS
    results['n_resamples_per_session'] = N_RESAMPLES
    results['n_shuffle_patterns_per_resample'] = N_SHUFFLES
    results['reliability_overlap_low'] = overlap_low
    results['reliability_overlap_high'] = overlap_high

    animals = (summary.groupby(['sensor', 'animal'], as_index=False)
               .agg(matched_top10_excess=('matched_top10_excess', 'median'),
                    matched_split_half_excess=(
                        'matched_split_half_excess', 'median'),
                    n_sessions=('recname', 'size')))

    source.to_csv(
        OUTPUT_ROOT / 'response_concentration_matching_inclusion.csv', index=False)
    resamples.to_csv(
        OUTPUT_ROOT / 'response_concentration_matching_resamples.csv', index=False)
    summary.to_csv(
        OUTPUT_ROOT / 'response_concentration_matching_sessions.csv', index=False)
    animals.to_csv(
        OUTPUT_ROOT / 'response_concentration_matching_animals.csv', index=False)
    results.to_csv(
        OUTPUT_ROOT / 'response_concentration_matching_summary.csv', index=False)
    plot_summary(summary, overlap_low, overlap_high)

    methods = f'''# Dimension-matched response-concentration analysis

The analysis recomputes the optogenetic top-10% contribution-minus-shuffle metric
from the stored trial-by-pixel arrays and ROI masks listed in
`roi_release_metrics.csv`. Sessions require at least {N_MATCHED_TRIALS} trials
and {N_MATCHED_ROIS} ROIs, which preserves all animals in both cohorts. Each
eligible session is subsampled without replacement to exactly
{N_MATCHED_TRIALS} trials and {N_MATCHED_ROIS} ROIs, with {N_RESAMPLES}
independent subsamples. Because exactly ten ROIs are retained, the largest
single ROI is the top 10%. Each subsample uses {N_SHUFFLES} within-session
ROI-identity shuffles per trial, preserving the recorded trial values, and
reports top-10% excess and split-half Spearman excess above the median shuffle.

The reliability-overlap analysis retains session medians between
`{overlap_low:.6f}` and `{overlap_high:.6f}`, the intersection of the sensors'
5th-to-95th-percentile matched split-half distributions. It therefore tests
the concentration contrast only where measured reliability overlaps; the
summary also reports the residual reliability contrast after restriction.
Session values are the median across subsamples. Animal medians, an exact
animal-level Mann-Whitney sensitivity test and an animal-then-session
bootstrap interval accompany the descriptive session-level comparison.

Random seed: `{RNG_SEED}`. The complete resample table is retained so every
session median can be reconstructed. This remains a sensitivity analysis,
not a replacement for an animal-balanced experiment.
'''
    (OUTPUT_ROOT / 'response_concentration_matching_methods.md').write_text(
        methods, encoding='utf-8')

    print(results.to_string(index=False))
    print(f'eligible sessions: {len(summary)}')
    print(f'reliability-overlap sessions: {len(overlap)}')


if __name__ == '__main__':
    main()
