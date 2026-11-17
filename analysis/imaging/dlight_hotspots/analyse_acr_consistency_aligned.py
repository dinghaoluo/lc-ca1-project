# -*- coding: utf-8 -*-
'''
Created on Wed Nov 11 15:14:01 2026

measure aligned single-trial dLight response consistency

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

import importlib

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from common_functions import mpl_formatting
mpl_formatting()
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

repo_root = Path(__file__).resolve().parents[3]
JINGYU_ROOT = pp.JINGYU_CODE_ROOT
RAW_ROOT = pp.DBH_DLIGHT_RAW_ROOT
OUTPUT_ROOT = pp.DLIGHT_HOTSPOTS_OUTPUT_ROOT
FIGURE_ROOT = pp.DLIGHT_HOTSPOTS_FIGURE_ROOT
if str(JINGYU_ROOT) not in sys.path:
    sys.path.insert(0, str(JINGYU_ROOT))
recnames = rec_list.recDLightDbh
from common.utils_basic import trace_filter  # noqa: E402
from grid_response_test import postprocess_dff, quantify_grid_response  # noqa: E402

THRESHOLDS = (90, 95, 99)
N_SHUFFLE = 1000
CHANCE_SEED = 42


#%% loaders and event frames
def load_behaviour(path):
    numpy_core = importlib.import_module('numpy.core')
    sys.modules.setdefault('numpy._core', numpy_core)
    for name in ('multiarray', 'numeric', '_multiarray_umath'):
        try:
            module = importlib.import_module(f'numpy.core.{name}')
        except ImportError:
            continue
        sys.modules.setdefault(f'numpy._core.{name}', module)
    return pd.read_pickle(path)


def load_recording(recname):
    dff_path = (
        RAW_ROOT / 'regression_res_grid_free_dilation' / recname
        / 'single_trial_regression' / 'dilation_k=0'
        / 'dff_corrected_dlight.npy'
    )
    profile_path = (
        RAW_ROOT / 'processed_dataframe_grid_free_dilation'
        / f'{recname}_profile_combined_dilation=0_pre(-1, 0)_post(0, 1)'
        '_ES=0.05_shuff95_0621.parquet'
    )
    behaviour = load_behaviour(RAW_ROOT / 'behaviour_profile' / f'{recname}.pkl')
    dff = postprocess_dff(np.load(dff_path, mmap_mode='r'), trace_filter)
    profile = pd.read_parquet(profile_path)
    valid = profile['dlight_valid'].astype(bool) & profile['red_valid'].astype(bool) & ~profile['edge'].astype(bool)
    profile = profile.loc[valid].copy()
    coords = np.asarray(profile['roi_id'].tolist(), dtype=int)
    return behaviour, profile, dff[coords[:, 0], coords[:, 1], :]


def event_frames(behaviour):
    run = np.asarray(behaviour['run_onset_frames'], dtype=float)
    reward = np.asarray(behaviour['reward_times'], dtype=float)
    non_stop = np.asarray(behaviour['non_stop_trials'])
    non_fullstop = np.asarray(behaviour['non_fullstop_trials'])
    n = min(len(run), len(reward), len(non_stop), len(non_fullstop))
    valid = (
        np.isfinite(run[:n])
        & (run[:n] != -1)
        & np.isfinite(reward[:n])
        & (non_stop[:n] == 0)
        & (non_fullstop[:n] == 0)
    )
    return run[:n][valid].astype(int)


#%% response consistency
def response_matrices(result, percentile):
    shuffle_single = result['shuffle_single_trial_amplitude']
    threshold = np.nanpercentile(
        shuffle_single.reshape(shuffle_single.shape[0], -1), percentile, axis=1)
    pooled_std = result['pooled_std']
    real_effect = result['trial_amplitude'] / pooled_std[:, None]
    real_response = (
        (result['trial_amplitude'] > threshold[:, None])
        & (real_effect > 0.05)
        & result['keep']
    )
    shuffle_effect = shuffle_single / pooled_std[:, None, None]
    shuffle_response = (
        (shuffle_single > threshold[:, None, None])
        & (shuffle_effect > 0.05)
        & result['keep'][:, :, None]
    )
    return threshold, real_response, shuffle_response


def jaccard_summary(real_response, shuffle_response, keep):
    real_sets = (real_response & keep).T.astype(np.float64)
    available = keep.T.astype(np.float64)
    intersection = real_sets @ real_sets.T
    counts_on_common = real_sets @ available.T
    union = counts_on_common + counts_on_common.T - intersection
    pair_mask = np.triu(union > 0, k=1)
    observed_values = intersection[pair_mask] / union[pair_mask]
    observed = float(np.mean(observed_values)) if len(observed_values) else np.nan
    null_per_shuffle = np.full(shuffle_response.shape[-1], np.nan)
    for start in range(0, shuffle_response.shape[-1], 32):
        stop = min(start + 32, shuffle_response.shape[-1])
        sets = (shuffle_response[:, :, start:stop] & keep[:, :, None])
        sets = sets.transpose(2, 1, 0).astype(np.float64)
        null_intersection = sets @ sets.transpose(0, 2, 1)
        null_counts_on_common = sets @ available.T
        null_union = (null_counts_on_common
                      + null_counts_on_common.transpose(0, 2, 1)
                      - null_intersection)
        for offset in range(stop - start):
            mask = np.triu(null_union[offset] > 0, k=1)
            if np.any(mask):
                null_per_shuffle[start + offset] = np.mean(
                    null_intersection[offset][mask] / null_union[offset][mask])
    return observed, null_per_shuffle, int(pair_mask.sum())


def group_trial_fraction(real_response, keep, group):
    values = real_response & keep
    if group is not None:
        values = values[group]
        keep = keep[group]
    fractions = []
    for trial in range(values.shape[1]):
        valid = keep[:, trial]
        if np.any(valid):
            fractions.append(float(values[valid, trial].mean()))
    return np.asarray(fractions, dtype=float)


def safe_wilcoxon_pair(first, second):
    first = np.asarray(first, dtype=float)
    second = np.asarray(second, dtype=float)
    valid = np.isfinite(first) & np.isfinite(second)
    first, second = first[valid], second[valid]
    if not len(first):
        return np.nan, np.nan, 0
    if np.allclose(first - second, 0):
        return 0.0, 1.0, len(first)
    result = wilcoxon(first, second, alternative='two-sided')
    return float(result.statistic), float(result.pvalue), len(first)


def safe_wilcoxon_constant(values, constant):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if not len(values):
        return np.nan, np.nan, 0
    if np.allclose(values - constant, 0):
        return 0.0, 1.0, len(values)
    result = wilcoxon(values - constant, alternative='two-sided')
    return float(result.statistic), float(result.pvalue), len(values)


#%% figures
def ecdf(values, grid):
    values = np.sort(np.asarray(values, dtype=float))
    return np.searchsorted(values, grid, side='right') / len(values)


def make_probability_figure(observed, null):
    observed = np.asarray(observed, dtype=float)
    null = np.asarray(null, dtype=float)
    grid = np.linspace(0, 1, 201)
    observed_ecdf = ecdf(observed, grid)
    null_ecdf = np.vstack([ecdf(null[:, i], grid) for i in range(null.shape[1])])
    median = np.nanmedian(null_ecdf, axis=0)
    low, high = np.nanpercentile(null_ecdf, [2.5, 97.5], axis=0)
    fig, ax = plt.subplots(figsize=(4.6, 3.8))
    ax.fill_between(grid, low, high, color='0.75', alpha=0.55, label='shuffle, 95% band')
    ax.plot(grid, median, color='0.35', linewidth=1.5, label='shuffle median')
    ax.plot(grid, observed_ecdf, color='#19a819', linewidth=1.5, label='observed')
    ax.set_xlabel('response probability per ACR')
    ax.set_ylabel('cumulative fraction of ACRs')
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.spines[['top', 'right']].set_visible(False)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIGURE_ROOT / 'single_trial_response_probability_95.png', dpi=200, facecolor='white')
    plt.close(fig)


def make_overlap_figure(per_recording):
    rows = per_recording.loc[
        per_recording['overlap_observed_95'].notna()
        & per_recording['overlap_null_median_95'].notna()
    ]
    fig, ax = plt.subplots(figsize=(5.2, 3.8))
    for _, row in rows.iterrows():
        ax.plot([0, 1], [row['overlap_observed_95'], row['overlap_null_median_95']], color='0.75', linewidth=0.7)
    ax.scatter(np.zeros(len(rows)), rows['overlap_observed_95'], color='#19a819', s=12, label='observed')
    ax.scatter(np.ones(len(rows)), rows['overlap_null_median_95'], color='0.35', s=12, label='shuffle median')
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['observed', 'shuffle median'])
    ax.set_xlabel('condition')
    ax.set_ylabel('mean Jaccard overlap')
    ax.spines[['top', 'right']].set_visible(False)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIGURE_ROOT / 'single_trial_overlap_95.png', dpi=200, facecolor='white')
    plt.close(fig)


#%% log
def write_log(per_recording, per_acr, all_observed, all_null, zero_n):
    lines = ['', '## 11 November 2026: corrected aligned single-trial response consistency', '', 'The validated CPU circular-shuffle test was applied to all valid trials behind Jingyu\'s saved classification, with her per-ACR NaN rule. Thresholds are the 90th, 95th and 99th percentiles of the pooled shuffled single-trial amplitudes for each ACR; real and shuffled responses also required single-trial effect size above 0.05.', '', f'Recordings: {per_recording.rec_id.nunique()} recordings, {per_recording.animal.nunique()} animals. Tables: `outputs/da_hotspots/aligned_acr_consistency_per_acr.parquet` and `outputs/da_hotspots/aligned_acr_consistency_per_recording.parquet`.']
    if zero_n:
        lines.append(f'ACRs with no kept trials were excluded from probability summaries: {zero_n} ACRs.')
    else:
        lines.append('Every valid ACR had at least one kept trial.')
    lines.append('')
    lines.append('single-trial responses, single-trial responding fractions:')
    lines.append('Fractions use the valid ACRs with a kept window on that trial; trials with no kept ACRs are omitted. Their run-onset frames, denominators and fractions are retained as aligned arrays in the per-recording table.')
    for percentile in THRESHOLDS:
        value = np.nanmedian(per_recording[f'single_trial_fraction_median_{percentile}'])
        pooled = np.concatenate(per_recording[f'single_trial_fractions_{percentile}'].to_numpy())
        lines.append(f'- {percentile}th threshold: pooled median across {len(pooled)} trial rows = {100 * np.nanmedian(pooled):.2f}%; median of per-recording trial medians = {100 * value:.2f}%; compare with the earlier 6 to 8% definition and the approximately 21% trial-averaged classification.')
    lines.append('')
    lines.append('response probability, response probability against the shuffle null:')
    lines.append(f'- The pooled 95th-threshold distribution contains {len(all_observed)} ACRs with kept trials; the null band uses the {N_SHUFFLE} aligned shuffle columns.')
    lines.append('')
    lines.append('response counts, empirical per-ACR shuffle p < 0.05 tested against 5%:')
    for percentile in THRESHOLDS:
        stat, p_value, n = safe_wilcoxon_constant(per_recording[f'excess_pct_{percentile}'], 5.0)
        lines.append(f'- {percentile}th threshold: median percentage={np.nanmedian(per_recording[f"excess_pct_{percentile}"]):.2f}%; W={stat:g}, p={p_value:.6g}, n={n} recordings, n={per_recording.animal.nunique()} animals.')
    lines.append('')
    lines.append('Up comparison, saved Up ACRs versus remaining valid ACRs, descriptive paired comparisons:')
    for percentile in THRESHOLDS:
        stat, p_value, n = safe_wilcoxon_pair(per_recording[f'up_response_probability_median_{percentile}'], per_recording[f'other_response_probability_median_{percentile}'])
        lines.append(f'- {percentile}th threshold, median response probability: W={stat:g}, p={p_value:.6g}, n={n} recordings, n={per_recording.animal.nunique()} animals.')
        stat, p_value, n = safe_wilcoxon_pair(per_recording[f'up_single_trial_fraction_mean_{percentile}'], per_recording[f'other_single_trial_fraction_mean_{percentile}'])
        lines.append(f'- {percentile}th threshold, mean single-trial responding fraction: W={stat:g}, p={p_value:.6g}, n={n} recordings, n={per_recording.animal.nunique()} animals.')
    lines.append('')
    lines.append('overlap, Jaccard overlap of responding ACR sets across trials:')
    lines.append('Each pair uses only ACRs kept on both trials, and non-empty unions are selected separately for the observed matrix and each shuffle; float64 accumulation avoids integer overflow.')
    for percentile in THRESHOLDS:
        stat, p_value, n = safe_wilcoxon_pair(per_recording[f'overlap_observed_{percentile}'], per_recording[f'overlap_null_median_{percentile}'])
        lines.append(f'- {percentile}th threshold: observed median={np.nanmedian(per_recording[f"overlap_observed_{percentile}"]):.4f}, null median={np.nanmedian(per_recording[f"overlap_null_median_{percentile}"]):.4f}; W={stat:g}, p={p_value:.6g}, n={n} recordings, n={per_recording.animal.nunique()} animals.')
    lines.extend(['', 'Figures: `figures/dlight_hotspots/single_trial_response_probability_95.png` and `figures/dlight_hotspots/single_trial_overlap_95.png`.', '', 'The saved behaviour profiles lack `reward_times_aligned`; no reconstruction was used because the validation trial set is defined by reward, stop flags and run-onset validity.'])
    (OUTPUT_ROOT / 'PROGRESS.md').open('a', encoding='utf-8').write('\n'.join(lines) + '\n')


#%% run
def main():
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    FIGURE_ROOT.mkdir(parents=True, exist_ok=True)
    per_acr_rows = []
    per_recording_rows = []
    all_observed = []
    all_null = []
    zero_n = 0
    for recname in list(recnames):
        behaviour, profile, traces = load_recording(recname)
        frames = event_frames(behaviour)
        result = quantify_grid_response(traces, frames, n_shuffle=N_SHUFFLE, seed=CHANCE_SEED, return_single_trial=True)
        keep = result['keep']
        n_trials = keep.sum(axis=1).astype(int)
        zero_n += int(np.sum(n_trials == 0))
        up = profile['Up'].to_numpy(dtype=bool)
        rec_acr = {'rec_id': recname, 'animal': recname.split('-')[0], 'roi_y': profile['roi_id'].map(lambda x: int(x[0])).to_numpy(), 'roi_x': profile['roi_id'].map(lambda x: int(x[1])).to_numpy(), 'Up': up, 'n': n_trials}
        rec_summary = {'rec_id': recname, 'animal': recname.split('-')[0], 'n_trials': len(frames), 'n_valid_acr': len(profile), 'n_zero_trial_acr': int(np.sum(n_trials == 0))}
        trial_has_acr = np.any(keep, axis=0)
        rec_summary['single_trial_frames'] = frames[trial_has_acr].tolist()
        rec_summary['single_trial_n_acr'] = keep.sum(axis=0)[trial_has_acr].tolist()
        for percentile in THRESHOLDS:
            threshold, real_response, shuffle_response = response_matrices(result, percentile)
            observed_k = real_response.sum(axis=1).astype(int)
            null_k = shuffle_response.sum(axis=1).astype(int)
            probability = np.divide(observed_k, n_trials, out=np.full(len(n_trials), np.nan), where=n_trials > 0)
            empirical_p = np.divide((null_k >= observed_k[:, None]).sum(axis=1), N_SHUFFLE, out=np.full(len(n_trials), np.nan), where=n_trials > 0)
            rec_acr[f'n_{percentile}'] = n_trials
            rec_acr[f'k_{percentile}'] = observed_k
            rec_acr[f'response_probability_{percentile}'] = probability
            rec_acr[f'empirical_p_{percentile}'] = empirical_p
            rec_acr[f'threshold_{percentile}'] = threshold
            rec_acr[f'null_k_median_{percentile}'] = np.nanmedian(null_k, axis=1)
            rec_acr[f'null_k_low_{percentile}'] = np.nanpercentile(null_k, 2.5, axis=1)
            rec_acr[f'null_k_high_{percentile}'] = np.nanpercentile(null_k, 97.5, axis=1)
            for label in ('median', 'low', 'high'):
                rec_acr[f'null_response_probability_{label}_{percentile}'] = np.divide(
                    rec_acr[f'null_k_{label}_{percentile}'], n_trials,
                    out=np.full(len(n_trials), np.nan), where=n_trials > 0)
            fractions = group_trial_fraction(real_response, keep, None)
            rec_summary[f'single_trial_fraction_median_{percentile}'] = float(np.nanmedian(fractions)) if len(fractions) else np.nan
            rec_summary[f'single_trial_fraction_n_{percentile}'] = len(fractions)
            rec_summary[f'single_trial_fractions_{percentile}'] = fractions.tolist()
            rec_summary[f'excess_pct_{percentile}'] = float(100 * np.nanmean(empirical_p < 0.05)) if np.any(np.isfinite(empirical_p)) else np.nan
            for group, label in ((True, 'up'), (False, 'other')):
                group_probability = probability[up == group]
                rec_summary[f'{label}_response_probability_median_{percentile}'] = float(np.nanmedian(group_probability)) if np.any(np.isfinite(group_probability)) else np.nan
                group_fractions = group_trial_fraction(real_response, keep, up == group)
                rec_summary[f'{label}_single_trial_fraction_mean_{percentile}'] = float(np.nanmean(group_fractions)) if len(group_fractions) else np.nan
            if percentile == 95:
                all_observed.extend(probability[np.isfinite(probability)].tolist())
                all_null.append(np.divide(null_k, n_trials[:, None], out=np.full_like(null_k, np.nan, dtype=float), where=n_trials[:, None] > 0))
            observed_overlap, null_overlap, pair_count = jaccard_summary(real_response, shuffle_response, keep)
            rec_summary[f'overlap_observed_{percentile}'] = observed_overlap
            rec_summary[f'overlap_null_median_{percentile}'] = float(np.nanmedian(null_overlap)) if np.any(np.isfinite(null_overlap)) else np.nan
            rec_summary[f'overlap_null_low_{percentile}'] = float(np.nanpercentile(null_overlap, 2.5)) if np.any(np.isfinite(null_overlap)) else np.nan
            rec_summary[f'overlap_null_high_{percentile}'] = float(np.nanpercentile(null_overlap, 97.5)) if np.any(np.isfinite(null_overlap)) else np.nan
            rec_summary[f'overlap_n_pairs_{percentile}'] = pair_count
        per_acr_rows.extend(pd.DataFrame(rec_acr).to_dict('records'))
        per_recording_rows.append(rec_summary)
        print(f'{recname}: {len(frames)} trials, {len(profile)} ACRs')
    per_acr = pd.DataFrame(per_acr_rows)
    per_recording = pd.DataFrame(per_recording_rows)
    per_acr.to_parquet(OUTPUT_ROOT / 'aligned_acr_consistency_per_acr.parquet', index=False)
    per_recording.to_parquet(OUTPUT_ROOT / 'aligned_acr_consistency_per_recording.parquet', index=False)
    all_observed = np.asarray(all_observed, dtype=float)
    all_null = np.concatenate(all_null, axis=0)
    make_probability_figure(all_observed, all_null)
    make_overlap_figure(per_recording)
    write_log(per_recording, per_acr, all_observed, all_null, zero_n)


if __name__ == '__main__':
    main()
