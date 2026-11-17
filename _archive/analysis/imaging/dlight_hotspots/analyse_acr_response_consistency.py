# -*- coding: utf-8 -*-
'''
Created on Tue Nov 10 21:09:53 2026

measure response consistency across valid dLight trials

@author: Dinghao Luo
'''

#%% imports
import sys
from pathlib import Path

repo_root = Path(__file__).resolve().parents[4]
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
from scipy.ndimage import gaussian_filter1d
from scipy.stats import binom, wilcoxon


#%% paths and parameters
repo_root = Path(__file__).resolve().parents[4]
JINGYU_ROOT = pp.JINGYU_CODE_ROOT
RAW_ROOT = pp.DBH_DLIGHT_RAW_ROOT
OUTPUT_ROOT = pp.DLIGHT_HOTSPOTS_OUTPUT_ROOT
FIGURE_ROOT = pp.DLIGHT_HOTSPOTS_FIGURE_ROOT

if str(JINGYU_ROOT) not in sys.path:
    sys.path.insert(0, str(JINGYU_ROOT))

recnames = rec_list.recDLightDbh
from common.utils_basic import trace_filter  # noqa: E402


SAMP_FREQ = 30
BASELINE_FRAMES = (-30, 0)
RESPONSE_FRAMES = (0, 30)
NULL_EARLY_FRAMES = (-90, -60)
NULL_LATE_FRAMES = (-60, -30)
NULL_PERCENTILES = (90, 95, 99)
N_CHANCE = 1000
CHANCE_RATE = 0.05
RNG_SEED = 0

PER_TRIAL_PATH = OUTPUT_ROOT / 'acr_response_consistency_per_trial.parquet'
PER_ACR_PATH = OUTPUT_ROOT / 'acr_response_consistency_per_acr.parquet'
PER_RECORDING_PATH = OUTPUT_ROOT / 'acr_response_consistency_per_recording.parquet'
RESPONSE_PROBABILITY_PATH = FIGURE_ROOT / 'acr_response_consistency_measure_1.png'
UP_COMPARISON_PATH = FIGURE_ROOT / 'acr_response_consistency_measure_3.png'


#%% loading and trial windows
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


def postprocess_dff(dff):
    dff = np.asarray(dff, dtype=np.float32)
    threshold = np.nanmean(dff) + 5 * np.nanstd(dff)
    filtered = np.vstack([
        trace_filter(row, fix_thresh=threshold)
        for row in dff.reshape(-1, dff.shape[-1])
    ]).reshape(dff.shape)
    return gaussian_filter1d(filtered, sigma=1, axis=-1, mode='reflect')


def valid_trial_mask(behaviour, n_frames):
    run_onset = np.asarray(behaviour['run_onset_frames'], dtype=float)
    reward = np.asarray(behaviour['reward_times'], dtype=float)
    non_stop = np.asarray(behaviour['non_stop_trials'])
    non_fullstop = np.asarray(behaviour['non_fullstop_trials'])
    n_trials = min(len(run_onset), len(reward), len(non_stop), len(non_fullstop))
    run_onset = run_onset[:n_trials]
    valid = (
        np.isfinite(run_onset)
        & (run_onset >= 0)
        & np.isfinite(reward[:n_trials])
        & (non_stop[:n_trials] == 0)
        & (non_fullstop[:n_trials] == 0)
    )
    valid &= run_onset - 90 >= 0
    valid &= run_onset + 120 <= n_frames
    return valid


def window_mean(trace, start, stop):
    window = trace[:, start:stop]
    finite = np.isfinite(window)
    count = finite.sum(axis=1)
    summed = np.nansum(window, axis=1)
    mean = np.full(trace.shape[0], np.nan, dtype=float)
    np.divide(summed, count, out=mean, where=count > 0)
    return mean


def load_profile(recname):
    path = (
        RAW_ROOT / 'processed_dataframe_grid_free_dilation'
        / f'{recname}_profile_combined_dilation=0_pre(-1, 0)_post(0, 1)'
        '_ES=0.05_shuff95_0621.parquet'
    )
    profile = pd.read_parquet(path)
    valid = (
        profile['dlight_valid'].astype(bool)
        & profile['red_valid'].astype(bool)
        & ~profile['edge'].astype(bool)
    )
    profile = profile.loc[valid].copy()
    coords = np.asarray(profile['roi_id'].tolist(), dtype=int)
    if coords.ndim != 2 or coords.shape[1] != 2:
        raise ValueError(f'{recname}: invalid valid-ACR coordinates')
    if len(np.unique(coords, axis=0)) != len(coords):
        raise ValueError(f'{recname}: duplicate valid-ACR coordinates')
    profile['Up'] = profile['Up'].astype(bool).to_numpy()
    return profile, coords


def load_recording(recname):
    dff_path = (
        RAW_ROOT / 'regression_res_grid_free_dilation' / recname
        / 'single_trial_regression' / 'dilation_k=0'
        / 'dff_corrected_dlight.npy'
    )
    behaviour_path = RAW_ROOT / 'behaviour_profile' / f'{recname}.pkl'
    dff = postprocess_dff(np.load(dff_path, mmap_mode='r'))
    behaviour = load_behaviour(behaviour_path)
    profile, coords = load_profile(recname)
    return dff, behaviour, profile, coords


def response_arrays(dff, coords, trial_indices):
    trace = dff[coords[:, 0], coords[:, 1], :]
    delta = np.full((len(trial_indices), len(coords)), np.nan)
    null_delta = np.full_like(delta, np.nan)
    for row, trial in enumerate(trial_indices):
        onset = int(trial)
        baseline = window_mean(trace, onset + BASELINE_FRAMES[0], onset + BASELINE_FRAMES[1])
        response = window_mean(trace, onset + RESPONSE_FRAMES[0], onset + RESPONSE_FRAMES[1])
        null_early = window_mean(trace, onset + NULL_EARLY_FRAMES[0], onset + NULL_EARLY_FRAMES[1])
        null_late = window_mean(trace, onset + NULL_LATE_FRAMES[0], onset + NULL_LATE_FRAMES[1])
        delta[row] = response - baseline
        null_delta[row] = null_late - null_early
    return delta, null_delta


def null_thresholds(null_delta):
    thresholds = {}
    for percentile in NULL_PERCENTILES:
        values = np.full(null_delta.shape[1], np.nan)
        for roi in range(null_delta.shape[1]):
            sample = null_delta[np.isfinite(null_delta[:, roi]), roi]
            if len(sample):
                values[roi] = np.nanpercentile(sample, percentile)
        thresholds[percentile] = values
    return thresholds


def build_trial_rows(recname, animal, profile, trial_indices, delta, null_delta, thresholds):
    rows = []
    coords = profile['roi_id'].tolist()
    up = profile['Up'].to_numpy(dtype=bool)
    for trial_row, trial_index in enumerate(trial_indices):
        for roi, (y, x) in enumerate(coords):
            d_value = delta[trial_row, roi]
            n_value = null_delta[trial_row, roi]
            response_valid = bool(np.isfinite(d_value) and np.isfinite(n_value))
            item = {
                'rec_id': recname,
                'animal': animal,
                'trial_index': int(trial_index),
                'roi_y': int(y),
                'roi_x': int(x),
                'Up': bool(up[roi]),
                'delta': d_value,
                'null_delta': n_value,
                'response_valid': response_valid,
            }
            for percentile in NULL_PERCENTILES:
                threshold = thresholds[percentile][roi]
                item[f'null_threshold_{percentile}'] = threshold
                item[f'responding_{percentile}'] = bool(
                    response_valid and np.isfinite(threshold) and d_value > threshold
                )
            rows.append(item)
    return rows


#%% summaries and statistics
def median_or_nan(values):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    return float(np.nanmedian(values)) if len(values) else np.nan


def per_acr_summary(table):
    group = table.groupby(['rec_id', 'animal', 'roi_y', 'roi_x', 'Up'])
    rows = group.agg(
        n=('response_valid', 'sum'),
        k_90=('responding_90', 'sum'),
        k_95=('responding_95', 'sum'),
        k_99=('responding_99', 'sum'),
    ).reset_index()
    for percentile in NULL_PERCENTILES:
        rows[f'response_probability_{percentile}'] = np.where(
            rows['n'] > 0, rows[f'k_{percentile}'] / rows['n'], np.nan)
        rows[f'excess_p_{percentile}'] = np.where(
            rows['n'] > 0,
            binom.sf(rows[f'k_{percentile}'] - 1, rows['n'], CHANCE_RATE),
            np.nan,
        )
    return rows


def single_trial_fractions(table, percentile, group=None):
    values = table.loc[table['response_valid']]
    if group is not None:
        values = values.loc[values['Up'].eq(group)]
    if values.empty:
        return np.array([], dtype=float)
    column = f'responding_{percentile}'
    fractions = values.groupby('trial_index')[column].mean().to_numpy(dtype=float)
    return fractions[np.isfinite(fractions)]


def recording_summary(recname, animal, table, per_acr, n_trials):
    row = {
        'rec_id': recname,
        'animal': animal,
        'n_trials': int(n_trials),
        'n_valid_acr': int(len(per_acr)),
        'n_up': int(per_acr['Up'].sum()),
    }
    for percentile in NULL_PERCENTILES:
        fractions = single_trial_fractions(table, percentile)
        row[f'single_trial_fraction_median_{percentile}'] = median_or_nan(fractions)
        row[f'single_trial_fraction_n_{percentile}'] = int(len(fractions))
        probability = per_acr[f'response_probability_{percentile}']
        excess = per_acr[f'excess_p_{percentile}']
        row[f'n_response_acr_{percentile}'] = int(probability.notna().sum())
        row[f'n_zero_response_acr_{percentile}'] = int(probability.isna().sum())
        row[f'excess_pct_{percentile}'] = (
            float(100 * (excess.dropna() < 0.05).mean())
            if excess.notna().any() else np.nan
        )
        for group, label in ((True, 'up'), (False, 'other')):
            group_probability = probability.loc[per_acr['Up'].eq(group)]
            row[f'{label}_response_probability_median_{percentile}'] = median_or_nan(group_probability)
            group_fractions = single_trial_fractions(table, percentile, group)
            row[f'{label}_single_trial_fraction_mean_{percentile}'] = (
                float(np.nanmean(group_fractions)) if len(group_fractions) else np.nan
            )
    return row


def wilcoxon_pair(first, second):
    first = np.asarray(first, dtype=float)
    second = np.asarray(second, dtype=float)
    keep = np.isfinite(first) & np.isfinite(second)
    first, second = first[keep], second[keep]
    if len(first) == 0:
        return 0.0, np.nan, 0
    difference = first - second
    if np.allclose(difference, 0):
        return 0.0, 1.0, len(first)
    result = wilcoxon(first, second, alternative='two-sided')
    return float(result.statistic), float(result.pvalue), len(first)


def wilcoxon_constant(values, constant):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if len(values) == 0:
        return 0.0, np.nan, 0
    difference = values - constant
    if np.allclose(difference, 0):
        return 0.0, 1.0, len(values)
    result = wilcoxon(difference, alternative='two-sided')
    return float(result.statistic), float(result.pvalue), len(values)


#%% figures
def ecdf(values, x_grid):
    values = np.sort(np.asarray(values, dtype=float))
    return np.searchsorted(values, x_grid, side='right') / len(values)


def plot_measure_1(per_acr):
    observed_rows = per_acr.loc[per_acr['response_probability_95'].notna()]
    values = observed_rows['response_probability_95'].to_numpy(dtype=float)
    counts = observed_rows['n'].to_numpy(dtype=int)
    x_grid = np.linspace(0, 1, 201)
    observed = ecdf(values, x_grid)
    rng = np.random.default_rng(RNG_SEED)
    simulated = rng.binomial(
        counts[None, :], CHANCE_RATE, size=(N_CHANCE, len(counts))) / counts[None, :]
    simulated_ecdf = np.vstack([ecdf(row, x_grid) for row in simulated])
    median = np.nanmedian(simulated_ecdf, axis=0)
    low, high = np.nanpercentile(simulated_ecdf, [2.5, 97.5], axis=0)
    fig, ax = plt.subplots(figsize=(4.6, 3.8))
    ax.fill_between(x_grid, low, high, color='0.75', alpha=0.55, label='chance, 95% band')
    ax.plot(x_grid, median, color='0.35', linewidth=1.5, label='chance median')
    ax.plot(x_grid, observed, color='#19a819', linewidth=1.5, label='observed')
    ax.set_xlabel('response probability per ACR')
    ax.set_ylabel('cumulative fraction of ACRs')
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.spines[['top', 'right']].set_visible(False)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(RESPONSE_PROBABILITY_PATH, dpi=200, facecolor='white')
    plt.close(fig)


def plot_measure_3(per_recording):
    rows = per_recording.loc[
        per_recording['up_response_probability_median_95'].notna()
        & per_recording['other_response_probability_median_95'].notna()
    ].copy()
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.5), sharex=True)
    pairs = [
        ('up_response_probability_median_95', 'other_response_probability_median_95', 'median response probability'),
        ('up_single_trial_fraction_mean_95', 'other_single_trial_fraction_mean_95', 'mean single-trial responding fraction'),
    ]
    for ax, (up_column, other_column, ylabel) in zip(axes, pairs):
        for _, row in rows.iterrows():
            ax.plot([0, 1], [row[up_column], row[other_column]], color='0.75', linewidth=0.7)
        ax.scatter(np.zeros(len(rows)), rows[up_column], color='#19a819', s=10, label='Up ACRs')
        ax.scatter(np.ones(len(rows)), rows[other_column], color='0.35', s=10, label='other ACRs')
        ax.set_xticks([0, 1])
        ax.set_xticklabels(['Up', 'other'])
        ax.set_xlabel('ACR group')
        ax.set_ylabel(ylabel)
        ax.spines[['top', 'right']].set_visible(False)
    axes[0].legend(frameon=False)
    fig.tight_layout()
    fig.savefig(UP_COMPARISON_PATH, dpi=200, facecolor='white')
    plt.close(fig)


#%% progress log
def write_progress(per_trial, per_acr, per_recording, selector_missing, length_issues):
    lines = [
        '',
        '## 10 November 2026: dLight ACR response consistency analysis',
        '',
        'The analysis used all trials passing the permissive `valid_trial_mask` from `compare_early_late_dlight.py`, with no first-lick restriction. The saved behaviour profiles do not contain `reward_times_aligned`, so Jingyu\'s exact selector could not be reproduced; that field was not reconstructed.',
        '',
        f'Recordings: {per_recording.rec_id.nunique()} recordings, {per_recording.animal.nunique()} animals.',
        'Tables: `outputs/da_hotspots/acr_response_consistency_per_trial.parquet`, `outputs/da_hotspots/acr_response_consistency_per_acr.parquet`, and `outputs/da_hotspots/acr_response_consistency_per_recording.parquet`.',
    ]
    if length_issues:
        lines.append('Trial-array length ambiguities were recorded:')
        lines.extend(f'- `{recname}`: {detail}.' for recname, detail in length_issues)
    else:
        lines.append('No trial-array length ambiguity was found.')
    lines.append('')
    lines.append('single-trial responses, single-trial responding fractions:')
    for percentile in NULL_PERCENTILES:
        column = f'single_trial_fraction_median_{percentile}'
        pooled = np.nanmedian(per_recording[column].to_numpy(dtype=float))
        lines.append(f'- {percentile}th threshold: median of per-recording trial medians = {100 * pooled:.2f}%; pooled trial rows = {int(per_recording[f"single_trial_fraction_n_{percentile}"].sum())}.')
    lines.append('')
    lines.append('response probability, response probability against the 5% independent-response null:')
    lines.append(f'- The pooled 95th-threshold distribution contains {int(per_acr["response_probability_95"].notna().sum())} ACRs with at least one valid response; chance distributions used {N_CHANCE} binomial draws with seed {RNG_SEED}, retaining each ACR\'s observed valid-trial count.')
    lines.append('')
    lines.append('response counts, ACRs with one-sided binomial p < 0.05, tested against 5%:')
    for percentile in NULL_PERCENTILES:
        stat, p_value, n = wilcoxon_constant(per_recording[f'excess_pct_{percentile}'], 5.0)
        median = np.nanmedian(per_recording[f'excess_pct_{percentile}'])
        lines.append(f'- {percentile}th threshold: median percentage = {median:.2f}%; W={stat:g}, p={p_value:.6g}, n={n} recordings, n={per_recording.animal.nunique()} animals.')
    lines.append('')
    lines.append('Up comparison, Up ACRs versus the remaining valid ACRs:')
    for percentile in NULL_PERCENTILES:
        first = per_recording[f'up_response_probability_median_{percentile}']
        second = per_recording[f'other_response_probability_median_{percentile}']
        stat, p_value, n = wilcoxon_pair(first, second)
        lines.append(f'- {percentile}th threshold, median response probability: Up median={np.nanmedian(first):.4f}, other median={np.nanmedian(second):.4f}; W={stat:g}, p={p_value:.6g}, n={n} recordings, n={per_recording.animal.nunique()} animals.')
        first = per_recording[f'up_single_trial_fraction_mean_{percentile}']
        second = per_recording[f'other_single_trial_fraction_mean_{percentile}']
        stat, p_value, n = wilcoxon_pair(first, second)
        lines.append(f'- {percentile}th threshold, mean single-trial responding fraction: Up mean={np.nanmean(first):.4f}, other mean={np.nanmean(second):.4f}; W={stat:g}, p={p_value:.6g}, n={n} recordings, n={per_recording.animal.nunique()} animals.')
    lines.extend([
        '',
        f'Figures: `figures/dlight_hotspots/{RESPONSE_PROBABILITY_PATH.name}` and `figures/dlight_hotspots/{UP_COMPARISON_PATH.name}`.',
        '',
        f'The missing selector field affected {len(selector_missing)} recordings; no missing input was resolved.',
    ])
    (OUTPUT_ROOT / 'PROGRESS.md').open('a', encoding='utf-8').write('\n'.join(lines) + '\n')


#%% run
def main():
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    FIGURE_ROOT.mkdir(parents=True, exist_ok=True)
    trial_rows = []
    recording_inputs = []
    selector_missing = []
    length_issues = []
    recordings = list(recnames)
    if len(recordings) != 57:
        raise ValueError(f'Expected 57 recordings, found {len(recordings)}')
    for recname in recordings:
        dff, behaviour, profile, coords = load_recording(recname)
        lengths = {key: len(np.asarray(behaviour[key])) for key in (
            'run_onset_frames', 'reward_times', 'non_stop_trials', 'non_fullstop_trials')}
        if len(set(lengths.values())) != 1:
            length_issues.append((recname, str(lengths)))
        mask = valid_trial_mask(behaviour, dff.shape[-1])
        trial_indices = np.flatnonzero(mask).tolist()
        if 'reward_times_aligned' not in behaviour:
            selector_missing.append(recname)
        delta, null_delta = response_arrays(dff, coords, trial_indices)
        thresholds = null_thresholds(null_delta)
        trial_rows.extend(build_trial_rows(
            recname, recname.split('-')[0], profile, trial_indices,
            delta, null_delta, thresholds))
        recording_inputs.append((recname, len(trial_indices), len(coords), int(profile['Up'].sum())))
        print(f'{recname}: {len(trial_indices)} valid trials, {len(coords)} valid ACRs')
    per_trial = pd.DataFrame(trial_rows)
    per_acr = per_acr_summary(per_trial)
    recording_rows = []
    for recname, n_trials, n_acr, n_up in recording_inputs:
        rec_trial = per_trial.loc[per_trial['rec_id'].eq(recname)]
        rec_acr = per_acr.loc[per_acr['rec_id'].eq(recname)]
        recording_rows.append(recording_summary(
            recname, recname.split('-')[0], rec_trial, rec_acr, n_trials))
    per_recording = pd.DataFrame(recording_rows)
    per_trial.to_parquet(PER_TRIAL_PATH, index=False)
    per_acr.to_parquet(PER_ACR_PATH, index=False)
    per_recording.to_parquet(PER_RECORDING_PATH, index=False)
    plot_measure_1(per_acr)
    plot_measure_3(per_recording)
    write_progress(per_trial, per_acr, per_recording, selector_missing, length_issues)


if __name__ == '__main__':
    main()
