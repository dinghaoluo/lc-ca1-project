# -*- coding: utf-8 -*-
'''
Created on Fri Nov  6 18:28:54 2026

compare early and late first-lick dLight hotspot responses

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

import argparse
import importlib

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from common_functions import mpl_formatting
mpl_formatting()
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d
from scipy.stats import wilcoxon


#%% paths and parameters
repo_root = Path(__file__).resolve().parents[3]
JINGYU_ROOT = pp.JINGYU_CODE_ROOT
RAW_ROOT = pp.DBH_DLIGHT_RAW_ROOT
OUTPUT_ROOT = pp.DLIGHT_HOTSPOTS_OUTPUT_ROOT
FIGURE_ROOT = pp.DLIGHT_HOTSPOTS_FIGURE_ROOT

if str(JINGYU_ROOT) not in sys.path:
    sys.path.insert(0, str(JINGYU_ROOT))

recnames = rec_list.recDLightDbh
from common.utils_basic import trace_filter  # noqa: E402


SAMP_FREQ = 30
GRID_SIZE = 32
BASELINE_FRAMES = (-30, 0)
RESPONSE_FRAMES = (0, 30)
NULL_EARLY_FRAMES = (-90, -60)
NULL_LATE_FRAMES = (-60, -30)
FIRST_LICK_MIN_MS = 1000
EARLY_MAX_MS = 2500
LATE_MAX_MS = 3500
N_SPEED_BINS = 7
SPEED_BIN_MS = 500
MIN_MATCHED_TRIALS = 10
NULL_PERCENTILES = (90, 95, 99)
TRIAL_COLUMNS = [
    'rec_id', 'animal', 'speed_mode', 'condition', 'trial_index',
    'roi_y', 'roi_x', 'delta', 'null_delta', 'response_valid',
    'null_threshold_90', 'null_threshold_95', 'null_threshold_99',
    'responding_90', 'responding_95', 'responding_99',
]


#%% loading
def postprocess_dff(dff):
    dff = np.asarray(dff, dtype=np.float32)
    threshold = np.nanmean(dff) + 5 * np.nanstd(dff)
    safe = np.vstack([
        trace_filter(row, fix_thresh=threshold)
        for row in dff.reshape(-1, dff.shape[-1])
    ]).reshape(dff.shape)
    return gaussian_filter1d(safe, sigma=1, axis=-1, mode='reflect')


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


def load_recording_inputs(recname):
    dff_path = (
        RAW_ROOT / 'regression_res_grid_free_dilation' / recname
        / 'single_trial_regression' / 'dilation_k=0'
        / 'dff_corrected_dlight.npy'
    )
    behaviour_path = RAW_ROOT / 'behaviour_profile' / f'{recname}.pkl'
    profile_path = (
        RAW_ROOT / 'processed_dataframe_grid_free_dilation'
        / f'{recname}_profile_combined_dilation=0_pre(-1, 0)_post(0, 1)'
        '_ES=0.05_shuff95_0621.parquet'
    )

    dff = postprocess_dff(np.load(dff_path, mmap_mode='r'))
    behaviour = load_behaviour(behaviour_path)
    profile = pd.read_parquet(profile_path)
    profile = profile.loc[
        profile['dlight_valid'].astype(bool)
        & profile['red_valid'].astype(bool)
        & ~profile['edge'].astype(bool)
    ].copy()
    coords = np.asarray(profile['roi_id'].tolist(), dtype=int)
    if coords.ndim != 2 or coords.shape[1] != 2:
        raise ValueError(f'{recname}: ROI coordinates are not an n-by-2 array')
    if np.any(coords < 0) or np.any(coords >= GRID_SIZE):
        raise ValueError(f'{recname}: ROI coordinates fall outside the 32-by-32 grid')
    if len(np.unique(coords, axis=0)) != len(coords):
        raise ValueError(f'{recname}: duplicate ROI coordinates in the saved profile')
    return dff, behaviour, profile, coords


#%% trial selection and speed matching
def first_lick_times(behaviour):
    first_lick = np.full(len(behaviour['lick_times_aligned']), np.nan)
    for trial, licks in enumerate(behaviour['lick_times_aligned']):
        values = np.asarray(licks, dtype=float).reshape(-1)
        values = values[np.isfinite(values) & (values >= FIRST_LICK_MIN_MS)]
        if values.size:
            first_lick[trial] = values[0]
    return first_lick


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


def speed_bins(behaviour, trial, run_onset_ms):
    trace = np.asarray(behaviour['speed_times_aligned'][trial], dtype=float)
    if trace.ndim != 2 or trace.shape[1] < 2:
        return None
    relative_ms = trace[:, 0] - run_onset_ms
    values = np.full(N_SPEED_BINS, np.nan)
    for bin_idx in range(N_SPEED_BINS):
        start = bin_idx * SPEED_BIN_MS
        stop = start + SPEED_BIN_MS
        in_bin = (relative_ms >= start) & (relative_ms < stop)
        speeds = trace[in_bin, 1]
        if len(speeds) < SPEED_BIN_MS or not np.all(np.isfinite(speeds)):
            return None
        values[bin_idx] = np.mean(speeds)
    return values


def speed_match_trials(behaviour, early_trials, late_trials):
    run_onsets = np.asarray(behaviour['run_onsets'], dtype=float)
    all_trials = list(early_trials) + list(late_trials)
    binned = {
        trial: speed_bins(behaviour, trial, run_onsets[trial])
        for trial in all_trials
        if trial < len(run_onsets) and np.isfinite(run_onsets[trial])
    }
    early_available = [trial for trial in early_trials if binned.get(trial) is not None]
    late_available = [trial for trial in late_trials if binned.get(trial) is not None]
    if not early_available or not late_available:
        return [], [], len(early_available), len(late_available)

    early_speed = np.vstack([binned[trial] for trial in early_available])
    late_speed = np.vstack([binned[trial] for trial in late_available])
    early_mean = np.mean(early_speed, axis=0)
    early_std = np.std(early_speed, axis=0)
    late_mean = np.mean(late_speed, axis=0)
    late_std = np.std(late_speed, axis=0)
    early_keep = np.all(
        (early_speed >= late_mean - late_std)
        & (early_speed <= late_mean + late_std),
        axis=1,
    )
    late_keep = np.all(
        (late_speed >= early_mean - early_std)
        & (late_speed <= early_mean + early_std),
        axis=1,
    )
    early_matched = [trial for trial, keep in zip(early_available, early_keep) if keep]
    late_matched = [trial for trial, keep in zip(late_available, late_keep) if keep]
    return early_matched, late_matched, len(early_available), len(late_available)


#%% dLight response calculation
def window_mean(trace, start, stop):
    window = trace[:, start:stop]
    finite = np.isfinite(window)
    count = finite.sum(axis=1)
    summed = np.nansum(window, axis=1)
    mean = np.full(trace.shape[0], np.nan, dtype=float)
    np.divide(summed, count, out=mean, where=count > 0)
    return mean


def trial_response_arrays(dff, coords, trial_indices, n_frames):
    trace = dff[coords[:, 0], coords[:, 1], :]
    delta = np.full((len(trial_indices), len(coords)), np.nan)
    null_delta = np.full_like(delta, np.nan)
    valid_trials = []
    for row, trial in enumerate(trial_indices):
        onset = int(trial['run_onset_frame'])
        if onset + RESPONSE_FRAMES[1] > n_frames:
            continue
        baseline = window_mean(trace, onset + BASELINE_FRAMES[0], onset + BASELINE_FRAMES[1])
        response = window_mean(trace, onset + RESPONSE_FRAMES[0], onset + RESPONSE_FRAMES[1])
        null_early = window_mean(trace, onset + NULL_EARLY_FRAMES[0], onset + NULL_EARLY_FRAMES[1])
        null_late = window_mean(trace, onset + NULL_LATE_FRAMES[0], onset + NULL_LATE_FRAMES[1])
        delta[row] = response - baseline
        null_delta[row] = null_late - null_early
        valid_trials.append(row)
    return delta, null_delta, valid_trials


def recording_trial_info(behaviour, n_frames):
    valid = valid_trial_mask(behaviour, n_frames)
    first_lick = first_lick_times(behaviour)
    n_trials = len(valid)
    run_onsets = np.asarray(behaviour['run_onset_frames'], dtype=float)[:n_trials]
    first_lick = first_lick[:n_trials]
    early = np.flatnonzero(valid & (first_lick < EARLY_MAX_MS)).tolist()
    late = np.flatnonzero(
        valid & (first_lick >= EARLY_MAX_MS) & (first_lick <= LATE_MAX_MS)
    ).tolist()
    trial_indices = [
        {'trial_index': trial, 'run_onset_frame': int(run_onsets[trial])}
        for trial in np.flatnonzero(valid)
    ]
    return trial_indices, early, late


def null_thresholds(null_delta):
    thresholds = {}
    for percentile in NULL_PERCENTILES:
        thresholds[percentile] = np.full(null_delta.shape[1], np.nan)
        for roi in range(null_delta.shape[1]):
            values = null_delta[np.isfinite(null_delta[:, roi]), roi]
            if values.size:
                thresholds[percentile][roi] = np.nanpercentile(values, percentile)
    return thresholds


def make_trial_table(recname, animal, coords, all_trials, delta, null_delta,
                     thresholds, selected_trials, condition, speed_mode):
    positions = {item['trial_index']: row for row, item in enumerate(all_trials)}
    rows = []
    for trial in selected_trials:
        row = positions[trial]
        for roi, (y, x) in enumerate(coords):
            d_value = delta[row, roi]
            n_value = null_delta[row, roi]
            item = {
                'rec_id': recname,
                'animal': animal,
                'speed_mode': speed_mode,
                'condition': condition,
                'trial_index': trial,
                'roi_y': int(y),
                'roi_x': int(x),
                'delta': d_value,
                'null_delta': n_value,
                'response_valid': bool(np.isfinite(d_value) and np.isfinite(n_value)),
            }
            for percentile in NULL_PERCENTILES:
                threshold = thresholds[percentile][roi]
                item[f'null_threshold_{percentile}'] = threshold
                item[f'responding_{percentile}'] = bool(
                    np.isfinite(d_value)
                    and np.isfinite(n_value)
                    and np.isfinite(threshold)
                    and d_value > threshold
                )
            rows.append(item)
    return rows


#%% recording summaries and statistics
def response_summary(trial_table, condition, percentile):
    subset = trial_table.loc[trial_table['condition'].eq(condition)]
    trial_rows = []
    for trial, group in subset.groupby('trial_index', sort=True):
        valid = group['response_valid'].to_numpy(dtype=bool)
        n_valid = int(valid.sum())
        if n_valid == 0:
            continue
        responding = group.loc[valid, f'responding_{percentile}'].to_numpy(dtype=bool)
        trial_rows.append(100 * responding.mean())
    values = np.asarray(trial_rows, dtype=float)
    if values.size == 0:
        return np.nan, np.nan, np.nan, 0
    median = float(np.nanmedian(values))
    sd = float(np.nanstd(values, ddof=1)) if values.size > 1 else np.nan
    cv = float(sd / np.nanmean(values)) if values.size > 1 and np.nanmean(values) else np.nan
    return median, sd, cv, int(values.size)


def recording_summary(recname, animal, speed_mode, included, early_trials, late_trials,
                      trial_table, n_acr):
    row = {
        'rec_id': recname,
        'animal': animal,
        'speed_mode': speed_mode,
        'included': included,
        'n_early_trials': len(early_trials),
        'n_late_trials': len(late_trials),
        'n_acr': n_acr,
    }
    condition_medians = {}
    for condition in ('early', 'late'):
        subset = trial_table.loc[trial_table['condition'].eq(condition)]
        per_acr = subset.groupby(['roi_y', 'roi_x'])['delta'].mean()
        per_acr = pd.to_numeric(per_acr, errors='coerce').dropna()
        median = float(np.nanmedian(per_acr)) if len(per_acr) else np.nan
        condition_medians[condition] = median
        row[f'{condition}_acr_delta_median'] = median
        row[f'{condition}_acr_delta_n'] = int(len(per_acr))
        for percentile in NULL_PERCENTILES:
            median_pct, sd_pct, cv_pct, n_pct_trials = response_summary(
                trial_table, condition, percentile
            )
            row[f'{condition}_responding_pct_median_{percentile}'] = median_pct
            row[f'{condition}_responding_pct_sd_{percentile}'] = sd_pct
            row[f'{condition}_responding_pct_cv_{percentile}'] = cv_pct
            row[f'{condition}_responding_pct_n_trials_{percentile}'] = n_pct_trials
    row['late_minus_early_acr_delta_median'] = (
        condition_medians['late'] - condition_medians['early']
        if np.isfinite(condition_medians['late']) and np.isfinite(condition_medians['early'])
        else np.nan
    )
    for percentile in NULL_PERCENTILES:
        for metric in ('median', 'sd', 'cv'):
            early_value = row[f'early_responding_pct_{metric}_{percentile}']
            late_value = row[f'late_responding_pct_{metric}_{percentile}']
            row[f'late_minus_early_responding_pct_{metric}_{percentile}'] = (
                late_value - early_value
                if np.isfinite(early_value) and np.isfinite(late_value)
                else np.nan
            )
    return row


def paired_test(table, metric, speed_mode):
    subset = table.loc[table['speed_mode'].eq(speed_mode) & table['included']].copy()
    early = subset[f'early_{metric}'].to_numpy(dtype=float)
    late = subset[f'late_{metric}'].to_numpy(dtype=float)
    keep = np.isfinite(early) & np.isfinite(late)
    early, late = early[keep], late[keep]
    if len(early) < 2:
        return {'speed_mode': speed_mode, 'metric': metric, 'n_recordings': len(early),
                'n_animals': subset.loc[keep, 'animal'].nunique(), 'statistic': np.nan,
                'p_value': np.nan}
    try:
        result = wilcoxon(late, early, alternative='two-sided')
        statistic, p_value = float(result.statistic), float(result.pvalue)
    except ValueError:
        statistic, p_value = 0.0, 1.0
    return {
        'speed_mode': speed_mode,
        'metric': metric,
        'n_recordings': len(early),
        'n_animals': subset.loc[keep, 'animal'].nunique(),
        'statistic': statistic,
        'p_value': p_value,
    }


#%% figures
def paired_plot(table, metric, ylabel, path):
    fig, axs = plt.subplots(1, 2, figsize=(5.4, 2.4), sharey=True)
    for ax, speed_mode, label in zip(
            axs,
            ('without_speed', 'with_speed'),
            ('without speed matching', 'with speed matching')):
        subset = table.loc[table['speed_mode'].eq(speed_mode) & table['included']]
        for _, row in subset.iterrows():
            ax.plot(
                [0, 1],
                [row[f'early_{metric}'], row[f'late_{metric}']],
                color='0.65', linewidth=0.6,
            )
        ax.scatter(
            np.zeros(len(subset)), subset[f'early_{metric}'],
            color='0.25', s=8, label='early', zorder=3,
        )
        ax.scatter(
            np.ones(len(subset)), subset[f'late_{metric}'],
            color='tab:green', s=8, label='late', zorder=3,
        )
        ax.set_xticks([0, 1], ['early', 'late'])
        ax.set_xlabel(label)
        ax.spines[['top', 'right']].set_visible(False)
    axs[0].set_ylabel(ylabel)
    axs[1].legend(frameon=False, fontsize=7)
    fig.tight_layout()
    fig.savefig(path, dpi=300)
    plt.close(fig)


def example_plot(trial_table, recname, path):
    subset = trial_table.loc[trial_table['rec_id'].eq(recname)]
    fig, axs = plt.subplots(1, 2, figsize=(5.2, 2.4))
    distributions = [
        subset.loc[subset['condition'].eq(condition), 'delta'].dropna().to_numpy()
        for condition in ('early', 'late')
    ]
    axs[0].violinplot(distributions, positions=[0, 1], showmedians=True)
    axs[0].set_xticks([0, 1], ['early', 'late'])
    axs[0].set_ylabel('dLight Δ')
    axs[0].spines[['top', 'right']].set_visible(False)
    for condition, colour in (('early', '0.25'), ('late', 'tab:green')):
        values = []
        for _, group in subset.loc[subset['condition'].eq(condition)].groupby('trial_index'):
            valid = group['response_valid'].to_numpy(dtype=bool)
            if valid.any():
                values.append(100 * group.loc[valid, 'responding_95'].mean())
        axs[1].plot(values, color=colour, linewidth=0.8, label=condition)
    axs[1].set_xlabel('trial')
    axs[1].set_ylabel('responding ACRs (%)')
    axs[1].legend(frameon=False, fontsize=7)
    axs[1].spines[['top', 'right']].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=300)
    plt.close(fig)


#%% run
def analyse_recording(recname):
    dff, behaviour, profile, coords = load_recording_inputs(recname)
    n_frames = dff.shape[-1]
    animal = recname.split('-')[0]
    all_trials, early_trials, late_trials = recording_trial_info(behaviour, n_frames)
    delta, null_delta, _ = trial_response_arrays(dff, coords, all_trials, n_frames)
    thresholds = null_thresholds(null_delta)

    tables = []
    summaries = []
    mode_trials = {'without_speed': (early_trials, late_trials, True)}
    early_matched, late_matched, early_available, late_available = speed_match_trials(
        behaviour, early_trials, late_trials
    )
    mode_trials['with_speed'] = (
        early_matched,
        late_matched,
        len(early_matched) >= MIN_MATCHED_TRIALS
        and len(late_matched) >= MIN_MATCHED_TRIALS,
    )
    for speed_mode, (early, late, included) in mode_trials.items():
        rows = []
        rows.extend(make_trial_table(
            recname, animal, coords, all_trials, delta, null_delta, thresholds,
            early, 'early', speed_mode,
        ))
        rows.extend(make_trial_table(
            recname, animal, coords, all_trials, delta, null_delta, thresholds,
            late, 'late', speed_mode,
        ))
        table = pd.DataFrame(rows, columns=TRIAL_COLUMNS)
        tables.append(table)
        summaries.append(recording_summary(
            recname, animal, speed_mode,
            included and len(early) > 0 and len(late) > 0,
            early, late, table, len(coords),
        ))
    return (
        pd.concat(tables, ignore_index=True),
        pd.DataFrame(summaries),
        {
            'n_profile_rois': len(profile),
            'n_early_trials': len(early_trials),
            'n_late_trials': len(late_trials),
            'n_early_speed_available': early_available,
            'n_late_speed_available': late_available,
        },
    )


def main(argv=None):
    parser = argparse.ArgumentParser(description='compare early and late first-lick dLight responses')
    parser.add_argument('--recording', help='optional recording substring')
    args = parser.parse_args(argv)

    recordings = list(recnames)
    if args.recording:
        recordings = [recname for recname in recordings if args.recording in recname]
    if not recordings:
        raise ValueError('no recordings matched the requested filter')

    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    FIGURE_ROOT.mkdir(parents=True, exist_ok=True)
    trial_tables = []
    summary_tables = []
    for recname in recordings:
        print(f'loading {recname}', flush=True)
        trial_table, summary_table, rec_metadata = analyse_recording(recname)
        trial_tables.append(trial_table)
        summary_tables.append(summary_table)
        print(
            f'  early={rec_metadata["n_early_trials"]}, '
            f'late={rec_metadata["n_late_trials"]}, '
            f'valid ACRs={rec_metadata["n_profile_rois"]}',
            flush=True,
        )

    per_trial = pd.concat(trial_tables, ignore_index=True)
    per_recording = pd.concat(summary_tables, ignore_index=True)
    per_trial_path = OUTPUT_ROOT / 'early_late_per_acr_trial.parquet'
    per_recording_path = OUTPUT_ROOT / 'early_late_per_recording.parquet'
    per_trial.to_parquet(per_trial_path, index=False)
    per_recording.to_parquet(per_recording_path, index=False)

    stats = []
    for speed_mode in ('without_speed', 'with_speed'):
        stats.append(paired_test(per_recording, 'acr_delta_median', speed_mode))
        for percentile in NULL_PERCENTILES:
            for metric in ('responding_pct_median', 'responding_pct_sd', 'responding_pct_cv'):
                stats.append(paired_test(
                    per_recording,
                    f'{metric}_{percentile}',
                    speed_mode,
                ))

    paired_plot(
        per_recording, 'acr_delta_median', 'dLight Δ',
        FIGURE_ROOT / 'early_late_amplitude_paired.png',
    )
    paired_plot(
        per_recording, 'responding_pct_median_95', 'responding ACRs (%)',
        FIGURE_ROOT / 'early_late_responding_paired_95.png',
    )
    example_candidates = per_recording.loc[
        per_recording['speed_mode'].eq('without_speed') & per_recording['included'],
        'rec_id',
    ]
    if len(example_candidates):
        example_plot(
            per_trial, example_candidates.iloc[0],
            FIGURE_ROOT / 'early_late_example_recording.png',
        )

    stats = pd.DataFrame(stats)
    print(per_recording.to_string(index=False))
    print(stats.to_string(index=False))


if __name__ == '__main__':
    main()
