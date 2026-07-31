# -*- coding: utf-8 -*-
'''
Created on 6 Jul 2026

support functions for CA1 place-cell, time-cell, and sequence analyses

@author: Dinghao Luo
'''

from __future__ import annotations

#%% imports
from dataclasses import dataclass

import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.stats import spearmanr

from common_functions import get_trialtype_indices_from_stim_conds
import first_lick_analysis_functions as flaf


#%% constants
SAMP_FREQ       = 1250
RUN_ONSET_BIN   = 3750
TRACK_LENGTH_CM = 200
DIST_BIN_CM     = 0.1


#%% dataclasses
@dataclass
class TrialGroups:
    baseline: list[int]
    stim: list[int]
    ctrl: list[int]
    early: list[int]
    late: list[int]
    matched_early: list[int]
    matched_late: list[int]
    first_lick_s: np.ndarray


def safe_pearson(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 3:
        return np.nan
    if np.nanstd(x[mask]) == 0 or np.nanstd(y[mask]) == 0:
        return np.nan
    return float(np.corrcoef(x[mask], y[mask])[0, 1])


def safe_spearman(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    mask = np.isfinite(x) & np.isfinite(y)
    if mask.sum() < 3:
        return np.nan
    if np.nanstd(x[mask]) == 0 or np.nanstd(y[mask]) == 0:
        return np.nan
    r, _ = spearmanr(x[mask], y[mask])
    return float(r)


def smooth_nan_1d(values, sigma):
    values = np.asarray(values, dtype=float)
    if sigma is None or sigma <= 0:
        return values
    if np.all(~np.isfinite(values)):
        return np.full_like(values, np.nan)

    finite = np.isfinite(values).astype(float)
    filled = np.nan_to_num(values, nan=0.0)
    smooth_values = gaussian_filter1d(filled, sigma=sigma, mode='nearest')
    smooth_weights = gaussian_filter1d(finite, sigma=sigma, mode='nearest')
    with np.errstate(divide='ignore', invalid='ignore'):
        out = smooth_values / smooth_weights
    out[smooth_weights == 0] = np.nan
    return out


def downsample_1d(values, factor):
    values = np.asarray(values, dtype=float)
    if factor <= 1:
        return values
    n_bins = values.size // factor
    if n_bins == 0:
        return np.array([], dtype=float)
    trimmed = values[:n_bins * factor]
    return np.nanmean(trimmed.reshape(n_bins, factor), axis=1)


#%% trial grouping
def valid_behaviour_trials(beh, n_trials):
    bad_trials = beh.get('bad_trials', [])
    valid = []
    for trial in range(n_trials):
        beh_trial = trial + 1
        if beh_trial >= len(bad_trials):
            valid.append(False)
            continue
        bad = bad_trials[beh_trial]
        valid.append(bool(bad is False or bad == 0))
    return np.asarray(valid, dtype=bool)


def first_lick_times_s(beh, n_trials, min_lick_ms=500):
    first = np.full(n_trials, np.nan, dtype=float)
    lick_times = beh.get('lick_times_aligned', [])

    for trial in range(n_trials):
        beh_trial = trial + 1
        if beh_trial >= len(lick_times):
            continue
        licks = lick_times[beh_trial]
        if licks is None or isinstance(licks, float):
            continue
        licks = np.asarray(licks, dtype=float)
        licks = licks[np.isfinite(licks) & (licks > min_lick_ms)]
        if licks.size:
            first[trial] = licks[0] / 1000

    return first


def stim_ctrl_trial_indices(beh, n_trials, exclude_bad=True):
    trial_statements = beh.get('trial_statements', [])
    stim_conds = [
        trial[15]
        for trial in trial_statements[1:n_trials + 1]
        if len(trial) > 15
    ]
    baseline, stim, ctrl = get_trialtype_indices_from_stim_conds(stim_conds)
    baseline = [
        int(trial)
        for trial in baseline
        if trial is not None and 0 <= int(trial) < n_trials
    ]
    stim = [
        int(trial)
        for trial in stim
        if trial is not None and 0 <= int(trial) < n_trials
    ]
    ctrl = [
        int(trial)
        for trial in ctrl
        if trial is not None and 0 <= int(trial) < n_trials
    ]

    if not exclude_bad:
        return baseline, stim, ctrl

    valid = valid_behaviour_trials(beh, n_trials)
    baseline = [trial for trial in baseline if valid[trial]]
    stim = [trial for trial in stim if valid[trial]]
    ctrl = [trial for trial in ctrl if valid[trial]]
    return baseline, stim, ctrl


def early_late_trial_indices(
        beh,
        n_trials,
        early_cutoff_s=2.5,
        late_window_s=(2.5, 3.5),
        exclude_stim=True,
        exclude_previous_stim=True
        ):
    first = first_lick_times_s(beh, n_trials)
    valid = valid_behaviour_trials(beh, n_trials)
    _, stim, _ = stim_ctrl_trial_indices(beh, n_trials, exclude_bad=False)
    stim_set = set(stim)

    early = []
    late = []
    for trial, lick_s in enumerate(first):
        if not valid[trial] or not np.isfinite(lick_s):
            continue
        if exclude_stim and trial in stim_set:
            continue
        if exclude_previous_stim and trial - 1 in stim_set:
            continue
        if lick_s < early_cutoff_s:
            early.append(trial)
        elif late_window_s[0] < lick_s < late_window_s[1]:
            late.append(trial)

    return early, late, first


def speed_matched_early_late_trials(
        beh,
        early,
        late,
        n_bins=7,
        bin_size=500,
        k=1.5,
        min_trials=10
        ):
    speed_times = beh.get('speed_times_aligned', [])[1:]
    early_bins, early_valid = flaf.compute_binned_speed_matrix(
        early,
        speed_times,
        n_bins=n_bins,
        bin_size=bin_size,
    )
    late_bins, late_valid = flaf.compute_binned_speed_matrix(
        late,
        speed_times,
        n_bins=n_bins,
        bin_size=bin_size,
    )
    if len(early_bins) == 0 or len(late_bins) == 0:
        return [], []

    early_mu = early_bins.mean(axis=0)
    early_sd = early_bins.std(axis=0, ddof=0)
    late_mu = late_bins.mean(axis=0)
    late_sd = late_bins.std(axis=0, ddof=0)

    late_mask = np.all(
        (late_bins >= early_mu - k * early_sd)
        & (late_bins <= early_mu + k * early_sd),
        axis=1,
    )
    early_mask = np.all(
        (early_bins >= late_mu - k * late_sd)
        & (early_bins <= late_mu + k * late_sd),
        axis=1,
    )

    matched_early = [early_valid[i] for i in np.where(early_mask)[0]]
    matched_late = [late_valid[i] for i in np.where(late_mask)[0]]
    if len(matched_early) < min_trials or len(matched_late) < min_trials:
        return [], []
    return matched_early, matched_late


def build_trial_groups(beh, n_trials, min_trials=10, has_stim=True):
    if has_stim:
        baseline, stim, ctrl = stim_ctrl_trial_indices(beh, n_trials)
        early, late, first = early_late_trial_indices(beh, n_trials)
    else:
        valid = valid_behaviour_trials(beh, n_trials)
        baseline = np.where(valid)[0].astype(int).tolist()
        stim = []
        ctrl = []
        early, late, first = early_late_trial_indices(
            beh,
            n_trials,
            exclude_stim=False,
            exclude_previous_stim=False,
        )
    matched_early, matched_late = speed_matched_early_late_trials(
        beh,
        early,
        late,
        min_trials=min_trials,
    )
    return TrialGroups(
        baseline=baseline,
        stim=stim,
        ctrl=ctrl,
        early=early,
        late=late,
        matched_early=matched_early,
        matched_late=matched_late,
        first_lick_s=first,
    )


#%% maps and information
def information_from_rate_map(rate_map, occupancy=None):
    rate_map = np.asarray(rate_map, dtype=float)
    valid = np.isfinite(rate_map) & (rate_map >= 0)
    if occupancy is not None:
        occupancy = np.asarray(occupancy, dtype=float)
        valid &= np.isfinite(occupancy) & (occupancy > 0)

    if valid.sum() < 3:
        return np.nan

    rates = rate_map[valid]
    if occupancy is None:
        p = np.ones(rates.size, dtype=float) / rates.size
    else:
        occ = occupancy[valid]
        p = occ / np.nansum(occ)

    mean_rate = np.nansum(p * rates)
    if not np.isfinite(mean_rate) or mean_rate <= 0:
        return np.nan

    with np.errstate(divide='ignore', invalid='ignore'):
        ratio = rates / mean_rate
        info = np.nansum(p * ratio * np.log2(ratio))
    return float(info)


def make_spatial_rate_maps(
        trains_dist_cell,
        occupancy,
        trials,
        smooth_sigma_bins=20,
        min_occupancy_s=0.02
        ):
    trial_maps = []
    sum_spikes = None
    sum_occ = None

    for trial in trials:
        if trial >= len(trains_dist_cell) or trial >= len(occupancy):
            continue
        spikes = trains_dist_cell[trial]
        occ = occupancy[trial]
        if spikes is None or occ is None:
            continue

        spikes = np.asarray(spikes, dtype=float)
        occ = np.asarray(occ, dtype=float)
        n = min(spikes.size, occ.size)
        if n == 0:
            continue
        spikes = spikes[:n]
        occ = occ[:n]

        with np.errstate(divide='ignore', invalid='ignore'):
            trial_map = spikes / occ
        trial_map[occ < min_occupancy_s] = np.nan
        trial_map = smooth_nan_1d(trial_map, smooth_sigma_bins)
        trial_maps.append(trial_map)

        if sum_spikes is None:
            sum_spikes = np.zeros(n, dtype=float)
            sum_occ = np.zeros(n, dtype=float)
        n_acc = min(n, sum_spikes.size)
        sum_spikes[:n_acc] += np.nan_to_num(spikes[:n_acc], nan=0.0)
        sum_occ[:n_acc] += np.nan_to_num(occ[:n_acc], nan=0.0)

    if not trial_maps or sum_spikes is None:
        return np.empty((0, 0)), np.array([]), np.array([])

    min_len = min(len(trial_map) for trial_map in trial_maps)
    trial_maps = np.vstack([trial_map[:min_len] for trial_map in trial_maps])
    sum_spikes = sum_spikes[:min_len]
    sum_occ = sum_occ[:min_len]
    with np.errstate(divide='ignore', invalid='ignore'):
        mean_map = sum_spikes / sum_occ
    mean_map[sum_occ < min_occupancy_s] = np.nan
    mean_map = smooth_nan_1d(mean_map, smooth_sigma_bins)
    return trial_maps, mean_map, sum_occ


def make_temporal_rate_maps(
        trains_cell,
        trials,
        start_s=-0.5,
        end_s=3.5,
        bin_size_s=0.05,
        run_onset_bin=RUN_ONSET_BIN,
        samp_freq=SAMP_FREQ
        ):
    start = int(run_onset_bin + start_s * samp_freq)
    stop = int(run_onset_bin + end_s * samp_freq)
    bin_samples = max(1, int(bin_size_s * samp_freq))
    trial_maps = []

    for trial in trials:
        if trial >= len(trains_cell):
            continue
        trace = np.asarray(trains_cell[trial], dtype=float)
        if trace.size < stop:
            continue
        window = trace[start:stop]
        trial_maps.append(downsample_1d(window, bin_samples))

    if not trial_maps:
        return np.empty((0, 0)), np.array([])

    min_len = min(len(trial_map) for trial_map in trial_maps)
    trial_maps = np.vstack([trial_map[:min_len] for trial_map in trial_maps])
    mean_map = np.nanmean(trial_maps, axis=0)
    return trial_maps, mean_map


def split_half_reliability(trial_maps):
    trial_maps = np.asarray(trial_maps, dtype=float)
    if trial_maps.ndim != 2 or trial_maps.shape[0] < 4:
        return np.nan
    odd = np.nanmean(trial_maps[::2], axis=0)
    even = np.nanmean(trial_maps[1::2], axis=0)
    return safe_pearson(odd, even)


def circular_shift_information_null(
        trial_maps,
        occupancy=None,
        n_shuf=200,
        rng=None
        ):
    trial_maps = np.asarray(trial_maps, dtype=float)
    if rng is None:
        rng = np.random.default_rng()
    if trial_maps.ndim != 2 or trial_maps.shape[0] == 0 or trial_maps.shape[1] < 4:
        return np.full(n_shuf, np.nan)

    n_bins = trial_maps.shape[1]
    null = np.full(n_shuf, np.nan)
    for shuf in range(n_shuf):
        shifted = np.empty_like(trial_maps)
        for trial in range(trial_maps.shape[0]):
            shift = int(rng.integers(1, n_bins))
            shifted[trial] = np.roll(trial_maps[trial], shift)
        null[shuf] = information_from_rate_map(
            np.nanmean(shifted, axis=0),
            occupancy=occupancy,
        )
    return null


def field_width(rate_map, bin_size, threshold_frac=0.2):
    rate_map = np.asarray(rate_map, dtype=float)
    if rate_map.size == 0 or np.all(~np.isfinite(rate_map)):
        return np.nan
    peak_idx = int(np.nanargmax(rate_map))
    peak = rate_map[peak_idx]
    if not np.isfinite(peak) or peak <= 0:
        return np.nan

    threshold = peak * threshold_frac
    left = peak_idx
    right = peak_idx
    while left > 0 and np.isfinite(rate_map[left - 1]) and rate_map[left - 1] >= threshold:
        left -= 1
    while right < rate_map.size - 1 and np.isfinite(rate_map[right + 1]) and rate_map[right + 1] >= threshold:
        right += 1
    return float((right - left + 1) * bin_size)


def field_summary(
        trial_maps,
        mean_map,
        bin_size,
        occupancy=None,
        axis_start=0.0,
        n_shuf=200,
        rng=None,
        alpha=0.05,
        min_peak_rate=1.0,
        min_split_r=0.2,
        min_width=0.0,
        max_width=np.inf
        ):
    mean_map = np.asarray(mean_map, dtype=float)
    if mean_map.size == 0 or np.all(~np.isfinite(mean_map)):
        return {
            'info': np.nan,
            'info_p': np.nan,
            'split_r': np.nan,
            'peak': np.nan,
            'peak_rate': np.nan,
            'field_width': np.nan,
            'is_field': False,
        }

    info = information_from_rate_map(mean_map, occupancy=occupancy)
    null = circular_shift_information_null(
        trial_maps,
        occupancy=occupancy,
        n_shuf=n_shuf,
        rng=rng,
    )
    if np.all(~np.isfinite(null)) or not np.isfinite(info):
        info_p = np.nan
    else:
        info_p = (np.sum(null >= info) + 1) / (np.sum(np.isfinite(null)) + 1)

    split_r = split_half_reliability(trial_maps)
    peak_idx = int(np.nanargmax(mean_map))
    peak_rate = float(np.nanmax(mean_map))
    width = field_width(mean_map, bin_size)
    peak = axis_start + peak_idx * bin_size
    is_field = (
        np.isfinite(info_p)
        and info_p <= alpha
        and np.isfinite(split_r)
        and split_r >= min_split_r
        and np.isfinite(peak_rate)
        and peak_rate >= min_peak_rate
        and np.isfinite(width)
        and min_width <= width <= max_width
    )

    return {
        'info': float(info) if np.isfinite(info) else np.nan,
        'info_p': float(info_p) if np.isfinite(info_p) else np.nan,
        'split_r': float(split_r) if np.isfinite(split_r) else np.nan,
        'peak': float(peak) if np.isfinite(peak) else np.nan,
        'peak_rate': peak_rate,
        'field_width': width,
        'is_field': bool(is_field),
    }


#%% sequence summaries
def sequence_peak_summary(reference_peaks, condition_peaks):
    reference_peaks = np.asarray(reference_peaks, dtype=float)
    condition_peaks = np.asarray(condition_peaks, dtype=float)
    mask = np.isfinite(reference_peaks) & np.isfinite(condition_peaks)
    if mask.sum() < 3:
        return {
            'n': int(mask.sum()),
            'rank_r': np.nan,
            'median_shift': np.nan,
            'median_abs_shift': np.nan,
        }

    shift = condition_peaks[mask] - reference_peaks[mask]
    return {
        'n': int(mask.sum()),
        'rank_r': safe_spearman(reference_peaks[mask], condition_peaks[mask]),
        'median_shift': float(np.nanmedian(shift)),
        'median_abs_shift': float(np.nanmedian(np.abs(shift))),
    }


def population_vector_correlation(matrix_a, matrix_b):
    matrix_a = np.asarray(matrix_a, dtype=float)
    matrix_b = np.asarray(matrix_b, dtype=float)
    if matrix_a.ndim != 2 or matrix_b.ndim != 2:
        return np.nan
    n_bins = min(matrix_a.shape[1], matrix_b.shape[1])
    if n_bins == 0:
        return np.nan
    corrs = [
        safe_pearson(matrix_a[:, idx], matrix_b[:, idx])
        for idx in range(n_bins)
    ]
    return float(np.nanmean(corrs)) if np.any(np.isfinite(corrs)) else np.nan


def order_by_peak(maps, peaks):
    maps = np.asarray(maps, dtype=float)
    peaks = np.asarray(peaks, dtype=float)
    valid = np.isfinite(peaks)
    if maps.ndim != 2 or valid.sum() == 0:
        return np.array([], dtype=int), np.empty((0, 0))
    valid_idx = np.where(valid)[0]
    order = valid_idx[np.argsort(peaks[valid_idx])]
    return order, maps[order]
