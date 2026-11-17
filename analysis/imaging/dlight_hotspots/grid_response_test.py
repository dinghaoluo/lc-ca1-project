# -*- coding: utf-8 -*-
'''
Created on Tue Nov 10 22:13:57 2026

cpu response and shuffle test for grid-aligned dLight traces

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


import numpy as np


#%% trace preparation and alignment
def postprocess_dff(dff, trace_filter):
    from scipy.ndimage import gaussian_filter1d
    dff = np.asarray(dff, dtype=np.float32)
    threshold = np.nanmean(dff) + 5 * np.nanstd(dff)
    filtered = np.vstack([
        trace_filter(row, fix_thresh=threshold)
        for row in dff.reshape(-1, dff.shape[-1])
    ]).reshape(dff.shape)
    return gaussian_filter1d(filtered, sigma=1, axis=-1, mode='reflect')


def align_trials(traces, event_frames, bef=2, aft=4, fs=30):
    traces = np.asarray(traces, dtype=np.float64)
    if traces.ndim != 2:
        raise ValueError('traces must have shape (n_acr, n_frames)')
    length = int((bef + aft) * fs)
    aligned = np.full((traces.shape[0], len(event_frames), length), np.nan)
    for trial, event in enumerate(np.asarray(event_frames, dtype=int)):
        if event == 0:
            continue
        start = event - int(bef * fs)
        stop = event + int(aft * fs)
        current = traces[:, start:stop]
        if current.shape[1] == length:
            aligned[:, trial, :] = current
    return aligned


def keep_aligned_trials(aligned):
    keep = np.all(np.isfinite(aligned), axis=2)
    cleaned = aligned.copy()
    cleaned[~keep] = np.nan
    return cleaned, keep


def pooled_std(baseline, response):
    n_baseline = baseline.shape[-1]
    n_response = response.shape[-1]
    std_baseline = np.nanstd(baseline.reshape(baseline.shape[0], -1), axis=1, ddof=1)
    std_response = np.nanstd(response.reshape(response.shape[0], -1), axis=1, ddof=1)
    return np.sqrt(
        ((n_baseline - 1) * std_baseline ** 2
         + (n_response - 1) * std_response ** 2)
        / (n_baseline + n_response - 2)
    )


def trial_amplitudes(aligned):
    baseline = aligned[:, :, 30:60]
    response = aligned[:, :, 60:90]
    baseline_mean = np.nanmean(baseline, axis=2)
    response_mean = np.nanmean(response, axis=2)
    return response_mean - baseline_mean, baseline, response


#%% shuffle response test
def shuffle_amplitudes(aligned, n_shuffle=1000, seed=42,
                       return_single_trial=False):
    n_acr, n_trial, length = aligned.shape
    rng = np.random.RandomState(seed)
    shuffled_mean = np.empty((n_acr, n_shuffle), dtype=np.float64)
    shuffled_single = (
        np.empty((n_acr, n_trial, n_shuffle), dtype=np.float64)
        if return_single_trial else None
    )
    finite_trial = np.isfinite(aligned[:, :, 0])
    safe = np.nan_to_num(aligned, nan=0.0)
    doubled = np.concatenate((safe, safe), axis=2)
    prefix = np.concatenate(
        (np.zeros((n_acr, n_trial, 1), dtype=np.float64),
         np.cumsum(doubled, axis=2)), axis=2)
    batch_size = 64
    for start in range(0, n_shuffle, batch_size):
        stop = min(start + batch_size, n_shuffle)
        batch = stop - start
        shifts = rng.randint(1, length, size=(batch, n_trial))
        baseline_start = (30 - shifts) % length
        response_start = (60 - shifts) % length
        prefix_batch = prefix[:, None, :, :]
        baseline_sum = (
            np.take_along_axis(
                prefix_batch,
                (baseline_start + 30)[None, :, :, None], axis=3)[:, :, :, 0]
            - np.take_along_axis(
                prefix_batch,
                baseline_start[None, :, :, None], axis=3)[:, :, :, 0]
        )
        response_sum = (
            np.take_along_axis(
                prefix_batch,
                (response_start + 30)[None, :, :, None], axis=3)[:, :, :, 0]
            - np.take_along_axis(
                prefix_batch,
                response_start[None, :, :, None], axis=3)[:, :, :, 0]
        )
        single = response_sum / 30 - baseline_sum / 30
        single = np.where(finite_trial[:, None, :], single, np.nan)
        shuffled_mean[:, start:stop] = np.nanmean(single, axis=2)
        if return_single_trial:
            shuffled_single[:, :, start:stop] = np.transpose(single, (0, 2, 1))
    return shuffled_mean, shuffled_single

def quantify_grid_response(traces, event_frames, n_shuffle=1000, seed=42,
                           return_single_trial=False):
    aligned = align_trials(traces, event_frames)
    cleaned, keep = keep_aligned_trials(aligned)
    amplitude_trial, baseline, response = trial_amplitudes(cleaned)
    amplitude = np.nanmean(amplitude_trial, axis=1)
    std = pooled_std(baseline, response)
    effect_size = amplitude / std
    mean_profile = np.nanmean(cleaned, axis=1)
    shuffled_mean, shuffled_single = shuffle_amplitudes(
        cleaned, n_shuffle=n_shuffle, seed=seed,
        return_single_trial=return_single_trial)
    return {
        'aligned': cleaned,
        'keep': keep,
        'trial_amplitude': amplitude_trial,
        'amplitude': amplitude,
        'pooled_std': std,
        'effect_size': effect_size,
        'mean_profile': mean_profile,
        'shuffle_amplitude': shuffled_mean,
        'shuffle_single_trial_amplitude': shuffled_single,
    }


def classify_up(amplitude, effect_size, shuffle_amplitude, edge):
    threshold = np.nanpercentile(shuffle_amplitude, 95, axis=1)
    up = (
        ~np.asarray(edge, dtype=bool)
        & (amplitude > threshold)
        & (effect_size > 0.05)
    )
    return up, threshold


def explicit_shift_amplitude(aligned_trial, shift):
    shifted = np.roll(aligned_trial, shift=int(shift), axis=-1)
    return np.nanmean(shifted[..., 60:90], axis=-1) - np.nanmean(shifted[..., 30:60], axis=-1)
