# -*- coding: utf-8 -*-
'''
Created on Fri Apr 10 2026

speed matching and profile comparisons used across the LC, HPC and axon first-lick analyses

@author: Dinghao Luo
'''

#%% imports
import numpy as np
from scipy.stats import sem, ttest_ind, ttest_rel, wilcoxon, ranksums

from console_formatting import (
    print_binwise_header,
    print_binwise_row,
)
from common_functions import smooth_convolve


#%% speed functions
def compute_binned_speed_matrix(trial_indices, speed_times, n_bins=7, bin_size=500):
    '''
    build a trial-by-bin mean-speed matrix for first-lick matching analyses.
    '''
    means = []
    valid = []
    total_len = n_bins * bin_size

    for trial in trial_indices:
        speeds = [pt[1] for pt in speed_times[trial]]
        if len(speeds) < total_len:
            continue
        speed_array = np.asarray(speeds[:total_len], dtype=float)
        means.append(speed_array.reshape(n_bins, bin_size).mean(axis=1))
        valid.append(trial)

    if not means:
        return np.empty((0, n_bins)), []
    return np.vstack(means), valid

def compute_session_mean_speed(trial_indices, speed_times, n=4000):
    '''
    average the leading n ms of speed across a set of trials.
    '''
    arrs = []
    for trial in trial_indices:
        speeds = [pt[1] for pt in speed_times[trial]]
        if len(speeds) >= n:
            arrs.append(np.asarray(speeds[:n], dtype=float))
    if not arrs:
        return None
    return np.nanmean(np.vstack(arrs), axis=0)


#%% profile extraction
def extract_run_onset_profiles(trains, trials, run_onset_bin, samp_freq, bef, aft,
                               rate_window=None, smoothing_sigma=None):
    '''
    extract peri-run-onset profiles and, optionally, scalar mean activity in a target window.
    '''
    profiles = []
    rates = []

    for trial in trials:
        curr_train = trains[trial]
        if smoothing_sigma is not None:
            curr_train = smooth_convolve(curr_train, smoothing_sigma)
        profile = curr_train[run_onset_bin - bef * samp_freq : run_onset_bin + aft * samp_freq]
        profiles.append(profile)
        if rate_window is not None:
            rates.append(np.mean(curr_train[rate_window[0]:rate_window[1]]))

    if rate_window is None:
        return profiles
    return profiles, rates

#%% binwise statistics
def annotate_binwise_test_suite(ax, p_ind, p_rs, p_rel, p_wil,
                                start=-0.5, bin_size=.5, fontsize=2.5, label_fontsize=2.5):
    '''
    annotate four parallel binwise test series above a first-lick profile plot.
    '''
    ymin, ymax = ax.get_ylim()
    yr = ymax - ymin if ymax > ymin else 1.0

    base = ymax + 0.04 * yr
    dy = 0.030 * yr

    label_x = -0.6
    ax.text(label_x, base + dy * 3, 'ind', transform=ax.transData, ha='right', va='bottom', fontsize=label_fontsize)
    ax.text(label_x, base + dy * 2, 'rs', transform=ax.transData, ha='right', va='bottom', fontsize=label_fontsize)
    ax.text(label_x, base + dy * 1, 'rel', transform=ax.transData, ha='right', va='bottom', fontsize=label_fontsize)
    ax.text(label_x, base + dy * 0, 'wil', transform=ax.transData, ha='right', va='bottom', fontsize=label_fontsize)

    for i in range(len(p_ind)):
        lo = start + i * bin_size
        hi = lo + bin_size
        mid = (lo + hi) / 2

        ax.text(mid, base + dy * 3, f'{p_ind[i]:.2g}', ha='center', va='bottom', fontsize=fontsize)
        ax.text(mid, base + dy * 2, f'{p_rs[i]:.2g}', ha='center', va='bottom', fontsize=fontsize)
        ax.text(mid, base + dy * 1, f'{p_rel[i]:.2g}', ha='center', va='bottom', fontsize=fontsize)
        ax.text(mid, base + dy * 0, f'{p_wil[i]:.2g}', ha='center', va='bottom', fontsize=fontsize)

def compute_binwise_test_suite(early_profiles, late_profiles, samp_freq=1250, bef=1,
                               start=-0.5, end=3.5, bin_size=.5, verbose=False,
                               label='profiles'):
    '''
    compute four complementary binwise tests for early-vs-late first-lick profiles.
    '''
    bins = np.arange(start, end, bin_size)
    n_bins = len(bins)

    p_ind = []
    p_rs = []
    p_rel = []
    p_wil = []

    if verbose:
        print_binwise_header(label, 'early', 'late', ['ind', 'rs', 'rel', 'wil'])

    for b in range(n_bins):
        lo, hi = bins[b], bins[b] + bin_size
        lo_idx = int((lo + bef) * samp_freq)
        hi_idx = int((hi + bef) * samp_freq)

        e_means = np.nanmean(early_profiles[:, lo_idx:hi_idx], axis=1)
        l_means = np.nanmean(late_profiles[:, lo_idx:hi_idx], axis=1)

        mask = np.isfinite(e_means) & np.isfinite(l_means)
        e = e_means[mask]
        l = l_means[mask]

        e_mean = np.mean(e) if e.size else np.nan
        l_mean = np.mean(l) if l.size else np.nan
        e_sem = sem(e) if e.size > 1 else np.nan
        l_sem = sem(l) if l.size > 1 else np.nan

        if e.size:
            _, p_i = ttest_ind(e, l, nan_policy='omit')
            _, p_r = ranksums(e, l, nan_policy='omit')
            _, p_p = ttest_rel(e, l, nan_policy='omit')
            if np.all(e == l):
                p_w = np.nan
            else:
                _, p_w = wilcoxon(e, l, nan_policy='omit')
        else:
            p_i = p_r = p_p = p_w = np.nan

        p_ind.append(p_i)
        p_rs.append(p_r)
        p_rel.append(p_p)
        p_wil.append(p_w)

        if verbose:
            print_binwise_row(
                f'{lo:g} to {hi:g} s',
                e_mean,
                e_sem,
                l_mean,
                l_sem,
                [p_i, p_r, p_p, p_w],
            )

    return (np.array(p_ind), np.array(p_rs), np.array(p_rel), np.array(p_wil))
