# -*- coding: utf-8 -*-
'''
Created on Mar 31 2026

LC-DA-CA1 model; generates figure 6 panels in
'figures/fig_6_lc_da_ca1_model/'

@author: Dinghao Luo
'''


#%% imports
from __future__ import annotations

import argparse
import os
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
from matplotlib.ticker import MaxNLocator
import numpy as np
from scipy.stats import ranksums, sem, ttest_ind, ttest_rel, wilcoxon

MODEL_DIR = Path(__file__).resolve().parent
REPO_ROOT = MODEL_DIR.parent
DEFAULT_OUTPUT_DIR = REPO_ROOT / 'figures' / 'fig_6_lc_da_ca1_model'

plt.rcParams.update(
    {
        'font.family': 'Arial',
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
        'font.size': 10,
        'axes.titlesize': 11,
        'axes.labelsize': 10,
        'xtick.labelsize': 9,
        'ytick.labelsize': 9,
        'legend.fontsize': 9,
    }
)

FIREBRICK = (178 / 255, 34 / 255, 34 / 255)
ROYALBLUE = (65 / 255, 105 / 255, 225 / 255)
PURPLE = (128 / 255, 0 / 255, 128 / 255)
GREY = (128 / 255, 128 / 255, 128 / 255)
DARKGREY = (0.35, 0.35, 0.35)
MIDGREY = (0.45, 0.45, 0.45)
LIGHTCORAL = (240 / 255, 128 / 255, 128 / 255)
PLUM = (221 / 255, 160 / 255, 221 / 255)
GOLDENROD = (218 / 255, 165 / 255, 32 / 255)
MAGENTA = (1.0, 0.0, 1.0)
DARKGREEN = (0.0, 100 / 255, 0.0)
RED = (1.0, 0.0, 0.0)
CLASS_COLORS = {
    'is_up': FIREBRICK,
    'is_down': PURPLE,
}

CONDITION_COLORS = {
    'baseline': DARKGREY,
    'lc': ROYALBLUE,
    'blocked': CLASS_COLORS['is_down'],
    'da_targeted': CLASS_COLORS['is_up'],
    'not_targeted': MIDGREY,
}

@dataclass
class Params:
    '''Selected model parameters.'''

    # simulation grid
    dt: float = 0.01
    t_pre: float = 1.00
    t_post: float = 6.00

    # bootstrap controls
    n_bootstrap: int = 50
    seed_start: int = 0
    lc_activation_fold: float = 2.30

    # synthetic population
    n_cells: int = 1000

    # population priors
    baseline_mean: float = 0.80
    baseline_sd: float = 0.44
    wR_mean: float = 0.88
    wR_sd: float = 1.35
    wW_mean: float = 0.97
    wW_sd: float = 0.98

    # DA targeting and state-dependent DA gain
    frac_da_targ: float = 0.30
    da_half_rate: float = 2.10
    da_rate_slope: float = 0.08
    wDA_global: float = 0.00
    da_block_scale: float = 0.40

    # cell-intrinsic recovery
    intrinsic_tau_mean: float = 0.77
    baseline_tau_coupling: float = 0.17
    intrinsic_tau_sd: float = 0.04
    intrinsic_tau_min: float = 0.05
    intrinsic_tau_max: float = 3.00

    # output nonlinearity / output limits
    softplus_beta: float = 2.10
    max_rate: float = 20.00

    # run-drive shape
    run_on_mid: float = 0.0
    run_off_mid: float = 1.60
    run_rise_scale: float = 0.09
    run_fall_scale: float = 0.57

    # reward/task-drive shape
    reward_on_mid: float = 0.18
    reward_off_mid: float = 1.31
    reward_rise_scale: float = 0.05
    reward_fall_scale: float = 0.69

    # LC -> DA conversion
    lc_baseline: float = 1.00
    lc_amp: float = 1.50
    lc_mu: float = 0.00
    lc_sigma: float = 0.20
    lc_to_da_gain: float = 5.60
    da_kernel_tau: float = 1.60
    da_ca1_delay: float = 0.34

    # classification
    pre_window: tuple[float, float] = (-1.00, -0.50)
    post_window: tuple[float, float] = (0.50, 1.50)
    up_thresh: float = 1.50
    down_thresh: float = 2.00 / 3.00

def sigmoid(x: np.ndarray | float, midpoint: float, scale: float) -> np.ndarray:
    return 1.0 / (1.0 + np.exp(-(np.asarray(x, dtype=float) - midpoint) / scale))

def softplus(x: np.ndarray, beta: float) -> np.ndarray:
    xb = beta * np.asarray(x, dtype=float)
    return (np.log1p(np.exp(-np.abs(xb))) + np.maximum(xb, 0.0)) / beta

def smooth_pulse(
    tt: np.ndarray,
    onset_mid: float,
    offset_mid: float,
    rise_scale: float,
    fall_scale: float,
) -> np.ndarray:
    pulse = sigmoid(tt, onset_mid, rise_scale) * (1.0 - sigmoid(tt, offset_mid, fall_scale))
    return pulse / np.max(pulse)

def exponential_kernel(dt: float, tau: float, t_max: float) -> np.ndarray:
    tk = np.arange(0.0, t_max + dt, dt)
    kernel = np.exp(-tk / tau) / tau
    kernel[0] = 0.0
    kernel /= np.sum(kernel) * dt
    return kernel

def make_drives(t: np.ndarray, p: Params) -> dict[str, np.ndarray]:
    '''build the selected rate-form shared drives.

    R(t) is a 0-1 run drive. W(t) is the selected rate-form task drive: the
    older suppressive pulse reparameterised as 1 - pulse, so it decays from
    high to low instead of becoming negative. D(t) is the baseline delayed DA
    waveform normalised to peak at 1; LC stimulation scales this same basis
    above 1.
    '''

    run_shape = smooth_pulse(t, p.run_on_mid, p.run_off_mid, p.run_rise_scale, p.run_fall_scale)
    reward_shape = smooth_pulse(t, p.reward_on_mid, p.reward_off_mid, p.reward_rise_scale, p.reward_fall_scale)
    R = run_shape
    W = 1.0 - reward_shape

    L = p.lc_baseline + p.lc_amp * np.exp(-0.5 * ((t - p.lc_mu) / p.lc_sigma) ** 2)
    kernel = exponential_kernel(p.dt, p.da_kernel_tau, p.t_pre + p.t_post)
    lc_excess = np.maximum(L - p.lc_baseline, 0.0)
    D_release = p.lc_to_da_gain * np.convolve(lc_excess, kernel, mode='full')[: len(t)] * p.dt
    if p.da_ca1_delay <= 0:
        D = np.asarray(D_release, dtype=float).copy()
    else:
        D = np.interp(t - p.da_ca1_delay, t, D_release, left=0.0, right=0.0)

    peak = float(np.max(D))
    D_release = D_release / peak
    D = D / peak

    return {
        'R': R,
        'W': W,
        'D': D,
        'L': L,
        'D_release': D_release,
        'run_shape': run_shape,
        'reward_shape': reward_shape,
    }

def response_strength(rates: np.ndarray, t: np.ndarray, p: Params) -> np.ndarray:
    pre_mask = (t >= p.pre_window[0]) & (t < p.pre_window[1])
    post_mask = (t >= p.post_window[0]) & (t < p.post_window[1])
    pre_mean = np.mean(rates[:, pre_mask], axis=1)
    post_mean = np.mean(rates[:, post_mask], axis=1)
    return post_mean / pre_mean

def match_paired_pre_baseline(
    left_traces: np.ndarray,
    right_traces: np.ndarray,
    t: np.ndarray,
    window: tuple[float, float],
) -> tuple[np.ndarray, np.ndarray]:
    left = np.asarray(left_traces, dtype=float).copy()
    right = np.asarray(right_traces, dtype=float).copy()
    mask = (t >= window[0]) & (t < window[1])
    left_pre = np.nanmean(left[:, mask], axis=1, keepdims=True)
    right_pre = np.nanmean(right[:, mask], axis=1, keepdims=True)
    shared_pre = 0.5 * (left_pre + right_pre)
    return left + (shared_pre - left_pre), right + (shared_pre - right_pre)

def make_population(p: Params, rng: np.random.Generator) -> dict[str, np.ndarray]:
    b = rng.normal(p.baseline_mean, p.baseline_sd, p.n_cells)
    baseline_z = (p.baseline_mean - b) / p.baseline_sd
    wR = rng.normal(p.wR_mean, p.wR_sd, p.n_cells)
    wW = rng.normal(p.wW_mean, p.wW_sd, p.n_cells)
    tau_center = p.intrinsic_tau_mean + p.baseline_tau_coupling * baseline_z
    tau = rng.normal(tau_center, p.intrinsic_tau_sd, p.n_cells)

    da_targ_strength = rng.beta(0.35, 0.35, p.n_cells).astype(float)
    threshold = np.quantile(da_targ_strength, 1.0 - p.frac_da_targ)
    da_targ = da_targ_strength >= threshold

    return {
        'b': b,
        'wR': wR,
        'wW': wW,
        'tau_intr': np.clip(tau, p.intrinsic_tau_min, p.intrinsic_tau_max),
        'da_targ': da_targ,
        'da_targ_strength': da_targ_strength,
    }

def simulate_population_condition(
    t: np.ndarray,
    p: Params,
    pop: dict[str, np.ndarray],
    drives: dict[str, np.ndarray],
    da_scale: float = 1.0,
) -> dict[str, Any]:
    x_no_da = (
        pop['b'][:, None]
        + pop['wR'][:, None] * drives['R'][None, :]
        + pop['wW'][:, None] * drives['W'][None, :]
    )

    n_cells, n_t = x_no_da.shape
    alpha = np.clip(p.dt / pop['tau_intr'], 0.0, 1.0)
    latent_state = np.zeros((n_cells, n_t), dtype=float)
    rates = np.zeros((n_cells, n_t), dtype=float)

    rate_prev = np.clip(softplus(x_no_da[:, 0], p.softplus_beta), 0.0, p.max_rate)

    for k_idx in range(n_t):
        wDA_extra = sigmoid(rate_prev, p.da_half_rate, p.da_rate_slope)
        wDA_t = p.wDA_global + pop['da_targ_strength'] * wDA_extra
        da_term = da_scale * wDA_t * drives['D'][k_idx]
        state_target = x_no_da[:, k_idx] + da_term

        if k_idx == 0:
            latent_state[:, k_idx] = state_target
        else:
            latent_state[:, k_idx] = latent_state[:, k_idx - 1] + alpha * (state_target - latent_state[:, k_idx - 1])

        rates[:, k_idx] = np.clip(softplus(latent_state[:, k_idx], p.softplus_beta), 0.0, p.max_rate)
        rate_prev = rates[:, k_idx]

    resp = response_strength(rates, t, p)
    is_up = resp >= p.up_thresh
    is_down = resp <= p.down_thresh
    return {
        'rates': rates,
        'classes': {
            'is_up': is_up,
            'is_down': is_down,
        },
        'mean_traces': {
            'is_up': np.mean(rates[is_up], axis=0),
            'is_down': np.mean(rates[is_down], axis=0),
        },
    }

def run_bootstrap_suite(p: Params) -> dict[str, Any]:
    t = np.arange(-p.t_pre, p.t_post, p.dt)
    drives = make_drives(t, p)

    trace_keys = [
        'base_up',
        'lc_up',
        'block_up',
        'base_down',
        'lc_down',
        'block_down',
        'exp1_base_up_same',
        'exp1_lc_up_same',
        'exp1_base_down_same',
        'exp1_lc_down_same',
        'exp3_base_up_same',
        'exp3_block_up_same',
        'exp3_base_down_same',
        'exp3_block_down_same',
        'da_up',
        'non_da_up',
        'da_down',
        'non_da_down',
    ]
    traces = {key: [] for key in trace_keys}
    stats = {
        'base_up_pct': [],
        'base_down_pct': [],
        'lc_up_pct': [],
        'lc_down_pct': [],
        'block_up_pct': [],
        'block_down_pct': [],
        'p_up_da_targeted': [],
        'p_up_not_targeted': [],
        'p_down_da_targeted': [],
        'p_down_not_targeted': [],
        'base_up_pre': [],
        'base_down_pre': [],
        'base_up_04': [],
        'base_down_04': [],
        'base_up_late': [],
        'base_down_late': [],
        'base_up_wW': [],
        'base_down_wW': [],
        'base_up_tau': [],
        'base_down_tau': [],
    }

    pre_mask = (t >= p.pre_window[0]) & (t < p.pre_window[1])
    win04_mask = (t >= 0.0) & (t < 4.0)
    late_mask = (t >= 3.0) & (t < 4.0)

    for seed in range(p.seed_start, p.seed_start + p.n_bootstrap):
        rng = np.random.default_rng(seed)
        pop = make_population(p, rng)

        base = simulate_population_condition(t, p, pop, drives, da_scale=1.0)
        lc = simulate_population_condition(t, p, pop, drives, da_scale=p.lc_activation_fold)
        block = simulate_population_condition(t, p, pop, drives, da_scale=p.da_block_scale)

        for condition, result in (('base', base), ('lc', lc), ('block', block)):
            traces[f'{condition}_up'].append(result['mean_traces']['is_up'])
            traces[f'{condition}_down'].append(result['mean_traces']['is_down'])

        da_mask = pop['da_targ']
        not_da_mask = ~da_mask
        base_classes = base['classes']
        lc_classes = lc['classes']
        block_classes = block['classes']

        exp1_up_same = base_classes['is_up'] & lc_classes['is_up']
        exp1_down_same = base_classes['is_down'] & lc_classes['is_down']
        exp3_up_same = base_classes['is_up'] & block_classes['is_up']
        exp3_down_same = base_classes['is_down'] & block_classes['is_down']

        traces['exp1_base_up_same'].append(np.mean(base['rates'][exp1_up_same], axis=0))
        traces['exp1_lc_up_same'].append(np.mean(lc['rates'][exp1_up_same], axis=0))
        traces['exp1_base_down_same'].append(np.mean(base['rates'][exp1_down_same], axis=0))
        traces['exp1_lc_down_same'].append(np.mean(lc['rates'][exp1_down_same], axis=0))
        traces['exp3_base_up_same'].append(np.mean(base['rates'][exp3_up_same], axis=0))
        traces['exp3_block_up_same'].append(np.mean(block['rates'][exp3_up_same], axis=0))
        traces['exp3_base_down_same'].append(np.mean(base['rates'][exp3_down_same], axis=0))
        traces['exp3_block_down_same'].append(np.mean(block['rates'][exp3_down_same], axis=0))

        da_up_mask = da_mask & base_classes['is_up']
        non_da_up_mask = not_da_mask & base_classes['is_up']
        da_down_mask = da_mask & base_classes['is_down']
        non_da_down_mask = not_da_mask & base_classes['is_down']

        traces['da_up'].append(np.mean(base['rates'][da_up_mask], axis=0))
        traces['non_da_up'].append(np.mean(base['rates'][non_da_up_mask], axis=0))
        traces['da_down'].append(np.mean(base['rates'][da_down_mask], axis=0))
        traces['non_da_down'].append(np.mean(base['rates'][non_da_down_mask], axis=0))

        base_up = base_classes['is_up']
        base_down = base_classes['is_down']
        stats['base_up_pct'].append(100.0 * np.mean(base_up))
        stats['base_down_pct'].append(100.0 * np.mean(base_down))
        stats['lc_up_pct'].append(100.0 * np.mean(lc_classes['is_up']))
        stats['lc_down_pct'].append(100.0 * np.mean(lc_classes['is_down']))
        stats['block_up_pct'].append(100.0 * np.mean(block_classes['is_up']))
        stats['block_down_pct'].append(100.0 * np.mean(block_classes['is_down']))
        stats['p_up_da_targeted'].append(100.0 * np.mean(base_classes['is_up'][da_mask]))
        stats['p_up_not_targeted'].append(100.0 * np.mean(base_classes['is_up'][not_da_mask]))
        stats['p_down_da_targeted'].append(100.0 * np.mean(base_classes['is_down'][da_mask]))
        stats['p_down_not_targeted'].append(100.0 * np.mean(base_classes['is_down'][not_da_mask]))

        for label, mask in (('base_up', base_up), ('base_down', base_down)):
            stats[f'{label}_pre'].append(float(np.mean(base['rates'][mask][:, pre_mask])))
            stats[f'{label}_04'].append(float(np.mean(base['rates'][mask][:, win04_mask])))
            stats[f'{label}_late'].append(float(np.mean(base['rates'][mask][:, late_mask])))
            stats[f'{label}_wW'].append(float(np.mean(pop['wW'][mask])))
            stats[f'{label}_tau'].append(float(np.mean(pop['tau_intr'][mask])))

    return {
        't': t,
        'params': p,
        'drives': drives,
        'traces': {key: np.asarray(value, dtype=float) for key, value in traces.items()},
        'stats': {key: np.asarray(value, dtype=float) for key, value in stats.items()},
    }

def paired_test(left: np.ndarray, right: np.ndarray) -> dict[str, float]:
    left = np.asarray(left, dtype=float)
    right = np.asarray(right, dtype=float)
    mask = np.isfinite(left) & np.isfinite(right)
    left = left[mask]
    right = right[mask]
    t_res = ttest_rel(right, left, nan_policy='omit')
    w_res = wilcoxon(right, left, zero_method='wilcox', method='approx')
    p_w = float(w_res.pvalue)
    return {
        'delta': float(np.mean(right - left)),
        'p_t': float(t_res.pvalue),
        'p_w': p_w,
        'n': float(len(left)),
    }

def independent_test(left: np.ndarray, right: np.ndarray) -> dict[str, float]:
    left = np.asarray(left, dtype=float)
    right = np.asarray(right, dtype=float)
    left = left[np.isfinite(left)]
    right = right[np.isfinite(right)]
    t_res = ttest_ind(right, left, equal_var=False, nan_policy='omit')
    rank_res = ranksums(right, left)
    p_rank = float(rank_res.pvalue)
    return {
        'delta': float(np.mean(right) - np.mean(left)),
        'p_t': float(t_res.pvalue),
        'p_w': p_rank,
        'n': float(min(len(left), len(right))),
    }

def add_paired_test_fields(summary: dict[str, float | str], prefix: str, left: np.ndarray, right: np.ndarray) -> None:
    test = paired_test(left, right)
    summary[f'{prefix}_delta_test'] = test['delta']
    summary[f'{prefix}_p_t'] = test['p_t']
    summary[f'{prefix}_p_w'] = test['p_w']
    summary[f'{prefix}_n'] = test['n']

def add_independent_test_fields(summary: dict[str, float | str], prefix: str, left: np.ndarray, right: np.ndarray) -> None:
    test = independent_test(left, right)
    summary[f'{prefix}_delta_test'] = test['delta']
    summary[f'{prefix}_p_t'] = test['p_t']
    summary[f'{prefix}_p_w'] = test['p_w']
    summary[f'{prefix}_n'] = test['n']

def summarize_result(result: dict[str, Any]) -> dict[str, float | str]:
    t = result['t']
    traces = result['traces']
    stats = result['stats']
    drives = result['drives']
    p = result['params']

    da_up_matched, non_da_up_matched = match_paired_pre_baseline(traces['da_up'], traces['non_da_up'], t, p.pre_window)
    da_down_matched, non_da_down_matched = match_paired_pre_baseline(traces['da_down'], traces['non_da_down'], t, p.pre_window)

    window_14 = (t >= 1.0) & (t < 4.0)
    window_04 = (t >= 0.0) & (t < 4.0)
    window_34 = (t >= 3.0) & (t < 4.0)
    window_354 = (t >= 3.5) & (t < 4.0)
    window_054 = (t >= 0.5) & (t < 4.0)
    exp1_up = np.nanmean((traces['exp1_lc_up_same'] - traces['exp1_base_up_same'])[:, window_14], axis=1)
    exp1_down = np.nanmean((traces['exp1_lc_down_same'] - traces['exp1_base_down_same'])[:, window_04], axis=1)
    exp2_up = np.nanmean((da_up_matched - non_da_up_matched)[:, window_14], axis=1)
    exp2_down = np.nanmean((da_down_matched - non_da_down_matched)[:, window_04], axis=1)
    exp3_up = np.nanmean((traces['exp3_base_up_same'] - traces['exp3_block_up_same'])[:, window_14], axis=1)
    exp3_down = np.nanmean((traces['exp3_base_down_same'] - traces['exp3_block_down_same'])[:, window_04], axis=1)
    exp1_up_late = np.nanmean((traces['exp1_lc_up_same'] - traces['exp1_base_up_same'])[:, window_34], axis=1)
    exp2_up_late = np.nanmean((da_up_matched - non_da_up_matched)[:, window_34], axis=1)
    exp3_up_late = np.nanmean((traces['exp3_base_up_same'] - traces['exp3_block_up_same'])[:, window_34], axis=1)
    exp1_up_tail = np.nanmean((traces['exp1_lc_up_same'] - traces['exp1_base_up_same'])[:, window_354], axis=1)
    exp2_up_tail = np.nanmean((da_up_matched - non_da_up_matched)[:, window_354], axis=1)
    exp3_up_tail = np.nanmean((traces['exp3_base_up_same'] - traces['exp3_block_up_same'])[:, window_354], axis=1)
    exp1_down_profile = np.nanmax(np.abs((traces['exp1_lc_down_same'] - traces['exp1_base_down_same'])[:, window_054]), axis=1)
    exp2_down_profile = np.nanmax(np.abs((da_down_matched - non_da_down_matched)[:, window_054]), axis=1)
    exp3_down_profile = np.nanmax(np.abs((traces['exp3_base_down_same'] - traces['exp3_block_down_same'])[:, window_054]), axis=1)

    d = np.asarray(drives['D'], dtype=float)
    d_release = np.asarray(drives['D_release'], dtype=float)
    base_up_mean = np.nanmean(traces['base_up'], axis=0)
    base_down_mean = np.nanmean(traces['base_down'], axis=0)
    base_up_mask = window_04 & np.isfinite(base_up_mean)
    base_down_mask = window_04 & np.isfinite(base_down_mean)
    base_up_values = base_up_mean[base_up_mask]
    base_down_values = base_down_mean[base_down_mask]
    base_up_min = float(np.min(base_up_values))
    base_up_max = float(np.max(base_up_values))
    base_down_min = float(np.min(base_down_values))
    base_down_max = float(np.max(base_down_values))
    base_up_peak_time = float(
        t[np.flatnonzero(base_up_mask)[int(np.argmax(base_up_values))]]
    )
    base_down_trough_time = float(
        t[np.flatnonzero(base_down_mask)[int(np.argmin(base_down_values))]]
    )
    summary: dict[str, float | str] = {
        'name': 'selected_model',
        'drive_form': 'rate',
        'D_peak_time': float(t[int(np.nanargmax(d))]),
        'D_release_peak_time': float(t[int(np.nanargmax(d_release))]),
        'D_at_0p5': float(np.interp(0.5, t, d)),
        'D_at_0p8': float(np.interp(0.8, t, d)),
        'D_at_1p0': float(np.interp(1.0, t, d)),
        'D_at_3p0': float(np.interp(3.0, t, d)),
        'D_at_4p0': float(np.interp(4.0, t, d)),
        'base_up_pct': float(np.nanmean(stats['base_up_pct'])),
        'base_down_pct': float(np.nanmean(stats['base_down_pct'])),
        'exp1_delta_up_pct': float(np.nanmean(stats['lc_up_pct'] - stats['base_up_pct'])),
        'exp1_delta_down_pct': float(np.nanmean(stats['lc_down_pct'] - stats['base_down_pct'])),
        'exp2_delta_up_pct': float(np.nanmean(stats['p_up_da_targeted'] - stats['p_up_not_targeted'])),
        'exp2_delta_down_pct': float(np.nanmean(stats['p_down_da_targeted'] - stats['p_down_not_targeted'])),
        'exp3_delta_up_pct': float(np.nanmean(stats['block_up_pct'] - stats['base_up_pct'])),
        'exp3_delta_down_pct': float(np.nanmean(stats['block_down_pct'] - stats['base_down_pct'])),
        'base_up_pre': float(np.nanmean(stats['base_up_pre'])),
        'base_down_pre': float(np.nanmean(stats['base_down_pre'])),
        'base_up_04': float(np.nanmean(stats['base_up_04'])),
        'base_down_04': float(np.nanmean(stats['base_down_04'])),
        'base_rate_gap_04': float(np.nanmean(stats['base_up_04'] - stats['base_down_04'])),
        'base_up_late': float(np.nanmean(stats['base_up_late'])),
        'base_down_late': float(np.nanmean(stats['base_down_late'])),
        'base_up_wW': float(np.nanmean(stats['base_up_wW'])),
        'base_down_wW': float(np.nanmean(stats['base_down_wW'])),
        'base_up_tau': float(np.nanmean(stats['base_up_tau'])),
        'base_down_tau': float(np.nanmean(stats['base_down_tau'])),
        'base_up_min_0_4': base_up_min,
        'base_up_max_0_4': base_up_max,
        'base_down_min_0_4': base_down_min,
        'base_down_max_0_4': base_down_max,
        'base_up_peak_time_0_4': base_up_peak_time,
        'base_down_trough_time_0_4': base_down_trough_time,
        'exp1_up_fr_1_4': float(np.nanmean(exp1_up)),
        'exp1_down_fr_0_4': float(np.nanmean(exp1_down)),
        'exp2_up_fr_1_4': float(np.nanmean(exp2_up)),
        'exp2_down_fr_0_4': float(np.nanmean(exp2_down)),
        'exp3_up_fr_1_4': float(np.nanmean(exp3_up)),
        'exp3_down_fr_0_4': float(np.nanmean(exp3_down)),
        'exp1_up_fr_3_4': float(np.nanmean(exp1_up_late)),
        'exp2_up_fr_3_4': float(np.nanmean(exp2_up_late)),
        'exp3_up_fr_3_4': float(np.nanmean(exp3_up_late)),
        'exp1_up_fr_3p5_4': float(np.nanmean(exp1_up_tail)),
        'exp2_up_fr_3p5_4': float(np.nanmean(exp2_up_tail)),
        'exp3_up_fr_3p5_4': float(np.nanmean(exp3_up_tail)),
        'exp1_down_profile_max_0p5_4': float(np.nanmean(exp1_down_profile)),
        'exp2_down_profile_max_0p5_4': float(np.nanmean(exp2_down_profile)),
        'exp3_down_profile_max_0p5_4': float(np.nanmean(exp3_down_profile)),
    }

    add_paired_test_fields(summary, 'exp1_up_pct', stats['base_up_pct'], stats['lc_up_pct'])
    add_paired_test_fields(summary, 'exp1_down_pct', stats['base_down_pct'], stats['lc_down_pct'])
    add_paired_test_fields(summary, 'exp2_up_pct', stats['p_up_not_targeted'], stats['p_up_da_targeted'])
    add_paired_test_fields(summary, 'exp2_down_pct', stats['p_down_not_targeted'], stats['p_down_da_targeted'])
    add_paired_test_fields(summary, 'exp3_up_pct', stats['base_up_pct'], stats['block_up_pct'])
    add_paired_test_fields(summary, 'exp3_down_pct', stats['base_down_pct'], stats['block_down_pct'])

    for label, left, right in (
        ('exp1_up', traces['exp1_base_up_same'], traces['exp1_lc_up_same']),
        ('exp2_up', non_da_up_matched, da_up_matched),
        ('exp3_up', traces['exp3_block_up_same'], traces['exp3_base_up_same']),
        ('exp1_down', traces['exp1_base_down_same'], traces['exp1_lc_down_same']),
        ('exp2_down', non_da_down_matched, da_down_matched),
        ('exp3_down', traces['exp3_block_down_same'], traces['exp3_base_down_same']),
    ):
        add_independent_test_fields(
            summary,
            f'{label}_fr_1_4',
            np.nanmean(left[:, window_14], axis=1),
            np.nanmean(right[:, window_14], axis=1),
        )
        add_independent_test_fields(
            summary,
            f'{label}_fr_3_4',
            np.nanmean(left[:, window_34], axis=1),
            np.nanmean(right[:, window_34], axis=1),
        )
        add_independent_test_fields(
            summary,
            f'{label}_fr_3p5_4',
            np.nanmean(left[:, window_354], axis=1),
            np.nanmean(right[:, window_354], axis=1),
        )
        add_independent_test_fields(
            summary,
            f'{label}_fr_0_4',
            np.nanmean(left[:, window_04], axis=1),
            np.nanmean(right[:, window_04], axis=1),
        )

    return summary

def plot_mean_sem(
    ax: plt.Axes,
    t: np.ndarray,
    traces: np.ndarray,
    color: str,
    label: str,
    linestyle: str = '-',
    alpha_fill: float = 0.18,
) -> None:
    traces = np.asarray(traces, dtype=float)
    mean_trace = np.nanmean(traces, axis=0)
    sem_trace = sem(traces, axis=0, nan_policy='omit')
    ax.plot(t, mean_trace, color=color, linewidth=2, linestyle=linestyle, label=label)
    ax.fill_between(t, mean_trace - sem_trace, mean_trace + sem_trace, color=color, alpha=alpha_fill, linewidth=0)

def trace_ylim(
    *trace_sets: np.ndarray,
    pad_frac: float = 0.08,
    lower_floor: float = 0.0,
) -> tuple[float, float]:
    lows: list[float] = []
    highs: list[float] = []
    for traces in trace_sets:
        traces = np.asarray(traces, dtype=float)
        mean_trace = np.nanmean(traces, axis=0)
        sem_trace = sem(traces, axis=0, nan_policy='omit')
        lows.append(float(np.nanmin(mean_trace - sem_trace)))
        highs.append(float(np.nanmax(mean_trace + sem_trace)))

    ymin = min(lows)
    ymax = max(highs)
    span = ymax - ymin
    pad = span * pad_frac
    return max(lower_floor, ymin - pad), ymax + pad

def paired_ylim(
    *arrays: np.ndarray,
    lower_floor: float = 0.0,
    pad_frac: float = 0.35,
    min_pad: float = 2.5,
    anchor_floor: float | None = None,
) -> tuple[float, float]:
    values = np.concatenate([
        np.asarray(arr, dtype=float).ravel()
        for arr in arrays
    ])
    vmin = float(np.min(values))
    vmax = float(np.max(values))
    span = vmax - vmin
    pad = max(span * pad_frac, min_pad)
    ymin = vmin - pad
    if anchor_floor is not None:
        ymin = min(float(anchor_floor), ymin)
    return max(lower_floor, ymin), vmax + pad

def clean_axis(ax: plt.Axes) -> None:
    ax.set_facecolor('white')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=4, integer=True, min_n_ticks=1))

def summary_stat_text(
    summary: dict[str, float | str],
    prefix: str,
    comparison: str,
    unit: str,
    test_label: str,
    precision: int = 3,
) -> str:
    delta = float(summary[f'{prefix}_delta_test'])
    p_value = float(summary[f'{prefix}_p_w'])
    delta_text = f'Delta={delta:+.{precision}f} {unit}'
    if p_value < 1e-4:
        p_text = 'p<1e-4'
    else:
        p_text = f'p={p_value:.3g}'
    return f'{comparison}\n{delta_text}; {test_label} {p_text}'

def add_trace_statistics(ax: plt.Axes, text: str) -> None:
    ax.text(
        0.98,
        0.96,
        text,
        transform=ax.transAxes,
        ha='right',
        va='top',
        fontsize=8,
        bbox={'facecolor': 'white', 'edgecolor': 'none', 'alpha': 0.78, 'pad': 2.0},
        zorder=10,
    )

def save_figure_bundle(fig: plt.Figure, output_dir: Path, stem: str) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    for ext in ('.png', '.pdf'):
        fig.savefig(output_dir / f'{stem}{ext}', dpi=300, bbox_inches='tight', facecolor='white')

def build_experiment_axes(figsize: tuple[float, float] = (5.6, 4.8)):
    fig = plt.figure(figsize=figsize, constrained_layout=True)
    fig.set_facecolor('white')
    outer = fig.add_gridspec(
        2,
        2,
        width_ratios=[2.4, 1.3],
        height_ratios=[1, 1],
        wspace=0.45,
        hspace=0.35,
    )
    ax_up = fig.add_subplot(outer[0, 0])
    ax_down = fig.add_subplot(outer[1, 0])
    ax_up_bar = fig.add_subplot(outer[0, 1])
    ax_down_bar = fig.add_subplot(outer[1, 1])
    for ax in (ax_up, ax_down, ax_up_bar, ax_down_bar):
        ax.set_facecolor('white')
    return fig, ax_up, ax_down, ax_up_bar, ax_down_bar

def plot_trace_pair(
    ax: plt.Axes,
    t: np.ndarray,
    left_traces: np.ndarray,
    right_traces: np.ndarray,
    left_color: str,
    right_color: str,
    left_label: str,
    right_label: str,
    title: str,
    show_xlabel: bool = False,
    lower_floor: float = 0.0,
    ylim: tuple[float, float] | None = None,
) -> None:
    plot_mean_sem(ax, t, left_traces, left_color, left_label)
    plot_mean_sem(ax, t, right_traces, right_color, right_label)
    ax.axvline(0, linestyle='--', color=RED, linewidth=1)
    ax.set_xlim([-1.0, 4.0])
    ax.set_ylabel('Firing rate (Hz)')
    ax.set_title(title)
    if ylim is None:
        ylim = trace_ylim(left_traces, right_traces, pad_frac=0.10, lower_floor=lower_floor)
    ax.set_ylim(ylim)
    ax.legend(frameon=False, fontsize=10)
    clean_axis(ax)
    if show_xlabel:
        ax.set_xlabel('Time from run onset (s)')
    else:
        ax.set_xticklabels([])

def plot_paired_bar_panel(
    ax: plt.Axes,
    left_values: np.ndarray,
    right_values: np.ndarray,
    left_color: str,
    right_color: str,
    title: str,
    xticklabels: tuple[str, str],
    anchor_floor: float,
    stats_text: str | None = None,
    ylim: tuple[float, float] | None = None,
) -> None:
    left_values = np.asarray(left_values, dtype=float)
    right_values = np.asarray(right_values, dtype=float)
    if ylim is None:
        ylim = paired_ylim(left_values, right_values, min_pad=3.0, anchor_floor=anchor_floor)
    left = left_values
    right = right_values

    means = [np.nanmean(left), np.nanmean(right)]
    errs = [float(sem(left, nan_policy='omit')), float(sem(right, nan_policy='omit'))]
    colors = (left_color, right_color)
    ax.bar(
        [0, 1],
        means,
        width=0.8,
        color=[to_rgba(color, 0.6) for color in colors],
        edgecolor='black',
        linewidth=1.0,
        zorder=2,
    )
    for y0, y1 in zip(left, right):
        ax.plot([0, 1], [y0, y1], lw=0.8, color='k', alpha=0.3, zorder=3)
    ax.scatter(np.zeros(len(left)), left, s=10, color=left_color, edgecolor='none', alpha=0.5, zorder=4)
    ax.scatter(np.ones(len(right)), right, s=10, color=right_color, edgecolor='none', alpha=0.5, zorder=4)
    ax.errorbar(
        [0, 1],
        means,
        yerr=errs,
        fmt='none',
        ecolor='k',
        elinewidth=1.2,
        capsize=3,
        capthick=1.0,
        zorder=6,
    )

    ax.set(
        xticks=[0, 1],
        xticklabels=xticklabels,
        ylabel='Proportion (%)',
        title=title,
        xlim=(-0.5, 1.5),
        ylim=ylim,
    )
    clean_axis(ax)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=4, integer=True, min_n_ticks=3))
    if stats_text:
        y0, y1 = ax.get_ylim()
        span = y1 - y0
        ax.text(
            0.5,
            y1 - 0.12 * span,
            stats_text,
            ha='center',
            va='top',
            fontsize=7.5,
            linespacing=1.05,
        )

def plot_overview(result: dict[str, Any], summary: dict[str, float | str], output_dir: Path) -> None:
    t = np.asarray(result['t'], dtype=float)
    p = result['params']
    traces = result['traces']
    drives = result['drives']

    r_grid = np.linspace(0.0, 8.0, 300)
    da_gate = sigmoid(r_grid, p.da_half_rate, p.da_rate_slope)

    fig = plt.figure(figsize=(7.4, 7.2), constrained_layout=True)
    fig.set_facecolor('white')
    gs = fig.add_gridspec(3, 2, height_ratios=[1.0, 1.0, 1.05], hspace=0.20, wspace=0.25)

    ax = fig.add_subplot(gs[0, 0])
    for key, color, label in (
        ('R', GOLDENROD, 'Run R'),
        ('W', MAGENTA, 'Task W'),
        ('D_release', DARKGREEN, 'DA release'),
        ('D', ROYALBLUE, 'CA1 DA D'),
    ):
        ax.plot(t, drives[key], color=color, linewidth=2.3, label=label)
    ax.axvline(0, linestyle='--', color=RED, linewidth=1)
    ax.set_xlim([-1.0, 4.0])
    ax.set_ylabel('Drive value')
    ax.set_title('Shared drives')
    ax.legend(frameon=False, fontsize=9, ncol=2)
    clean_axis(ax)

    ax = fig.add_subplot(gs[0, 1])
    ax.plot(r_grid, da_gate, color=DARKGREEN, linewidth=2.4)
    ax.axvline(p.da_half_rate, linestyle='--', color='0.55', linewidth=1)
    ax.set_xlim([0.0, 4.0])
    ax.set_ylim([0.0, 1.05])
    ax.set_xlabel('Previous firing rate (Hz)')
    ax.set_ylabel('Extra DA weight')
    ax.set_title('Activity-dependent targeted DA weight')
    clean_axis(ax)

    ax = fig.add_subplot(gs[1, 0])
    plot_mean_sem(ax, t, traces['base_up'], CLASS_COLORS['is_up'], 'PyrUp')
    plot_mean_sem(ax, t, traces['base_down'], CLASS_COLORS['is_down'], 'PyrDown')
    ax.axvline(0, linestyle='--', color=RED, linewidth=1)
    ax.set_xlim([-1.0, 4.0])
    overview_ylim = trace_ylim(traces['base_up'], traces['base_down'], pad_frac=0.10, lower_floor=0.0)
    ax.set_ylim(overview_ylim)
    ax.set_xlabel('Time from run onset (s)')
    ax.set_ylabel('Firing rate (Hz)')
    ax.set_title('Baseline class profiles')
    ax.legend(frameon=False, fontsize=10, loc='upper right')
    clean_axis(ax)

    ax = fig.add_subplot(gs[1, 1])
    stats = result['stats']
    means = [np.nanmean(stats['base_up_pct']), np.nanmean(stats['base_down_pct'])]
    errs = [sem(stats['base_up_pct'], axis=0, nan_policy='omit'), sem(stats['base_down_pct'], axis=0, nan_policy='omit')]
    ax.bar(
        [0, 1],
        means,
        yerr=errs,
        capsize=2,
        color=[to_rgba(CLASS_COLORS['is_up'], 0.65), to_rgba(CLASS_COLORS['is_down'], 0.65)],
        edgecolor='black',
        linewidth=1.0,
    )
    ax.set_xticks([0, 1])
    ax.set_xticklabels(['PyrUp', 'PyrDown'])
    ax.set_ylabel('Cells (%)')
    ax.set_title('Baseline class proportions')
    ax.set_ylim(paired_ylim(stats['base_up_pct'], stats['base_down_pct'], min_pad=4.0, anchor_floor=0.0))
    ax.text(
        0.5,
        0.96,
        f'PyrUp {means[0]:.2f}+/-{errs[0]:.2f}%\nPyrDown {means[1]:.2f}+/-{errs[1]:.2f}%',
        transform=ax.transAxes,
        ha='center',
        va='top',
        fontsize=8,
        bbox={'facecolor': 'white', 'edgecolor': 'none', 'alpha': 0.78, 'pad': 2.0},
    )
    clean_axis(ax)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=4, integer=True, min_n_ticks=3))

    ax = fig.add_subplot(gs[2, :])
    text = (
        f"Baseline PyrUp/PyrDown: {summary['base_up_pct']:.2f}% / {summary['base_down_pct']:.2f}%\n"
        'PyrUp late FR effects, 3-4 s: '
        f"Exp1 {summary['exp1_up_fr_3_4']:+.3f}, "
        f"Exp2 {summary['exp2_up_fr_3_4']:+.3f}, "
        f"Exp3 {summary['exp3_up_fr_3_4']:+.3f} Hz\n"
        'PyrDown 0-4 s rank-sum p: '
        f"Exp1 {summary['exp1_down_fr_0_4_p_w']:.3g}, "
        f"Exp2 {summary['exp2_down_fr_0_4_p_w']:.3g}, "
        f"Exp3 {summary['exp3_down_fr_0_4_p_w']:.3g}"
    )
    ax.text(0.02, 0.80, text, transform=ax.transAxes, family='monospace', fontsize=10, va='top')
    ax.set_axis_off()

    save_figure_bundle(fig, output_dir, 'fig_6_lc_da_ca1_model_overview')
    plt.close(fig)

def experiment_shared_ylims(result: dict[str, Any]) -> dict[str, tuple[float, float]]:
    t = result['t']
    traces = result['traces']
    stats = result['stats']
    p = result['params']

    exp1_base_up, exp1_lc_up = match_paired_pre_baseline(
        traces['exp1_base_up_same'], traces['exp1_lc_up_same'], t, p.pre_window)
    exp1_base_down, exp1_lc_down = match_paired_pre_baseline(
        traces['exp1_base_down_same'], traces['exp1_lc_down_same'], t, p.pre_window)
    exp2_da_up, exp2_non_da_up = match_paired_pre_baseline(
        traces['da_up'], traces['non_da_up'], t, p.pre_window)
    exp2_da_down, exp2_non_da_down = match_paired_pre_baseline(
        traces['da_down'], traces['non_da_down'], t, p.pre_window)

    return {
        'profile': trace_ylim(
            exp1_base_up,
            exp1_lc_up,
            exp1_base_down,
            exp1_lc_down,
            exp2_da_up,
            exp2_non_da_up,
            exp2_da_down,
            exp2_non_da_down,
            traces['exp3_base_up_same'],
            traces['exp3_block_up_same'],
            traces['exp3_base_down_same'],
            traces['exp3_block_down_same'],
            pad_frac=0.10,
            lower_floor=0.0,
        ),
        'proportion': paired_ylim(
            stats['base_up_pct'],
            stats['lc_up_pct'],
            stats['base_down_pct'],
            stats['lc_down_pct'],
            stats['p_up_not_targeted'],
            stats['p_up_da_targeted'],
            stats['p_down_not_targeted'],
            stats['p_down_da_targeted'],
            stats['block_up_pct'],
            stats['block_down_pct'],
            min_pad=3.0,
            anchor_floor=0.0,
        ),
    }

def plot_experiment_1(
    result: dict[str, Any],
    summary: dict[str, float | str],
    output_dir: Path,
    ylims: dict[str, tuple[float, float]] | None = None,
) -> None:
    t = result['t']
    traces = result['traces']
    stats = result['stats']
    p = result['params']
    fig, ax_up, ax_down, ax_up_bar, ax_down_bar = build_experiment_axes()

    base_up_plot, lc_up_plot = match_paired_pre_baseline(traces['exp1_base_up_same'], traces['exp1_lc_up_same'], t, p.pre_window)
    base_down_plot, lc_down_plot = match_paired_pre_baseline(
        traces['exp1_base_down_same'], traces['exp1_lc_down_same'], t, p.pre_window
    )
    profile_ylim = trace_ylim(
        base_up_plot,
        lc_up_plot,
        base_down_plot,
        lc_down_plot,
        pad_frac=0.10,
        lower_floor=0.0,
    )
    proportion_ylim = paired_ylim(
        stats['base_up_pct'],
        stats['lc_up_pct'],
        stats['base_down_pct'],
        stats['lc_down_pct'],
        min_pad=3.0,
        anchor_floor=10.0,
    )
    if ylims is not None:
        profile_ylim = ylims['profile']
        proportion_ylim = ylims['proportion']

    plot_trace_pair(
        ax_up,
        t,
        base_up_plot,
        lc_up_plot,
        CLASS_COLORS['is_up'],
        CONDITION_COLORS['lc'],
        'Ctrl.',
        'Stim.',
        'PyrUp',
        ylim=profile_ylim,
    )
    add_trace_statistics(ax_up, summary_stat_text(summary, 'exp1_up_fr_1_4', '1-4 s Stim.-Ctrl.', 'Hz', 'rank-sum'))
    plot_trace_pair(
        ax_down,
        t,
        base_down_plot,
        lc_down_plot,
        CLASS_COLORS['is_down'],
        CONDITION_COLORS['lc'],
        'Ctrl.',
        'Stim.',
        'PyrDown',
        show_xlabel=True,
        ylim=profile_ylim,
    )
    add_trace_statistics(ax_down, summary_stat_text(summary, 'exp1_down_fr_0_4', '0-4 s Stim.-Ctrl.', 'Hz', 'rank-sum'))
    plot_paired_bar_panel(
        ax_up_bar,
        stats['base_up_pct'],
        stats['lc_up_pct'],
        CLASS_COLORS['is_up'],
        CONDITION_COLORS['lc'],
        'PyrUp proportion',
        ('Ctrl.', 'Stim.'),
        10.0,
        stats_text=summary_stat_text(summary, 'exp1_up_pct', 'Stim.-Ctrl.', 'pp', 'Wilcoxon', precision=2),
        ylim=proportion_ylim,
    )
    plot_paired_bar_panel(
        ax_down_bar,
        stats['base_down_pct'],
        stats['lc_down_pct'],
        CLASS_COLORS['is_down'],
        CONDITION_COLORS['lc'],
        'PyrDown proportion',
        ('Ctrl.', 'Stim.'),
        10.0,
        stats_text=summary_stat_text(summary, 'exp1_down_pct', 'Stim.-Ctrl.', 'pp', 'Wilcoxon', precision=2),
        ylim=proportion_ylim,
    )
    save_figure_bundle(fig, output_dir, 'fig_6_lc_da_ca1_model_exp_1_lc_activation')
    plt.close(fig)

def plot_experiment_2(
    result: dict[str, Any],
    summary: dict[str, float | str],
    output_dir: Path,
    ylims: dict[str, tuple[float, float]] | None = None,
) -> None:
    t = result['t']
    traces = result['traces']
    stats = result['stats']
    p = result['params']
    fig, ax_up, ax_down, ax_up_bar, ax_down_bar = build_experiment_axes()

    da_up_plot, non_da_up_plot = match_paired_pre_baseline(traces['da_up'], traces['non_da_up'], t, p.pre_window)
    da_down_plot, non_da_down_plot = match_paired_pre_baseline(traces['da_down'], traces['non_da_down'], t, p.pre_window)
    profile_ylim = trace_ylim(
        da_up_plot,
        non_da_up_plot,
        da_down_plot,
        non_da_down_plot,
        pad_frac=0.10,
        lower_floor=0.0,
    )
    proportion_ylim = paired_ylim(
        stats['p_up_not_targeted'],
        stats['p_up_da_targeted'],
        stats['p_down_not_targeted'],
        stats['p_down_da_targeted'],
        min_pad=3.0,
        anchor_floor=0.0,
    )
    if ylims is not None:
        profile_ylim = ylims['profile']
        proportion_ylim = ylims['proportion']

    up_non_targ_color = LIGHTCORAL
    down_non_targ_color = PLUM
    plot_trace_pair(
        ax_up,
        t,
        da_up_plot,
        non_da_up_plot,
        CLASS_COLORS['is_up'],
        up_non_targ_color,
        'DA-Up',
        'non-DA-Up',
        'PyrUp',
        ylim=profile_ylim,
    )
    add_trace_statistics(ax_up, summary_stat_text(summary, 'exp2_up_fr_1_4', '1-4 s DA-Up - non-DA-Up', 'Hz', 'rank-sum'))
    plot_trace_pair(
        ax_down,
        t,
        da_down_plot,
        non_da_down_plot,
        CLASS_COLORS['is_down'],
        down_non_targ_color,
        'DA-Up',
        'non-DA-Up',
        'PyrDown',
        show_xlabel=True,
        ylim=profile_ylim,
    )
    add_trace_statistics(ax_down, summary_stat_text(summary, 'exp2_down_fr_0_4', '0-4 s DA-Up - non-DA-Up', 'Hz', 'rank-sum'))
    plot_paired_bar_panel(
        ax_up_bar,
        stats['p_up_not_targeted'],
        stats['p_up_da_targeted'],
        up_non_targ_color,
        CLASS_COLORS['is_up'],
        'PyrUp proportion',
        ('non-DA-Up', 'DA-Up'),
        10.0,
        stats_text=summary_stat_text(summary, 'exp2_up_pct', 'DA-Up - non-DA-Up', 'pp', 'Wilcoxon', precision=2),
        ylim=proportion_ylim,
    )
    plot_paired_bar_panel(
        ax_down_bar,
        stats['p_down_not_targeted'],
        stats['p_down_da_targeted'],
        down_non_targ_color,
        CLASS_COLORS['is_down'],
        'PyrDown proportion',
        ('non-DA-Up', 'DA-Up'),
        0.0,
        stats_text=summary_stat_text(summary, 'exp2_down_pct', 'DA-Up - non-DA-Up', 'pp', 'Wilcoxon', precision=2),
        ylim=proportion_ylim,
    )
    save_figure_bundle(fig, output_dir, 'fig_6_lc_da_ca1_model_exp_2_da_targeted_subsets')
    plt.close(fig)

def plot_experiment_3(
    result: dict[str, Any],
    summary: dict[str, float | str],
    output_dir: Path,
    ylims: dict[str, tuple[float, float]] | None = None,
) -> None:
    t = result['t']
    traces = result['traces']
    stats = result['stats']
    fig, ax_up, ax_down, ax_up_bar, ax_down_bar = build_experiment_axes()
    profile_ylim = trace_ylim(
        traces['exp3_base_up_same'],
        traces['exp3_block_up_same'],
        traces['exp3_base_down_same'],
        traces['exp3_block_down_same'],
        pad_frac=0.10,
        lower_floor=0.0,
    )
    proportion_ylim = paired_ylim(
        stats['base_up_pct'],
        stats['block_up_pct'],
        stats['base_down_pct'],
        stats['block_down_pct'],
        min_pad=3.0,
        anchor_floor=10.0,
    )
    if ylims is not None:
        profile_ylim = ylims['profile']
        proportion_ylim = ylims['proportion']

    plot_trace_pair(
        ax_up,
        t,
        traces['exp3_base_up_same'],
        traces['exp3_block_up_same'],
        CLASS_COLORS['is_up'],
        MIDGREY,
        'Baseline',
        'Blockade',
        'PyrUp',
        ylim=profile_ylim,
    )
    add_trace_statistics(ax_up, summary_stat_text(summary, 'exp3_up_fr_1_4', '1-4 s Baseline - Blockade', 'Hz', 'rank-sum'))
    plot_trace_pair(
        ax_down,
        t,
        traces['exp3_base_down_same'],
        traces['exp3_block_down_same'],
        CLASS_COLORS['is_down'],
        MIDGREY,
        'Baseline',
        'Blockade',
        'PyrDown',
        show_xlabel=True,
        ylim=profile_ylim,
    )
    add_trace_statistics(ax_down, summary_stat_text(summary, 'exp3_down_fr_0_4', '0-4 s Baseline - Blockade', 'Hz', 'rank-sum'))
    plot_paired_bar_panel(
        ax_up_bar,
        stats['base_up_pct'],
        stats['block_up_pct'],
        CLASS_COLORS['is_up'],
        MIDGREY,
        'PyrUp proportion',
        ('Baseline', 'Blockade'),
        10.0,
        stats_text=summary_stat_text(summary, 'exp3_up_pct', 'Blockade - Baseline', 'pp', 'Wilcoxon', precision=2),
        ylim=proportion_ylim,
    )
    plot_paired_bar_panel(
        ax_down_bar,
        stats['base_down_pct'],
        stats['block_down_pct'],
        CLASS_COLORS['is_down'],
        MIDGREY,
        'PyrDown proportion',
        ('Baseline', 'Blockade'),
        10.0,
        stats_text=summary_stat_text(summary, 'exp3_down_pct', 'Blockade - Baseline', 'pp', 'Wilcoxon', precision=2),
        ylim=proportion_ylim,
    )
    save_figure_bundle(fig, output_dir, 'fig_6_lc_da_ca1_model_exp_3_partial_da_block')
    plt.close(fig)

def quartile_text(values: np.ndarray, precision: int = 3) -> str:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    q1, med, q3 = np.percentile(values, [25, 50, 75])
    return f'[{q1:.{precision}f} {med:.{precision}f} {q3:.{precision}f}]'

def print_summary(result: dict[str, Any], summary: dict[str, float | str], p: Params) -> None:
    t = result['t']
    traces = result['traces']
    stats = result['stats']

    da_up_matched, non_da_up_matched = match_paired_pre_baseline(
        traces['da_up'], traces['non_da_up'], t, p.pre_window)
    da_down_matched, non_da_down_matched = match_paired_pre_baseline(
        traces['da_down'], traces['non_da_down'], t, p.pre_window)
    window_04 = (t >= 0.0) & (t < 4.0)

    fr_readout = [
        (
            'Exp. 1 LC stim',
            'Ctrl.',
            np.nanmean(traces['exp1_base_up_same'][:, window_04], axis=1),
            np.nanmean(traces['exp1_base_down_same'][:, window_04], axis=1),
            'Stim.',
            np.nanmean(traces['exp1_lc_up_same'][:, window_04], axis=1),
            np.nanmean(traces['exp1_lc_down_same'][:, window_04], axis=1),
            'exp1',
        ),
        (
            'Exp. 2 DA-targeted',
            'non-DA',
            np.nanmean(non_da_up_matched[:, window_04], axis=1),
            np.nanmean(non_da_down_matched[:, window_04], axis=1),
            'DA',
            np.nanmean(da_up_matched[:, window_04], axis=1),
            np.nanmean(da_down_matched[:, window_04], axis=1),
            'exp2',
        ),
        (
            'Exp. 3 DA blockade',
            'Baseline',
            np.nanmean(traces['exp3_base_up_same'][:, window_04], axis=1),
            np.nanmean(traces['exp3_base_down_same'][:, window_04], axis=1),
            'Blockade',
            np.nanmean(traces['exp3_block_up_same'][:, window_04], axis=1),
            np.nanmean(traces['exp3_block_down_same'][:, window_04], axis=1),
            'exp3',
        ),
    ]

    prop_readout = [
        (
            'Exp. 1 LC stim',
            'Ctrl.',
            stats['base_up_pct'],
            stats['base_down_pct'],
            'Stim.',
            stats['lc_up_pct'],
            stats['lc_down_pct'],
            'exp1',
        ),
        (
            'Exp. 2 DA-targeted',
            'non-DA',
            stats['p_up_not_targeted'],
            stats['p_down_not_targeted'],
            'DA',
            stats['p_up_da_targeted'],
            stats['p_down_da_targeted'],
            'exp2',
        ),
        (
            'Exp. 3 DA blockade',
            'Baseline',
            stats['base_up_pct'],
            stats['base_down_pct'],
            'Blockade',
            stats['block_up_pct'],
            stats['block_down_pct'],
            'exp3',
        ),
    ]

    print('LC-DA-CA1 model')
    print(f'  bootstrapping = {p.n_bootstrap}')
    print(f'  cells = {p.n_cells}')
    print()

    print('0-4 s firing rate [Q1 median Q3]')
    for exp_name, left_name, left_up, left_down, right_name, right_up, right_down, prefix in fr_readout:
        print(f'  {exp_name}')
        for label, left_values, right_values, suffix in [
                ('PyrUp', left_up, right_up, 'up'),
                ('PyrDown', left_down, right_down, 'down'),
                ]:
            p_rank = float(summary[f'{prefix}_{suffix}_fr_0_4_p_w'])
            print(
                f'    {label}: '
                f'{left_name} {quartile_text(left_values)} Hz, '
                f'{right_name} {quartile_text(right_values)} Hz, '
                f'rank-sum p={p_rank:.3g}'
                )

    print()
    print('Cell proportions [Q1 median Q3]')
    for exp_name, left_name, left_up, left_down, right_name, right_up, right_down, prefix in prop_readout:
        print(f'  {exp_name}')
        for label, left_values, right_values, suffix in [
                ('PyrUp', left_up, right_up, 'up'),
                ('PyrDown', left_down, right_down, 'down'),
                ]:
            p_wilcoxon = float(summary[f'{prefix}_{suffix}_pct_p_w'])
            print(
                f'    {label}: '
                f'{left_name} {quartile_text(left_values, precision=2)}%, '
                f'{right_name} {quartile_text(right_values, precision=2)}%, '
                f'Wilcoxon p={p_wilcoxon:.3g}'
                )

def main() -> None:
    parser = argparse.ArgumentParser(description='Run the LC-DA-CA1 model.')
    parser.add_argument(
        '--output-dir',
        type=Path,
        default=Path(os.getenv('LC_DA_CA1_OUTPUT_DIR', DEFAULT_OUTPUT_DIR)),
        help='Directory for generated figures.',
    )
    parser.add_argument('--n-bootstrap', type=int, default=50, help='Paired bootstrap replicates.')
    parser.add_argument('--n-cells', type=int, default=1000, help='Synthetic cells per replicate.')
    parser.add_argument('--seed-start', type=int, default=0, help='First bootstrap RNG seed.')
    parser.add_argument('--skip-plots', action='store_true', help='Compute and print the model without writing figures.')
    args = parser.parse_args()
    output_dir = args.output_dir.expanduser().resolve()
    p = Params(n_bootstrap=args.n_bootstrap, n_cells=args.n_cells, seed_start=args.seed_start)

    result = run_bootstrap_suite(p)
    summary = summarize_result(result)

    print_summary(result, summary, p)

    if not args.skip_plots:
        experiment_ylims = experiment_shared_ylims(result)
        plot_overview(result, summary, output_dir)
        plot_experiment_1(result, summary, output_dir, experiment_ylims)
        plot_experiment_2(result, summary, output_dir, experiment_ylims)
        plot_experiment_3(result, summary, output_dir, experiment_ylims)

if __name__ == '__main__':
    main()
