# -*- coding: utf-8 -*-
'''
Created on Wed Nov 11 15:14:01 2026

compare early and late lick trials with Jingyu's response test

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
from scipy.stats import ttest_ind

repo_root = Path(__file__).resolve().parents[3]
JINGYU_ROOT = pp.JINGYU_CODE_ROOT
RAW_ROOT = pp.DBH_DLIGHT_RAW_ROOT
OUTPUT_ROOT = pp.DLIGHT_HOTSPOTS_OUTPUT_ROOT
FIGURE_ROOT = pp.DLIGHT_HOTSPOTS_FIGURE_ROOT
if str(JINGYU_ROOT) not in sys.path:
    sys.path.insert(0, str(JINGYU_ROOT))
recnames = rec_list.recDLightDbh
from common.utils_basic import trace_filter  # noqa: E402
from grid_response_test import postprocess_dff, quantify_grid_response, classify_up  # noqa: E402

N_DRAWS = 20
RNG_SEED = 42
EARLY_LIMIT = 2500
LATE_LIMIT = 4000
SPEED_SAMPLES = 4000
SPEED_BINS = 8
SPEED_BIN_LENGTH = 500
SPEED_TOLERANCE = 1.5


#%% loading and trial selection
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
    traces = dff[coords[:, 0], coords[:, 1], :]
    return behaviour, profile, traces


def valid_trials(behaviour):
    reward = np.asarray(behaviour['reward_times'], dtype=float)
    non_stop = np.asarray(behaviour['non_stop_trials'])
    non_fullstop = np.asarray(behaviour['non_fullstop_trials'])
    run = np.asarray(behaviour['run_onset_frames'], dtype=float)
    n = min(len(reward), len(non_stop), len(non_fullstop), len(run))
    valid = (
        np.isfinite(reward[:n])
        & (non_stop[:n] == 0)
        & (non_fullstop[:n] == 0)
    )
    first_lick = np.full(n, np.nan)
    for trial, licks in enumerate(behaviour['lick_times_aligned'][:n]):
        values = np.asarray(licks, dtype=float).reshape(-1)
        values = values[np.isfinite(values) & (values > 500)]
        if len(values):
            first_lick[trial] = values[0]
    early = np.flatnonzero(valid & (first_lick < EARLY_LIMIT)).tolist()
    late = np.flatnonzero(valid & (first_lick > EARLY_LIMIT) & (first_lick < LATE_LIMIT)).tolist()
    return valid, early, late, first_lick


def speed_vectors(behaviour, n):
    vectors = {}
    for trial, speed in enumerate(behaviour['speed_times_aligned'][:n]):
        values = np.asarray(speed, dtype=float)
        if values.ndim != 2 or values.shape[1] < 2 or len(values) < SPEED_SAMPLES:
            continue
        vectors[trial] = values[:SPEED_SAMPLES, 1].reshape(SPEED_BINS, SPEED_BIN_LENGTH).mean(axis=1)
    return vectors


def speed_match(early, late, vectors):
    early_available = [trial for trial in early if trial in vectors]
    late_available = [trial for trial in late if trial in vectors]
    if not early_available or not late_available:
        return [], [], 'no speed-complete trials'
    early_speed = np.vstack([vectors[trial] for trial in early_available])
    late_speed = np.vstack([vectors[trial] for trial in late_available])
    early_mean, early_sd = np.mean(early_speed, axis=0), np.std(early_speed, axis=0)
    late_mean, late_sd = np.mean(late_speed, axis=0), np.std(late_speed, axis=0)
    early_keep = np.all(
        (early_speed >= late_mean - SPEED_TOLERANCE * late_sd)
        & (early_speed <= late_mean + SPEED_TOLERANCE * late_sd), axis=1)
    late_keep = np.all(
        (late_speed >= early_mean - SPEED_TOLERANCE * early_sd)
        & (late_speed <= early_mean + SPEED_TOLERANCE * early_sd), axis=1)
    early_matched = [trial for trial, keep in zip(early_available, early_keep) if keep]
    late_matched = [trial for trial, keep in zip(late_available, late_keep) if keep]
    if len(early_matched) < 10 or len(late_matched) < 10:
        return early_matched, late_matched, 'fewer than 10 matched trials in one group'
    early_test = np.vstack([vectors[trial] for trial in early_matched])
    late_test = np.vstack([vectors[trial] for trial in late_matched])
    p_values = np.array([ttest_ind(early_test[:, i], late_test[:, i]).pvalue for i in range(SPEED_BINS)])
    if np.any(p_values < 0.05):
        return early_matched, late_matched, 'speed-bin t-test p < 0.05'
    return early_matched, late_matched, ''


#%% condition test
def condition_draws(traces, profile, run_frames, early, late, rng):
    n = min(len(early), len(late))
    if n == 0:
        return []
    draws = []
    for draw in range(N_DRAWS):
        if len(early) <= len(late):
            early_draw = list(early)
            late_draw = list(rng.choice(late, size=n, replace=False))
        else:
            early_draw = list(rng.choice(early, size=n, replace=False))
            late_draw = list(late)
        values = {}
        for condition, trials in (('early', early_draw), ('late', late_draw)):
            frames = np.asarray(run_frames, dtype=int)[trials]
            result = quantify_grid_response(traces, frames, return_single_trial=False)
            up, _ = classify_up(
                result['amplitude'], result['effect_size'],
                result['shuffle_amplitude'], profile['edge'].to_numpy())
            values[f'{condition}_pct_up'] = 100 * float(np.mean(up))
            values[f'{condition}_median_amplitude'] = float(np.nanmedian(result['amplitude']))
        values['draw'] = draw
        draws.append(values)
    return draws


#%% figure
def make_figure(summary):
    fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.5), sharey=True)
    for ax, variant, label in zip(axes, ('speed', 'unmatched'), ('speed matched', 'unmatched')):
        included = summary.loc[summary[f'{variant}_included']]
        for _, row in included.iterrows():
            ax.plot([0, 1], [row[f'{variant}_early_pct_up'], row[f'{variant}_late_pct_up']], color='0.75', linewidth=0.7)
        ax.scatter(np.zeros(len(included)), included[f'{variant}_early_pct_up'], color='0.35', s=10, label='early')
        ax.scatter(np.ones(len(included)), included[f'{variant}_late_pct_up'], color='#19a819', s=10, label='late')
        ax.set_xticks([0, 1])
        ax.set_xticklabels(['early', 'late'])
        ax.set_xlabel(label)
        ax.set_ylabel('Up ACRs (%)')
        ax.spines[['top', 'right']].set_visible(False)
    axes[0].legend(frameon=False)
    fig.tight_layout()
    path = FIGURE_ROOT / 'early_late_pct_up.png'
    fig.savefig(path, dpi=200, facecolor='white')
    plt.close(fig)


#%% run
def main():
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    FIGURE_ROOT.mkdir(parents=True, exist_ok=True)
    rng = np.random.RandomState(RNG_SEED)
    summary_rows = []
    draw_rows = []
    exclusions = {'speed': {}, 'unmatched': {}}
    for recname in list(recnames):
        behaviour, profile, traces = load_recording(recname)
        valid, early, late, first_lick = valid_trials(behaviour)
        vectors = speed_vectors(behaviour, len(valid))
        early_speed, late_speed, speed_reason = speed_match(early, late, vectors)
        if not early or not late:
            speed_reason = 'no early or late trials'
        if speed_reason:
            exclusions['speed'][speed_reason] = exclusions['speed'].get(speed_reason, 0) + 1
        variants = {
            'speed': (early_speed, late_speed, not bool(speed_reason), speed_reason),
            'unmatched': (early, late, bool(early and late), 'no early or late trials' if not early or not late else ''),
        }
        row = {'rec_id': recname, 'animal': recname.split('-')[0], 'n_early': len(early), 'n_late': len(late), 'n_early_speed': len(early_speed), 'n_late_speed': len(late_speed), 'speed_included': variants['speed'][2], 'unmatched_included': variants['unmatched'][2], 'speed_exclusion': variants['speed'][3], 'unmatched_exclusion': variants['unmatched'][3]}
        for variant, (early_group, late_group, included, reason) in variants.items():
            if not included:
                if variant == 'unmatched':
                    exclusions['unmatched'][reason] = exclusions['unmatched'].get(reason, 0) + 1
                for key in ('early_pct_up', 'late_pct_up', 'early_median_amplitude', 'late_median_amplitude'):
                    row[f'{variant}_{key}'] = np.nan
                continue
            draws = condition_draws(traces, profile, behaviour['run_onset_frames'], early_group, late_group, rng)
            draw_table = pd.DataFrame(draws)
            for _, draw in draw_table.iterrows():
                draw_rows.append({'rec_id': recname, 'animal': recname.split('-')[0], 'variant': variant, **draw.to_dict()})
            for key in ('early_pct_up', 'late_pct_up', 'early_median_amplitude', 'late_median_amplitude'):
                row[f'{variant}_{key}'] = float(draw_table[key].mean())
        summary_rows.append(row)
        print(f'{recname}: early={len(early)}, late={len(late)}, speed={len(early_speed)}/{len(late_speed)}')
    summary = pd.DataFrame(summary_rows)
    draws = pd.DataFrame(draw_rows)
    summary.to_parquet(OUTPUT_ROOT / 'aligned_early_late_per_recording.parquet', index=False)
    draws.to_parquet(OUTPUT_ROOT / 'aligned_early_late_draws.parquet', index=False)
    make_figure(summary)


if __name__ == '__main__':
    main()
