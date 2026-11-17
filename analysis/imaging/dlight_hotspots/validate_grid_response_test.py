# -*- coding: utf-8 -*-
'''
Created on Tue Nov 10 22:13:57 2026

validate the cpu grid response test against saved dLight profiles

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

import numpy as np
import pandas as pd

repo_root = Path(__file__).resolve().parents[3]
JINGYU_ROOT = pp.JINGYU_CODE_ROOT
RAW_ROOT = pp.DBH_DLIGHT_RAW_ROOT
OUTPUT_ROOT = pp.DLIGHT_HOTSPOTS_OUTPUT_ROOT
if str(JINGYU_ROOT) not in sys.path:
    sys.path.insert(0, str(JINGYU_ROOT))

recnames = rec_list.recDLightDbh
from common.utils_basic import trace_filter  # noqa: E402
from grid_response_test import postprocess_dff, quantify_grid_response, classify_up, explicit_shift_amplitude  # noqa: E402


#%% loaders
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


def profile_path(recname):
    return (
        RAW_ROOT / 'processed_dataframe_grid_free_dilation'
        / f'{recname}_profile_combined_dilation=0_pre(-1, 0)_post(0, 1)'
        '_ES=0.05_shuff95_0621.parquet'
    )


def load_inputs(recname):
    dff_path = (
        RAW_ROOT / 'regression_res_grid_free_dilation' / recname
        / 'single_trial_regression' / 'dilation_k=0'
        / 'dff_corrected_dlight.npy'
    )
    dff = postprocess_dff(np.load(dff_path, mmap_mode='r'), trace_filter)
    behaviour = load_behaviour(RAW_ROOT / 'behaviour_profile' / f'{recname}.pkl')
    profile = pd.read_parquet(profile_path(recname))
    valid = (
        profile['dlight_valid'].astype(bool)
        & profile['red_valid'].astype(bool)
        & ~profile['edge'].astype(bool)
    )
    profile = profile.loc[valid].copy()
    coords = np.asarray(profile['roi_id'].tolist(), dtype=int)
    traces = dff[coords[:, 0], coords[:, 1], :]
    return traces, behaviour, profile, coords


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
    return run[:n][valid].astype(int), n, (len(run), len(reward), len(non_stop), len(non_fullstop))


def max_relative_difference(observed, expected):
    observed = np.asarray(observed, dtype=float)
    expected = np.asarray(expected, dtype=float)
    denominator = np.maximum(np.abs(observed), 1e-12)
    return float(np.nanmax(np.abs(expected - observed) / denominator))


#%% validation
def main():
    recordings = list(recnames)
    if len(recordings) != 57:
        raise ValueError(f'Expected 57 recordings, found {len(recordings)}')
    recording_rows = []
    acr_rows = []
    hard_failures = []
    shuffle_diffs = []
    explicit_shift_checked = False
    for recname in recordings:
        traces, behaviour, profile, coords = load_inputs(recname)
        frames, n_used, lengths = event_frames(behaviour)
        if len(set(lengths)) != 1:
            print(f'{recname}: trial-array lengths differ: {lengths}')
        result = quantify_grid_response(traces, frames, return_single_trial=False)
        up_cpu, threshold_cpu = classify_up(
            result['amplitude'], result['effect_size'],
            result['shuffle_amplitude'], profile['edge'].to_numpy())
        saved_amp = profile['response_amplitude'].to_numpy(dtype=float)
        saved_effect = profile['effect_size'].to_numpy(dtype=float)
        saved_keep = profile['n_keep_trial'].to_numpy(dtype=float)
        saved_profile = np.stack(profile['mean_profile'].to_numpy())
        saved_threshold = profile['shuffle_amps_thresh_up'].to_numpy(dtype=float)
        amp_ok = np.allclose(saved_amp, result['amplitude'], rtol=1e-4, atol=1e-8, equal_nan=True)
        effect_ok = np.allclose(saved_effect, result['effect_size'], rtol=1e-4, atol=1e-8, equal_nan=True)
        keep_ok = np.array_equal(saved_keep.astype(int), result['keep'].sum(axis=1).astype(int))
        profile_ok = np.allclose(saved_profile, result['mean_profile'], rtol=1e-4, atol=1e-8, equal_nan=True)
        if not (keep_ok and amp_ok and effect_ok and profile_ok):
            hard_failures.append({
                'rec_id': recname,
                'n_keep_ok': keep_ok,
                'amplitude_ok': amp_ok,
                'effect_size_ok': effect_ok,
                'mean_profile_ok': profile_ok,
                'amplitude_max_relative_difference': max_relative_difference(saved_amp, result['amplitude']),
                'effect_size_max_relative_difference': max_relative_difference(saved_effect, result['effect_size']),
            })
        threshold_relative = np.abs(threshold_cpu - saved_threshold) / np.maximum(np.abs(saved_threshold), 1e-12)
        shuffle_diffs.extend(threshold_relative[np.isfinite(threshold_relative)].tolist())
        up_agreement = np.mean(up_cpu == profile['Up'].to_numpy(dtype=bool))
        recording_rows.append({
            'rec_id': recname,
            'animal': recname.split('-')[0],
            'n_events': len(frames),
            'n_valid_acr': len(profile),
            'n_keep_identical': keep_ok,
            'response_amplitude_identical': amp_ok,
            'effect_size_identical': effect_ok,
            'mean_profile_identical': profile_ok,
            'shuffle_threshold_median_relative_difference': float(np.nanmedian(threshold_relative)),
            'shuffle_threshold_max_relative_difference': float(np.nanmax(threshold_relative)),
            'up_agreement': float(up_agreement),
            'saved_up_pct': float(100 * profile['Up'].mean()),
            'cpu_up_pct': float(100 * up_cpu.mean()),
        })
        for roi, (y, x) in enumerate(coords):
            acr_rows.append({
                'rec_id': recname,
                'animal': recname.split('-')[0],
                'roi_y': int(y),
                'roi_x': int(x),
                'saved_n_keep_trial': int(saved_keep[roi]),
                'cpu_n_keep_trial': int(result['keep'][roi].sum()),
                'saved_response_amplitude': saved_amp[roi],
                'cpu_response_amplitude': result['amplitude'][roi],
                'saved_effect_size': saved_effect[roi],
                'cpu_effect_size': result['effect_size'][roi],
                'saved_shuffle_amps_thresh_up': saved_threshold[roi],
                'cpu_shuffle_amps_thresh_up': threshold_cpu[roi],
                'saved_up': bool(profile['Up'].iloc[roi]),
                'cpu_up': bool(up_cpu[roi]),
            })
        if not explicit_shift_checked:
            aligned = result['aligned']
            checked = quantify_grid_response(traces, frames, n_shuffle=5,
                                             return_single_trial=True)
            rng = np.random.RandomState(42)
            shifts = rng.randint(1, aligned.shape[-1], size=(5, len(frames)))
            explicit = np.empty_like(checked['shuffle_single_trial_amplitude'])
            for shuffle in range(5):
                for trial in range(len(frames)):
                    explicit[:, trial, shuffle] = explicit_shift_amplitude(
                        aligned[:, trial, :], shifts[shuffle, trial])
            if not np.allclose(explicit, checked['shuffle_single_trial_amplitude'],
                               rtol=1e-10, atol=1e-12, equal_nan=True):
                raise RuntimeError('Batched single-trial circular-shift verification failed')
            if not np.allclose(np.nanmean(explicit, axis=1),
                               checked['shuffle_amplitude'],
                               rtol=1e-10, atol=1e-12, equal_nan=True):
                raise RuntimeError('Batched trial-mean circular-shift verification failed')
            explicit_shift_checked = True
        print(f'{recname}: keep={keep_ok}, amplitude={amp_ok}, effect={effect_ok}, profile={profile_ok}, Up agreement={up_agreement:.4f}')
    recording_table = pd.DataFrame(recording_rows)
    acr_table = pd.DataFrame(acr_rows)
    recording_table.to_parquet(OUTPUT_ROOT / 'aligned_validation_per_recording.parquet', index=False)
    acr_table.to_parquet(OUTPUT_ROOT / 'aligned_validation_per_acr.parquet', index=False)
    if hard_failures:
        lines = ['', '## 11 November 2026: aligned CPU validation recheck stopped', '', 'The CPU counterpart was compared with Jingyu\'s saved `_0621` dataframes on the valid ACRs. The required hard-stop checks failed; the condition and consistency analyses were not run.']
        for failure in hard_failures:
            lines.append(f"- `{failure['rec_id']}`: n_keep identical={failure['n_keep_ok']}, amplitude identical={failure['amplitude_ok']} (max relative difference {failure['amplitude_max_relative_difference']:.6g}), effect size identical={failure['effect_size_ok']} (max relative difference {failure['effect_size_max_relative_difference']:.6g}), mean profile identical={failure['mean_profile_ok']}.")
        lines.append('Validation tables were written; no variants were tried.')
        (OUTPUT_ROOT / 'PROGRESS.md').open('a', encoding='utf-8').write('\n'.join(lines) + '\n')
        raise SystemExit(2)
    pooled_agreement = np.mean(acr_table['saved_up'] == acr_table['cpu_up'])
    lines = ['', '## 11 November 2026: aligned CPU validation recheck', '', 'The CPU response test matched Jingyu\'s saved `_0621` dataframes for n_keep_trial, response amplitude, effect size and mean profile within the specified tolerance for all 57 recordings. The batched single-trial and trial-mean shuffle amplitudes matched explicit circular shifts for five complete shuffle columns on the first recording (rtol=1e-10, atol=1e-12).', '', f"Shuffle-threshold relative differences across valid ACRs: median={np.nanmedian(shuffle_diffs):.6g}, maximum={np.nanmax(shuffle_diffs):.6g}; exact CPU/GPU agreement was not expected.", f"Up agreement: pooled={pooled_agreement:.6f}; per-recording minimum={recording_table['up_agreement'].min():.6f}, median={recording_table['up_agreement'].median():.6f}.", '', 'Validation tables: `outputs/da_hotspots/aligned_validation_per_recording.parquet` and `outputs/da_hotspots/aligned_validation_per_acr.parquet`.', '', 'Validation passed; the condition and consistency analyses may proceed.']
    (OUTPUT_ROOT / 'PROGRESS.md').open('a', encoding='utf-8').write('\n'.join(lines) + '\n')


if __name__ == '__main__':
    main()
