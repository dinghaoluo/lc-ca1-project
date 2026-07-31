'''
Created on 6 Jul 2026
Modified on 24 June 2026

build CA1 place-cell, time-cell, and sequence metrics; compare stimulation
with both clean flanking controls on the metric scale

@author: Dinghao Luo
'''

from __future__ import annotations

#%% imports
import argparse
from pathlib import Path
import pickle
import sys

import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d
from scipy.stats import norm, spearmanr
import statsmodels.formula.api as smf
from tqdm import tqdm

repo_root = Path(__file__).resolve().parents[3]
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))

from console_formatting import print_files_saved, print_session, print_status
from common_functions import get_trialtype_indices_from_stim_conds
import first_lick_analysis_functions as flaf
import hpc_ephys_support as support
import project_paths as pp

import rec_list


#%% paths and parameters
OUTPUT_STEM = pp.DATA_ROOT / 'analysis' / 'hpc' / 'sequence_analysis'
DEFAULT_PAYLOAD_PATH = OUTPUT_STEM / 'hpc_sequence_place_time_metrics.pkl'
DEFAULT_CELL_CSV_PATH = OUTPUT_STEM / 'hpc_sequence_place_time_cell_metrics.csv'
DEFAULT_SESSION_CSV_PATH = OUTPUT_STEM / 'hpc_sequence_metrics_by_session.csv'
DEFAULT_MIXED_CSV_PATH = OUTPUT_STEM / 'hpc_sequence_stim_mixed_effects.csv'

SPATIAL_BIN_SIZE_CM = 0.1
SPATIAL_PAYLOAD_DOWNSAMPLE = 10
TEMPORAL_BIN_SIZE_S = 0.05
TEMPORAL_START_S = -0.5
TEMPORAL_END_S = 3.5
OPTO_DATASETS = ('HPCLC', 'HPCLCterm')


#%% constants
SAMP_FREQ       = 1250
RUN_ONSET_BIN   = 3750
TRACK_LENGTH_CM = 200
DIST_BIN_CM     = 0.1

def smooth_nan_1d(values, sigma):
    values = np.asarray(values, dtype=float)
    finite = np.isfinite(values).astype(float)
    filled = np.nan_to_num(values, nan=0.0)
    smooth_values = gaussian_filter1d(filled, sigma=sigma, mode='nearest')
    smooth_weights = gaussian_filter1d(finite, sigma=sigma, mode='nearest')
    out = smooth_values / smooth_weights
    out[smooth_weights == 0] = np.nan
    return out

def downsample_1d(values, factor):
    values = np.asarray(values, dtype=float)
    n_bins = values.size // factor
    trimmed = values[:n_bins * factor]
    return np.nanmean(trimmed.reshape(n_bins, factor), axis=1)


#%% speed matching
def speed_matched_early_late_trials(
        beh,
        early,
        late,
        n_bins=7,
        bin_size=500,
        k=1.5,
        min_trials=10
        ):
    speed_times = beh['speed_times_aligned'][1:]
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
        spikes = np.asarray(trains_dist_cell[trial], dtype=float)
        occ = np.asarray(occupancy[trial], dtype=float)
        n = min(spikes.size, occ.size)
        spikes = spikes[:n]
        occ = occ[:n]

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

    min_len = min(len(trial_map) for trial_map in trial_maps)
    trial_maps = np.vstack([trial_map[:min_len] for trial_map in trial_maps])
    sum_spikes = sum_spikes[:min_len]
    sum_occ = sum_occ[:min_len]
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
    bin_samples = int(bin_size_s * samp_freq)
    trial_maps = []

    for trial in trials:
        trace = np.asarray(trains_cell[trial], dtype=float)
        if trace.size < stop:
            continue
        window = trace[start:stop]
        trial_maps.append(downsample_1d(window, bin_samples))

    min_len = min(len(trial_map) for trial_map in trial_maps)
    trial_maps = np.vstack([trial_map[:min_len] for trial_map in trial_maps])
    mean_map = np.nanmean(trial_maps, axis=0)
    return trial_maps, mean_map

def split_half_reliability(trial_maps, block_split=False):
    trial_maps = np.asarray(trial_maps, dtype=float)
    if trial_maps.ndim != 2 or trial_maps.shape[0] < 4:
        return np.nan
    if block_split:
        middle = trial_maps.shape[0] // 2
        # leave one cycle between the two contiguous session blocks
        half_a = np.nanmean(trial_maps[:middle], axis=0)
        half_b = np.nanmean(trial_maps[middle + 1:], axis=0)
    else:
        half_a = np.nanmean(trial_maps[::2], axis=0)
        half_b = np.nanmean(trial_maps[1::2], axis=0)
    finite = np.isfinite(half_a) & np.isfinite(half_b)
    return np.corrcoef(half_a[finite], half_b[finite])[0, 1]

def circular_shift_information_null(
        trial_maps,
        rng,
        occupancy=None,
        n_shuf=200
        ):
    trial_maps = np.asarray(trial_maps, dtype=float)
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
        rng,
        occupancy=None,
        axis_start=0.0,
        n_shuf=200,
        alpha=0.05,
        min_peak_rate=1.0,
        min_split_r=0.2,
        min_width=0.0,
        max_width=np.inf,
        block_split=False
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

    split_r = split_half_reliability(
        trial_maps,
        block_split=block_split,
    )
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

    reference_peaks = reference_peaks[mask]
    condition_peaks = condition_peaks[mask]
    shift = condition_peaks - reference_peaks
    rank_r, _ = spearmanr(reference_peaks, condition_peaks)

    return {
        'n': int(mask.sum()),
        'rank_r': float(rank_r),
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
    corrs = []
    for idx in range(n_bins):
        values_a = matrix_a[:, idx]
        values_b = matrix_b[:, idx]
        finite = np.isfinite(values_a) & np.isfinite(values_b)
        corrs.append(np.corrcoef(values_a[finite], values_b[finite])[0, 1])

    return float(np.nanmean(corrs)) if np.any(np.isfinite(corrs)) else np.nan


#%% metric helpers
def build_cell_condition_metrics(
        cluname,
        row,
        recname,
        dataset,
        trains_cell,
        trains_dist_cell,
        occupancy,
        conditions,
        args,
        rng
        ):
    records = []
    maps = {}

    for condition, trials in conditions.items():
        cycle_condition = (
            dataset in OPTO_DATASETS
            and condition in ['ctrl_pre', 'ctrl_post', 'stim']
        )
        spatial_trial_maps, spatial_map, spatial_occ = make_spatial_rate_maps(
            trains_dist_cell,
            occupancy,
            trials,
        )
        temporal_trial_maps, temporal_map = make_temporal_rate_maps(
            trains_cell,
            trials,
            start_s=TEMPORAL_START_S,
            end_s=TEMPORAL_END_S,
            bin_size_s=TEMPORAL_BIN_SIZE_S,
        )

        spatial = field_summary(
            spatial_trial_maps,
            spatial_map,
            bin_size=SPATIAL_BIN_SIZE_CM,
            occupancy=spatial_occ,
            n_shuf=args.n_shuf,
            rng=rng,
            min_peak_rate=1.0,
            min_split_r=0.2,
            min_width=5.0,
            max_width=120.0,
            block_split=cycle_condition,
        )
        temporal = field_summary(
            temporal_trial_maps,
            temporal_map,
            bin_size=TEMPORAL_BIN_SIZE_S,
            axis_start=TEMPORAL_START_S,
            n_shuf=args.n_shuf,
            rng=rng,
            min_peak_rate=0.5,
            min_split_r=0.2,
            min_width=0.1,
            max_width=2.5,
            block_split=cycle_condition,
        )

        record = {
            'animal': recname.split('-')[0],
            'recname': recname,
            'dataset': dataset,
            'cluname': cluname,
            'profile_source': row['profile_source'],
            'cell_identity': row['cell_identity'],
            'spike_rate': row['spike_rate'],
            'place_cell_matlab': row['place_cell'],
            'run_onset_class': row['class'],
            'run_onset_class_ctrl': row['class_ctrl'],
            'run_onset_class_stim': row['class_stim'],
        }
        record.update({
            'condition': condition,
            'is_stim': int(condition == 'stim'),
            'n_trials': len(trials),
            'n_sandwich_cycles': (
                len(trials)
                if dataset in OPTO_DATASETS
                and condition in ['ctrl_pre', 'ctrl_post', 'stim']
                else np.nan
            ),
            'spatial_info': spatial['info'],
            'spatial_info_p': spatial['info_p'],
            'spatial_split_r': spatial['split_r'],
            'spatial_peak_cm': spatial['peak'],
            'spatial_peak_rate': spatial['peak_rate'],
            'spatial_field_width_cm': spatial['field_width'],
            'is_place_cell': spatial['is_field'],
            'temporal_info': temporal['info'],
            'temporal_info_p': temporal['info_p'],
            'temporal_split_r': temporal['split_r'],
            'temporal_peak_s': temporal['peak'],
            'temporal_peak_rate': temporal['peak_rate'],
            'temporal_field_width_s': temporal['field_width'],
            'is_time_cell': temporal['is_field'],
        })
        records.append(record)

        maps[condition] = {
            'spatial_map': spatial_map,
            'temporal_map': temporal_map,
            'spatial_peak_cm': spatial['peak'],
            'temporal_peak_s': temporal['peak'],
            'is_place_cell': spatial['is_field'],
            'is_time_cell': temporal['is_field'],
        }

    return records, maps

def build_sequence_metrics(recname, dataset, cell_metrics, cell_maps):
    comparisons = []
    payload = {}
    specs = [
        ('early_late', 'baseline', 'early', 'late'),
    ]
    if dataset in OPTO_DATASETS:
        specs.extend([
            ('stim_ctrl_pre', 'baseline', 'ctrl_pre', 'stim'),
            ('stim_ctrl_post', 'baseline', 'ctrl_post', 'stim'),
        ])

    for comparison, ref_condition, condition_a, condition_b in specs:
        if ref_condition not in set(cell_metrics['condition']):
            continue
        for domain in ['spatial', 'temporal']:
            class_col = 'is_place_cell' if domain == 'spatial' else 'is_time_cell'
            reference_cells = cell_metrics[
                (cell_metrics['condition'] == ref_condition)
                & (cell_metrics[class_col])
            ]
            selected = []
            for cluname in reference_cells['cluname']:
                cell_rows = cell_metrics[cell_metrics['cluname'] == cluname]
                conditions = set(cell_rows['condition'])
                if condition_a in conditions and condition_b in conditions:
                    selected.append(cluname)
            if len(selected) < 3:
                continue

            conditions = [ref_condition, condition_a, condition_b]
            peak_key = (
                'spatial_peak_cm' if domain == 'spatial'
                else 'temporal_peak_s'
            )
            ref_peaks, peaks_a, peaks_b = [
                np.asarray([
                    cell_maps[cluname][curr_condition][peak_key]
                    for cluname in selected
                ], dtype=float)
                for curr_condition in conditions
            ]

            map_key = 'spatial_map' if domain == 'spatial' else 'temporal_map'
            matrices = []
            for curr_condition in conditions:
                maps = [
                    np.asarray(
                        cell_maps[cluname][curr_condition][map_key],
                        dtype=float,
                    )
                    for cluname in selected
                ]
                min_len = min(map_.size for map_ in maps)
                matrices.append(np.vstack([
                    map_[:min_len] for map_ in maps
                ]))
            matrix_ref, matrix_a, matrix_b = matrices

            valid_idx = np.where(np.isfinite(ref_peaks))[0]
            order = valid_idx[np.argsort(ref_peaks[valid_idx])]
            ordered_ref = matrix_ref[order]
            ordered_a = matrix_a[order]
            ordered_b = matrix_b[order]
            if domain == 'spatial':
                ordered_ref = np.vstack([
                    downsample_1d(row, SPATIAL_PAYLOAD_DOWNSAMPLE)
                    for row in ordered_ref
                ])
                ordered_a = np.vstack([
                    downsample_1d(row, SPATIAL_PAYLOAD_DOWNSAMPLE)
                    for row in ordered_a
                ])
                ordered_b = np.vstack([
                    downsample_1d(row, SPATIAL_PAYLOAD_DOWNSAMPLE)
                    for row in ordered_b
                ])

            summary_a = sequence_peak_summary(ref_peaks, peaks_a)
            summary_b = sequence_peak_summary(ref_peaks, peaks_b)
            direct = sequence_peak_summary(peaks_a, peaks_b)
            comparisons.append({
                'recname': recname,
                'dataset': dataset,
                'comparison': comparison,
                'domain': domain,
                'ref_condition': ref_condition,
                'condition_a': condition_a,
                'condition_b': condition_b,
                'n_cells': len(selected),
                'condition_a_rank_r': summary_a['rank_r'],
                'condition_b_rank_r': summary_b['rank_r'],
                'condition_a_median_abs_shift': summary_a['median_abs_shift'],
                'condition_b_median_abs_shift': summary_b['median_abs_shift'],
                'direct_rank_r': direct['rank_r'],
                'direct_median_shift': direct['median_shift'],
                'direct_median_abs_shift': direct['median_abs_shift'],
                'condition_a_population_vector_r': (
                    population_vector_correlation(matrix_ref, matrix_a)
                ),
                'condition_b_population_vector_r': (
                    population_vector_correlation(matrix_ref, matrix_b)
                ),
                'population_vector_r': population_vector_correlation(matrix_a, matrix_b),
            })

            payload.setdefault(comparison, {})[domain] = {
                'clunames': [selected[idx] for idx in order],
                'ref_condition': ref_condition,
                'condition_a': condition_a,
                'condition_b': condition_b,
                'ref_matrix': ordered_ref,
                'condition_a_matrix': ordered_a,
                'condition_b_matrix': ordered_b,
            }

    return comparisons, payload


#%% mixed effects
def fit_stim_mixed_models(cell_df):
    if cell_df.empty:
        return pd.DataFrame()

    records = []
    metrics = [
        'spatial_info',
        'spatial_split_r',
        'temporal_info',
        'temporal_split_r',
    ]
    contrast_specs = {
        'stim_vs_pre': np.array([1, 0], dtype=float),
        'stim_vs_post': np.array([1, -1], dtype=float),
        'stim_vs_sandwich': np.array([1, -1/3], dtype=float),
        'post_vs_pre': np.array([0, 1], dtype=float),
    }

    for dataset in OPTO_DATASETS:
        for metric in metrics:
            data = cell_df[
                (cell_df['dataset'] == dataset)
                & (cell_df['condition'].isin([
                    'ctrl_pre', 'ctrl_post', 'stim'
                ]))
            ][['recname', 'cluname', 'condition', metric]].dropna().copy()
            data['cell_id'] = (
                data['recname'] + ':' + data['cluname'].astype(str)
            )
            complete = data.groupby('cell_id')['condition'].nunique()
            data = data[
                data['cell_id'].isin(complete[complete == 3].index)
            ].copy()
            data['is_stim'] = (data['condition'] == 'stim').astype(int)
            data['is_post'] = (data['condition'] == 'ctrl_post').astype(int)

            # recording-level checks cannot estimate a session random effect
            if data['recname'].nunique() < 2:
                continue

            model = smf.mixedlm(
                f'{metric} ~ is_stim + is_post',
                data,
                groups=data['recname'],
                re_formula='1',
                vc_formula={'cell': '0 + C(cell_id)'},
            )
            fit = model.fit(reml=True, method='lbfgs')
            coef_names = ['is_stim', 'is_post']
            coefficients = fit.fe_params[coef_names].to_numpy()
            covariance = fit.cov_params().loc[
                coef_names, coef_names
            ].to_numpy()

            # pre is the reference, so beta_stim - beta_post/3 compares stim
            # with 2/3 pre + 1/3 post after fitting the nonlinear metric
            for contrast, vector in contrast_specs.items():
                estimate = vector @ coefficients
                standard_error = np.sqrt(vector @ covariance @ vector)
                z = estimate / standard_error
                p = 2 * norm.sf(abs(z))
                records.append({
                    'dataset': dataset,
                    'outcome': metric,
                    'contrast': contrast,
                    'n_sessions': data['recname'].nunique(),
                    'n_cells': data['cell_id'].nunique(),
                    'estimate': estimate,
                    'standard_error': standard_error,
                    'ci_low': estimate - 1.96*standard_error,
                    'ci_high': estimate + 1.96*standard_error,
                    'z': z,
                    'p': p,
                    'session_variance': fit.cov_re.iloc[0, 0],
                    'cell_variance': fit.vcomp[0],
                    'converged': fit.converged,
                })

                if contrast == 'stim_vs_sandwich':
                    print(
                        f'{dataset} {metric}: beta={estimate:.4g}, '
                        f'p={p:.4g}, converged={fit.converged}'
                    )

    return pd.DataFrame(records)


#%% session loop
def process_session(recname, dataset, raw_path, session_profiles, args, rng):
    beh_path = pp.behaviour_experiment_stem(dataset) / f'{recname}.pkl'
    with open(beh_path, 'rb') as f:
        beh = pickle.load(f)

    if dataset == 'HPCRaphi':
        train_path = pp.HPC_EPHYS_STEM / 'all_sessions_raphi' / recname / f'{recname}_all_trains_run.npy'
    else:
        train_path = pp.HPC_EPHYS_STEM / 'all_sessions' / recname / f'{recname}_all_trains_run.npy'
    trains = np.load(train_path, allow_pickle=True).item()

    dist_candidates = []
    if dataset != 'HPCRaphi':
        dist_candidates.append(raw_path / (
            f'{recname}_DataStructure_mazeSection1_TrialType1_'
            'convSpikesDistAligned_msess1_Run0.mat'
        ))
    for run in [0, 1]:
        for maze_sess in range(10):
            dist_candidates.append(raw_path / (
                f'{recname}_DataStructure_mazeSection1_TrialType1_'
                f'convSpikesDistAligned_msess{maze_sess}_Run{run}.mat'
            ))
    dist_paths = [
        pp.resolve_matlab_pipeline_file(path, recname)
        for path in dist_candidates
    ]
    dist_path = next((path for path in dist_paths if path.exists()), None)
    if dist_path is None:
        raise FileNotFoundError(
            'missing distance-aligned spike file:\n' +
            '\n'.join(str(path) for path in dist_paths)
        )
    trains_dist = support.load_dist_spike_array(dist_path)
    speeds = support.load_speeds(beh)
    distance_bins = np.arange(
        0,
        TRACK_LENGTH_CM + DIST_BIN_CM,
        DIST_BIN_CM,
    )
    occupancy = [
        support.calculate_occupancy(speed_row, dt=.02, distance_bins=distance_bins)
        for speed_row in speeds
    ]

    n_trials = min(
        len(next(iter(trains.values()))),
        len(occupancy),
    )

    valid_trials = np.asarray([
        beh['bad_trials'][trial + 1] is False
        or beh['bad_trials'][trial + 1] == 0
        for trial in range(n_trials)
    ], dtype=bool)
    stim_conds = [
        beh['trial_statements'][trial + 1][15]
        for trial in range(n_trials)
    ]
    baseline_all, stim_all, _ = get_trialtype_indices_from_stim_conds(stim_conds)
    baseline_all = [int(trial) for trial in baseline_all]
    stim_all = [int(trial) for trial in stim_all]

    if dataset in OPTO_DATASETS:
        baseline = [trial for trial in baseline_all if valid_trials[trial]]
        cycles = [
            (trial, trial - 1, trial + 2)
            for trial in stim_all
            if (
                trial >= 3
                and trial + 2 < n_trials
                and stim_conds[trial - 3] != '0'
                and all(stim_conds[idx] == '0'
                        for idx in [trial - 2, trial - 1, trial + 1, trial + 2])
                and valid_trials[trial]
                and valid_trials[trial - 1]
                and valid_trials[trial + 2]
            )
        ]
        stim = [cycle[0] for cycle in cycles]
        ctrl_pre = [cycle[1] for cycle in cycles]
        ctrl_post = [cycle[2] for cycle in cycles]
    else:
        baseline = np.where(valid_trials)[0].astype(int).tolist()
        cycles = []
        stim = []
        ctrl_pre = []
        ctrl_post = []

    first_lick_s = np.full(n_trials, np.nan, dtype=float)
    for trial in range(n_trials):
        licks = beh['lick_times_aligned'][trial + 1]
        if licks is None or isinstance(licks, float):
            continue
        licks = np.asarray(licks, dtype=float)
        licks = licks[np.isfinite(licks) & (licks > 500)]
        if licks.size:
            first_lick_s[trial] = licks[0] / 1000

    early = []
    late = []
    stim_set = set(stim_all)
    for trial, lick_s in enumerate(first_lick_s):
        if not valid_trials[trial] or not np.isfinite(lick_s):
            continue
        if dataset in OPTO_DATASETS and (
                trial in stim_set or trial - 1 in stim_set
                ):
            continue
        if lick_s < 2.5:
            early.append(trial)
        elif 2.5 < lick_s < 3.5:
            late.append(trial)

    matched_early, matched_late = speed_matched_early_late_trials(
        beh,
        early,
        late,
        min_trials=args.min_trials,
    )
    conditions = {
        'baseline': baseline,
        'ctrl_pre': ctrl_pre,
        'ctrl_post': ctrl_post,
        'stim': stim,
        'early': matched_early,
        'late': matched_late,
    }
    conditions = {
        key: trials
        for key, trials in conditions.items()
        if len(trials) >= args.min_trials
    }

    session_info = {
        'recname': recname,
        'dataset': dataset,
        'train_path': str(train_path),
        'distance_path': str(dist_path),
        'n_trials': n_trials,
        'n_cells': session_profiles.shape[0],
        'n_baseline_trials': len(baseline),
        'n_ctrl_trials': len(ctrl_pre),
        'n_stim_trials': len(stim),
        'n_sandwich_cycles': len(cycles),
        'n_early_trials': len(early),
        'n_late_trials': len(late),
        'n_matched_early_trials': len(matched_early),
        'n_matched_late_trials': len(matched_late),
        'analysed': bool(conditions),
    }
    if not conditions:
        return [], [], {}, session_info

    cluname_to_idx = {
        cluname: idx
        for idx, cluname in enumerate(trains.keys())
    }

    cell_records = []
    cell_maps = {}
    for cluname, row in tqdm(
            session_profiles.iterrows(),
            total=session_profiles.shape[0],
            desc=f'{recname} cells',
            leave=False,
            ):
        clu_idx = cluname_to_idx[cluname]

        records, maps = build_cell_condition_metrics(
            cluname,
            row,
            recname,
            dataset,
            trains[cluname],
            trains_dist[clu_idx],
            occupancy,
            conditions,
            args,
            rng,
        )
        cell_records.extend(records)
        cell_maps[cluname] = maps

    cell_metrics = pd.DataFrame(cell_records)
    sequence_records, sequence_payload = build_sequence_metrics(
        recname,
        dataset,
        cell_metrics,
        cell_maps,
    )

    return cell_records, sequence_records, sequence_payload, session_info


#%% main
def main(argv=None):
    parser = argparse.ArgumentParser(
        description='build CA1 place/time-cell and sequence metrics'
    )
    parser.add_argument('--n-shuf', type=int, default=200,
                        help='number of circular-shift shuffles for each field call')
    parser.add_argument('--min-trials', type=int, default=10,
                        help='minimum trials per condition')
    parser.add_argument('--recording', dest='recording_filter',
                        help='optional recname substring filter')
    parser.add_argument('--max-recordings', type=int,
                        help='optional cap for quick checks')
    parser.add_argument('--no-raphi', action='store_true',
                        help='skip Raphi/passive HPC recordings')
    parser.add_argument('--seed', type=int, default=6102026,
                        help='random seed for shuffle nulls')
    parser.add_argument('--payload', type=Path, default=DEFAULT_PAYLOAD_PATH,
                        help='pickle output path')
    parser.add_argument('--cell-csv', type=Path, default=DEFAULT_CELL_CSV_PATH,
                        help='cell-level CSV output path')
    parser.add_argument('--session-csv', type=Path, default=DEFAULT_SESSION_CSV_PATH,
                        help='session-level CSV output path')
    parser.add_argument('--mixed-csv', type=Path, default=DEFAULT_MIXED_CSV_PATH,
                        help='mixed-effects summary CSV output path')
    args = parser.parse_args(argv)
    rng = np.random.default_rng(args.seed)

    bad_behaviour_sessions = {
        Path(path).name for path in rec_list.pathHPCbadbeh
    }
    soma_paths = [
        path for path in rec_list.pathHPCLCopt
        if Path(path).name not in bad_behaviour_sessions
    ]
    lookup = {}
    for dataset, paths in [
            ('HPCLC', soma_paths),
            ('HPCLCterm', rec_list.pathHPCLCtermopt),
            ('HPCRaphi', rec_list.pathHPC_Raphi),
            ]:
        for path in paths:
            lookup[Path(path).name] = {
                'dataset': dataset,
                'path': Path(path),
            }

    profile_path = pp.HPC_EPHYS_STEM / 'hpc_all_profiles.pkl'
    profiles = pd.read_pickle(profile_path).copy()
    profiles['profile_source'] = profile_path.name
    if not args.no_raphi:
        raphi_path = pp.HPC_EPHYS_STEM / 'hpc_all_profiles_raphi.pkl'
        raphi_profiles = pd.read_pickle(raphi_path).copy()
        raphi_profiles['profile_source'] = raphi_path.name
        raphi_profiles['spike_rate'] = raphi_profiles['firing_rate']
        profiles = pd.concat([profiles, raphi_profiles], axis=0, sort=False)
    profiles = profiles[profiles['cell_identity'] == 'pyr'].copy()
    profiles = profiles[profiles['recname'].isin(lookup)].copy()

    recnames = []
    for recname in sorted(pd.unique(profiles['recname']).tolist()):
        if args.recording_filter and args.recording_filter not in recname:
            continue
        recnames.append(recname)
        if args.max_recordings and len(recnames) >= args.max_recordings:
            break

    all_cell_records = []
    all_sequence_records = []
    session_payload = {}
    session_infos = []

    for recname in recnames:
        dataset = lookup[recname]['dataset']
        raw_path = lookup[recname]['path']
        session_profiles = profiles[profiles['recname'] == recname]
        print_session(recname)
        (
            cell_records,
            sequence_records,
            sequence_payload,
            session_info,
        ) = process_session(
            recname,
            dataset,
            raw_path,
            session_profiles,
            args,
            rng,
        )

        all_cell_records.extend(cell_records)
        all_sequence_records.extend(sequence_records)
        session_payload[recname] = sequence_payload
        session_infos.append(session_info)
        print_status(
            'done' if session_info['analysed'] else 'not enough valid trials'
        )

    cell_df = pd.DataFrame(all_cell_records)
    session_df = pd.DataFrame(all_sequence_records)
    session_info_df = pd.DataFrame(session_infos)
    mixed_df = fit_stim_mixed_models(cell_df)

    payload = {
        'metadata': {
            'n_shuf': args.n_shuf,
            'min_trials': args.min_trials,
            'seed': args.seed,
            'temporal_window_s': (TEMPORAL_START_S, TEMPORAL_END_S),
            'temporal_bin_size_s': TEMPORAL_BIN_SIZE_S,
            'spatial_bin_size_cm': SPATIAL_BIN_SIZE_CM,
            'spatial_payload_bin_size_cm': (
                SPATIAL_BIN_SIZE_CM * SPATIAL_PAYLOAD_DOWNSAMPLE
            ),
            'control_definition': (
                'metric-scale 2/3 stim-1 + 1/3 stim+2'
            ),
            'excluded_bad_behaviour_sessions': sorted(bad_behaviour_sessions),
        },
        'cell_metrics': cell_df,
        'session_metrics': session_df,
        'session_info': session_info_df,
        'mixed_effects': mixed_df,
        'session_payload': session_payload,
    }

    args.payload.parent.mkdir(parents=True, exist_ok=True)
    with open(args.payload, 'wb') as f:
        pickle.dump(payload, f, protocol=pickle.HIGHEST_PROTOCOL)
    cell_df.to_csv(args.cell_csv, index=False)
    session_df.to_csv(args.session_csv, index=False)
    session_info_df.to_csv(
        args.session_csv.with_name('hpc_sequence_session_info.csv'),
        index=False,
    )
    mixed_df.to_csv(args.mixed_csv, index=False)

    print_files_saved([
        ('analysis data', args.payload),
        ('cell summary', args.cell_csv),
        ('session summary', args.session_csv),
        ('session info', args.session_csv.with_name('hpc_sequence_session_info.csv')),
        ('mixed-effects summary', args.mixed_csv),
    ])

if __name__ == '__main__':
    main()
