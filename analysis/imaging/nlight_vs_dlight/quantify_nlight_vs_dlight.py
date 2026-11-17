# -*- coding: utf-8 -*-
'''
Created on Sat Aug 29 02:24:22 2026

quantify natural and optogenetic nlight responses against chapter 2 dlight analyses

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

from argparse import ArgumentParser
import warnings

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from common_functions import mpl_formatting
mpl_formatting()
import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d
from scipy.optimize import curve_fit
from scipy.stats import mannwhitneyu, spearmanr, wilcoxon

# Some stored arrays were written by NumPy 2, while this analysis environment
# uses NumPy 1.  The alias changes pickle lookup only; it does not alter arrays.
sys.modules.setdefault('numpy._core', np.core)
sys.modules.setdefault('numpy._core.multiarray', np.core.multiarray)


#%% paths and parameters
repo_root = Path(__file__).resolve().parents[3]
JINGYU_CODE = pp.JINGYU_CODE_ROOT
NLIGHT_ROOT = pp.NLIGHT_RAW_ROOT
DLIGHT_NATURAL_ROOT = pp.DBH_DLIGHT_RAW_ROOT

NATURAL_DLIGHT_TABLE = (
    DLIGHT_NATURAL_ROOT / 'processed_dataframe_grid_free_dilation'
    / 'df_population_profile_pooled_dilation=0_pre(-1, 0)_post(0, 1)_ES=0.05_shuff95_0621.parquet'
)
NATURAL_NLIGHT_TABLE_ROOT = NLIGHT_ROOT / 'processed_dataframe_grid_free_dilation'
NATURAL_DLIGHT_DILATION_ROOT = DLIGHT_NATURAL_ROOT / 'processed_dataframe_grid_free_dilation'
NATURAL_DLIGHT_REGRESSION_ROOT = DLIGHT_NATURAL_ROOT / 'regression_res_grid_free_dilation'
NATURAL_NLIGHT_REGRESSION_ROOT = NLIGHT_ROOT / 'regression_res_grid_free_dilation'
NATURAL_DLIGHT_BEHAVIOUR_ROOT = DLIGHT_NATURAL_ROOT / 'behaviour_profile'
NATURAL_NLIGHT_BEHAVIOUR_ROOT = (
    repo_root / 'data' / 'behaviour' / 'all_experiments' / 'hpc_nlight_lc_opto'
)

OPTO_TABLE_ROOT = repo_root / 'data' / 'catecholamine_comparison'
OUTPUT_ROOT = OPTO_TABLE_ROOT / 'nlight_vs_dlight_single_state'
FIGURE_ROOT = repo_root / 'figures' / 'nlight_vs_dlight_single_state'

IMAGING_RATE             = 30
PROFILE_BEFORE_S         = 2
PROFILE_AFTER_S          = 4
TEMPORAL_SMOOTH_SIGMA    = 5
TEMPORAL_R2_MIN          = 0.7
TEMPORAL_PEAK_RANGE_S    = (0, 3)
SPATIAL_R2_MIN           = 0.3
PRIMARY_SHUFFLE_PERCENT  = 95
PRIMARY_EFFECT_SIZE      = 0.05
SHUFFLE_SENSITIVITY      = (90, 95, 99)
EFFECT_SIZE_SENSITIVITY  = (0, 0.05, 0.1)

SENSOR_COLOURS = {'dLight': '#4f8f62', 'nLight': '#8b5ea7'}


#%% small helpers
def require_inputs(paths):
    missing = [str(path) for path in paths if not path.exists()]
    if missing:
        raise FileNotFoundError('Missing required inputs:\n' + '\n'.join(missing))


def profile_is_valid(profile):
    if profile is None:
        return False
    values = np.asarray(profile)
    return values.size > 0 and bool(np.isfinite(values).any())


def finite_summary(values):
    values = pd.Series(values, dtype=float).replace([np.inf, -np.inf], np.nan).dropna()
    if values.empty:
        return {'n': 0, 'median': np.nan, 'q1': np.nan, 'q3': np.nan}
    return {
        'n': int(len(values)),
        'median': float(values.median()),
        'q1': float(values.quantile(0.25)),
        'q3': float(values.quantile(0.75)),
    }


def format_number(value, digits=3):
    if not np.isfinite(value):
        return 'unavailable'
    return f'{value:.{digits}f}'


def format_p(value):
    if not np.isfinite(value):
        return 'unavailable'
    if value < 0.001:
        return f'{value:.2e}'
    return f'{value:.3f}'






def fit_temporal_profile(profile):
    y = gaussian_filter1d(np.asarray(profile, dtype=float) * 100, sigma=TEMPORAL_SMOOTH_SIGMA)
    if y.size < 5 or not np.isfinite(y).all():
        return None
    peak_index = int(np.nanargmax(y))
    t = np.arange(y.size, dtype=float)
    t_fit = t[peak_index:] - t[peak_index]
    y_fit = y[peak_index:]
    if y_fit.size < 5:
        return None
    baseline = float(np.nanmin(y_fit))

    try:
        fitted, _ = curve_fit(
            lambda x, amplitude, tau: amplitude * np.exp(-(x - x[0]) / tau) + baseline,
            t_fit,
            y_fit,
            maxfev=10000,
        )
    except (RuntimeError, ValueError, FloatingPointError):
        return None

    amplitude, tau_frames = fitted
    prediction = amplitude * np.exp(-(t_fit - t_fit[0]) / tau_frames) + baseline
    ss_res = float(np.sum((y_fit - prediction) ** 2))
    ss_tot = float(np.sum((y_fit - np.mean(y_fit)) ** 2))
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
    peak_s = peak_index / IMAGING_RATE - PROFILE_BEFORE_S
    tau_s = tau_frames / IMAGING_RATE
    included = (
        np.isfinite(tau_s)
        and tau_s > 0
        and np.isfinite(r2)
        and r2 > TEMPORAL_R2_MIN
        and TEMPORAL_PEAK_RANGE_S[0] < peak_s < TEMPORAL_PEAK_RANGE_S[1]
    )
    return {
        'peak_index': peak_index,
        'peak_s': float(peak_s),
        'tau_s': float(tau_s),
        'r2': float(r2),
        'amplitude': float(amplitude),
        'baseline': baseline,
        'profile_percent': y,
        'fit_t_s': t[peak_index:] / IMAGING_RATE - PROFILE_BEFORE_S,
        'fit_y_percent': prediction,
        'included': bool(included),
    }


def fit_spatial_profile(distance, response, r2_min=SPATIAL_R2_MIN):
    distance = np.asarray(distance, dtype=float)
    response = np.asarray(response, dtype=float)
    valid = np.isfinite(distance) & np.isfinite(response)
    distance = distance[valid]
    response = response[valid]
    if response.size < 6:
        return None
    baseline = float(np.min(response))
    try:
        fitted, _ = curve_fit(
            lambda x, amplitude, tau: amplitude * np.exp(-(x - x[0]) / tau) + baseline,
            distance,
            response,
            maxfev=10000,
        )
    except (RuntimeError, ValueError, FloatingPointError):
        return None
    amplitude, tau_px = fitted
    prediction = amplitude * np.exp(-(distance - distance[0]) / tau_px) + baseline
    ss_res = float(np.sum((response - prediction) ** 2))
    ss_tot = float(np.sum((response - np.mean(response)) ** 2))
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else np.nan
    return {
        'tau_px': float(tau_px),
        'r2': float(r2),
        'amplitude': float(amplitude),
        'baseline': baseline,
        'distance_px': distance,
        'response': response,
        'fit_response': prediction,
        'included': bool(np.isfinite(tau_px) and tau_px > 0 and r2 > r2_min),
    }


#%% natural-run tables and fits
NATURAL_COLUMNS = [
    'roi_id', 'n_valid_event', 'n_keep_trial', 'response_amplitude', 'effect_size',
    'shuff_response_amplitude', 'mean_profile', 'edge', 'rec_id', 'anm',
    'dlight_valid', 'red_valid', 'Up',
]


def selected_nlight_sessions():
    pattern = '*_profile_dilation=0_stat_dlight_pre(-1, 0)_post(0, 1).parquet'
    return sorted(path.name.split('_profile_dilation')[0] for path in NATURAL_NLIGHT_TABLE_ROOT.glob(pattern))


def rebuild_current_nlight_table(recording):
    dlight_path = NATURAL_NLIGHT_TABLE_ROOT / (
        f'{recording}_profile_dilation=0_stat_dlight_pre(-1, 0)_post(0, 1).parquet'
    )
    red_path = NATURAL_NLIGHT_TABLE_ROOT / (
        f'{recording}_profile_dilation=0_stat_red_pre(-1, 0)_post(0, 1).parquet'
    )
    regression_dir = (
        NATURAL_NLIGHT_REGRESSION_ROOT / recording / 'single_trial_regression'
        / 'dilation_k=0'
    )
    masks_dir = NATURAL_NLIGHT_REGRESSION_ROOT / recording / 'masks'
    require_inputs([
        dlight_path,
        red_path,
        regression_dir / 'baseline_corrected_dlight.npy',
        regression_dir / 'baseline_red.npy',
        masks_dir / 'dilated_global_axon_k=0.npy',
        masks_dir / 'global_dlight_mask_enhanced.npy',
    ])

    stats = pd.read_parquet(dlight_path).copy()
    red = pd.read_parquet(red_path)
    if len(stats) != len(red):
        raise ValueError(f'nLight/red row mismatch for {recording}')
    stats['mean_profile_red'] = red['mean_profile'].to_numpy()
    stats['rec_id'] = recording
    stats['anm'] = recording.split('-')[0]

    baseline_dlight = np.load(regression_dir / 'baseline_corrected_dlight.npy')
    baseline_red = np.load(regression_dir / 'baseline_red.npy')
    coords = np.asarray(stats['roi_id'].tolist(), dtype=int)
    ys, xs = coords[:, 0], coords[:, 1]
    stats['baseline_dlight_min'] = np.nanmin(baseline_dlight[ys, xs], axis=-1)
    stats['baseline_red_min'] = np.nanmin(baseline_red[ys, xs], axis=-1)

    axon = np.load(masks_dir / 'dilated_global_axon_k=0.npy')
    dlight = np.load(masks_dir / 'global_dlight_mask_enhanced.npy')
    grid_rows, grid_cols = baseline_dlight.shape[:2]
    if axon.ndim == 2:
        image_rows, image_cols = axon.shape
        tile_rows, tile_cols = image_rows // grid_rows, image_cols // grid_cols
        reshape = (grid_rows, tile_rows, grid_cols, tile_cols)
        stats['n_pixels_axon'] = axon.reshape(reshape).sum((1, 3))[ys, xs]
        stats['n_pixels_dlight'] = dlight.reshape(reshape).sum((1, 3))[ys, xs]
        stats['n_pixels_axon_and_dlight'] = (
            (axon & dlight).reshape(reshape).sum((1, 3))[ys, xs]
        )
    elif axon.ndim == 4:
        stats['n_pixels_axon'] = axon[ys, xs].sum((1, 2))
        stats['n_pixels_axon_and_dlight'] = (axon[ys, xs] & dlight).sum((1, 2))
        tile_rows = dlight.shape[0] // grid_rows
        tile_cols = dlight.shape[1] // grid_cols
        dlight_counts = dlight.reshape(
            grid_rows, tile_rows, grid_cols, tile_cols
        ).sum((1, 3))
        stats['n_pixels_dlight'] = dlight_counts[ys, xs]
    else:
        raise ValueError(f'Unexpected axon mask shape for {recording}: {axon.shape}')
    stats['edge'] = (
        (ys == 0) | (ys == grid_rows - 1)
        | (xs == 0) | (xs == grid_cols - 1)
    )
    stats['dlight_valid'] = (
        (stats['baseline_dlight_min'] > 2)
        & stats['mean_profile'].apply(profile_is_valid)
    )
    stats['red_valid'] = (
        (stats['baseline_red_min'] > 1)
        & stats['mean_profile_red'].apply(profile_is_valid)
    )
    stats['shuffle_amps_thresh_up'] = stats['shuff_response_amplitude'].apply(
        lambda values: np.nanpercentile(values, PRIMARY_SHUFFLE_PERCENT)
    )
    stats['Up'] = (
        ~stats['edge']
        & (stats['response_amplitude'] > stats['shuffle_amps_thresh_up'])
        & (stats['effect_size'] > PRIMARY_EFFECT_SIZE)
    )
    return stats


def load_natural_tables():
    require_inputs([NATURAL_DLIGHT_TABLE, NATURAL_NLIGHT_TABLE_ROOT])
    dlight = pd.read_parquet(NATURAL_DLIGHT_TABLE, columns=NATURAL_COLUMNS)
    nlight = pd.concat(
        [rebuild_current_nlight_table(recording) for recording in selected_nlight_sessions()],
        ignore_index=True,
    )[NATURAL_COLUMNS]
    return {'dLight': dlight, 'nLight': nlight}


def natural_metrics(tables):
    session_rows = []
    roi_fit_rows = []
    session_profiles = {'dLight': {}, 'nLight': {}}

    for sensor, table in tables.items():
        valid = table['dlight_valid'] & table['red_valid'] & ~table['edge']
        table = table.copy()
        table['valid_primary'] = valid
        table['up_primary'] = valid & table['Up']

        for row in table.loc[table['up_primary']].itertuples(index=False):
            fit = fit_temporal_profile(row.mean_profile)
            roi_y, roi_x = row.roi_id
            roi_fit_rows.append({
                'sensor': sensor,
                'event': 'natural_run_onset',
                'session': row.rec_id,
                'animal': row.anm,
                'roi_y': int(roi_y),
                'roi_x': int(roi_x),
                'temporal_peak_s': fit['peak_s'] if fit else np.nan,
                'temporal_tau_s': fit['tau_s'] if fit else np.nan,
                'temporal_r2': fit['r2'] if fit else np.nan,
                'temporal_included': bool(fit and fit['included']),
            })

        fits = pd.DataFrame(roi_fit_rows)
        fits = fits[(fits['sensor'] == sensor) & fits['temporal_included']]

        for recording, group in table.groupby('rec_id', sort=True):
            animal = str(group['anm'].iloc[0])
            group_valid = group[group['valid_primary']]
            group_up = group[group['up_primary']]
            group_fits = fits[fits['session'] == recording]
            if len(group_up):
                session_profiles[sensor][recording] = np.nanmean(
                    np.stack(group_up['mean_profile'].to_numpy()), axis=0
                )
            session_rows.append({
                'sensor': sensor,
                'transmitter': 'dopamine' if sensor == 'dLight' else 'norepinephrine',
                'condition': 'untreated',
                'event': 'natural_run_onset',
                'animal': animal,
                'session': recording,
                'included': True,
                'exclusion_reason': '',
                'n_valid_roi': int(len(group_valid)),
                'n_up_roi': int(len(group_up)),
                'up_fraction': len(group_up) / len(group_valid) if len(group_valid) else np.nan,
                'response_amplitude_all_dff': float(group_valid['response_amplitude'].median()) if len(group_valid) else np.nan,
                'response_amplitude_up_dff': float(group_up['response_amplitude'].median()) if len(group_up) else np.nan,
                'n_temporal_fit_roi': int(len(group_fits)),
                'temporal_tau_s': float(group_fits['temporal_tau_s'].median()) if len(group_fits) else np.nan,
                'rise_time_s': float(group_fits['temporal_peak_s'].median()) if len(group_fits) else np.nan,
                'spatial_tau_px': np.nan,
                'n_spatial_fit_roi': 0,
                'first_lick_rho': np.nan,
                'first_lick_p': np.nan,
                'n_first_lick_trials': 0,
            })

    return pd.DataFrame(session_rows), pd.DataFrame(roi_fit_rows), session_profiles


#%% natural spatial fits, threshold sensitivity and lick timing
def natural_dlight_spatial(session_metrics, roi_metrics):
    recordings = sorted(
        session_metrics.loc[
            (session_metrics['sensor'] == 'dLight')
            & (session_metrics['event'] == 'natural_run_onset'),
            'session',
        ]
    )
    columns = [
        'roi_id', 'dilation_k', 'response_amplitude', 'rec_id',
        'dlight_valid', 'red_valid', 'edge', 'Up',
    ]
    frames = []
    for recording in recordings:
        for dilation in range(11):
            path = NATURAL_DLIGHT_DILATION_ROOT / (
                f'{recording}_profile_combined_dilation={dilation}_pre(-1, 0)_post(0, 1)_'
                'ES=0.05_shuff95.parquet'
            )
            require_inputs([path])
            frames.append(pd.read_parquet(path, columns=columns))

    dilation = pd.concat(frames, ignore_index=True)
    dilation['roi_tuple'] = dilation['roi_id'].apply(tuple)
    valid = dilation['dlight_valid'] & ~dilation['edge']
    selected = dilation.loc[
        (dilation['dilation_k'] == 0)
        & valid
        & dilation['red_valid']
        & dilation['Up'],
        ['rec_id', 'roi_tuple'],
    ].drop_duplicates()
    fit_source = dilation.loc[valid].merge(selected, on=['rec_id', 'roi_tuple'])

    rows = []
    curve_rows = []
    for (recording, roi), group in fit_source.groupby(['rec_id', 'roi_tuple']):
        group = group.sort_values('dilation_k')
        fit = fit_spatial_profile(group['dilation_k'], group['response_amplitude'])
        rows.append({
            'sensor': 'dLight',
            'event': 'natural_run_onset',
            'session': recording,
            'animal': recording.split('-')[0],
            'roi_y': int(roi[0]),
            'roi_x': int(roi[1]),
            'spatial_tau_px': fit['tau_px'] if fit else np.nan,
            'spatial_r2': fit['r2'] if fit else np.nan,
            'spatial_included': bool(fit and fit['included']),
        })
        if fit and fit['included']:
            for distance, response in zip(fit['distance_px'], fit['response']):
                curve_rows.append({
                    'sensor': 'dLight',
                    'event': 'natural_run_onset',
                    'session': recording,
                    'roi_y': int(roi[0]),
                    'roi_x': int(roi[1]),
                    'distance_px': float(distance),
                    'response': float(response),
                })

    spatial_roi = pd.DataFrame(rows)
    included = spatial_roi[spatial_roi['spatial_included']]
    for recording, group in included.groupby('session'):
        mask = (
            session_metrics['sensor'].eq('dLight')
            & session_metrics['event'].eq('natural_run_onset')
            & session_metrics['session'].eq(recording)
        )
        session_metrics.loc[mask, 'spatial_tau_px'] = group['spatial_tau_px'].median()
        session_metrics.loc[mask, 'n_spatial_fit_roi'] = len(group)
        session_metrics.loc[mask, 'spatial_exclusion_reason'] = ''

    no_fit = (
        session_metrics['sensor'].eq('dLight')
        & session_metrics['event'].eq('natural_run_onset')
        & session_metrics['spatial_tau_px'].isna()
    )
    session_metrics.loc[no_fit, 'spatial_exclusion_reason'] = 'no ACR passed spatial R2 and positive-tau QC'
    nlight = session_metrics['sensor'].eq('nLight') & session_metrics['event'].eq('natural_run_onset')
    session_metrics.loc[nlight, 'spatial_exclusion_reason'] = (
        'unavailable: saved nLight natural outputs contain dilation k=0 only'
    )

    roi_metrics = roi_metrics.merge(
        spatial_roi,
        on=['sensor', 'event', 'session', 'animal', 'roi_y', 'roi_x'],
        how='outer',
    )
    return session_metrics, roi_metrics, pd.DataFrame(curve_rows)


def threshold_sensitivity(tables):
    rows = []
    for sensor, table in tables.items():
        valid = table['dlight_valid'] & table['red_valid'] & ~table['edge']
        source = table.loc[valid].copy()
        shuffled = source['shuff_response_amplitude'].to_numpy()
        for percentile in SHUFFLE_SENSITIVITY:
            shuffle_threshold = np.array([
                np.nanpercentile(np.asarray(values, dtype=float), percentile)
                for values in shuffled
            ])
            for effect_size in EFFECT_SIZE_SENSITIVITY:
                classified = (
                    source['response_amplitude'].to_numpy() > shuffle_threshold
                ) & (source['effect_size'].to_numpy() > effect_size)
                work = pd.DataFrame({
                    'session': source['rec_id'].to_numpy(),
                    'classified': classified,
                })
                fractions = work.groupby('session')['classified'].mean()
                summary = finite_summary(fractions)
                rows.append({
                    'sensor': sensor,
                    'shuffle_percentile': percentile,
                    'effect_size_threshold': effect_size,
                    'n_sessions': summary['n'],
                    'median_up_fraction': summary['median'],
                    'q1_up_fraction': summary['q1'],
                    'q3_up_fraction': summary['q3'],
                    'pooled_up_fraction': float(np.mean(classified)),
                })
    return pd.DataFrame(rows)


def response_trial_indices(behaviour):
    onsets = np.asarray(behaviour['run_onset_frames'], dtype=float)
    selected = np.isfinite(onsets) & (onsets >= 0)
    for key in ['reward_times', 'non_stop_trials', 'non_fullstop_trials']:
        if key not in behaviour:
            return np.array([], dtype=int), np.array([], dtype=int)
    selected &= np.isfinite(np.asarray(behaviour['reward_times'], dtype=float))
    selected &= np.asarray(behaviour['non_stop_trials']) == 0
    selected &= np.asarray(behaviour['non_fullstop_trials']) == 0
    if 'trial_statements' in behaviour:
        pulse_mode = np.asarray([statement[15] for statement in behaviour['trial_statements']], dtype=float)
        selected &= np.isfinite(pulse_mode) & (pulse_mode == 0)
    indices = np.where(selected)[0]
    return indices, onsets[selected].astype(int)


def add_first_lick_correlations(session_metrics, natural_tables):
    for sensor, table in natural_tables.items():
        valid_up = table['dlight_valid'] & table['red_valid'] & ~table['edge'] & table['Up']
        for recording, group in table.loc[valid_up].groupby('rec_id'):
            if sensor == 'dLight':
                behaviour_path = NATURAL_DLIGHT_BEHAVIOUR_ROOT / f'{recording}.pkl'
                regression_root = NATURAL_DLIGHT_REGRESSION_ROOT
            else:
                behaviour_path = NATURAL_NLIGHT_BEHAVIOUR_ROOT / f'{recording}.pkl'
                regression_root = NATURAL_NLIGHT_REGRESSION_ROOT
            dff_path = (
                regression_root / recording / 'single_trial_regression' / 'dilation_k=0'
                / 'dff_corrected_dlight.npy'
            )
            if not behaviour_path.exists() or not dff_path.exists():
                continue
            behaviour = pd.read_pickle(behaviour_path)
            trial_indices, onsets = response_trial_indices(behaviour)
            if len(trial_indices) != len(onsets):
                continue
            dff = np.load(dff_path, mmap_mode='r')
            coordinates = np.asarray(group['roi_id'].tolist(), dtype=int)
            ys, xs = coordinates[:, 0], coordinates[:, 1]
            amplitudes = []
            first_licks = []
            for trial_index, onset in zip(trial_indices, onsets):
                if onset < IMAGING_RATE or onset + IMAGING_RATE > dff.shape[-1]:
                    continue
                licks = np.asarray(behaviour['lick_times_aligned'][trial_index], dtype=float)
                licks = licks[np.isfinite(licks) & (licks > 500)]
                if not len(licks):
                    continue
                traces = np.asarray(dff[ys, xs, onset - IMAGING_RATE:onset + IMAGING_RATE], dtype=float)
                complete = np.isfinite(traces).all(axis=1)
                if not complete.any():
                    continue
                traces = traces[complete]
                amplitude = np.mean(traces[:, IMAGING_RATE:]) - np.mean(traces[:, :IMAGING_RATE])
                if np.isfinite(amplitude):
                    amplitudes.append(float(amplitude))
                    first_licks.append(float(licks[0] / 1000))
            if len(amplitudes) < 8 or len(np.unique(first_licks)) < 2:
                continue
            rho, p_value = spearmanr(first_licks, amplitudes)
            mask = (
                session_metrics['sensor'].eq(sensor)
                & session_metrics['event'].eq('natural_run_onset')
                & session_metrics['session'].eq(recording)
            )
            session_metrics.loc[mask, 'first_lick_rho'] = rho
            session_metrics.loc[mask, 'first_lick_p'] = p_value
            session_metrics.loc[mask, 'n_first_lick_trials'] = len(amplitudes)
    return session_metrics


#%% optogenetic tables
def load_opto_metrics(selected_sessions):
    time = pd.read_csv(OPTO_TABLE_ROOT / 'roi_neuropil_timecourse_metrics.csv')
    spatial = pd.read_csv(OPTO_TABLE_ROOT / 'spatial_tau_metrics.csv')
    release = pd.read_csv(OPTO_TABLE_ROOT / 'roi_release_metrics.csv')
    curves = pd.read_csv(OPTO_TABLE_ROOT / 'spatial_dilation_curves.csv')
    spatial = spatial[spatial['window'] == '0-1 s'].copy()
    curves = curves[curves['window'] == '0-1 s'].copy()

    allowed = {'dLight': None, 'nLight': set(selected_sessions)}
    rows = []
    for sensor in ['dLight', 'nLight']:
        base = release[release['sensor'] == sensor].copy().set_index('recname')
        sensor_time_frame = time[time['sensor'] == sensor].copy()
        sensor_spatial_frame = spatial[spatial['sensor'] == sensor].copy()
        if allowed[sensor] is not None:
            base = base[base.index.isin(allowed[sensor])]
            sensor_time_frame = sensor_time_frame[sensor_time_frame['recname'].isin(allowed[sensor])]
            sensor_spatial_frame = sensor_spatial_frame[sensor_spatial_frame['recname'].isin(allowed[sensor])]
        sensor_time = sensor_time_frame.set_index('recname')
        sensor_spatial = sensor_spatial_frame.set_index('recname')
        recordings = sorted(set(base.index) | set(sensor_time.index) | set(sensor_spatial.index))
        for recording in recordings:
            release_row = base.loc[recording] if recording in base.index else None
            time_row = sensor_time.loc[recording] if recording in sensor_time.index else None
            spatial_row = sensor_spatial.loc[recording] if recording in sensor_spatial.index else None
            tau_s = float(time_row['roi_tau_s']) if time_row is not None else np.nan
            tau_px = float(spatial_row['spatial_tau_px']) if spatial_row is not None else np.nan
            if release_row is not None:
                animal = release_row['animal']
                n_valid_roi = int(release_row['n_roi'])
                n_up_roi = int(round(release_row['prop_releasing'] * release_row['n_roi']))
                up_fraction = float(release_row['prop_releasing'])
                response_amplitude = float(release_row['mean_RI'])
            else:
                source_row = time_row if time_row is not None else spatial_row
                animal = source_row['animal']
                n_valid_roi = int(source_row['n_roi'] if time_row is not None else source_row['n_roi_total'])
                n_up_roi = 0
                up_fraction = np.nan
                response_amplitude = np.nan
            rows.append({
                'sensor': sensor,
                'transmitter': 'dopamine' if sensor == 'dLight' else 'norepinephrine',
                'condition': 'untreated',
                'event': 'opto_lc_stimulation',
                'animal': animal,
                'session': recording,
                'included': True,
                'exclusion_reason': '',
                'n_valid_roi': n_valid_roi,
                'n_up_roi': n_up_roi,
                'up_fraction': up_fraction,
                'response_amplitude_all_dff': response_amplitude,
                'response_amplitude_up_dff': np.nan,
                'n_temporal_fit_roi': int(time_row['n_roi']) if time_row is not None and np.isfinite(tau_s) else 0,
                'temporal_tau_s': tau_s,
                'rise_time_s': np.nan,
                'spatial_tau_px': tau_px,
                'n_spatial_fit_roi': int(spatial_row['n_roi_releasing']) if spatial_row is not None and np.isfinite(tau_px) else 0,
                'first_lick_rho': np.nan,
                'first_lick_p': np.nan,
                'n_first_lick_trials': 0,
                'spatial_exclusion_reason': '' if np.isfinite(tau_px) else 'no session spatial fit passed Chapter 2 QC',
                'temporal_exclusion_reason': '' if np.isfinite(tau_s) else 'no session temporal fit passed Chapter 2 QC',
            })

    filtered_curves = []
    for sensor in ['dLight', 'nLight']:
        curr = curves[curves['sensor'] == sensor].copy()
        if allowed[sensor] is not None:
            curr = curr[curr['recname'].isin(allowed[sensor])]
        filtered_curves.append(curr)
    return pd.DataFrame(rows), pd.concat(filtered_curves, ignore_index=True)


#%% optogenetic recruitment sensitivity and formal comparisons
def add_wilcoxon_release(session_metrics):
    np_core = np.core
    sys.modules.setdefault('numpy._core', np_core)
    sys.modules.setdefault('numpy._core.multiarray', np_core.multiarray)
    roots = {
        'dLight': repo_root / 'data' / 'hpc_dlight' / 'lc_opto',
        'nLight': repo_root / 'data' / 'hpc_nlight' / 'lc_opto',
    }
    opto = session_metrics[session_metrics['event'] == 'opto_lc_stimulation']
    for row in opto.itertuples(index=False):
        processed = roots[row.sensor] / 'all_sessions' / row.session / 'processed_data'
        response_path = processed / f'{row.session}_pixel_RI_stim.npy'
        roi_path = processed / f'{row.session}_ROI_dict.npy'
        if not response_path.exists() or not roi_path.exists():
            continue
        response = np.load(response_path, mmap_mode='r')
        roi_dict = np.load(roi_path, allow_pickle=True).item()
        n_releasing = 0
        for roi in roi_dict.values():
            means = np.nanmean(response[roi['ypix'], roi['xpix'], :], axis=0)
            means = means[np.isfinite(means)]
            if len(means) <= 2:
                continue
            try:
                p_value = wilcoxon(means, alternative='greater').pvalue
            except ValueError:
                p_value = 1
            if p_value < 0.05 and np.mean(means) > 0.1:
                n_releasing += 1
        fraction = n_releasing / len(roi_dict) if len(roi_dict) else np.nan
        mask = (
            session_metrics['event'].eq('opto_lc_stimulation')
            & session_metrics['sensor'].eq(row.sensor)
            & session_metrics['session'].eq(row.session)
        )
        session_metrics.loc[mask, 'up_fraction_wilcoxon'] = fraction
        session_metrics.loc[mask, 'n_up_roi_wilcoxon'] = n_releasing
    return session_metrics


def mann_comparison(session_metrics, event, metric, label):
    current = session_metrics[session_metrics['event'] == event]
    nlight = current.loc[current['sensor'] == 'nLight', metric].dropna().astype(float)
    dlight = current.loc[current['sensor'] == 'dLight', metric].dropna().astype(float)
    if len(nlight) and len(dlight):
        result = mannwhitneyu(nlight, dlight, alternative='two-sided', method='auto')
        rank_biserial = 2 * result.statistic / (len(nlight) * len(dlight)) - 1
        u_value = float(result.statistic)
        p_value = float(result.pvalue)
    else:
        u_value = p_value = rank_biserial = np.nan

    def sample_fields(sensor, values):
        group = current[(current['sensor'] == sensor) & current[metric].notna()]
        if metric in ['temporal_tau_s', 'rise_time_s']:
            n_roi = int(group['n_temporal_fit_roi'].sum())
        elif metric == 'spatial_tau_px':
            n_roi = int(group['n_spatial_fit_roi'].sum())
        elif metric == 'first_lick_rho':
            n_roi = int(group['n_first_lick_trials'].sum())
        elif metric == 'response_amplitude_up_dff':
            n_roi = int(group['n_up_roi'].sum())
        else:
            n_roi = int(group['n_valid_roi'].sum())
        summary = finite_summary(values)
        return summary, int(group['animal'].nunique()), n_roi

    n_summary, n_mice, n_roi = sample_fields('nLight', nlight)
    d_summary, d_mice, d_roi = sample_fields('dLight', dlight)
    return {
        'event': event,
        'metric': metric,
        'label': label,
        'nlight_n_sessions': n_summary['n'],
        'nlight_n_mice': n_mice,
        'nlight_n_roi_or_trials': n_roi,
        'nlight_median': n_summary['median'],
        'nlight_q1': n_summary['q1'],
        'nlight_q3': n_summary['q3'],
        'dlight_n_sessions': d_summary['n'],
        'dlight_n_mice': d_mice,
        'dlight_n_roi_or_trials': d_roi,
        'dlight_median': d_summary['median'],
        'dlight_q1': d_summary['q1'],
        'dlight_q3': d_summary['q3'],
        'mann_whitney_u': u_value,
        'p_value': p_value,
        'rank_biserial_nlight_minus_dlight': rank_biserial,
    }


def formal_comparisons(session_metrics):
    specifications = [
        ('natural_run_onset', 'temporal_tau_s', 'decay tau (s)'),
        ('natural_run_onset', 'rise_time_s', 'rise time to peak (s)'),
        ('natural_run_onset', 'up_fraction', 'DA-Up ACR fraction'),
        ('natural_run_onset', 'response_amplitude_all_dff', 'all-ACR amplitude (dF/F)'),
        ('natural_run_onset', 'response_amplitude_up_dff', 'DA-Up amplitude (dF/F)'),
        ('natural_run_onset', 'first_lick_rho', 'first-lick Spearman rho'),
        ('opto_lc_stimulation', 'temporal_tau_s', 'decay tau (s)'),
        ('opto_lc_stimulation', 'up_fraction', 'releasing ACR fraction, t-test'),
        ('opto_lc_stimulation', 'up_fraction_wilcoxon', 'releasing ACR fraction, Wilcoxon'),
        ('opto_lc_stimulation', 'response_amplitude_all_dff', 'mean response index'),
    ]
    rows = [mann_comparison(session_metrics, *specification) for specification in specifications]
    comparisons = pd.DataFrame(rows)
    comparisons['q_value_bh'] = np.nan
    finite = comparisons['p_value'].notna()
    if finite.any():
        ordered = comparisons.loc[finite, 'p_value'].sort_values()
        ranks = np.arange(1, len(ordered) + 1)
        adjusted = ordered.to_numpy() * len(ordered) / ranks
        adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
        comparisons.loc[ordered.index, 'q_value_bh'] = np.minimum(adjusted, 1)
    return comparisons


#%% figures
def comparison_row(comparisons, event, metric):
    rows = comparisons[(comparisons['event'] == event) & (comparisons['metric'] == metric)]
    return rows.iloc[0] if len(rows) else None


def dot_panel(ax, session_metrics, comparisons, event, metric, ylabel, scale=1):
    rng = np.random.default_rng(42)
    current = session_metrics[session_metrics['event'] == event]
    for index, sensor in enumerate(['dLight', 'nLight']):
        values = current.loc[current['sensor'] == sensor, metric].dropna().to_numpy(dtype=float) * scale
        if not len(values):
            continue
        jitter = rng.uniform(-0.10, 0.10, len(values))
        ax.scatter(
            np.full(len(values), index) + jitter,
            values,
            s=19,
            color=SENSOR_COLOURS[sensor],
            edgecolor='white',
            linewidth=0.35,
            alpha=0.82,
            zorder=2,
        )
        q1, median, q3 = np.percentile(values, [25, 50, 75])
        ax.vlines(index, q1, q3, color='black', linewidth=1.2, zorder=3)
        ax.hlines(median, index - 0.18, index + 0.18, color='black', linewidth=1.5, zorder=3)
    result = comparison_row(comparisons, event, metric)
    if result is not None and np.isfinite(result['p_value']):
        u_value = result['mann_whitney_u']
        p_value = result['p_value']
        rank_biserial = result['rank_biserial_nlight_minus_dlight']
        ax.set_title(
            f'U={u_value:.1f}, p={format_p(p_value)}\n'
            f'rank-biserial={rank_biserial:.2f}',
            fontsize=8,
        )
    ax.set_xticks([0, 1], ['dLight', 'nLight'])
    ax.set_ylabel(ylabel)
    ax.spines[['top', 'right']].set_visible(False)


def plot_natural_summary(session_metrics, comparisons, session_profiles):
    fig, axes = plt.subplots(2, 3, figsize=(10, 6.2))
    dot_panel(axes[0, 0], session_metrics, comparisons, 'natural_run_onset', 'temporal_tau_s', 'decay tau (s)')
    dot_panel(axes[0, 1], session_metrics, comparisons, 'natural_run_onset', 'rise_time_s', 'time to peak (s)')
    dot_panel(axes[0, 2], session_metrics, comparisons, 'natural_run_onset', 'up_fraction', 'DA-Up ACRs (%)', 100)
    dot_panel(axes[1, 0], session_metrics, comparisons, 'natural_run_onset', 'response_amplitude_up_dff', 'DA-Up amplitude (% dF/F)', 100)
    dot_panel(axes[1, 1], session_metrics, comparisons, 'natural_run_onset', 'first_lick_rho', 'response versus first lick, rho')

    ax = axes[1, 2]
    xaxis = np.arange(IMAGING_RATE * (PROFILE_BEFORE_S + PROFILE_AFTER_S)) / IMAGING_RATE - PROFILE_BEFORE_S
    for sensor in ['dLight', 'nLight']:
        profiles = session_profiles[sensor]
        if not profiles:
            continue
        matrix = np.stack(list(profiles.values())) * 100
        median_profile = np.nanmedian(matrix, axis=0)
        q1 = np.nanpercentile(matrix, 25, axis=0)
        q3 = np.nanpercentile(matrix, 75, axis=0)
        colour = SENSOR_COLOURS[sensor]
        ax.plot(xaxis, median_profile, color=colour, linewidth=1.4, label=f'{sensor}, {len(matrix)} sessions')
        ax.fill_between(xaxis, q1, q3, color=colour, alpha=0.14, linewidth=0)
        fit = fit_temporal_profile(median_profile / 100)
        if fit:
            ax.plot(fit['fit_t_s'], fit['fit_y_percent'], color=colour, linestyle='--', linewidth=1)
    ax.axvline(0, color='0.6', linewidth=0.8)
    ax.set(xlabel='time from run onset (s)', ylabel='DA-Up signal (% dF/F)', title='session profiles and exponential fits')
    ax.legend(frameon=False, fontsize=7)
    ax.spines[['top', 'right']].set_visible(False)
    fig.tight_layout()
    return fig


def plot_opto_summary(session_metrics, comparisons):
    fig, axes = plt.subplots(1, 3, figsize=(9.2, 3))
    dot_panel(axes[0], session_metrics, comparisons, 'opto_lc_stimulation', 'temporal_tau_s', 'decay tau (s)')
    dot_panel(axes[1], session_metrics, comparisons, 'opto_lc_stimulation', 'up_fraction', 'releasing ACRs (%)', 100)
    dot_panel(axes[2], session_metrics, comparisons, 'opto_lc_stimulation', 'response_amplitude_all_dff', 'mean response index')
    fig.tight_layout()
    return fig


def bounded_spatial_fit(distance, response):
    distance = np.asarray(distance, dtype=float)
    response = np.asarray(response, dtype=float)
    valid = np.isfinite(distance) & np.isfinite(response)
    distance, response = distance[valid], response[valid]
    if len(response) < 6:
        return None
    baseline = float(np.min(response))
    shifted = response - baseline
    try:
        fitted, _ = curve_fit(
            lambda x, amplitude, tau: amplitude * np.exp(-x / tau),
            distance,
            shifted,
            bounds=([0, 0.5], [np.inf, 50]),
        )
    except (RuntimeError, ValueError):
        return None
    return baseline, fitted[1], fitted[0] * np.exp(-distance / fitted[1]) + baseline


def plot_spatial_curves(natural_curves, opto_curves, selected_sessions):
    fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.2))
    ax = axes[0]
    if len(natural_curves):
        session_curve = natural_curves.groupby(['session', 'distance_px'])['response'].mean().reset_index()
        population = session_curve.groupby('distance_px')['response'].median().reset_index()
        fit = fit_spatial_profile(population['distance_px'], population['response'])
        ax.plot(population['distance_px'], population['response'], 'o', color=SENSOR_COLOURS['dLight'], label='dLight data')
        if fit:
            tau_px = fit['tau_px']
            ax.plot(fit['distance_px'], fit['fit_response'], '--', color=SENSOR_COLOURS['dLight'], label=f'fit, tau={tau_px:.2f} px')
    ax.text(0.98, 0.95, 'nLight unavailable\nonly k=0 was saved', transform=ax.transAxes, ha='right', va='top', fontsize=8)
    ax.set(title='natural run onset', xlabel='distance from axon (pixels)', ylabel='response amplitude')
    ax.legend(frameon=False, fontsize=7, loc='lower left')
    ax.spines[['top', 'right']].set_visible(False)

    ax = axes[1]
    for sensor in ['dLight', 'nLight']:
        current = opto_curves[opto_curves['sensor'] == sensor].copy()
        if sensor == 'nLight':
            current = current[current['recname'].isin(selected_sessions)]
        session_curve = current.groupby(['recname', 'dilation_px'])['mean_RI'].mean().reset_index()
        population = session_curve.groupby('dilation_px')['mean_RI'].median().reset_index()
        fit = bounded_spatial_fit(population['dilation_px'], population['mean_RI'])
        colour = SENSOR_COLOURS[sensor]
        ax.plot(population['dilation_px'], population['mean_RI'], 'o', color=colour, label=f'{sensor} data')
        if fit:
            baseline, tau, prediction = fit
            ax.plot(population['dilation_px'], prediction, '--', color=colour, label=f'{sensor} fit, tau={tau:.2f} px')
    ax.set(title='LC optogenetic stimulation, 0 to 1 s', xlabel='distance from releasing ACRs (pixels)', ylabel='mean response index')
    ax.legend(frameon=False, fontsize=7)
    ax.spines[['top', 'right']].set_visible(False)
    fig.tight_layout()
    return fig


def plot_threshold_sensitivity(sensitivity):
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3))
    for ax, sensor in zip(axes, ['dLight', 'nLight']):
        current = sensitivity[sensitivity['sensor'] == sensor]
        for effect_size, group in current.groupby('effect_size_threshold'):
            group = group.sort_values('shuffle_percentile')
            ax.plot(group['shuffle_percentile'], group['median_up_fraction'] * 100, marker='o', label=f'effect size > {effect_size:g}')
        ax.set(title=sensor, xlabel='shuffle percentile', ylabel='median DA-Up ACRs (%)', xticks=SHUFFLE_SENSITIVITY)
        ax.legend(frameon=False, fontsize=7)
        ax.spines[['top', 'right']].set_visible(False)
    fig.tight_layout()
    return fig


def save_figure(fig, stem):
    for suffix in ['.png', '.svg']:
        fig.savefig(stem.with_suffix(suffix), dpi=300, bbox_inches='tight')
    plt.close(fig)


#%% inclusion tables and reports
def inclusion_table(session_metrics):
    rows = []
    for row in session_metrics.itertuples(index=False):
        common = {
            'sensor': row.sensor,
            'event': row.event,
            'animal': row.animal,
            'session': row.session,
        }
        rows.append({
            **common,
            'metric': 'recruitment_and_amplitude',
            'included': bool(row.n_valid_roi > 0),
            'reason': '' if row.n_valid_roi > 0 else 'no valid ACR',
            'n_roi_or_trials': int(row.n_valid_roi),
        })
        rows.append({
            **common,
            'metric': 'temporal_kinetics',
            'included': bool(np.isfinite(row.temporal_tau_s)),
            'reason': '' if np.isfinite(row.temporal_tau_s) else getattr(row, 'temporal_exclusion_reason', '') or 'no ACR or session fit passed temporal QC',
            'n_roi_or_trials': int(row.n_temporal_fit_roi),
        })
        rows.append({
            **common,
            'metric': 'spatial_kinetics',
            'included': bool(np.isfinite(row.spatial_tau_px)),
            'reason': '' if np.isfinite(row.spatial_tau_px) else getattr(row, 'spatial_exclusion_reason', '') or 'no spatial fit passed QC',
            'n_roi_or_trials': int(row.n_spatial_fit_roi),
        })
        if row.event == 'natural_run_onset':
            rows.append({
                **common,
                'metric': 'first_lick_correlation',
                'included': bool(np.isfinite(row.first_lick_rho)),
                'reason': '' if np.isfinite(row.first_lick_rho) else 'fewer than eight valid first-lick trials or no DA-Up ACR',
                'n_roi_or_trials': int(row.n_first_lick_trials),
            })
    return pd.DataFrame(rows)


def interpretation(row):
    if not np.isfinite(row.p_value):
        return 'not testable from the available cohort'
    if row.q_value_bh < 0.05:
        direction = 'higher' if row.rank_biserial_nlight_minus_dlight > 0 else 'lower'
        return f'nLight was {direction}'
    return 'no evidence of a sensor difference; this is not an equivalence result'


def comparison_markdown_row(row, percent=False):
    scale = 100 if percent else 1
    nlight = (
        f'{format_number(row.nlight_median * scale)} '
        f'[{format_number(row.nlight_q1 * scale)}, {format_number(row.nlight_q3 * scale)}]; '
        f'{row.nlight_n_sessions} sessions, {row.nlight_n_mice} mice, '
        f'{row.nlight_n_roi_or_trials} ACRs/trials'
    )
    dlight = (
        f'{format_number(row.dlight_median * scale)} '
        f'[{format_number(row.dlight_q1 * scale)}, {format_number(row.dlight_q3 * scale)}]; '
        f'{row.dlight_n_sessions} sessions, {row.dlight_n_mice} mice, '
        f'{row.dlight_n_roi_or_trials} ACRs/trials'
    )
    test = (
        f'U={format_number(row.mann_whitney_u, 1)}, p={format_p(row.p_value)}, '
        f'BH q={format_p(row.q_value_bh)}, rank-biserial='
        f'{format_number(row.rank_biserial_nlight_minus_dlight, 2)}'
    )
    return f'| {row.label} | {nlight} | {dlight} | {test} | {interpretation(row)} |'


def summary_text(session_metrics, roi_metrics, comparisons, sensitivity, selected_sessions):
    natural = session_metrics[session_metrics['event'] == 'natural_run_onset']
    natural_roi = roi_metrics[roi_metrics['event'] == 'natural_run_onset']
    dlight_temporal = natural_roi[
        (natural_roi['sensor'] == 'dLight') & natural_roi['temporal_included'].fillna(False)
    ]
    dlight_fraction = natural.loc[natural['sensor'] == 'dLight', 'up_fraction']
    dlight_natural_tau = finite_summary(dlight_temporal['temporal_tau_s'])

    opto = session_metrics[session_metrics['event'] == 'opto_lc_stimulation']
    dlight_opto_tau = finite_summary(opto.loc[opto['sensor'] == 'dLight', 'temporal_tau_s'])
    dlight_opto_release = finite_summary(opto.loc[opto['sensor'] == 'dLight', 'up_fraction'])

    all_time = pd.read_csv(OPTO_TABLE_ROOT / 'roi_neuropil_timecourse_metrics.csv')
    broad_time = finite_summary(all_time.loc[all_time['sensor'].eq('nLight'), 'roi_tau_s'])
    nlight_natural = natural[natural['sensor'].eq('nLight')]
    nlight_opto = opto[opto['sensor'].eq('nLight')]
    nlight_up_total = int(nlight_natural['n_up_roi'].sum())
    nlight_valid_total = int(nlight_natural['n_valid_roi'].sum())
    nlight_fraction_median = float(nlight_natural['up_fraction'].median())
    nlight_opto_temporal_n = int(nlight_opto['temporal_tau_s'].notna().sum())
    dlight_natural_tau_median = dlight_natural_tau['median']
    dlight_natural_tau_n = dlight_natural_tau['n']
    dlight_opto_tau_median = dlight_opto_tau['median']
    dlight_opto_tau_n = dlight_opto_tau['n']
    dlight_opto_release_median = dlight_opto_release['median']
    broad_time_n = broad_time['n']
    broad_time_median = broad_time['median']

    lines = [
        '# nLight versus dLight summary',
        '',
        '## Main result',
        '',
        (
            'Natural run-onset nLight responses were selectively recruited in the prespecified cohort: '
            f'{nlight_up_total} of {nlight_valid_total} valid ACRs '
            f'were Up and the median session recruited {100 * nlight_fraction_median:.2f}% of '
            'valid ACRs. The session-level Mann-Whitney results below describe a lower nLight '
            'recruitment and amplitude in this table, while the kinetic estimate applies only to '
            'the subset of sessions with at least one DA-Up ACR whose fit passed quality control. '
            'These session-level comparisons are not the Chapter 3 biological inference.'
        ),
        '',
        (
            'The strict optogenetic Chapter 2 temporal fitting procedure produced a usable nLight time '
            f'constant in only {nlight_opto_temporal_n} '
            f'of the {len(selected_sessions)} prespecified sessions. That value is reported, but a '
            'one-session group is not a stable population estimate. Across every metric, a '
            'non-significant result is treated as '
            'absence of evidence rather than evidence of equivalence.'
        ),
        '',
        '## Session-level comparisons',
        '',
        'Values are median [IQR]. The original comparison generator treats session as the '
        'independent unit; the rank-biserial sign is positive when nLight is higher. BH q values '
        'correct that generator family. For Chapter 3, session tests are cohort descriptions and '
        'animal-level exact tests and animal-then-session bootstrap intervals govern biological '
        'inference.',
        '',
        '| Metric | nLight | dLight | Test | Reading |',
        '|---|---|---|---|---|',
    ]
    percent_metrics = {'up_fraction', 'up_fraction_wilcoxon'}
    for row in comparisons.itertuples(index=False):
        lines.append(comparison_markdown_row(row, row.metric in percent_metrics))

    lines.extend([
        '',
        '## Reconciliation with the Chapter 2 reference values',
        '',
        (
            f'- Natural dLight decay: pooled DA-Up ACR median {dlight_natural_tau_median:.4f} s '
            f'(n={dlight_natural_tau_n} fits), which reproduces 1.42 s after rounding.'
        ),
        (
            f'- Natural dLight recruitment: session mean {100 * dlight_fraction.mean():.2f}% '
            f'and SEM {100 * dlight_fraction.sem():.2f}% across {len(dlight_fraction)} sessions, '
            'which reproduces 21.4 +/- 1.9%.'
        ),
        (
            f'- Optogenetic dLight decay: session median {dlight_opto_tau_median:.4f} s '
            f'(n={dlight_opto_tau_n}), reproducing 1.37 s.'
        ),
        (
            f'- Optogenetic dLight recruitment: the current untreated t-test implementation gives '
            f'{100 * dlight_opto_release_median:.2f}% median releasing ACRs, not 87.5%. The '
            'historical Wilcoxon histogram script also mixed untreated and nepicastat sessions and '
            'discarded zero-release sessions; 87.5% could not be reproduced from the present outputs.'
        ),
        '',
        '## Metrics that the data do not support',
        '',
        '- Natural nLight spatial tau is unavailable because the saved analysis products contain only '
        'dilation k=0. Reconstructing k=1 through 10 would require raw movie re-extraction, not a fit to '
        'the existing ACR time series.',
        '- Spatial fits are retained in native image pixels in the supporting tables and diagnostic '
        'figure. No physical-distance result is reported because the image calibration remains unresolved.',
        '- Axonal calcium lag is unavailable because these nLight sessions contain sensor, red structural '
        'and neuropil/regressor channels, but no simultaneous axonal GCaMP channel. The red regressor was '
        'not relabelled as calcium.',
        '- No explicitly labelled nLight nepicastat session was found. The repository contains a separate '
        'dLight nepicastat cohort, which is not an nLight drug comparison.',
        '',
        '## Additional context and sensitivity',
        '',
        (
            f'The broader repository nLight optogenetic table, which includes sessions outside Jingyu\'s '
            f'{len(selected_sessions)}-session cohort, contains {broad_time_n} temporal fits with median '
            f'{broad_time_median:.3f} s. This is useful context, but those sessions were not substituted '
            'for the prespecified cohort.'
        ),
        '',
        'The threshold-sensitivity table and figure show the joint effect of shuffle percentiles 90, 95 '
        'and 99 with effect-size cut-offs 0, 0.05 and 0.1. The first-lick analysis is secondary and uses '
        'within-session Spearman correlations on trials with a first lick later than 0.5 s. Raw dF/F or '
        'response-index amplitudes are not molar calibrations, so their cross-sensor comparison should be '
        'read as an empirical signal comparison rather than a concentration ratio.',
        '',
        '## Files',
        '',
        '- `per_session_metrics.csv`: one row per sensor, event and session.',
        '- `per_roi_fits.csv`: natural-run temporal and native-pixel dLight spatial ACR fits with QC flags.',
        '- `comparisons.csv`: U tests, effect sizes and BH q values.',
        '- `inclusion_exclusion.csv`: explicit metric-wise inclusion and reasons.',
        '- `threshold_sensitivity.csv`: recruitment sensitivity grid.',
        '- `methods.md`: paths and parameters.',
        '- `figures/`: PNG and editable SVG sources.',
        '',
    ])
    return '\n'.join(lines)


def methods_text(session_metrics, selected_sessions):
    nlight = session_metrics[session_metrics['sensor'] == 'nLight']
    dlight = session_metrics[session_metrics['sensor'] == 'dLight']
    return f'''# Methods and provenance

## Cohorts and inputs

The prespecified nLight cohort contains {len(selected_sessions)} sessions from {nlight['animal'].nunique()} mice: {', '.join(selected_sessions)}. Natural-run inputs came from `{NLIGHT_ROOT}` and the matched behavioural pickles under `{NATURAL_NLIGHT_BEHAVIOUR_ROOT}`. For each nLight session, the analysis rebuilt the combined table in memory from the current green and red response tables, baseline arrays and axon masks, then reapplied the current pooling and classification rules. It did not consume the saved combined tables because 16 A193 component tables dated 20 August 2026 were newer than their combined tables dated 10 August 2026. The earlier `nlight_vs_dlight_complete` output is therefore retained only as a stale-input provenance record. The dLight natural comparison contains {dlight.loc[dlight.event.eq('natural_run_onset'), 'session'].nunique()} sessions from {dlight.loc[dlight.event.eq('natural_run_onset'), 'animal'].nunique()} mice, selected by the stored lick, speed, valid-trial, trial-length and motion-correlation rules in Jingyu's `Dbh_dlight` recording list. Its pooled table was verified against all 57 constituent combined tables, with 11,502 rows and no mismatches. Optogenetic metrics came from `{OPTO_TABLE_ROOT}` and their processed session arrays under this repository's `data/hpc_dlight/lc_opto` and `data/hpc_nlight/lc_opto` roots.

An ACR is a 16 by 16 pixel bin. Spatial diagnostics retain native image pixels because the physical image calibration remains unresolved; no physical-distance conversion is applied here.

## Natural run onset

Events were rewarded, full-stop, non-stimulation run onsets with a complete -2 to +4 s window outside shutter closure. The green sensor trace was corrected with neuropil coefficient 0.2 and trial-wise regression against the green anatomical regressor and red channel, requiring at least 80 shared finite frames. dF/F used a rolling 20th percentile over 9,000 frames, the trace filter used global mean plus five standard deviations, and smoothing before response quantification used sigma one frame.

Valid ACRs were non-edge, had finite green and red profiles, minimum green baseline above 2 and minimum red baseline above 1. DA-Up required response amplitude from 0 to 1 s minus -1 to 0 s above the ACR-specific 95th percentile of 1,000 circular shifts and effect size above 0.05. Sensitivity used shuffle percentiles {SHUFFLE_SENSITIVITY} crossed with effect-size thresholds {EFFECT_SIZE_SENSITIVITY}.

Temporal fits followed `compute_tau_with_qc`: each DA-Up mean profile was Gaussian-smoothed with sigma {TEMPORAL_SMOOTH_SIGMA} frames, its maximum was located, and `A exp(-t/tau) + C` was fitted from that peak through +4 s with C fixed to the minimum fitted value. Fits required positive tau, R squared above {TEMPORAL_R2_MIN} and a peak strictly between {TEMPORAL_PEAK_RANGE_S[0]} and {TEMPORAL_PEAK_RANGE_S[1]} s. One median across included ACR fits represented each session.

dLight natural spatial fits used non-overlapping grid-free one-pixel dilation rings k=0 through 10. ACRs were selected as DA-Up at k=0, each ACR response-amplitude series was fitted with the same fixed-minimum exponential convention, and positive-tau fits required at least six points and R squared above {SPATIAL_R2_MIN}. nLight natural spatial tau was not calculated because only k=0 was saved.

First-lick correlations used the mean DA-Up ACR response on each retained trial, defined as 0 to 1 s minus -1 to 0 s, and the first lick later than 0.5 s. Spearman's rho required at least eight valid trials per session.

## LC optogenetic stimulation

The Chapter 2 implementations identified releasing ACRs with a one-sample greater-than-zero response test at alpha 0.05 plus mean response index above 0.1. The recruitment sensitivity repeated the dedicated Wilcoxon branch while retaining zero-release sessions. Temporal session-mask fits required R squared above 0.7. Spatial fits used 0 to 1 s response-index rings from zero through 20 pixels, subtracted the minimum response across distance, fitted a bounded exponential with tau 0.5 to 50 pixels, and excluded sessions with more than 50% expression-mask pixel reduction.

## Inference

The original comparison generator used session as its independent unit. nLight and dLight used different animals, so its available comparisons used a two-sided Mann-Whitney U test rather than a paired test. Rank-biserial effect size was `2U/(n_nLight*n_dLight)-1`; positive values mean higher nLight. Benjamini-Hochberg q values cover that generator's family of testable metrics. ROI-level fits are retained in the source tables but are not counted as independent animals.

That source-table convention is retained for reproducibility, but it is not the Chapter 3 biological inference rule. For the Chapter 3 analysis, animal is the biological inference unit and session is a repeated observation; Chapter 3 uses animal-median exact tests and animal-then-session bootstrap intervals, with session-level tests shown only as cohort descriptions.

'''


#%% run and write
def validate_results(session_metrics, roi_metrics, comparisons):
    dlight_natural = roi_metrics[
        (roi_metrics['sensor'] == 'dLight')
        & roi_metrics['event'].eq('natural_run_onset')
        & roi_metrics['temporal_included'].fillna(False)
    ]
    if not np.isclose(dlight_natural['temporal_tau_s'].median(), 1.42, atol=0.01):
        raise AssertionError('Natural dLight temporal benchmark did not reproduce 1.42 s')
    dlight_opto = session_metrics[
        session_metrics['sensor'].eq('dLight')
        & session_metrics['event'].eq('opto_lc_stimulation')
    ]
    if not np.isclose(dlight_opto['temporal_tau_s'].median(), 1.37, atol=0.01):
        raise AssertionError('Optogenetic dLight temporal benchmark did not reproduce 1.37 s')
    if comparisons.empty or session_metrics.empty or roi_metrics.empty:
        raise AssertionError('Analysis produced an empty required table')


def main(write=False):
    warnings.filterwarnings('ignore', category=RuntimeWarning)
    required = [
        OPTO_TABLE_ROOT / 'roi_neuropil_timecourse_metrics.csv',
        OPTO_TABLE_ROOT / 'spatial_tau_metrics.csv',
        OPTO_TABLE_ROOT / 'roi_release_metrics.csv',
        OPTO_TABLE_ROOT / 'spatial_dilation_curves.csv',
    ]
    require_inputs(required)
    selected_sessions = selected_nlight_sessions()
    if len(selected_sessions) != 26:
        raise AssertionError(f'Expected 26 prespecified nLight sessions, found {len(selected_sessions)}')

    natural_tables = load_natural_tables()
    natural_sessions, roi_metrics, session_profiles = natural_metrics(natural_tables)
    natural_sessions, roi_metrics, natural_curves = natural_dlight_spatial(natural_sessions, roi_metrics)
    natural_sessions = add_first_lick_correlations(natural_sessions, natural_tables)
    sensitivity = threshold_sensitivity(natural_tables)
    opto_sessions, opto_curves = load_opto_metrics(selected_sessions)
    session_metrics = pd.concat([natural_sessions, opto_sessions], ignore_index=True, sort=False)
    session_metrics = add_wilcoxon_release(session_metrics)
    comparisons = formal_comparisons(session_metrics)
    inclusion = inclusion_table(session_metrics)
    validate_results(session_metrics, roi_metrics, comparisons)

    if not write:
        print(comparisons.to_string(index=False))
        return
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    FIGURE_ROOT.mkdir(parents=True, exist_ok=True)

    session_metrics.sort_values(['event', 'sensor', 'session']).to_csv(
        OUTPUT_ROOT / 'per_session_metrics.csv', index=False
    )
    roi_metrics.sort_values(['event', 'sensor', 'session', 'roi_y', 'roi_x']).to_csv(
        OUTPUT_ROOT / 'per_roi_fits.csv', index=False
    )
    comparisons.to_csv(OUTPUT_ROOT / 'comparisons.csv', index=False)
    inclusion.to_csv(OUTPUT_ROOT / 'inclusion_exclusion.csv', index=False)
    sensitivity.to_csv(OUTPUT_ROOT / 'threshold_sensitivity.csv', index=False)

    save_figure(
        plot_natural_summary(session_metrics, comparisons, session_profiles),
        FIGURE_ROOT / 'natural_run_summary',
    )
    save_figure(
        plot_opto_summary(session_metrics, comparisons),
        FIGURE_ROOT / 'opto_summary',
    )
    save_figure(
        plot_spatial_curves(natural_curves, opto_curves, selected_sessions),
        FIGURE_ROOT / 'spatial_fit_curves',
    )
    save_figure(
        plot_threshold_sensitivity(sensitivity),
        FIGURE_ROOT / 'threshold_sensitivity',
    )

    (OUTPUT_ROOT / 'nlight_vs_dlight_summary.md').write_text(
        summary_text(session_metrics, roi_metrics, comparisons, sensitivity, selected_sessions),
        encoding='utf-8',
    )
    (OUTPUT_ROOT / 'methods.md').write_text(
        methods_text(session_metrics, selected_sessions), encoding='utf-8'
    )
    print(f'Wrote analysis to {OUTPUT_ROOT}')
    print(f'Wrote figures to {FIGURE_ROOT}')


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('--write', action='store_true', help='write new tables, reports and figures')
    arguments = parser.parse_args()
    main(write=arguments.write)
