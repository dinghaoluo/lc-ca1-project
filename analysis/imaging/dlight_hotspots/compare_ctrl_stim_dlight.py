# -*- coding: utf-8 -*-
'''
Created on Fri Nov  6 18:42:11 2026

compare run-onset LC stimulation with interleaved control trials in dLight ACRs

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
OUTPUT_ROOT = pp.DLIGHT_HOTSPOTS_OUTPUT_ROOT
FIGURE_ROOT = pp.DLIGHT_HOTSPOTS_FIGURE_ROOT
INVENTORY_PATH = OUTPUT_ROOT / 'ctrl_stim_inventory.csv'
EVENT_PATH = OUTPUT_ROOT / 'ctrl_stim_events.csv'

if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))
if str(JINGYU_ROOT) not in sys.path:
    sys.path.insert(0, str(JINGYU_ROOT))

import behaviour_functions as bf  # noqa: E402
import project_paths as pp  # noqa: E402
from common.mask.generate_masks import generate_axon_dlight_masks_wrap  # noqa: E402
from common.utils_basic import trace_filter  # noqa: E402
from dlight_imaging.regression.trace_extraction import (  # noqa: E402
    TraceExtractionBatch,
)
from dlight_imaging.regression.utils_regression import (  # noqa: E402
    single_trial_regression_parallel,
)


SAMP_FREQ = 30
GRID_SIZE = 16
GRID_N = 32
PRE_FRAMES = 90
POST_FRAMES = 120
NULL_PERCENTILES = (90, 95, 99)
MIN_FINITE_FRAMES = 100
TRIAL_COLUMNS = [
    'rec_id', 'animal', 'condition', 'event_index', 'trial_index',
    'onset_frame', 'offset_frame', 'duration_frames',
    'response_start_frame', 'response_end_frame',
    'grid_y', 'grid_x', 'delta', 'null_delta', 'red_delta',
    'response_valid', 'red_valid',
    'null_threshold_90', 'null_threshold_95', 'null_threshold_99',
    'responding_90', 'responding_95', 'responding_99',
]


#%% movie masking and input tables
class MaskedMovie:
    def __init__(self, movie, frame_mask):
        self.movie = movie
        self.frame_mask = np.asarray(frame_mask, dtype=bool)

    @property
    def shape(self):
        return self.movie.shape

    def __getitem__(self, index):
        values = np.asarray(self.movie[index], dtype=np.float64).copy()
        if isinstance(index, slice):
            frame_indices = np.arange(*index.indices(self.shape[0]))
            values[self.frame_mask[frame_indices], ...] = np.nan
        elif np.isscalar(index) and self.frame_mask[int(index)]:
            values[...] = np.nan
        return values


def final_sessions():
    inventory = pd.read_csv(INVENTORY_PATH)
    selected = inventory.loc[
        inventory['qualifies_ctrl_stim'].astype(bool)
        & inventory['session_lists'].eq('pathdLightLCOpto')
    ].copy()
    if len(selected) != 40:
        raise ValueError(f'expected 40 comparison-2 sessions, found {len(selected)}')
    return selected


def load_events(recname):
    events = pd.read_csv(EVENT_PATH)
    events = events.loc[
        events['session'].eq(recname)
        & events['qualifies_ctrl_stim'].astype(bool)
    ].copy()
    if events.empty:
        raise ValueError(f'{recname}: no detected stimulation events')
    events['onset_frame'] = events['onset_frame'].astype(int)
    events['offset_frame'] = events['offset_frame'].astype(int)
    return events


def load_behaviour(row):
    txt_path = Path(row['behaviour_path'])
    behaviour = bf.process_behavioural_data_imaging(txt_path)
    if len(behaviour['run_onset_frames']) == 0:
        raise ValueError(f'{row.session}: behaviour parser returned no run onsets')
    return behaviour


def trial_statements(behaviour):
    return [str(statement[15]) for statement in behaviour['trial_statements']]


def valid_run_onsets(behaviour, n_frames):
    run_onsets = np.asarray(behaviour['run_onset_frames'], dtype=float)
    methods = trial_statements(behaviour)
    n_trials = min(len(run_onsets), len(methods))
    out = []
    for trial in range(n_trials):
        onset = run_onsets[trial]
        if (
                np.isfinite(onset)
                and onset - PRE_FRAMES >= 0
                and onset + POST_FRAMES <= n_frames
        ):
            out.append((trial, int(onset), methods[trial]))
    return out


def build_stim_mask(events, n_frames):
    frame_mask = np.zeros(n_frames, dtype=bool)
    for _, event in events.iterrows():
        onset = int(event['onset_frame'])
        offset = int(event['offset_frame'])
        if onset < 0 or offset >= n_frames or offset < onset:
            raise ValueError('detected stimulation edge lies outside the movie')
        frame_mask[onset:offset + 1] = True
    return frame_mask


#%% masks, extraction and regression
def ensure_masks(recname, processed_path, mask_path):
    expected = [
        mask_path / 'dilated_global_axon_k=0.npy',
        mask_path / 'global_dlight_mask_enhanced.npy',
        mask_path / 'dlight_regressor_fiber_dilation_k=0.npy',
        mask_path / 'fiber_neuropil_masks_dilation_k=0.npy',
    ]
    if not all(path.exists() for path in expected):
        red_ref = np.load(processed_path / f'{recname}_ref_mat_ch2.npy')
        green_ref = np.load(processed_path / f'{recname}_ref_mat_ch1.npy')
        mask_path.mkdir(parents=True, exist_ok=True)
        generate_axon_dlight_masks_wrap(
            red_ref,
            green_ref,
            mask_path,
            p_fiber_masks=processed_path,
            fiber_mask_name=f'{recname}_ROI_dict.npy',
            dilation_steps=(0,),
            grid_size=GRID_SIZE,
            constrain_to_grid=False,
            visualize_neu=False,
            visualize_dilation=False,
            visualize_regressor=False,
        )


def extract_and_regress(recname, row, events, behaviour, movie_mask, output_path):
    plane = Path(row['registration_plane0'])
    ops = np.load(plane / 'ops.npy', allow_pickle=True).item()
    n_frames = int(ops['nframes'])
    height, width = int(ops['Ly']), int(ops['Lx'])
    processed_path = Path(row['roi_path']).parent
    mask_path = output_path / 'masks'
    ensure_masks(recname, processed_path, mask_path)

    raw_path = output_path / f'{recname}_raw_traces.npz'
    if not raw_path.exists():
        green_movie = np.memmap(
            plane / 'data.bin', mode='r', dtype='int16',
            shape=(n_frames, height, width),
        )
        red_movie = np.memmap(
            plane / 'data_chan2.bin', mode='r', dtype='int16',
            shape=(n_frames, height, width),
        )
        batch = TraceExtractionBatch((height, width))
        axon = np.load(mask_path / 'dilated_global_axon_k=0.npy')
        dlight = np.load(mask_path / 'global_dlight_mask_enhanced.npy')
        regressor = np.load(mask_path / 'dlight_regressor_fiber_dilation_k=0.npy')
        neuropil = np.load(mask_path / 'fiber_neuropil_masks_dilation_k=0.npy')
        batch.add_grid(
            raw_path,
            axon,
            axon,
            dlight,
            regressor,
            neuropil,
            grid_free=True,
            grid_size=GRID_SIZE,
        )
        batch.run(
            MaskedMovie(green_movie, movie_mask),
            MaskedMovie(red_movie, movie_mask),
            batch_size=256,
            n_jobs=1,
        )

    regression_path = output_path / 'single_trial_regression'
    corrected_path = regression_path / 'corrected_dlight_trace.npy'
    if not corrected_path.exists():
        with np.load(raw_path) as raw:
            original = raw['original_dlight']
            neuropil = raw['neuropil_dlight']
            regressor = raw['dlight_regressor']
            red = raw['red_trace']
        corrected_input = original - 0.2 * neuropil
        run_events = [
            onset for _, onset, _ in valid_run_onsets(behaviour, n_frames)
        ]
        if not run_events:
            raise ValueError(f'{recname}: no complete run-onset regression windows')
        regression_path.mkdir(parents=True, exist_ok=True)
        single_trial_regression_parallel(
            original_green_traces=corrected_input,
            green_regressor_traces=regressor,
            red_regressor_traces=red,
            event_frames=np.asarray(run_events, dtype=int),
            out_dir=regression_path,
            pre_frames=PRE_FRAMES,
            post_frames=POST_FRAMES,
            min_finite_frames=MIN_FINITE_FRAMES,
            n_jobs=1,
        )
        np.save(regression_path / 'dlight_regressor_trace.npy', regressor)
        np.save(regression_path / 'red_trace.npy', red)

    with np.load(raw_path) as raw:
        red = raw['red_trace']
    corrected = np.load(corrected_path, mmap_mode='r')
    corrected = np.asarray(corrected).reshape(GRID_N * GRID_N, n_frames)
    red = np.asarray(red).reshape(GRID_N * GRID_N, n_frames)
    dff_path = regression_path / 'dff_corrected_dlight.npy'
    red_dff_path = regression_path / 'dff_red.npy'
    if dff_path.exists() and red_dff_path.exists():
        dff = np.load(dff_path, mmap_mode='r')
        red_dff = np.load(red_dff_path, mmap_mode='r')
    else:
        valid = np.isfinite(corrected).any(axis=1) & np.isfinite(red).any(axis=1)
        dff = np.full(corrected.shape, np.nan, dtype=np.float32)
        red_dff = np.full(red.shape, np.nan, dtype=np.float32)
        dff[valid] = percentile_dff_cpu(corrected[valid])
        red_dff[valid] = percentile_dff_cpu(red[valid])
        np.save(dff_path, dff.astype(np.float32))
        np.save(red_dff_path, red_dff.astype(np.float32))
    return (
        dff,
        red_dff,
        np.isfinite(dff).any(axis=1) & np.isfinite(red_dff).any(axis=1),
        n_frames,
    )


def percentile_dff_cpu(trace):
    finite = np.isfinite(trace)
    filled = np.where(finite, trace, np.inf)
    half_window = 9000 // 2
    baseline = np.vstack([
        pd.Series(np.pad(row, (half_window, half_window), mode='reflect'))
        .rolling(9000, center=True, min_periods=1)
        .quantile(0.2)
        .to_numpy()[half_window:half_window + trace.shape[1]]
        for row in filled
    ])
    dff = (trace - baseline) / (baseline + 1e-8)
    dff[~finite] = np.nan
    threshold = np.nanmean(dff) + 5 * np.nanstd(dff)
    dff = np.vstack([
        trace_filter(row, fix_thresh=threshold)
        for row in dff
    ])
    masked = ~np.isfinite(dff)
    smoothed = gaussian_filter1d(
        np.nan_to_num(dff, nan=0.0),
        sigma=1,
        axis=1,
        mode='reflect',
    )
    smoothed[masked] = np.nan
    return smoothed.astype(np.float32)


#%% windows and measures
def nanmean_window(trace, start, stop):
    values = trace[:, start:stop]
    count = np.isfinite(values).sum(axis=1)
    total = np.nansum(values, axis=1)
    output = np.full(trace.shape[0], np.nan, dtype=float)
    np.divide(total, count, out=output, where=count > 0)
    return output


def event_values(dff, red_dff, coords, specs, max_duration_frames):
    rows = []
    for spec in specs:
        onset = spec['onset_frame']
        if spec['condition'] == 'ctrl':
            response_start = onset
        else:
            response_start = onset + max_duration_frames + 1
        response_end = response_start + SAMP_FREQ
        baseline_start = onset - SAMP_FREQ
        baseline_end = onset
        null_early_start = onset - 3 * SAMP_FREQ
        null_early_end = onset - 2 * SAMP_FREQ
        null_late_start = onset - 2 * SAMP_FREQ
        null_late_end = onset - SAMP_FREQ
        delta = (
            nanmean_window(dff, response_start, response_end)
            - nanmean_window(dff, baseline_start, baseline_end)
        )
        null_delta = (
            nanmean_window(dff, null_late_start, null_late_end)
            - nanmean_window(dff, null_early_start, null_early_end)
        )
        red_delta = (
            nanmean_window(red_dff, response_start, response_end)
            - nanmean_window(red_dff, baseline_start, baseline_end)
        )
        for grid_y, grid_x in coords:
            grid = int(grid_y) * GRID_N + int(grid_x)
            rows.append({
                'rec_id': spec['rec_id'],
                'animal': spec['animal'],
                'condition': spec['condition'],
                'event_index': spec['event_index'],
                'trial_index': spec['trial_index'],
                'onset_frame': onset,
                'offset_frame': spec['offset_frame'],
                'duration_frames': spec['duration_frames'],
                'response_start_frame': response_start,
                'response_end_frame': response_end,
                'grid_y': int(grid_y),
                'grid_x': int(grid_x),
                'delta': delta[grid],
                'null_delta': null_delta[grid],
                'red_delta': red_delta[grid],
                'response_valid': False,
                'red_valid': False,
            })
    return pd.DataFrame(rows, columns=TRIAL_COLUMNS)


def add_thresholds(table):
    for percentile in NULL_PERCENTILES:
        threshold_column = f'null_threshold_{percentile}'
        if threshold_column in table:
            table = table.drop(columns=threshold_column)
        if f'responding_{percentile}' in table:
            table = table.drop(columns=f'responding_{percentile}')
        thresholds = (
            table.groupby(['grid_y', 'grid_x'])['null_delta']
            .apply(lambda values: np.nanpercentile(values.dropna(), percentile)
                   if values.notna().any() else np.nan)
            .rename(f'null_threshold_{percentile}')
            .reset_index()
        )
        table = table.merge(thresholds, on=['grid_y', 'grid_x'], how='left')
        table[f'responding_{percentile}'] = (
            table['delta'].notna()
            & table['null_delta'].notna()
            & table[f'null_threshold_{percentile}'].notna()
            & (table['delta'] > table[f'null_threshold_{percentile}'])
        )
    table['response_valid'] = table['delta'].notna() & table['null_delta'].notna()
    table['red_valid'] = table['red_delta'].notna()
    return table


def response_summary(table, condition, percentile):
    subset = table.loc[table['condition'].eq(condition)]
    values = []
    for _, group in subset.groupby('event_index'):
        valid = group['response_valid'].to_numpy(dtype=bool)
        if valid.any():
            values.append(100 * group.loc[valid, f'responding_{percentile}'].mean())
    values = np.asarray(values, dtype=float)
    if values.size == 0:
        return np.nan, np.nan, np.nan, 0
    median = float(np.nanmedian(values))
    sd = float(np.nanstd(values, ddof=1)) if values.size > 1 else np.nan
    mean = float(np.nanmean(values))
    cv = float(sd / mean) if values.size > 1 and mean else np.nan
    return median, sd, cv, int(values.size)


def session_summary(table, recname, animal, n_events_ctrl, n_events_stim, n_grids):
    row = {
        'rec_id': recname,
        'animal': animal,
        'n_ctrl_trials': n_events_ctrl,
        'n_stim_events': n_events_stim,
        'n_grids': n_grids,
    }
    medians = {}
    red_medians = {}
    for condition in ('ctrl', 'stim'):
        subset = table.loc[table['condition'].eq(condition)]
        per_grid = subset.groupby(['grid_y', 'grid_x'])['delta'].mean()
        per_grid = pd.to_numeric(per_grid, errors='coerce').dropna()
        medians[condition] = float(np.nanmedian(per_grid)) if len(per_grid) else np.nan
        row[f'{condition}_acr_delta_median'] = medians[condition]
        row[f'{condition}_acr_delta_n'] = int(len(per_grid))
        red_grid = subset.groupby(['grid_y', 'grid_x'])['red_delta'].mean()
        red_grid = pd.to_numeric(red_grid, errors='coerce').dropna()
        red_medians[condition] = (
            float(np.nanmedian(red_grid)) if len(red_grid) else np.nan
        )
        row[f'{condition}_red_delta_median'] = red_medians[condition]
        for percentile in NULL_PERCENTILES:
            median, sd, cv, n_trials = response_summary(table, condition, percentile)
            row[f'{condition}_responding_pct_median_{percentile}'] = median
            row[f'{condition}_responding_pct_sd_{percentile}'] = sd
            row[f'{condition}_responding_pct_cv_{percentile}'] = cv
            row[f'{condition}_responding_pct_n_trials_{percentile}'] = n_trials
    row['stim_minus_ctrl_acr_delta_median'] = (
        medians['stim'] - medians['ctrl']
        if np.isfinite(medians['stim']) and np.isfinite(medians['ctrl'])
        else np.nan
    )
    row['stim_minus_ctrl_red_delta_median'] = (
        red_medians['stim'] - red_medians['ctrl']
        if np.isfinite(red_medians['stim']) and np.isfinite(red_medians['ctrl'])
        else np.nan
    )
    for percentile in NULL_PERCENTILES:
        for metric in ('median', 'sd', 'cv'):
            ctrl = row[f'ctrl_responding_pct_{metric}_{percentile}']
            stim = row[f'stim_responding_pct_{metric}_{percentile}']
            row[f'stim_minus_ctrl_responding_pct_{metric}_{percentile}'] = (
                stim - ctrl if np.isfinite(stim) and np.isfinite(ctrl) else np.nan
            )
    return row


def paired_test(table, ctrl_column, stim_column):
    valid = table[[ctrl_column, stim_column]].dropna()
    if len(valid) < 2:
        return {'metric': stim_column, 'n_recordings': len(valid),
                'n_animals': table.loc[valid.index, 'animal'].nunique(),
                'statistic': np.nan, 'p_value': np.nan}
    result = wilcoxon(
        valid[stim_column].to_numpy(dtype=float),
        valid[ctrl_column].to_numpy(dtype=float),
        alternative='two-sided',
    )
    return {
        'metric': stim_column,
        'n_recordings': len(valid),
        'n_animals': table.loc[valid.index, 'animal'].nunique(),
        'statistic': float(result.statistic),
        'p_value': float(result.pvalue),
    }


#%% figures
def paired_plot(table, ctrl_column, stim_column, ylabel, path):
    fig, ax = plt.subplots(figsize=(3.0, 2.6))
    valid = table[[ctrl_column, stim_column]].dropna()
    for _, row in valid.iterrows():
        ax.plot([0, 1], [row[ctrl_column], row[stim_column]],
                color='0.65', linewidth=0.6)
    ax.scatter(np.zeros(len(valid)), valid[ctrl_column], s=10, color='0.25',
               label='ctrl')
    ax.scatter(np.ones(len(valid)), valid[stim_column], s=10, color='tab:blue',
               label='stim')
    ax.set_xticks([0, 1], ['ctrl', 'stim'])
    ax.set_ylabel(ylabel)
    ax.legend(frameon=False, fontsize=7)
    ax.spines[['top', 'right']].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=300)
    plt.close(fig)


def example_plot(table, recname, path):
    subset = table.loc[table['rec_id'].eq(recname)]
    fig, ax = plt.subplots(figsize=(3.2, 2.6))
    for condition, colour in (('ctrl', '0.25'), ('stim', 'tab:blue')):
        values = []
        for _, group in subset.loc[subset['condition'].eq(condition)].groupby('event_index'):
            valid = group['response_valid'].to_numpy(dtype=bool)
            if valid.any():
                values.append(100 * group.loc[valid, 'responding_95'].mean())
        ax.plot(values, linewidth=0.8, color=colour, label=condition)
    ax.set_xlabel('event') 
    ax.set_ylabel('responding ACRs (%)')
    ax.legend(frameon=False, fontsize=7)
    ax.spines[['top', 'right']].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=300)
    plt.close(fig)


def diagnostic_pass(summary):
    ctrl = summary['ctrl_red_delta_median'].dropna()
    stim = summary['stim_red_delta_median'].dropna()
    if len(ctrl) < 2:
        return False, np.nan, np.nan
    result = wilcoxon(stim.to_numpy(), ctrl.to_numpy(), alternative='two-sided')
    pooled = np.nanmedian(np.r_[ctrl.to_numpy(), stim.to_numpy()])
    return bool(abs(pooled) <= 0.01 and result.pvalue >= 0.05), float(result.statistic), float(result.pvalue)


def analyse_recording(recname, row):
    output_path = OUTPUT_ROOT / 'ctrl_stim_sessions' / recname
    output_path.mkdir(parents=True, exist_ok=True)
    behaviour = load_behaviour(row)
    events = load_events(recname)
    plane = Path(row['registration_plane0'])
    n_frames = int(np.load(plane / 'ops.npy', allow_pickle=True).item()['nframes'])
    movie_mask = build_stim_mask(events, n_frames)
    max_duration_frames = int(round(float(row['max_stim_duration_s']) * SAMP_FREQ))
    valid_runs = valid_run_onsets(behaviour, n_frames)
    ctrl_specs = [
        {
            'rec_id': recname,
            'animal': row['animal'],
            'condition': 'ctrl',
            'event_index': trial,
            'trial_index': trial,
            'onset_frame': onset,
            'offset_frame': onset,
            'duration_frames': 0,
        }
        for trial, onset, method in valid_runs if method == '0'
    ]
    stim_specs = [
        {
            'rec_id': recname,
            'animal': row['animal'],
            'condition': 'stim',
            'event_index': int(event['stimulation_index_zero_based']),
            'trial_index': int(event['preceding_trial_index_zero_based'])
            if pd.notna(event['preceding_trial_index_zero_based']) else -1,
            'onset_frame': int(event['onset_frame']),
            'offset_frame': int(event['offset_frame']),
            'duration_frames': int(event['duration_frames']),
        }
        for _, event in events.iterrows()
    ]
    for spec in stim_specs:
        if (
                spec['onset_frame'] - PRE_FRAMES < 0
                or spec['onset_frame'] + max_duration_frames + SAMP_FREQ + 1 > n_frames
        ):
            raise ValueError(f'{recname}: detected stimulation lacks a complete analysis window')
    dff, red_dff, valid_grid, n_frames = extract_and_regress(
        recname, row, events, behaviour, movie_mask, output_path,
    )
    coords = np.argwhere(valid_grid.reshape(GRID_N, GRID_N))
    specs = ctrl_specs + stim_specs
    table = event_values(
        dff, red_dff, coords, specs, max_duration_frames,
    )
    table = add_thresholds(table)
    table_path = output_path / f'{recname}_per_acr_trial.parquet'
    table.to_parquet(table_path, index=False)
    summary = session_summary(
        table,
        recname,
        row['animal'],
        len(ctrl_specs),
        len(stim_specs),
        len(coords),
    )
    return table, summary


def main(argv=None):
    parser = argparse.ArgumentParser(description='compare run-onset LC stimulation with ctrl trials')
    parser.add_argument('--recording', help='optional recording substring')
    args = parser.parse_args(argv)

    sessions = final_sessions()
    if args.recording:
        sessions = sessions.loc[sessions['session'].str.contains(args.recording)]
    if sessions.empty:
        raise ValueError('no sessions matched the requested filter')

    trial_tables = []
    summary_tables = []
    for _, row in sessions.iterrows():
        recname = row['session']
        print(f'loading {recname}', flush=True)
        table, summary = analyse_recording(recname, row)
        trial_tables.append(table)
        summary_tables.append(pd.DataFrame([summary]))
        print(
            f'  ctrl={summary["n_ctrl_trials"]}, '
            f'stim={summary["n_stim_events"]}, '
            f'grids={summary["n_grids"]}',
            flush=True,
        )

    per_trial = pd.concat(trial_tables, ignore_index=True)
    per_recording = pd.concat(summary_tables, ignore_index=True)
    per_trial.to_parquet(
        OUTPUT_ROOT / 'ctrl_stim_per_acr_trial.parquet',
        index=False,
    )
    per_recording.to_parquet(
        OUTPUT_ROOT / 'ctrl_stim_per_recording.parquet',
        index=False,
    )

    stats = [
        paired_test(
            per_recording,
            'ctrl_acr_delta_median',
            'stim_acr_delta_median',
        )
    ]
    for percentile in NULL_PERCENTILES:
        for metric in ('median', 'sd', 'cv'):
            stats.append(paired_test(
                per_recording,
                f'ctrl_responding_pct_{metric}_{percentile}',
                f'stim_responding_pct_{metric}_{percentile}',
            ))
    if len(per_recording) < 2:
        red_ok, red_statistic, red_pvalue = True, np.nan, np.nan
    else:
        red_ok, red_statistic, red_pvalue = diagnostic_pass(per_recording)
        if not red_ok:
            raise RuntimeError(
                f'tdTomato diagnostic failed: W={red_statistic}, p={red_pvalue}'
            )

    paired_plot(
        per_recording,
        'ctrl_acr_delta_median',
        'stim_acr_delta_median',
        'dLight Δ',
        FIGURE_ROOT / 'ctrl_stim_amplitude_paired.png',
    )
    paired_plot(
        per_recording,
        'ctrl_responding_pct_median_95',
        'stim_responding_pct_median_95',
        'responding ACRs (%)',
        FIGURE_ROOT / 'ctrl_stim_responding_paired_95.png',
    )
    example_plot(
        per_trial,
        per_recording.iloc[0]['rec_id'],
        FIGURE_ROOT / 'ctrl_stim_example_recording.png',
    )

    print(per_recording.to_string(index=False))
    print(pd.DataFrame(stats).to_string(index=False))


if __name__ == '__main__':
    main()
