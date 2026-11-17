# -*- coding: utf-8 -*-
'''
Created on Sun Nov  8 10:32:07 2026

diagnose the red control signal during LC stimulation

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
from scipy.stats import spearmanr, wilcoxon


#%% paths
repo_root = Path(__file__).resolve().parents[3]
OUTPUT_ROOT = pp.DLIGHT_HOTSPOTS_OUTPUT_ROOT
FIGURE_ROOT = pp.DLIGHT_HOTSPOTS_FIGURE_ROOT
INVENTORY_PATH = OUTPUT_ROOT / 'ctrl_stim_inventory.csv'
EVENT_PATH = OUTPUT_ROOT / 'ctrl_stim_events.csv'
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))
import behaviour_functions as bf
SESSION_ROOT = OUTPUT_ROOT / 'ctrl_stim_sessions'


#%% helpers
def load_behaviour(path):
    return bf.process_behavioural_data_imaging(path)


def mean_window(trace, start, stop):
    values = trace[start:stop]
    return float(np.nanmean(values)) if np.isfinite(values).any() else np.nan


def event_speed(behaviour, event_time_ms, max_duration_s):
    times = np.asarray(behaviour['upsampled_timestamps_ms'], dtype=float)
    speed = np.asarray(behaviour['upsampled_speed_cm_s'], dtype=float)
    response_start = event_time_ms + max_duration_s * 1000
    response_end = response_start + 1000
    baseline = speed[(times >= event_time_ms - 1000) & (times < event_time_ms)]
    response = speed[(times >= response_start) & (times < response_end)]
    return {
        'speed_baseline': np.nanmean(baseline) if baseline.size else np.nan,
        'speed_response': np.nanmean(response) if response.size else np.nan,
        'speed_change': (
            np.nanmean(response) - np.nanmean(baseline)
            if baseline.size and response.size else np.nan
        ),
        'speed_response_fraction_above_10': (
            np.mean(response >= 10) if response.size else np.nan
        ),
    }


def event_motion(xoff, yoff, onset, max_duration_frames):
    response_start = onset + max_duration_frames
    response_end = response_start + 30
    baseline_start = max(0, onset - 30)
    baseline_end = onset
    x = np.asarray(xoff, dtype=float)
    y = np.asarray(yoff, dtype=float)
    if response_end > len(x):
        return {
            'motion_baseline': np.nan,
            'motion_response': np.nan,
            'motion_change': np.nan,
            'motion_step_response': np.nan,
        }
    baseline = np.sqrt(x[baseline_start:baseline_end] ** 2
                       + y[baseline_start:baseline_end] ** 2)
    response = np.sqrt(x[response_start:response_end] ** 2
                       + y[response_start:response_end] ** 2)
    steps = np.sqrt(np.diff(x[response_start:response_end]) ** 2
                    + np.diff(y[response_start:response_end]) ** 2)
    return {
        'motion_baseline': np.nanmean(baseline),
        'motion_response': np.nanmean(response),
        'motion_change': np.nanmean(response) - np.nanmean(baseline),
        'motion_step_response': np.nanmean(steps) if steps.size else np.nan,
    }


def session_rows(row, events):
    recname = row['session']
    behaviour = load_behaviour(row['behaviour_path'])
    frame_times = np.asarray(behaviour['frame_times'], dtype=float)
    plane = Path(row['registration_plane0'])
    ops = np.load(plane / 'ops.npy', allow_pickle=True).item()
    xoff = np.asarray(ops['xoff'], dtype=float)
    yoff = np.asarray(ops['yoff'], dtype=float)
    max_duration_s = float(row['max_stim_duration_s'])
    max_duration_frames = int(round(max_duration_s * 30))

    dff_red = np.load(
        SESSION_ROOT / recname / 'single_trial_regression' / 'dff_red.npy',
        mmap_mode='r',
    ).reshape(1024, -1)
    whole_red = np.nanmean(dff_red, axis=0)
    ctrl_stim = pd.read_parquet(
        SESSION_ROOT / recname / f'{recname}_per_acr_trial.parquet'
    )
    stim_red = (
        ctrl_stim.loc[ctrl_stim['condition'].eq('stim')]
        .groupby('event_index')['red_delta']
        .median()
    )

    rows = []
    for _, event in events.iterrows():
        event_index = int(event['stimulation_index_zero_based'])
        onset = int(event['onset_frame'])
        if onset >= len(frame_times):
            continue
        event_time = frame_times[onset]
        speed_values = event_speed(behaviour, event_time, max_duration_s)
        motion_values = event_motion(
            xoff, yoff, onset, max_duration_frames
        )
        response_start = onset + max_duration_frames
        response_end = response_start + 30
        whole_red_baseline = mean_window(whole_red, onset - 30, onset)
        whole_red_response = mean_window(
            whole_red, response_start, response_end
        )
        rows.append({
            'rec_id': recname,
            'animal': row['animal'],
            'event_index': event_index,
            'onset_frame': onset,
            'red_delta_grid_median': stim_red.get(event_index, np.nan),
            'red_dff_wholefield_delta': (
                whole_red_response - whole_red_baseline
                if np.isfinite(whole_red_response)
                and np.isfinite(whole_red_baseline) else np.nan
            ),
            **speed_values,
            **motion_values,
        })
    return pd.DataFrame(rows)


#%% run
def main():
    inventory = pd.read_csv(INVENTORY_PATH)
    inventory = inventory.loc[
        inventory['qualifies_ctrl_stim'].astype(bool)
        & inventory['session_lists'].eq('pathdLightLCOpto')
    ].copy()
    events = pd.read_csv(EVENT_PATH)
    events = events.loc[events['qualifies_ctrl_stim'].astype(bool)].copy()
    event_tables = []
    for _, row in inventory.iterrows():
        print(row['session'], flush=True)
        event_tables.append(
            session_rows(
                row,
                events.loc[events['session'].eq(row['session'])],
            )
        )
    event_table = pd.concat(event_tables, ignore_index=True)
    event_table.to_parquet(
        OUTPUT_ROOT / 'red_control_diagnostic_events.parquet',
        index=False,
    )

    summary_rows = []
    ctrl_stim_summary = pd.read_parquet(
        OUTPUT_ROOT / 'ctrl_stim_per_recording.parquet'
    )
    for recname, group in event_table.groupby('rec_id'):
        ctrl_stim = ctrl_stim_summary.loc[ctrl_stim_summary['rec_id'].eq(recname)].iloc[0]
        summary_rows.append({
            'rec_id': recname,
            'animal': group['animal'].iloc[0],
            'red_grid_median': group['red_delta_grid_median'].median(),
            'red_wholefield_median': group['red_dff_wholefield_delta'].median(),
            'speed_change_median': group['speed_change'].median(),
            'motion_change_median': group['motion_change'].median(),
            'motion_step_response_median': group['motion_step_response'].median(),
            'ctrl_stim_red_stim_minus_ctrl': ctrl_stim['stim_minus_ctrl_red_delta_median'],
        })
    summary = pd.DataFrame(summary_rows)
    summary.to_parquet(
        OUTPUT_ROOT / 'red_control_diagnostic_sessions.parquet',
        index=False,
    )

    for x, y in [
            ('red_grid_median', 'speed_change_median'),
            ('red_grid_median', 'motion_change_median'),
            ('red_grid_median', 'motion_step_response_median'),
            ('red_grid_median', 'red_wholefield_median'),
    ]:
        values = summary[[x, y]].dropna()
        result = spearmanr(values[x], values[y]) if len(values) >= 3 else None
        print(
            x, y,
            'n=', len(values),
            'rho=', getattr(result, 'statistic', np.nan),
            'p=', getattr(result, 'pvalue', np.nan),
        )

    for column in [
            'red_grid_median',
            'red_wholefield_median',
            'speed_change_median',
            'motion_change_median',
            'motion_step_response_median',
    ]:
        values = summary[column].dropna()
        print(
            column,
            'median=', values.median(),
            'q25=', values.quantile(.25),
            'q75=', values.quantile(.75),
        )

    fig, axes = plt.subplots(1, 3, figsize=(8.2, 2.6))
    for axis, x, xlabel in zip(
            axes,
            ('speed_change_median', 'motion_change_median',
             'motion_step_response_median'),
            ('speed change', 'registration shift change',
             'registration step response')):
        axis.scatter(
            summary[x],
            summary['red_grid_median'],
            s=14,
            color='tab:red',
            alpha=.7,
        )
        axis.set_xlabel(xlabel)
        axis.set_ylabel('red dF/F change')
        axis.spines[['top', 'right']].set_visible(False)
    fig.tight_layout()
    fig.savefig(FIGURE_ROOT / 'red_control_diagnostic_correlations.png', dpi=300)
    plt.close(fig)


if __name__ == '__main__':
    main()
