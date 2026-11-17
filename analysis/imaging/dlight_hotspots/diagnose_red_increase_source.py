# -*- coding: utf-8 -*-
'''
Created on Tue Nov 17 09:18:21 2026

characterise registration, green-to-red coupling and red smoothing in ctrl and stim trials

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

import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d
from scipy.stats import spearmanr
from threadpoolctl import threadpool_limits

import compare_ctrl_stim_dlight as comparison
import diagnose_red_ctrl_motion as motion
from compare_red_coefficients import session_medians
from dlight_imaging.regression.utils_regression import (
    _prepare_shared_regression_screening,
    _single_trial_regression_worker_nan_aware)

OUTPUT_ROOT = pp.DLIGHT_HOTSPOTS_OUTPUT_ROOT
SESSION_ROOT = OUTPUT_ROOT / 'ctrl_stim_sessions'


#%% event summaries
def acr_median(table, condition, column):
    values = table.loc[table.condition.eq(condition)].groupby(
        ['grid_y', 'grid_x'])[column].mean()
    return float(values.median())


def event_specs(table):
    columns = ['rec_id', 'animal', 'condition', 'event_index', 'trial_index',
               'onset_frame', 'offset_frame', 'duration_frames']
    return table[columns].drop_duplicates().to_dict('records')


def finite_dff(traces):
    valid = np.isfinite(traces).any(axis=1)
    dff = np.full(traces.shape, np.nan, dtype=np.float32)
    dff[valid] = comparison.percentile_dff_cpu(traces[valid])
    return dff


def green_red_slope(table, condition):
    values = table.loc[table.condition.eq(condition),
                       ['raw_green_change', 'red_delta']].to_numpy(dtype=float)
    values = values[np.isfinite(values).all(axis=1)]
    if len(values) < 2 or np.ptp(values[:, 0]) == 0:
        return np.nan, len(values)
    design = np.column_stack((values[:, 0], np.ones(len(values))))
    slope = np.linalg.lstsq(design, values[:, 1], rcond=None)[0][0]
    return float(slope), len(values)


#%% smoothed-red regression
def smooth_red(red):
    finite = np.isfinite(red)
    values = gaussian_filter1d(np.where(finite, red, 0.0), sigma=3,
                               axis=-1, mode='nearest')
    weights = gaussian_filter1d(finite.astype(float), sigma=3,
                                axis=-1, mode='nearest')
    smoothed = np.divide(values, weights, out=np.full(red.shape, np.nan),
                          where=weights > 0)
    smoothed[~finite] = np.nan
    return smoothed


def regress_smoothed_red(target, green, red, run_onset_frames):
    red = smooth_red(red)
    valid_acr, finite_frames, windows = _prepare_shared_regression_screening(
        target, green, red, run_onset_frames, comparison.PRE_FRAMES,
        comparison.POST_FRAMES, comparison.MIN_FINITE_FRAMES)
    corrected = np.full(target.shape, np.nan)
    records = []
    for index in np.flatnonzero(valid_acr.ravel()):
        result = _single_trial_regression_worker_nan_aware(
            int(index), target, green, red, windows, valid_acr, finite_frames)
        if result is not None:
            trace, coefficients = result
            corrected.reshape(-1, corrected.shape[-1])[index] = trace
            records.append(coefficients)
    return corrected, pd.DataFrame(records)


#%% recording diagnostics
def analyse_session(row, diagnostics, previous_coefficients):
    recname = row.session
    folder = SESSION_ROOT / recname
    table = pd.read_parquet(folder / f'{recname}_per_acr_trial.parquet')
    behaviour = comparison.load_behaviour(pd.Series(row._asdict()))
    result = dict(rec_id=recname, animal=row.animal, status='ok')
    ctrl_events = motion.control_rows(pd.Series(row._asdict()), table)
    stim_events = diagnostics.loc[diagnostics.rec_id.eq(recname)].copy()
    stim_events['condition'] = 'stim'
    ctrl_ids = table.loc[table.condition.eq('ctrl'), 'event_index'].unique()
    stim_ids = table.loc[table.condition.eq('stim'), 'event_index'].unique()
    ctrl_events = ctrl_events.loc[ctrl_events.event_index.isin(ctrl_ids)]
    stim_events = stim_events.loc[stim_events.event_index.isin(stim_ids)]
    assert not ctrl_events.event_index.duplicated().any()
    assert not stim_events.event_index.duplicated().any()
    result.update(n_ctrl_events=len(ctrl_ids), n_stim_events=len(stim_ids),
                  n_ctrl_motion=int(ctrl_events.motion_step_response.notna().sum()),
                  n_stim_motion=int(stim_events.motion_step_response.notna().sum()))
    for condition, events in (('ctrl', ctrl_events), ('stim', stim_events)):
        result[f'{condition}_registration_step'] = events.motion_step_response.median()
    stim_red = table.loc[table.condition.eq('stim')].groupby('event_index').red_delta.median()
    paired_events = stim_events[['event_index', 'motion_step_response']].merge(
        stim_red.rename('red_change'), on='event_index', validate='one_to_one').dropna()
    result['n_registration_red_pairs'] = len(paired_events)
    result['registration_red_rho'] = np.nan
    if (len(paired_events) >= 3 and paired_events.motion_step_response.nunique() > 1
            and paired_events.red_change.nunique() > 1):
        result['registration_red_rho'] = float(spearmanr(
            paired_events.motion_step_response, paired_events.red_change)[0])
    cutoff = ctrl_events.motion_step_response.quantile(.95)
    kept_ids = stim_events.loc[stim_events.motion_step_response.le(cutoff), 'event_index']
    filtered = table.loc[table.condition.eq('ctrl') | table.event_index.isin(kept_ids)
                         & table.condition.eq('stim')]
    result.update(ctrl_registration_q95=cutoff, n_stim_kept=len(kept_ids),
                  fraction_stim_kept=len(kept_ids) / len(stim_ids) if len(stim_ids) else np.nan)
    for condition in ('ctrl', 'stim'):
        for label, column in (('red', 'red_delta'), ('dlight', 'delta')):
            result[f'filtered_{condition}_{label}_change'] = acr_median(filtered, condition, column)
            result[f'unsmoothed_{condition}_{label}_change'] = acr_median(table, condition, column)

    with np.load(folder / f'{recname}_raw_traces.npz') as raw:
        target = raw['original_dlight'] - .2 * raw['neuropil_dlight']
        green = raw['dlight_regressor']
        red = raw['red_trace']
    assert target.shape == green.shape == red.shape
    n_frames = target.shape[-1]
    raw_green_dff = finite_dff(target.reshape(-1, n_frames))
    specs = event_specs(table)
    coords = table[['grid_y', 'grid_x']].drop_duplicates().to_numpy(dtype=int)
    max_duration_frames = int(round(float(row.max_stim_duration_s) * comparison.SAMP_FREQ))
    red_dff = np.load(folder / 'single_trial_regression' / 'dff_red.npy', mmap_mode='r').reshape(-1, n_frames)
    raw_values = comparison.event_values(raw_green_dff, red_dff, coords, specs, max_duration_frames)
    raw_values = raw_values.rename(columns={'delta': 'raw_green_change'})
    del raw_green_dff
    for condition in ('ctrl', 'stim'):
        slope, count = green_red_slope(raw_values, condition)
        result[f'{condition}_red_on_green_slope'] = slope
        result[f'{condition}_n_green_red_pairs'] = count
        result[f'{condition}_raw_green_change'] = acr_median(raw_values, condition, 'raw_green_change')
    result['observed_red_difference'] = (result['unsmoothed_stim_red_change']
                                          - result['unsmoothed_ctrl_red_change'])
    result['raw_green_difference'] = (result['stim_raw_green_change'] - result['ctrl_raw_green_change'])
    result['predicted_red_difference'] = result['ctrl_red_on_green_slope'] * result['raw_green_difference']
    del raw_values

    run_onset_frames = np.asarray([onset for _, onset, _ in
                                    comparison.valid_run_onsets(behaviour, n_frames)], dtype=int)
    corrected, coefficients = regress_smoothed_red(target, green, red, run_onset_frames)
    del target, green, red
    smoothed_dff = finite_dff(corrected.reshape(-1, n_frames))
    del corrected
    smoothed_values = comparison.event_values(smoothed_dff, red_dff, coords, specs, max_duration_frames)
    del smoothed_dff
    coefficient_summary = session_medians(row, coefficients)
    previous = previous_coefficients.loc[previous_coefficients.rec_id.eq(recname)].iloc[0]
    for condition in ('ctrl', 'stim'):
        result[f'smoothed_{condition}_red_coef'] = coefficient_summary[f'{condition}_red_coefficient_median']
        result[f'unsmoothed_{condition}_red_coef'] = previous[f'{condition}_red_coefficient_median']
        result[f'smoothed_{condition}_dlight_change'] = acr_median(smoothed_values, condition, 'delta')
    for variant in ('smoothed', 'unsmoothed'):
        result[f'{variant}_dlight_difference'] = (result[f'{variant}_stim_dlight_change']
                                                  - result[f'{variant}_ctrl_dlight_change'])
    return result


#%% tests and progress


def main():
    output = OUTPUT_ROOT / 'red_increase_source.csv'
    if output.exists():
        raise FileExistsError(output)
    sessions = comparison.final_sessions()
    diagnostics = pd.read_parquet(OUTPUT_ROOT / 'red_control_diagnostic_events.parquet')
    previous_coefficients = pd.read_csv(OUTPUT_ROOT / 'red_coefficient_ctrl_stim.csv')
    rows = []
    for row in sessions.itertuples():
        result = analyse_session(row, diagnostics, previous_coefficients)
        rows.append(result)
        print(f'{row.session}: kept={result["fraction_stim_kept"]:.3f}; '
              f'red difference={result["observed_red_difference"]:.6g}; '
              f'predicted={result["predicted_red_difference"]:.6g}; '
              f'smoothed dLight difference={result["smoothed_dlight_difference"]:.6g}', flush=True)
    table = pd.DataFrame(rows)
    table.to_csv(output, index=False)
    print(f'Completed {len(table)} sessions from {table.animal.nunique()} animals.', flush=True)


#%% run
if __name__ == '__main__':
    with threadpool_limits(limits=1):
        main()
