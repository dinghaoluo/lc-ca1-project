# -*- coding: utf-8 -*-
'''
Created on Tue Nov 17 09:18:21 2026

characterise registration, green-to-red coupling and red smoothing in ctrl and stim trials

@author: Dinghao Luo
'''

#%% imports
import sys
from pathlib import Path
from contextlib import redirect_stdout

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
from scipy.stats import spearmanr, wilcoxon
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


def regress_smoothed_red(target, green, red, run_onset_frames, log):
    red = smooth_red(red)
    with redirect_stdout(log):
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
def analyse_session(row, diagnostics, previous_coefficients, log):
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
    log.write(f'{recname}: ctrl registration q95={cutoff:.9g}; kept {len(kept_ids)}/{len(stim_ids)} stim events; '
              f'{len(paired_events)} finite registration/red event pairs; '
              f'missing motion ctrl={len(ctrl_ids)-len(ctrl_events)}, stim={len(stim_ids)-len(stim_events)}.\n')
    log.flush()

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
    corrected, coefficients = regress_smoothed_red(target, green, red, run_onset_frames, log)
    del target, green, red
    smoothed_dff = finite_dff(corrected.reshape(-1, n_frames))
    del corrected
    smoothed_values = comparison.event_values(smoothed_dff, red_dff, coords, specs, max_duration_frames)
    del smoothed_dff
    coefficient_summary = session_medians(row, coefficients, log)
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
def paired_test(table, left, right, label, log):
    selected = table.loc[np.isfinite(table[left])].copy()
    if right is not None:
        selected = selected.loc[np.isfinite(selected[right])]
        differences = selected[left].to_numpy() - selected[right].to_numpy()
    else:
        differences = selected[left].to_numpy()
    if len(differences):
        test = wilcoxon(differences, alternative='two-sided') if np.any(differences != 0) else None
        statistic, pvalue = (float(test.statistic), float(test.pvalue)) if test else (0.0, 1.0)
        q25, median, q75 = np.percentile(differences, [25, 50, 75])
    else:
        statistic, pvalue, q25, median, q75 = (np.nan,) * 5
    log.write(f'{label}: median difference={median:.9g}, IQR [{q25:.9g}, {q75:.9g}]; '
              f'two-sided paired Wilcoxon W={statistic:.9g}, p={pvalue:.12g}, '
              f'n={len(selected)} sessions, n={selected.animal.nunique()} animals.\n')
    return len(selected), pvalue


def write_results(table, log):
    tests = [
        ('stim_registration_step', 'ctrl_registration_step', 'Registration step, stim minus ctrl'),
        ('registration_red_rho', None, 'Within-session stim registration/red Spearman rho against 0'),
        ('filtered_stim_red_change', 'filtered_ctrl_red_change', 'After registration exclusion, red stim minus ctrl'),
        ('filtered_stim_dlight_change', 'filtered_ctrl_dlight_change', 'After registration exclusion, dLight stim minus ctrl'),
        ('stim_red_on_green_slope', 'ctrl_red_on_green_slope', 'Red-on-green OLS slope, stim minus ctrl'),
        ('observed_red_difference', 'predicted_red_difference', 'Observed red difference minus ctrl-slope prediction'),
        ('unsmoothed_stim_red_coef', 'unsmoothed_ctrl_red_coef', 'Unsmoothed red coefficient, stim minus ctrl'),
        ('smoothed_stim_red_coef', 'smoothed_ctrl_red_coef', 'Smoothed red coefficient, stim minus ctrl'),
        ('unsmoothed_dlight_difference', None, 'Unsmoothed dLight stim minus ctrl'),
        ('smoothed_dlight_difference', None, 'Smoothed-red dLight stim minus ctrl')]
    for left, right, label in tests:
        paired_test(table, left, right, label, log)
    log.write('\nPooled session summaries give each session equal weight; IQRs describe dispersion across sessions.\n\n')
    for column in table.select_dtypes(include='number'):
        finite = table.loc[np.isfinite(table[column])]
        if not len(finite):
            continue
        q25, median, q75 = np.percentile(finite[column], [25, 50, 75])
        log.write(f'{column}: median={median:.9g}, IQR [{q25:.9g}, {q75:.9g}], '
                  f'n={len(finite)} sessions, n={finite.animal.nunique()} animals.\n')
    log.write('\nPer-session results:\n\n' + table.to_csv(index=False) + '\n')


def main():
    output = OUTPUT_ROOT / 'red_increase_source.csv'
    if output.exists():
        raise FileExistsError(output)
    sessions = comparison.final_sessions()
    diagnostics = pd.read_parquet(OUTPUT_ROOT / 'red_control_diagnostic_events.parquet')
    previous_coefficients = pd.read_csv(OUTPUT_ROOT / 'red_coefficient_ctrl_stim.csv')
    with (OUTPUT_ROOT / 'PROGRESS.md').open('a', encoding='utf-8') as log:
        log.write('\n## 17 November 2026: sources of the red increase\n\n'
                  'I use the 40 saved ctrl/stim sessions, saved ACR-event tables and saved stim-event diagnostics. '
                  'The combined ctrl/stim diagnostic file is missing, so control rows are recomputed with '
                  'diagnose_red_ctrl_motion.control_rows unchanged and retained in memory only, restricted to saved comparison events. '
                  'The existing registration diagnostic uses onset plus the rounded maximum duration, '
                  'whereas fluorescence uses one additional frame for stim; both existing definitions are retained. '
                  'Session fluorescence is the median across ACRs of each ACR\'s mean over condition events, as in session_summary. '
                  'The motion filter keeps stim steps <= the ctrl 95th percentile; non-finite steps cannot pass. '
                  'Spearman rho pairs stim steps with that event\'s median ACR red change.\n\n'
                  'Raw green is original_dlight minus 0.2 neuropil_dlight before regression, converted with '
                  'percentile_dff_cpu unchanged and evaluated with event_values unchanged. OLS includes an intercept '
                  'and pools finite ACR-event pairs within each condition; prediction is the ctrl slope multiplied '
                  'by the stim-minus-ctrl difference of the median per-ACR condition means of raw green. '
                  'The observed red difference uses the identical aggregation.\n\n'
                  'Smoothed-red regression calls the original shared screening and per-ACR worker with the same '
                  'regressors, windows, 100-frame minimum and subtraction rule, retaining arrays and QC in memory. '
                  'Gaussian sigma=3, nearest mode was unspecified at missing shutter frames: I normalise the '
                  'convolution by finite weights and restore every original NaN, so smoothing does not expand '
                  'the excluded frame mask. The green regressor and target are unchanged. Coefficients are '
                  'mapped to ctrl/stim trials by session_medians unchanged; repeated events on one trial do '
                  'not duplicate its coefficient. Saved unsmoothed coefficients and dLight changes are shown beside it. '
                  'All tests are two-sided Wilcoxon across sessions, zero_method=wilcox and method=auto; '
                  'animal counts accompany every test. No animal-level test is added.\n\n')
        log.flush()
        rows = []
        for row in sessions.itertuples():
            try:
                result = analyse_session(row, diagnostics, previous_coefficients, log)
            except Exception as error:
                log.write(f'{row.session}: unavailable input or failed calculation: {type(error).__name__}: {error}.\n')
                log.flush()
                raise
            rows.append(result)
            print(f'{row.session}: kept={result["fraction_stim_kept"]:.3f}; '
                  f'red difference={result["observed_red_difference"]:.6g}; '
                  f'predicted={result["predicted_red_difference"]:.6g}; '
                  f'smoothed dLight difference={result["smoothed_dlight_difference"]:.6g}', flush=True)
            log.write(f'{row.session}: completed all red-source calculations.\n')
            log.flush()
        table = pd.DataFrame(rows)
        table.to_csv(output, index=False)
        write_results(table, log)
        log.write('\nTable: outputs/da_hotspots/red_increase_source.csv. No arrays, regression folders or earlier numerical outputs were written.\n')
        print(f'Completed {len(table)} sessions from {table.animal.nunique()} animals.', flush=True)


#%% run
if __name__ == '__main__':
    with threadpool_limits(limits=1):
        main()
