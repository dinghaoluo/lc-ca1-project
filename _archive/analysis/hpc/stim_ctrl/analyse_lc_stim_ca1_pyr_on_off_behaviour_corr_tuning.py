# -*- coding: utf-8 -*-
'''
temporary amplitude-window tuning for HPCLC PyrUp/PyrDown behaviour correlations.

@author: Dinghao Luo
'''

#%% imports
from pathlib import Path
import argparse
import pickle
import sys

import numpy as np
import pandas as pd
from scipy.stats import linregress

repo_root = Path(__file__).resolve().parent
while not (repo_root / 'utils').exists() and repo_root != repo_root.parent:
    repo_root = repo_root.parent
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))

import project_paths as pp
import rec_list


#%% command-line options
parser = argparse.ArgumentParser(
    description='Tune HPCLC pyramidal ON/OFF behaviour-correlation windows.')
parser.add_argument(
    '--recompute',
    action='store_true',
    help='recompute session mean profiles instead of reusing the saved tuning cache',
)
args = parser.parse_args()


#%% parameters
SAMP_FREQ   = 1250
RUN_ONSET   = 3 * SAMP_FREQ
RUN_ONSET_ON  = 'run-onset ON'
RUN_ONSET_OFF = 'run-onset OFF'

WINDOW_SPECS = [
    ('0.5-1.5 s', 0.50, 1.50),
    ('0.5-2.5 s', 0.50, 2.50),
    ('0.5-3.5 s', 0.50, 3.50),
    ('0.75-1.25 s', 0.75, 1.25),
]

ON_TIME_AMP_BOUNDS    = (0.0, 2.0)
ON_TIME_UPPER_BOUNDS  = (None, 2.0)
ON_EXTRA_LOWER_BOUNDS = (-0.5,)
OFF_TIME_AMP_BOUNDS   = (-2.0, 0.5)
OFF_TIME_DELTA_BOUNDS = (-400.0, 1000.0)

all_sess_stem = pp.HPC_EPHYS_STEM / 'all_sessions'
all_beh_stem  = pp.BEHAVIOUR_EXPERIMENTS_STEM
output_stem   = repo_root / 'data' / 'analysis' / 'hpc' / 'stim_ctrl'
output_stem.mkdir(parents=True, exist_ok=True)

csv_output_path = output_stem / 'lc_stim_ca1_pyr_on_off_amp_window_tuning.csv'
pkl_output_path = output_stem / 'lc_stim_ca1_pyr_on_off_amp_window_tuning.pkl'
loo_csv_output_path = output_stem / 'lc_stim_ca1_pyr_on_off_unbounded_ON_drop_one_session.csv'
loo_pkl_output_path = output_stem / 'lc_stim_ca1_pyr_on_off_unbounded_ON_drop_one_session.pkl'


#%% helpers
def _window_slice(start_s, end_s):
    start = RUN_ONSET + int(round(start_s * SAMP_FREQ))
    end = RUN_ONSET + int(round(end_s * SAMP_FREQ))
    return slice(start, end)


def _load_behaviour(recname):
    beh_path = all_beh_stem / 'HPCLC' / f'{recname}.pkl'
    if not beh_path.exists():
        beh_path = all_beh_stem / 'HPCLCterm' / f'{recname}.pkl'

    with open(beh_path, 'rb') as f:
        return pickle.load(f)


def _behaviour_deltas(beh):
    first_lick_times = [t[0][0] - s for t, s
                        in zip(beh['lick_times'], beh['run_onsets']) if t]
    stim_conds = [t[15] for t in beh['trial_statements']][1:]
    stim_idx = [trial for trial, cond in enumerate(stim_conds) if cond != '0']
    ctrl_idx = [trial+2 for trial in stim_idx]

    lick_time_delta = np.median(
        [t for i, t in enumerate(first_lick_times) if i in stim_idx]
        ) - np.median(
        [t for i, t in enumerate(first_lick_times) if i in ctrl_idx]
        )

    return stim_idx, ctrl_idx, lick_time_delta


def _mean_profile_or_none(traces):
    if not traces:
        return None
    return np.mean(traces, axis=0)


def _profile_window_delta(ctrl_profile, stim_profile, amp_window):
    if ctrl_profile is None or stim_profile is None:
        return np.nan

    return np.mean(stim_profile[amp_window]) - np.mean(ctrl_profile[amp_window])


def _within_bounds(value, bounds):
    if value is None or np.isnan(value):
        return False
    low, high = bounds
    if low is not None and value < low:
        return False
    if high is not None and value > high:
        return False
    return True


def _fit_records(records, amp_key, amp_bounds, delta_bounds=(None, None)):
    filtered = [
        rec for rec in records
        if _within_bounds(rec[amp_key], amp_bounds)
        and _within_bounds(rec['lick_time_delta'], delta_bounds)
    ]

    if len(filtered) < 2:
        return {
            'n': len(filtered),
            'slope': np.nan,
            'intercept': np.nan,
            'r': np.nan,
            'p': np.nan,
        }

    x_values = np.array([rec[amp_key] for rec in filtered], dtype=float)
    y_values = np.array([rec['lick_time_delta'] for rec in filtered], dtype=float)
    slope, intercept, r, p, _ = linregress(x_values, y_values)

    return {
        'n': len(filtered),
        'slope': slope,
        'intercept': intercept,
        'r': r,
        'p': p,
    }


def _add_result(rows, window_label, start_s, end_s, analysis_label,
                records, amp_key, amp_bounds, delta_bounds):
    fit = _fit_records(records, amp_key, amp_bounds, delta_bounds)
    rows.append({
        'window_label': window_label,
        'window_start_s': start_s,
        'window_end_s': end_s,
        'analysis_label': analysis_label,
        'amp_key': amp_key,
        'amp_bound_low': amp_bounds[0],
        'amp_bound_high': amp_bounds[1],
        'delta_bound_low': delta_bounds[0],
        'delta_bound_high': delta_bounds[1],
        **fit,
        'abs_r': abs(fit['r']) if np.isfinite(fit['r']) else np.nan,
    })


def _add_baseline_deltas(result_df):
    result_df = result_df.copy()
    result_df['delta_r_from_0.5-1.5_s'] = np.nan
    result_df['delta_abs_r_from_0.5-1.5_s'] = np.nan

    for analysis_label in result_df['analysis_label'].unique():
        baseline = result_df[
            (result_df['analysis_label'] == analysis_label)
            & (result_df['window_label'] == '0.5-1.5 s')
        ]
        if baseline.empty:
            continue

        base_r = baseline.iloc[0]['r']
        base_abs_r = baseline.iloc[0]['abs_r']
        idx = result_df['analysis_label'] == analysis_label
        result_df.loc[idx, 'delta_r_from_0.5-1.5_s'] = result_df.loc[idx, 'r'] - base_r
        result_df.loc[idx, 'delta_abs_r_from_0.5-1.5_s'] = (
            result_df.loc[idx, 'abs_r'] - base_abs_r
        )

    return result_df


def _integer_lower_bound_sweep(records, amp_key, upper_bound):
    values = np.array([
        rec[amp_key] for rec in records
        if rec[amp_key] is not None
        and np.isfinite(rec[amp_key])
        and rec[amp_key] <= upper_bound
    ], dtype=float)

    if values.size == 0:
        return []

    lowest_integer = int(np.floor(np.min(values)))
    lower_bounds = []
    for lower_bound in ON_EXTRA_LOWER_BOUNDS:
        if lower_bound < 0 and np.min(values) <= lower_bound:
            lower_bounds.append(float(lower_bound))

    if lowest_integer < 0:
        lower_bounds.extend(
            float(lower) for lower in range(-1, lowest_integer - 1, -1))

    return sorted(set(lower_bounds), reverse=True)


def _load_cached_tuning_records():
    if args.recompute or not pkl_output_path.exists():
        return None, None

    with open(pkl_output_path, 'rb') as f:
        payload = pickle.load(f)

    cached_records = payload.get('base_session_records')
    cached_window_records = payload.get('window_records')

    if cached_window_records:
        print(f'Reusing cached window records from: {pkl_output_path}')
    elif cached_records:
        print(f'Reusing cached session profiles from: {pkl_output_path}')

    return cached_records, cached_window_records


def _drop_one_session_results(window_records):
    rows = []

    for window_label, records in window_records.items():
        full_fit = _fit_records(
            records,
            'amp_remain_ON_delta_mean',
            (None, None),
            (None, None),
        )

        for dropped in records:
            kept_records = [
                rec for rec in records
                if rec['recname'] != dropped['recname']
            ]
            drop_fit = _fit_records(
                kept_records,
                'amp_remain_ON_delta_mean',
                (None, None),
                (None, None),
            )

            full_abs_r = abs(full_fit['r']) if np.isfinite(full_fit['r']) else np.nan
            drop_abs_r = abs(drop_fit['r']) if np.isfinite(drop_fit['r']) else np.nan

            rows.append({
                'window_label': window_label,
                'dropped_recname': dropped['recname'],
                'full_n': full_fit['n'],
                'full_r': full_fit['r'],
                'full_p': full_fit['p'],
                'drop_n': drop_fit['n'],
                'drop_r': drop_fit['r'],
                'drop_p': drop_fit['p'],
                'delta_r_from_full': drop_fit['r'] - full_fit['r'],
                'delta_abs_r_from_full': drop_abs_r - full_abs_r,
                'dropped_amp_remain_ON_delta_mean': dropped['amp_remain_ON_delta_mean'],
                'dropped_lick_time_delta': dropped['lick_time_delta'],
                'dropped_n_remain_ON': dropped['n_remain_ON'],
            })

    return pd.DataFrame(rows)


#%% load data
base_session_records, cached_window_records = _load_cached_tuning_records()


#%% tuning
window_records = cached_window_records if cached_window_records is not None else {}
rows = []

if base_session_records is None and cached_window_records is None:
    print('Loading dataframes...')
    cell_profiles = pd.read_pickle(pp.HPC_EPHYS_STEM / 'hpc_all_profiles.pkl')
    df_pyr = cell_profiles[cell_profiles['cell_identity'] == 'pyr']

    base_session_records = []
    print('Loading per-session ON/OFF profiles...')
    for path in rec_list.pathHPCLCopt:
        recname = Path(path).name
        print(recname)

        train_path = all_sess_stem / recname / f'{recname}_all_trains_run.npy'
        trains = np.load(train_path, allow_pickle=True).item()

        beh = _load_behaviour(recname)
        stim_idx, ctrl_idx, ctrl_stim_lick_time_delta = _behaviour_deltas(beh)

        curr_df_pyr = df_pyr[df_pyr['recname'] == recname]
        sess_ctrl_remain_ON, sess_stim_remain_ON = [], []
        sess_ctrl_remain_OFF, sess_stim_remain_OFF = [], []

        for cluname, session in curr_df_pyr.iterrows():
            if session['class_ctrl'] == RUN_ONSET_ON and session['class_stim'] == RUN_ONSET_ON:
                sess_ctrl_remain_ON.append(np.mean(trains[cluname][ctrl_idx], axis=0))
                sess_stim_remain_ON.append(np.mean(trains[cluname][stim_idx], axis=0))

            if session['class_ctrl'] == RUN_ONSET_OFF and session['class_stim'] == RUN_ONSET_OFF:
                sess_ctrl_remain_OFF.append(np.mean(trains[cluname][ctrl_idx], axis=0))
                sess_stim_remain_OFF.append(np.mean(trains[cluname][stim_idx], axis=0))

        base_session_records.append({
            'recname': recname,
            'lick_time_delta': ctrl_stim_lick_time_delta,
            'ctrl_remain_ON_profile': _mean_profile_or_none(sess_ctrl_remain_ON),
            'stim_remain_ON_profile': _mean_profile_or_none(sess_stim_remain_ON),
            'ctrl_remain_OFF_profile': _mean_profile_or_none(sess_ctrl_remain_OFF),
            'stim_remain_OFF_profile': _mean_profile_or_none(sess_stim_remain_OFF),
            'n_remain_ON': len(sess_ctrl_remain_ON),
            'n_remain_OFF': len(sess_ctrl_remain_OFF),
        })

for window_label, start_s, end_s in WINDOW_SPECS:
    print(f'\nWindow: {window_label}')

    if cached_window_records is None:
        amp_window = _window_slice(start_s, end_s)
        session_delta_records = []

        for session_record in base_session_records:
            session_delta_records.append({
                'recname': session_record['recname'],
                'window_label': window_label,
                'window_start_s': start_s,
                'window_end_s': end_s,
                'amp_remain_ON_delta_mean': _profile_window_delta(
                    session_record['ctrl_remain_ON_profile'],
                    session_record['stim_remain_ON_profile'],
                    amp_window),
                'amp_remain_OFF_delta_mean': _profile_window_delta(
                    session_record['ctrl_remain_OFF_profile'],
                    session_record['stim_remain_OFF_profile'],
                    amp_window),
                'lick_time_delta': session_record['lick_time_delta'],
                'n_remain_ON': session_record['n_remain_ON'],
                'n_remain_OFF': session_record['n_remain_OFF'],
            })

        window_records[window_label] = session_delta_records
    else:
        session_delta_records = window_records[window_label]

    _add_result(
        rows, window_label, start_s, end_s,
        'ON_all_finite',
        session_delta_records,
        'amp_remain_ON_delta_mean',
        (None, None),
        (None, None),
    )
    _add_result(
        rows, window_label, start_s, end_s,
        'ON_current_0_to_2_Hz_bound',
        session_delta_records,
        'amp_remain_ON_delta_mean',
        ON_TIME_AMP_BOUNDS,
        (None, None),
    )
    _add_result(
        rows, window_label, start_s, end_s,
        'ON_upper_le_2_Hz_bound',
        session_delta_records,
        'amp_remain_ON_delta_mean',
        ON_TIME_UPPER_BOUNDS,
        (None, None),
    )
    for lower_bound in _integer_lower_bound_sweep(
            session_delta_records, 'amp_remain_ON_delta_mean', 2.0):
        _add_result(
            rows, window_label, start_s, end_s,
            f'ON_ge_{lower_bound:g}_le_2_Hz_bound',
            session_delta_records,
            'amp_remain_ON_delta_mean',
            (lower_bound, 2.0),
            (None, None),
        )
    _add_result(
        rows, window_label, start_s, end_s,
        'OFF_current_bounds',
        session_delta_records,
        'amp_remain_OFF_delta_mean',
        OFF_TIME_AMP_BOUNDS,
        OFF_TIME_DELTA_BOUNDS,
    )


result_df = _add_baseline_deltas(pd.DataFrame(rows))
loo_result_df = _drop_one_session_results(window_records)
result_df.to_csv(csv_output_path, index=False)
loo_result_df.to_csv(loo_csv_output_path, index=False)

with open(pkl_output_path, 'wb') as f:
    pickle.dump({
        'base_session_records': base_session_records,
        'window_records': window_records,
        'result_table': result_df,
        'window_specs': WINDOW_SPECS,
        'sample_frequency_hz': SAMP_FREQ,
        'run_onset_sample': RUN_ONSET,
    }, f)

with open(loo_pkl_output_path, 'wb') as f:
    pickle.dump({
        'drop_one_session_table': loo_result_df,
        'window_records': window_records,
        'analysis_label': 'ON_all_finite',
        'amp_bounds': (None, None),
        'delta_bounds': (None, None),
    }, f)

print('\n-----------\nSTATISTICS\n-----------')
print('amplitude-window tuning:')
print(result_df[[
    'window_label',
    'analysis_label',
    'n',
    'r',
    'p',
    'delta_r_from_0.5-1.5_s',
    'delta_abs_r_from_0.5-1.5_s',
]].to_string(index=False))

print('\nunbounded ON drop-one-session-out: largest positive r shifts:')
print(loo_result_df.sort_values(
    'delta_r_from_full',
    ascending=False,
)[[
    'window_label',
    'dropped_recname',
    'full_r',
    'drop_r',
    'drop_p',
    'delta_r_from_full',
    'delta_abs_r_from_full',
]].head(12).to_string(index=False))

print('\nunbounded ON drop-one-session-out: largest |r| shifts:')
print(loo_result_df.sort_values(
    'delta_abs_r_from_full',
    ascending=False,
)[[
    'window_label',
    'dropped_recname',
    'full_r',
    'drop_r',
    'drop_p',
    'delta_r_from_full',
    'delta_abs_r_from_full',
]].head(12).to_string(index=False))

print(f'\nSaved tuning table to: {csv_output_path}')
print(f'Saved tuning data to: {pkl_output_path}')
print(f'Saved drop-one-session table to: {loo_csv_output_path}')
print(f'Saved drop-one-session data to: {loo_pkl_output_path}')
