# -*- coding: utf-8 -*-
'''
Created on Sun Nov  8 10:48:00 2026

compare red, movement and speed signals in ctrl and stimulation trials

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
from scipy.stats import wilcoxon

repo_root = Path(__file__).resolve().parents[3]
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))
sys.path.insert(0, str(Path(__file__).resolve().parent))

import diagnose_red_control as diag


#%% ctrl extraction
def control_rows(row, ctrl_stim):
    recname = row['session']
    behaviour = diag.load_behaviour(row['behaviour_path'])
    frame_times = np.asarray(behaviour['frame_times'], dtype=float)
    run_onsets = np.asarray(behaviour['run_onset_frames'], dtype=float)
    methods = [str(statement[15]) for statement in behaviour['trial_statements']]
    plane = Path(row['registration_plane0'])
    ops = np.load(plane / 'ops.npy', allow_pickle=True).item()
    xoff = np.asarray(ops['xoff'], dtype=float)
    yoff = np.asarray(ops['yoff'], dtype=float)
    dff_red = np.load(
        diag.SESSION_ROOT / recname / 'single_trial_regression' / 'dff_red.npy',
        mmap_mode='r',
    ).reshape(1024, -1)
    whole_red = np.nanmean(dff_red, axis=0)
    ctrl_red = (
        ctrl_stim.loc[ctrl_stim['condition'].eq('ctrl')]
        .groupby('event_index')['red_delta']
        .median()
    )
    rows = []
    for trial, onset in enumerate(run_onsets):
        if (
                trial >= len(methods)
                or methods[trial] != '0'
                or not np.isfinite(onset)
        ):
            continue
        onset = int(onset)
        if onset < 30 or onset + 30 > len(frame_times):
            continue
        event_time = frame_times[onset]
        speed_values = diag.event_speed(behaviour, event_time, 0)
        motion_values = diag.event_motion(xoff, yoff, onset, 0)
        whole_baseline = diag.mean_window(whole_red, onset - 30, onset)
        whole_response = diag.mean_window(whole_red, onset, onset + 30)
        rows.append({
            'rec_id': recname,
            'animal': row['animal'],
            'condition': 'ctrl',
            'event_index': trial,
            'red_delta_grid_median': ctrl_red.get(trial, np.nan),
            'red_dff_wholefield_delta': (
                whole_response - whole_baseline
                if np.isfinite(whole_response)
                and np.isfinite(whole_baseline) else np.nan
            ),
            **speed_values,
            **motion_values,
        })
    return pd.DataFrame(rows)


#%% run
def main():
    inventory = pd.read_csv(diag.INVENTORY_PATH)
    inventory = inventory.loc[
        inventory['qualifies_ctrl_stim'].astype(bool)
        & inventory['session_lists'].eq('pathdLightLCOpto')
    ].copy()
    events = pd.read_csv(diag.EVENT_PATH)
    ctrl_stim = pd.read_parquet(
        diag.OUTPUT_ROOT / 'ctrl_stim_per_acr_trial.parquet'
    )
    rows = []
    for _, row in inventory.iterrows():
        recname = row['session']
        stim = diag.session_rows(
            row, events.loc[events['session'].eq(recname)]
        )
        stim['condition'] = 'stim'
        ctrl = control_rows(
            row,
            ctrl_stim.loc[ctrl_stim['rec_id'].eq(recname)],
        )
        rows.extend([stim, ctrl])
    table = pd.concat(rows, ignore_index=True)
    table.to_parquet(
        diag.OUTPUT_ROOT / 'red_control_diagnostic_ctrl_stim_events.parquet',
        index=False,
    )
    summary = (
        table.groupby(['rec_id', 'animal', 'condition'], as_index=False)
        .median(numeric_only=True)
    )
    summary.to_parquet(
        diag.OUTPUT_ROOT / 'red_control_diagnostic_ctrl_stim_sessions.parquet',
        index=False,
    )
    wide = summary.pivot(index=['rec_id', 'animal'], columns='condition')
    for metric in [
            'red_delta_grid_median',
            'red_dff_wholefield_delta',
            'speed_change',
            'motion_change',
            'motion_step_response',
    ]:
        ctrl = wide[(metric, 'ctrl')]
        stim = wide[(metric, 'stim')]
        keep = ctrl.notna() & stim.notna()
        result = wilcoxon(stim[keep], ctrl[keep])
        print(
            metric,
            'n=', int(keep.sum()),
            'ctrl=', float(ctrl[keep].median()),
            'stim=', float(stim[keep].median()),
            'W=', float(result.statistic),
            'p=', float(result.pvalue),
        )


if __name__ == '__main__':
    main()
