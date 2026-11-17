# -*- coding: utf-8 -*-
'''
Created on Fri Nov 13 00:16:16 2026

compare saved single-trial red slopes between ctrl and stim

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

from compare_ctrl_stim_dlight import final_sessions, load_behaviour

repo_root    = Path(__file__).resolve().parents[3]
OUTPUT_ROOT  = pp.DLIGHT_HOTSPOTS_OUTPUT_ROOT
SESSION_ROOT = OUTPUT_ROOT / 'ctrl_stim_sessions'


#%% saved files and coefficient identification
def inspect_folders(sessions):
    saved = {}
    for recname in sessions.session:
        folder = SESSION_ROOT / recname / 'single_trial_regression'
        path = folder / 'regression_qc.npy'
        if not path.exists():
            continue
        qc = pd.read_pickle(path)
        saved[recname] = qc
    return saved


#%% trial mapping and medians
def session_medians(row, qc):
    recname = row.session
    behaviour = load_behaviour(pd.Series(row._asdict()))
    run_onset_frames = np.asarray(behaviour['run_onset_frames'], dtype=float)
    path = SESSION_ROOT / recname / f'{recname}_per_acr_trial.parquet'
    labels = pd.read_parquet(path, columns=['condition', 'trial_index', 'event_index',
                                          'grid_y', 'grid_x'])
    coords = labels[['grid_y', 'grid_x']].drop_duplicates()
    n_saved = len(qc)
    qc = qc.merge(coords, on=['grid_y', 'grid_x'], validate='one_to_one')
    records = []
    for acr in qc.itertuples():
        assert len(acr.b2) == len(acr.event_frame) == len(acr.start_frame) == len(acr.end_frame)
        records.extend((acr.grid_y, acr.grid_x, int(frame), float(beta))
                       for frame, beta in zip(acr.event_frame, acr.b2))
    coefficients = pd.DataFrame(records, columns=['grid_y', 'grid_x', 'run_frame', 'b2'])
    duplicate = coefficients.duplicated(['grid_y', 'grid_x', 'run_frame'], keep=False)
    if duplicate.any():
        coefficients = coefficients.loc[~duplicate]
    ctrl_ids = set(labels.loc[labels.condition.eq('ctrl'), 'trial_index'])
    stim_ids = set(labels.loc[labels.condition.eq('stim'), 'trial_index'])
    overlap = sorted(ctrl_ids & stim_ids)
    result = dict(rec_id=recname, animal=row.animal, n_saved_acr=n_saved, n_acr=len(qc),
                  n_ctrl_stim_overlap_trials=len(overlap),
                  n_ambiguous_coefficient_rows=int(duplicate.sum()), status='ok')
    for condition in ('ctrl', 'stim'):
        subset = labels.loc[labels.condition.eq(condition)]
        trials = subset[['trial_index', 'grid_y', 'grid_x']].drop_duplicates().copy()
        valid = (trials.trial_index >= 0) & (trials.trial_index < len(run_onset_frames))
        trials['run_frame'] = np.nan
        trials.loc[valid, 'run_frame'] = run_onset_frames[trials.loc[valid, 'trial_index'].to_numpy(dtype=int)]
        joined = trials.merge(coefficients, on=['grid_y', 'grid_x', 'run_frame'],
                              how='left', validate='many_to_one', indicator=True)
        finite = np.isfinite(joined.b2)
        values = joined.loc[finite, 'b2'].to_numpy()
        low, median, high = np.percentile(values, [25, 50, 75]) if len(values) else (np.nan, np.nan, np.nan)
        matched = joined.loc[joined['_merge'].eq('both')]
        reused = int(matched.duplicated(['grid_y', 'grid_x', 'run_frame']).sum())
        nonfinite = int((~finite & joined['_merge'].eq('both')).sum())
        result.update({f'{condition}_red_coefficient_median': median,
                       f'{condition}_red_coefficient_q25': low,
                       f'{condition}_red_coefficient_q75': high,
                       f'{condition}_n_events': subset.event_index.nunique(),
                       f'{condition}_n_trials': subset.trial_index.nunique(),
                       f'{condition}_n_trials_with_coefficient': joined.loc[finite, 'trial_index'].nunique(),
                       f'{condition}_n_coefficients': len(values),
                       f'{condition}_n_unmatched_acr_trials': int(joined['_merge'].eq('left_only').sum()),
                       f'{condition}_n_nonfinite_coefficients': nonfinite,
                       f'{condition}_n_reused_run_frame_pairs': reused})
    result['stim_minus_ctrl_median'] = (result['stim_red_coefficient_median']
                                        - result['ctrl_red_coefficient_median'])
    return result


#%% comparison
def main():
    output = OUTPUT_ROOT / 'red_coefficient_ctrl_stim.csv'
    if output.exists():
        raise FileExistsError(f'Existing output left unchanged: {output}')
    sessions = final_sessions()
    folders = {path.name for path in SESSION_ROOT.iterdir() if path.is_dir()}
    assert folders == set(sessions.session)
    saved = inspect_folders(sessions)
    rows = []
    for row in sessions.itertuples():
        if row.session not in saved:
            rows.append(dict(rec_id=row.session, animal=row.animal, status='missing saved coefficients'))
            continue
        result = session_medians(row, saved[row.session])
        rows.append(result)
        print(f'{row.session}: ctrl={result["ctrl_red_coefficient_median"]:.9g}, '
              f'stim={result["stim_red_coefficient_median"]:.9g}', flush=True)
    table = pd.DataFrame(rows)
    table.to_csv(output, index=False)


if __name__ == '__main__':
    main()
