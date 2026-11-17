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
from scipy.stats import wilcoxon

from compare_ctrl_stim_dlight import final_sessions, load_behaviour

repo_root    = Path(__file__).resolve().parents[3]
OUTPUT_ROOT  = pp.DLIGHT_HOTSPOTS_OUTPUT_ROOT
SESSION_ROOT = OUTPUT_ROOT / 'ctrl_stim_sessions'
LOG_PATH     = OUTPUT_ROOT / 'PROGRESS.md'


#%% saved files and coefficient identification
def inspect_folders(sessions, log):
    saved = {}
    for recname in sessions.session:
        folder = SESSION_ROOT / recname / 'single_trial_regression'
        entries = []
        for path in sorted(folder.rglob('*')):
            if not path.is_file():
                continue
            entry = f'{path.relative_to(folder)} ({path.stat().st_size} bytes)'
            if path.suffix == '.npy' and path.name != 'regression_qc.npy':
                values = np.load(path, mmap_mode='r')
                entry += f', shape {values.shape}, dtype {values.dtype}'
            entries.append(entry)
        log.write(f'{recname}: ' + '; '.join(entries) + '.\n')
        path = folder / 'regression_qc.npy'
        if not path.exists():
            log.write(f'{recname}: QC missing; no coefficient inferred from traces.\n')
            continue
        qc = pd.read_pickle(path)
        log.write(f'{recname}: pickled dataframe, {len(qc)} ACR rows, columns '
                  f'{list(qc.columns)}, b2 list lengths '
                  f'{sorted(qc.b2.map(len).unique().tolist())}.\n')
        saved[recname] = qc
    log.write('\nCoefficient identification before summarising: in '
              'Z:/Jingyu/code_mpfi_Jingyu/dlight_imaging/regression/utils_regression.py, '
              '_single_trial_regression_worker_nan_aware constructs '
              'column_stack((green, red)), fits LinearRegression with its intercept, '
              'assigns b2=model.coef_[1], and appends b2 alongside event_frame, '
              'start_frame and end_frame. single_trial_regression_parallel uses '
              'to_pickle to save regression_qc.npy, despite the extension. Thus '
              'b2 is the unstandardised tdTomato slope, in target green fluorescence '
              'units per red fluorescence unit; b1 is the green-regressor slope '
              'and r2 is fit quality. compare_ctrl_stim_dlight.py supplied '
              'original_dlight minus 0.2 times neuropil_dlight as the target. '
              'Fit windows extend from the previous run event plus 120 frames '
              '(frame 0 for the first event) to the current run onset plus 120 '
              'frames, using shared finite frames. This slope was not fitted '
              'only within the run-onset response window. No regression is rerun.\n\n')
    log.flush()
    return saved


#%% trial mapping and medians
def session_medians(row, qc, log):
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
        log.write(f'{recname}: {int(duplicate.sum())} ambiguous duplicate ACR/run-frame '
                  'coefficient rows excluded rather than selecting a fit.\n')
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
        missing = sorted(joined.loc[joined['_merge'].eq('left_only'), 'trial_index'].unique().tolist())
        log.write(f'{recname}, {condition}: {subset.event_index.nunique()} events, '
                  f'{subset.trial_index.nunique()} distinct trials, {len(values)} '
                  f'finite coefficients; unmatched trial indices (zero-based) {missing}; '
                  f'reused ACR/run-frame coefficients from tied behaviour onsets {reused}; '
                  f'non-finite matched coefficients {nonfinite}.\n')
    result['stim_minus_ctrl_median'] = (result['stim_red_coefficient_median']
                                        - result['ctrl_red_coefficient_median'])
    log.write(f'{recname}: trial indices labelled both ctrl and stim (zero-based) {overlap}; '
              f'comparison ACRs with QC {len(qc)} / {len(coords)}, saved QC ACRs {n_saved}.\n')
    log.flush()
    return result


#%% comparison and progress
def main():
    output = OUTPUT_ROOT / 'red_coefficient_ctrl_stim.csv'
    if output.exists():
        raise FileExistsError(f'Existing output left unchanged: {output}')
    sessions = final_sessions()
    folders = {path.name for path in SESSION_ROOT.iterdir() if path.is_dir()}
    assert folders == set(sessions.session)
    with LOG_PATH.open('a', encoding='utf-8') as log:
        log.write('\n## 13 November 2026: saved red coefficients, ctrl versus stim\n\n'
                  'I inspected all 40 single_trial_regression folders before using '
                  'their coefficients; the contents and identification follow.\n\n')
        saved = inspect_folders(sessions, log)
        log.write('Trial labels come from each saved per_acr_trial.parquet written '
                  'by compare_ctrl_stim_dlight.py: ctrl is its method-0 valid run-onset '
                  'trial, and stim is its detected PMT event assigned to '
                  'preceding_trial_index_zero_based. Stim events are not relabelled '
                  'by behavioural method. The requested unit is a trial, so repeated '
                  'events assigned to one trial contribute its coefficient once per '
                  'ACR and condition. Trials labelled both ctrl and stim retain both '
                  'original labels, with overlaps logged. Behaviour run-onset frames '
                  'map exactly to QC event_frame, without nearest-frame matching or '
                  'using stimulation-onset frames as regression keys. Only ACRs in '
                  'the saved comparison table are retained, without response or Up '
                  'selection.\n\n')
        log.flush()
        rows = []
        for row in sessions.itertuples():
            if row.session not in saved:
                rows.append(dict(rec_id=row.session, animal=row.animal, status='missing saved coefficients'))
                continue
            result = session_medians(row, saved[row.session], log)
            rows.append(result)
            print(f'{row.session}: ctrl={result["ctrl_red_coefficient_median"]:.9g}, '
                  f'stim={result["stim_red_coefficient_median"]:.9g}', flush=True)
        table = pd.DataFrame(rows)
        table.to_csv(output, index=False)
        ctrl, stim = 'ctrl_red_coefficient_median', 'stim_red_coefficient_median'
        paired = table.loc[np.isfinite(table[ctrl]) & np.isfinite(table[stim])]
        test = wilcoxon(paired[stim], paired[ctrl], alternative='two-sided')
        lines = ['', 'Each session median pools all finite ACR-trial b2 values within '
                 'its condition. Pooled summaries are the median and IQR across these '
                 'session medians, giving each session equal weight; IQRs describe '
                 'session dispersion. The paired two-sided Wilcoxon uses session '
                 'medians, stim against ctrl, with default zero_method=wilcox and '
                 'method=auto. Animals are counted without an additional animal-level test.', '']
        for condition in ('ctrl', 'stim'):
            group = table.loc[np.isfinite(table[f'{condition}_red_coefficient_median'])]
            low, median, high = np.percentile(group[f'{condition}_red_coefficient_median'], [25, 50, 75])
            count = int(group[f'{condition}_n_coefficients'].sum())
            lines.append(f'{condition}: pooled session-median red coefficient={median:.9g}, '
                         f'IQR [{low:.9g}, {high:.9g}], n={len(group)} sessions, '
                         f'n={group.animal.nunique()} animals, {count} finite ACR-trial coefficients.')
        low, median, high = np.percentile(paired.stim_minus_ctrl_median, [25, 50, 75])
        lines += [f'Paired stim minus ctrl: median difference={median:.9g}, '
                  f'IQR [{low:.9g}, {high:.9g}]; two-sided Wilcoxon signed-rank '
                  f'W={test.statistic:.9g}, p={test.pvalue:.12g}, n={len(paired)} sessions, '
                  f'n={paired.animal.nunique()} animals.', '',
                  '| session | animal | ACRs | ctrl median | ctrl IQR | stim median | stim IQR | stim minus ctrl |',
                  '| --- | --- | ---: | ---: | --- | ---: | --- | ---: |']
        for row in table.itertuples():
            if row.status != 'ok':
                lines.append(f'| {row.rec_id} | {row.animal} | unavailable | | | | | |')
                continue
            lines.append(f'| {row.rec_id} | {row.animal} | {row.n_acr} | '
                         f'{row.ctrl_red_coefficient_median:.9g} | '
                         f'[{row.ctrl_red_coefficient_q25:.9g}, {row.ctrl_red_coefficient_q75:.9g}] | '
                         f'{row.stim_red_coefficient_median:.9g} | '
                         f'[{row.stim_red_coefficient_q25:.9g}, {row.stim_red_coefficient_q75:.9g}] | '
                         f'{row.stim_minus_ctrl_median:.9g} |')
        lines += ['', f'Table: outputs/da_hotspots/{output.name}. No figure or other '
                  'analysis was produced; earlier scripts, regression files and '
                  'comparison outputs were left unchanged.']
        log.write('\n'.join(lines) + '\n')
        print('\n'.join(lines[3:6]), flush=True)


if __name__ == '__main__':
    main()
