# -*- coding: utf-8 -*-
'''
Created on Tue Mar 18 16:39:49 2025

Collect off-target run-bout tables for behaviour sessions that have a
run-bout companion file.

@author: Dinghao Luo
'''

#%% imports
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from tqdm import tqdm

repo_root = Path(__file__).resolve().parents[2]
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))

import project_paths as pp
from console_formatting import print_files_saved


#%% path stems
BEHAVIOUR_SESSION_ROOT = pp.BEHAVIOUR_EXPERIMENTS_STEM
RUN_BOUT_OUTPUT_STEM = pp.RUN_BOUTS_STEM

RUN_BOUT_TABLE_PATH = RUN_BOUT_OUTPUT_STEM / 'off_target_run_bouts.csv'
RUN_BOUT_PAYLOAD_PATH = RUN_BOUT_OUTPUT_STEM / 'off_target_run_bouts_payload.npy'

CSV_SUFFIX = '_run_bouts_py.csv'
MAT_SUFFIX = '_run_bouts.mat'


#%% analysis
session_infos = [
    {
        'dataset': session_path.parent.name,
        'rec_name': session_path.stem,
        'behaviour_path': session_path,
    }
    for session_path in sorted(BEHAVIOUR_SESSION_ROOT.glob('*/*.pkl'))
]
if not session_infos:
    raise FileNotFoundError(
        f'no behaviour session pickles found in {BEHAVIOUR_SESSION_ROOT}'
    )

matched_sessions = []
skipped_sessions = []
for session_info in session_infos:
    rec_name = session_info['rec_name']
    csv_path = RUN_BOUT_OUTPUT_STEM / f'{rec_name}{CSV_SUFFIX}'
    mat_path = RUN_BOUT_OUTPUT_STEM / f'{rec_name}{MAT_SUFFIX}'

    if csv_path.exists():
        matched = dict(session_info)
        matched['run_bout_csv'] = csv_path
        matched_sessions.append(matched)
        continue

    reason = 'no run-bout companion file'
    if mat_path.exists():
        reason = 'run-bout MAT found but CSV companion missing'
    skipped_sessions.append(
        {
            'dataset': session_info['dataset'],
            'rec_name': session_info['rec_name'],
            'behaviour_path': str(session_info['behaviour_path']),
            'reason': reason,
        }
    )

if not matched_sessions:
    raise FileNotFoundError('no behaviour sessions had a run-bout CSV companion')

print(f'found {len(session_infos)} behaviour sessions in {BEHAVIOUR_SESSION_ROOT}')
print(f'processing {len(matched_sessions)} sessions with run-bout CSV companions')
if skipped_sessions:
    print(f'skipping {len(skipped_sessions)} sessions without CSV companions')

run_bout_tables = []
for session_info in tqdm(matched_sessions, desc='run-bout tables'):
    table = pd.read_csv(session_info['run_bout_csv'])
    if 'rec_name' not in table.columns:
        table.insert(0, 'rec_name', session_info['rec_name'])
    table.insert(0, 'dataset', session_info['dataset'])
    table['behaviour_path'] = str(session_info['behaviour_path'])
    table['run_bout_csv'] = str(session_info['run_bout_csv'])
    run_bout_tables.append(table)

run_bout_df = pd.concat(run_bout_tables, ignore_index=True)
RUN_BOUT_OUTPUT_STEM.mkdir(parents=True, exist_ok=True)
run_bout_df.to_csv(RUN_BOUT_TABLE_PATH, index=False)

session_tables = {
    f'{dataset}/{rec_name}': session_df.reset_index(drop=True).to_dict('list')
    for (dataset, rec_name), session_df in run_bout_df.groupby(
        ['dataset', 'rec_name']
    )
}
payload = {
    'run_bouts': run_bout_df.to_dict('list'),
    'session_run_bouts': session_tables,
    'sessions_processed': [
        {
            'dataset': info['dataset'],
            'rec_name': info['rec_name'],
            'behaviour_path': str(info['behaviour_path']),
            'run_bout_csv': str(info['run_bout_csv']),
        }
        for info in matched_sessions
    ],
    'sessions_skipped': skipped_sessions,
    'copied_session_csvs': [],
    'behaviour_session_root': BEHAVIOUR_SESSION_ROOT,
    'run_bout_companion_stems': [str(RUN_BOUT_OUTPUT_STEM)],
}
np.save(RUN_BOUT_PAYLOAD_PATH, payload)

print_files_saved([
    ('run-bout table', RUN_BOUT_TABLE_PATH),
    ('run-bout plotting data', RUN_BOUT_PAYLOAD_PATH),
])
