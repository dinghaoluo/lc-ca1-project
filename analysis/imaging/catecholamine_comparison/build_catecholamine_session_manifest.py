# -*- coding: utf-8 -*-
'''
Created on Sat Aug 29 12:33:23 2026

build one row per analysed session and event with metric-specific inclusion status

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


#%% paths
repo_root = Path(__file__).resolve().parents[3]
COMPARISON_ROOT = repo_root / 'data' / 'catecholamine_comparison'
NATURAL_ROOT = COMPARISON_ROOT / 'nlight_vs_dlight_single_state'
OUTPUT_ROOT = repo_root / 'data' / 'catecholamine_comparison' / 'evidence'
OUTPUT_PATH = OUTPUT_ROOT / 'session_manifest.csv'


#%% helpers




def natural_rows():
    path = NATURAL_ROOT / 'per_session_metrics.csv'
    natural = pd.read_csv(path)
    natural = natural.loc[
        natural['event'].eq('natural_run_onset')
        & natural['condition'].eq('untreated')
    ].copy()
    if natural['session'].duplicated().any():
        raise ValueError('natural manifest contains duplicate session IDs')
    return pd.DataFrame({
        'domain': 'natural_run_onset',
        'sensor': natural['sensor'],
        'transmitter': natural['transmitter'],
        'animal': natural['animal'],
        'session_id': natural['session'],
        'event': natural['event'],
        'condition': natural['condition'],
        'included': natural['included'],
        'exclusion_reason': natural['exclusion_reason'],
        'n_valid_acr_or_roi': natural['n_valid_roi'],
        'n_trials': np.nan,
        'n_first_lick_trials': natural['n_first_lick_trials'],
        'n_up_or_releasing': natural['n_up_roi'],
        'n_temporal_fit_roi': natural['n_temporal_fit_roi'],
        'n_spatial_fit_roi': natural['n_spatial_fit_roi'],
        'response_concentration_eligible': np.nan,
        'response_concentration_exclusion_reason': '',
        'source_table': str(path.relative_to(repo_root)),
        })


def opto_rows():
    path = COMPARISON_ROOT / 'roi_release_metrics.csv'
    release = pd.read_csv(path)
    spatial_path = COMPARISON_ROOT / 'spatial_tau_metrics.csv'
    temporal_path = COMPARISON_ROOT / 'roi_neuropil_timecourse_metrics.csv'
    spatial = pd.read_csv(spatial_path)
    spatial = spatial[spatial['window'].eq('0-1 s')]
    spatial_counts = spatial.groupby('recname')['n_roi_releasing'].sum()
    if temporal_path.exists():
        temporal_counts = pd.read_csv(temporal_path).groupby('recname').size()
    else:
        temporal_counts = pd.Series(dtype=float)

    matched = (release['n_roi'] >= 10) & (release['n_trials'] >= 10)
    reasons = np.where(matched, '', 'fewer than 10 trials or 10 ROIs')
    return pd.DataFrame({
        'domain': 'lc_stimulation',
        'sensor': release['sensor'],
        'transmitter': release['transmitter'],
        'animal': release['animal'],
        'session_id': release['recname'],
        'event': 'lc_opto_stimulation',
        'condition': 'stimulation',
        'included': True,
        'exclusion_reason': '',
        'n_valid_acr_or_roi': release['n_roi'],
        'n_trials': release['n_trials'],
        'n_first_lick_trials': np.nan,
        'n_up_or_releasing': np.round(
            release['prop_releasing'] * release['n_roi']).astype(int),
        'n_temporal_fit_roi': release['recname'].map(temporal_counts),
        'n_spatial_fit_roi': release['recname'].map(spatial_counts),
        'response_concentration_eligible': matched,
        'response_concentration_exclusion_reason': reasons,
        'source_table': str(path.relative_to(repo_root)),
        'spatial_source_table': str(spatial_path.relative_to(repo_root)),
        'temporal_source_table': str(temporal_path.relative_to(repo_root)),
        })


#%% run
def main():
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    manifest = pd.concat([natural_rows(), opto_rows()], ignore_index=True)
    manifest = manifest.sort_values(
        ['domain', 'sensor', 'animal', 'session_id']).reset_index(drop=True)
    manifest.to_csv(OUTPUT_PATH, index=False)
    print(f'wrote {len(manifest)} rows to {OUTPUT_PATH}')
    print(manifest.groupby(['domain', 'sensor']).size().to_string())


if __name__ == '__main__':
    main()
