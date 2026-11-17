# -*- coding: utf-8 -*-
'''
Created on Sat Aug 29 13:02:05 2026

summarise fit success and censoring from the catecholamine comparison tables

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
CORRECTED_ROOT = COMPARISON_ROOT / 'nlight_vs_dlight_single_state'
OUTPUT_ROOT = repo_root / 'data' / 'catecholamine_comparison' / 'evidence'


#%% helpers
def finite_median(values):
    values = pd.to_numeric(values, errors='coerce').replace([np.inf, -np.inf], np.nan).dropna()
    return float(values.median()) if len(values) else np.nan


def reason_summary(frame):
    reasons = frame.loc[~frame['included'].astype(bool), 'reason'].fillna('')
    reasons = reasons.replace('', 'unspecified exclusion')
    return '; '.join(f'{reason} (n={count})'
                     for reason, count in reasons.value_counts().items())


def build_rows():
    inclusion = pd.read_csv(CORRECTED_ROOT / 'inclusion_exclusion.csv')
    session = pd.read_csv(CORRECTED_ROOT / 'per_session_metrics.csv')
    roi = pd.read_csv(CORRECTED_ROOT / 'per_roi_fits.csv')
    rows = []
    for (sensor, event, metric), group in inclusion.groupby(
            ['sensor', 'event', 'metric'], sort=True):
        if metric not in ('temporal_kinetics', 'spatial_kinetics'):
            continue
        recruitment = session.loc[
            session['sensor'].eq(sensor) & session['event'].eq(event)
        ]
        included = group['included'].astype(bool)
        units = int(group.loc[included, 'n_roi_or_trials'].sum())
        fit_rows = roi.loc[
            roi['sensor'].eq(sensor) & roi['event'].eq(event)
        ]
        if metric == 'temporal_kinetics':
            fit_rows = fit_rows.loc[fit_rows['temporal_included'].astype(bool)]
            r2 = finite_median(fit_rows['temporal_r2'])
        else:
            fit_rows = fit_rows.loc[fit_rows['spatial_included'].astype(bool)]
            r2 = finite_median(fit_rows['spatial_r2'])
        rows.append({
            'sensor': sensor,
            'event': event,
            'metric': metric,
            'n_eligible_sessions': int(len(recruitment)),
            'n_included_sessions': int(included.sum()),
            'n_excluded_sessions': int((~included).sum()),
            'session_success_fraction': float(included.mean()),
            'n_included_fit_units_recorded': units,
            'median_included_fit_r2': r2,
            'exclusion_reasons': reason_summary(group),
            'pre_event_noise_fields_available': False,
            'status': 'diagnostic',
        })
    return pd.DataFrame(rows)


#%% run
def main():
    output = build_rows()
    output.to_csv(OUTPUT_ROOT / 'fit_success_censoring.csv', index=False)
    print(f'wrote {len(output)} fit-success rows')


if __name__ == '__main__':
    main()
