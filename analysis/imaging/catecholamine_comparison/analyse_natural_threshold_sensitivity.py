# -*- coding: utf-8 -*-
'''
Created on Sat Aug 29 12:54:25 2026

build per-session natural classifier sensitivity tables for the catecholamine comparison

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

from importlib.util import module_from_spec, spec_from_file_location

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu


#%% paths and parameters
repo_root = Path(__file__).resolve().parents[3]
OUTPUT_ROOT = repo_root / 'data' / 'catecholamine_comparison' / 'evidence'
SOURCE_SCRIPT = (
    repo_root / 'analysis' / 'imaging' / 'nlight_vs_dlight'
    / 'quantify_nlight_vs_dlight.py'
)
SHUFFLE_PERCENTILES = (90, 95, 99)
EFFECT_SIZE_THRESHOLDS = (0, 0.05, 0.1)
BOOTSTRAP_SEED = 20260829
N_BOOTSTRAPS = 5000


#%% helpers
def load_source_module():
    specification = spec_from_file_location('catecholamine_nlight_dlight', SOURCE_SCRIPT)
    module = module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


def session_rows(table, sensor, percentile, effect_size):
    valid = table['dlight_valid'] & table['red_valid'] & ~table['edge']
    source = table.loc[valid].copy()
    thresholds = np.array([
        np.nanpercentile(np.asarray(values, dtype=float), percentile)
        for values in source['shuff_response_amplitude']
    ])
    classified = (
        source['response_amplitude'].to_numpy(dtype=float) > thresholds
    ) & (source['effect_size'].to_numpy(dtype=float) > effect_size)
    source['classified'] = classified
    rows = []
    for session, group in source.groupby('rec_id', sort=True):
        rows.append({
            'sensor': sensor,
            'animal': str(group['anm'].iloc[0]),
            'session': session,
            'shuffle_percentile': percentile,
            'effect_size_threshold': effect_size,
            'n_valid_acr': int(len(group)),
            'n_classified_acr': int(group['classified'].sum()),
            'up_fraction': float(group['classified'].mean()),
        })
    return rows


def animal_bootstrap(session_table, rng):
    groups = {
        sensor: [group['up_fraction'].to_numpy(dtype=float)
                 for _, group in session_table.loc[
                     session_table['sensor'].eq(sensor)
                 ].groupby('animal')]
        for sensor in ['dLight', 'nLight']
    }
    differences = np.empty(N_BOOTSTRAPS)
    for index in range(N_BOOTSTRAPS):
        medians = {}
        for sensor, animal_groups in groups.items():
            sampled_animals = rng.integers(0, len(animal_groups), len(animal_groups))
            animal_medians = [
                np.median(rng.choice(
                    animal_groups[animal], len(animal_groups[animal]), replace=True
                ))
                for animal in sampled_animals
            ]
            medians[sensor] = float(np.median(animal_medians))
        differences[index] = medians['nLight'] - medians['dLight']
    return {
        'bootstrap_difference_median': float(np.median(differences)),
        'bootstrap_ci_low': float(np.percentile(differences, 2.5)),
        'bootstrap_ci_high': float(np.percentile(differences, 97.5)),
    }


def comparison_rows(session_table):
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    rows = []
    for (percentile, effect_size), group in session_table.groupby(
            ['shuffle_percentile', 'effect_size_threshold'], sort=True):
        nlight = group.loc[group['sensor'].eq('nLight')]
        dlight = group.loc[group['sensor'].eq('dLight')]
        session_test = mannwhitneyu(
            nlight['up_fraction'], dlight['up_fraction'],
            alternative='two-sided', method='auto',
        )
        animal = group.groupby(['sensor', 'animal'])['up_fraction'].median().reset_index()
        nlight_animal = animal.loc[animal['sensor'].eq('nLight'), 'up_fraction']
        dlight_animal = animal.loc[animal['sensor'].eq('dLight'), 'up_fraction']
        animal_test = mannwhitneyu(
            nlight_animal, dlight_animal,
            alternative='two-sided', method='exact',
        )
        session_effect = 2 * session_test.statistic / (len(nlight) * len(dlight)) - 1
        rows.append({
            'shuffle_percentile': int(percentile),
            'effect_size_threshold': float(effect_size),
            'nlight_n_sessions': int(len(nlight)),
            'dlight_n_sessions': int(len(dlight)),
            'nlight_n_animals': int(nlight['animal'].nunique()),
            'dlight_n_animals': int(dlight['animal'].nunique()),
            'nlight_session_median': float(nlight['up_fraction'].median()),
            'dlight_session_median': float(dlight['up_fraction'].median()),
            'session_mann_whitney_u': float(session_test.statistic),
            'session_p_value': float(session_test.pvalue),
            'session_rank_biserial_nlight_minus_dlight': float(session_effect),
            'nlight_animal_median': float(nlight_animal.median()),
            'dlight_animal_median': float(dlight_animal.median()),
            'animal_mann_whitney_u': float(animal_test.statistic),
            'animal_exact_p_value': float(animal_test.pvalue),
            **animal_bootstrap(group, rng),
        })
    return pd.DataFrame(rows)


#%% run
def main():
    source = load_source_module()
    tables = source.load_natural_tables()
    rows = []
    for sensor, table in tables.items():
        for percentile in SHUFFLE_PERCENTILES:
            for effect_size in EFFECT_SIZE_THRESHOLDS:
                rows.extend(session_rows(table, sensor, percentile, effect_size))
    session_table = pd.DataFrame(rows).sort_values(
        ['shuffle_percentile', 'effect_size_threshold', 'sensor', 'animal', 'session']
    ).reset_index(drop=True)
    animal_table = session_table.groupby(
        ['sensor', 'animal', 'shuffle_percentile', 'effect_size_threshold'],
        as_index=False,
    ).agg(
        n_sessions=('session', 'size'),
        animal_median_up_fraction=('up_fraction', 'median'),
    )
    comparisons = comparison_rows(session_table)
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    session_table.to_csv(OUTPUT_ROOT / 'natural_threshold_session_metrics.csv', index=False)
    animal_table.to_csv(OUTPUT_ROOT / 'natural_threshold_animal_metrics.csv', index=False)
    comparisons.to_csv(OUTPUT_ROOT / 'natural_threshold_comparisons.csv', index=False)
    print(f'wrote {len(session_table)} session-setting rows')
    print(f'wrote {len(animal_table)} animal-setting rows')
    print(f'wrote {len(comparisons)} comparison rows')


if __name__ == '__main__':
    main()
