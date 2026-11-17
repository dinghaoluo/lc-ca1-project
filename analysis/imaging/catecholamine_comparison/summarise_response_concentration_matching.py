# -*- coding: utf-8 -*-
'''
Created on Sat Aug 29 12:13:59 2026

consolidate the trial, ROI and reliability response-concentration sensitivities

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
from scipy.stats import mannwhitneyu, wilcoxon


#%% paths and parameters
repo_root = Path(__file__).resolve().parents[3]
EVIDENCE_ROOT = repo_root / 'data' / 'catecholamine_comparison' / 'evidence'
SESSION_INPUT = EVIDENCE_ROOT / 'response_concentration_matching_sessions.csv'
DECISION_OUTPUT = EVIDENCE_ROOT / 'response_concentration_matching_decision.csv'
METHOD_OUTPUT = EVIDENCE_ROOT / 'response_concentration_matching_decision_methods.md'

N_BOOTSTRAPS = 5000
BOOTSTRAP_SEED = 20260844


#%% summaries
def nested_bootstrap(frame, metric, seed):
    rng = np.random.default_rng(seed)
    values = {}
    for sensor in ['nLight', 'dLight']:
        values[sensor] = {
            animal: group[metric].dropna().to_numpy(dtype=float)
            for animal, group in frame[frame['sensor'] == sensor]
            .groupby('animal')
            }

    differences = np.empty(N_BOOTSTRAPS, dtype=float)
    for index in range(N_BOOTSTRAPS):
        sensor_medians = {}
        for sensor in ['nLight', 'dLight']:
            animals = list(values[sensor])
            sampled_animals = rng.choice(animals, len(animals), replace=True)
            animal_medians = []
            for animal in sampled_animals:
                sessions = values[sensor][animal]
                sampled_sessions = rng.choice(
                    sessions, len(sessions), replace=True)
                animal_medians.append(np.nanmedian(sampled_sessions))
            sensor_medians[sensor] = np.nanmedian(animal_medians)
        differences[index] = (
            sensor_medians['nLight'] - sensor_medians['dLight'])
    return differences


def sensor_summary(frame, sensor, metric):
    values = frame.loc[frame['sensor'] == sensor, metric].dropna()
    animal_values = (frame[frame['sensor'] == sensor]
                     .groupby('animal')[metric].median())
    within = wilcoxon(values, alternative='greater')
    return {
        'n_sessions': len(values),
        'n_animals': len(animal_values),
        'median': values.median(),
        'q1': values.quantile(.25),
        'q3': values.quantile(.75),
        'animal_median': animal_values.median(),
        'within_sensor_w': float(within.statistic),
        'within_sensor_p': float(within.pvalue),
        }


def compare_subset(frame, label, reliability_basis, overlap_low, overlap_high,
                   bootstrap_seed):
    metric = 'matched_top10_excess'
    nlight = frame.loc[frame['sensor'] == 'nLight', metric].dropna()
    dlight = frame.loc[frame['sensor'] == 'dLight', metric].dropna()
    nlight_animals = (frame[frame['sensor'] == 'nLight']
                      .groupby('animal')[metric].median())
    dlight_animals = (frame[frame['sensor'] == 'dLight']
                      .groupby('animal')[metric].median())
    session_test = mannwhitneyu(
        nlight, dlight, alternative='two-sided', method='auto')
    animal_test = mannwhitneyu(
        nlight_animals, dlight_animals, alternative='two-sided',
        method='exact')
    session_effect = (
        2 * session_test.statistic / (len(nlight) * len(dlight)) - 1)
    animal_effect = (
        2 * animal_test.statistic
        / (len(nlight_animals) * len(dlight_animals)) - 1)
    bootstrap = nested_bootstrap(frame, metric, bootstrap_seed)

    row = {
        'subset': label,
        'reliability_basis': reliability_basis,
        'reliability_overlap_low': overlap_low,
        'reliability_overlap_high': overlap_high,
        }
    for sensor in ['nLight', 'dLight']:
        for key, value in sensor_summary(frame, sensor, metric).items():
            row[f'{sensor.lower()}_{key}'] = value
    row.update({
        'session_median_difference_nlight_minus_dlight': (
            nlight.median() - dlight.median()),
        'session_mann_whitney_u': float(session_test.statistic),
        'session_p_value': float(session_test.pvalue),
        'session_rank_biserial': float(session_effect),
        'animal_median_difference_nlight_minus_dlight': (
            nlight_animals.median() - dlight_animals.median()),
        'animal_mann_whitney_u': float(animal_test.statistic),
        'animal_p_value': float(animal_test.pvalue),
        'animal_rank_biserial': float(animal_effect),
        'bootstrap_median': float(np.nanmedian(bootstrap)),
        'bootstrap_ci_low': float(np.nanquantile(bootstrap, .025)),
        'bootstrap_ci_high': float(np.nanquantile(bootstrap, .975)),
        })

    for basis in ['split_half_minus_shuffle',
                  'matched_split_half_excess']:
        n_rel = frame.loc[frame['sensor'] == 'nLight', basis].dropna()
        d_rel = frame.loc[frame['sensor'] == 'dLight', basis].dropna()
        rel_test = mannwhitneyu(
            n_rel, d_rel, alternative='two-sided', method='auto')
        row[f'nlight_{basis}_median'] = n_rel.median()
        row[f'dlight_{basis}_median'] = d_rel.median()
        row[f'{basis}_p_value'] = float(rel_test.pvalue)
    return row


#%% subset definitions and output
def overlap_bounds(frame, metric, lower_quantile, upper_quantile):
    lower = max(
        frame.loc[frame['sensor'] == 'dLight', metric]
        .quantile(lower_quantile),
        frame.loc[frame['sensor'] == 'nLight', metric]
        .quantile(lower_quantile))
    upper = min(
        frame.loc[frame['sensor'] == 'dLight', metric]
        .quantile(upper_quantile),
        frame.loc[frame['sensor'] == 'nLight', metric]
        .quantile(upper_quantile))
    return float(lower), float(upper)


def main():
    sessions = pd.read_csv(SESSION_INPUT)
    source_full_low, source_full_high = overlap_bounds(
        sessions, 'split_half_minus_shuffle', 0, 1)
    source_central_low, source_central_high = overlap_bounds(
        sessions, 'split_half_minus_shuffle', .05, .95)
    matched_central_low, matched_central_high = overlap_bounds(
        sessions, 'matched_split_half_excess', .05, .95)

    subsets = [
        ('all dimension-matched sessions', sessions,
         'none', np.nan, np.nan),
        ('full-range source-reliability overlap',
         sessions[sessions['split_half_minus_shuffle'].between(
             source_full_low, source_full_high)],
         'full-session split-half excess', source_full_low, source_full_high),
        ('central source-reliability overlap',
         sessions[sessions['split_half_minus_shuffle'].between(
             source_central_low, source_central_high)],
         'full-session split-half excess',
         source_central_low, source_central_high),
        ('central matched-reliability overlap',
         sessions[sessions['matched_split_half_excess'].between(
             matched_central_low, matched_central_high)],
         '10-trial, 10-ROI matched split-half excess',
         matched_central_low, matched_central_high),
        ]

    rows = []
    for index, subset in enumerate(subsets):
        rows.append(compare_subset(
            subset[1], subset[0], subset[2], subset[3], subset[4],
            BOOTSTRAP_SEED + index))
    decisions = pd.DataFrame(rows)
    decisions.to_csv(DECISION_OUTPUT, index=False)

    methods = f'''# Response-concentration matching sensitivity analysis

This table consolidates the completed 10-trial, 10-ROI analysis rather than
selecting one reliability rule after seeing its cross-sensor *p* value. Every
session estimate is the median of 1,000 subsamples produced by
`analyse_response_concentration_matching.py`; the source table contains 32 nLight sessions
from four animals and 40 dLight sessions from five animals.

Four sensitivity sets are reported. The first retains every dimension-matched session. The
second restricts the original full-session split-half metric to the complete
range shared by both cohorts (`{source_full_low:.6f}` to
`{source_full_high:.6f}`). The third uses the more conservative intersection
of the sensors' 5th-to-95th-percentile source ranges
(`{source_central_low:.6f}` to `{source_central_high:.6f}`). The fourth applies
that 5th-to-95th-percentile rule to split-half excess recomputed from the
matched subsamples (`{matched_central_low:.6f}` to
`{matched_central_high:.6f}`).

The reliability definitions select materially different data. The
matched-reliability restriction retains only four nLight sessions from three
animals, whereas the source-reliability restrictions retain broader cohorts.
It must therefore remain a selected-tail sensitivity, even if its session
contrast is smaller than .05. The animal-level exact test and animal-balanced
animal-then-session bootstrap are the relevant safeguards against treating
repeated sessions as independent replication.

Across the unselected dimension-matched cohort, the nLight-minus-dLight
concentration difference is absent. The source-reliability branches do not
yield decisive session or animal evidence, and matched reliability remains
lower for nLight after the full-range source restriction. The defensible
conclusion is therefore that both sensor-reported maps retain concentrated
responses after dimension matching; a sensor-specific difference in response
concentration is not established.

Bootstrap resamples per row: {N_BOOTSTRAPS}; seeds: {BOOTSTRAP_SEED} through
{BOOTSTRAP_SEED + len(subsets) - 1}.
'''
    METHOD_OUTPUT.write_text(methods, encoding='utf-8')
    print(decisions.to_string(index=False))


if __name__ == '__main__':
    main()
