# -*- coding: utf-8 -*-
'''
Created on Sat Aug 29 15:23:31 2026

test local spatial organisation of optogenetic reporter maps with geometry-preserving nulls

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

from _sha256 import sha256

import numpy as np
import pandas as pd
from scipy.stats import binomtest, mannwhitneyu, wilcoxon

np_core = np.core
sys.modules.setdefault('numpy._core', np_core)
sys.modules.setdefault('numpy._core.multiarray', np_core.multiarray)
sys.modules.setdefault('numpy._core.numeric', np_core.numeric)
try:
    sys.modules.setdefault('numpy._core._multiarray_umath',
                           np_core._multiarray_umath)
except:
    pass

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from common_functions import mpl_formatting
mpl_formatting()


#%% paths and parameters
repo_root = Path(__file__).resolve().parents[3]
INPUT_ROOT = repo_root / 'data' / 'catecholamine_comparison'
OUTPUT_ROOT = repo_root / 'data' / 'catecholamine_comparison' / 'evidence'
FIGURE_ROOT = repo_root / 'figures' / 'catecholamine_comparison' / 'evidence'

PROCESSED_ROOTS = {
    'dLight': repo_root / 'data' / 'hpc_dlight' / 'lc_opto' / 'all_sessions',
    'nLight': repo_root / 'data' / 'hpc_nlight' / 'lc_opto' / 'all_sessions',
}

N_PERMUTATIONS = 999
N_MATCHED_PERMUTATIONS = 199
N_MATCHED_RESAMPLES = 1000
N_MATCHED_TRIALS = 10
N_MATCHED_ROIS = 10
N_BOOTSTRAPS = 5000
RNG_SEED = 20260829

SENSOR_ORDER = ['dLight', 'nLight']
SENSOR_COLOURS = {'dLight': '#4f8f62', 'nLight': '#8b5ea7'}


#%% map statistics
def load_roi_trial_matrix(sensor, recname):
    processed = PROCESSED_ROOTS[sensor] / recname / 'processed_data'
    response_path = processed / f'{recname}_pixel_RI_stim.npy'
    roi_path = processed / f'{recname}_ROI_dict.npy'
    response = np.load(response_path, mmap_mode='r')
    roi_dict = np.load(roi_path, allow_pickle=True).item()

    values = []
    centres = []
    for roi in roi_dict.values():
        ypix = np.asarray(roi['ypix'], dtype=int)
        xpix = np.asarray(roi['xpix'], dtype=int)
        values.append(np.nanmean(response[ypix, xpix, :], axis=0))
        centres.append([np.mean(xpix), np.mean(ypix)])

    matrix = np.asarray(values, dtype=float)
    centres = np.asarray(centres, dtype=float)
    keep = np.sum(np.isfinite(matrix), axis=1) >= 5
    return matrix[keep], centres[keep], response_path, roi_path


def spatial_weights(centres):
    distances = np.sqrt(
        np.sum((centres[:, None, :] - centres[None, :, :]) ** 2, axis=2)
    )
    np.fill_diagonal(distances, np.inf)
    bandwidth = float(np.median(np.min(distances, axis=1)))
    if not np.isfinite(bandwidth) or bandwidth <= 0:
        raise ValueError('ROI centroids do not define a positive bandwidth')
    weights = np.exp(-0.5 * (distances / bandwidth) ** 2)
    np.fill_diagonal(weights, 0)
    return weights, bandwidth


def moran_i(values, weights):
    values = np.asarray(values, dtype=float)
    if not np.isfinite(values).all():
        return np.nan
    centred = values - np.mean(values)
    denominator = np.sum(centred ** 2)
    weight_sum = np.sum(weights)
    if denominator <= 0 or weight_sum <= 0:
        return np.nan
    return float(
        len(values) / weight_sum
        * (centred @ weights @ centred)
        / denominator
    )


def detrend_field(values, centres):
    design = np.column_stack([np.ones(len(centres)), centres])
    coefficients = np.linalg.lstsq(design, values, rcond=None)[0]
    return values - design @ coefficients


def locality_summary(values, centres, rng, n_permutations):
    values = np.asarray(values, dtype=float)
    valid = np.isfinite(values)
    values = values[valid]
    centres = np.asarray(centres, dtype=float)[valid]
    if len(values) < 5:
        return {key: np.nan for key in [
            'moran_i', 'null_median', 'null_q05', 'null_q95',
            'excess', 'empirical_p', 'bandwidth_px',
        ]}

    weights, bandwidth = spatial_weights(centres)
    observed = moran_i(values, weights)
    orders = np.asarray(
        [rng.permutation(len(values)) for _ in range(n_permutations)]
    )
    permuted = values[orders]
    centred = permuted - np.mean(values)
    denominator = np.sum((values - np.mean(values)) ** 2)
    null_values = (
        len(values) / np.sum(weights)
        * np.einsum('bi,ij,bj->b', centred, weights, centred)
        / denominator
    )
    return {
        'moran_i': observed,
        'null_median': float(np.median(null_values)),
        'null_q05': float(np.quantile(null_values, .05)),
        'null_q95': float(np.quantile(null_values, .95)),
        'excess': float(observed - np.median(null_values)),
        'empirical_p': float(
            (1 + np.sum(null_values >= observed)) / (n_permutations + 1)
        ),
        'bandwidth_px': bandwidth,
    }


def session_seed(sensor, recname):
    token = f'{RNG_SEED}:{sensor}:{recname}'.encode('utf-8')
    return int.from_bytes(sha256(token).digest()[:4], 'little')


def analyse_full_session(matrix, centres, seed):
    rng = np.random.default_rng(seed)
    response = np.nanmean(matrix, axis=1)
    valid = np.isfinite(response)
    raw = locality_summary(response, centres, rng, N_PERMUTATIONS)
    detrended_response = detrend_field(response[valid], centres[valid])
    detrended = locality_summary(
        detrended_response, centres[valid], rng, N_PERMUTATIONS
    )
    return raw, detrended


def analyse_matched_session(matrix, centres, seed):
    rng = np.random.default_rng(seed)
    records = []
    for _ in range(N_MATCHED_RESAMPLES):
        roi_indices = rng.choice(
            len(matrix), size=N_MATCHED_ROIS, replace=False
        )
        trial_indices = rng.choice(
            matrix.shape[1], size=N_MATCHED_TRIALS, replace=False
        )
        response = np.nanmean(
            matrix[np.ix_(roi_indices, trial_indices)], axis=1
        )
        valid = np.isfinite(response)
        if valid.sum() < N_MATCHED_ROIS:
            continue
        sample_centres = centres[roi_indices][valid]
        response = response[valid]
        raw = locality_summary(
            response, sample_centres, rng, N_MATCHED_PERMUTATIONS
        )
        records.append({
            'raw_moran_i': raw['moran_i'],
            'raw_null_median': raw['null_median'],
            'raw_excess': raw['excess'],
            'raw_empirical_p': raw['empirical_p'],
            'detrended_moran_i': np.nan,
            'detrended_null_median': np.nan,
            'detrended_excess': np.nan,
            'detrended_empirical_p': np.nan,
            'bandwidth_px': raw['bandwidth_px'],
        })
    return pd.DataFrame(records)


#%% summaries and comparisons
def equal_animal_bootstrap(frame, metric, seed):
    rng = np.random.default_rng(seed)
    values = {
        sensor: {
            animal: group[metric].dropna().to_numpy(dtype=float)
            for animal, group in frame[frame['sensor'].eq(sensor)].groupby('animal')
        }
        for sensor in SENSOR_ORDER
    }
    if not values['dLight'] or not values['nLight']:
        return {'median': np.nan, 'ci_low': np.nan, 'ci_high': np.nan}

    differences = np.empty(N_BOOTSTRAPS)
    for index in range(N_BOOTSTRAPS):
        sensor_medians = {}
        for sensor in SENSOR_ORDER:
            animals = list(values[sensor])
            sampled = rng.choice(animals, len(animals), replace=True)
            animal_medians = [
                np.median(rng.choice(
                    values[sensor][animal],
                    len(values[sensor][animal]),
                    replace=True,
                ))
                for animal in sampled
            ]
            sensor_medians[sensor] = np.median(animal_medians)
        differences[index] = (
            sensor_medians['nLight'] - sensor_medians['dLight']
        )
    return {
        'median': float(np.median(differences)),
        'ci_low': float(np.quantile(differences, .025)),
        'ci_high': float(np.quantile(differences, .975)),
    }


def comparison_row(frame, analysis, metric, seed):
    current = frame[frame['analysis'].eq(analysis)].dropna(subset=[metric])
    nlight = current.loc[current['sensor'].eq('nLight'), metric].to_numpy()
    dlight = current.loc[current['sensor'].eq('dLight'), metric].to_numpy()
    session_test = mannwhitneyu(
        nlight, dlight, alternative='two-sided', method='auto'
    )
    animal = current.groupby(['sensor', 'animal'], as_index=False)[metric].median()
    nlight_animal = animal.loc[animal['sensor'].eq('nLight'), metric].to_numpy()
    dlight_animal = animal.loc[animal['sensor'].eq('dLight'), metric].to_numpy()
    animal_test = mannwhitneyu(
        nlight_animal, dlight_animal, alternative='two-sided', method='exact'
    )
    return {
        'analysis': analysis,
        'metric': metric,
        'nlight_n_sessions': len(nlight),
        'nlight_n_animals': len(nlight_animal),
        'nlight_median': float(np.median(nlight)),
        'nlight_q1': float(np.quantile(nlight, .25)),
        'nlight_q3': float(np.quantile(nlight, .75)),
        'dlight_n_sessions': len(dlight),
        'dlight_n_animals': len(dlight_animal),
        'dlight_median': float(np.median(dlight)),
        'dlight_q1': float(np.quantile(dlight, .25)),
        'dlight_q3': float(np.quantile(dlight, .75)),
        'session_mann_whitney_u': float(session_test.statistic),
        'session_p_value': float(session_test.pvalue),
        'session_rank_biserial': float(
            2 * session_test.statistic / (len(nlight) * len(dlight)) - 1
        ),
        'nlight_animal_median': float(np.median(nlight_animal)),
        'dlight_animal_median': float(np.median(dlight_animal)),
        'animal_mann_whitney_u': float(animal_test.statistic),
        'animal_exact_p_value': float(animal_test.pvalue),
        'animal_rank_biserial': float(
            2 * animal_test.statistic
            / (len(nlight_animal) * len(dlight_animal))
            - 1
        ),
        **{
            key: value
            for key, value in zip(
                ['bootstrap_difference_median', 'bootstrap_ci_low',
                 'bootstrap_ci_high'],
                equal_animal_bootstrap(current, metric, seed).values(),
            )
        },
        'n_bootstraps': N_BOOTSTRAPS,
        'bootstrap_seed': seed,
        'nlight_session_wilcoxon_p': float(
            wilcoxon(nlight, alternative='greater', method='auto').pvalue
        ),
        'dlight_session_wilcoxon_p': float(
            wilcoxon(dlight, alternative='greater', method='auto').pvalue
        ),
        'nlight_n_positive_animals': int(
            (nlight_animal > 0).sum()
        ),
        'dlight_n_positive_animals': int(
            (dlight_animal > 0).sum()
        ),
        'nlight_animal_sign_p_greater': float(
            binomtest((nlight_animal > 0).sum(), len(nlight_animal),
                      .5, alternative='greater').pvalue
        ),
        'dlight_animal_sign_p_greater': float(
            binomtest((dlight_animal > 0).sum(), len(dlight_animal),
                      .5, alternative='greater').pvalue
        ),
    }


def select_representatives(sessions):
    rows = []
    for sensor in SENSOR_ORDER:
        current = sessions[
            sessions['sensor'].eq(sensor)
            & sessions['analysis'].eq('dimension_matched')
        ].dropna(subset=['raw_excess'])
        target = current['raw_excess'].median()
        selected = current.iloc[
            np.argmin(np.abs(current['raw_excess'] - target))
        ]
        rows.append({
            'sensor': sensor,
            'session': selected['session'],
            'animal': selected['animal'],
            'selection_rule': 'nearest to sensor median matched raw excess',
        })
    return pd.DataFrame(rows)


#%% figure
def plot_results(sessions, release, representatives):
    figure, axes = plt.subplots(
        3, 2, figsize=(8.2, 10.0), constrained_layout=True
    )
    for column, row in enumerate(representatives.itertuples(index=False)):
        matrix, centres, _, _ = load_roi_trial_matrix(row.sensor, row.session)
        response = np.nanmean(matrix, axis=1)
        limits = np.nanpercentile(response, [5, 95])
        axes[0, column].scatter(
            centres[:, 0], centres[:, 1], c=response, cmap='viridis',
            vmin=limits[0], vmax=limits[1], s=42, edgecolor='white',
            linewidth=.35,
        )
        axes[0, column].set_title(
            f'{chr(65 + column)}  {row.sensor} representative map',
            loc='left', color=SENSOR_COLOURS[row.sensor],
        )
        axes[0, column].set_aspect('equal')
        axes[0, column].invert_yaxis()
        axes[0, column].set_xlabel('ROI x (pixels)')
        axes[0, column].set_ylabel('ROI y (pixels)')

    for column, analysis in enumerate(['full_session', 'dimension_matched']):
        ax = axes[1, column]
        for xpos, sensor in enumerate(SENSOR_ORDER):
            values = sessions.loc[
                sessions['sensor'].eq(sensor)
                & sessions['analysis'].eq(analysis),
                'detrended_excess',
            ].dropna().to_numpy()
            jitter = np.random.default_rng(RNG_SEED + column + xpos).uniform(
                -.1, .1, len(values)
            )
            ax.scatter(
                np.full(len(values), xpos) + jitter, values,
                color=SENSOR_COLOURS[sensor], alpha=.55, s=18,
            )
            animal_values = sessions.loc[
                sessions['sensor'].eq(sensor)
                & sessions['analysis'].eq(analysis)
            ].groupby('animal')['detrended_excess'].median().dropna()
            ax.scatter(
                np.full(len(animal_values), xpos), animal_values,
                marker='o', facecolor='white', edgecolor=SENSOR_COLOURS[sensor],
                s=42, zorder=3,
            )
            ax.hlines(
                np.median(values), xpos - .16, xpos + .16,
                color='black', linewidth=1.5,
            )
        ax.set_xticks([0, 1], SENSOR_ORDER)
        ax.set_ylabel('Detrended Moran excess')
        ax.set_title(
            f'{chr(67 + column)}  {analysis.replace("_", " ")}',
            loc='left',
        )

    ax = axes[2, 0]
    for sensor in SENSOR_ORDER:
        current = sessions[
            sessions['sensor'].eq(sensor)
            & sessions['analysis'].eq('dimension_matched')
        ]
        ax.scatter(
            current['raw_excess'], current['detrended_excess'],
            color=SENSOR_COLOURS[sensor], alpha=.65, s=20, label=sensor,
        )
    ax.set_xlabel('Raw Moran excess')
    ax.set_ylabel('Detrended Moran excess')
    ax.set_title('E  Field-gradient sensitivity', loc='left')
    ax.legend(frameon=False, fontsize=8)

    ax = axes[2, 1]
    merged = sessions[
        sessions['analysis'].eq('dimension_matched')
    ].merge(
        release[['sensor', 'recname', 'split_half_minus_shuffle']],
        left_on=['sensor', 'session'],
        right_on=['sensor', 'recname'],
        how='left',
        validate='one_to_one',
    )
    for sensor in SENSOR_ORDER:
        current = merged[merged['sensor'].eq(sensor)]
        ax.scatter(
            current['split_half_minus_shuffle'],
            current['detrended_excess'],
            color=SENSOR_COLOURS[sensor], alpha=.65, s=20, label=sensor,
        )
    ax.set_xlabel('Original split-half excess')
    ax.set_ylabel('Matched detrended Moran excess')
    ax.set_title('F  Locality and map reliability', loc='left')
    ax.legend(frameon=False, fontsize=8)
    for ax in axes.ravel():
        ax.spines[['top', 'right']].set_visible(False)
    return figure


#%% run
def main():
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    FIGURE_ROOT.mkdir(parents=True, exist_ok=True)
    source = pd.read_csv(INPUT_ROOT / 'roi_release_metrics.csv')
    source = source[source['sensor'].isin(SENSOR_ORDER)].copy()
    rows = []
    for source_row in source.itertuples(index=False):
        if source_row.n_roi < N_MATCHED_ROIS or source_row.n_trials < N_MATCHED_TRIALS:
            continue
        matrix, centres, response_path, roi_path = load_roi_trial_matrix(
            source_row.sensor, source_row.recname
        )
        if matrix.shape != (source_row.n_roi, source_row.n_trials):
            raise ValueError(
                f'{source_row.recname}: source and loaded dimensions disagree'
            )
        seed = session_seed(source_row.sensor, source_row.recname)
        raw, detrended = analyse_full_session(matrix, centres, seed)
        matched = analyse_matched_session(matrix, centres, seed + 1)
        if matched.empty:
            raise ValueError(f'{source_row.recname}: no matched resample survived')
        for analysis, values in [
            ('full_session', {
                'raw_moran_i': raw['moran_i'],
                'raw_null_median': raw['null_median'],
                'raw_excess': raw['excess'],
                'raw_empirical_p': raw['empirical_p'],
                'detrended_moran_i': detrended['moran_i'],
                'detrended_null_median': detrended['null_median'],
                'detrended_excess': detrended['excess'],
                'detrended_empirical_p': detrended['empirical_p'],
                'bandwidth_px': detrended['bandwidth_px'],
                'matched_resamples': np.nan,
            }),
            ('dimension_matched', {
                'raw_moran_i': matched['raw_moran_i'].median(),
                'raw_null_median': matched['raw_null_median'].median(),
                'raw_excess': matched['raw_excess'].median(),
                'raw_empirical_p': matched['raw_empirical_p'].median(),
                'detrended_moran_i': np.nan,
                'detrended_null_median': np.nan,
                'detrended_excess': np.nan,
                'detrended_empirical_p': np.nan,
                'bandwidth_px': matched['bandwidth_px'].median(),
                'matched_resamples': len(matched),
            }),
        ]:
            rows.append({
                'sensor': source_row.sensor,
                'transmitter': source_row.transmitter,
                'session': source_row.recname,
                'animal': source_row.animal,
                'analysis': analysis,
                'n_roi': matrix.shape[0],
                'n_trials': matrix.shape[1],
                **values,
                'response_source': str(response_path.relative_to(repo_root)),
                'roi_source': str(roi_path.relative_to(repo_root)),
            })
        print(f'{source_row.sensor} {source_row.recname}', flush=True)

    sessions = pd.DataFrame(rows)
    comparisons = pd.DataFrame([
        comparison_row(sessions, analysis, metric, RNG_SEED + index)
        for index, (analysis, metric) in enumerate([
            ('full_session', 'raw_excess'),
            ('full_session', 'detrended_excess'),
            ('dimension_matched', 'raw_excess'),
        ])
    ])
    animals = (
        sessions.groupby(['analysis', 'sensor', 'animal'], as_index=False)
        .agg(
            n_sessions=('session', 'size'),
            raw_excess=('raw_excess', 'median'),
            detrended_excess=('detrended_excess', 'median'),
        )
    )
    representatives = select_representatives(sessions)
    sessions.to_csv(OUTPUT_ROOT / 'spatial_locality_sessions.csv', index=False)
    animals.to_csv(OUTPUT_ROOT / 'spatial_locality_animals.csv', index=False)
    comparisons.to_csv(OUTPUT_ROOT / 'spatial_locality_comparisons.csv', index=False)
    representatives.to_csv(
        OUTPUT_ROOT / 'spatial_locality_representatives.csv', index=False
    )
    figure = plot_results(sessions, source, representatives)
    figure.savefig(
        FIGURE_ROOT / 'spatial_locality_sensitivity.png',
        dpi=300,
        bbox_inches='tight',
    )
    figure.savefig(
        FIGURE_ROOT / 'spatial_locality_sensitivity.svg',
        bbox_inches='tight',
    )
    plt.close(figure)

    print(f'wrote {len(sessions)} session-analysis rows')
    print(comparisons.to_string(index=False))


if __name__ == '__main__':
    main()
