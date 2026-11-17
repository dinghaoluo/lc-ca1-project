# -*- coding: utf-8 -*-
'''
Created on Sat Aug 29 10:37:54 2026

build the supporting statistical analyses for catecholamine-reporter maps

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

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from common_functions import mpl_formatting
mpl_formatting()
import numpy as np
import pandas as pd
from scipy.stats import binomtest, mannwhitneyu, wilcoxon


#%% paths and parameters
repo_root       = Path(__file__).resolve().parents[3]
COMPARISON_ROOT = repo_root / 'data' / 'catecholamine_comparison'
CORRECTED_ROOT  = COMPARISON_ROOT / 'nlight_vs_dlight_single_state'
OUTPUT_ROOT     = repo_root / 'data' / 'catecholamine_comparison' / 'evidence'
FIGURE_ROOT     = repo_root / 'figures' / 'catecholamine_comparison' / 'evidence'

BOOTSTRAP_SEED = 20260829
N_BOOTSTRAPS   = 5000

SENSOR_ORDER   = ['dLight', 'nLight']
SENSOR_COLOURS = {'dLight': '#2E5D3E', 'nLight': '#2050A0'}


#%% statistical helpers
def summarise(values):
    values = pd.Series(values, dtype=float).replace([np.inf, -np.inf], np.nan).dropna()
    return {
        'n': int(len(values)),
        'median': float(values.median()) if len(values) else np.nan,
        'q1': float(values.quantile(0.25)) if len(values) else np.nan,
        'q3': float(values.quantile(0.75)) if len(values) else np.nan,
    }


def values_by_animal(frame, metric):
    return [group[metric].dropna().to_numpy(dtype=float)
            for _, group in frame.groupby('animal')]


def resample_sensor_median(groups, rng):
    sampled = rng.integers(0, len(groups), size=len(groups))
    animal_medians = [
        np.median(rng.choice(groups[index], size=len(groups[index]), replace=True))
        for index in sampled
    ]
    return float(np.median(animal_medians))


def hierarchical_difference(frame, metric):
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    dlight = values_by_animal(frame[frame['sensor'] == 'dLight'], metric)
    nlight = values_by_animal(frame[frame['sensor'] == 'nLight'], metric)
    differences = np.empty(N_BOOTSTRAPS)
    for index in range(N_BOOTSTRAPS):
        differences[index] = (
            resample_sensor_median(nlight, rng)
            - resample_sensor_median(dlight, rng)
        )
    return {
        'bootstrap_difference_median': float(np.median(differences)),
        'bootstrap_ci_low': float(np.percentile(differences, 2.5)),
        'bootstrap_ci_high': float(np.percentile(differences, 97.5)),
    }


def compare_sensors(frame, metric, label, domain, family, unit, source, status, caveat):
    current = frame[['sensor', 'animal', 'session', metric]].dropna().copy()
    nlight = current.loc[current['sensor'] == 'nLight', metric].astype(float)
    dlight = current.loc[current['sensor'] == 'dLight', metric].astype(float)
    test = mannwhitneyu(nlight, dlight, alternative='two-sided', method='auto')
    effect = 2 * test.statistic / (len(nlight) * len(dlight)) - 1
    n_summary = summarise(nlight)
    d_summary = summarise(dlight)
    animal = current.groupby(['sensor', 'animal'], as_index=False)[metric].median()
    nlight_animal = animal.loc[animal['sensor'] == 'nLight', metric].to_numpy(float)
    dlight_animal = animal.loc[animal['sensor'] == 'dLight', metric].to_numpy(float)
    animal_test = mannwhitneyu(
        nlight_animal, dlight_animal, alternative='two-sided', method='exact'
    )
    return {
        'domain': domain,
        'family': family,
        'metric': metric,
        'label': label,
        'unit': unit,
        'source': source,
        'status': status,
        'caveat': caveat,
        'nlight_n_sessions': n_summary['n'],
        'nlight_n_animals': current.loc[current['sensor'] == 'nLight', 'animal'].nunique(),
        'nlight_median': n_summary['median'],
        'nlight_q1': n_summary['q1'],
        'nlight_q3': n_summary['q3'],
        'dlight_n_sessions': d_summary['n'],
        'dlight_n_animals': current.loc[current['sensor'] == 'dLight', 'animal'].nunique(),
        'dlight_median': d_summary['median'],
        'dlight_q1': d_summary['q1'],
        'dlight_q3': d_summary['q3'],
        'mann_whitney_u': float(test.statistic),
        'p_value': float(test.pvalue),
        'rank_biserial_nlight_minus_dlight': float(effect),
        'nlight_animal_median': float(np.median(nlight_animal)),
        'dlight_animal_median': float(np.median(dlight_animal)),
        'animal_mann_whitney_u': float(animal_test.statistic),
        'animal_exact_p_value': float(animal_test.pvalue),
        'animal_rank_biserial_nlight_minus_dlight': float(
            2 * animal_test.statistic
            / (len(nlight_animal) * len(dlight_animal))
            - 1
        ),
        **hierarchical_difference(current, metric),
    }


def add_holm_correction(results):
    results = results.copy()
    results['holm_p'] = np.nan
    for indices in results.groupby('family').groups.values():
        ordered = results.loc[indices, 'p_value'].sort_values()
        adjusted = []
        running = 0
        for rank, value in enumerate(ordered, start=1):
            candidate = min(1, value * (len(ordered) - rank + 1))
            running = max(running, candidate)
            adjusted.append(running)
        results.loc[ordered.index, 'holm_p'] = adjusted
    return results


def mouse_level_result(frame, metric, label, domain):
    medians = frame.groupby(['sensor', 'animal'], as_index=False)[metric].median().dropna()
    nlight = medians.loc[medians['sensor'] == 'nLight', metric]
    dlight = medians.loc[medians['sensor'] == 'dLight', metric]
    test = mannwhitneyu(nlight, dlight, alternative='two-sided', method='auto')
    return {
        'domain': domain,
        'metric': metric,
        'label': label,
        'nlight_n_animals': len(nlight),
        'nlight_animal_median': nlight.median(),
        'dlight_n_animals': len(dlight),
        'dlight_animal_median': dlight.median(),
        'mann_whitney_u': float(test.statistic),
        'p_value': float(test.pvalue),
    }


#%% data assembly
def load_inputs():
    natural = pd.read_csv(CORRECTED_ROOT / 'per_session_metrics.csv')
    natural = natural[natural['event'] == 'natural_run_onset'].copy()
    release = pd.read_csv(COMPARISON_ROOT / 'roi_release_metrics.csv').rename(
        columns={'recname': 'session'}
    )
    temporal = pd.read_csv(COMPARISON_ROOT / 'roi_neuropil_timecourse_metrics.csv').rename(
        columns={'recname': 'session'}
    )
    normalised = pd.read_csv(COMPARISON_ROOT / 'sensor_normalised_session_metrics.csv').rename(
        columns={'recname': 'session'}
    )
    ratio = normalised['roi_tau_s'] / normalised['sensor_tau_off_s']
    difference = (ratio - normalised['roi_tau_sensor_ratio']).abs().dropna()
    if len(difference) and difference.max() > 1e-12:
        raise AssertionError('roi_tau_sensor_ratio does not equal roi_tau_s / sensor_tau_off_s')

    spatial = pd.read_csv(COMPARISON_ROOT / 'spatial_tau_metrics.csv').rename(
        columns={'recname': 'session'}
    )
    spatial = spatial[spatial['window'].eq('0-1 s')].copy()
    spatial['response_qualified_fraction'] = (
        spatial['n_roi_releasing'] / spatial['n_roi_total']
    )
    return natural, release, temporal, normalised, spatial


def load_spatial_geometry_diagnostic():
    geometry = pd.read_csv(
        COMPARISON_ROOT / 'opto_session_spatial_calibration.csv',
        )
    common = geometry['common_zoom3_geometry'].astype(bool)
    excluded = set(geometry.loc[~common, 'session'])
    expected_excluded = {
        'A152i-20251113-01',
        'A152i-20251114-01',
        }
    if len(geometry) != 47 or common.sum() != 45:
        raise AssertionError(
            'expected 45/47 common 512-by-512, zoom-three fields'
            )
    if excluded != expected_excluded:
        raise AssertionError(
            f'unexpected common-geometry exclusions: {sorted(excluded)}'
            )
    return geometry


def build_results(natural, release, temporal, normalised, spatial):
    specifications = [
        (natural, 'up_fraction', 'Run-onset Up-region fraction', 'Natural run onset',
         'natural_primary', 'Fraction', 'Corrected single-state natural table', 'Exploratory',
         'Threshold-dependent; only two nLight animals'),
        (natural, 'response_amplitude_all_dff', 'All-region response', 'Natural run onset',
         'natural_primary', 'dF/F', 'Corrected single-state natural table', 'Descriptive',
         'Raw fluorescence scales differ between sensors'),
        (natural, 'response_amplitude_up_dff', 'Up-region response', 'Natural run onset',
         'natural_primary', 'dF/F', 'Corrected single-state natural table', 'Descriptive',
         'Conditioned on passing the response criterion'),
        (natural, 'temporal_tau_s', 'Run-onset decay time constant', 'Natural run onset',
         'natural_kinetics', 's', 'Corrected single-state natural table', 'Secondary',
         'Fit-qualified responder sessions only'),
        (natural, 'rise_time_s', 'Run-onset rise time', 'Natural run onset',
         'natural_kinetics', 's', 'Corrected single-state natural table', 'Secondary',
         'Fit-qualified responder sessions only'),
        (natural, 'first_lick_rho', 'Response versus first-lick latency', 'Natural run onset',
         'natural_behaviour', 'Spearman rho', 'Corrected single-state natural table', 'Exploratory',
         'Responder sessions and valid post-onset licks only'),
        (release, 'gini_positive_RI', 'Positive-response Gini coefficient', 'LC stimulation',
         'opto_topology', 'Coefficient', 'ROI release metrics', 'Descriptive',
         'Concentration depends partly on signal-to-noise ratio'),
        (release, 'top10_minus_shuffle', 'Top-10% excess contribution', 'LC stimulation',
         'opto_topology', 'Fraction', 'ROI release metrics', 'Primary within reporter',
         'Observed contribution minus within-session ROI-identity shuffle'),
        (release, 'split_half_minus_shuffle', 'Split-half map stability above shuffle',
         'LC stimulation', 'opto_topology', 'Spearman rho', 'ROI release metrics', 'Validation',
         'nLight remains reproducible but is less stable than dLight'),
        (temporal, 'roi_tau_s', 'Stimulus-evoked decay time constant', 'LC stimulation',
         'opto_temporal_fit', 's', 'ROI-neuropil time-course metrics',
         'Primary kinetic comparison',
         'Reporter trace includes sensor kinetics and fit selection'),
        (release, 'mean_positive_prob', 'Positive-trial probability', 'LC stimulation',
         'opto_scale_free', 'Probability', 'ROI release metrics', 'Secondary',
         'Sign-based endpoint, not transmitter amount'),
        (release, 'mean_RI', 'Mean response index', 'LC stimulation',
         'opto_raw_scale', 'Response index', 'ROI release metrics', 'Descriptive',
         'Raw fluorescence scales are not molar measurements'),
        (normalised, 'roi_mean_C_over_EC50', 'Sensor-normalised apparent C/EC50',
         'LC stimulation', 'sensor_normalisation', 'Ratio', 'Sensor-normalised session metrics',
         'Sensitivity only', 'Published in vitro constants are not an in vivo calibration'),
        (normalised, 'roi_tau_sensor_ratio', 'Decay tau normalised to sensor off-kinetics',
         'LC stimulation', 'opto_temporal_fit', 'Ratio', 'Sensor-normalised session metrics',
         'Secondary kinetic comparison',
         'Published off-kinetics scale the measured trace; they do not remove expression or fit selection'),
        (spatial, 'spatial_tau_px', 'Chapter 2-style spatial decay constant',
         'LC stimulation', 'opto_spatial', 'pixels', 'Spatial tau metrics, 0-1 s window',
         'Secondary spatial comparison',
         'Native pixel scale; response-qualified and fit-qualified session subset'),
        (spatial, 'response_qualified_fraction',
         'Fraction of LC-axon ROIs meeting the response criterion',
         'LC stimulation', 'opto_spatial', 'Fraction', 'Spatial tau metrics, 0-1 s window',
         'Secondary spatial comparison',
         'Chapter 2 criterion: one-sample greater-than-zero test and mean response index above 0.1'),
    ]
    return add_holm_correction(pd.DataFrame([
        compare_sensors(*specification) for specification in specifications
    ]))


def build_null_tests(release):
    rows = []
    for label, observed, shuffled in [
        ('Top-10% contribution', 'top10_contribution', 'top10_shuf_median'),
        ('Split-half map correlation', 'split_half_r', 'split_half_shuf_median'),
    ]:
        for sensor in SENSOR_ORDER:
            current = release[release['sensor'] == sensor].dropna(subset=[observed, shuffled])
            test = wilcoxon(current[observed], current[shuffled], alternative='greater')
            current = current.copy()
            current['observed_minus_shuffle'] = current[observed] - current[shuffled]
            animal_excess = current.groupby('animal')['observed_minus_shuffle'].median()
            n_positive_animals = int((animal_excess > 0).sum())
            animal_sign_p = float(
                binomtest(n_positive_animals, len(animal_excess), 0.5,
                          alternative='greater').pvalue
            ) if len(animal_excess) else np.nan
            rows.append({
                'metric': label,
                'sensor': sensor,
                'n_sessions': len(current),
                'n_animals': current['animal'].nunique(),
                'observed_median': current[observed].median(),
                'shuffle_median': current[shuffled].median(),
                'median_observed_minus_shuffle': (current[observed] - current[shuffled]).median(),
                'wilcoxon_w': float(test.statistic),
                'p_value_greater': float(test.pvalue),
                'n_positive_animals': n_positive_animals,
                'animal_median_excess_min': float(animal_excess.min()) if len(animal_excess) else np.nan,
                'animal_median_excess_max': float(animal_excess.max()) if len(animal_excess) else np.nan,
                'animal_sign_p_greater': animal_sign_p,
            })
    return pd.DataFrame(rows)


#%% plotting
def result_for(results, metric):
    return results.loc[results['metric'] == metric].iloc[0]


def format_p(value):
    return f'{value:.1e}' if value < 0.001 else f'{value:.3f}'


def dot_panel(ax, frame, metric, ylabel, result, scale=1):
    rng = np.random.default_rng(42)
    for x, sensor in enumerate(SENSOR_ORDER):
        current = frame[frame['sensor'] == sensor].dropna(subset=[metric])
        values = current[metric].to_numpy(dtype=float) * scale
        jitter = rng.uniform(-0.11, 0.11, len(values))
        ax.scatter(np.full(len(values), x) + jitter, values, s=19,
                   color=SENSOR_COLOURS[sensor], edgecolor='white', linewidth=0.35,
                   alpha=0.72, zorder=2)
        animal_medians = current.groupby('animal')[metric].median().to_numpy() * scale
        ax.scatter(np.full(len(animal_medians), x), animal_medians, s=42, marker='o',
                   facecolor='white', edgecolor=SENSOR_COLOURS[sensor], linewidth=1.2,
                   zorder=3)
        q1, median, q3 = np.percentile(values, [25, 50, 75])
        ax.vlines(x, q1, q3, color='black', linewidth=1.2, zorder=4)
        ax.hlines(median, x - 0.14, x + 0.14, color='black', linewidth=1.5, zorder=4)
    ax.set_xticks([0, 1], SENSOR_ORDER)
    ax.set_ylabel(ylabel)
    ax.text(0.03, 0.97,
            f'Session P (descriptive) = {format_p(result.p_value)}\n'
            f'Animal exact P = {format_p(result.animal_exact_p_value)}\n'
            f'Animal rank-biserial = '
            f'{result.animal_rank_biserial_nlight_minus_dlight:.2f}',
            transform=ax.transAxes, va='top', fontsize=7,
            bbox={'facecolor': 'white', 'edgecolor': 'none', 'alpha': 0.78, 'pad': 1.5})
    ax.spines[['top', 'right']].set_visible(False)


def panel_figure(panels, results):
    figure, axes = plt.subplots(1, 4, figsize=(10, 2.5), constrained_layout=True)
    for letter, ax, (frame, metric, ylabel, scale) in zip('ABCD', axes, panels):
        dot_panel(ax, frame, metric, ylabel, result_for(results, metric), scale)
        title = ylabel.split(' (')[0]
        ax.set_title(f'{letter}  {title}', loc='left', fontsize=9)
    for ax in axes[len(panels):]:
        ax.set_visible(False)
    return figure


def save_figure(figure, stem):
    figure.savefig(stem.with_suffix('.png'), dpi=300, bbox_inches='tight')
    figure.savefig(stem.with_suffix('.svg'), bbox_inches='tight')
    plt.close(figure)


#%% provisional EM analysis
def em_profile_table():
    return pd.DataFrame([
        ['P1', 32, 0.1745, 0.100, 106, 607, 675, 52.0, 64.0, 83.5, 43.9, 500.7],
        ['P2', 15, 0.0862, 0.262, 76, 882, 1195, 39.1, 53.0, 73.9, 35.4, 440.0],
        ['P3', 28, 0.5081, 0.407, 100, 197, 332, 57.2, 90.3, 122.0, 50.0, 105.5],
        ['P4', 41, 0.2126, 0.395, 22, 103, 171, 56.7, 99.9, 192.0, 44.4, 83.8],
    ], columns=[
        'profile', 'n_sections', 'axon_volume_um3', 'mitochondrial_fraction',
        'n_vesicles', 'vesicle_density_total_per_um3', 'vesicle_density_cytoplasm_per_um3',
        'observed_nn_median_nm', 'csr_nn_95_low_nm', 'csr_nn_95_high_nm',
        'membrane_distance_median_nm', 'mitochondrion_distance_median_nm',
    ])


def plot_em(em):
    figure, ax = plt.subplots(figsize=(4.2, 2.8), constrained_layout=True)
    y = np.arange(len(em))
    ax.hlines(y, em['csr_nn_95_low_nm'], em['csr_nn_95_high_nm'],
              color='#4c78a8', linewidth=5, alpha=0.45,
              label='Geometry-matched null 95% interval')
    ax.scatter(em['observed_nn_median_nm'], y, color='black', s=35,
               zorder=3, label='Observed vesicles')
    ax.set_yticks(y, em['profile'])
    ax.invert_yaxis()
    ax.set_xlabel('Median nearest-neighbour distance (nm)')
    ax.set_ylabel('Axonal profile')
    ax.spines[['top', 'right']].set_visible(False)
    ax.legend(frameon=False, fontsize=7, loc='upper right')
    return figure


#%% provenance and metric catalogue




def metric_catalogue():
    return pd.DataFrame([
        ['Figure 18c', 'Representative whole-field stimulation traces', 'Descriptive',
         'Selected field', 'Trials in the two selected fields',
         'Mean and SEM across trials within each field'],
        ['Figure 18e', 'Stimulus-evoked decay time constant',
         'Primary kinetic comparison', 'Animal',
         'Fit-qualified stimulation sessions',
         'Session Mann-Whitney U; exact test on animal medians'],
        ['Figure 18f', 'Fit-independent persistence before and after neuropil subtraction',
         'Control', 'Session', 'Included stimulation time courses',
         'Paired Wilcoxon within each reporter cohort'],
        ['Figure 18g', 'Decay tau normalised to sensor off-kinetics',
         'Secondary kinetic comparison', 'Animal',
         'Fit-qualified stimulation sessions',
         'Session Mann-Whitney U; exact test on animal medians'],
        ['Figure 19c', 'Response index by distance from labelled LC axons',
         'Descriptive spatial profile', 'Session',
         'Included stimulation sessions with one-pixel distance rings',
         'Cohort mean and SEM across session summaries'],
        ['Figure 19d', 'Chapter 2-style spatial decay constant',
         'Secondary spatial comparison', 'Animal',
         'Sessions meeting the Chapter 2 spatial-fit criterion',
         'Session Mann-Whitney U; exact test on animal medians'],
        ['Figure 19e', 'Fraction of LC-axon ROIs meeting the response criterion',
         'Secondary spatial comparison', 'Animal',
         'Included stimulation sessions',
         'Session Mann-Whitney U; exact test on animal medians'],
        ['Figure 20c', 'Run-onset-responsive fraction across classifier settings',
         'Classifier sensitivity', 'Animal', 'Valid non-edge axonal candidate regions',
         'Animal traces; session Mann-Whitney U at each setting'],
        ['Figure 20d', 'Run-onset decay time constant', 'Secondary', 'Animal',
         'Fit-qualified responsive candidate regions',
         'Session Mann-Whitney U; exact test on animal medians'],
        ['Figure 20e', 'Run-onset rise time', 'Secondary', 'Animal',
         'Fit-qualified responsive candidate regions',
         'Session Mann-Whitney U; exact test on animal medians'],
        ['Figure 20f', 'Natural-event measurement controls', 'Control', 'Animal',
         'Included natural-event sessions', 'Session distributions'],
        ['Figure 21b', 'Vesicle nearest-neighbour clustering within cytoplasm',
         'Profile-level anatomical result', 'Axonal profile',
         'All annotated vesicle centres; mitochondrial volume excluded',
         'One observed profile statistic against 500 cytoplasmic null patterns'],
        ['Figure 21c', 'Vesicle distance to nearest nearby synapse centre',
         'Profile-level anatomical result', 'Axonal profile',
         'All annotated vesicle centres',
         'One observed profile statistic against 500 cytoplasmic null patterns'],
        ['Figure 21d', 'Variation in synapse-assigned vesicle counts',
         'Profile-level anatomical result', 'Axonal profile',
         'Nearby synapse centres within the fixed inclusion distance',
         'One observed profile statistic against 500 cytoplasmic null patterns'],
        ['Figure 21e', 'Reserved EM image: LC axons near a local CA1 synapse',
         'Composition placeholder', 'Author-supplied image',
         'To be inserted in the final chapter figure', 'No analysis'],
    ], columns=[
        'panel', 'metric', 'priority', 'inference_unit', 'denominator', 'statistical_test',
    ])


#%% run
def main():
    OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)
    FIGURE_ROOT.mkdir(parents=True, exist_ok=True)
    natural, release, temporal, normalised, spatial = load_inputs()
    geometry = load_spatial_geometry_diagnostic()
    results = build_results(natural, release, temporal, normalised, spatial)

    mouse_specs = [
        (natural, 'up_fraction', 'Run-onset-responsive fraction', 'Natural run onset'),
        (natural, 'response_amplitude_all_dff', 'All-region response', 'Natural run onset'),
        (natural, 'response_amplitude_up_dff', 'Responsive-region response', 'Natural run onset'),
        (release, 'top10_minus_shuffle', 'Top-10% excess contribution', 'LC stimulation'),
        (release, 'split_half_minus_shuffle', 'Split-half map stability', 'LC stimulation'),
        (temporal, 'roi_tau_s', 'Stimulus-evoked decay time constant', 'LC stimulation'),
        (normalised, 'roi_tau_sensor_ratio', 'Decay tau normalised to sensor off-kinetics',
         'LC stimulation'),
        (spatial, 'spatial_tau_px', 'Chapter 2-style spatial decay constant',
         'LC stimulation'),
        (spatial, 'response_qualified_fraction',
         'Fraction of LC-axon ROIs meeting the response criterion', 'LC stimulation'),
    ]
    mouse_results = pd.DataFrame([mouse_level_result(*spec) for spec in mouse_specs])
    em = em_profile_table()

    results.to_csv(OUTPUT_ROOT / 'statistical_results.csv', index=False)
    build_null_tests(release).to_csv(OUTPUT_ROOT / 'within_session_null_tests.csv', index=False)
    mouse_results.to_csv(OUTPUT_ROOT / 'mouse_level_sensitivity.csv', index=False)
    em.to_csv(OUTPUT_ROOT / 'em_profile_metrics.csv', index=False)
    metric_catalogue().to_csv(OUTPUT_ROOT / 'metric_catalogue.csv', index=False)

    natural_panels = [
        (natural, 'up_fraction', 'Responsive regions (%)', 100),
        (natural, 'response_amplitude_all_dff', 'All-region response (% ΔF/F)', 100),
        (natural, 'response_amplitude_up_dff', 'Responsive-region response (% ΔF/F)', 100),
        (natural, 'temporal_tau_s', 'Decay time constant (s)', 1),
    ]
    opto_panels = [
        (release, 'top10_minus_shuffle', 'Top-10% excess contribution (%)', 100),
        (release, 'gini_positive_RI', 'Positive-response Gini coefficient', 1),
        (release, 'split_half_minus_shuffle', 'Split-half stability above shuffle', 1),
    ]
    control_panels = [
        (temporal, 'roi_tau_s', 'Decay time constant (s)', 1),
        (release, 'mean_positive_prob', 'Positive-trial probability (%)', 100),
        (release, 'mean_RI', 'Mean response index', 1),
        (normalised, 'roi_mean_C_over_EC50', 'Apparent C/EC50', 1),
    ]
    save_figure(panel_figure(natural_panels, results), FIGURE_ROOT / 'natural_comparison_first')
    save_figure(panel_figure(opto_panels, results), FIGURE_ROOT / 'opto_map_organisation')
    save_figure(panel_figure(control_panels, results),
                FIGURE_ROOT / 'opto_response_kinetics_and_scale_metrics')
    save_figure(plot_em(em), FIGURE_ROOT / 'em_vesicle_nearest_neighbour_analysis')
    print(f'Wrote evidence tables to {OUTPUT_ROOT}')
    print(f'Wrote comparison-first figures to {FIGURE_ROOT}')
    print(
        f'Pixel-geometry diagnostic: '
        f'{geometry.common_zoom3_geometry.astype(bool).sum()}/{len(geometry)} '
        'fields have common 512-by-512, zoom-three geometry; '
        'A152i-20251113-01 and A152i-20251114-01 are excluded. '
        'No physical conversion is made.'
        )


if __name__ == '__main__':
    main()
