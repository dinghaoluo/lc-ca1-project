'''
Created on 25 Jun 2026

parameter sensitivity diagnostic for the LC-DA-CA1 model

@author: Dinghao Luo
'''


#%% imports
from __future__ import annotations

import argparse
import csv
import os
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
import numpy as np

from run_model import (
    CLASS_COLORS,
    DARKGREEN,
    DEFAULT_OUTPUT_DIR,
    GREY,
    Params,
    clean_axis,
    run_bootstrap_suite,
    save_figure_bundle,
    sigmoid,
    summarize_result,
)


#%% constants
OUTPUT_DIR = DEFAULT_OUTPUT_DIR / 'parameter_sensitivity'

plt.rcParams.update(
    {
        'font.family': 'Arial',
        'pdf.fonttype': 42,
        'ps.fonttype': 42,
        'font.size': 9,
        'axes.titlesize': 10,
        'axes.labelsize': 9,
        'xtick.labelsize': 8,
        'ytick.labelsize': 8,
        'legend.fontsize': 8,
    }
)

@dataclass
class SensitivityRow:
    '''one point in the sigmoid parameter grid'''

    da_half_rate: float
    da_rate_slope: float
    base_up_pct: float
    base_down_pct: float
    exp1_delta_up_pct: float
    exp1_delta_down_pct: float
    exp2_delta_up_pct: float
    exp2_delta_down_pct: float
    exp3_delta_up_pct: float
    exp3_delta_down_pct: float
    exp1_up_fr_3_4: float
    exp1_down_fr_0_4: float
    exp2_up_fr_3_4: float
    exp2_down_fr_0_4: float
    exp3_up_fr_3_4: float
    exp3_down_fr_0_4: float
    expected_signs: bool


#%% command line

#%% sigmoid sweep
def run_sensitivity_grid(
        args: argparse.Namespace,
        half_rates: np.ndarray,
        slopes: np.ndarray,
) -> list[SensitivityRow]:
    p_ref = Params(n_bootstrap=args.n_bootstrap, n_cells=args.n_cells, seed_start=args.seed_start)
    rows: list[SensitivityRow] = []
    total = len(half_rates) * len(slopes)
    counter = 0

    for slope in slopes:
        for half_rate in half_rates:
            counter += 1
            print(
                f'[{counter:02d}/{total:02d}] '
                f'da_half_rate={half_rate:.3f}, da_rate_slope={slope:.3f}'
            )
            p = replace(
                p_ref,
                da_half_rate=float(half_rate),
                da_rate_slope=float(slope),
            )
            result = run_bootstrap_suite(p)
            summary = summarize_result(result)
            expected_signs = (
                float(summary['base_up_pct']) > 5.0
                and float(summary['base_down_pct']) > 5.0
                and float(summary['exp1_delta_up_pct']) > 0.0
                and float(summary['exp1_delta_down_pct']) < 0.0
                and float(summary['exp2_delta_up_pct']) > 0.0
                and float(summary['exp3_delta_up_pct']) < 0.0
                and float(summary['exp3_delta_down_pct']) > 0.0
            )
            rows.append(SensitivityRow(
                da_half_rate=float(half_rate),
                da_rate_slope=float(slope),
                base_up_pct=float(summary['base_up_pct']),
                base_down_pct=float(summary['base_down_pct']),
                exp1_delta_up_pct=float(summary['exp1_delta_up_pct']),
                exp1_delta_down_pct=float(summary['exp1_delta_down_pct']),
                exp2_delta_up_pct=float(summary['exp2_delta_up_pct']),
                exp2_delta_down_pct=float(summary['exp2_delta_down_pct']),
                exp3_delta_up_pct=float(summary['exp3_delta_up_pct']),
                exp3_delta_down_pct=float(summary['exp3_delta_down_pct']),
                exp1_up_fr_3_4=float(summary['exp1_up_fr_3_4']),
                exp1_down_fr_0_4=float(summary['exp1_down_fr_0_4']),
                exp2_up_fr_3_4=float(summary['exp2_up_fr_3_4']),
                exp2_down_fr_0_4=float(summary['exp2_down_fr_0_4']),
                exp3_up_fr_3_4=float(summary['exp3_up_fr_3_4']),
                exp3_down_fr_0_4=float(summary['exp3_down_fr_0_4']),
                expected_signs=expected_signs,
            ))

    return rows

def heatmap_axis(
        ax: plt.Axes,
        mat: np.ndarray,
        half_rates: np.ndarray,
        slopes: np.ndarray,
        p_ref: Params,
        title: str,
        cmap: str,
        label: str,
        norm: TwoSlopeNorm | None = None,
) -> None:
    im = ax.imshow(
        mat,
        origin='lower',
        aspect='auto',
        extent=[half_rates[0], half_rates[-1], slopes[0], slopes[-1]],
        cmap=cmap,
        norm=norm,
    )
    ax.scatter(p_ref.da_half_rate, p_ref.da_rate_slope, marker='x', s=55, color='black', linewidths=1.4)
    ax.set_title(title)
    ax.set_xlabel('Inflection point (Hz)')
    ax.set_ylabel('Slope (Hz)')
    ax.set_xticks(half_rates)
    ax.set_yticks(slopes)
    ax.tick_params(axis='x', rotation=45)
    cbar = ax.figure.colorbar(im, ax=ax, shrink=0.82, pad=0.015)
    cbar.set_label(label)
    clean_axis(ax)

def plot_sensitivity_heatmaps(
        rows: list[SensitivityRow],
        half_rates: np.ndarray,
        slopes: np.ndarray,
        output_dir: Path,
) -> None:
    p_ref = Params()
    mats = {
        'base_up_pct': np.asarray([row.base_up_pct for row in rows], dtype=float).reshape(len(slopes), len(half_rates)),
        'base_down_pct': np.asarray([row.base_down_pct for row in rows], dtype=float).reshape(len(slopes), len(half_rates)),
        'exp1_delta_up_pct': np.asarray([row.exp1_delta_up_pct for row in rows], dtype=float).reshape(len(slopes), len(half_rates)),
        'exp1_delta_down_pct': np.asarray([row.exp1_delta_down_pct for row in rows], dtype=float).reshape(len(slopes), len(half_rates)),
        'exp2_delta_up_pct': np.asarray([row.exp2_delta_up_pct for row in rows], dtype=float).reshape(len(slopes), len(half_rates)),
        'exp2_delta_down_pct': np.asarray([row.exp2_delta_down_pct for row in rows], dtype=float).reshape(len(slopes), len(half_rates)),
        'exp3_delta_up_pct': np.asarray([row.exp3_delta_up_pct for row in rows], dtype=float).reshape(len(slopes), len(half_rates)),
        'exp3_delta_down_pct': np.asarray([row.exp3_delta_down_pct for row in rows], dtype=float).reshape(len(slopes), len(half_rates)),
    }
    abs_max = np.nanmax([
        np.nanmax(np.abs(mat))
        for mat in (
            mats['exp1_delta_up_pct'],
            mats['exp1_delta_down_pct'],
            mats['exp2_delta_up_pct'],
            mats['exp2_delta_down_pct'],
            mats['exp3_delta_up_pct'],
            mats['exp3_delta_down_pct'],
        )
    ])
    abs_max = max(float(abs_max), 1.0)
    prop_norm = TwoSlopeNorm(vmin=-abs_max, vcenter=0.0, vmax=abs_max)

    fig = plt.figure(figsize=(11.2, 8.5), constrained_layout=True)
    fig.set_facecolor('white')
    gs = fig.add_gridspec(3, 3, hspace=0.22, wspace=0.25)

    ax = fig.add_subplot(gs[0, 0])
    r_grid = np.linspace(0.0, 4.0, 240)
    for half_rate in half_rates:
        for slope in slopes:
            ax.plot(r_grid, sigmoid(r_grid, half_rate, slope), color=GREY, alpha=0.16, linewidth=0.8)
    ax.plot(
        r_grid,
        sigmoid(r_grid, p_ref.da_half_rate, p_ref.da_rate_slope),
        color=DARKGREEN,
        linewidth=2.4,
        label='selected',
    )
    ax.axvline(p_ref.da_half_rate, linestyle='--', color='0.45', linewidth=1.0)
    ax.set_xlim([0.0, 4.0])
    ax.set_ylim([0.0, 1.05])
    ax.set_title('Sigmoid sweep')
    ax.set_xlabel('Previous firing rate (Hz)')
    ax.set_ylabel('Extra DA weight')
    ax.legend(frameon=False, loc='lower right')
    clean_axis(ax)

    heatmap_axis(
        fig.add_subplot(gs[0, 1]),
        mats['base_up_pct'],
        half_rates,
        slopes,
        p_ref,
        'Baseline PyrUp',
        'Reds',
        'Cells (%)',
    )
    heatmap_axis(
        fig.add_subplot(gs[0, 2]),
        mats['base_down_pct'],
        half_rates,
        slopes,
        p_ref,
        'Baseline PyrDown',
        'Purples',
        'Cells (%)',
    )
    heatmap_axis(
        fig.add_subplot(gs[1, 0]),
        mats['exp1_delta_up_pct'],
        half_rates,
        slopes,
        p_ref,
        'Exp. 1 PyrUp, Stim. - Ctrl.',
        'RdBu_r',
        'Delta cells (pp)',
        norm=prop_norm,
    )
    heatmap_axis(
        fig.add_subplot(gs[1, 1]),
        mats['exp1_delta_down_pct'],
        half_rates,
        slopes,
        p_ref,
        'Exp. 1 PyrDown, Stim. - Ctrl.',
        'RdBu_r',
        'Delta cells (pp)',
        norm=prop_norm,
    )
    heatmap_axis(
        fig.add_subplot(gs[1, 2]),
        mats['exp2_delta_up_pct'],
        half_rates,
        slopes,
        p_ref,
        'Exp. 2 PyrUp, DA - non-DA',
        'RdBu_r',
        'Delta cells (pp)',
        norm=prop_norm,
    )
    heatmap_axis(
        fig.add_subplot(gs[2, 0]),
        mats['exp2_delta_down_pct'],
        half_rates,
        slopes,
        p_ref,
        'Exp. 2 PyrDown, DA - non-DA',
        'RdBu_r',
        'Delta cells (pp)',
        norm=prop_norm,
    )
    heatmap_axis(
        fig.add_subplot(gs[2, 1]),
        mats['exp3_delta_up_pct'],
        half_rates,
        slopes,
        p_ref,
        'Exp. 3 PyrUp, Blockade - Baseline',
        'RdBu_r',
        'Delta cells (pp)',
        norm=prop_norm,
    )
    heatmap_axis(
        fig.add_subplot(gs[2, 2]),
        mats['exp3_delta_down_pct'],
        half_rates,
        slopes,
        p_ref,
        'Exp. 3 PyrDown, Blockade - Baseline',
        'RdBu_r',
        'Delta cells (pp)',
        norm=prop_norm,
    )

    save_figure_bundle(fig, output_dir, 'model_parameter_sensitivity_heatmaps')
    plt.close(fig)

def print_grid_summary(rows: list[SensitivityRow], csv_path: Path, output_dir: Path) -> None:
    n_expected = sum(row.expected_signs for row in rows)
    n_total = len(rows)
    print()
    print('LC-DA-CA1 sigmoid parameter sensitivity')
    print(f'  grid points with expected sign pattern: {n_expected}/{n_total}')
    print(f'  baseline PyrUp range: {min(row.base_up_pct for row in rows):.2f}-{max(row.base_up_pct for row in rows):.2f}%')
    print(f'  baseline PyrDown range: {min(row.base_down_pct for row in rows):.2f}-{max(row.base_down_pct for row in rows):.2f}%')
    print(
        '  Exp. 1 PyrUp Stim.-Ctrl. range: '
        f'{min(row.exp1_delta_up_pct for row in rows):+.2f} to {max(row.exp1_delta_up_pct for row in rows):+.2f} pp'
    )
    print(
        '  Exp. 3 PyrDown Blockade-Baseline range: '
        f'{min(row.exp3_delta_down_pct for row in rows):+.2f} to {max(row.exp3_delta_down_pct for row in rows):+.2f} pp'
    )
    if n_expected < n_total:
        print('  grid points outside the expected sign pattern:')
        for row in rows:
            if not row.expected_signs:
                print(f'    h={row.da_half_rate:.3f}, s={row.da_rate_slope:.3f}')
    print(f'  table: {csv_path}')
    print(f'  figures: {output_dir}')

#%% broader parameter diagnostics
@dataclass(frozen=True)
class MetricSpec:
    key: str
    label: str
    unit: str

@dataclass(frozen=True)
class ParameterSweep:
    name: str
    lower: float
    upper: float
    label: str
    family: str

METRIC_SPECS = (
    MetricSpec('base_up_pct', 'baseline PyrUp', '%'),
    MetricSpec('base_down_pct', 'baseline PyrDown', '%'),
    MetricSpec('exp1_delta_up_pct', 'Exp. 1 PyrUp proportion, Stim. - Ctrl.', 'pp'),
    MetricSpec('exp1_delta_down_pct', 'Exp. 1 PyrDown proportion, Stim. - Ctrl.', 'pp'),
    MetricSpec('exp2_delta_up_pct', 'Exp. 2 PyrUp proportion, DA - non-DA', 'pp'),
    MetricSpec('exp2_delta_down_pct', 'Exp. 2 PyrDown proportion, DA - non-DA', 'pp'),
    MetricSpec('exp3_delta_up_pct', 'Exp. 3 PyrUp proportion, Blockade - Baseline', 'pp'),
    MetricSpec('exp3_delta_down_pct', 'Exp. 3 PyrDown proportion, Blockade - Baseline', 'pp'),
    MetricSpec('exp1_up_fr_0_4', 'Exp. 1 PyrUp FR, Stim. - Ctrl., 0-4 s', 'Hz'),
    MetricSpec('exp1_down_fr_0_4', 'Exp. 1 PyrDown FR, Stim. - Ctrl., 0-4 s', 'Hz'),
    MetricSpec('exp2_up_fr_0_4', 'Exp. 2 PyrUp FR, DA - non-DA, 0-4 s', 'Hz'),
    MetricSpec('exp2_down_fr_0_4', 'Exp. 2 PyrDown FR, DA - non-DA, 0-4 s', 'Hz'),
    MetricSpec('exp3_up_fr_0_4', 'Exp. 3 PyrUp FR, Baseline - Blockade, 0-4 s', 'Hz'),
    MetricSpec('exp3_down_fr_0_4', 'Exp. 3 PyrDown FR, Baseline - Blockade, 0-4 s', 'Hz'),
    MetricSpec('exp1_up_fr_3p5_4', 'Exp. 1 PyrUp FR, Stim. - Ctrl., 3.5-4 s', 'Hz'),
    MetricSpec('exp1_down_fr_3p5_4', 'Exp. 1 PyrDown FR, Stim. - Ctrl., 3.5-4 s', 'Hz'),
    MetricSpec('exp2_up_fr_3p5_4', 'Exp. 2 PyrUp FR, DA - non-DA, 3.5-4 s', 'Hz'),
    MetricSpec('exp2_down_fr_3p5_4', 'Exp. 2 PyrDown FR, DA - non-DA, 3.5-4 s', 'Hz'),
    MetricSpec('exp3_up_fr_3p5_4', 'Exp. 3 PyrUp FR, Baseline - Blockade, 3.5-4 s', 'Hz'),
    MetricSpec('exp3_down_fr_3p5_4', 'Exp. 3 PyrDown FR, Baseline - Blockade, 3.5-4 s', 'Hz'),
    MetricSpec('exp1_up_retention', 'Exp. 1 PyrUp tail / 0-4 s', 'ratio'),
    MetricSpec('exp2_up_retention', 'Exp. 2 PyrUp tail / 0-4 s', 'ratio'),
    MetricSpec('exp3_up_retention', 'Exp. 3 PyrUp tail / 0-4 s', 'ratio'),
)

METRIC_LABELS = {spec.key: spec.label for spec in METRIC_SPECS}
METRIC_UNITS = {spec.key: spec.unit for spec in METRIC_SPECS}
PRIMARY_METRICS = tuple(spec.key for spec in METRIC_SPECS)

PARAMETER_SWEEPS = (
    ParameterSweep('baseline_mean', 0.45, 1.20, 'baseline mean', 'population'),
    ParameterSweep('baseline_sd', 0.20, 0.75, 'baseline SD', 'population'),
    ParameterSweep('wR_mean', 0.25, 1.60, 'run weight mean', 'population'),
    ParameterSweep('wR_sd', 0.40, 2.30, 'run weight SD', 'population'),
    ParameterSweep('wW_mean', 0.35, 1.70, 'task weight mean', 'population'),
    ParameterSweep('wW_sd', 0.35, 1.80, 'task weight SD', 'population'),
    ParameterSweep('frac_da_targ', 0.05, 0.65, 'DA-targeted fraction', 'DA targeting'),
    ParameterSweep('da_half_rate', 0.80, 3.80, 'DA sigmoid inflection', 'DA targeting'),
    ParameterSweep('da_rate_slope', 0.02, 0.32, 'DA sigmoid slope', 'DA targeting'),
    ParameterSweep('wDA_global', -0.20, 0.45, 'global DA weight', 'DA targeting'),
    ParameterSweep('lc_activation_fold', 0.60, 4.20, 'LC activation fold', 'DA drive'),
    ParameterSweep('da_block_scale', 0.00, 1.10, 'DA blockade scale', 'DA drive'),
    ParameterSweep('intrinsic_tau_mean', 0.20, 1.70, 'intrinsic tau mean', 'cell dynamics'),
    ParameterSweep('baseline_tau_coupling', 0.00, 0.45, 'baseline-tau coupling', 'cell dynamics'),
    ParameterSweep('intrinsic_tau_sd', 0.00, 0.16, 'intrinsic tau SD', 'cell dynamics'),
    ParameterSweep('softplus_beta', 0.80, 4.20, 'softplus beta', 'output'),
    ParameterSweep('max_rate', 6.00, 35.00, 'maximum rate', 'output'),
    ParameterSweep('run_on_mid', -0.30, 0.30, 'run onset midpoint', 'run drive'),
    ParameterSweep('run_off_mid', 0.90, 2.60, 'run offset midpoint', 'run drive'),
    ParameterSweep('run_rise_scale', 0.03, 0.25, 'run rise scale', 'run drive'),
    ParameterSweep('run_fall_scale', 0.20, 1.20, 'run fall scale', 'run drive'),
    ParameterSweep('reward_on_mid', -0.10, 0.55, 'task onset midpoint', 'task drive'),
    ParameterSweep('reward_off_mid', 0.60, 2.30, 'task offset midpoint', 'task drive'),
    ParameterSweep('reward_rise_scale', 0.02, 0.18, 'task rise scale', 'task drive'),
    ParameterSweep('reward_fall_scale', 0.20, 1.40, 'task fall scale', 'task drive'),
    ParameterSweep('lc_mu', -0.30, 0.40, 'LC pulse centre', 'DA waveform'),
    ParameterSweep('lc_sigma', 0.08, 0.45, 'LC pulse width', 'DA waveform'),
    ParameterSweep('da_kernel_tau', 0.50, 3.40, 'DA kernel tau', 'DA waveform'),
    ParameterSweep('da_ca1_delay', 0.00, 0.95, 'DA CA1 delay', 'DA waveform'),
    ParameterSweep('up_thresh', 1.05, 2.10, 'PyrUp threshold', 'classification'),
    ParameterSweep('down_thresh', 0.40, 0.95, 'PyrDown threshold', 'classification'),
)

def collect_extended_metrics(result: dict[str, Any]) -> dict[str, float]:
    summary = summarize_result(result)
    metrics = {
        'base_up_pct': float(summary['base_up_pct']),
        'base_down_pct': float(summary['base_down_pct']),
        'exp1_delta_up_pct': float(summary['exp1_delta_up_pct']),
        'exp1_delta_down_pct': float(summary['exp1_delta_down_pct']),
        'exp2_delta_up_pct': float(summary['exp2_delta_up_pct']),
        'exp2_delta_down_pct': float(summary['exp2_delta_down_pct']),
        'exp3_delta_up_pct': float(summary['exp3_delta_up_pct']),
        'exp3_delta_down_pct': float(summary['exp3_delta_down_pct']),
        'exp1_up_fr_0_4': float(summary['exp1_up_fr_0_4_delta_test']),
        'exp1_down_fr_0_4': float(summary['exp1_down_fr_0_4_delta_test']),
        'exp2_up_fr_0_4': float(summary['exp2_up_fr_0_4_delta_test']),
        'exp2_down_fr_0_4': float(summary['exp2_down_fr_0_4_delta_test']),
        'exp3_up_fr_0_4': float(summary['exp3_up_fr_0_4_delta_test']),
        'exp3_down_fr_0_4': float(summary['exp3_down_fr_0_4_delta_test']),
        'exp1_up_fr_3p5_4': float(summary['exp1_up_fr_3p5_4_delta_test']),
        'exp1_down_fr_3p5_4': float(summary['exp1_down_fr_3p5_4_delta_test']),
        'exp2_up_fr_3p5_4': float(summary['exp2_up_fr_3p5_4_delta_test']),
        'exp2_down_fr_3p5_4': float(summary['exp2_down_fr_3p5_4_delta_test']),
        'exp3_up_fr_3p5_4': float(summary['exp3_up_fr_3p5_4_delta_test']),
        'exp3_down_fr_3p5_4': float(summary['exp3_down_fr_3p5_4_delta_test']),
    }
    metrics['exp1_up_retention'] = metrics['exp1_up_fr_3p5_4'] / metrics['exp1_up_fr_0_4']
    metrics['exp2_up_retention'] = metrics['exp2_up_fr_3p5_4'] / metrics['exp2_up_fr_0_4']
    metrics['exp3_up_retention'] = metrics['exp3_up_fr_3p5_4'] / metrics['exp3_up_fr_0_4']
    return metrics

def core_failures(metrics: dict[str, float]) -> list[str]:
    core_metrics = (
        'base_up_pct',
        'base_down_pct',
        'exp1_delta_up_pct',
        'exp1_delta_down_pct',
        'exp2_delta_up_pct',
        'exp3_delta_up_pct',
        'exp3_delta_down_pct',
        'exp1_up_fr_0_4',
        'exp2_up_fr_0_4',
        'exp3_up_fr_0_4',
        'exp1_up_fr_3p5_4',
        'exp2_up_fr_3p5_4',
        'exp3_up_fr_3p5_4',
    )
    invalid = [
        key for key in core_metrics
        if key not in metrics or not np.isfinite(metrics[key])
    ]
    if invalid:
        return [f'{METRIC_LABELS[key]} is missing or non-finite' for key in invalid]

    failures: list[str] = []
    if metrics['base_up_pct'] < 10.0:
        failures.append('baseline PyrUp <10%')
    if metrics['base_down_pct'] < 10.0:
        failures.append('baseline PyrDown <10%')
    if metrics['exp1_delta_up_pct'] <= 0.0:
        failures.append('Exp. 1 PyrUp proportion does not increase')
    if metrics['exp1_delta_down_pct'] >= 0.0:
        failures.append('Exp. 1 PyrDown proportion does not decrease')
    if metrics['exp2_delta_up_pct'] <= 0.0:
        failures.append('Exp. 2 DA-targeted PyrUp enrichment lost')
    if metrics['exp3_delta_up_pct'] >= 0.0:
        failures.append('Exp. 3 PyrUp blockade reduction lost')
    if metrics['exp3_delta_down_pct'] <= 0.0:
        failures.append('Exp. 3 PyrDown blockade increase lost')
    if metrics['exp1_up_fr_0_4'] <= 0.0:
        failures.append('Exp. 1 PyrUp FR effect <=0')
    if metrics['exp2_up_fr_0_4'] <= 0.0:
        failures.append('Exp. 2 PyrUp FR effect <=0')
    if metrics['exp3_up_fr_0_4'] <= 0.0:
        failures.append('Exp. 3 PyrUp FR effect <=0')
    if metrics['exp1_up_fr_3p5_4'] <= 0.0:
        failures.append('Exp. 1 PyrUp tail effect <=0')
    if metrics['exp2_up_fr_3p5_4'] <= 0.0:
        failures.append('Exp. 2 PyrUp tail effect <=0')
    if metrics['exp3_up_fr_3p5_4'] <= 0.0:
        failures.append('Exp. 3 PyrUp tail effect <=0')
    return failures

def relative_delta_pct(value: float, default: float) -> float:
    if default == 0:
        return np.nan
    return 100.0 * (value - default) / abs(default)

def write_rows(rows: list[dict[str, Any]], output_dir: Path, filename: str) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / filename
    with path.open('w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
    return path

def read_rows(path: Path) -> list[dict[str, Any]]:
    with path.open('r', newline='', encoding='utf-8') as f:
        rows: list[dict[str, Any]] = []
        for row in csv.DictReader(f):
            parsed: dict[str, Any] = {}
            for key, value in row.items():
                if value == '':
                    parsed[key] = ''
                elif value == 'True':
                    parsed[key] = True
                elif value == 'False':
                    parsed[key] = False
                else:
                    try:
                        parsed[key] = float(value)
                    except ValueError:
                        parsed[key] = value
            rows.append(parsed)
        return rows

def redraw_from_tables(output_dir: Path) -> None:
    sigmoid_rows = []
    for row in read_rows(output_dir / 'model_parameter_sensitivity_summary.csv'):
        row['expected_signs'] = bool(row['expected_signs'])
        sigmoid_rows.append(SensitivityRow(**row))
    half_rates = np.asarray(sorted({row.da_half_rate for row in sigmoid_rows}), dtype=float)
    slopes = np.asarray(sorted({row.da_rate_slope for row in sigmoid_rows}), dtype=float)
    default_rows = read_rows(output_dir / 'model_parameter_sensitivity_default_metrics.csv')
    default_metrics = {str(row['metric']): float(row['value']) for row in default_rows}
    sweep_rows = read_rows(output_dir / 'model_parameter_sensitivity_parameter_sweep.csv')
    breakpoint_rows = read_rows(output_dir / 'model_parameter_sensitivity_breakpoints.csv')
    importance_rows = read_rows(output_dir / 'model_parameter_sensitivity_importance.csv')

    plot_sensitivity_heatmaps(sigmoid_rows, half_rates, slopes, output_dir)
    plot_parameter_importance(importance_rows, output_dir)
    plot_breakpoints(breakpoint_rows, output_dir)
    plot_top_metric_curves(sweep_rows, importance_rows, default_metrics, output_dir)
    write_text_report(output_dir, default_metrics, importance_rows, breakpoint_rows, sigmoid_rows)

def run_parameter_sweeps(args: argparse.Namespace, default_metrics: dict[str, float]) -> list[dict[str, Any]]:
    p_ref = Params(n_bootstrap=args.n_bootstrap, n_cells=args.n_cells, seed_start=args.seed_start)
    rows: list[dict[str, Any]] = []
    sweep_values = []
    for sweep in PARAMETER_SWEEPS:
        default_value = float(getattr(p_ref, sweep.name))
        values = np.linspace(sweep.lower, sweep.upper, args.param_grid_size)
        if not np.any(np.isclose(values, default_value, rtol=1e-9, atol=1e-9)):
            values = np.sort(np.append(values, default_value))
        sweep_values.append((sweep, default_value, values))

    total = sum(len(values) for _, _, values in sweep_values)
    counter = 0

    for sweep, default_value, values in sweep_values:
        for value in values:
            counter += 1
            value = float(value)
            print(f'[parameter {counter:03d}/{total:03d}] {sweep.name}={value:.5g}')
            if np.isclose(value, default_value, rtol=1e-10, atol=1e-10):
                metrics = dict(default_metrics)
            else:
                p = replace(p_ref, **{sweep.name: value})
                metrics = collect_extended_metrics(run_bootstrap_suite(p))
            failures = core_failures(metrics)
            row: dict[str, Any] = {
                'parameter': sweep.name,
                'label': sweep.label,
                'family': sweep.family,
                'value': value,
                'default_value': default_value,
                'relative_delta_pct': relative_delta_pct(value, default_value),
                'is_default': bool(np.isclose(value, default_value, rtol=1e-10, atol=1e-10)),
                'passes_core': len(failures) == 0,
                'failed_criteria': '; '.join(failures),
            }
            row.update(metrics)
            rows.append(row)
    return rows

def parameter_breakpoints(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    out: list[dict[str, Any]] = []
    for sweep in PARAMETER_SWEEPS:
        p_rows = sorted([row for row in rows if row['parameter'] == sweep.name], key=lambda row: row['value'])
        default_value = float(p_rows[0]['default_value'])
        lower_rows = sorted(
            [row for row in p_rows if float(row['value']) < default_value],
            key=lambda row: abs(float(row['value']) - default_value),
        )
        upper_rows = sorted(
            [row for row in p_rows if float(row['value']) > default_value],
            key=lambda row: abs(float(row['value']) - default_value),
        )
        lower_fail = next((row for row in lower_rows if not row['passes_core']), None)
        upper_fail = next((row for row in upper_rows if not row['passes_core']), None)
        lower_edge = lower_rows[-1] if lower_rows else None
        upper_edge = upper_rows[-1] if upper_rows else None
        lower_pct = relative_delta_pct(lower_fail['value'], default_value) if lower_fail else np.nan
        upper_pct = relative_delta_pct(upper_fail['value'], default_value) if upper_fail else np.nan
        nearest_vals = np.asarray([abs(lower_pct), abs(upper_pct)], dtype=float)
        nearest = float(np.nanmin(nearest_vals)) if np.any(np.isfinite(nearest_vals)) else np.nan
        out.append(
            {
                'parameter': sweep.name,
                'label': sweep.label,
                'family': sweep.family,
                'default_value': default_value,
                'lower_failure_value': lower_fail['value'] if lower_fail else np.nan,
                'lower_failure_delta_pct': lower_pct,
                'lower_failure_criteria': lower_fail['failed_criteria'] if lower_fail else '',
                'lower_sampled_edge_value': lower_edge['value'] if lower_edge else np.nan,
                'lower_sampled_edge_delta_pct': relative_delta_pct(lower_edge['value'], default_value) if lower_edge else np.nan,
                'upper_failure_value': upper_fail['value'] if upper_fail else np.nan,
                'upper_failure_delta_pct': upper_pct,
                'upper_failure_criteria': upper_fail['failed_criteria'] if upper_fail else '',
                'upper_sampled_edge_value': upper_edge['value'] if upper_edge else np.nan,
                'upper_sampled_edge_delta_pct': relative_delta_pct(upper_edge['value'], default_value) if upper_edge else np.nan,
                'nearest_failure_delta_pct': nearest,
                'breaks_in_sampled_range': bool(lower_fail is not None or upper_fail is not None),
            }
        )
    return out

def parameter_importance(
        rows: list[dict[str, Any]],
        default_metrics: dict[str, float],
        breakpoints: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    scales = {}
    for metric in PRIMARY_METRICS:
        values = np.asarray([float(row[metric]) for row in rows], dtype=float)
        q1, q3 = np.percentile(values[np.isfinite(values)], [25, 75])
        scales[metric] = (q3 - q1) / 1.349

    bp_lookup = {row['parameter']: row for row in breakpoints}
    out: list[dict[str, Any]] = []
    for sweep in PARAMETER_SWEEPS:
        p_rows = [row for row in rows if row['parameter'] == sweep.name]
        metric_scores = []
        for metric in PRIMARY_METRICS:
            vals = np.asarray([float(row[metric]) for row in p_rows], dtype=float)
            deltas = vals - float(default_metrics[metric])
            idx = int(np.nanargmax(np.abs(deltas)))
            max_delta = float(deltas[idx])
            score = abs(max_delta) / scales[metric]
            value_at_max = float(p_rows[idx]['value'])
            metric_scores.append((metric, score, max_delta, value_at_max))
        largest_metric, largest_score, largest_delta, largest_value = max(
            metric_scores,
            key=lambda item: -np.inf if not np.isfinite(item[1]) else item[1],
        )
        bp = bp_lookup[sweep.name]
        out.append(
            {
                'parameter': sweep.name,
                'label': sweep.label,
                'family': sweep.family,
                'default_value': float(getattr(Params(), sweep.name)),
                'mean_abs_scaled_change': float(np.nanmean([item[1] for item in metric_scores])),
                'max_abs_scaled_change': float(largest_score),
                'largest_metric': largest_metric,
                'largest_metric_label': METRIC_LABELS[largest_metric],
                'largest_metric_delta': largest_delta,
                'largest_metric_value_at_max': largest_value,
                'nearest_failure_delta_pct': bp['nearest_failure_delta_pct'],
                'breaks_in_sampled_range': bp['breaks_in_sampled_range'],
            }
        )
    return sorted(out, key=lambda row: row['max_abs_scaled_change'], reverse=True)

def plot_parameter_importance(importance_rows: list[dict[str, Any]], output_dir: Path, top_n: int = 16) -> None:
    rows = list(reversed(importance_rows[:top_n]))
    fig, ax = plt.subplots(figsize=(7.0, 7.0), constrained_layout=True)
    fig.set_facecolor('white')
    y = np.arange(len(rows))
    scores = [row['max_abs_scaled_change'] for row in rows]
    ax.barh(y, scores, color=DARKGREEN, alpha=0.72, edgecolor='black', linewidth=0.7)
    ax.set_yticks(y)
    ax.set_yticklabels([row['label'] for row in rows])
    ax.tick_params(axis='y', length=0)
    ax.set_xlabel('Sensitivity score (higher = more sensitive)')
    ax.set_title('One-parameter sensitivity, top 16')
    xmax = max(scores)
    ax.set_xlim([0.0, xmax * 1.08])
    ax.set_ylim([-0.6, len(rows) - 0.4])
    ax.grid(axis='x', color='0.88', linewidth=0.6)
    ax.set_axisbelow(True)
    ax.set_facecolor('white')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    fig.text(
        0.01,
        0.01,
        'Score = largest absolute change in any diagnostic readout, divided by the sampled scale for that readout; each bar varies one parameter while the rest stay fixed.',
        ha='left',
        va='bottom',
        fontsize=7,
        color='0.25',
    )
    save_figure_bundle(fig, output_dir, 'model_parameter_sensitivity_parameter_importance')
    plt.close(fig)

def plot_breakpoints(breakpoint_rows: list[dict[str, Any]], output_dir: Path) -> None:
    rows = sorted(
        breakpoint_rows,
        key=lambda row: row['nearest_failure_delta_pct'] if np.isfinite(row['nearest_failure_delta_pct']) else np.inf,
    )
    fig, ax = plt.subplots(figsize=(8.5, 8.8), constrained_layout=True)
    fig.set_facecolor('white')
    y = np.arange(len(rows))
    ax.axvline(0.0, color='black', linewidth=0.8)
    for yy, row in zip(y, rows):
        lower_x = row['lower_failure_delta_pct'] if np.isfinite(row['lower_failure_delta_pct']) else row['lower_sampled_edge_delta_pct']
        upper_x = row['upper_failure_delta_pct'] if np.isfinite(row['upper_failure_delta_pct']) else row['upper_sampled_edge_delta_pct']
        lower_failed = np.isfinite(row['lower_failure_delta_pct'])
        upper_failed = np.isfinite(row['upper_failure_delta_pct'])
        if np.isfinite(lower_x):
            ax.scatter(
                lower_x,
                yy,
                s=34,
                facecolors=CLASS_COLORS['is_down'] if lower_failed else 'white',
                edgecolors=CLASS_COLORS['is_down'] if lower_failed else '0.45',
                linewidths=1.0,
                zorder=3,
            )
        if np.isfinite(upper_x):
            ax.scatter(
                upper_x,
                yy,
                s=34,
                facecolors=CLASS_COLORS['is_up'] if upper_failed else 'white',
                edgecolors=CLASS_COLORS['is_up'] if upper_failed else '0.45',
                linewidths=1.0,
                zorder=3,
            )
        if np.isfinite(lower_x) and np.isfinite(upper_x):
            ax.plot([lower_x, upper_x], [yy, yy], color='0.82', linewidth=0.8, zorder=1)
    ax.set_yticks(y)
    ax.set_yticklabels([row['label'] for row in rows])
    ax.set_xlabel('Parameter change from selected value (%)')
    ax.set_title('First sampled parameter value that breaks core criteria')
    ax.text(0.01, 0.02, 'filled = first failed point; open = sampled edge without failure', transform=ax.transAxes, fontsize=7)
    ax.set_ylim([-0.6, len(rows) - 0.4])
    ax.set_facecolor('white')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    save_figure_bundle(fig, output_dir, 'model_parameter_sensitivity_breakpoints')
    plt.close(fig)

def plot_top_metric_curves(
        sweep_rows: list[dict[str, Any]],
        importance_rows: list[dict[str, Any]],
        default_metrics: dict[str, float],
        output_dir: Path,
        top_n: int = 8,
) -> None:
    rows = importance_rows[:top_n]
    n_cols = 4
    n_rows = int(np.ceil(len(rows) / n_cols))
    fig, axs = plt.subplots(n_rows, n_cols, figsize=(11.5, 2.6 * n_rows), constrained_layout=True)
    fig.set_facecolor('white')
    axs = np.asarray(axs).reshape(-1)
    for ax, imp in zip(axs, rows):
        p_rows = sorted([row for row in sweep_rows if row['parameter'] == imp['parameter']], key=lambda row: row['value'])
        metric = imp['largest_metric']
        x = np.asarray([row['relative_delta_pct'] for row in p_rows], dtype=float)
        y = np.asarray([row[metric] for row in p_rows], dtype=float)
        if np.all(~np.isfinite(x)):
            x = np.asarray([row['value'] - row['default_value'] for row in p_rows], dtype=float)
            xlabel = 'Absolute parameter change'
        else:
            xlabel = 'Change from selected value (%)'
        ax.plot(x, y, color=DARKGREEN, linewidth=1.8)
        ax.scatter(x, y, s=18, color=DARKGREEN, edgecolor='none', alpha=0.75)
        ax.axhline(default_metrics[metric], color='0.45', linestyle='--', linewidth=0.9)
        ax.axvline(0.0, color='black', linewidth=0.8)
        ax.set_title(imp['label'])
        ax.set_xlabel(xlabel)
        ax.set_ylabel(f'{METRIC_LABELS[metric]} ({METRIC_UNITS[metric]})')
        clean_axis(ax)
    for ax in axs[len(rows):]:
        ax.set_axis_off()
    save_figure_bundle(fig, output_dir, 'model_parameter_sensitivity_top_metric_curves')
    plt.close(fig)

def write_text_report(
        output_dir: Path,
        default_metrics: dict[str, float],
        importance_rows: list[dict[str, Any]],
        breakpoint_rows: list[dict[str, Any]],
        sigmoid_rows: list[SensitivityRow],
) -> Path:
    path = output_dir / 'model_parameter_sensitivity_report.txt'
    sigmoid_pass = sum(row.expected_signs for row in sigmoid_rows)
    fragile = [row for row in sorted(
        breakpoint_rows,
        key=lambda row: row['nearest_failure_delta_pct'] if np.isfinite(row['nearest_failure_delta_pct']) else np.inf,
    ) if np.isfinite(row['nearest_failure_delta_pct'])]
    lines = [
        'LC-DA-CA1 model parameter sensitivity',
        '',
        f'default baseline PyrUp/PyrDown = {default_metrics["base_up_pct"]:.2f}% / {default_metrics["base_down_pct"]:.2f}%',
        f'default Exp. 1 proportion deltas = {default_metrics["exp1_delta_up_pct"]:+.2f} / {default_metrics["exp1_delta_down_pct"]:+.2f} pp',
        f'default Exp. 2 PyrUp DA enrichment = {default_metrics["exp2_delta_up_pct"]:+.2f} pp',
        f'default Exp. 3 proportion deltas = {default_metrics["exp3_delta_up_pct"]:+.2f} / {default_metrics["exp3_delta_down_pct"]:+.2f} pp',
        f'sigmoid grid pass count = {sigmoid_pass}/{len(sigmoid_rows)}',
        '',
        'top one-parameter sensitivities',
    ]
    for row in importance_rows[:10]:
        lines.append(f'{row["label"]}: score {row["max_abs_scaled_change"]:.2f}, largest metric = {row["largest_metric_label"]}')
    lines.extend(['', 'nearest sampled failures'])
    for row in fragile[:10]:
        lines.append(f'{row["label"]}: {row["nearest_failure_delta_pct"]:+.1f}% from selected value')
    if not fragile:
        lines.append('no sampled parameter value broke the core criteria')
    path.write_text('\n'.join(lines) + '\n', encoding='utf-8')
    return path

def print_extended_summary(default_metrics: dict[str, float], importance_rows: list[dict[str, Any]], breakpoint_rows: list[dict[str, Any]], output_dir: Path) -> None:
    print('  top one-parameter sensitivities:')
    for row in importance_rows[:8]:
        print(f'    {row["label"]}: {row["max_abs_scaled_change"]:.2f} ({row["largest_metric_label"]})')
    print('  nearest sampled failures:')
    fragile = [row for row in sorted(
        breakpoint_rows,
        key=lambda row: row['nearest_failure_delta_pct'] if np.isfinite(row['nearest_failure_delta_pct']) else np.inf,
    ) if np.isfinite(row['nearest_failure_delta_pct'])]
    for row in fragile[:8]:
        print(f'    {row["label"]}: {row["nearest_failure_delta_pct"]:+.1f}%')
    print(f'  outputs: {output_dir}')

#%% main
def main() -> None:
    parser = argparse.ArgumentParser(
        description='Sweep the activity-dependent DA sigmoid used by run_model.py.'
    )
    parser.add_argument(
        '--output-dir',
        type=Path,
        default=Path(os.getenv('LC_DA_CA1_SENSITIVITY_DIR', OUTPUT_DIR)),
        help='Directory for generated diagnostic figures and table.',
    )
    parser.add_argument('--n-bootstrap', type=int, default=12, help='Bootstrap replicates per grid point.')
    parser.add_argument('--n-cells', type=int, default=1000, help='Synthetic cells per replicate.')
    parser.add_argument('--seed-start', type=int, default=0, help='First bootstrap RNG seed.')
    parser.add_argument('--grid-size', type=int, default=7, help='Number of values per sigmoid parameter.')
    parser.add_argument('--param-grid-size', type=int, default=11, help='Number of values per one-parameter sweep.')
    parser.add_argument('--half-rate-min', type=float, default=1.70, help='Lowest sigmoid inflection point in Hz.')
    parser.add_argument('--half-rate-max', type=float, default=2.50, help='Highest sigmoid inflection point in Hz.')
    parser.add_argument('--slope-min', type=float, default=0.04, help='Lowest sigmoid slope in Hz.')
    parser.add_argument('--slope-max', type=float, default=0.16, help='Highest sigmoid slope in Hz.')
    parser.add_argument(
        '--plot-only',
        action='store_true',
        help='Redraw figures from existing sensitivity CSVs without rerunning simulations.',
    )
    args = parser.parse_args()
    output_dir = args.output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    if args.plot_only:
        redraw_from_tables(output_dir)
        print(f'[plot-only] redrew figures from existing tables in {output_dir}')
        return

    p_ref = Params(n_bootstrap=args.n_bootstrap, n_cells=args.n_cells, seed_start=args.seed_start)

    print('[default] selected parameter set')
    default_result = run_bootstrap_suite(p_ref)
    default_metrics = collect_extended_metrics(default_result)
    default_rows = [
        {'metric': spec.key, 'label': spec.label, 'value': default_metrics[spec.key], 'unit': spec.unit}
        for spec in METRIC_SPECS
    ]
    write_rows(default_rows, output_dir, 'model_parameter_sensitivity_default_metrics.csv')

    half_rates = np.linspace(args.half_rate_min, args.half_rate_max, args.grid_size)
    slopes = np.linspace(args.slope_min, args.slope_max, args.grid_size)
    if not np.any(np.isclose(half_rates, p_ref.da_half_rate, rtol=1e-9, atol=1e-9)):
        half_rates = np.sort(np.append(half_rates, p_ref.da_half_rate))
    if not np.any(np.isclose(slopes, p_ref.da_rate_slope, rtol=1e-9, atol=1e-9)):
        slopes = np.sort(np.append(slopes, p_ref.da_rate_slope))
    rows = run_sensitivity_grid(args, half_rates, slopes)
    csv_path = write_rows(
        [asdict(row) for row in rows],
        output_dir,
        'model_parameter_sensitivity_summary.csv',
    )
    plot_sensitivity_heatmaps(rows, half_rates, slopes, output_dir)
    print_grid_summary(rows, csv_path, output_dir)

    sweep_rows = run_parameter_sweeps(args, default_metrics)
    write_rows(sweep_rows, output_dir, 'model_parameter_sensitivity_parameter_sweep.csv')
    breakpoint_rows = parameter_breakpoints(sweep_rows)
    write_rows(breakpoint_rows, output_dir, 'model_parameter_sensitivity_breakpoints.csv')
    importance_rows = parameter_importance(sweep_rows, default_metrics, breakpoint_rows)
    write_rows(importance_rows, output_dir, 'model_parameter_sensitivity_importance.csv')
    plot_parameter_importance(importance_rows, output_dir)
    plot_breakpoints(breakpoint_rows, output_dir)
    plot_top_metric_curves(sweep_rows, importance_rows, default_metrics, output_dir)
    write_text_report(output_dir, default_metrics, importance_rows, breakpoint_rows, rows)
    print_extended_summary(default_metrics, importance_rows, breakpoint_rows, output_dir)

if __name__ == '__main__':
    main()
