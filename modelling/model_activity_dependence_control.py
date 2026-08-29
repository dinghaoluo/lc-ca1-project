'''
Created on Aug 28 2026

test whether fixed DA gain can reproduce selective PyrUp enhancement

@author: Dinghao Luo
'''


#%% imports
from __future__ import annotations

import argparse
import csv
from dataclasses import replace
from pathlib import Path

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import sem, wilcoxon

from run_model import (
    CLASS_COLORS,
    DEFAULT_OUTPUT_DIR,
    Params,
    experiment_shared_ylims,
    make_drives,
    make_population,
    plot_experiment_1,
    plot_experiment_2,
    plot_experiment_3,
    run_bootstrap_suite,
    sigmoid,
    simulate_population_condition,
    summarize_result,
)


#%% parameters
CONTROL_OUTPUT_DIR = DEFAULT_OUTPUT_DIR / 'activity_dependence_control'

EXPERIMENT_LABELS = {
    'lc_activation': 'LC activation',
    'da_targeted': 'DA-targeted cells',
    'da_blockade': 'DA blockade',
}


#%% paired simulation
def population_effects(
    base: dict,
    lc: dict,
    block: dict,
    pop: dict[str, np.ndarray],
    analysis_masks: dict[str, dict[str, np.ndarray]],
    t: np.ndarray,
    p: Params,
) -> dict[str, dict[str, float]]:
    '''read out every control using the selected model's original cell sets.'''

    pre = (t >= p.pre_window[0]) & (t < p.pre_window[1])
    window = (t >= 1.0) & (t < 4.0)
    base_change = np.nanmean(base['rates'][:, window], axis=1) - np.nanmean(
        base['rates'][:, pre], axis=1)

    effects = {
        experiment: {'up': np.nan, 'down': np.nan}
        for experiment in EXPERIMENT_LABELS
    }
    for group in ('up', 'down'):
        lc_mask = analysis_masks['lc_activation'][group]
        targeted_mask = analysis_masks['da_targeted'][group]
        block_mask = analysis_masks['da_blockade'][group]
        da_mask = targeted_mask & pop['da_targ']
        non_da_mask = targeted_mask & ~pop['da_targ']

        effects['lc_activation'][group] = float(np.nanmean(
            (lc['rates'][lc_mask] - base['rates'][lc_mask])[:, window]
        ))
        effects['da_targeted'][group] = float(
            np.nanmean(base_change[da_mask]) - np.nanmean(base_change[non_da_mask])
        )
        effects['da_blockade'][group] = float(np.nanmean(
            (base['rates'][block_mask] - block['rates'][block_mask])[:, window]
        ))

    return effects


def append_effects(
    target: dict[str, dict[str, list[float]]],
    source: dict[str, dict[str, float]],
) -> None:
    for experiment in EXPERIMENT_LABELS:
        for group in ('up', 'down'):
            target[experiment][group].append(source[experiment][group])


def as_arrays(
    effects: dict[str, dict[str, list[float]]],
) -> dict[str, dict[str, np.ndarray]]:
    return {
        experiment: {
            group: np.asarray(values, dtype=float)
            for group, values in groups.items()
        }
        for experiment, groups in effects.items()
    }


def run_activity_dependent_reference(
    p: Params,
) -> tuple[
    np.ndarray,
    dict[str, np.ndarray],
    list[tuple[dict[str, np.ndarray], dict[str, dict[str, np.ndarray]]]],
    dict[str, dict[str, np.ndarray]],
    float,
]:
    '''run the selected model and retain its baseline cell classes for every control.'''

    t = np.arange(-p.t_pre, p.t_post, p.dt)
    drives = make_drives(t, p)
    records = []
    effects = {
        experiment: {'up': [], 'down': []}
        for experiment in EXPERIMENT_LABELS
    }
    dose_numerator = 0.0
    dose_denominator = 0.0

    for seed in range(p.seed_start, p.seed_start + p.n_bootstrap):
        pop = make_population(p, np.random.default_rng(seed))
        base = simulate_population_condition(t, p, pop, drives, da_scale=1.0)
        lc = simulate_population_condition(t, p, pop, drives, da_scale=p.lc_activation_fold)
        block = simulate_population_condition(t, p, pop, drives, da_scale=p.da_block_scale)
        analysis_masks = {
            'lc_activation': {
                'up': base['classes']['is_up'] & lc['classes']['is_up'],
                'down': base['classes']['is_down'] & lc['classes']['is_down'],
            },
            'da_targeted': {
                'up': base['classes']['is_up'],
                'down': base['classes']['is_down'],
            },
            'da_blockade': {
                'up': base['classes']['is_up'] & block['classes']['is_up'],
                'down': base['classes']['is_down'] & block['classes']['is_down'],
            },
        }
        append_effects(
            effects,
            population_effects(base, lc, block, pop, analysis_masks, t, p),
        )
        records.append((pop, analysis_masks))

        previous_rates = np.column_stack((base['rates'][:, 0], base['rates'][:, :-1]))
        da_gain = sigmoid(previous_rates, p.da_half_rate, p.da_rate_slope)
        da_weights = pop['da_targ_strength'][:, None] * drives['D'][None, :]
        dose_numerator += float(np.sum(da_weights * da_gain))
        dose_denominator += float(np.sum(da_weights))

    return t, drives, records, as_arrays(effects), dose_numerator / dose_denominator


def run_fixed_control(
    p: Params,
    t: np.ndarray,
    drives: dict[str, np.ndarray],
    records: list[tuple[dict[str, np.ndarray], dict[str, dict[str, np.ndarray]]]],
    gain: float,
) -> dict[str, dict[str, np.ndarray]]:
    fixed_p = replace(p, da_gain_mode='fixed', da_fixed_gain=float(gain))
    effects = {
        experiment: {'up': [], 'down': []}
        for experiment in EXPERIMENT_LABELS
    }

    for pop, analysis_masks in records:
        base = simulate_population_condition(t, fixed_p, pop, drives, da_scale=1.0)
        lc = simulate_population_condition(
            t, fixed_p, pop, drives, da_scale=fixed_p.lc_activation_fold)
        block = simulate_population_condition(
            t, fixed_p, pop, drives, da_scale=fixed_p.da_block_scale)
        append_effects(
            effects,
            population_effects(base, lc, block, pop, analysis_masks, t, fixed_p),
        )

    return as_arrays(effects)


def signed_rank_p(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if len(values) == 0 or np.allclose(values, 0.0):
        return 1.0
    return float(wilcoxon(values, zero_method='wilcox', method='approx').pvalue)


def one_row(
    mode: str,
    gain: float,
    experiment: str,
    effects: dict[str, np.ndarray],
    active_effects: dict[str, np.ndarray],
) -> dict[str, float | str]:
    up = np.asarray(effects['up'], dtype=float)
    down = np.asarray(effects['down'], dtype=float)
    selectivity = up - down
    active_selectivity = active_effects['up'] - active_effects['down']

    return {
        'mode': mode,
        'fixed_gain': gain,
        'experiment': experiment,
        'pyrup_effect_hz': float(np.nanmean(up)),
        'pyrup_sem_hz': float(sem(up, nan_policy='omit')),
        'pyrdown_effect_hz': float(np.nanmean(down)),
        'pyrdown_sem_hz': float(sem(down, nan_policy='omit')),
        'selectivity_hz': float(np.nanmean(selectivity)),
        'selectivity_sem_hz': float(sem(selectivity, nan_policy='omit')),
        'pyrup_p': signed_rank_p(up),
        'pyrdown_p': signed_rank_p(down),
        'selectivity_p': signed_rank_p(selectivity),
        'selectivity_vs_activity_dependent_p': (
            np.nan if mode == 'activity_dependent'
            else signed_rank_p(selectivity - active_selectivity)
        ),
    }


#%% output
def write_table(rows: list[dict[str, float | str]], output_dir: Path) -> Path:
    output_path = output_dir / 'activity_dependence_fixed_gain_sweep.csv'
    with output_path.open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return output_path


def plot_control(
    rows: list[dict[str, float | str]],
    matched_gain: float,
    output_dir: Path,
) -> tuple[Path, Path]:
    fig, axs = plt.subplots(1, 3, figsize=(10.8, 3.35), sharey=True)

    for ax, (experiment, title) in zip(axs, EXPERIMENT_LABELS.items()):
        active = next(
            row for row in rows
            if row['mode'] == 'activity_dependent' and row['experiment'] == experiment
        )
        fixed = [
            row for row in rows
            if row['mode'] == 'fixed' and row['experiment'] == experiment
        ]
        gains = np.asarray([row['fixed_gain'] for row in fixed], dtype=float)

        for group, label in (('pyrup', 'PyrUp'), ('pyrdown', 'PyrDown')):
            means = np.asarray([row[f'{group}_effect_hz'] for row in fixed], dtype=float)
            errors = np.asarray([row[f'{group}_sem_hz'] for row in fixed], dtype=float)
            color = CLASS_COLORS['is_up' if group == 'pyrup' else 'is_down']
            ax.plot(gains, means, color=color, linewidth=1.8, label=f'fixed, {label}')
            ax.fill_between(gains, means - errors, means + errors, color=color, alpha=0.16, linewidth=0)
            ax.axhline(
                active[f'{group}_effect_hz'],
                color=color,
                linestyle='--',
                linewidth=1.1,
                label=f'activity-dependent, {label}',
            )

        ax.axvline(matched_gain, color='0.35', linestyle=':', linewidth=1.2)
        ax.axhline(0.0, color='0.75', linewidth=0.8)
        ax.set_title(title)
        ax.set_xlabel('Fixed DA gain')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    axs[0].set_ylabel('DA-dependent firing-rate effect (Hz)')
    handles, labels = axs[-1].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', bbox_to_anchor=(0.5, 1.03), ncol=2, frameon=False)
    fig.suptitle(
        f'Activity-independent DA control; dotted line is input-dose matched (gain={matched_gain:.3f})',
        y=1.14,
    )
    fig.tight_layout()

    png_path = output_dir / 'activity_dependence_fixed_gain_sweep.png'
    pdf_path = output_dir / 'activity_dependence_fixed_gain_sweep.pdf'
    fig.savefig(png_path, dpi=300, bbox_inches='tight')
    fig.savefig(pdf_path, bbox_inches='tight')
    plt.close(fig)
    return png_path, pdf_path


def plot_fixed_gain_presentation(
    p: Params,
    matched_gain: float,
    output_dir: Path,
) -> Path:
    '''render the dose-matched control with the regular Figure 6 panels.'''

    presentation_dir = output_dir / 'fixed_gain_dose_matched'
    fixed_p = replace(p, da_gain_mode='fixed', da_fixed_gain=matched_gain)
    result = run_bootstrap_suite(fixed_p)
    summary = summarize_result(result)
    ylims = experiment_shared_ylims(result)
    plot_experiment_1(result, summary, presentation_dir, ylims)
    plot_experiment_2(result, summary, presentation_dir, ylims)
    plot_experiment_3(result, summary, presentation_dir, ylims)
    return presentation_dir


def print_matched_control(
    rows: list[dict[str, float | str]],
    matched_gain: float,
) -> None:
    print('Activity-independent DA control')
    print(f'  input-dose-matched fixed gain = {matched_gain:.4f}')

    for experiment, label in EXPERIMENT_LABELS.items():
        active = next(
            row for row in rows
            if row['mode'] == 'activity_dependent' and row['experiment'] == experiment
        )
        matched = min(
            (
                row for row in rows
                if row['mode'] == 'fixed' and row['experiment'] == experiment
            ),
            key=lambda row: abs(float(row['fixed_gain']) - matched_gain),
        )
        active_up = float(active['pyrup_effect_hz'])
        active_down = float(active['pyrdown_effect_hz'])
        matched_up = float(matched['pyrup_effect_hz'])
        matched_down = float(matched['pyrdown_effect_hz'])
        comparison_p = float(matched['selectivity_vs_activity_dependent_p'])
        fixed_rows = sorted(
            (
                row for row in rows
                if row['mode'] == 'fixed' and row['experiment'] == experiment
            ),
            key=lambda row: float(row['fixed_gain']),
        )
        fixed_gains = np.asarray([row['fixed_gain'] for row in fixed_rows], dtype=float)
        fixed_up = np.asarray([row['pyrup_effect_hz'] for row in fixed_rows], dtype=float)
        fixed_down = np.asarray([row['pyrdown_effect_hz'] for row in fixed_rows], dtype=float)
        up_matched_gain = float(np.interp(active_up, fixed_up, fixed_gains))
        up_matched_down = float(np.interp(up_matched_gain, fixed_gains, fixed_down))
        print(f'  {label}')
        print(
            '    activity-dependent: '
            f'PyrUp={active_up:+.4f} Hz, '
            f'PyrDown={active_down:+.4f} Hz'
        )
        print(
            '    fixed, dose-matched: '
            f'PyrUp={matched_up:+.4f} Hz, '
            f'PyrDown={matched_down:+.4f} Hz, '
            f'selectivity comparison p={comparison_p:.3g}'
        )
        print(
            f'    fixed, PyrUp-matched gain={up_matched_gain:.3f}: '
            f'PyrDown={up_matched_down:+.4f} Hz '
            f'(activity-dependent PyrDown={active_down:+.4f} Hz)'
        )


#%% main
def main() -> None:
    parser = argparse.ArgumentParser(
        description='Compare activity-dependent DA gain with fixed-gain controls.'
    )
    parser.add_argument('--output-dir', type=Path, default=CONTROL_OUTPUT_DIR)
    parser.add_argument('--n-bootstrap', type=int, default=50)
    parser.add_argument('--n-cells', type=int, default=1000)
    parser.add_argument('--seed-start', type=int, default=0)
    parser.add_argument('--gain-min', type=float, default=0.0)
    parser.add_argument('--gain-max', type=float, default=1.0)
    parser.add_argument('--grid-size', type=int, default=11)
    args = parser.parse_args()

    output_dir = args.output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    p = Params(
        n_bootstrap=args.n_bootstrap,
        n_cells=args.n_cells,
        seed_start=args.seed_start,
    )
    t, drives, records, active_effects, matched_gain = run_activity_dependent_reference(p)

    gains = np.linspace(args.gain_min, args.gain_max, args.grid_size)
    if not np.any(np.isclose(gains, matched_gain, rtol=1e-9, atol=1e-9)):
        gains = np.sort(np.append(gains, matched_gain))

    rows = [
        one_row('activity_dependent', np.nan, experiment, effects, effects)
        for experiment, effects in active_effects.items()
    ]

    for gain in gains:
        print(f'running fixed gain {gain:.4f}')
        fixed_effects = run_fixed_control(p, t, drives, records, float(gain))
        rows.extend(
            one_row('fixed', float(gain), experiment, effects, active_effects[experiment])
            for experiment, effects in fixed_effects.items()
        )

    table_path = write_table(rows, output_dir)
    png_path, pdf_path = plot_control(rows, matched_gain, output_dir)
    presentation_dir = plot_fixed_gain_presentation(p, matched_gain, output_dir)
    print_matched_control(rows, matched_gain)
    print(f'  table: {table_path}')
    print(f'  figure: {png_path}')
    print(f'  vector figure: {pdf_path}')
    print(f'  regular-format control panels: {presentation_dir}')


if __name__ == '__main__':
    main()
