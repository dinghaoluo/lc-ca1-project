'''
Created on Aug 28 2026

test phasic DA withdrawal and sign-reversed modulation in the LC-DA-CA1 model

@author: Dinghao Luo
'''


#%% imports
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
from matplotlib.ticker import MaxNLocator
import numpy as np
from scipy.stats import sem, wilcoxon

from run_model import (
    CLASS_COLORS,
    DEFAULT_OUTPUT_DIR,
    MIDGREY,
    Params,
    RED,
    ROYALBLUE,
    clean_axis,
    make_drives,
    make_population,
    plot_mean_sem,
    save_figure_bundle,
    simulate_population_condition,
    trace_ylim,
)


#%% parameters
OUTPUT_DIR = DEFAULT_OUTPUT_DIR / 'phasic_inhibition_diagnostic'


#%% simulation
def signed_rank_p(values: np.ndarray) -> float:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if len(values) == 0 or np.allclose(values, 0.0):
        return 1.0
    return float(wilcoxon(values, zero_method='wilcox', method='approx').pvalue)


def class_mean_trace(rates: np.ndarray, mask: np.ndarray) -> np.ndarray:
    return np.nanmean(rates[mask], axis=0)


def run_diagnostic(p: Params, amplitudes: np.ndarray) -> dict:
    '''compare normal phasic DA, complete withdrawal, and conditional sign reversal.'''

    t = np.arange(-p.t_pre, p.t_post, p.dt)
    drives = make_drives(t, p)
    window = (t >= 1.0) & (t < 4.0)
    max_amplitude = float(np.max(amplitudes))

    trace_keys = (
        'normal_up', 'zero_up', 'reversed_up',
        'normal_down', 'zero_down', 'reversed_down',
    )
    traces = {key: [] for key in trace_keys}
    proportions = {
        'normal_up': [], 'zero_up': [], 'reversed_up': [],
        'normal_down': [], 'zero_down': [], 'reversed_down': [],
    }
    positive = {'up': [], 'down': []}
    withdrawal = {'up': [], 'down': []}
    negative = {
        float(amplitude): {'up': [], 'down': []}
        for amplitude in amplitudes
    }

    for seed in range(p.seed_start, p.seed_start + p.n_bootstrap):
        pop = make_population(p, np.random.default_rng(seed))
        normal = simulate_population_condition(t, p, pop, drives, da_scale=1.0)
        zero = simulate_population_condition(t, p, pop, drives, da_scale=0.0)
        selected_masks = normal['classes']

        for group, class_key in (('up', 'is_up'), ('down', 'is_down')):
            mask = selected_masks[class_key]
            traces[f'normal_{group}'].append(class_mean_trace(normal['rates'], mask))
            traces[f'zero_{group}'].append(class_mean_trace(zero['rates'], mask))
            positive[group].append(float(np.nanmean(
                (normal['rates'][mask] - zero['rates'][mask])[:, window]
            )))
            withdrawal[group].append(float(np.nanmean(
                (zero['rates'][mask] - normal['rates'][mask])[:, window]
            )))

        proportions['normal_up'].append(100.0 * np.mean(normal['classes']['is_up']))
        proportions['normal_down'].append(100.0 * np.mean(normal['classes']['is_down']))
        proportions['zero_up'].append(100.0 * np.mean(zero['classes']['is_up']))
        proportions['zero_down'].append(100.0 * np.mean(zero['classes']['is_down']))

        for amplitude in amplitudes:
            amplitude = float(amplitude)
            if np.isclose(amplitude, 0.0):
                reversed_result = zero
            else:
                reversed_result = simulate_population_condition(
                    t, p, pop, drives, da_scale=-amplitude)

            for group, class_key in (('up', 'is_up'), ('down', 'is_down')):
                mask = selected_masks[class_key]
                negative[amplitude][group].append(float(np.nanmean(
                    (reversed_result['rates'][mask] - zero['rates'][mask])[:, window]
                )))
                if np.isclose(amplitude, max_amplitude):
                    traces[f'reversed_{group}'].append(
                        class_mean_trace(reversed_result['rates'], mask)
                    )

            if np.isclose(amplitude, max_amplitude):
                proportions['reversed_up'].append(
                    100.0 * np.mean(reversed_result['classes']['is_up']))
                proportions['reversed_down'].append(
                    100.0 * np.mean(reversed_result['classes']['is_down']))

    return {
        't': t,
        'params': p,
        'amplitudes': amplitudes,
        'traces': {
            key: np.asarray(values, dtype=float)
            for key, values in traces.items()
        },
        'proportions': {
            key: np.asarray(values, dtype=float)
            for key, values in proportions.items()
        },
        'positive': {
            group: np.asarray(values, dtype=float)
            for group, values in positive.items()
        },
        'withdrawal': {
            group: np.asarray(values, dtype=float)
            for group, values in withdrawal.items()
        },
        'negative': {
            amplitude: {
                group: np.asarray(values, dtype=float)
                for group, values in groups.items()
            }
            for amplitude, groups in negative.items()
        },
    }


#%% output
def plot_three_condition_bar(
    ax: plt.Axes,
    values: tuple[np.ndarray, np.ndarray, np.ndarray],
    colors: tuple,
    title: str,
) -> None:
    x = np.arange(3)
    matrix = np.column_stack(values)
    means = np.nanmean(matrix, axis=0)
    errors = sem(matrix, axis=0, nan_policy='omit')

    ax.bar(
        x,
        means,
        width=0.8,
        color=[to_rgba(color, 0.6) for color in colors],
        edgecolor='black',
        linewidth=1.0,
        zorder=2,
    )
    for row in matrix:
        ax.plot(x, row, color='black', alpha=0.25, linewidth=0.7, zorder=3)
    for index, color in enumerate(colors):
        ax.scatter(
            np.full(len(matrix), index), matrix[:, index],
            s=9, color=color, edgecolor='none', alpha=0.48, zorder=4,
        )
    ax.errorbar(
        x, means, yerr=errors, fmt='none', ecolor='black',
        elinewidth=1.1, capsize=3, zorder=5,
    )
    ax.set_xticks(x)
    ax.set_xticklabels(['Normal', 'No phasic', 'Sign-reversed'], rotation=20, ha='right')
    ax.set_ylabel('Proportion (%)')
    ax.set_title(title)
    clean_axis(ax)
    ax.yaxis.set_major_locator(MaxNLocator(nbins=4, integer=True, min_n_ticks=3))


def effect_text(withdrawal: np.ndarray, reversed_effect: np.ndarray) -> str:
    withdrawal_mean = float(np.nanmean(withdrawal))
    reversed_mean = float(np.nanmean(reversed_effect))
    return (
        f'No phasic - normal: {withdrawal_mean:+.3f} Hz; '
        f'p={signed_rank_p(withdrawal):.3g}\n'
        f'Sign-reversed - no phasic: {reversed_mean:+.3f} Hz; '
        f'p={signed_rank_p(reversed_effect):.3g}'
    )


def plot_presentation(result: dict, output_dir: Path) -> None:
    t = result['t']
    traces = result['traces']
    proportions = result['proportions']
    amplitudes = result['amplitudes']
    max_amplitude = float(np.max(amplitudes))
    negative = result['negative'][max_amplitude]

    fig = plt.figure(figsize=(6.2, 5.0), constrained_layout=True)
    grid = fig.add_gridspec(
        2, 2, width_ratios=[2.4, 1.4], height_ratios=[1, 1],
        wspace=0.42, hspace=0.35,
    )
    ax_up = fig.add_subplot(grid[0, 0])
    ax_down = fig.add_subplot(grid[1, 0])
    ax_up_bar = fig.add_subplot(grid[0, 1])
    ax_down_bar = fig.add_subplot(grid[1, 1])

    ylim = trace_ylim(*traces.values(), pad_frac=0.10, lower_floor=0.0)
    for ax, group, title in (
        (ax_up, 'up', 'PyrUp'),
        (ax_down, 'down', 'PyrDown'),
    ):
        plot_mean_sem(ax, t, traces[f'normal_{group}'], CLASS_COLORS[f'is_{group}'], 'Normal phasic DA')
        plot_mean_sem(ax, t, traces[f'zero_{group}'], MIDGREY, 'No phasic DA')
        plot_mean_sem(
            ax, t, traces[f'reversed_{group}'], ROYALBLUE,
            f'Sign-reversed ({max_amplitude:.1f}x)',
        )
        ax.axvline(0.0, linestyle='--', color=RED, linewidth=1.0)
        ax.set_xlim([-1.0, 4.0])
        ax.set_ylim(ylim)
        ax.set_ylabel('Firing rate (Hz)')
        ax.set_title(title)
        ax.legend(frameon=False, fontsize=8)
        ax.text(
            0.98, 0.97,
            effect_text(result['withdrawal'][group], negative[group]),
            transform=ax.transAxes, ha='right', va='top', fontsize=7.4,
            bbox={'facecolor': 'white', 'edgecolor': 'none', 'alpha': 0.78, 'pad': 2.0},
        )
        clean_axis(ax)

    ax_up.set_xticklabels([])
    ax_down.set_xlabel('Time from run onset (s)')

    plot_three_condition_bar(
        ax_up_bar,
        (
            proportions['normal_up'],
            proportions['zero_up'],
            proportions['reversed_up'],
        ),
        (CLASS_COLORS['is_up'], MIDGREY, ROYALBLUE),
        'PyrUp proportion',
    )
    plot_three_condition_bar(
        ax_down_bar,
        (
            proportions['normal_down'],
            proportions['zero_down'],
            proportions['reversed_down'],
        ),
        (CLASS_COLORS['is_down'], MIDGREY, ROYALBLUE),
        'PyrDown proportion',
    )

    fig.suptitle(
        'Phasic DA withdrawal and conditional sign reversal',
        fontsize=11,
    )
    save_figure_bundle(fig, output_dir, 'fig_6_lc_da_ca1_model_phasic_inhibition_diagnostic')
    plt.close(fig)


def plot_amplitude_sweep(result: dict, output_dir: Path) -> None:
    amplitudes = np.asarray(result['amplitudes'], dtype=float)
    fig, ax = plt.subplots(figsize=(5.6, 3.8))

    for group, label in (('up', 'PyrUp'), ('down', 'PyrDown')):
        values = np.vstack([
            result['negative'][float(amplitude)][group]
            for amplitude in amplitudes
        ])
        means = np.nanmean(values, axis=1)
        errors = sem(values, axis=1, nan_policy='omit')
        color = CLASS_COLORS[f'is_{group}']
        ax.plot(amplitudes, means, color=color, linewidth=2.0, label=label)
        ax.fill_between(
            amplitudes, means - errors, means + errors,
            color=color, alpha=0.16, linewidth=0,
        )
        ax.axhline(
            -float(np.nanmean(result['positive'][group])),
            color=color, linestyle='--', linewidth=1.0,
        )

    ax.axhline(0.0, color='0.65', linewidth=0.8)
    ax.set_xlabel('Sign-reversed amplitude / selected positive pulse')
    ax.set_ylabel('Firing-rate effect relative to no phasic DA (Hz)')
    ax.set_title('Conditional phasic-inhibition amplitude sweep')
    ax.legend(frameon=False)
    clean_axis(ax)
    save_figure_bundle(fig, output_dir, 'fig_6_lc_da_ca1_model_phasic_inhibition_sweep')
    plt.close(fig)


def write_table(result: dict, output_dir: Path) -> Path:
    path = output_dir / 'phasic_inhibition_amplitude_sweep.csv'
    with path.open('w', newline='', encoding='utf-8') as stream:
        fields = [
            'amplitude_relative_to_positive_pulse',
            'pyrup_effect_hz', 'pyrup_sem_hz', 'pyrup_p',
            'pyrdown_effect_hz', 'pyrdown_sem_hz', 'pyrdown_p',
        ]
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for amplitude in result['amplitudes']:
            groups = result['negative'][float(amplitude)]
            writer.writerow({
                'amplitude_relative_to_positive_pulse': float(amplitude),
                'pyrup_effect_hz': float(np.nanmean(groups['up'])),
                'pyrup_sem_hz': float(sem(groups['up'], nan_policy='omit')),
                'pyrup_p': signed_rank_p(groups['up']),
                'pyrdown_effect_hz': float(np.nanmean(groups['down'])),
                'pyrdown_sem_hz': float(sem(groups['down'], nan_policy='omit')),
                'pyrdown_p': signed_rank_p(groups['down']),
            })
    return path


def amplitude_matching_positive(result: dict, group: str) -> tuple[float, float, bool]:
    amplitudes = np.asarray(result['amplitudes'], dtype=float)
    suppression = np.asarray([
        -np.nanmean(result['negative'][float(amplitude)][group])
        for amplitude in amplitudes
    ])
    target = float(np.nanmean(result['positive'][group]))
    if target > np.max(suppression):
        return float(amplitudes[-1]), target, False
    matched_amplitude = float(np.interp(target, suppression, amplitudes))
    return matched_amplitude, target, True


def print_summary(result: dict, table_path: Path, output_dir: Path) -> None:
    max_amplitude = float(np.max(result['amplitudes']))
    p = result['params']
    print('Phasic DA withdrawal and sign-reversal diagnostic')
    print('  negative amplitudes are conditional on tonic modulatory reserve')
    print(f'  bootstrapping = {p.n_bootstrap}')
    print(f'  cells = {p.n_cells}')

    for group, label in (('up', 'PyrUp'), ('down', 'PyrDown')):
        positive = result['positive'][group]
        withdrawal = result['withdrawal'][group]
        reversed_effect = result['negative'][max_amplitude][group]
        matched_amplitude, target, target_reached = amplitude_matching_positive(result, group)
        print(f'  {label}')
        print(f'    normal - no phasic = {np.nanmean(positive):+.4f} Hz')
        print(f'    no phasic - normal = {np.nanmean(withdrawal):+.4f} Hz')
        print(
            f'    sign-reversed {max_amplitude:.1f}x - no phasic = '
            f'{np.nanmean(reversed_effect):+.4f} Hz'
        )
        if target_reached:
            print(
                f'    amplitude matching the {target:.4f} Hz positive effect = '
                f'{matched_amplitude:.3f}x'
            )
        else:
            print(
                f'    {target:.4f} Hz positive effect not matched by '
                f'{matched_amplitude:.3f}x sign reversal'
            )

    print(f'  table: {table_path}')
    print(f'  figures: {output_dir}')


#%% main
def main() -> None:
    parser = argparse.ArgumentParser(
        description='Test phasic DA withdrawal and conditional sign-reversed modulation.'
    )
    parser.add_argument('--output-dir', type=Path, default=OUTPUT_DIR)
    parser.add_argument('--n-bootstrap', type=int, default=50)
    parser.add_argument('--n-cells', type=int, default=1000)
    parser.add_argument('--seed-start', type=int, default=0)
    parser.add_argument('--amplitude-max', type=float, default=1.0)
    parser.add_argument('--grid-size', type=int, default=11)
    args = parser.parse_args()

    output_dir = args.output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    p = Params(
        n_bootstrap=args.n_bootstrap,
        n_cells=args.n_cells,
        seed_start=args.seed_start,
    )
    amplitudes = np.linspace(0.0, args.amplitude_max, args.grid_size)
    result = run_diagnostic(p, amplitudes)
    plot_presentation(result, output_dir)
    plot_amplitude_sweep(result, output_dir)
    table_path = write_table(result, output_dir)
    print_summary(result, table_path, output_dir)


if __name__ == '__main__':
    main()
