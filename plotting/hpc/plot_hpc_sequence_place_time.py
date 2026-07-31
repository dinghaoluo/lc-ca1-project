'''
Created on 6 Jul 2026
Modified on 24 June 2026

plot CA1 place-cell, time-cell, and sequence metrics, including the
equal-noise stimulation contrasts

@author: Dinghao Luo
'''


#%% imports
import argparse
from pathlib import Path
import pickle
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

repo_root = Path(__file__).resolve().parents[2]
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))

from common_functions import mpl_formatting, normalise
from console_formatting import print_files_saved
import project_paths as pp

mpl_formatting()


#%% paths
DEFAULT_PAYLOAD_PATH = (
    pp.DATA_ROOT
    / 'analysis'
    / 'hpc'
    / 'sequence_analysis'
    / 'hpc_sequence_place_time_metrics.pkl'
)
DEFAULT_OUTPUT_DIR = pp.HPC_EPHYS_FIGURES_STEM / 'sequence_analysis'

OPTO_DATASETS = ('HPCLC', 'HPCLCterm')
RAPHI_DATASETS = ('HPCRaphi',)

SUMMARY_METRICS = [
    'direct_median_abs_shift',
    'population_vector_r',
    'direct_rank_r',
]


#%% plotting


def save_figure(fig, stem):
    saved = []
    for ext in ['.png', '.pdf']:
        path = stem.with_suffix(ext)
        path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(path, dpi=300, bbox_inches='tight')
        saved.append(('figure', path))
    plt.close(fig)
    return saved

def plot_metric_strip(session_df, domain, metric, ylabel, output_dir):
    df = session_df[session_df['domain'] == domain].copy()
    df = df[np.isfinite(df[metric])]

    comparison_values = []
    for comparison in sorted(df['comparison'].unique()):
        values = np.asarray(df[df['comparison'] == comparison][metric], dtype=float)
        values = values[np.isfinite(values)]
        comparison_values.append((comparison, values))

    fig, ax = plt.subplots(figsize=(2.4, 2.0))
    positions = np.arange(1, len(comparison_values) + 1)
    ax.boxplot(
        [values for _, values in comparison_values],
        positions=positions,
        widths=0.45,
        showfliers=False,
        medianprops={'color': 'k', 'linewidth': 1},
        boxprops={'linewidth': 0.8},
        whiskerprops={'linewidth': 0.8},
        capprops={'linewidth': 0.8},
    )

    for x, (_, vals) in zip(positions, comparison_values):
        jitter = np.linspace(-0.10, 0.10, vals.size) if vals.size > 1 else [0]
        ax.scatter(
            np.asarray(jitter) + x,
            vals,
            s=8,
            color='0.25',
            alpha=0.75,
            linewidths=0,
        )

    ax.set(
        ylabel=ylabel,
        xticks=positions,
        xticklabels=[
            comparison.replace('_', ' ')
            for comparison, _ in comparison_values
        ],
        title=domain,
    )
    ax.tick_params(axis='x', rotation=30)
    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)
    fig.tight_layout()
    return save_figure(fig, output_dir / f'{domain}_{metric}')

def plot_condition_rank_pairs(session_df, domain, output_dir):
    df = session_df[session_df['domain'] == domain].copy()
    df = df[
        np.isfinite(df['condition_a_rank_r'])
        & np.isfinite(df['condition_b_rank_r'])
    ]

    saved = []
    for comparison, comp_df in df.groupby('comparison'):
        comparison_text = comparison.replace('_', ' ')
        fig, ax = plt.subplots(figsize=(1.8, 2.0))
        for _, row in comp_df.iterrows():
            ax.plot(
                [0, 1],
                [row['condition_a_rank_r'], row['condition_b_rank_r']],
                color='0.65',
                linewidth=0.8,
                zorder=1,
            )
            ax.scatter(
                [0, 1],
                [row['condition_a_rank_r'], row['condition_b_rank_r']],
                color='0.2',
                s=9,
                zorder=2,
            )

        labels = [
            comp_df.iloc[0]['condition_a'],
            comp_df.iloc[0]['condition_b'],
        ]
        ax.set(
            xticks=[0, 1],
            xticklabels=labels,
            ylabel='rank r to reference',
            ylim=(-1.05, 1.05),
            title=f'{domain}, {comparison_text}',
        )
        for spine in ['top', 'right']:
            ax.spines[spine].set_visible(False)
        fig.tight_layout()
        saved.extend(save_figure(
            fig,
            output_dir / f'{domain}_{comparison}_rank_pairs',
        ))
    return saved

def plot_heatmap_triplet(recname, comparison, domain, data, metadata, output_dir):
    ref = np.vstack([normalise(row) for row in data['ref_matrix']])
    a = np.vstack([normalise(row) for row in data['condition_a_matrix']])
    b = np.vstack([normalise(row) for row in data['condition_b_matrix']])

    comparison_text = comparison.replace('_', ' ')
    fig, axs = plt.subplots(1, 3, figsize=(5.4, 2.4), sharey=True)
    matrices = [ref, a, b]
    titles = [
        data['ref_condition'],
        data['condition_a'],
        data['condition_b'],
    ]
    for ax, matrix, title in zip(axs, matrices, titles):
        if domain == 'spatial':
            bin_size = metadata['spatial_payload_bin_size_cm']
            extent = [0, matrix.shape[1] * bin_size, matrix.shape[0], 0]
        else:
            start, stop = metadata['temporal_window_s']
            extent = [start, stop, matrix.shape[0], 0]
        ax.imshow(
            matrix,
            aspect='auto',
            cmap='turbo',
            interpolation='nearest',
            extent=extent,
        )
        ax.set(title=title)
        if domain == 'spatial':
            ax.set(xlabel='distance (cm)')
        else:
            ax.set(xlabel='time from run onset (s)')

    axs[0].set(ylabel='cells sorted by reference peak')
    fig.suptitle(f'{recname}, {comparison_text}, {domain}')
    fig.tight_layout()
    return save_figure(
        fig,
        output_dir / 'examples' / f'{recname}_{comparison}_{domain}',
    )

def plot_example_heatmaps(
        payload,
        output_dir,
        max_examples,
        allowed_recnames=None,
        valid_pairs=None
        ):
    metadata = payload['metadata']
    session_payload = payload['session_payload']
    if allowed_recnames is not None:
        allowed_recnames = set(allowed_recnames)
    if valid_pairs is not None:
        valid_pairs = set(valid_pairs)
    counts = {}
    saved = []

    for recname, rec_payload in session_payload.items():
        if allowed_recnames is not None and recname not in allowed_recnames:
            continue
        for comparison, comp_payload in rec_payload.items():
            for domain, data in comp_payload.items():
                if valid_pairs is not None and (recname, comparison, domain) not in valid_pairs:
                    continue
                key = (comparison, domain)
                if key not in counts:
                    counts[key] = 0
                if counts[key] >= max_examples:
                    continue
                saved.extend(plot_heatmap_triplet(
                    recname,
                    comparison,
                    domain,
                    data,
                    metadata,
                    output_dir,
                ))
                counts[key] += 1
    return saved

def plot_summary_set(session_df, output_dir):
    saved = []
    for domain in ['spatial', 'temporal']:
        saved.extend(plot_metric_strip(
            session_df,
            domain,
            'direct_median_abs_shift',
            'median abs. peak shift',
            output_dir,
        ))
        saved.extend(plot_metric_strip(
            session_df,
            domain,
            'population_vector_r',
            'population-vector r',
            output_dir,
        ))
        saved.extend(plot_condition_rank_pairs(
            session_df,
            domain,
            output_dir,
        ))
    return saved

def plot_stim_mixed_effects(mixed_df, output_dir):
    outcome_groups = [
        (['spatial_info', 'temporal_info'], 'information'),
        (['spatial_split_r', 'temporal_split_r'], 'split-half r'),
    ]
    contrast_styles = [
        ('stim_vs_pre', -.12, '^', 'white', 'stim - pre ctrl.'),
        ('stim_vs_sandwich', 0, 'o', 'royalblue', 'stim - sandwich ctrl.'),
        ('stim_vs_post', .12, 's', 'white', 'stim - post ctrl.'),
    ]
    saved = []

    for dataset, dataset_df in mixed_df.groupby('dataset'):
        fig, axs = plt.subplots(1, 2, figsize=(5.2, 2.2))
        for ax, (outcomes, title) in zip(axs, outcome_groups):
            for y, outcome in zip([1, 0], outcomes):
                outcome_df = dataset_df[dataset_df['outcome'] == outcome]
                for contrast, offset, marker, face, label in contrast_styles:
                    row = outcome_df[
                        outcome_df['contrast'] == contrast
                    ].iloc[0]
                    estimate = row['estimate']
                    ax.errorbar(
                        estimate,
                        y + offset,
                        xerr=np.array([[
                            estimate - row['ci_low']
                        ], [
                            row['ci_high'] - estimate
                        ]]),
                        fmt=marker,
                        color='royalblue',
                        markerfacecolor=face,
                        markeredgecolor='royalblue',
                        capsize=2,
                        label=label if y == 1 else None,
                    )

            ax.axvline(0, color='grey', lw=.8, ls='--')
            ax.set(
                yticks=[0, 1],
                yticklabels=['temporal', 'spatial'],
                ylim=(-.45, 1.45),
                xlabel='Stimulation contrast',
            )
            ax.set_title(title, fontsize=10)
            for spine in ['top', 'right']:
                ax.spines[spine].set_visible(False)

        axs[1].legend(
            frameon=False,
            fontsize=6,
            loc='center left',
            bbox_to_anchor=(1.02, .5),
        )
        fig.suptitle(dataset, fontsize=11)
        fig.tight_layout()
        saved.extend(save_figure(
            fig,
            output_dir / f'{dataset}_stim_mixed_effects',
        ))

    return saved


#%% main
def main(argv=None):
    parser = argparse.ArgumentParser(
        description='plot CA1 place/time-cell sequence metrics'
    )
    parser.add_argument(
        '--payload',
        type=Path,
        default=DEFAULT_PAYLOAD_PATH,
        help='analysis data from build_hpc_sequence_place_time_metrics.py',
    )
    parser.add_argument(
        '--output-dir',
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help='figure output folder',
    )
    parser.add_argument(
        '--max-examples',
        type=int,
        default=8,
        help='maximum sorted heatmap examples per comparison/domain',
    )
    args = parser.parse_args(argv)
    with open(args.payload, 'rb') as f:
        payload = pickle.load(f)
    session_df = payload['session_metrics'].copy()
    session_df['dataset_group'] = 'other'
    session_df.loc[session_df['dataset'].isin(OPTO_DATASETS), 'dataset_group'] = 'opto'
    session_df.loc[session_df['dataset'].isin(RAPHI_DATASETS), 'dataset_group'] = 'raphi'
    session_df = session_df[
        (session_df['comparison'] != 'stim_ctrl')
        | session_df['dataset'].isin(OPTO_DATASETS)
    ].copy()
    valid_pairs = set(zip(
        session_df['recname'],
        session_df['comparison'],
        session_df['domain'],
    ))

    grouped_frames = [
        session_df.assign(group_type='all', group_label='all'),
        session_df.assign(group_type='dataset_group', group_label=session_df['dataset_group']),
        session_df.assign(group_type='dataset', group_label=session_df['dataset']),
    ]
    summary_source = pd.concat(grouped_frames, axis=0, ignore_index=True)

    rows = []
    for keys, group in summary_source.groupby([
            'group_type',
            'group_label',
            'comparison',
            'domain',
            ]):
        group_type, group_label, comparison, domain = keys
        row = {
            'group_type': group_type,
            'group_label': group_label,
            'comparison': comparison,
            'domain': domain,
            'n_sessions': group['recname'].nunique(),
            'n_rows': group.shape[0],
            'median_n_cells': group['n_cells'].median(),
        }
        for metric in SUMMARY_METRICS:
            values = pd.to_numeric(group[metric], errors='coerce').dropna()
            row[f'{metric}_median'] = values.median() if not values.empty else np.nan
            row[f'{metric}_q25'] = values.quantile(.25) if not values.empty else np.nan
            row[f'{metric}_q75'] = values.quantile(.75) if not values.empty else np.nan
        rows.append(row)

    summary_df = pd.DataFrame(rows).sort_values([
        'group_type',
        'group_label',
        'comparison',
        'domain',
    ])
    path = args.payload.with_name('hpc_sequence_summary_by_group.csv')
    summary_df.to_csv(path, index=False)
    saved = [('summary table', path)]
    saved.extend(plot_summary_set(session_df, args.output_dir))
    saved.extend(plot_stim_mixed_effects(
        payload['mixed_effects'],
        args.output_dir,
    ))

    for group_label in ['opto', 'raphi']:
        group_df = session_df[session_df['dataset_group'] == group_label].copy()

        saved.extend(plot_summary_set(
            group_df,
            args.output_dir / f'{group_label}_only',
        ))

    saved.extend(plot_example_heatmaps(
        payload,
        args.output_dir,
        args.max_examples,
        valid_pairs=valid_pairs,
    ))
    for group_label in ['opto', 'raphi']:
        group_df = session_df[session_df['dataset_group'] == group_label]

        saved.extend(plot_example_heatmaps(
            payload,
            args.output_dir / f'{group_label}_only',
            args.max_examples,
            allowed_recnames=group_df['recname'].unique(),
            valid_pairs=valid_pairs,
        ))

    print_files_saved(saved)

if __name__ == '__main__':
    main()
