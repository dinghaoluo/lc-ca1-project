'''
Created on 30 April 2026
Modified on 24 June 2026

plot activation and inhibition latency distributions from the saved HPCLC
stimulation analysis, including the probabilistic latency estimate

@author: Dinghao Luo
'''

#%% imports
from pathlib import Path
import argparse
import pickle
import sys

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import ranksums, sem, ttest_ind

repo_root = Path(__file__).resolve().parents[3]
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))

from common_functions import mpl_formatting
from console_formatting import print_files_saved
import project_paths as pp
mpl_formatting()


#%% parameters
ACT_COLOR        = 'indianred'
INH_HIST_COLOR   = 'lightsteelblue'
INH_STRONG_COLOR = 'steelblue'

default_input_path = (
    repo_root / 'data' / 'analysis' / 'hpc' / 'stim_ctrl'
    / 'lc_stim_ca1_effect_latency.pkl')
default_output_stem = pp.LC_STIM_CA1_EFFECT_FIGURES_STEM


#%% plotting
def _plot_latency_histograms(summary, output_stem):
    activation = np.asarray(summary['activation_times'], dtype=float)
    activation = activation[np.isfinite(activation)]
    inhibition = np.asarray(summary['inhibition_times'], dtype=float)
    inhibition = inhibition[np.isfinite(inhibition)]

    max_latency = max(max(activation), max(inhibition))
    hist_specs = [
        (0.1, (3, 3.4), (0, max_latency), 'HPCLC_act_inh_full_hist'),
        (0.05, (2.6, 3), (0, 1), 'HPCLC_act_inh_0_1_hist'),
    ]

    for bin_width, figsize, xlim, filename in hist_specs:
        bin_edges = np.arange(0.0, max_latency, bin_width)
        fig, axs = plt.subplots(2, 1, figsize=figsize, sharex=True)

        axs[0].hist(
            activation,
            bins=bin_edges,
            density=True,
            color=ACT_COLOR,
            edgecolor='k',
            label='activation',
        )
        axs[1].hist(
            inhibition,
            bins=bin_edges,
            density=True,
            color=INH_HIST_COLOR,
            edgecolor='k',
            label='inhibition',
        )

        axs[0].axvline(summary['activation_median'], color='k', lw=1)
        axs[1].axvline(summary['inhibition_median'], color='k', lw=1)

        act_q25, act_q75 = summary['activation_iqr']
        inh_q25, inh_q75 = summary['inhibition_iqr']
        axs[0].set(
            title=(
                f'Activation\n'
                f'med={summary["activation_median"]:.3g}s, '
                f'MAD/sqrt(n)={summary["activation_mad_sem"]:.3g}s\n'
                f'IQR=[{act_q25:.3g}, {act_q75:.3g}]s'
            )
        )
        axs[1].set(
            title=(
                f'Inhibition\n'
                f'med={summary["inhibition_median"]:.3g}s, '
                f'MAD/sqrt(n)={summary["inhibition_mad_sem"]:.3g}s\n'
                f'IQR=[{inh_q25:.3g}, {inh_q75:.3g}]s'
            ),
            xlabel='Time from run/stim. onset (s)',
        )

        for ax in axs:
            ax.set(xlim=xlim, yticks=[0, 1], ylim=(0, 1.2),
                   ylabel='Density')

        fig.tight_layout()
        for ext in ['.png', '.pdf']:
            fig.savefig(
                output_stem / f'{filename}{ext}',
                dpi=300,
                bbox_inches='tight',
            )
        plt.close(fig)

def _plot_latency_ecdf(summary, output_stem):
    activation = np.asarray(summary['activation_times'], dtype=float)
    activation = np.sort(activation[np.isfinite(activation)])
    inhibition = np.asarray(summary['inhibition_times'], dtype=float)
    inhibition = np.sort(inhibition[np.isfinite(inhibition)])

    x0 = np.concatenate([[activation[0]], activation])
    x1 = np.concatenate([[inhibition[0]], inhibition])
    y0 = np.concatenate([[0], np.arange(1, len(activation) + 1) / len(activation)])
    y1 = np.concatenate([[0], np.arange(1, len(inhibition) + 1) / len(inhibition)])

    fig, ax = plt.subplots(figsize=(2.4, 3.0))
    ax.plot(x0, y0, label='activation', color=ACT_COLOR)
    ax.plot(x1, y1, label='inhibition', color=INH_STRONG_COLOR)

    ax.set_xlabel('Time from run/stim. onset (s)')
    ax.set_ylabel('Cumulative probability')
    ax.set_title('Latency to sustained effect')
    ax.legend(frameon=False, fontsize=6)

    ax.text(
        0.97,
        0.03,
        f'KS D={summary["ks_stat"]:.3f}\np={summary["ks_p"]:.3g}',
        transform=ax.transAxes,
        ha='right',
        va='bottom',
        fontsize=7,
    )

    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)

    fig.tight_layout()
    for ext in ['.png', '.pdf']:
        fig.savefig(
            output_stem / f'HPCLC_latency_ecdf{ext}',
            dpi=300,
            bbox_inches='tight',
        )
    plt.close(fig)

def _plot_latency_violin(summary, output_stem):
    activation = np.asarray(summary['activation_times'], dtype=float)
    activation = activation[np.isfinite(activation)]
    inhibition = np.asarray(summary['inhibition_times'], dtype=float)
    inhibition = inhibition[np.isfinite(inhibition)]

    fig, ax = plt.subplots(figsize=(1.8, 2.2))

    vp = ax.violinplot(
        [activation, inhibition],
        positions=[1.1, 1.9],
        showmeans=False,
        showmedians=True,
        showextrema=False,
    )

    vp['bodies'][0].set_color(ACT_COLOR)
    vp['bodies'][1].set_color(INH_STRONG_COLOR)
    vp['cmedians'].set_color('k')
    vp['cmedians'].set_linewidth(2)

    ax.scatter(1.25, np.median(activation), s=30, c=ACT_COLOR,
               ec='none', lw=.5, zorder=2)
    ax.scatter(1.75, np.median(inhibition), s=30, c=INH_STRONG_COLOR,
               ec='none', lw=.5, zorder=2)
    ax.plot([1.25, 1.75], [np.median(activation), np.median(inhibition)],
            color='k', linewidth=2, zorder=1)

    for i in [0, 1]:
        vp['bodies'][i].set_edgecolor('none')
        vp['bodies'][i].set_alpha(.75)
        body = vp['bodies'][i]
        mid_x = np.mean(body.get_paths()[0].vertices[:, 0])
        if i == 0:
            body.get_paths()[0].vertices[:, 0] = np.clip(
                body.get_paths()[0].vertices[:, 0], -np.inf, mid_x)
        else:
            body.get_paths()[0].vertices[:, 0] = np.clip(
                body.get_paths()[0].vertices[:, 0], mid_x, np.inf)

    ax.scatter([1.25] * len(activation), activation, s=10, c=ACT_COLOR,
               ec='none', lw=.5, alpha=.25)
    ax.scatter([1.75] * len(inhibition), inhibition, s=10,
               c=INH_STRONG_COLOR, ec='none', lw=.5, alpha=.25)

    mean0 = np.nanmean(activation)
    mean1 = np.nanmean(inhibition)
    sem0 = sem(activation, nan_policy='omit')
    sem1 = sem(inhibition, nan_policy='omit')
    med0 = np.nanmedian(activation)
    med1 = np.nanmedian(inhibition)
    q25_0, q75_0 = np.percentile(activation, [25, 75])
    q25_1, q75_1 = np.percentile(inhibition, [25, 75])

    y_max = max(np.nanmax(activation), np.nanmax(inhibition))
    y_min = min(np.nanmin(activation), np.nanmin(inhibition))
    y_span = y_max - y_min if y_max > y_min else 1
    y_offset = 0.05 * y_span

    ax.text(
        1.1,
        y_max + y_offset,
        f'Med = {med0:.2f}\n'
        f'IQR = [{q25_0:.2f}, {q75_0:.2f}]\n'
        f'{mean0:.2f} +/- {sem0:.2f}',
        ha='center',
        va='bottom',
        fontsize=7,
        color=ACT_COLOR,
    )
    ax.text(
        1.9,
        y_max + y_offset,
        f'Med = {med1:.2f}\n'
        f'IQR = [{q25_1:.2f}, {q75_1:.2f}]\n'
        f'{mean1:.2f} +/- {sem1:.2f}',
        ha='center',
        va='bottom',
        fontsize=7,
        color=INH_STRONG_COLOR,
    )

    y_top = max(np.max(activation), np.max(inhibition))
    y_bottom = min(np.min(activation), np.min(inhibition))
    y_range_tot = y_top - y_bottom if y_top > y_bottom else 1

    _, rank_p = ranksums(activation, inhibition)
    _, ttest_p = ttest_ind(activation, inhibition)
    stat_y = y_top + 0.05 * y_range_tot
    ax.plot([1.1, 1.9], [stat_y, stat_y], c='k', lw=.5)
    ax.text(
        1.5,
        stat_y,
        f'ranksums_p={rank_p:.2e}\nttest_p={ttest_p:.2e}',
        ha='center',
        va='bottom',
        color='k',
        fontsize=8,
    )

    ax.set(xticks=[1, 2],
           xticklabels=['Act.', 'Inh.'],
           ylabel='Time from stim. onset (s)',
           title=' ')
    ax.set(xlim=(.5, 2.5))

    for spine in ['top', 'right', 'bottom']:
        ax.spines[spine].set_visible(False)

    fig.tight_layout()
    plt.grid(False)
    for ext in ['.png', '.pdf']:
        fig.savefig(
            output_stem / f'HPCLC_latency_violin{ext}',
            dpi=300,
            bbox_inches='tight',
        )
    plt.close(fig)

def _plot_probabilistic_latency(latency, output_stem):
    onset = latency['onset_median_s'].to_numpy()
    activation_weight = latency['p_activation'].to_numpy()
    inhibition_weight = latency['p_inhibition'].to_numpy()
    no_response_weight = latency['p_no_response'].to_numpy()
    bins = np.arange(0, 3.6, .1)

    # hard labels throw away the main gain of the model, so each cell contributes
    # to both curves according to its posterior sign probability
    fig, ax = plt.subplots(figsize=(3.2, 2.8))
    ax.hist(
        onset, bins=bins,
        weights=activation_weight / len(latency),
        histtype='step', linewidth=1.3, color=ACT_COLOR,
        label=f'activation (mean p={np.mean(activation_weight):.2f})',
    )
    ax.hist(
        onset, bins=bins,
        weights=inhibition_weight / len(latency),
        histtype='step', linewidth=1.3, color=INH_STRONG_COLOR,
        label=f'inhibition (mean p={np.mean(inhibition_weight):.2f})',
    )
    ax.set(
        xlim=(0, 3.5), xlabel='Time from stim. onset (s)',
        ylabel='Posterior cell fraction / 100 ms',
    )
    ax.set_title('Probabilistic response latency', fontsize=10)
    ax.text(
        .97, .74,
        f'mean p(no response)={np.mean(no_response_weight):.2f}',
        transform=ax.transAxes, ha='right', va='top', fontsize=7,
    )
    ax.legend(frameon=False, fontsize=7)
    fig.tight_layout()
    for ext in ['.png', '.pdf']:
        fig.savefig(
            output_stem / f'HPCLC_probabilistic_latency{ext}',
            dpi=300, bbox_inches='tight',
        )
    plt.close(fig)

def main():
    parser = argparse.ArgumentParser(
        description='plot LC-stimulation CA1 effect latencies.'
    )
    parser.add_argument(
        '--input', type=Path, default=default_input_path,
        help='analysis data produced by estimate_lc_stim_ca1_effect_latency.py',
    )
    parser.add_argument(
        '--output-dir', type=Path, default=default_output_stem,
        help='directory for exported figure files',
    )
    args = parser.parse_args()
    with open(args.input, 'rb') as f:
        payload = pickle.load(f)

    output_stem = args.output_dir
    output_stem.mkdir(parents=True, exist_ok=True)
    summary = payload['latency_summary']
    _plot_latency_histograms(summary, output_stem)
    _plot_latency_ecdf(summary, output_stem)
    _plot_latency_violin(summary, output_stem)
    _plot_probabilistic_latency(payload['bayesian_latency'], output_stem)
    print_files_saved([
        ('figures', output_stem),
    ], gap=1)

if __name__ == '__main__':
    main()
