# -*- coding: utf-8 -*-
'''
Created on 20 May 2026

plot LC lick-sensitivity and peak-detection figures from the saved profile table

@author: Dinghao Luo
'''

#%% imports
from pathlib import Path
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

repo_root = Path(__file__).resolve().parents[2]
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))

from common_functions import mpl_formatting
import peak_detection_functions as pdf
import project_paths as pp
mpl_formatting()


#%% paths
payload_path = pp.LC_EPHYS_STEM / 'LC_all_profiles_plot_payload.npy'


#%% lick-sensitivity plots
def plot_lick_sensitivity(payload, samp_freq=1250):
    cluname = payload['cluname']
    suffix = payload['suffix']
    identity = payload['identity']
    lick_sensitive_signif = payload['lick_sensitive_signif']
    aligned_rasters = payload['aligned_rasters']
    aligned_prof_mean = payload['aligned_prof_mean']
    aligned_prof_sem = payload['aligned_prof_sem']
    shuf_ratios = payload['shuf_ratios']
    true_ratio = payload['true_ratio']

    xaxis = np.arange(6 * samp_freq) / samp_freq - 3

    fig = plt.figure(figsize=(3.5, 1.9))
    gs = fig.add_gridspec(1, 2, width_ratios=[2.5, 1])
    axs = [fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[0, 1])]

    for i, trial in enumerate(aligned_rasters):
        axs[0].scatter(trial, [i + 1] * len(trial), color='grey', alpha=.25, s=.6)

    axs[0].set(
        xlabel='time to 1st lick (s)',
        xlim=(-3, 3),
        xticks=[-3, 0, 3],
        ylabel='trial #',
        title=cluname + suffix
    )
    axs[0].title.set_fontsize(10)

    ax_twin = axs[0].twinx()
    ax_twin.plot(xaxis, aligned_prof_mean, color='k')
    ax_twin.fill_between(
        xaxis,
        aligned_prof_mean + aligned_prof_sem,
        aligned_prof_mean - aligned_prof_sem,
        color='k', alpha=.25, edgecolor='none'
    )
    ax_twin.set(ylabel='spike rate (Hz)')

    axs[1].plot([-1, 1], [shuf_ratios[2], shuf_ratios[2]], color='grey')
    axs[1].plot([-1, 1], [shuf_ratios[3], shuf_ratios[3]], color='grey')
    axs[1].plot([-1, 1], [shuf_ratios[-3], shuf_ratios[-3]], color='grey')
    axs[1].plot([-1, 1], [true_ratio, true_ratio], color='red')
    axs[1].set(
        xlim=(-2, 2),
        xticks=[],
        xticklabels=[],
        ylabel='post-pre ratio',
        title=lick_sensitive_signif
    )
    for spine in ['top', 'right', 'bottom']:
        axs[1].spines[spine].set_visible(False)
    axs[0].spines['top'].set_visible(False)
    ax_twin.spines['top'].set_visible(False)

    fig.tight_layout()

    idstring = f'{cluname} {identity} {suffix}'
    output_stem = (
        pp.LC_EPHYS_FIGURES_STEM
        / 'lick_sensitivity'
        / 'rasters_first_lick_aligned'
        / idstring
        )
    output_stem.parent.mkdir(parents=True, exist_ok=True)
    for ext in ['.png', '.pdf']:
        fig.savefig(
            str(output_stem) + ext,
            dpi=300,
            bbox_inches='tight'
            )
    plt.close(fig)


#%% load data
payload = np.load(payload_path, allow_pickle=True).item()

for peak_payload in payload['peak_detection'].values():
    cluname = peak_payload['cluname']
    identity = peak_payload['identity']
    peak = peak_payload['peak']
    output_stem = (
        pp.LC_EPHYS_FIGURES_STEM
        / 'peak_detection'
        / f'{cluname} {identity} {peak}'
        )
    output_stem.parent.mkdir(parents=True, exist_ok=True)
    pdf.plot_peak_v_shuf(
        cluname,
        peak_payload['mean_prof'],
        peak_payload['shuf_prof'],
        peak,
        peak_width=2,
        savepath=output_stem
        )

for lick_payload in payload['lick_sensitivity'].values():
    plot_lick_sensitivity(lick_payload)
