# -*- coding: utf-8 -*-
'''
Created on Wed Feb 12 10:30:59 2025

a stud script
plot traces for Ch1 and Ch2 for whole-field signals

@author: Dinghao Luo
'''

#%% imports
import sys
from pathlib import Path

from scipy.stats import pearsonr
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

plot_first = 2000  # frames

repo_root = Path(__file__).resolve().parents[2]
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))

from console_formatting import print_session
import project_paths as pp


#%% recording lists
import rec_list
pathGRABNE = rec_list.pathHPCGRABNE
output_stem = pp.GRABNE_FIGURES_STEM / 'single_session_whole_field_ch1_ch2_correlation'


#%% main
all_r = []
all_p = []

for rec_path in pathGRABNE:
    grid_stem = Path(f'{rec_path}_grid_extract')
    whole_field_ch1 = np.squeeze(np.load(
        grid_stem / 'grid_traces_496.npy',
        allow_pickle=True))
    whole_field_ch2 = np.squeeze(np.load(
        grid_stem / 'grid_traces_496_ch2.npy',
        allow_pickle=True))

    recname = Path(rec_path).name
    print_session(recname)

    r, p = pearsonr(whole_field_ch1, whole_field_ch2)
    all_r.append(r); all_p.append(p)

    # plotting
    fig, ax = plt.subplots(figsize=(5,2.4))

    ax.plot(whole_field_ch1[:plot_first], alpha=.8, linewidth=1)
    ax.plot(whole_field_ch2[:plot_first], alpha=.8, linewidth=1)

    ax.set(title=f'R = {r}\np = {p:.10e}')

    fig.tight_layout()

    output_stem.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        output_stem / f'{recname}_first{plot_first}.png',
        dpi=300,
        bbox_inches='tight')

    plt.close(fig)


#%% main dFF
all_r_dFF = []
all_p_dFF = []

for rec_path in pathGRABNE:
    grid_stem = Path(f'{rec_path}_grid_extract')
    whole_field_ch1 = np.squeeze(np.load(
        grid_stem / 'grid_traces_dFF_496.npy',
        allow_pickle=True))
    whole_field_ch2 = np.squeeze(np.load(
        grid_stem / 'grid_traces_dFF_496_ch2.npy',
        allow_pickle=True))

    recname = Path(rec_path).name
    print_session(recname)

    r, p = pearsonr(whole_field_ch1, whole_field_ch2)
    all_r_dFF.append(r); all_p_dFF.append(p)

    s1 = whole_field_ch1 - np.mean(whole_field_ch1)
    s2 = whole_field_ch2 - np.mean(whole_field_ch2)
    corr = np.correlate(s1, s2, mode='full')
    lags = np.arange(-len(s1) + 1, len(s1))

    # plotting
    fig, ax = plt.subplots(figsize=(5,2.4))

    ax.plot(whole_field_ch1[:plot_first], alpha=.8, linewidth=1)
    ax.plot(whole_field_ch2[:plot_first], alpha=.8, linewidth=1)

    ax.set(title=f'R = {r}\np = {p:.10e}')

    fig.tight_layout()

    output_stem.mkdir(parents=True, exist_ok=True)
    fig.savefig(
        output_stem / f'{recname}_dFF_first{plot_first}.png',
        dpi=300,
        bbox_inches='tight')

    plt.close(fig)

    # plot correlagram based on dFF
    fig, ax = plt.subplots(figsize=(3,3))

    ax.plot(lags[int(len(lags)/2-600):int(len(lags)/2+600)], corr[int(len(lags)/2-600):int(len(lags)/2+600)])

    ax.set(xlabel='lag',
           ylabel='cross-correlation',
           title=f'peak @ {np.argmax(corr)-lags[-1]}')

    fig.tight_layout()

    fig.savefig(
        output_stem / f'{recname}_correlagram.png',
        dpi=300,
        bbox_inches='tight')

    plt.close(fig)
