# -*- coding: utf-8 -*-
'''
Created on Fri Jul 12 17:32:01 2024

plot run-onset- and reward-aligned sorted heatmaps of each session

@author: Dinghao Luo
'''

#%% imports
import sys
from pathlib import Path
import numpy as np
import pandas as pd

# plotting parameters
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = 'Arial'
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

repo_root = Path(__file__).resolve().parents[2]
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))
from common_functions import smooth_convolve, normalise
from console_formatting import print_session
import project_paths as pp

import rec_list
pathGRABNE = rec_list.pathHPCGRABNE


#%% significant activity only?
sig_act_only = True

if sig_act_only:
    df = pd.read_pickle(pp.GRABNE_STEM / 'significant_activity_roi.pkl')

heatmap_stem = pp.GRABNE_FIGURES_STEM / 'single_session_heatmaps'


#%% main loop
for path in pathGRABNE:  # start from the first good recording animal (A093i)
    recname = Path(path).name
    print_session(recname)

    aligned_run_path = Path(f'{path}_roi_extract') / 'suite2pROI_run_dFF_aligned.npy'
    aligned_rew_path = Path(f'{path}_roi_extract') / 'suite2pROI_rew_dFF_aligned.npy'
    aligned_run_dff = np.load(aligned_run_path, allow_pickle=True)
    aligned_rew_dff = np.load(aligned_rew_path, allow_pickle=True)
    tot_roi, tot_trial, tot_frame = aligned_run_dff.shape

    # calculate mean trace for each grid
    run_trace_means = []
    rew_trace_means = []
    for g in range(tot_roi):
        # run aligned
        run_traces = np.zeros(aligned_run_dff.shape[1:])
        for trial in range(aligned_run_dff.shape[1]):
            run_traces[trial, :] = smooth_convolve(aligned_run_dff[g, trial, :])
        run_trace_means.append(normalise(np.mean(run_traces, axis=0)))

        # rew aligned
        rew_traces = np.zeros(aligned_rew_dff.shape[1:])
        for trial in range(aligned_rew_dff.shape[1]):
            rew_traces[trial, :] = smooth_convolve(aligned_rew_dff[g, trial, :])
        rew_trace_means.append(normalise(np.mean(rew_traces, axis=0)))

    run_curr_max_pt = {}  # argmax for all mean trace
    for g in range(tot_roi):
        run_curr_max_pt[g] = np.argmax(run_trace_means[g])
    run_ord_ind = sorted(np.arange(tot_roi), key=run_curr_max_pt.get)  # ordered indices
    run_trace_means_ordered = np.zeros((tot_roi, 5*30))  # ordered heatmap
    for i, grid in enumerate(run_ord_ind):
        run_trace_means_ordered[i,:] = run_trace_means[grid]

    rew_curr_max_pt = {}  # argmax for all mean trace
    for g in range(tot_roi):
        rew_curr_max_pt[g] = np.argmax(rew_trace_means[g])
    rew_ord_ind = sorted(np.arange(tot_roi), key=rew_curr_max_pt.get)  # ordered indices
    rew_trace_means_ordered = np.zeros((tot_roi, 5*30))  # ordered heatmap
    for i, grid in enumerate(rew_ord_ind):
        rew_trace_means_ordered[i,:] = rew_trace_means[grid]

    # plot run
    fig, axs = plt.subplots(1,2, figsize=(5.5,3.3))

    axs[0].imshow(run_trace_means, aspect='auto', extent=[-1,4,0,tot_roi], cmap='viridis', interpolation='none')
    axs[1].imshow(run_trace_means_ordered, aspect='auto', extent=[-1,4,0,tot_roi], cmap='viridis', interpolation='none')

    axs[0].set(title='mean traces\n(ROIs, RO-aligned)')
    axs[1].set(title='mean traces\n(ROIs, RO-aligned, sorted)')
    for p in [0,1]:
        axs[p].set(xlabel='time (s)', ylabel='grid #')

    fig.suptitle(recname)

    fig.tight_layout()

    ro_stem = heatmap_stem / 'sorted_heatmaps_RO_aligned_rois'
    ro_stem.mkdir(parents=True, exist_ok=True)
    fig.savefig(ro_stem / f'{recname}.png', dpi=500, bbox_inches='tight')
    fig.savefig(ro_stem / f'{recname}.pdf', bbox_inches='tight')

    plt.close(fig)

    # plot rew
    fig, axs = plt.subplots(1,2, figsize=(5.5,3.3))

    axs[0].imshow(rew_trace_means, aspect='auto', extent=[-1,4,0,tot_roi], cmap='viridis', interpolation='none')
    axs[1].imshow(rew_trace_means_ordered, aspect='auto', extent=[-1,4,0,tot_roi], cmap='viridis', interpolation='none')

    axs[0].set(title='mean traces\n(grid, rew-aligned)')
    axs[1].set(title='mean traces\n(grid, rew-aligned, sorted)')
    for p in [0,1]:
        axs[p].set(xlabel='time (s)', ylabel='grid #')

    fig.suptitle(recname)

    fig.tight_layout()

    rew_stem = heatmap_stem / 'sorted_heatmaps_rew_aligned_rois'
    rew_stem.mkdir(parents=True, exist_ok=True)
    fig.savefig(rew_stem / f'{recname}.png', dpi=120, bbox_inches='tight')
    fig.savefig(rew_stem / f'{recname}.pdf', bbox_inches='tight')
