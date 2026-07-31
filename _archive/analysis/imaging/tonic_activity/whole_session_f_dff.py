# -*- coding: utf-8 -*-
"""
Created on Wed Jul 24 15:18:22 2024

plots F and dF/F of single ROIs throughout the session for long-timescale changes

@author: Dinghao Luo
"""

#%% imports
import numpy as np
import matplotlib.pyplot as plt
import sys
from pathlib import Path

repo_root = Path(__file__).resolve().parent
while not (repo_root / 'utils').exists() and repo_root != repo_root.parent:
    repo_root = repo_root.parent
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))
import imaging_pipeline_functions as ipf
import project_paths as pp
from common_functions import smooth_convolve

# plotting parameters
import matplotlib
plt.rcParams['font.family'] = 'Arial'
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42


#%% recording list
import rec_list
pathGRABNE = rec_list.pathHPCGRABNE


#%% main loop
for path in pathGRABNE:
    recname = path[-17:]
    print(recname)

    outpath = pp.GRABNE_FIGURES_STEM / 'tonic_activity' / recname
    outpath.mkdir(parents=True, exist_ok=True)

    plane_stem = pp.resolve_suite2p_plane_stem(path)
    roi_traces = np.load(plane_stem / 'F.npy', allow_pickle=True)
    roi_traces = ipf.filter_outlier(roi_traces)
    # roi_traces_dff = ipf.calculate_dFF(roi_traces)  # dFF
    tot_roi = roi_traces.shape[0]
    tot_frame = roi_traces.shape[1]

    xaxis = np.arange(tot_frame)/30

    p_per_plot = 10
    n_plot = int(np.ceil(tot_roi/p_per_plot))  # ROIs per plot


    # raw
    curr_roi = 0
    roi_range = [1, 10]
    for plot in range(n_plot):
        fig = plt.figure(1, figsize=(24,13))

        for p in range(p_per_plot):
            ax = fig.add_subplot(10,1, p+1)
            ax.plot(xaxis, roi_traces[curr_roi], c='darkgreen', lw=.5)
            ax.set(xlim=(0,xaxis[-1]),
                   title='ROI {}'.format(curr_roi+1))

            curr_roi+=1
            if curr_roi>=tot_roi:  # halt plotting if all ROIs have been plotted
                break

        ax.set(xlabel='time (s)')
        fig.tight_layout()
        plt.show(fig)

        fig.savefig(outpath / f'rois_{roi_range[0]}_{roi_range[1]}.png',
                    dpi=120, bbox_inches='tight')
        fig.savefig(outpath / f'rois_{roi_range[0]}_{roi_range[1]}.pdf',
                    bbox_inches='tight')

        roi_range = [r+10 for r in roi_range]
        if roi_range[1]>tot_roi:
            roi_range[1] = tot_roi

    # smoothed
    curr_roi = 0
    roi_range = [1, 10]
    for plot in range(n_plot):
        fig = plt.figure(1, figsize=(25,12.5))

        for p in range(p_per_plot):
            ax = fig.add_subplot(10,1, p+1)
            ax.plot(xaxis, smooth_convolve(roi_traces[curr_roi],sigma=10), c='darkgreen', lw=.5)
            ax.set(xlim=(0,xaxis[-1]),
                   title='ROI {}'.format(curr_roi+1))

            curr_roi+=1
            if curr_roi>=tot_roi:  # halt plotting if all ROIs have been plotted
                break

        ax.set(xlabel='time (s)')
        fig.tight_layout()
        plt.show(fig)

        fig.savefig(outpath / f'rois_{roi_range[0]}_{roi_range[1]}_smoothed.png',
                    dpi=120, bbox_inches='tight')
        fig.savefig(outpath / f'rois_{roi_range[0]}_{roi_range[1]}_smoothed.pdf',
                    bbox_inches='tight')

        roi_range = [r+10 for r in roi_range]
        if roi_range[1]>tot_roi:
            roi_range[1] = tot_roi
