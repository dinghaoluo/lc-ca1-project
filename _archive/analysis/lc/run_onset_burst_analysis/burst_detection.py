# -*- coding: utf-8 -*-
"""
Created on Thu Mar 13 11:53:25 2025

burst analysis for LC cells

@author: Dinghao Luo

- extracted generic burst-detection helpers into utils/burst_detection_functions.py
"""

#%% imports
import numpy as np
import sys
import pandas as pd
import matplotlib.pyplot as plt
import os
from pathlib import Path

repo_root = Path(__file__).resolve().parent
while not (repo_root / 'utils').exists() and repo_root != repo_root.parent:
    repo_root = repo_root.parent
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))
from common_functions import mpl_formatting
from burst_detection_functions import detect_bursts_from_isi, percentile_based_thresholds
import project_paths as pp
mpl_formatting()

import rec_list
paths = rec_list.pathLC


#%% load dataframe
df = pd.read_pickle(
    pp.LC_EPHYS_STEM / 'LC_all_cell_profiles.pkl'
    )


#%% functions
def plot_bursts(spike_times,
                bursts,
                cluname,
                save_dir,
                window_size=3,
                max_bursts_per_fig=50):
    """
    plots spikes with bursts highlighted in black and tonic spikes in grey, centering each burst on its own subplot.
    if a cell has too many bursts, splits them into multiple figures.

    parameters:
    - spike_times: list or numpy array of spike timestamps in seconds
    - bursts: list of tuples (burst_start_index, burst_end_index) from detect_bursts_from_isi
    - cluname: string, name of the cell (used for figure title)
    - save_dir: directory where figures will be saved
    - window_size: time window around each burst center (default: 3 seconds)
    - max_bursts_per_fig: maximum number of bursts per figure (default: 50)

    returns:
    - None (saves the figure)
    """
    if len(spike_times) == 0 or len(bursts) == 0:
        print(f'no spikes or bursts to plot for {cluname}.')
        return

    # save directory
    save_dir = rf'{save_dir}'

    spike_times = np.array(spike_times)

    num_total_bursts = len(bursts)
    num_figs = int(np.ceil(num_total_bursts / max_bursts_per_fig))

    for fig_idx in range(num_figs):
        start_burst = fig_idx * max_bursts_per_fig
        end_burst = min((fig_idx + 1) * max_bursts_per_fig, num_total_bursts)
        bursts_subset = bursts[start_burst:end_burst]

        num_bursts = len(bursts_subset)
        fig, axes = plt.subplots(num_bursts, 1, figsize=(5, 0.5 * num_bursts), sharex=True)

        if num_bursts == 1:
            axes = [axes]

        for i, (start_idx, end_idx) in enumerate(bursts_subset):
            ax = axes[i]

            burst_center = spike_times[(start_idx + end_idx) // 2]
            start_time = burst_center - window_size / 2
            end_time = burst_center + window_size / 2

            spikes_in_window = spike_times[(spike_times >= start_time) & (spike_times < end_time)]
            burst_spikes = spike_times[start_idx:end_idx + 1]
            tonic_spikes = np.setdiff1d(spikes_in_window, burst_spikes)

            ax.scatter(tonic_spikes - burst_center, np.ones_like(tonic_spikes), color='grey', s=10, alpha=0.6)
            ax.scatter(burst_spikes - burst_center, np.ones_like(burst_spikes), color='black', s=10)

            ax.set_xlim(-window_size / 2, window_size / 2)
            ax.set_yticks([])

            for spine in ['top', 'right', 'left']:
                ax.spines[spine].set_visible(False)

            if i == num_bursts - 1:
                ax.set_xlabel('time from burst (s)', fontsize=8)
            else:
                ax.set_xticklabels([])

            ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, _: f'{x:.1f}'))
            ax.tick_params(axis='x', labelsize=8)

        fig.suptitle(f'{cluname} (part {fig_idx + 1})', fontsize=10)
        plt.tight_layout()

        # generate and save figure
        plt.savefig(save_dir, dpi=300, bbox_inches='tight')
        plt.close()


#%% main
for path in paths:
    recname = path[-17:]
    print(f'processing {recname}')

    sess_folder = pp.LC_EPHYS_STEM / 'all_sessions' / recname
    burst_folder = pp.LC_EPHYS_FIGURES_STEM / 'single_cell_bursts'
    os.makedirs(burst_folder, exist_ok=True)

    ISI_dict = np.load(sess_folder / f'{recname}_all_ISIs.npy', allow_pickle=True).item()
    spike_dict = np.load(sess_folder / f'{recname}_all_spikes.npy', allow_pickle=True).item()

    burst_dict = {}

    for cluname, ISIs in ISI_dict.items():
        identity = df.loc[cluname]['identity']
        if identity != 'other':
            burst_thresh, end_thresh = percentile_based_thresholds(ISIs)
            bursts = detect_bursts_from_isi(ISIs, burst_thresh, end_thresh)

            spike_times = [t/20000 for t in spike_dict[cluname]]
            save_path = os.path.join(burst_folder, f'{cluname}_{identity}.png')

            plot_bursts(spike_times, bursts, cluname, save_path)

            burst_dict[cluname] = bursts

    # save burst dictionary for session
    np.save(sess_folder / f'{recname}_all_bursts.npy',
            burst_dict)
