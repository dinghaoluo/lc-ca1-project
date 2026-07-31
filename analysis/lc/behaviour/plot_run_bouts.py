# -*- coding: utf-8 -*-
'''
Created on Thu Mar  2 14:33:49 2023

plot run bouts

@author: Dinghao Luo

'''

#%% imports
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.io as sio
import sys
from pathlib import Path
import mat73

repo_root = Path(__file__).resolve().parents[3]
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))
from common_functions import mpl_formatting, smooth_convolve, gaussian_kernel_unity
from console_formatting import print_session
import project_paths as pp
mpl_formatting()

import rec_list
paths = rec_list.pathLC


#%% gaussian kernel for speed smoothing
SAMP_FREQ = 1250  # Hz
gaus_speed = gaussian_kernel_unity(sigma=SAMP_FREQ*0.03)  # same as spike


#%% load cell profiles
print('loading cell profiles...')
cell_profiles = pd.read_pickle(pp.LC_EPHYS_STEM / 'LC_all_cell_profiles.pkl')


#%% main
run_bout_dir = pp.RUN_BOUTS_STEM
save_path_base = pp.RUN_BOUTS_FIGURES_STEM / 'fsa_run_bouts_plots_python'

for path in paths:
    recname = path[-17:]

    rec_stem = pp.MICEEXP_ROOT / f'ANMD{recname[1:5]}' / recname[:14] / recname
    run_bout_path = run_bout_dir / f'{recname}_run_bouts_py.csv'
    alignedRun_path = rec_stem / f'{recname}_DataStructure_mazeSection1_TrialType1_alignRun_msess1.mat'
    alignedCue_path = rec_stem / f'{recname}_DataStructure_mazeSection1_TrialType1_alignCue_msess1.mat'
    alignedRew_path = rec_stem / f'{recname}_DataStructure_mazeSection1_TrialType1_alignRew_msess1.mat'
    behave_lfp_path = rec_stem / f'{recname}_BehavElectrDataLFP.mat'
    clu_path = rec_stem / f'{recname}.clu.1'
    res_path = rec_stem / f'{recname}.res.1'
    print_session(recname)

    # load data
    run_bout_table = pd.read_csv(run_bout_path)

    beh_lfp = mat73.loadmat(behave_lfp_path)
    tracks = beh_lfp['Track']
    laps = beh_lfp['Laps']

    aligned = sio.loadmat(alignedRun_path)['trialsRun'][0][0]
    alignedCue = sio.loadmat(alignedCue_path)['trialsCue'][0][0]  # for cue marking
    alignedRew = sio.loadmat(alignedRew_path)['trialsRew'][0][0]  # for cue marking

    # read cues
    cueLfpInd = alignedCue['startLfpInd'].flatten()

    # read rewards
    rewLfpInd = alignedRew['startLfpInd'].flatten()

    # spike reading
    clusters = np.loadtxt(clu_path, dtype=int, skiprows=1)  # first line = number of clusters
    spike_times = np.loadtxt(res_path, dtype=int) / (20_000 / 1_250)  # convert to behavioural time scale
    spike_times = spike_times.astype(int)  # integer indices for indexing

    unique_clus = [clu for clu in np.unique(clusters) if clu not in [0, 1]]
    clu_to_row = {clu: i for i, clu in enumerate(unique_clus)}  # map cluster ID to row index

    max_time = spike_times.max() + 1  # +1 to make sure last time index is included
    spike_map = np.zeros((len(unique_clus), max_time), dtype=int)
    spike_array = np.zeros((len(unique_clus), max_time))

    for time, clu in zip(spike_times, clusters):
        if clu in [0, 1]:
            continue  # skip noise
        row = clu_to_row[clu]
        spike_map[row, time] = 1  # set spike bin to 1

    for i in range(len(unique_clus)):
        spike_array[i, :] = smooth_convolve(spike_map[i, :], sigma=int(SAMP_FREQ*0.1))

    # read licks
    lickLfp = laps['lickLfpInd']
    lickLfp_flat = []
    for trial in range(len(lickLfp)):
        if isinstance(lickLfp[trial][0], np.ndarray):
            for i in range(len(lickLfp[trial][0])):
                lickLfp_flat.append(int(lickLfp[trial][0][i]))
    lickLfp_flat = np.array(lickLfp_flat)
    speed_MMsec = tracks['speed_MMsecAll']
    for tbin in range(len(speed_MMsec)):
        if speed_MMsec[tbin]<0:
            speed_MMsec[tbin] = (speed_MMsec[tbin-1]+speed_MMsec[tbin+1])/2
    speed_MMsec = np.convolve(speed_MMsec, gaus_speed, mode='same')/10  # /10 for cm
    startLfpInd = aligned['startLfpInd'][0]
    endLfpInd = aligned['endLfpInd'][0]

    # identify cells of interest
    selected_indices = []
    selected_names = []

    for i, clu in enumerate(unique_clus):
        cluname = f'{recname} clu{clu}'
        profile = cell_profiles.loc[cluname]
        if profile['identity'] in ['tagged', 'putative'] and profile['run_onset_peak'] is True:
            selected_indices.append(i)
            selected_names.append(cluname)

    if not selected_indices:
        print(f'no tagged/putative run-onset peak cells found for {recname}')
        continue

    print(f'{len(selected_indices)} tagged/putative run-onset peak cells selected for {recname}')

    ## plotting
    save_path_sess = save_path_base / recname

    # plot mean firing rate of selected cells across trial windows
    for t in np.arange(2, len(endLfpInd)-2, 3):
        lfp_indices_t = np.arange(startLfpInd[t]-1250, min(endLfpInd[t+2]+1250, len(speed_MMsec)))
        lap_start = lfp_indices_t[0]
        xaxis = np.arange(0, len(lfp_indices_t)) / SAMP_FREQ

        fig, ax = plt.subplots(figsize=(len(lfp_indices_t)/5000, 1.2))
        ax.set(xlabel='time (s)', ylabel='speed (cm/s)',
               ylim=(0, 1.2 * max(speed_MMsec[lfp_indices_t])),
               xlim=(0, len(lfp_indices_t) / SAMP_FREQ),
               title=f'{recname} trials {t-1} to {t+1}')

        ax.plot(xaxis, speed_MMsec[lfp_indices_t], color='royalblue', label='speed')

        ax_spk = ax.twinx()
        ax_spk.set_ylabel('firing rate (Hz)', color='black')
        ax_spk.spines['right'].set_color('black')
        ax_spk.tick_params(axis='y', colors='black', labelcolor='black')

        # compute mean trace
        spike_subset = spike_array[selected_indices, :][:, lfp_indices_t] * SAMP_FREQ
        mean_trace = np.mean(spike_subset, axis=0)
        mean_trace = np.clip(mean_trace, 0, np.percentile(mean_trace, 99.5))

        ax_spk.plot(xaxis, mean_trace, color='black', linewidth=1.2, label='mean spike')
        ax_spk.set_ylim(0, np.max(mean_trace) * 1.1)

        cueLfpInd_t = cueLfpInd[np.isin(cueLfpInd, lfp_indices_t)]
        ax.vlines((cueLfpInd_t - lap_start)/SAMP_FREQ, 0, ax.get_ylim()[1], 'darkgrey', zorder=10)

        rewLfpInd_t = rewLfpInd[np.isin(rewLfpInd, lfp_indices_t)]
        ax.vlines((rewLfpInd_t - lap_start)/SAMP_FREQ, ax.get_ylim()[1], ax.get_ylim()[1]*.95, 'forestgreen', linewidth=1.5, zorder=10)

        startLfpInd_t = startLfpInd[np.isin(startLfpInd, lfp_indices_t)]
        ax.vlines((startLfpInd_t - lap_start)/SAMP_FREQ, 0, ax.get_ylim()[1], 'red', linestyle='dashed', zorder=10)

        run_bout_t = run_bout_table.iloc[:,1][np.isin(run_bout_table.iloc[:,1], lfp_indices_t)]
        ax.vlines((run_bout_t - lap_start)/SAMP_FREQ, 0, ax.get_ylim()[1], 'green', linestyle='dashed', zorder=10)

        licks_t = lickLfp_flat[np.isin(lickLfp_flat, lfp_indices_t)]
        ax.vlines((licks_t - lap_start)/SAMP_FREQ, ax.get_ylim()[1], ax.get_ylim()[1] * 0.96, 'magenta')

        ax.spines['top'].set_visible(False)
        ax_spk.spines['top'].set_visible(False)

        save_path_sess.mkdir(parents=True, exist_ok=True)
        rb_tag = ' rb' if len(run_bout_t) > 0 else ''
        save_path = save_path_sess / f'trials_{t-1}_to_{t+1}{rb_tag}'

        for ext in ['.pdf', '.png']:
            fig.savefig(f'{save_path}{ext}', dpi=300, bbox_inches='tight')
        plt.close(fig)
