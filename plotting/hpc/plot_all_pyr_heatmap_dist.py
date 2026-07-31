# -*- coding: utf-8 -*-
'''
Created on Tue 26 Sep 14:05:31 2023

pyramidal cell heatmap for each session

@author: Dinghao Luo
'''

#%% imports
import numpy as np
import matplotlib.pyplot as plt
import mat73
import scipy.io as sio
import sys
from pathlib import Path


#%% load local helpers and paths
repo_root = Path(__file__).resolve().parents[2]
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))
from common_functions import normalise
import project_paths as pp

import rec_list
pathHPC = rec_list.pathHPCLCopt
OUTPUT_DIR = pp.HPC_ALL_FIGURES_STEM / 'HPC_cont_stim_sequence_heatmap'


#%% main
for pathname in pathHPC:
    pathname = Path(pathname)
    recname = pathname.name

    # spikes by distance
    spikes_dist_path = (
        pathname / f'{recname}_DataStructure_mazeSection1_TrialType1_convSpikesDist20mm_Run1.mat'
        )
    spikes_dist = mat73.loadmat(str(spikes_dist_path))['filteredSpikeArray']
    for trial_spikes in spikes_dist:
        if trial_spikes is not None:
            tot_dist = trial_spikes.shape[1]  # in mm
            break
    else:
        raise ValueError('no valid distance-spike trials')
    tot_trial = len(spikes_dist)

    # if each neurones is an interneurone or a pyramidal cell
    info = sio.loadmat(
        pathname / f'{recname}_DataStructure_mazeSection1_TrialType1_Info.mat'
        )
    rec_info = info['rec'][0][0]
    intern_id = rec_info['isIntern'][0]
    pyr_id = [not(clu) for clu in intern_id]

    # behaviour parameters
    beh_info = info['beh'][0][0]
    stim_ind = np.where(beh_info['pulseMethod']!=0)[1]
    stim_start = stim_ind[0]; stim_end = stim_ind[-1]
    stim_trials = np.arange(stim_start, stim_end)

    tot_clu = len(pyr_id)
    tot_pyr = sum(pyr_id)  # how many pyramidal cells are in this recording
    pyr_mat = np.zeros((tot_pyr, tot_dist))
    pyr_mat_cont = np.zeros((tot_pyr, tot_dist))
    pyr_mat_stim = np.zeros((tot_pyr, tot_dist))

    pyr_counter = 0
    argmax_pyr = []; argmax_pyr_cont = []; argmax_pyr_stim = []
    for i in range(tot_clu):
        if pyr_id[i]==True:
            temp = np.zeros((tot_trial, tot_dist))  # temporary to contain all trials of one clu
            temp_cont = np.zeros((stim_start, tot_dist))
            temp_stim = np.zeros((len(stim_trials), tot_dist))
            for trial in range(tot_trial):
                trial_spikes = spikes_dist[trial]
                if trial_spikes is not None:
                    temp[trial, :] = trial_spikes[i, :]
            for trial in range(stim_start):
                trial_spikes = spikes_dist[trial]
                if trial_spikes is not None:
                    temp_cont[trial, :] = trial_spikes[i, :]
            for ind, trial in enumerate(stim_trials):
                trial_spikes = spikes_dist[trial]
                if trial_spikes is not None:
                    temp_stim[ind, :] = trial_spikes[i, :]
            pyr_mat[pyr_counter, :] = normalise(np.mean(temp, axis=0))
            pyr_mat_cont[pyr_counter, :] = normalise(np.mean(temp_cont, axis=0))
            pyr_mat_stim[pyr_counter, :] = normalise(np.mean(temp_stim, axis=0))
            argmax_pyr.append(np.argmax(pyr_mat[pyr_counter, :]))
            argmax_pyr_cont.append(np.argmax(pyr_mat_cont[pyr_counter, :]))
            argmax_pyr_stim.append(np.argmax(pyr_mat_stim[pyr_counter, :]))

            pyr_counter+=1

    temp = list(np.arange(tot_pyr))
    peak_ordered, pyr_id_ordered = zip(*sorted(zip(argmax_pyr, temp)))
    pyr_mat_ordered = np.zeros((tot_pyr, tot_dist))
    for i in range(tot_pyr):
        curr_id = pyr_id_ordered[i]
        pyr_mat_ordered[i,:] = pyr_mat[curr_id,:]

    peak_ordered_cont, pyr_id_ordered_cont = zip(*sorted(zip(argmax_pyr_cont, temp)))
    pyr_mat_ordered_cont = np.zeros((tot_pyr, tot_dist))
    for i in range(tot_pyr):
        curr_id = pyr_id_ordered_cont[i]  # use the same order
        pyr_mat_ordered_cont[i,:] = pyr_mat_cont[curr_id,:]

    peak_ordered_stim, pyr_id_ordered_stim = zip(*sorted(zip(argmax_pyr_stim, temp)))
    pyr_mat_ordered_stim = np.zeros((tot_pyr, tot_dist))
    for i in range(tot_pyr):
        curr_id = pyr_id_ordered_stim[i]  # use the same order
        pyr_mat_ordered_stim[i,:] = pyr_mat_stim[curr_id,:]

    # plot heatmap
    fig, axs = plt.subplots(1,2, figsize=(8,4))
    image_cont = axs[0].imshow(pyr_mat_ordered_cont, aspect='auto', cmap='turbo',
                               extent=[0, 180.1, 0, tot_pyr])
    plt.colorbar(image_cont)
    axs[0].set(title='control')

    image_stim = axs[1].imshow(pyr_mat_ordered_stim, aspect='auto', cmap='turbo',
                               extent=[0, 180.1, 0, tot_pyr])
    plt.colorbar(image_stim)
    axs[1].set(title='stim')

    for i in range(2):
        axs[i].set(xlabel='distance (cm)', ylabel='pyr cell #')

    fig.suptitle(recname)

    plt.tight_layout()
    plt.show()

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT_DIR / f'{recname}.png',
                dpi=500, bbox_inches='tight')

    plt.close(fig)
