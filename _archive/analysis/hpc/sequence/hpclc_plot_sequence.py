# -*- coding: utf-8 -*-
"""
Created on Sat 5 Aug 14:23:46 2023

plot sequence given firing rate profiles and place cell classification (from MATLAB pipeline)

@author: Dinghao Luo


"""

#%% imports
import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio
import sys
from pathlib import Path

# plotting parameters
import matplotlib
plt.rcParams['font.family'] = 'Arial'
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

repo_root = Path(__file__).resolve().parent
while not (repo_root / 'utils').exists() and repo_root != repo_root.parent:
    repo_root = repo_root.parent
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))
from common_functions import normalise
import project_paths as pp

HPC_ALL_STEM = pp.HPC_ALL_STEM
HPC_ALL_FIGURES_STEM = pp.HPC_ALL_FIGURES_STEM


def save_figure(fig, filepath, **kwargs):
    filepath.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(filepath, **kwargs)


#%% load paths to recordings
import rec_list
pathHPC = rec_list.pathHPCLCopt


#%% main
# Use firing-rate profiles from the Python pipeline (HPC all train)
# and classification results from the MATLAB preprocessing pipeline.
for pathname in pathHPC:
    recname = pathname[-17:]
    print(recname)

    classification = sio.loadmat('{}\{}_DataStructure_mazeSection1_TrialType1_FieldSpCorrAligned_Run1_Run0.mat'.format(pathname, recname))
    place_cells = classification['fieldSpCorrSessNonStimGood'][0][0]['indNeuron'][0]
    tot_pc = len(place_cells)
    if tot_pc == 0:
        print('session has no detected place cells under current criteria\n')
        continue
    print('session has {} detected place cells'.format(tot_pc))

    trains = list(np.load(HPC_ALL_STEM / recname / f'HPC_all_info_{recname}.npy',
                          allow_pickle=True).item().values())
    tot_trial = len(trains[0])

    # behaviour parameters
    info = sio.loadmat('{}\{}_DataStructure_mazeSection1_TrialType1_Info.mat'.format(pathname, recname))
    beh_info = info['beh'][0][0]
    behPar = sio.loadmat('{}\{}_DataStructure_mazeSection1_TrialType1_behPar_msess1.mat'.format(pathname, recname))
    stimOn = behPar['behPar']['stimOn'][0][0][0][1:]
    stim_trials = np.where(stimOn!=0)[0]+1
    cont_trials = stim_trials+2

    profile_cont = np.zeros((tot_pc, 5*1250))
    profile_stim = np.zeros((tot_pc, 5*1250))
    for i, cell in enumerate(place_cells):
        # take average
        temp_cont = np.zeros((len(cont_trials), 5*1250))
        temp_stim = np.zeros((len(stim_trials), 5*1250))
        for ind, trial in enumerate(cont_trials):
            trial_length = len(trains[cell-2][trial])-2500
            if trial_length<5*1250 and trial_length>0:
                temp_cont[ind, :trial_length] = trains[cell-2][trial][2500:2500+1250*5]
            elif trial_length>0:
                temp_cont[ind, :] = trains[cell-2][trial][2500:2500+1250*5]
        for ind, trial in enumerate(stim_trials):
            trial_length = len(trains[cell-2][trial])-2500
            if trial_length<5*1250 and trial_length>0:
                temp_stim[ind, :trial_length] = trains[cell-2][trial][2500:2500+5*1250]
            elif trial_length>0:
                temp_stim[ind, :] = trains[cell-2][trial][2500:2500+5*1250]

        profile_cont[i,:] = normalise(np.mean(temp_cont, axis=0))
        profile_stim[i,:] = normalise(np.mean(temp_stim, axis=0))

    # order stuff by argmax
    max_pt = {}  # argmax for conts for all pyrs
    for i in range(tot_pc):
        max_pt[i] = np.argmax(profile_cont[i,:])
    def helper(x):
        return max_pt[x]
    ord_ind = sorted(np.arange(tot_pc), key=helper)

    im_mat_cont = np.zeros((tot_pc, 5*1250))
    im_mat_stim = np.zeros((tot_pc, 5*1250))
    for i, ind in enumerate(ord_ind):
        im_mat_cont[i,:] = profile_cont[ind,:]
        im_mat_stim[i,:] = profile_stim[ind,:]

    # yticks
    if tot_pc>10:
        ytks = np.arange(0, tot_pc, 10)
        ytks[0] = 1
    elif tot_pc==10:
        ytks = [1,5,10]
    else:
        ytks = [1,5]

    # stimcont sequence
    fig, ax = plt.subplots(figsize=(3,2.2))
    image_cont = ax.imshow(im_mat_cont,
                           aspect='auto', cmap='Greys', interpolation='none',
                           extent=(-1, 4, 0, tot_pc))
    plt.colorbar(image_cont, shrink=.5)
    ax.set(yticks=ytks,
           ylabel='cell #', xlabel='time (s)',
           title='{} ctrl.'.format(recname))

    plt.show()

    # save figure
    outdirroot = HPC_ALL_FIGURES_STEM / recname
    save_figure(fig, outdirroot / 'sequence_pyr_stimcont.png',
                dpi=500, bbox_inches='tight')
    save_figure(fig, outdirroot / 'sequence_pyr_stimcont.pdf',
                bbox_inches='tight')

    plt.close(fig)

    # stim sequence ordered by stimcont
    fig, ax = plt.subplots(figsize=(3,2.2))
    image_stim = ax.imshow(im_mat_stim,
                           aspect='auto', cmap='Greys', interpolation='none',
                           extent=(-1, 4, 0, tot_pc))
    plt.colorbar(image_stim, shrink=.5)
    ax.set(yticks=ytks,
           ylabel='cell #', xlabel='time (s)',
           title='{} stim.'.format(recname))

    plt.show()

    # save figure
    outdirroot = HPC_ALL_FIGURES_STEM / recname
    save_figure(fig, outdirroot / 'sequence_pyr_stim_by_stimcont.png',
                dpi=500, bbox_inches='tight')
    save_figure(fig, outdirroot / 'sequence_pyr_stim_by_stimcont.pdf',
                bbox_inches='tight')

    plt.close(fig)
