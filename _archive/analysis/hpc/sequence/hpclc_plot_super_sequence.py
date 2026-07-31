# -*- coding: utf-8 -*-
"""
Created on Sat 5 Aug 14:23:46 2023
Modified to plot super-sequences

plot sequence given firing rate profiles and place cell classification (from MATLAB pipeline)
supersequence

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
pathHPC = rec_list.pathHPCLCopt+rec_list.pathHPCLCtermopt


#%% main
# Use firing-rate profiles from the Python pipeline (HPC all train)
# and classification results from the MATLAB preprocessing pipeline.
all_profiles = np.zeros((1, 5*1250))
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
    stim_start = np.where(stimOn!=0)[0][0]

    profile = np.zeros((tot_pc, 5*1250))
    for i, cell in enumerate(place_cells):
        # take average
        temp = np.zeros((stim_start, 5*1250))
        for trial in range(stim_start):
            trial_length = len(trains[cell-2][trial])-2500
            if trial_length<5*1250 and trial_length>0:
                temp[trial, :trial_length] = trains[cell-2][trial][2500:2500+1250*5]
            elif trial_length>0:
                temp[trial, :] = trains[cell-2][trial][2500:2500+1250*5]

        profile[i,:] = normalise(np.mean(temp, axis=0))

    all_profiles = np.vstack((all_profiles, profile))


#%% remove the first row which is 0's
all_profiles = np.delete(all_profiles, (0), axis=0)


#%% plotting
# order stuff by argmax
max_pt = {}  # argmax for conts for all pyrs
for i in range(all_profiles.shape[0]):
    max_pt[i] = np.argmax(all_profiles[i,1250:])  # 1250: because we want to look at from 0 to 4
def helper(x):
    return max_pt[x]
ord_ind = sorted(np.arange(all_profiles.shape[0]), key=helper)

im_mat = np.zeros(all_profiles.shape)
for i, ind in enumerate(ord_ind):
    im_mat[i,:] = all_profiles[ind,:]

# yticks
ytks = np.arange(0, all_profiles.shape[0], 100)
ytks[0] = 1

# stimcont sequence
fig, ax = plt.subplots(figsize=(2.1,1.6))
image = ax.imshow(im_mat,
                  aspect='auto', cmap='Greys', interpolation='none', origin='lower',
                  extent=(-1, 4, 0, all_profiles.shape[0]))
plt.colorbar(image, shrink=.5)
ax.set(yticks=ytks,
       ylabel='cell #', xlabel='time (s)',
       title='supersequence_all')

plt.show()

# save figure
save_figure(fig, HPC_ALL_FIGURES_STEM / 'HPCLC_HPCLCterm_supersequence.png',
            dpi=500, bbox_inches='tight')
save_figure(fig, HPC_ALL_FIGURES_STEM / 'HPCLC_HPCLCterm_supersequence.pdf',
            bbox_inches='tight')

plt.close(fig)
