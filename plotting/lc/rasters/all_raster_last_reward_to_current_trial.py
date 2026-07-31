# -*- coding: utf-8 -*-
'''
Created on Wed May 24 17:09:18 2023
Originally named all_raster_cue_rew_run_lasttocurtr.py

save all rasters into a big .npy file

@author: Dinghao Luo
'''

#%% imports
import numpy as np
import scipy.io as sio
import h5py
import sys
from pathlib import Path

repo_root = Path(__file__).resolve().parents[3]
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))
from console_formatting import print_session
import project_paths as pp

import rec_list
pathLC = rec_list.pathLC
pp.LC_ALL_STEM.mkdir(parents=True, exist_ok=True)


#%% create rasters
all_rasters_cue = {}

for pathname in pathLC:
    recname = Path(pathname).name
    print_session(recname)

    filename = str(Path(pathname) / f'{recname}_DataStructure_mazeSection1_TrialType1')

    # check if this is an opto-stim session
    beh_info_file = sio.loadmat(filename+'_Info.mat')['beh']
    stim_trial = np.squeeze(np.where(beh_info_file['pulseMethod'][0][0][0]!=0))
    if stim_trial.size!=0:
        stimtype = beh_info_file['pulseMethod'][0][0][0][stim_trial[0]]
        stimwind = [stim_trial[0], stim_trial[-1]]
    else:
        stimtype = 'NA'
        stimwind = ['NA', 'NA']

    spike_time_file = h5py.File(filename + '_alignedSpikesPerNPerT_msess1_Run0.mat')['trialsCueSpikes']
    spike_struct = spike_time_file['Time_LasttoCurTr']

    tot_clu = spike_struct.shape[1]
    tot_trial = spike_struct.shape[0]  # trial 1 is empty but tot_trial includes it for now
    all_id = list(np.arange(2, tot_clu+2))

    # spike reading
    spike_time = np.empty(shape=(tot_clu, tot_trial - 1), dtype='object')
    for i in range(tot_clu):
        i = int(i)-2
        for j in range(1, tot_trial):  # trial 0 is an empty trial
            spike_time[i,j-1] = spike_time_file[spike_struct[j,i]][0]

    max_length = 16250  # 3 s before, 10 s after
    spike_train_all = np.empty(shape=(tot_clu, tot_trial-1), dtype='object')
    for i in range(tot_clu):
        for trial in range(tot_trial-1):
            clu_id = int(all_id[i])-2
            spikes = spike_time[clu_id][trial].reshape(-1)
            spike_train_all[i, trial] = spikes

    # save into all_rasters
    i = 0
    for clu in all_id:
        cluname = recname+' clu'+str(clu)
        all_rasters_cue[cluname] = spike_train_all[i]
        i+=1

# at this point all_rasters should have rasters of all tagged cells from only
# the opto-stim sessions


#%% save all rasters
print('\nWriting to LC_all_rasters_cue.npy\n')
np.save(pp.LC_ALL_STEM / 'LC_all_rasters_cue.npy',
        all_rasters_cue)
print('success!')


#%% create rasters
all_rasters_run = {}

for pathname in pathLC:
    recname = Path(pathname).name
    print_session(recname)

    filename = str(Path(pathname) / f'{recname}_DataStructure_mazeSection1_TrialType1')

    # check if this is an opto-stim session
    beh_info_file = sio.loadmat(filename+'_Info.mat')['beh']
    stim_trial = np.squeeze(np.where(beh_info_file['pulseMethod'][0][0][0]!=0))
    if stim_trial.size!=0:
        stimtype = beh_info_file['pulseMethod'][0][0][0][stim_trial[0]]
        stimwind = [stim_trial[0], stim_trial[-1]]
    else:
        stimtype = 'NA'
        stimwind = ['NA', 'NA']

    spike_time_file = h5py.File(filename + '_alignedSpikesPerNPerT_msess1_Run0.mat')['trialsRunSpikes']
    spike_struct = spike_time_file['Time_LasttoCurTr']

    tot_clu = spike_struct.shape[1]
    tot_trial = spike_struct.shape[0]  # trial 1 is empty but tot_trial includes it for now
    all_id = list(np.arange(2, tot_clu+2))

    # spike reading
    spike_time = np.empty(shape=(tot_clu, tot_trial - 1), dtype='object')
    for i in range(tot_clu):
        i = int(i)-2
        for j in range(1, tot_trial):  # trial 0 is an empty trial
            spike_time[i,j-1] = spike_time_file[spike_struct[j,i]][0]

    max_length = 16250  # 3 s before, 10 s after
    spike_train_all = np.empty(shape=(tot_clu, tot_trial-1), dtype='object')
    for i in range(tot_clu):
        for trial in range(tot_trial-1):
            clu_id = int(all_id[i])-2
            spikes = spike_time[clu_id][trial].reshape(-1)
            spike_train_all[i, trial] = spikes

    # save into all_rasters
    i = 0
    for clu in all_id:
        cluname = recname+' clu'+str(clu)
        all_rasters_run[cluname] = spike_train_all[i]
        i+=1

# at this point all_rasters should have rasters of all tagged cells from only
# the opto-stim sessions


#%% save all rasters
print('\nWriting to LC_all_rasters_run.npy\n')
np.save(pp.LC_ALL_STEM / 'LC_all_rasters_run.npy',
        all_rasters_run)
print('success!')


#%% create rasters
all_rasters_rew = {}

for pathname in pathLC:
    recname = Path(pathname).name
    print_session(recname)

    filename = str(Path(pathname) / f'{recname}_DataStructure_mazeSection1_TrialType1')

    # check if this is an opto-stim session
    beh_info_file = sio.loadmat(filename+'_Info.mat')['beh']
    stim_trial = np.squeeze(np.where(beh_info_file['pulseMethod'][0][0][0]!=0))
    if stim_trial.size!=0:
        stimtype = beh_info_file['pulseMethod'][0][0][0][stim_trial[0]]
        stimwind = [stim_trial[0], stim_trial[-1]]
    else:
        stimtype = 'NA'
        stimwind = ['NA', 'NA']

    spike_time_file = h5py.File(filename + '_alignedSpikesPerNPerT_msess1_Run0.mat')['trialsRewSpikes']
    spike_struct = spike_time_file['Time_LasttoCurTr']

    tot_clu = spike_struct.shape[1]
    tot_trial = spike_struct.shape[0]  # trial 1 is empty but tot_trial includes it for now
    all_id = list(np.arange(2, tot_clu+2))

    # spike reading
    spike_time = np.empty(shape=(tot_clu, tot_trial - 1), dtype='object')
    for i in range(tot_clu):
        i = int(i)-2
        for j in range(1, tot_trial):  # trial 0 is an empty trial
            spike_time[i,j-1] = spike_time_file[spike_struct[j,i]][0]

    max_length = 16250  # 3 s before, 10 s after
    spike_train_all = np.empty(shape=(tot_clu, tot_trial-1), dtype='object')
    for i in range(tot_clu):
        for trial in range(tot_trial-1):
            clu_id = int(all_id[i])-2
            spikes = spike_time[clu_id][trial].reshape(-1)
            spike_train_all[i, trial] = spikes

    # save into all_rasters
    i = 0
    for clu in all_id:
        cluname = recname+' clu'+str(clu)
        all_rasters_rew[cluname] = spike_train_all[i]
        i+=1

# at this point all_rasters should have rasters of all tagged cells from only
# the opto-stim sessions


#%% save all rasters
print('\nWriting to LC_all_rasters_rew.npy\n')
np.save(pp.LC_ALL_STEM / 'LC_all_rasters_rew.npy',
        all_rasters_rew)
print('success!')
