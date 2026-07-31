# -*- coding: utf-8 -*-
'''
Created on Mon Aug  1 18:43:47 2022
update:
    25 Jan 2023, input rec_list
    11 Feb 2025, modified to include all cells

criterion: 0.33 response rate and <20 Hz

saves the average and tagged waveforms of all recordings in rec_list, pathLC
@author: Dinghao Luo


'''

#%% imports
import sys
from pathlib import Path
import numpy as np
from random import sample
import scipy.io as sio
import mat73
import os
from scipy.stats import sem

repo_root = Path(__file__).resolve().parents[2]
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))
from common_functions import normalise
from spike_text_io import param2array, get_clu
import project_paths as pp

import rec_list
paths = rec_list.pathLC

# parameters
NUM_CH = 32  # 32 channels on probe
NUM_SPK_SAMP = 32  # 32 samples per spikes
NUM_RAND_SAMP = 1000  # how many random spikes to read for waveforms
NUM_TAG_PULSE = 60  # 60 pulses for tagging
TAG_THRESHOLD = .33


#%% main function
def spk_w_sem(fspk, clu, nth_clu,
              spikes_to_load=None, n_chan=NUM_CH, n_spk_samp=NUM_SPK_SAMP):

    # get ID of every single spike of clu
    clu_n_id = np.asarray(get_clu(nth_clu, clu)).astype(int).ravel().tolist()

    rnd_samp_size = NUM_RAND_SAMP
    if spikes_to_load is None:
        if len(clu_n_id) >= rnd_samp_size:
            clu_n_id = sample(clu_n_id, rnd_samp_size)
    else:
        clu_n_id = spikes_to_load

    # load spikes (could cut out to be a separate function)
    tot_spks = len(clu_n_id)
    spks_wfs = []  # wfs for all spikes

    for i in range(tot_spks):  # reading .spk in binary ***might be finicky
        fspk.seek(clu_n_id[i]*n_chan*n_spk_samp*2)  # go to correct part of file

        spk = fspk.read(2048)  # 32*32 bts for a spike, 32*32 bts for valence of each point
        spk_fmtd = np.zeros([n_chan, n_spk_samp])  # spk but formatted as a 32x32 matrix
        for j in range(n_spk_samp):
            for k in range(n_chan):
                spk_fmtd[k, j] = spk[k*2+j*64]
                if spk[k*2+j*64+1] == 255:  # byte following value signifies valence (255 means negative)
                    spk_fmtd[k, j] = spk_fmtd[k, j] - 256  # flip sign, negative values work as subtracting from 256

        spks_wfs.append(spk_fmtd)

    spks_wfs = np.asarray(spks_wfs)

    # average spike waveforms
    av_spks = np.zeros([tot_spks, n_spk_samp])

    for i in range(tot_spks):
        spk_single = spks_wfs[i, :, :]
        av_spks[i, :] = np.nanmean(spk_single, axis=0)  # wf of averaged amplitude (channels)

    norm_spks = np.zeros([tot_spks, n_spk_samp])
    for i in range(tot_spks):
        norm_spks[i, :] = normalise(av_spks[i, :])  # normalisation

    av_spk = norm_spks.mean(0)  # 1d vector for the average tagged spike wf

    # sem calculation
    spk_sem = np.zeros(n_spk_samp)

    for i in range(n_spk_samp):
        spk_sem[i] = sem(norm_spks[:, i])

    return av_spk, spk_sem


#%% analysis
for pathname in paths:
    recname = pathname[-17:]
    print(f'\n\nProcessing {recname}')

    sess_folder = pp.LC_EPHYS_STEM / 'all_sessions' / recname
    os.makedirs(sess_folder, exist_ok=True)

    # load .mat
    mat_BTDT = sio.loadmat(
        rf'{pathname}/{recname}BTDT.mat'
        )
    behEvents = mat_BTDT['behEventsTdt']
    spInfo = sio.loadmat(
        rf'{pathname}/{recname}_DataStructure_mazeSection1_TrialType1_SpInfo_Run0'
        )
    spike_rate = spInfo['spatialInfoSess'][0][0]['meanFR'][0][0][0]

    clu = param2array(
        rf'{pathname}/{recname}.clu.1'
        )  # load .clu
    res = param2array(
        rf'{pathname}/{recname}.res.1'
        )  # load .res

    clu = np.delete(clu, 0)  # delete 1st element (noise clusters)
    all_clu = np.delete(np.unique(clu), [0, 1])
    all_clu = np.array([int(x) for x in all_clu])
    all_clu = all_clu[all_clu>=2]
    tot_clu = len(all_clu)

    with open(rf'{pathname}\{recname}.spk.1', 'rb') as fspk:  # load .spk into a byte bufferedreader

        # tagged
        stim_tp = np.zeros([NUM_TAG_PULSE, 1])  # hard-coded for LC stim protocol
        if recname=='A060r-20230530-02':
            stim_tp = np.zeros([120, 1])  # in this session there was an accidental extra 60 pulses
        tag_id = 0
        for i in range(behEvents['stimPulse'][0, 0].shape[0]):
            i = int(i)
            if behEvents['stimPulse'][0, 0][i, 3]<10:  # ~5ms tagged pulses
                temp = behEvents['stimPulse'][0, 0][i, 0]
                stim_tp[tag_id] = temp
                tag_id += 1
        if tag_id not in [NUM_TAG_PULSE, 120]:
            raise Exception('not enough tag pulses (expected 60 or 120)')

        tag_rate = np.zeros(tot_clu)
        if_tagged_spks = np.zeros([tot_clu, NUM_TAG_PULSE])
        tagged = np.zeros([tot_clu, 2])

        for iclu in range(tot_clu):
            nth_clu = iclu + 2
            clu_n_id = np.asarray(get_clu(nth_clu, clu)).astype(int).ravel().tolist()

            tagged[iclu, 0] = nth_clu

            for i in range(NUM_TAG_PULSE):  # hard-coded
                t_0 = stim_tp[i, 0]  # stim time point
                t_1 = stim_tp[i, 0] + 200  # stim time point +10ms (stricter than Takeuchi et al.)
                spks_in_range = (
                    x for x in clu_n_id if res[x].strip() and t_0 <= int(res[x]) <= t_1
                    )
                if_tagged_spks[iclu, i] = next(spks_in_range, 0)  # 1st spike in range
            tag_rate[iclu] = round(len([x for x in if_tagged_spks[iclu, :] if x > 0])/len(if_tagged_spks[iclu, :]), 3)

            # spike rate upper bound added 26 Jan 2023 to filter out non-principal cells
            if tag_rate[iclu] > TAG_THRESHOLD and spike_rate[iclu] < 20:
                tagged[iclu, 1] = 1
                print('%s%s%s%s%s' % ('clu ', nth_clu, ' tag rate = ', tag_rate[iclu], ', tagged'))
            else:
                print('%s%s%s%s' % ('clu ', nth_clu, ' tag rate = ', tag_rate[iclu]))

        waveforms = []
        for iclu in range(tot_clu):
            if tagged[iclu, 1]:
                tagged_clu = iclu+2
                spont_mean, _ = spk_w_sem(fspk, clu, tagged_clu)  # be careful that here is tagged_clu, which is iclu+2
                waveforms.append(spont_mean)

            else:
                nontagged_clu = iclu+2  # corresponds to tagged_clu above

                spont_mean, _ = spk_w_sem(fspk, clu, nontagged_clu)
                waveforms.append(spont_mean)

    # keys for saving dictionaries
    keys = [f'{recname} clu{nth_clu}' for nth_clu in range(2, tot_clu+2)]

    # save waveforms
    waveforms_dict = {keys[i]: waveforms[i] for i in range(len(keys))}
    np.save(rf'{sess_folder}\{recname}_all_waveforms.npy',
            waveforms_dict)

    # get and save ACGs
    ccg_file = mat73.loadmat(
        rf'{pathname}\{pathname[-17:]}_DataStructure_mazeSection1_TrialType1_CCG_Ctrl_Run0_mazeSess1.mat'
    )
    CCGs = ccg_file['CCGSessCtrl']['ccgVal']
    print(f'CCGs: {CCGs.shape[1]}')
    print(f'keys: {len(keys)}')
    # diagonal control CCGs are each cell's autocorrelogram
    ACGs_dict = {keys[i]: CCGs[:,i,i] for i in range(len(keys))}
    np.save(rf'{sess_folder}\{recname}_all_ACGs.npy',
            ACGs_dict)

    # save tagged identities
    tagged_dict = {keys[i]: int(tagged[i,1]) for i in range(len(keys))}  # tagged[n,0] is just the clu index
    np.save(rf'{sess_folder}\{recname}_all_identities.npy',
            tagged_dict)
