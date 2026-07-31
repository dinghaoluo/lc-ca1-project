# -*- coding: utf-8 -*-
'''
Created on Tue Sep  9 13:12:32 2025

stimulated dopamine release vs endogenous dopamine release

@author: Dinghao Luo
'''

#%% imports
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

repo_root = Path(__file__).resolve().parents[3]
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))

from common_functions import mpl_formatting
from console_formatting import print_files_saved, print_session, print_status
mpl_formatting()

import imaging_pipeline_functions as ipf

import behaviour_functions as bf
import project_paths as pp

import rec_list
paths = rec_list.pathdLightLCOpto


#%% params
SAMP_FREQ = 30
BEF = 1
AFT = 4
TAXIS = np.arange(-BEF * SAMP_FREQ, AFT * SAMP_FREQ) / SAMP_FREQ
WIN_LEN = int((BEF + AFT) * SAMP_FREQ)
PMT_BUFFER_FRAMES = 10

all_sess_stem = pp.HPC_DLIGHT_LC_OPTO_STEM / 'all_sessions'
all_miceexp_stem = pp.MICEEXP_ROOT


#%% iterate sessions
for path in paths:
    recname = Path(path).name
    print_session(recname)

    animal = f'ANMD{recname[1:4]}'
    recname_txt = recname.replace('i', '')

    binpath = Path(path) / 'suite2p' / 'plane0' / 'data.bin'
    opspath = Path(path) / 'suite2p' / 'plane0' / 'ops.npy'
    txtpath = all_miceexp_stem / animal / f'{recname_txt}T.txt'

    savepath = all_sess_stem / recname
    (savepath / 'processed_data').mkdir(parents=True, exist_ok=True)

    out_fig = (
        pp.HPC_DLIGHT_LC_OPTO_ALL_SESSIONS_FIGURES_STEM
        / recname
        / f'{recname}_dLight_profiles_stim_ctrl'
        )
    out_fig.parent.mkdir(parents=True, exist_ok=True)

    # behaviour check
    beh = bf.process_behavioural_data_imaging(txtpath)
    run_onsets = beh['run_onset_frames']
    stim_methods = [t[15] for t in beh['trial_statements']]

    if len(run_onsets) == 0:
        print_status('skipped', 'no behaviour in session')
        continue
    elif '2' not in stim_methods:
        print_status('skipped', 'no valid stim in session')
        continue
    else:
        print_status('analysing', 'behaviour session')

    # load data
    print_status('loading movie')
    ops = np.load(opspath, allow_pickle=True).item()
    tot_frames = int(ops['nframes'])
    Ly, Lx = int(ops['Ly']), int(ops['Lx'])
    shape = (tot_frames, Ly, Lx)
    mov = np.memmap(binpath, mode='r', dtype='int16', shape=shape)

    # stim info
    print_status('extracting stim info')
    frame_times = beh['frame_times']
    pulse_times = beh['pulse_times']
    pulse_params = [ls for ls in beh['pulse_descriptions'] if ls]

    pulse_width    = float(pulse_params[-1][4]) / 1_000_000

    diffs = np.diff(pulse_times)
    split_idx = np.where(diffs >= 1000)[0] + 1
    pulse_trains = np.split(pulse_times, split_idx)

    pulse_frames = [
        [ipf.find_nearest(p, frame_times) for p in train]
        for train in pulse_trains
    ]

    stim_mask = np.zeros(tot_frames, dtype=bool)
    for train in pulse_frames:
        start = train[0] - 1  # one-frame pre-buffer
        end = train[-1] + int(pulse_width * SAMP_FREQ) + PMT_BUFFER_FRAMES
        if start < 0 or end > tot_frames:
            raise ValueError(f'Pulse train extends beyond movie frames in {recname}')
        stim_mask[start:end] = True

    # split trials
    ctrl_trials = [i for i, f in enumerate(stim_methods) if f == '0']
    stim_trials = [i for i, f in enumerate(stim_methods) if f != '0']

    # only complete run-onset windows enter the comparison
    ctrl_frames = []
    for i in ctrl_trials:
        f = run_onsets[i]
        if f is not None and not np.isnan(f):
            start = int(f) - int(BEF * SAMP_FREQ)
            end   = int(f) + int(AFT * SAMP_FREQ)
            if start >= 0 and end <= tot_frames:
                ctrl_frames.append(int(f))

    stim_frames = []
    for i in stim_trials:
        f = run_onsets[i]
        if f is not None and not np.isnan(f):
            start = int(f) - int(BEF * SAMP_FREQ)
            end   = int(f) + int(AFT * SAMP_FREQ)
            if start >= 0 and end <= tot_frames:
                stim_frames.append(int(f))

    if not ctrl_frames or not stim_frames:
        print_status('skipped', 'no complete ctrl. or stim. run-onset windows')
        continue

    # extract data
    print_status('extracting fluorescence data')
    F = np.sum(mov, axis=(1, 2)).astype(np.float32)
    dFF = ipf.calculate_dFF(F, t_axis=0)
    dFF[stim_mask] = np.nan  # mask out the stim frames

    ctrl_profiles = np.zeros((len(ctrl_frames), WIN_LEN), dtype=np.float32)
    stim_profiles = np.zeros((len(stim_frames), WIN_LEN), dtype=np.float32)
    for j, f in enumerate(ctrl_frames):
        start = f - int(BEF * SAMP_FREQ)
        end   = f + int(AFT * SAMP_FREQ)
        ctrl_profiles[j, :] = dFF[start : end]
    for j, f in enumerate(stim_frames):
        start = f - int(BEF * SAMP_FREQ)
        end   = f + int(AFT * SAMP_FREQ)
        stim_profiles[j, :] = dFF[start : end]

    # single session plot
    ctrl_mean = np.nanmean(ctrl_profiles, axis=0)
    stim_mean = np.nanmean(stim_profiles, axis=0)

    print_status('plotting')
    fig, ax = plt.subplots(figsize=(2.4, 2))
    ax.plot(TAXIS, ctrl_mean,
            c='grey', label='ctrl.', linewidth=1)
    ax.plot(TAXIS, stim_mean,
            c='royalblue', label='stim.', linewidth=1)
    ax.axvline(0, linestyle='--', linewidth=1, color='red')
    ax.set(xlabel='time from run onset (s)',
           ylabel='ΔF/F',
           title=recname)
    ax.legend(frameon=False)

    for s in ['top', 'right']:
        ax.spines[s].set_visible(False)
    fig.tight_layout()

    saved_paths = []
    for ext in ['.png', '.pdf']:
        curr_path = f'{out_fig}{ext}'
        fig.savefig(curr_path,
                    dpi=300, bbox_inches='tight')
        saved_paths.append((ext[1:], curr_path))
    print_files_saved(saved_paths)
    plt.close(fig)
