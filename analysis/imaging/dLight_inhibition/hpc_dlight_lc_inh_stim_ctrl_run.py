# -*- coding: utf-8 -*-
'''
Created on Wed Aug 13 16:25:06 2025

focussing on the inhibition experiments, investigate whether single-grid run
    onset aligned average activity changes between ctrl and stim trials

@author: Dinghao Luo
'''

#%% imports
import os
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
from common_functions import mpl_formatting, get_GPU_availability
from console_formatting import print_files_saved, print_session
import project_paths as pp
mpl_formatting()

import imaging_pipeline_functions as ipf

import behaviour_functions as bf

import rec_list
paths = rec_list.pathdLightLCOptoInh


#%% parameters
SAMP_FREQ = 30  # Hz
BORDER = 6  # in pixels
CELL_SIZE = 50  # in pixels

TAXIS = np.arange(-3*SAMP_FREQ, 5*SAMP_FREQ) / SAMP_FREQ


#%% GPU acceleration
_, GPU_AVAILABLE, gpu_name = get_GPU_availability()
if GPU_AVAILABLE:
    print(f'using GPU-acceleration with {gpu_name} and CuPy')


#%% main
for path in paths[-1:]:
    recname = Path(path).name
    print_session(recname)

    plane_stem = pp.resolve_suite2p_session_stem(path) / 'plane0'
    binpath = plane_stem / 'data.bin'
    bin2path = plane_stem / 'data_chan2.bin'
    opspath = plane_stem / 'ops.npy'
    txtpath = (
        pp.MICEEXP_ROOT
        / f'ANMD{recname[1:4]}'
        / f'{recname[:4]}{recname[5:]}T.txt'
    )

    savepath = pp.HPC_DLIGHT_LC_OPTO_ALL_SESSIONS_FIGURES_STEM / recname
    os.makedirs(savepath, exist_ok=True)

    # load data
    ops = np.load(opspath, allow_pickle=True).item()
    tot_frames = ops['nframes']
    shape = tot_frames, ops['Ly'], ops['Lx']

    print('loading movies...')
    mov = np.memmap(binpath, mode='r', dtype='int16', shape=shape).astype(np.float32)
    mov2 = np.memmap(bin2path, mode='r', dtype='int16', shape=shape).astype(np.float32)

    tot_frames = mov.shape[0]

    # crop the movies
    mov_c = mov[:, BORDER : ops['Ly']-BORDER, BORDER : ops['Lx']-BORDER]
    mov2_c = mov2[:, BORDER : ops['Ly']-BORDER, BORDER : ops['Lx']-BORDER]

    # grid-wise sums
    print('computing grid-wise sums...')
    # (time, 500, 500) -> (100, time), preserving row-major grid order
    m5 = mov_c.reshape(tot_frames, 10, CELL_SIZE, 10, CELL_SIZE)
    grid_sums = m5.sum(axis=(2, 4)).transpose(1, 2, 0).reshape(100, tot_frames)
    m5 = mov2_c.reshape(tot_frames, 10, CELL_SIZE, 10, CELL_SIZE)
    grid2_sums = m5.sum(axis=(2, 4)).transpose(1, 2, 0).reshape(100, tot_frames)

    # compute dFF
    print('computing dFF traces...')
    # grids_dFF = ipf.calculate_dFF(grid_sums, sigma=300, t_axis=1,
    #                               GPU_AVAILABLE=GPU_AVAILABLE)
    # grids2_dFF = ipf.calculate_dFF(grid2_sums, sigma=300, t_axis=1,
    #                                GPU_AVAILABLE=GPU_AVAILABLE)
    grids_dFF = grid_sums
    grids2_dFF = grid2_sums

    # behaviour
    print('processing .txt file...')
    txt = bf.process_behavioural_data_imaging(txtpath)

    run_onsets = txt['run_onset_frames']

    # since we are stimming at the reward, the actual trial to look at is the next trial
    stim_conds = [t[15] for t in txt['trial_statements']]
    stim_idx = [trial for trial, cond in enumerate(stim_conds)
                if cond!='0']
    stim_idx_actual = [trial + 1 for trial in stim_idx]

    frame_times = txt['frame_times']
    pulse_times = txt['pulse_times']
    pulse_parameters = [pp for pp in txt['pulse_descriptions'] if pp]

    # extract pulse parameters
    pulse_width = float(pulse_parameters[-1][3])/1000000  # in s
    pulse_times = np.array(pulse_times)
    diffs = np.diff(pulse_times)
    split_idx = np.where(diffs >= 1000)[0] + 1
    pulse_trains = np.split(pulse_times, split_idx)

    pulse_frames = [
        [ipf.find_nearest(p, frame_times) for p in train]
        for train in pulse_trains
    ]
    # checks $FM against tot_frame
    if tot_frames < len(frame_times) - 3 or tot_frames > len(frame_times):
        raise ValueError(f'check $FM; frame count mismatch for {recname}')

    # block out opto artefact periods
    pulse_period_frames = np.concatenate(
        [np.arange(
            pulse_train[0]-3,
            pulse_train[-1]+int(pulse_width)*SAMP_FREQ+40
            )
        for pulse_train in pulse_frames]
        )
    grids_dFF[:, pulse_period_frames] = np.nan
    grids2_dFF[:, pulse_period_frames] = np.nan

    run_onsets_ctrl = [f for trial, f in enumerate(run_onsets)
                      if not np.isnan(f)
                      and trial not in stim_idx and trial not in stim_idx_actual
                      and f > 3 * SAMP_FREQ
                      and f < tot_frames - 5 * SAMP_FREQ]
    run_onsets_stim = [f for trial, f in enumerate(run_onsets)
                      if not np.isnan(f)
                      and trial in stim_idx_actual
                      and f > 3 * SAMP_FREQ
                      and f < tot_frames - 5 * SAMP_FREQ]

    # alignmeet
    aligned_ctrl = np.zeros((grids_dFF.shape[0],
                             len(run_onsets_ctrl),
                             8 * SAMP_FREQ))
    aligned_stim = np.zeros((grids_dFF.shape[0],
                             len(run_onsets_stim),
                             8 * SAMP_FREQ))
    for i, f in enumerate(run_onsets_ctrl):
        aligned_ctrl[:, i, :] = grids_dFF[:, f - 3 * SAMP_FREQ : f + 5 * SAMP_FREQ]
    for i, f in enumerate(run_onsets_stim):
        aligned_stim[:, i, :] = grids_dFF[:, f - 3 * SAMP_FREQ : f + 5 * SAMP_FREQ]

    # plotting
    # mean across trials (ignoring NaNs)
    means_ctrl = np.nanmean(aligned_ctrl, axis=1)  # (100, T)
    means_stim = np.nanmean(aligned_stim, axis=1)  # (100, T)

    # global y-lims for fair comparison across grids
    y_min = np.nanmin([means_ctrl, means_stim])
    y_max = np.nanmax([means_ctrl, means_stim])
    pad = 0.05 * (y_max - y_min)
    y_lims = (y_min - pad, y_max + pad)

    fig, axes = plt.subplots(10, 10, figsize=(20, 20), sharex=True, sharey=True)
    fig.suptitle(f'{recname} - per-grid mean dF/F (ctrl vs stim), run-aligned', y=0.92)

    for g in range(100):
        r, c = divmod(g, 10)
        ax = axes[r, c]

        ax.plot(TAXIS, means_ctrl[g], lw=1, label='ctrl', color='grey')
        ax.plot(TAXIS, means_stim[g], lw=1, label='stim', color='royalblue')

        # vertical line at run onset (t=0)
        ax.axvline(0, ls='--', lw=0.8, alpha=0.6, color='r')

        # light annotation: row/col index + trial counts
        n_ctrl = np.sum(np.any(np.isfinite(aligned_ctrl[g]), axis=1))
        n_stim = np.sum(np.any(np.isfinite(aligned_stim[g]), axis=1))
        ax.text(0.02, 0.95, f'({r},{c})\nC:{n_ctrl} S:{n_stim}',
                transform=ax.transAxes, va='top', ha='left', fontsize=7)

        # axes cosmetics
        ax.set_xlim(TAXIS[0], TAXIS[-1])
        ax.set_ylim(*y_lims)

        if g == 0:
            ax.legend(frameon=False, fontsize=8, loc='lower right')

    # shared labels on outer figure
    fig.text(0.5, 0.06, 'time from run onset (s)', ha='center')
    fig.text(0.04, 0.5, 'dF/F', va='center', rotation='vertical')

    plt.tight_layout(rect=[0.04, 0.06, 0.995, 0.92])
    out_png = os.path.join(savepath, f'{recname}_grid_means_ctrl_stim_run.png')
    fig.savefig(out_png, dpi=300)
    plt.close(fig)
    print_files_saved([
        ('figure', out_png),
    ])
