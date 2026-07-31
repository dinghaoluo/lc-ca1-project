# -*- coding: utf-8 -*-
'''
Created on Fri Nov  7 14:25:26 2025

extract GRABNE signals locked to tone activation

@author: Dinghao Luo
'''

#%% imports
import argparse
import gc
from pathlib import Path
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

repo_root = Path(__file__).resolve().parents[3]
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))

import imaging_pipeline_functions as ipf
from common_functions import mpl_formatting
from console_formatting import print_session
import project_paths as pp
mpl_formatting()

import rec_list
paths = rec_list.pathGRABNETone + rec_list.pathGRABNEToneDbhBlock


#%% path stems
BEF = 10
AFT = 30
SAMP_FREQ = 30  # Hz

all_sess_stem = pp.GRABNE_STEM / 'all_sessions_tone'
all_sess_fig_stem = pp.GRABNE_FIGURES_STEM / 'all_sessions_tone'
mice_exp_stem = pp.MICEEXP_ROOT


#%% motion correction and plotting
def regress_out_motion(ch, ref, fit_mask=None, allow_lag=True, max_lag=5):
    '''
    regress out motion carried by ref from ch using linear regression with optional small lag search.

    parameters:
    - ch: 1d ndarray, target signal (green dff)
    - ref: 1d ndarray, motion reference (red dff)
    - fit_mask: 1d boolean ndarray marking frames used to fit beta (True = use). if None, use all
    - allow_lag: bool, whether to search integer lag between channels
    - max_lag: int, maximum absolute lag in frames to try

    returns:
    - resid: 1d ndarray, ch with fitted motion component removed (residual)
    - beta: tuple of (intercept, slope) from best fit
    - best_lag: int, lag (in frames) applied to ref (positive means ref shifted forward)
    - r2: float, variance explained by motion on the fit_mask
    '''
    ch = np.asarray(ch).astype(np.float64)
    ref = np.asarray(ref).astype(np.float64)
    n = ch.size

    if fit_mask is None:
        fit_mask = np.ones(n, dtype=bool)
    else:
        fit_mask = fit_mask.astype(bool)

    def fit_for_lag(lag):
        # shift ref by lag (np.roll), then zero out edges that wrap
        ref_shift = np.roll(ref, lag).copy()
        if lag > 0:
            ref_shift[:lag] = ref_shift[lag]  # avoid wrap artefact
        elif lag < 0:
            ref_shift[lag:] = ref_shift[lag-1]
        X = np.vstack([np.ones(n), ref_shift]).T
        m = fit_mask.copy()
        # exclude padded edge points from the fit
        if lag > 0:
            m[:lag] = False
        elif lag < 0:
            m[lag:] = False
        if m.sum() < 5:
            return np.inf, (0.0, 0.0), None, -np.inf
        beta = np.linalg.lstsq(X[m], ch[m], rcond=None)[0]  # [intercept, slope]
        pred = X @ beta
        resid = ch - pred
        # variance explained on the mask
        r2 = 1.0 - np.var(resid[m]) / np.var(ch[m])
        mse = np.mean((ch[m] - pred[m])**2)
        return mse, (float(beta[0]), float(beta[1])), resid, r2

    best = (np.inf, (0.0, 0.0), None, -np.inf, 0)
    lags = [0]
    if allow_lag and max_lag > 0:
        lags = list(range(-max_lag, max_lag + 1))

    for lag in lags:
        mse, beta, resid, r2 = fit_for_lag(lag)
        if mse < best[0]:
            best = (mse, beta, resid, r2, lag)

    _, beta, resid, r2, best_lag = best
    return resid, beta, best_lag, float(r2)

def save_tone_figures(savepath, recname, payload):
    savepath.mkdir(parents=True, exist_ok=True)

    dFF = payload['dFF']
    dFF2 = payload['dFF2']
    dFF_corr = payload['dFF_corr']
    buzzer_frames = payload['buzzer_frames']
    r2 = payload['r2']
    tone_aligned = payload['tone_aligned']
    tone_aligned2 = payload['tone_aligned2']
    tone_aligned_c = payload['tone_aligned_c']

    fig, ax = plt.subplots(figsize=(10, 3))
    ax.plot(dFF, color='darkgreen', alpha=.6, label='green dFF')
    ax.plot(dFF2, color='darkred', alpha=.3, label='red dFF (motion)')
    ax.plot(dFF_corr, color='k', lw=1.0, label='green corrected')
    for buzz in buzzer_frames:
        ax.axvspan(buzz-1, buzz+1, alpha=.1)
    ax.legend(loc='upper right')
    ax.set_title(f'{recname}  (R² motion explained = {r2:.3f})')
    fig.tight_layout()
    fig.savefig(savepath / f'{recname}_qc_motion_correction.png', dpi=200)
    plt.close(fig)

    # plot aligned means
    fig, ax = plt.subplots(figsize=(3, 2))
    ax.plot(tone_aligned.T, color='darkgreen', alpha=.02)
    ax.plot(np.mean(tone_aligned, axis=0), color='darkgreen', lw=1.2, label='green mean')
    ax.plot(np.mean(tone_aligned2, axis=0), color='darkred', alpha=.5, label='red mean')
    ax.plot(np.mean(tone_aligned_c, axis=0), color='k', lw=1.2, label='green corrected mean')
    ax.axvline(BEF * SAMP_FREQ, ls='--', lw=.8, color='k')
    ax.legend(frameon=False, fontsize=7)
    fig.tight_layout()
    fig.savefig(savepath / f'{recname}_tone_aligned_means.png', dpi=200)
    plt.close(fig)


#%% main
def process_session(path, make_plots=False):
    recname = Path(path).name
    print_session(recname)

    plane_stem = pp.resolve_suite2p_session_stem(path) / 'plane0'
    sessname = recname.replace('i', '')

    binpath = plane_stem / 'data.bin'
    bin2path = plane_stem / 'data_chan2.bin'
    opspath = plane_stem / 'ops.npy'
    txtpath = mice_exp_stem / f'ANMD{recname[1:4]}' / f'{sessname}T.txt'

    savepath = all_sess_stem / f'{recname}'
    savepath.mkdir(parents=True, exist_ok=True)
    figurepath = all_sess_fig_stem / f'{recname}'
    if make_plots:
        figurepath.mkdir(parents=True, exist_ok=True)
    proc_data_path = savepath / 'processed_data'
    proc_data_path.mkdir(parents=True, exist_ok=True)

    # load data
    ops = np.load(opspath, allow_pickle=True).item()
    tot_frames = ops['nframes']
    shape = tot_frames, ops['Ly'], ops['Lx']

    if make_plots:
        print('loading movies, saving references, and generating reference figures...')
    else:
        print('loading movies and saving reference arrays...')
    mov  = np.memmap(binpath, mode='r', dtype='int16', shape=shape).astype(np.float32)
    mov2 = np.memmap(bin2path, mode='r', dtype='int16', shape=shape).astype(np.float32)

    tot_frames = mov.shape[0]  # once loaded, update tot_frames to be the max frame number

    ipf.plot_reference(
        mov, recname=recname, outpath=savepath, channel=1,
        save_figure=make_plots, figure_outpath=figurepath
        )
    ipf.plot_reference(
        mov2, recname=recname, outpath=savepath, channel=2,
        save_figure=make_plots, figure_outpath=figurepath
        )

    raw_trace  = np.sum(mov, axis=(1,2))
    raw_trace2 = np.sum(mov2, axis=(1,2))

    dFF  = ipf.calculate_dFF(raw_trace, t_axis=0)
    dFF2 = ipf.calculate_dFF(raw_trace2, t_axis=0)

    # get tone stamps
    txt = ipf.process_txt_nobeh(txtpath)
    frame_times   = txt['frame_times']
    buzzer_times  = txt['buzzer_times']
    buzzer_frames = [ipf.find_nearest(t, frame_times) for t in buzzer_times]

    # build fit mask: use all frames except tone-locked windows when fitting motion regression
    fit_mask = np.ones(tot_frames, dtype=bool)
    for buzz in buzzer_frames:
        lo = max(0, int(buzz - BEF * SAMP_FREQ))
        hi = min(tot_frames, int(buzz + AFT * SAMP_FREQ))
        fit_mask[lo:hi] = False

    # optionally also exclude frames with large suite2p shifts (if available)
    for key in ('xoff', 'yoff'):
        if key in ops:
            shifts = np.asarray(ops[key])
            thr = np.percentile(np.abs(shifts), 99.5)
            fit_mask[np.abs(shifts) > thr] = False

    # regress out motion from green using red as reference (search up to ±5-frame lag)
    dFF_corr, beta, lag, r2 = regress_out_motion(
        dFF, dFF2, fit_mask=fit_mask, allow_lag=True, max_lag=5
        )
    print(f'motion correction: beta0={beta[0]:.4f}, beta1={beta[1]:.4f}, lag={lag} frames, R²={r2:.3f}')

    # align around tones
    tot_tones = len(buzzer_frames)
    win = (BEF + AFT) * SAMP_FREQ
    tone_aligned   = np.zeros((tot_tones, win))
    tone_aligned2  = np.zeros((tot_tones, win))
    tone_aligned_c = np.zeros((tot_tones, win))
    for i, buzz in enumerate(buzzer_frames):
        lo = int(buzz - BEF * SAMP_FREQ)
        hi = int(buzz + AFT * SAMP_FREQ)
        if lo >= 0 and hi < tot_frames:
            tone_aligned[i, :]   = dFF[lo:hi]
            tone_aligned2[i, :]  = dFF2[lo:hi]
            tone_aligned_c[i, :] = dFF_corr[lo:hi]

    np.save(proc_data_path / f'{recname}_tone_aligned.npy', tone_aligned)
    np.save(proc_data_path / f'{recname}_tone_aligned_ch2.npy', tone_aligned2)
    np.save(proc_data_path / f'{recname}_tone_aligned_corrected.npy',
            tone_aligned_c)

    plot_payload = {
        'dFF': dFF,
        'dFF2': dFF2,
        'dFF_corr': dFF_corr,
        'buzzer_frames': np.asarray(buzzer_frames),
        'buzzer_times': np.asarray(buzzer_times),
        'motion_beta': beta,
        'motion_lag': lag,
        'r2': r2,
        'tone_aligned': tone_aligned,
        'tone_aligned2': tone_aligned2,
        'tone_aligned_c': tone_aligned_c,
        }
    np.save(proc_data_path / f'{recname}_tone_plot_payload.npy',
            plot_payload)

    if make_plots:
        save_tone_figures(figurepath, recname, plot_payload)


#%% main
def main(argv=None):
    parser = argparse.ArgumentParser(
        description='extract GRABNE tone-locked imaging data.'
    )
    parser.add_argument(
        '--recording',
        dest='recording_filter',
        help='optional substring filter for recname or recording path',
    )
    parser.add_argument(
        '--plots',
        dest='make_plots',
        action='store_true',
        help='also regenerate figure outputs during extraction; default is data only',
    )
    parser.add_argument(
        '--no-plots',
        dest='make_plots',
        action='store_false',
        help=argparse.SUPPRESS,
    )
    parser.set_defaults(make_plots=False)
    args = parser.parse_args(argv)

    session_paths = paths
    if args.recording_filter:
        session_paths = [
            path for path in paths
            if args.recording_filter in path
            or args.recording_filter in Path(path).name
        ]
        if not session_paths:
            raise ValueError(f'no recordings matched filter: {args.recording_filter}')

    for path in session_paths:
        process_session(path, make_plots=args.make_plots)
        print(f'clearing memory after {Path(path).name}...', flush=True)
        plt.close('all')
        gc.collect()

if __name__ == '__main__':
    main()
