# -*- coding: utf-8 -*-
'''
Created on Thu Aug  8 18:20:50 2024

pick out ROIs that have significant activity
- criteria:
    1. exceeds 99% shuffled and
    2. does not correlate with neuropil
    in the 2 seconds following the aligned landmark

@author: Dinghao Luo

- kept generated outputs on project_paths stems
'''

#%% imports
import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import pearsonr, sem


# import util functions
repo_root = Path(__file__).resolve().parent
while not (repo_root / 'utils').exists() and repo_root != repo_root.parent:
    repo_root = repo_root.parent
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))

import imaging_utility_functions as iuf
import imaging_pipeline_functions as ipf
from common_functions import smooth_convolve
import project_paths as pp

import rec_list
pathGRABNE = rec_list.pathHPCGRABNE


#%% output paths
PLOT_PAYLOAD_STEM = pp.GRABNE_STEM / 'single_session_avgs' / 'processed_data'
RUN_FIG_STEM = pp.GRABNE_FIGURES_STEM / 'single_session_avgs' / 'RO_aligned_roi_sig_only'
REW_FIG_STEM = pp.GRABNE_FIGURES_STEM / 'single_session_avgs' / 'rew_aligned_roi_sig_only'


#%% statistic parameters
alpha = .01
sig_window = 2  # in seconds


#%% command line
def build_arg_parser():
    parser = argparse.ArgumentParser(
        description='extract significant GRABNE ROIs from aligned ROI outputs.'
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
        help='also regenerate significant-ROI summary figures during extraction; default is data only',
        )
    parser.add_argument(
        '--no-plots',
        dest='make_plots',
        action='store_false',
        help=argparse.SUPPRESS,
        )
    parser.set_defaults(make_plots=False)
    return parser


def get_recording_paths(recording_filter=None):
    if not recording_filter:
        return pathGRABNE

    filtered_paths = [
        rec_path for rec_path in pathGRABNE
        if recording_filter in rec_path or recording_filter in rec_path[-17:]
        ]
    if not filtered_paths:
        raise ValueError(f'no recordings matched filter: {recording_filter}')
    return filtered_paths


#%% plotting
def save_sig_roi_plots(recname, payload):
    signif_act_roi_run = payload['signif_act_roi_run']
    signif_act_roi_rew = payload['signif_act_roi_rew']
    sig_means_run = payload['sig_means_run']
    sig_sems_run = payload['sig_sems_run']
    sig_shuf_run = payload['sig_shuf_run']
    sig_shuf_95_run = payload['sig_shuf_95_run']
    sig_shuf_5_run = payload['sig_shuf_5_run']
    sig_means_rew = payload['sig_means_rew']
    sig_sems_rew = payload['sig_sems_rew']
    sig_shuf_rew = payload['sig_shuf_rew']
    sig_shuf_95_rew = payload['sig_shuf_95_rew']
    sig_shuf_5_rew = payload['sig_shuf_5_rew']

    RUN_FIG_STEM.mkdir(parents=True, exist_ok=True)
    REW_FIG_STEM.mkdir(parents=True, exist_ok=True)

    tot_sig_roi_run = len(signif_act_roi_run)
    tot_sig_roi_rew = len(signif_act_roi_rew)
    xaxis = (np.arange(5*30)-30)/30

    n_row = int(tot_sig_roi_run/5)
    if n_row == 0:
        n_col = tot_sig_roi_rew
        n_row = 1
    else:
        n_col = int(np.ceil(tot_sig_roi_run/n_row))
    fig = plt.figure(1, figsize=(n_col*3, n_row*2.1))
    for p, r in enumerate(signif_act_roi_run):
        ax = fig.add_subplot(n_row, n_col, p+1)
        ax.set(xlim=(-1,4), xlabel='time (s)', xticks=[0,2,4],
               ylabel='dFF', title='roi {}'.format(r))
        ax.plot(xaxis, sig_shuf_run[p], color='grey', linewidth=.8)
        ax.fill_between(xaxis, sig_shuf_95_run[p],
                               sig_shuf_5_run[p],
                        color='grey', edgecolor='none', alpha=.2)
        ax.plot(xaxis, sig_means_run[p], color='darkgreen', linewidth=.8)
        ax.fill_between(xaxis, sig_means_run[p]+sig_sems_run[p],
                               sig_means_run[p]-sig_sems_run[p],
                        color='darkgreen', edgecolor='none', alpha=.2)
        ax.axvspan(0, 0, color='burlywood', alpha=.5,
                   linestyle='dashed', linewidth=1)
    fig.suptitle('run_aligned_sig_act_only')
    fig.tight_layout()
    fig.savefig(RUN_FIG_STEM / f'{recname}.png',
                dpi=120,
                bbox_inches='tight')
    plt.close(fig)

    n_row = int(tot_sig_roi_rew/5)
    if n_row == 0:
        n_col = tot_sig_roi_rew
        n_row = 1
    else:
        n_col = int(np.ceil(tot_sig_roi_rew/n_row))
    fig = plt.figure(1, figsize=(n_col*3, n_row*2.1))
    for p, r in enumerate(signif_act_roi_rew):
        ax = fig.add_subplot(n_row, n_col, p+1)
        ax.set(xlim=(-1,4), xlabel='time (s)', xticks=[0,2,4],
               ylabel='dFF', title='roi {}'.format(r))
        ax.plot(xaxis, sig_shuf_rew[p], color='grey', linewidth=.8)
        ax.fill_between(xaxis, sig_shuf_95_rew[p],
                               sig_shuf_5_rew[p],
                        color='grey', edgecolor='none', alpha=.2)
        ax.plot(xaxis, sig_means_rew[p], color='darkgreen', linewidth=.8)
        ax.fill_between(xaxis, sig_means_rew[p]+sig_sems_rew[p],
                               sig_means_rew[p]-sig_sems_rew[p],
                        color='darkgreen', edgecolor='none', alpha=.2)
        ax.axvspan(0, 0, color='burlywood', alpha=.5,
                   linestyle='dashed', linewidth=1)
    fig.suptitle('rew_aligned_sig_act_only')
    fig.tight_layout()
    fig.savefig(REW_FIG_STEM / f'{recname}.png',
                dpi=120,
                bbox_inches='tight')
    plt.close(fig)


#%% extraction
def process_recording(rec_path, processed_recs, make_plots=False):
    recname = rec_path[-17:]

    # if processed, keep the original printout but still recompute
    if recname in processed_recs:
        print(recname+' already processed; skipped')
    else:
        print(recname)

    ext_path = rec_path+r'_roi_extract'

    run_aligned = np.load(ext_path+r'\suite2pROI_run_dFF_aligned.npy', allow_pickle=True)
    run_aligned_neu = np.load(ext_path+r'\suite2pROI_run_dFF_aligned_neu.npy', allow_pickle=True)

    rew_aligned = np.load(ext_path+r'\suite2pROI_rew_dFF_aligned.npy', allow_pickle=True)
    rew_aligned_neu = np.load(ext_path+r'\suite2pROI_rew_dFF_aligned_neu.npy', allow_pickle=True)

    tot_roi, tot_trial_run, _ = run_aligned.shape
    tot_trial_rew = rew_aligned.shape[1]

    signif_act_roi_run = []
    signif_act_roi_rew = []

    sig_means_run = []
    sig_sems_run = []
    sig_shuf_run = []
    sig_shuf_95_run = []
    sig_shuf_5_run = []
    sig_means_rew = []
    sig_sems_rew = []
    sig_shuf_rew = []
    sig_shuf_95_rew = []
    sig_shuf_5_rew = []

    for r in range(tot_roi):
        curr_run_aligned = run_aligned[r,:,:]
        curr_run_aligned_neu = run_aligned_neu[r,:,:]
        curr_rew_aligned = rew_aligned[r,:,:]
        curr_rew_aligned_neu = rew_aligned_neu[r,:,:]

        # check for nan
        if ipf.sum_mat(np.isnan(curr_run_aligned)) != 0 or ipf.sum_mat(np.isnan(curr_rew_aligned)) != 0:
            continue

        # smoothing?
        for t in range(tot_trial_run):
            curr_run_aligned[t,:] = smooth_convolve(curr_run_aligned[t,:])
            curr_run_aligned_neu[t,:] = smooth_convolve(curr_run_aligned_neu[t,:])
        for t in range(tot_trial_rew):
            curr_rew_aligned[t,:] = smooth_convolve(curr_rew_aligned[t,:])
            curr_rew_aligned_neu[t,:] = smooth_convolve(curr_rew_aligned_neu[t,:])

        curr_run_mean = np.mean(curr_run_aligned, axis=0)
        curr_run_mean_neu = np.mean(curr_run_aligned_neu, axis=0)
        curr_rew_mean = np.mean(curr_rew_aligned, axis=0)
        curr_rew_mean_neu = np.mean(curr_rew_aligned_neu, axis=0)

        shuf_run, shuf_95_run, shuf_5_run = iuf.circ_shuffle(curr_run_aligned, alpha=alpha, num_shuf=500)
        shuf_rew, shuf_95_rew, shuf_5_rew = iuf.circ_shuffle(curr_rew_aligned, alpha=alpha, num_shuf=500)

        # test for data-shuff significance by looking at the 2 seconds after RO
        sig_95up_run = [f for f, [v, vs] in enumerate(zip(curr_run_mean[:30*sig_window], shuf_95_run[:30*sig_window])) if v>vs]
        sig_5down_run = [f for f, [v, vs] in enumerate(zip(curr_run_mean[:30*sig_window], shuf_5_run[:30*sig_window])) if v<vs]
        sig_95up_rew = [f for f, [v, vs] in enumerate(zip(curr_rew_mean[:30*sig_window], shuf_95_rew[:30*sig_window])) if v>vs]
        sig_5down_rew = [f for f, [v, vs] in enumerate(zip(curr_rew_mean[:30*sig_window], shuf_5_rew[:30*sig_window])) if v<vs]

        # test for correlation
        _, pval_run = pearsonr(curr_run_mean, curr_run_mean_neu)
        _, pval_rew = pearsonr(curr_rew_mean, curr_rew_mean_neu)

        if pval_run >= .05 and (len(sig_95up_run)>30 or len(sig_5down_run)>30):
            signif_act_roi_run.append(r)
        if pval_rew >= .05 and (len(sig_95up_rew)>30 or len(sig_5down_rew)>30):
            signif_act_roi_rew.append(r)

        sig_means_run.append(curr_run_mean)
        sig_sems_run.append(sem(curr_run_aligned, axis=0))
        sig_shuf_run.append(shuf_run)
        sig_shuf_95_run.append(shuf_95_run)
        sig_shuf_5_run.append(shuf_5_run)
        sig_means_rew.append(curr_rew_mean)
        sig_sems_rew.append(sem(curr_rew_aligned, axis=0))
        sig_shuf_rew.append(shuf_rew)
        sig_shuf_95_rew.append(shuf_95_rew)
        sig_shuf_5_rew.append(shuf_5_rew)

    plot_payload = {
        'signif_act_roi_run': signif_act_roi_run,
        'signif_act_roi_rew': signif_act_roi_rew,
        'sig_means_run': sig_means_run,
        'sig_sems_run': sig_sems_run,
        'sig_shuf_run': sig_shuf_run,
        'sig_shuf_95_run': sig_shuf_95_run,
        'sig_shuf_5_run': sig_shuf_5_run,
        'sig_means_rew': sig_means_rew,
        'sig_sems_rew': sig_sems_rew,
        'sig_shuf_rew': sig_shuf_rew,
        'sig_shuf_95_rew': sig_shuf_95_rew,
        'sig_shuf_5_rew': sig_shuf_5_rew,
        }
    PLOT_PAYLOAD_STEM.mkdir(parents=True, exist_ok=True)
    np.save(PLOT_PAYLOAD_STEM / f'{recname}_sig_roi_plot_payload.npy',
            plot_payload)

    if make_plots:
        save_sig_roi_plots(recname, plot_payload)

    return recname, signif_act_roi_run, signif_act_roi_rew


#%% main
def main(argv=None):
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    # read dataframe and get indices
    existing_df = pd.read_pickle(pp.GRABNE_STEM / 'significant_activity_roi.pkl')
    processed_recs = list(existing_df.index.values)

    # dataframe to contain all results
    profiles = {'sig_act_roi_run': [],
                'sig_act_roi_rew': []
                }
    df = pd.DataFrame(profiles, dtype=object)

    for rec_path in get_recording_paths(args.recording_filter):
        recname, signif_act_roi_run, signif_act_roi_rew = process_recording(
            rec_path,
            processed_recs,
            make_plots=args.make_plots
            )
        df.loc[recname] = pd.Series(
            [signif_act_roi_run, signif_act_roi_rew],
            index=['sig_act_roi_run', 'sig_act_roi_rew']
            )

    # save dataframe
    df.to_csv(pp.GRABNE_STEM / 'significant_activity_roi.csv')
    df.to_pickle(pp.GRABNE_STEM / 'significant_activity_roi.pkl')


if __name__ == '__main__':
    main()
