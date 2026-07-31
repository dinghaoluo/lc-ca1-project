# -*- coding: utf-8 -*-
'''
Created on Thu Jun 25 2026

test first-lick early/late split thresholds for HPC pyramidal-cell profiles.

@author: Dinghao Luo
'''

#%% imports
import argparse
from pathlib import Path
import pickle
import sys

repo_root = Path(__file__).resolve().parents[3]
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.io as sio
from scipy.stats import sem, ttest_ind, ttest_rel, wilcoxon

from common_functions import mpl_formatting
from console_formatting import print_files_saved, print_session, print_statistics_section
from plotting_functions import plot_violin_with_scatter
import first_lick_analysis_functions as flaf
import project_paths as pp
mpl_formatting()

import rec_list
bad_behs = rec_list.pathHPCbadbeh
paths = rec_list.pathHPCLCopt + rec_list.pathHPCLCtermopt + rec_list.pathHPC_Raphi
paths = [p for p in paths if p not in bad_behs]
recnames = [Path(path).name for path in paths]


#%% parameters
SAMP_FREQ = 1250
RUN_ONSET_BIN = 3750
BEF = 1
AFT = 4

DEFAULT_THRESHOLDS = (2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8)
LATE_MAX = 3.5
MIN_MATCHED = 10
MATCH_K = 1.5

BIN_SIZE_MS = 500
TOTAL_LEN_MS = 3500
N_BINS = TOTAL_LEN_MS // BIN_SIZE_MS
ACC_WIN_S = 1
SUMMARY_WINDOW = (-.5, .5)

XAXIS = np.arange((BEF + AFT) * SAMP_FREQ) / SAMP_FREQ - BEF
X_SEC_PLOT = np.arange(4000) / 1000.0
YLIM_SPEED = (0, 70)

PYRUP_CLASS = 'run-onset ON'
PYRDOWN_CLASS = 'run-onset OFF'
EARLY_PYRUP = 'early PyrUp'
LATE_PYRUP = 'late PyrUp'
EARLY_PYRDOWN = 'early PyrDown'
LATE_PYRDOWN = 'late PyrDown'

early_c = (168 / 255, 155 / 255, 202 / 255)
late_c  = (102 / 255, 83 / 255, 162 / 255)

ORIGINAL_HPC_STEM = Path('Z:/Dinghao/code_dinghao/HPC_ephys')
ORIGINAL_BEHAVIOUR_ROOT = Path(
    'Z:/Dinghao/code_dinghao/behaviour/all_experiments'
)


#%% command line
parser = argparse.ArgumentParser(
    description='test alternative first-lick split thresholds for HPC pyramidal-cell profiles.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)
parser.add_argument('--thresholds', nargs='+', type=float,
                    default=list(DEFAULT_THRESHOLDS),
                    help='candidate split thresholds in seconds')
parser.add_argument('--late-max', type=float, default=LATE_MAX,
                    help='upper bound for late first-lick trials, in seconds')
parser.add_argument('--min-matched', type=int, default=MIN_MATCHED,
                    help='minimum early and late speed-matched trials required per session')
parser.add_argument('--match-k', type=float, default=MATCH_K,
                    help='symmetric SD multiplier for 7-bin speed matching')
parser.add_argument('--hpc-stem', type=Path, default=ORIGINAL_HPC_STEM,
                    help='folder containing HPC profile pickles and per-session train arrays')
parser.add_argument('--behaviour-root', type=Path, default=ORIGINAL_BEHAVIOUR_ROOT,
                    help='folder containing HPCLC, HPCLCterm, and HPCRaphi behaviour pickles')
parser.add_argument('--repo-inputs', action='store_true',
                    help='use generated in-repo profile, train, and behaviour files')
parser.add_argument('--recnames', nargs='+', default=None,
                    help='recording names to include')
args = parser.parse_args()
thresholds = sorted(set(args.thresholds))
late_max = args.late_max
min_matched = args.min_matched
match_k = args.match_k

if args.recnames is not None:
    requested_recnames = set(args.recnames)
    recnames = [recname for recname in recnames if recname in requested_recnames]

if args.repo_inputs:
    HPC_stem = pp.HPC_EPHYS_STEM
    all_exp_stem = pp.behaviour_experiment_stem('HPCLC')
    all_exp_term_stem = pp.behaviour_experiment_stem('HPCLCterm')
    all_exp_raphi_stem = pp.behaviour_experiment_stem('HPCRaphi')
else:
    HPC_stem = args.hpc_stem
    all_exp_stem = args.behaviour_root / 'HPCLC'
    all_exp_term_stem = args.behaviour_root / 'HPCLCterm'
    all_exp_raphi_stem = args.behaviour_root / 'HPCRaphi'

first_lick_stem = (
    pp.HPC_EPHYS_FIGURES_STEM
    / 'first_lick_analysis_split_threshold_diagnostics'
)
first_lick_stem.mkdir(parents=True, exist_ok=True)


#%% speed matching and diagnostic plots
def _speed_match_trials(early_trials, late_trials, speed_times):
    E_bins, e_valid = flaf.compute_binned_speed_matrix(
        early_trials,
        speed_times,
        n_bins=N_BINS,
        bin_size=BIN_SIZE_MS,
    )
    L_bins, l_valid = flaf.compute_binned_speed_matrix(
        late_trials,
        speed_times,
        n_bins=N_BINS,
        bin_size=BIN_SIZE_MS,
    )

    if not len(E_bins) or not len(L_bins):
        return [], []

    e_mu = E_bins.mean(axis=0)
    e_sd = E_bins.std(axis=0, ddof=0)
    l_mu = L_bins.mean(axis=0)
    l_sd = L_bins.std(axis=0, ddof=0)

    e_low, e_high = e_mu - match_k * e_sd, e_mu + match_k * e_sd
    l_low, l_high = l_mu - match_k * l_sd, l_mu + match_k * l_sd

    late_mask = np.all((L_bins >= e_low) & (L_bins <= e_high), axis=1)
    early_mask = np.all((E_bins >= l_low) & (E_bins <= l_high), axis=1)

    matched_late = [l_valid[i] for i in np.where(late_mask)[0]]
    matched_early = [e_valid[i] for i in np.where(early_mask)[0]]
    return matched_early, matched_late

def _speed_scalar_and_accel(trial_list, speed_times):
    speed_scalars = []
    accel_scalars = []
    acc_n = int(ACC_WIN_S * SAMP_FREQ)

    for trial in trial_list:
        if trial >= len(speed_times):
            continue
        speeds = np.asarray([pt[1] for pt in speed_times[trial]], dtype=float)
        if len(speeds) < acc_n + 1 or not np.any(np.isfinite(speeds)):
            continue
        speed_scalars.append(np.nanmean(speeds))
        accel_scalars.append(np.nanmean(np.diff(speeds[:acc_n + 1]) * SAMP_FREQ))

    return np.asarray(speed_scalars, dtype=float), np.asarray(accel_scalars, dtype=float)

def _plot_speed_pre_matched(result, threshold, stem):
    threshold_text = f'{threshold:.2f}'.rstrip('0').rstrip('.')
    if not result['sess_early_speed_raw'] or not result['sess_late_speed_raw']:
        print(f'pre-matching speed skipped for t={threshold_text}: no traces')
        return

    E = np.vstack(result['sess_early_speed_raw'])
    L = np.vstack(result['sess_late_speed_raw'])
    E_mean = np.mean(E, axis=0)
    E_sem = sem(E, axis=0)
    L_mean = np.mean(L, axis=0)
    L_sem = sem(L, axis=0)

    fig, ax = plt.subplots(figsize=(2.1, 2.0))
    ax.plot(X_SEC_PLOT, E_mean, c=early_c, label=f'early (<{threshold_text} s)')
    ax.fill_between(
        X_SEC_PLOT,
        E_mean + E_sem,
        E_mean - E_sem,
        color=early_c,
        edgecolor='none',
        alpha=.25,
    )
    ax.plot(X_SEC_PLOT, L_mean, c=late_c, label=f'late ({threshold_text}-{late_max:g} s)')
    ax.fill_between(
        X_SEC_PLOT,
        L_mean + L_sem,
        L_mean - L_sem,
        color=late_c,
        edgecolor='none',
        alpha=.25,
    )

    ax.set(
        xlabel='Time from run onset (s)',
        ylabel='Speed (cm/s)',
        xlim=(0, 4),
        ylim=YLIM_SPEED,
        title='pre-matching speed',
    )
    ax.legend(frameon=False, fontsize=7)
    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)

    fig.tight_layout()
    for ext in ['.png', '.pdf']:
        fig.savefig(stem / f'speed_pre_matched{ext}', dpi=300, bbox_inches='tight')
    plt.close(fig)

def _plot_behaviour_raw_violins(result, threshold, stem):
    threshold_text = f'{threshold:.2f}'.rstrip('0').rstrip('.')
    raw_E = np.asarray(result['sess_early_speed_scalar_raw'], dtype=float)
    raw_L = np.asarray(result['sess_late_speed_scalar_raw'], dtype=float)
    valid_speed = np.isfinite(raw_E) & np.isfinite(raw_L)
    raw_E, raw_L = raw_E[valid_speed], raw_L[valid_speed]
    if len(raw_E) and len(raw_L):
        plot_violin_with_scatter(
            raw_E,
            raw_L,
            colour0=early_c,
            colour1=late_c,
            showscatter=False,
            xticklabels=['early', 'late'],
            title=f't = {threshold_text} s',
            stats_labels=[
                f'pre-matching speed, first lick <{threshold_text} s',
                f'pre-matching speed, first lick {threshold_text}-{late_max:g} s',
            ],
            print_statistics=True,
            save=True,
            savepath=stem / 'all_run_onset_speed_violin',
            show=False,
            close=True,
        )

    raw_ea = np.asarray(result['sess_early_accel_raw'], dtype=float)
    raw_la = np.asarray(result['sess_late_accel_raw'], dtype=float)
    valid_acceleration = np.isfinite(raw_ea) & np.isfinite(raw_la)
    raw_ea, raw_la = raw_ea[valid_acceleration], raw_la[valid_acceleration]
    if len(raw_ea) and len(raw_la):
        plot_violin_with_scatter(
            raw_ea,
            raw_la,
            colour0=early_c,
            colour1=late_c,
            showscatter=False,
            xticklabels=['early', 'late'],
            title=f't = {threshold_text} s',
            stats_labels=[
                f'pre-matching initial acceleration, first lick <{threshold_text} s',
                f'pre-matching initial acceleration, first lick {threshold_text}-{late_max:g} s',
            ],
            print_statistics=True,
            save=True,
            savepath=stem / 'all_run_onset_acc_boxplot',
            show=False,
            close=True,
        )

def _plot_speed_post_matched(result, threshold, stem):
    threshold_text = f'{threshold:.2f}'.rstrip('0').rstrip('.')
    if not result['sess_early_speed'] or not result['sess_late_speed']:
        print(f'post-matching speed skipped for t={threshold_text}: no traces')
        return

    E = np.vstack(result['sess_early_speed'])
    L = np.vstack(result['sess_late_speed'])
    E_mean = np.mean(E, axis=0)
    E_sem = sem(E, axis=0)
    L_mean = np.mean(L, axis=0)
    L_sem = sem(L, axis=0)

    fig, ax = plt.subplots(figsize=(2.1, 2.0))
    ax.plot(X_SEC_PLOT, E_mean, c=early_c, label=f'early (<{threshold_text} s)')
    ax.fill_between(
        X_SEC_PLOT,
        E_mean + E_sem,
        E_mean - E_sem,
        color=early_c,
        edgecolor='none',
        alpha=.25,
    )
    ax.plot(X_SEC_PLOT, L_mean, c=late_c, label=f'late ({threshold_text}-{late_max:g} s)')
    ax.fill_between(
        X_SEC_PLOT,
        L_mean + L_sem,
        L_mean - L_sem,
        color=late_c,
        edgecolor='none',
        alpha=.25,
    )

    E_bins = np.vstack([
        E[:, i * BIN_SIZE_MS:(i + 1) * BIN_SIZE_MS].mean(axis=1)
        for i in range(N_BINS)
    ]).T
    L_bins = np.vstack([
        L[:, i * BIN_SIZE_MS:(i + 1) * BIN_SIZE_MS].mean(axis=1)
        for i in range(N_BINS)
    ]).T

    pvals, tvals = [], []
    dz = []
    for i in range(N_BINS):
        t, p = ttest_ind(
            E_bins[:, i],
            L_bins[:, i],
            equal_var=False,
            nan_policy='omit',
        )
        tvals.append(t)
        pvals.append(p)
        d = E_bins[:, i] - L_bins[:, i]
        d = d[np.isfinite(d)]
        dz.append(np.mean(d) / np.std(d, ddof=1))
    pvals = np.asarray(pvals)

    ymax = max((E_mean + E_sem).max(), (L_mean + L_sem).max())
    ymin = min((E_mean - E_sem).min(), (L_mean - L_sem).min())
    yr = ymax - ymin if ymax > ymin else 1.0
    bar_y = ymax + 0.06 * yr
    text_y = ymax + 0.08 * yr

    for i in range(N_BINS):
        x_left = i * 0.5 + 0.1
        x_right = (i + 1) * 0.5 - 0.1
        ax.hlines(bar_y, x_left, x_right, color='k', lw=1)
        ax.text(
            (x_left + x_right) / 2,
            text_y,
            f'{pvals[i]:.3e}',
            ha='center',
            va='bottom',
            fontsize=3,
        )

    ax.set(
        xlabel='Time from run onset (s)',
        ylabel='Speed (cm/s)',
        xlim=(0, 4),
        ylim=YLIM_SPEED,
        title='post-matching speed',
    )
    ax.legend(frameon=False, fontsize=7)
    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)

    fig.tight_layout()
    print(f'\npost-matching speed binwise t-tests, t={label} s:')
    print('bin | t | p | dz')
    for i in range(N_BINS):
        lo = i * BIN_SIZE_MS
        hi = (i + 1) * BIN_SIZE_MS
        print(f'{lo:>4d} to {hi:<4d} ms | {tvals[i]:.3f} | {pvals[i]:.4f} | {dz[i]:.3f}')

    for ext in ['.png', '.pdf']:
        fig.savefig(stem / f'speed_post_matched{ext}', dpi=300, bbox_inches='tight')
    plt.close(fig)

    mean_E = np.nanmean(E, axis=1)
    mean_L = np.nanmean(L, axis=1)
    plot_violin_with_scatter(
        mean_E,
        mean_L,
        early_c,
        late_c,
        ylabel='Speed (cm/s)',
        xticklabels=['early', 'late'],
        title=f't = {threshold_text} s',
        stats_labels=[
            f'post-matching speed, first lick <{threshold_text} s',
            f'post-matching speed, first lick {threshold_text}-{late_max:g} s',
        ],
        print_statistics=True,
        save=True,
        savepath=stem / 'speed_post_matched_violin',
        show=False,
        close=True,
    )

def _plot_acceleration_post_matched(result, threshold, stem):
    threshold_text = f'{threshold:.2f}'.rstrip('0').rstrip('.')
    ea = np.asarray(result['sess_early_accel'], dtype=float)
    la = np.asarray(result['sess_late_accel'], dtype=float)
    valid_acceleration = np.isfinite(ea) & np.isfinite(la)
    ea, la = ea[valid_acceleration], la[valid_acceleration]
    if len(ea) and len(la):
        plot_violin_with_scatter(
            ea,
            la,
            early_c,
            late_c,
            ylabel='Acceleration (cm/s^2)',
            xticklabels=['early', 'late'],
            title=f't = {threshold_text} s',
            stats_labels=[
                f'post-matching initial acceleration, first lick <{threshold_text} s',
                f'post-matching initial acceleration, first lick {threshold_text}-{late_max:g} s',
            ],
            print_statistics=True,
            save=True,
            savepath=stem / 'accel_post_matched_violin',
            show=False,
            close=True,
        )

def _plot_one_cell_class(result, threshold, stem, label, early_key, late_key,
                         colour_early, colour_late):
    threshold_text = f'{threshold:.2f}'.rstrip('0').rstrip('.')
    early_profiles = result[early_key]
    late_profiles = result[late_key]
    if not early_profiles or not late_profiles:
        print(f'{label} plot skipped for t={threshold_text}: no cells')
        return

    early_arr = np.asarray(early_profiles)
    late_arr = np.asarray(late_profiles)
    early_mean = np.mean(early_arr, axis=0)
    early_sem = sem(early_arr, axis=0)
    late_mean = np.mean(late_arr, axis=0)
    late_sem = sem(late_arr, axis=0)

    fig, ax = plt.subplots(figsize=(2.3, 2.0))
    ax.plot(XAXIS, early_mean, c=colour_early, label=f'<{threshold_text}')
    ax.fill_between(
        XAXIS,
        early_mean + early_sem,
        early_mean - early_sem,
        color=colour_early,
        edgecolor='none',
        alpha=.25,
    )
    ax.plot(XAXIS, late_mean, c=colour_late, label=f'{threshold_text}-{late_max:g}')
    ax.fill_between(
        XAXIS,
        late_mean + late_sem,
        late_mean - late_sem,
        color=colour_late,
        edgecolor='none',
        alpha=.25,
    )

    p_ind, p_rs, p_rel, p_wil = flaf.compute_binwise_test_suite(
        early_arr,
        late_arr,
        samp_freq=SAMP_FREQ,
        bef=BEF,
        start=-0.5,
        end=3.5,
        bin_size=1,
        label=f'{label} cells, t={threshold_text} s',
        verbose=True,
    )
    flaf.annotate_binwise_test_suite(
        ax,
        p_ind,
        p_rs,
        p_rel,
        p_wil,
        start=-0.5,
        bin_size=1,
        fontsize=2.5,
    )

    ax.legend(fontsize=5, frameon=False)
    ax.set(
        xlabel='Time from run onset (s)',
        xlim=(-1, 4),
        ylabel='Firing rate (Hz)',
        ylim=(0.75, 3.15),
        yticks=[1, 2, 3],
        title=f't = {threshold_text} s',
    )
    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)

    fig.tight_layout()
    for ext in ['.png', '.pdf']:
        fig.savefig(
            stem / f'all_run_onset_{label}_mean_profiles{ext}',
            dpi=300,
            bbox_inches='tight',
        )
    plt.close(fig)

def _profile_delta(result, early_key, late_key):
    profile_start = int((SUMMARY_WINDOW[0] + BEF) * SAMP_FREQ)
    profile_stop = int((SUMMARY_WINDOW[1] + BEF) * SAMP_FREQ)
    early_profiles = np.asarray(result[early_key], dtype=float)
    late_profiles = np.asarray(result[late_key], dtype=float)
    if early_profiles.size == 0:
        early_vals = np.asarray([], dtype=float)
    else:
        early_vals = np.nanmean(early_profiles[:, profile_start:profile_stop], axis=1)
    if late_profiles.size == 0:
        late_vals = np.asarray([], dtype=float)
    else:
        late_vals = np.nanmean(late_profiles[:, profile_start:profile_stop], axis=1)
    valid_profiles = np.isfinite(early_vals) & np.isfinite(late_vals)
    early_vals, late_vals = early_vals[valid_profiles], late_vals[valid_profiles]
    if len(early_vals) == 0:
        return np.nan, np.nan, np.nan, np.nan
    delta = late_vals - early_vals
    delta_sem = sem(delta) if len(delta) > 1 else np.nan
    _, p_w = wilcoxon(early_vals, late_vals)
    _, p_t = ttest_rel(early_vals, late_vals)
    return np.nanmean(delta), delta_sem, p_w, p_t


#%% load data
print('loading dataframes...')

cell_profiles_path = HPC_stem / 'hpc_all_profiles.pkl'
cell_profiles = pd.read_pickle(cell_profiles_path)

cell_profiles_raphi_path = HPC_stem / 'hpc_all_profiles_raphi.pkl'
cell_profiles_raphi = pd.read_pickle(cell_profiles_raphi_path)

df_pyr = cell_profiles[cell_profiles['cell_identity'] == 'pyr']
df_pyr_raphi = cell_profiles_raphi[cell_profiles_raphi['cell_identity'] == 'pyr']

pyrup = pd.concat(
    [
        df_pyr[df_pyr['class'] == PYRUP_CLASS][['recname']],
        df_pyr_raphi[df_pyr_raphi['class'] == PYRUP_CLASS][['recname']],
    ],
    axis=0,
)
pyrdown = pd.concat(
    [
        df_pyr[df_pyr['class'] == PYRDOWN_CLASS][['recname']],
        df_pyr_raphi[df_pyr_raphi['class'] == PYRDOWN_CLASS][['recname']],
    ],
    axis=0,
)

active_recnames = set(recnames)
pyrup = pyrup[pyrup['recname'].isin(active_recnames)]
pyrdown = pyrdown[pyrdown['recname'].isin(active_recnames)]
pyrup_keys = set(pyrup.index)
pyrdown_keys = set(pyrdown.index)

all_valid_clunames = sorted(
    pyrup_keys | pyrdown_keys,
    key=lambda x: x.split(' ')[0],
)
clunames_by_rec = {}
for cluname in all_valid_clunames:
    recname = cluname.split(' ')[0]
    clunames_by_rec.setdefault(recname, []).append(cluname)


#%% main
results = {
    threshold: {
        'threshold': threshold,
        EARLY_PYRUP: [],
        LATE_PYRUP: [],
        EARLY_PYRDOWN: [],
        LATE_PYRDOWN: [],
        'sess_early_speed_raw': [],
        'sess_late_speed_raw': [],
        'sess_early_speed': [],
        'sess_late_speed': [],
        'sess_early_speed_scalar_raw': [],
        'sess_late_speed_scalar_raw': [],
        'sess_early_accel_raw': [],
        'sess_late_accel_raw': [],
        'sess_early_accel': [],
        'sess_late_accel': [],
        'raw_early_counts': [],
        'raw_late_counts': [],
        'matched_early_counts': [],
        'matched_late_counts': [],
        'n_sessions': 0,
        'n_pyrup': 0,
        'n_pyrdown': 0,
    }
    for threshold in thresholds
}

for recname in sorted(clunames_by_rec):
    print_session(recname)
    if recname not in active_recnames:
        print('skipped: not in active recording list')
        continue

    for maze_section in range(6):
        raphi_stem = (
            pp.RAPHAEL_ROOT / f'ANM{recname[1:4]}' / recname[:-3]
            / recname
        )
        align_path = (
            raphi_stem / f'{recname}_DataStructure_mazeSection1_'
            f'TrialType1_alignRun_msess{maze_section}.mat'
        )
        beh_path = (
            raphi_stem / f'{recname}_DataStructure_mazeSection1_'
            f'TrialType1_behPar_msess{maze_section}.mat'
        )
        if align_path.exists() and beh_path.exists():
            break
    else:
        mice_stem = (
            pp.MICEEXP_ROOT / f'ANMD{recname[1:5]}' / recname[:-3]
            / recname
        )
        align_path = pp.resolve_matlab_pipeline_file(
            mice_stem / f'{recname}_DataStructure_mazeSection1_'
            'TrialType1_alignRun_msess1.mat',
            recname,
        )
        beh_path = pp.resolve_matlab_pipeline_file(
            mice_stem / f'{recname}_DataStructure_mazeSection1_'
            'TrialType1_behPar_msess1.mat',
            recname,
        )
    alignRun = sio.loadmat(align_path)
    behPar = sio.loadmat(beh_path)

    licks = alignRun['trialsRun']['lickLfpInd'][0][0][0][1:]
    starts = alignRun['trialsRun']['startLfpInd'][0][0][0][1:]
    tot_trial = licks.shape[0]

    bad_idx = np.where(behPar['behPar'][0]['indTrBadBeh'][0] == 1)[1] - 1
    stim_idx = np.where(behPar['behPar'][0]['stimOn'][0] != 0)[1] - 1
    bad_set = set(bad_idx)
    stim_set = set(stim_idx)

    behaviour_paths = [
        stem / f'{recname}.pkl'
        for stem in [all_exp_stem, all_exp_term_stem, all_exp_raphi_stem]
    ]
    behaviour_path = next(
        (path for path in behaviour_paths if path.exists()),
        None,
    )
    if behaviour_path is None:
        raise FileNotFoundError(f'no behaviour file for {recname}')
    with open(behaviour_path, 'rb') as f:
        beh = pickle.load(f)

    speed_times = beh['speed_times_aligned'][1:]

    first_licks = []
    for trial in range(tot_trial):
        lk = [
            lick for lick in licks[trial]
            if (lick - starts[trial]).item() > .5 * SAMP_FREQ
        ]
        if len(lk) == 0:
            first_licks.append(np.nan)
        else:
            first_licks.append((lk[0] - starts[trial]).item() / SAMP_FREQ)

    matched_by_threshold = {}
    for threshold in thresholds:
        result = results[threshold]
        early_trials = []
        late_trials = []
        for trial, first_lick in enumerate(first_licks):
            if (
                trial in bad_set
                or trial in stim_set
                or trial - 1 in stim_set
                or np.isnan(first_lick)
            ):
                continue
            if first_lick < threshold:
                early_trials.append(trial)
            elif threshold < first_lick < late_max:
                late_trials.append(trial)
        result['raw_early_counts'].append(len(early_trials))
        result['raw_late_counts'].append(len(late_trials))

        if len(early_trials) < min_matched or len(late_trials) < min_matched:
            result['matched_early_counts'].append(0)
            result['matched_late_counts'].append(0)
            matched_by_threshold[threshold] = None
            continue

        raw_e_speed = flaf.compute_session_mean_speed(early_trials, speed_times, n=4000)
        raw_l_speed = flaf.compute_session_mean_speed(late_trials, speed_times, n=4000)
        if raw_e_speed is not None and raw_l_speed is not None:
            result['sess_early_speed_raw'].append(raw_e_speed)
            result['sess_late_speed_raw'].append(raw_l_speed)

        e_speed_scalar, e_accel_raw = _speed_scalar_and_accel(early_trials, speed_times)
        l_speed_scalar, l_accel_raw = _speed_scalar_and_accel(late_trials, speed_times)
        if len(e_speed_scalar) > min_matched and len(l_speed_scalar) > min_matched:
            result['sess_early_speed_scalar_raw'].append(np.nanmean(e_speed_scalar))
            result['sess_late_speed_scalar_raw'].append(np.nanmean(l_speed_scalar))
        if len(e_accel_raw) > min_matched and len(l_accel_raw) > min_matched:
            result['sess_early_accel_raw'].append(np.nanmean(e_accel_raw))
            result['sess_late_accel_raw'].append(np.nanmean(l_accel_raw))

        matched_early, matched_late = _speed_match_trials(
            early_trials,
            late_trials,
            speed_times,
        )
        result['matched_early_counts'].append(len(matched_early))
        result['matched_late_counts'].append(len(matched_late))

        if len(matched_early) < min_matched or len(matched_late) < min_matched:
            matched_by_threshold[threshold] = None
            continue

        matched_e_speed = flaf.compute_session_mean_speed(matched_early, speed_times, n=4000)
        matched_l_speed = flaf.compute_session_mean_speed(matched_late, speed_times, n=4000)
        if matched_e_speed is not None and matched_l_speed is not None:
            result['sess_early_speed'].append(matched_e_speed)
            result['sess_late_speed'].append(matched_l_speed)

        _, e_accel = _speed_scalar_and_accel(matched_early, speed_times)
        _, l_accel = _speed_scalar_and_accel(matched_late, speed_times)
        if len(e_accel) and len(l_accel):
            result['sess_early_accel'].append(np.nanmean(e_accel))
            result['sess_late_accel'].append(np.nanmean(l_accel))

        matched_by_threshold[threshold] = (matched_early, matched_late)
        threshold_text = f'{threshold:.2f}'.rstrip('0').rstrip('.')
        print(
            f't={threshold_text} s: '
            f'{len(early_trials)} early, {len(late_trials)} late; '
            f'{len(matched_early)} early, {len(matched_late)} late after matching'
        )

    if all(matched is None for matched in matched_by_threshold.values()):
        continue

    train_path = (
        HPC_stem / 'all_sessions'
        / recname / f'{recname}_all_trains_run.npy'
    )
    if not train_path.exists():
        train_path = (
            HPC_stem / 'all_sessions_raphi'
            / recname / f'{recname}_all_trains.npy'
        )
    all_trains = np.load(train_path, allow_pickle=True).item()

    profile_start = RUN_ONSET_BIN - BEF * SAMP_FREQ
    profile_stop = RUN_ONSET_BIN + AFT * SAMP_FREQ
    session_added = {threshold: False for threshold in thresholds}
    for cluname in clunames_by_rec[recname]:
        trains = all_trains[cluname]

        for threshold in thresholds:
            matched = matched_by_threshold[threshold]
            if matched is None:
                continue
            matched_early, matched_late = matched
            early_profile = np.nanmean(np.asarray(trains)[matched_early, profile_start:profile_stop], axis=0)
            late_profile = np.nanmean(np.asarray(trains)[matched_late, profile_start:profile_stop], axis=0)
            if int(recname[1:4]) > 40 and 'r' not in recname:
                early_profile = early_profile * SAMP_FREQ
                late_profile = late_profile * SAMP_FREQ

            result = results[threshold]
            if cluname in pyrup_keys:
                result[EARLY_PYRUP].append(early_profile)
                result[LATE_PYRUP].append(late_profile)
                result['n_pyrup'] += 1
                session_added[threshold] = True
            if cluname in pyrdown_keys:
                result[EARLY_PYRDOWN].append(early_profile)
                result[LATE_PYRDOWN].append(late_profile)
                result['n_pyrdown'] += 1
                session_added[threshold] = True

    for threshold, added in session_added.items():
        if added:
            results[threshold]['n_sessions'] += 1


#%% figures
summary_rows = []
saved_entries = []

for threshold in thresholds:
    threshold_text = f'{threshold:.2f}'.rstrip('0').rstrip('.')
    stem = first_lick_stem / f't_{threshold:.2f}'.replace('.', 'p')
    stem.mkdir(parents=True, exist_ok=True)
    result = results[threshold]

    print_statistics_section()
    print(f'threshold = {threshold_text} s')
    _plot_speed_pre_matched(result, threshold, stem)
    _plot_behaviour_raw_violins(result, threshold, stem)
    _plot_speed_post_matched(result, threshold, stem)
    _plot_acceleration_post_matched(result, threshold, stem)
    print(
        f'\nprofile summary for t={threshold_text} s: '
        f'{result["n_sessions"]} sessions, '
        f'{result["n_pyrup"]} PyrUp cells, {result["n_pyrdown"]} PyrDown cells'
    )
    _plot_one_cell_class(
        result,
        threshold,
        stem,
        'PyrUp',
        EARLY_PYRUP,
        LATE_PYRUP,
        'lightcoral',
        'firebrick',
    )
    _plot_one_cell_class(
        result,
        threshold,
        stem,
        'PyrDown',
        EARLY_PYRDOWN,
        LATE_PYRDOWN,
        'violet',
        'purple',
    )

    pyrup_delta, pyrup_sem, pyrup_w, pyrup_t = _profile_delta(
        result,
        EARLY_PYRUP,
        LATE_PYRUP,
    )
    pyrdown_delta, pyrdown_sem, pyrdown_w, pyrdown_t = _profile_delta(
        result,
        EARLY_PYRDOWN,
        LATE_PYRDOWN,
    )
    summary_rows.append(
        {
            'threshold_s': result['threshold'],
            'n_sessions': result['n_sessions'],
            'n_pyrup': result['n_pyrup'],
            'n_pyrdown': result['n_pyrdown'],
            'median_raw_early_trials': np.nanmedian(result['raw_early_counts']),
            'median_raw_late_trials': np.nanmedian(result['raw_late_counts']),
            'median_matched_early_trials': np.nanmedian(result['matched_early_counts']),
            'median_matched_late_trials': np.nanmedian(result['matched_late_counts']),
            'late_minus_early_pyrup_mean': pyrup_delta,
            'late_minus_early_pyrup_sem': pyrup_sem,
            'late_minus_early_pyrup_wilcoxon_p': pyrup_w,
            'late_minus_early_pyrup_ttest_p': pyrup_t,
            'late_minus_early_pyrdown_mean': pyrdown_delta,
            'late_minus_early_pyrdown_sem': pyrdown_sem,
            'late_minus_early_pyrdown_wilcoxon_p': pyrdown_w,
            'late_minus_early_pyrdown_ttest_p': pyrdown_t,
        }
    )
    saved_entries.append((f't={threshold_text} s', stem))

summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv(first_lick_stem / 'threshold_summary.csv', index=False)
for ext in ['.png', '.pdf']:
    summary_fig = first_lick_stem / f'threshold_summary{ext}'
    if summary_fig.exists():
        summary_fig.unlink()

saved_entries.append(('summary', first_lick_stem / 'threshold_summary.csv'))
print_files_saved(saved_entries)
