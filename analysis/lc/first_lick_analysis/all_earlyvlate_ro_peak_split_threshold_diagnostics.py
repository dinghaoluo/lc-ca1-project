# -*- coding: utf-8 -*-
'''
Created on Thu Jun 25 2026

test first-lick early/late split thresholds for LC run-onset cells.

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
from scipy.stats import sem, ttest_rel, wilcoxon

from common_functions import mpl_formatting
from console_formatting import print_files_saved, print_session, print_statistics_section
from plotting_functions import plot_violin_with_scatter
import first_lick_analysis_functions as flaf
import project_paths as pp
mpl_formatting()

import rec_list
paths = rec_list.pathLC
pathnames = [Path(p).name for p in paths]


#%% parameters
SAMP_FREQ = 1250
RUN_ONSET_BIN = 3750
BEF = 1
AFT = 4
WINDOW_HALF_SIZE = .5
RO_WINDOW = [
    int(RUN_ONSET_BIN - WINDOW_HALF_SIZE * SAMP_FREQ),
    int(RUN_ONSET_BIN + WINDOW_HALF_SIZE * SAMP_FREQ),
]

DEFAULT_THRESHOLDS = (2.2, 2.3, 2.4, 2.5, 2.6, 2.7, 2.8)
LATE_MAX = 3.5

BIN_SIZE_MS = 500
TOTAL_LEN_MS = 3500
N_BINS = TOTAL_LEN_MS // BIN_SIZE_MS
MATCH_K = 1
MIN_MATCHED = 5
ACC_WIN_MS = 1000

XAXIS = np.arange((BEF + AFT) * SAMP_FREQ) / SAMP_FREQ - BEF
X_SEC_PLOT = np.arange(4000) / 1000.0
YLIM_SPEED = (0, 70)

early_c = (168 / 255, 155 / 255, 202 / 255)
late_c  = (102 / 255, 83 / 255, 162 / 255)


#%% command line
parser = argparse.ArgumentParser(
    description='test alternative first-lick split thresholds for LC run-onset cells.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)
parser.add_argument(
    '--thresholds',
    nargs='+',
    type=float,
    default=list(DEFAULT_THRESHOLDS),
    help='candidate split thresholds in seconds',
)
parser.add_argument(
    '--late-max',
    type=float,
    default=LATE_MAX,
    help='upper bound for late first-lick trials, in seconds',
)
parser.add_argument(
    '--min-matched',
    type=int,
    default=MIN_MATCHED,
    help='minimum early and late speed-matched trials required per session',
)
args = parser.parse_args()
thresholds = sorted(set(args.thresholds))
late_max = args.late_max
min_matched = args.min_matched

first_lick_stem = (
    pp.LC_EPHYS_FIGURES_STEM
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

    e_low, e_high = e_mu - MATCH_K * e_sd, e_mu + MATCH_K * e_sd
    l_low, l_high = l_mu - MATCH_K * l_sd, l_mu + MATCH_K * l_sd

    late_mask = np.all((L_bins >= e_low) & (L_bins <= e_high), axis=1)
    early_mask = np.all((E_bins >= l_low) & (E_bins <= l_high), axis=1)

    matched_late = [l_valid[i] for i in np.where(late_mask)[0]]
    matched_early = [e_valid[i] for i in np.where(early_mask)[0]]
    return matched_early, matched_late

def _trial_mean_accel(trial_list, speed_times):
    out = []
    n_needed = ACC_WIN_MS + 1
    for trial in trial_list:
        speeds = [pt[1] for pt in speed_times[trial]]
        if len(speeds) < n_needed:
            continue
        s = np.asarray(speeds[:n_needed], dtype=float)
        out.append(np.nanmean(np.diff(s) * 1000.0))

    return np.asarray(out, dtype=float)

def _plot_speed_trace(result, threshold, stem, matched=False):
    key0 = 'sess_early_speed' if matched else 'sess_early_speed_raw'
    key1 = 'sess_late_speed' if matched else 'sess_late_speed_raw'
    title = 'post-matching speed' if matched else 'pre-matching speed'
    filename = 'speed_post_matched' if matched else 'speed_pre_matched'

    label = f'{threshold:.2f}'.rstrip('0').rstrip('.')
    if not result[key0] or not result[key1]:
        print(f'{title} skipped for t={label}: no speed traces')
        return

    E = np.vstack(result[key0])
    L = np.vstack(result[key1])
    E_mean = np.mean(E, axis=0)
    E_sem = sem(E, axis=0)
    L_mean = np.mean(L, axis=0)
    L_sem = sem(L, axis=0)

    fig, ax = plt.subplots(figsize=(2.1, 2.0))
    ax.plot(X_SEC_PLOT, E_mean, c=early_c, label=f'early (<{label} s)')
    ax.fill_between(
        X_SEC_PLOT,
        E_mean + E_sem,
        E_mean - E_sem,
        color=early_c,
        edgecolor='none',
        alpha=.25,
    )

    ax.plot(X_SEC_PLOT, L_mean, c=late_c, label=f'late ({label}-{late_max:g} s)')
    ax.fill_between(
        X_SEC_PLOT,
        L_mean + L_sem,
        L_mean - L_sem,
        color=late_c,
        edgecolor='none',
        alpha=.25,
    )

    if matched:
        E_binmean = np.nanmean(E, axis=1)
        L_binmean = np.nanmean(L, axis=1)
        e_mean = np.mean(E_binmean)
        e_sem = sem(E_binmean)
        l_mean = np.mean(L_binmean)
        l_sem = sem(L_binmean)
        e_med = np.median(E_binmean)
        l_med = np.median(L_binmean)
        e_q25, e_q75 = np.percentile(E_binmean, [25, 75])
        l_q25, l_q75 = np.percentile(L_binmean, [25, 75])
        _, p_t = ttest_rel(E_binmean, L_binmean)
        _, p_w = wilcoxon(E_binmean, L_binmean)
        p_t_text = f'{p_t:.3g}' if p_t < 0.01 else f'{p_t:.3f}'
        p_w_text = f'{p_w:.3g}' if p_w < 0.01 else f'{p_w:.3f}'

        stats_txt = (
            f'mean speed (0-4 s)\n'
            f'early: mean {e_mean:.2f} +/- {e_sem:.2f}\n'
            f'       med  {e_med:.2f} [{e_q25:.2f}, {e_q75:.2f}]\n'
            f'late:  mean {l_mean:.2f} +/- {l_sem:.2f}\n'
            f'       med  {l_med:.2f} [{l_q25:.2f}, {l_q75:.2f}]\n'
            f'ttest p = {p_t_text}\n'
            f'wilc  p = {p_w_text}'
        )
        ax.text(
            0.02,
            0.98,
            stats_txt,
            transform=ax.transAxes,
            ha='left',
            va='top',
            fontsize=6,
        )

    ax.set(
        xlabel='Time from run onset (s)',
        ylabel='Speed (cm/s)',
        xlim=(0, 4),
        ylim=YLIM_SPEED,
        title=title,
    )
    ax.legend(frameon=False, fontsize=7)
    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)

    fig.tight_layout()
    for ext in ['.png', '.pdf']:
        fig.savefig(stem / f'{filename}{ext}', dpi=300, bbox_inches='tight')
    plt.close(fig)
    return E, L

def _plot_profile(result, threshold, stem):
    early_profiles = np.asarray(result['early_profiles'])
    late_profiles = np.asarray(result['late_profiles'])
    label = f'{threshold:.2f}'.rstrip('0').rstrip('.')
    if early_profiles.size == 0 or late_profiles.size == 0:
        print(f'profile plot skipped for t={label}: no cells')
        return

    early_mean = np.mean(early_profiles, axis=0)
    early_sem = sem(early_profiles, axis=0)
    late_mean = np.mean(late_profiles, axis=0)
    late_sem = sem(late_profiles, axis=0)

    fig, ax = plt.subplots(figsize=(2.2, 2.1))
    ax.plot(XAXIS, early_mean, c=early_c, label=f'<{label} s')
    ax.fill_between(
        XAXIS,
        early_mean + early_sem,
        early_mean - early_sem,
        color=early_c,
        edgecolor='none',
        alpha=.25,
    )
    ax.plot(XAXIS, late_mean, c=late_c, label=f'{label}-{late_max:g} s')
    ax.fill_between(
        XAXIS,
        late_mean + late_sem,
        late_mean - late_sem,
        color=late_c,
        edgecolor='none',
        alpha=.25,
    )

    a = np.asarray(result['early_rates'], dtype=float)
    b = np.asarray(result['late_rates'], dtype=float)
    paired = np.isfinite(a) & np.isfinite(b)
    a, b = a[paired], b[paired]
    if len(a) and len(b):
        a_mean = np.mean(a)
        a_sem = sem(a) if len(a) > 1 else np.nan
        b_mean = np.mean(b)
        b_sem = sem(b) if len(b) > 1 else np.nan
        a_med = np.median(a)
        b_med = np.median(b)
        a_q25, a_q75 = np.percentile(a, [25, 75])
        b_q25, b_q75 = np.percentile(b, [25, 75])

        _, p_w = wilcoxon(a, b)
        _, p_t = ttest_rel(a, b)
        p_w_text = f'{p_w:.3g}' if p_w < 0.01 else f'{p_w:.3f}'
        p_t_text = f'{p_t:.3g}' if p_t < 0.01 else f'{p_t:.3f}'

        stats_txt = (
            f'early: mean {a_mean:.2f} +/- {a_sem:.2f}\n'
            f'       med  {a_med:.2f} [{a_q25:.2f}, {a_q75:.2f}]\n'
            f'late:  mean {b_mean:.2f} +/- {b_sem:.2f}\n'
            f'       med  {b_med:.2f} [{b_q25:.2f}, {b_q75:.2f}]\n'
            f'wilc  p = {p_w_text}\n'
            f'ttest p = {p_t_text}'
        )
        ax.text(
            0.02,
            0.98,
            stats_txt,
            transform=ax.transAxes,
            ha='left',
            va='top',
            fontsize=6,
        )

    ax.legend(fontsize=7, frameon=False)
    ax.set(
        xlabel='Time from run-onset (s)',
        xlim=(-1, 4),
        ylabel='Firing rate (Hz)',
        title=f't = {label} s',
    )
    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)

    fig.tight_layout()
    for ext in ['.png', '.pdf']:
        fig.savefig(
            stem / f'all_run_onset_mean_profiles_early_v_late{ext}',
            dpi=300,
            bbox_inches='tight',
        )
    plt.close(fig)

def _plot_accel(result, threshold, stem):
    ea = np.asarray(result['sess_early_accel'], dtype=float)
    la = np.asarray(result['sess_late_accel'], dtype=float)
    paired = np.isfinite(ea) & np.isfinite(la)
    ea, la = ea[paired], la[paired]
    label = f'{threshold:.2f}'.rstrip('0').rstrip('.')
    if len(ea) == 0 or len(la) == 0:
        print(f'acceleration plot skipped for t={label}: no paired sessions')
        return

    print_statistics_section()
    print(f'threshold = {label} s')
    plot_violin_with_scatter(
        data0=ea,
        data1=la,
        colour0=early_c,
        colour1=late_c,
        paired=True,
        alpha=.3,
        xticklabels=['early', 'late'],
        ylabel='Mean acceleration (cm/s$^2$)',
        title=f't = {label} s',
        showscatter=False,
        print_statistics=True,
        plot_statistics=True,
        figsize=(1.8, 2.2),
        save=True,
        savepath=stem / 'post_matching_acceleration_0to1s',
        dpi=300,
        show=False,
        close=True,
    )


#%% load cells
print('loading data...')
cell_prop = pd.read_pickle(pp.LC_EPHYS_STEM / 'LC_all_cell_profiles.pkl')

RO_keys = []
for clu in cell_prop.itertuples():
    if clu.identity in ['tagged', 'putative'] and clu.run_onset_peak:
        RO_keys.append(clu.Index)

RO_keys = [
    cluname for cluname in RO_keys
    if cluname.split(' ')[0] in pathnames
]

ro_keys_by_rec = {}
for cluname in RO_keys:
    recname = cluname.split(' ')[0]
    ro_keys_by_rec.setdefault(recname, []).append(cluname)

#%% main
results = {
    threshold: {
        'threshold': threshold,
        'early_profiles': [],
        'late_profiles': [],
        'early_rates': [],
        'late_rates': [],
        'sess_early_speed_raw': [],
        'sess_late_speed_raw': [],
        'sess_early_speed': [],
        'sess_late_speed': [],
        'sess_early_accel': [],
        'sess_late_accel': [],
        'raw_early_counts': [],
        'raw_late_counts': [],
        'matched_early_counts': [],
        'matched_late_counts': [],
        'n_sessions': 0,
        'n_cells': 0,
    }
    for threshold in thresholds
}

for recname in sorted(ro_keys_by_rec):
    print_session(recname)

    rec_stem = pp.MICEEXP_ROOT / f'ANMD{recname[1:5]}' / recname[:14] / recname
    alignRun = sio.loadmat(
        rec_stem / f'{recname}_DataStructure_mazeSection1_TrialType1_alignRun_msess1.mat'
    )
    licks = alignRun['trialsRun']['lickLfpInd'][0][0][0][1:]
    starts = alignRun['trialsRun']['startLfpInd'][0][0][0][1:]
    tot_trial = licks.shape[0]

    behPar = sio.loadmat(
        rec_stem / f'{recname}_DataStructure_mazeSection1_TrialType1_behPar_msess1.mat'
    )
    bad_idx = np.where(behPar['behPar'][0]['indTrBadBeh'][0] == 1)[1] - 1
    stim_idx = np.where(behPar['behPar'][0]['stimOn'][0] == 1)[1] - 1

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

    all_trains = np.load(
        pp.LC_EPHYS_STEM / 'all_sessions' / recname / f'{recname}_all_trains.npy',
        allow_pickle=True,
    ).item()
    beh_path = pp.behaviour_experiment_stem('LC') / f'{recname}.pkl'
    if not beh_path.exists():
        beh_path = pp.behaviour_experiment_stem('LCterm') / f'{recname}.pkl'
    with open(beh_path, 'rb') as f:
        beh = pickle.load(f)
    speed_times = beh['speed_times_aligned'][1:]
    bad_set = set(bad_idx)
    stim_set = set(stim_idx)

    for threshold in thresholds:
        label = f'{threshold:.2f}'.rstrip('0').rstrip('.')
        result = results[threshold]
        early_trials = []
        late_trials = []
        for trial, t in enumerate(first_licks):
            if trial in bad_set or trial in stim_set or np.isnan(t):
                continue
            if t < threshold:
                early_trials.append(trial)
            elif threshold < t < late_max:
                late_trials.append(trial)
        result['raw_early_counts'].append(len(early_trials))
        result['raw_late_counts'].append(len(late_trials))

        raw_e_speed = flaf.compute_session_mean_speed(early_trials, speed_times, n=4000)
        raw_l_speed = flaf.compute_session_mean_speed(late_trials, speed_times, n=4000)
        if raw_e_speed is not None and raw_l_speed is not None:
            result['sess_early_speed_raw'].append(raw_e_speed)
            result['sess_late_speed_raw'].append(raw_l_speed)

        matched_early, matched_late = _speed_match_trials(
            early_trials,
            late_trials,
            speed_times,
        )
        result['matched_early_counts'].append(len(matched_early))
        result['matched_late_counts'].append(len(matched_late))

        print(
            f't={label} s: '
            f'{len(early_trials)} early, {len(late_trials)} late; '
            f'{len(matched_early)} early, {len(matched_late)} late after matching'
        )

        # LC profiles follow the fixed first-lick split; speed matching stays diagnostic.
        if len(early_trials) and len(late_trials):
            added_cells = 0
            for cluname in ro_keys_by_rec[recname]:
                trains = all_trains[cluname]

                early_profiles, early_rates = flaf.extract_run_onset_profiles(
                    trains,
                    early_trials,
                    RUN_ONSET_BIN,
                    SAMP_FREQ,
                    BEF,
                    AFT,
                    rate_window=RO_WINDOW,
                )
                late_profiles, late_rates = flaf.extract_run_onset_profiles(
                    trains,
                    late_trials,
                    RUN_ONSET_BIN,
                    SAMP_FREQ,
                    BEF,
                    AFT,
                    rate_window=RO_WINDOW,
                )
                if not len(early_profiles) or not len(late_profiles):
                    continue

                result['early_profiles'].append(np.nanmean(early_profiles, axis=0))
                result['late_profiles'].append(np.nanmean(late_profiles, axis=0))
                result['early_rates'].append(np.nanmean(early_rates))
                result['late_rates'].append(np.nanmean(late_rates))
                added_cells += 1

            if added_cells:
                result['n_sessions'] += 1
                result['n_cells'] += added_cells

        if len(matched_early) < min_matched or len(matched_late) < min_matched:
            continue

        matched_e_speed = flaf.compute_session_mean_speed(matched_early, speed_times, n=4000)
        matched_l_speed = flaf.compute_session_mean_speed(matched_late, speed_times, n=4000)
        if matched_e_speed is not None and matched_l_speed is not None:
            result['sess_early_speed'].append(matched_e_speed)
            result['sess_late_speed'].append(matched_l_speed)

        e_acc = _trial_mean_accel(matched_early, speed_times)
        l_acc = _trial_mean_accel(matched_late, speed_times)
        if len(e_acc) and len(l_acc):
            result['sess_early_accel'].append(np.nanmean(e_acc))
            result['sess_late_accel'].append(np.nanmean(l_acc))


#%% figures
summary_rows = []
saved_entries = []
for threshold in thresholds:
    label = f'{threshold:.2f}'.rstrip('0').rstrip('.')
    stem = first_lick_stem / f't_{threshold:.2f}'.replace('.', 'p')
    stem.mkdir(parents=True, exist_ok=True)
    result = results[threshold]

    _plot_speed_trace(result, threshold, stem, matched=False)
    _plot_speed_trace(result, threshold, stem, matched=True)
    _plot_profile(result, threshold, stem)
    _plot_accel(result, threshold, stem)

    a = np.asarray(result['early_rates'], dtype=float)
    b = np.asarray(result['late_rates'], dtype=float)
    paired = np.isfinite(a) & np.isfinite(b)
    a, b = a[paired], b[paired]
    if len(a) and len(b):
        delta = b - a
        delta_mean = np.nanmean(delta)
        delta_sem = sem(delta) if len(delta) > 1 else np.nan
        _, p_w = wilcoxon(a, b)
        _, p_t = ttest_rel(a, b)
    else:
        delta_mean = np.nan
        delta_sem = np.nan
        p_w = np.nan
        p_t = np.nan

    summary_rows.append({
        'threshold_s': result['threshold'],
        'n_sessions': result['n_sessions'],
        'n_cells': result['n_cells'],
        'median_raw_early_trials': np.nanmedian(result['raw_early_counts']),
        'median_raw_late_trials': np.nanmedian(result['raw_late_counts']),
        'median_matched_early_trials': np.nanmedian(result['matched_early_counts']),
        'median_matched_late_trials': np.nanmedian(result['matched_late_counts']),
        'late_minus_early_rate_mean': delta_mean,
        'late_minus_early_rate_sem': delta_sem,
        'late_minus_early_wilcoxon_p': p_w,
        'late_minus_early_ttest_p': p_t,
    })
    saved_entries.append((f't={label} s', stem))

summary_df = pd.DataFrame(summary_rows)
summary_df.to_csv(first_lick_stem / 'threshold_summary.csv', index=False)
saved_entries.append(('summary', first_lick_stem / 'threshold_summary.csv'))
print_files_saved(saved_entries)
