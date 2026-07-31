# -*- coding: utf-8 -*-
'''
Created on Mon Oct 13 12:03:43 2025
Originally named amp_since_last_reward.py
Modified on Friday to get reward-aligned firing profiles

Compare LC run-onset amplitudes between trials with low and high 'time since
    last reward' according to GLM results

@author: Dinghao Luo
'''

#%% imports
from pathlib import Path
import sys

repo_root = Path(__file__).resolve().parents[3]

if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
utils_path = repo_root / 'utils'
if str(utils_path) not in sys.path:
    sys.path.insert(0, str(utils_path))

import numpy as np
import pandas as pd
import pickle
import scipy.io as sio
from scipy.stats import sem, linregress, ttest_1samp, wilcoxon
import matplotlib.pyplot as plt

from common_functions import mpl_formatting
from console_formatting import print_session
import glm_functions as gf
from plotting_functions import plot_violin_with_scatter
import project_paths as pp
mpl_formatting()

import rec_list
paths = rec_list.pathLC


#%% paths and parameters
all_sess_stem = pp.LC_EPHYS_STEM / 'all_sessions'
LC_beh_stem   = pp.behaviour_experiment_stem('LC')
GLM_stem      = pp.LC_EPHYS_FIGURES_STEM / 'GLM'
GLM_stem.mkdir(parents=True, exist_ok=True)
for fig_dir in [
    GLM_stem / 'rew_to_run_single_session_split',
    GLM_stem / 'rew_to_run_single_session',
]:
    fig_dir.mkdir(parents=True, exist_ok=True)

SAMP_FREQ = 1250
SAMP_FREQ_BEH = 1000
RUN_ONSET_IDX = 3 * SAMP_FREQ
BURST_WINDOW = (-.5, .5)  # for amplitude

PERMS = 1000  # permutate for 1000 times (per session) for signif test

# fine 0.1-s bins for time-since-last-reward traces
BIN_START = 1.0
BIN_END = 2.1
BIN_WIDTH = 0.1
N_BINS = int((BIN_END - BIN_START) / BIN_WIDTH)

save = True

BASELINE_OFFSET = 0.25
SMOOTH_WINDOW = 0.05
MIN_TRIALS_PER_BIN = 3


#%% session scatter plots
def _plot_session_scatter_pair(recname, trial_times, trial_rates, regression, suffix=''):
    fig, axes = plt.subplots(1, 2, figsize=(3.2, 2.2), sharex=True, sharey=True)

    axes[0].scatter(trial_times, trial_rates, s=10, color='forestgreen', ec='none', alpha=0.7)
    axes[0].plot(regression['xfit'], regression['yfit'], color='black', lw=1)
    axes[0].text(0.05, 0.95, f'r = {regression["r"]:.2f}\np = {regression["p"]:.3f}',
                 transform=axes[0].transAxes, ha='left', va='top',
                 fontsize=7, color='black')
    axes[0].set(xlabel='Time since rew. (s)', ylabel='Run-onset FR (Hz)')
    axes[0].spines[['top', 'right']].set_visible(False)

    axes[1].scatter(trial_times, regression['trial_rates_shuf'], s=10, color='gray', ec='none', alpha=0.6)
    axes[1].plot(regression['xfit'], regression['yfit_shuf'], color='black', lw=1)
    axes[1].text(0.05, 0.95, f'r = {regression["r_shuf"]:.2f}\np = {regression["p_shuf"]:.3f}',
                 transform=axes[1].transAxes, ha='left', va='top',
                 fontsize=7, color='black')
    axes[1].set(xlabel='Time since rew. (s)')
    axes[1].spines[['top', 'right']].set_visible(False)

    if suffix == '_lim':
        for ax in axes:
            ax.set(xlim=(.5, 1.5))

    fig.suptitle(recname)
    plt.tight_layout()

    if save:
        for ext in ['.pdf', '.png']:
            fig.savefig(
                GLM_stem / 'rew_to_run_single_session' / f'{recname}{suffix}{ext}',
                dpi=300,
                bbox_inches='tight',
            )
    plt.close(fig)


#%% analysis
cell_prop = pd.read_pickle(pp.LC_EPHYS_STEM / 'LC_all_cell_profiles.pkl')

all_t_since = []
low_prof_list = []
high_prof_list = []
low_prof_rew_list = []
high_prof_rew_list = []
regress_r = []
regress_shuf_r = []
all_prof_rew_list_bins = [[] for _ in range(N_BINS)]
all_cell_rate_dict = {}
all_cell_time_dict = {}
cell_bin_traces = {}

for path in paths:
    recname = Path(path).name
    print_session(recname)

    with open(LC_beh_stem / f'{recname}.pkl', 'rb') as f:
        beh = pickle.load(f)
    reward_times = beh['reward_times'][1:]
    run_onsets = beh['run_onsets'][1:]
    trials_sts = beh['trial_statements'][1:]
    opto_idx = [i for i, trial in enumerate(trials_sts) if trial[15] != '0']

    rec_stem = pp.MICEEXP_ROOT / f'ANMD{recname[1:5]}' / recname[:14] / recname
    # reward-aligned LFP samples index the full spike maps
    aligned_rew = sio.loadmat(
        rec_stem / f'{recname}_DataStructure_mazeSection1_TrialType1_alignRew_msess1.mat'
    )['trialsRew'][0][0]
    rew_spike = aligned_rew['startLfpInd'][0][1:]
    spike_maps = np.load(
        all_sess_stem / recname / f'{recname}_smoothed_spike_map.npy',
        allow_pickle=True,
    )
    all_trains = np.load(
        all_sess_stem / recname / f'{recname}_all_trains.npy',
        allow_pickle=True,
    ).item()
    curr_cell_prop = cell_prop[cell_prop['sessname'] == recname]
    eligible_cells = [
        cluname for cluname, row in curr_cell_prop.iterrows()
        if row['identity'] != 'other' and row['run_onset_peak']
    ]
    if not eligible_cells:
        print('no eligible cells; skipped')
        continue

    for cluname in eligible_cells:
        if cluname not in cell_bin_traces:
            cell_bin_traces[cluname] = [[] for _ in range(N_BINS)]

    prof_rew_list_bins = [[] for _ in range(N_BINS)]
    trial_counts = [0] * N_BINS
    # exclude optogenetic trials, the immediately following trial, and trials without a run onset
    valid_trials = [
        trial for trial, ro in enumerate(run_onsets[:-1])
        if (
            trial not in opto_idx
            and trial - 1 not in opto_idx
            and not np.isnan(ro)
        )
    ]

    trial_times = []
    trial_rates = []
    cell_rate_dict = {clu: [] for clu in eligible_cells}
    cell_time_dict = {clu: [] for clu in eligible_cells}

    for ti in valid_trials:
        onset_time = run_onsets[ti] / SAMP_FREQ_BEH
        t_since = gf.time_since_last_reward(reward_times, onset_time, ti)
        all_t_since.append(t_since)

        if np.isnan(t_since) or t_since < 0 or t_since > 8:
            continue

        bin_idx = int((t_since - BIN_START) // BIN_WIDTH) if BIN_START <= t_since < BIN_END else None
        if bin_idx is not None:
            trial_counts[bin_idx] += 1

        curr_rates = []
        cell_rates = {}
        cell_rew_full = {}
        cell_rew_short = {}
        for cluname in eligible_cells:
            clu_idx = int(cluname.split('clu')[-1]) - 2
            train = all_trains[cluname][ti]
            curr_rate = np.mean(
                train[
                    RUN_ONSET_IDX + int(BURST_WINDOW[0] * SAMP_FREQ):
                    RUN_ONSET_IDX + int(BURST_WINDOW[1] * SAMP_FREQ)
                ]
            )
            curr_rates.append(curr_rate)
            cell_rates[cluname] = curr_rate

            # fixed 10-s reward-aligned trace for the low/high interval comparison
            cell_rew_full[cluname] = spike_maps[clu_idx][
                rew_spike[ti] - 3 * SAMP_FREQ: rew_spike[ti] + 7 * SAMP_FREQ
            ]
            # shorter reward-aligned trace stops at the current run onset
            cell_rew_short[cluname] = spike_maps[clu_idx][
                rew_spike[ti] - 1 * SAMP_FREQ: rew_spike[ti] + int(t_since * SAMP_FREQ)
            ]
        mean_curr_rate = np.mean(curr_rates)

        for cluname in eligible_cells:
            cell_rate_dict[cluname].append(cell_rates[cluname])
            cell_time_dict[cluname].append(t_since)

            if t_since < 1.5:
                low_prof_list.append(all_trains[cluname][ti])
                low_prof_rew_list.append(cell_rew_full[cluname])
            else:
                high_prof_list.append(all_trains[cluname][ti])
                high_prof_rew_list.append(cell_rew_full[cluname])

            if bin_idx is not None:
                prof_rew_list_bins[bin_idx].append(cell_rew_short[cluname])
                cell_bin_traces[cluname][bin_idx].append(cell_rew_short[cluname])

        trial_times.append(t_since)
        trial_rates.append(mean_curr_rate)

    # require ten trials across the interval, then plot only bins with at least three
    if sum(trial_counts) >= 10:
        fig, ax = plt.subplots(figsize=(3, 2))
        for bi in range(N_BINS):
            if trial_counts[bi] < 3:
                continue
            bin_trials = prof_rew_list_bins[bi]
            bin_min = min(len(trial) for trial in bin_trials)
            xaxis_bin = np.arange(bin_min) / SAMP_FREQ - 1
            prof_rew_array_bin = np.array([trial[:bin_min] for trial in bin_trials])
            mean_prof_rew_bin = np.mean(prof_rew_array_bin, axis=0)
            colour = plt.cm.Greens(0.3 + 0.6 * bi / (N_BINS - 1))
            ax.plot(xaxis_bin, mean_prof_rew_bin, color=colour)

        fig.suptitle(recname)
        plt.tight_layout()
        if save:
            for ext in ['.pdf', '.png']:
                fig.savefig(
                    GLM_stem / 'rew_to_run_single_session_split' / f'{recname}{ext}',
                    dpi=300,
                    bbox_inches='tight',
                )
        plt.close(fig)
        for bi in range(N_BINS):
            all_prof_rew_list_bins[bi].extend(prof_rew_list_bins[bi])
    else:
        print('not enough trials; skipped reward-binned plot for this session')

    all_cell_rate_dict.update(cell_rate_dict)
    all_cell_time_dict.update(cell_time_dict)

    slope, intercept, r, p, _ = linregress(trial_times, trial_rates)
    xfit = np.linspace(min(trial_times), max(trial_times), 2)
    yfit = intercept + slope * xfit

    curr_r_shuf = []
    # permutations stay within session
    for _ in range(PERMS):
        trial_times_shuf = np.random.permutation(trial_times)
        trial_rates_shuf = np.random.permutation(trial_rates)
        slope_shuf, intercept_shuf, r_shuf, p_shuf, _ = linregress(
            trial_times_shuf,
            trial_rates_shuf,
        )
        curr_r_shuf.append(r_shuf)
        yfit_shuf = intercept_shuf + slope_shuf * xfit

    regression = {
        'r': r,
        'p': p,
        'xfit': xfit,
        'yfit': yfit,
        'trial_rates_shuf': trial_rates_shuf,
        'yfit_shuf': yfit_shuf,
        'r_shuf': r_shuf,
        'p_shuf': p_shuf,
    }
    regress_r.append(r)
    regress_shuf_r.append(np.mean(curr_r_shuf))
    _plot_session_scatter_pair(recname, trial_times, trial_rates, regression, suffix='')
    _plot_session_scatter_pair(recname, trial_times, trial_rates, regression, suffix='_lim')

low_prof_list = np.array(low_prof_list)
high_prof_list = np.array(high_prof_list)

cell_ramp_rates = {clu: np.full(N_BINS, np.nan) for clu in cell_bin_traces}
cell_baselines = {clu: np.full(N_BINS, np.nan) for clu in cell_bin_traces}
cell_ends = {clu: np.full(N_BINS, np.nan) for clu in cell_bin_traces}
for cluname, bin_lists in cell_bin_traces.items():
    for bi in range(N_BINS):
        trials = bin_lists[bi]
        if len(trials) < MIN_TRIALS_PER_BIN:
            continue

        bin_min = min(len(tr) for tr in trials)
        arr = np.stack([tr[:bin_min] for tr in trials])
        mean_prof = np.mean(arr, axis=0)
        xaxis = np.arange(bin_min) / SAMP_FREQ - 1

        t_end = xaxis[-1]
        t_base = t_end - BASELINE_OFFSET
        if t_base <= xaxis[0]:
            continue

        idx_base = np.where(
            (xaxis >= (t_base - SMOOTH_WINDOW)) &
            (xaxis <= (t_base + SMOOTH_WINDOW))
        )[0]
        idx_end = np.where(
            (xaxis >= (t_end - SMOOTH_WINDOW)) &
            (xaxis <= (t_end + SMOOTH_WINDOW))
        )[0]
        baseline = np.nanmean(mean_prof[idx_base])
        end = np.nanmean(mean_prof[idx_end])
        slope = (end - baseline) / BASELINE_OFFSET

        cell_baselines[cluname][bi] = baseline
        cell_ramp_rates[cluname][bi] = slope
        cell_ends[cluname][bi] = end

xaxis_bins = []
mean_prof_rew_bins = []
sem_prof_rew_bins = []
for bi in range(N_BINS):
    bin_trials = all_prof_rew_list_bins[bi]
    if len(bin_trials) == 0:
        xaxis_bins.append(None)
        mean_prof_rew_bins.append(None)
        sem_prof_rew_bins.append(None)
        continue

    bin_end = BIN_START + bi * BIN_WIDTH
    duration = bin_end + 1.0
    n_samples = int(duration * SAMP_FREQ)
    arr = []
    for trial in bin_trials:
        if len(trial) >= n_samples:
            arr.append(trial[:n_samples])
        else:
            arr.append(np.pad(trial, (0, n_samples - len(trial)), mode='edge'))
    arr = np.vstack(arr)
    mean_prof_rew_bins.append(np.mean(arr, axis=0))
    sem_prof_rew_bins.append(sem(arr, axis=0))
    xaxis_bins.append(np.arange(n_samples) / SAMP_FREQ - 1)

mean_low_prof = np.mean(low_prof_list, axis=0)
sem_low_prof = sem(low_prof_list, axis=0)
mean_high_prof = np.mean(high_prof_list, axis=0)
sem_high_prof = sem(high_prof_list, axis=0)

xaxis = np.arange(len(mean_low_prof)) / SAMP_FREQ - 3
fig, ax = plt.subplots(figsize=(2.4, 2.0))
ax.plot(xaxis, mean_low_prof, color='forestgreen', alpha=.5, label='Shorter time')
ax.fill_between(
    xaxis,
    mean_low_prof - sem_low_prof,
    mean_low_prof + sem_low_prof,
    color='forestgreen',
    alpha=.15,
)
ax.plot(xaxis, mean_high_prof, color='forestgreen', label='Longer time')
ax.fill_between(
    xaxis,
    mean_high_prof - sem_high_prof,
    mean_high_prof + sem_high_prof,
    color='forestgreen',
    alpha=.3,
)
ax.set(
    xlabel='Time from run onset (s)',
    xlim=(-1, 4),
    ylabel='Firing rate (Hz)',
    ylim=(1.6, 5),
)
ax.legend(frameon=False)
ax.spines[['top', 'right']].set_visible(False)
plt.tight_layout()
plt.show()
if save:
    for ext in ['.pdf', '.png']:
        fig.savefig(GLM_stem / f'rew_to_run_profiles{ext}', dpi=300, bbox_inches='tight')

tval, p_t = ttest_1samp(regress_r, 0)
wstat, p_w = wilcoxon(regress_r)
q25, q75 = np.percentile(regress_r, [25, 75])
shuf_mean = np.nanmean(regress_shuf_r)
shuf_std = np.nanstd(regress_shuf_r)
ci_low = shuf_mean - 1.96 * shuf_std
ci_high = shuf_mean + 1.96 * shuf_std

fig, ax = plt.subplots(figsize=(1.6, 2.2))
parts = ax.violinplot(
    regress_r,
    positions=[1],
    showmeans=False,
    showmedians=True,
    showextrema=False,
)
for pc in parts['bodies']:
    pc.set_facecolor('forestgreen')
    pc.set_edgecolor('none')
    pc.set_alpha(0.35)
parts['cmedians'].set_color('k')
parts['cmedians'].set_linewidth(1.2)
ax.scatter(
    np.ones(len(regress_r)),
    regress_r,
    color='forestgreen',
    ec='none',
    s=10,
    alpha=0.5,
    zorder=3,
)
ax.axhline(shuf_mean, color='gray', lw=1.2, ls='--')
ax.fill_between(
    [0.5, 1.5],
    [ci_low, ci_low],
    [ci_high, ci_high],
    color='gray',
    alpha=0.20,
    edgecolor='none',
)
mean_r = np.nanmean(regress_r)
sem_r = sem(regress_r)
ymax = np.max(regress_r)
ymin = np.min(regress_r)
ax.text(
    1,
    ymax + 0.05 * (ymax - ymin),
    f'Med = {np.median(regress_r):.2f}\n'
    f'IQR = [{q25:.2f}, {q75:.2f}]\n'
    f'{mean_r:.2f} +/- {sem_r:.2f}',
    ha='center',
    va='bottom',
    fontsize=7,
    color='forestgreen',
)
ax.text(
    1,
    ymin - 0.10 * (ymax - ymin),
    f't(1-samp)={tval:.2f}, p={p_t:.2e}\n'
    f'Wilcoxon={wstat:.2f}, p={p_w:.2e}',
    ha='center',
    va='top',
    fontsize=6.5,
    color='black',
)
ax.set(
    xlim=(0.5, 1.5),
    xticks=[1],
    xticklabels=['Real r'],
    ylabel='Correlation (r)',
    title='Across-sess. r',
)
ax.spines[['top', 'right', 'bottom']].set_visible(False)
plt.tight_layout()
plt.show()
if save:
    for ext in ['.pdf', '.png']:
        fig.savefig(GLM_stem / f'rew_to_run_r_violinplot{ext}', dpi=300, bbox_inches='tight')

low_rates_cells = []
high_rates_cells = []
for cluname in all_cell_rate_dict:
    rates = np.array(all_cell_rate_dict[cluname])
    times = np.array(all_cell_time_dict[cluname])
    curr_low = []
    curr_high = []
    for idx, t in enumerate(times):
        if t < 1.5:
            curr_low.append(rates[idx])
        else:
            curr_high.append(rates[idx])
    low_rates_cells.append(np.mean(curr_low))
    high_rates_cells.append(np.mean(curr_high))

plot_violin_with_scatter(
    low_rates_cells,
    high_rates_cells,
    'forestgreen',
    'darkgreen',
    paired=True,
    xticklabels=['Low', 'High'],
    ylabel='Firing rate (Hz)',
    save=True,
    savepath=GLM_stem / 'rew_to_run_FR_violinplot',
    print_statistics=True,
)
