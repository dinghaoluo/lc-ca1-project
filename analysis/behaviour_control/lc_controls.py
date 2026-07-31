# -*- coding: utf-8 -*-
'''
Created on Thu Jun  5 15:29:38 2025

controls for LC run-onset peaks

@author: Dinghao Luo
'''

#%% imports
from pathlib import Path
import sys

import numpy as np
import matplotlib.pyplot as plt
import pickle
import pandas as pd
from scipy.stats import sem, wilcoxon, ttest_1samp, linregress

repo_root = Path(__file__).resolve().parents[2]

if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
utils_path = repo_root / 'utils'
if str(utils_path) not in sys.path:
    sys.path.insert(0, str(utils_path))

import rec_list
paths = rec_list.pathLC

import plotting_functions as pf
from common_functions import mpl_formatting
from console_formatting import print_session, print_statistics_section, print_status
import project_paths as pp
mpl_formatting()


#%% paths and parameters
LC_stem = pp.LC_EPHYS_STEM
all_sess_stem = LC_stem / 'all_sessions'
LC_beh_stem = pp.behaviour_experiment_stem('LC')
SPEED_CONTROL_FIGURE_STEM = pp.LC_EPHYS_FIGURES_STEM / 'speed_controls'
SINGLE_SESSION_INIT_SPEED_FR_STEM = SPEED_CONTROL_FIGURE_STEM / 'single_session_init_speed_FR'
SINGLE_SESSION_INIT_ACCEL_FR_STEM = SPEED_CONTROL_FIGURE_STEM / 'single_session_init_accel_FR'

SAMP_FREQ = 1_250
SAMP_FREQ_BEH = 1_000
RUN_ONSET_BIN = SAMP_FREQ * 3
BEF = 1  # s, how much time before run-onset to get
AFT = 4  # same as above
WINDOW_HALF_SIZE = .5

RO_WINDOW = [
    int(RUN_ONSET_BIN - WINDOW_HALF_SIZE * SAMP_FREQ),
    int(RUN_ONSET_BIN + WINDOW_HALF_SIZE * SAMP_FREQ)
]  # window for spike summation, half a sec around run onsets

XAXIS_SPIKE_TIME = np.arange(SAMP_FREQ * (BEF + AFT)) / SAMP_FREQ - 1  # 5 seconds
XAXIS_SPEED_TIME = np.arange(SAMP_FREQ_BEH * 4) / SAMP_FREQ_BEH  # 4 seconds


#%% analysis
for stem in [
        SPEED_CONTROL_FIGURE_STEM,
        SINGLE_SESSION_INIT_SPEED_FR_STEM,
        SINGLE_SESSION_INIT_ACCEL_FR_STEM,
        ]:
    stem.mkdir(parents=True, exist_ok=True)

print('loading data...')
cell_prop_path = LC_stem / 'LC_all_cell_profiles.pkl'
cell_prop = pd.read_pickle(cell_prop_path)

tagged_ro_keys = []
putative_ro_keys = []
for clu in cell_prop.itertuples():
    if clu.identity == 'tagged' and clu.run_onset_peak:
        tagged_ro_keys.append(clu.Index)
    if clu.identity == 'putative' and clu.run_onset_peak:
        putative_ro_keys.append(clu.Index)
ro_keys = tagged_ro_keys + putative_ro_keys

all_high_speed_speed = []
all_low_speed_speed = []
all_high_init_speed_speed = []
all_low_init_speed_speed = []
all_high_accel_speed = []
all_low_accel_speed = []

all_high_speed_curve = []
all_low_speed_curve = []
all_high_init_speed_curve = []
all_low_init_speed_curve = []
all_high_accel_curve = []
all_low_accel_curve = []

all_high_speed_amp = []
all_low_speed_amp = []
all_high_init_speed_amp = []
all_low_init_speed_amp = []
all_high_accel_amp = []
all_low_accel_amp = []

all_session_init_speed_FR_r = []
all_session_init_accel_FR_r = []
all_init_speeds = []
all_init_accel = []
all_init_FR = []

animals = set()
for path in paths:
    recname = Path(path).name
    print_session(recname)

    trains_path = all_sess_stem / recname / f'{recname}_all_trains_run.npy'
    trains = np.load(trains_path, allow_pickle=True).item()

    if not [clu for clu in trains.keys() if clu in ro_keys]:
        print_status('skipped', 'no run-onset cell')
        continue

    anmname = recname.split('-')[0]
    animals.add(anmname)
    with open(LC_beh_stem / f'{recname}.pkl', 'rb') as f:
        beh = pickle.load(f)

    # trial filtering (ignore bad trials & stim trials)
    stim_conds = [t[15] for t in beh['trial_statements']][1:]
    stim_idx = [trial for trial, cond in enumerate(stim_conds) if cond != '0']
    bad_idx = [trial for trial, bad in enumerate(beh['bad_trials'][1:]) if bad]
    valid_trial_idx = [
        trial for trial in np.arange(len(stim_conds))
        if trial not in stim_idx and trial not in bad_idx
    ]

    speed_trials = [[t[1] for t in trial] for trial in beh['speed_times_aligned'][1:]]

    mean_speed_trials = [np.mean(trial) for trial in speed_trials]
    init_speed_trials = [np.mean(trial[:500]) for trial in speed_trials]
    acceleration_trials = [
        np.mean(np.diff(trial[:500])) * 1_000
        for trial in speed_trials
    ]

    high_speed_idx = [
        trial for trial, speed in enumerate(mean_speed_trials)
        if 45 < speed < 55
    ]
    low_speed_idx = [
        trial for trial, speed in enumerate(mean_speed_trials)
        if 35 < speed < 45
    ]
    high_init_speed_idx = [
        trial for trial, speed in enumerate(init_speed_trials)
        if 35 < speed < 45
    ]
    low_init_speed_idx = [
        trial for trial, speed in enumerate(init_speed_trials)
        if 25 < speed < 35
    ]
    high_accel_idx = [
        trial for trial, accel in enumerate(acceleration_trials)
        if accel > 80
    ]
    low_accel_idx = [
        trial for trial, accel in enumerate(acceleration_trials)
        if 60 < accel < 80
    ]

    high_speed_trials = [speed_trials[trial] for trial in valid_trial_idx if trial in high_speed_idx]
    low_speed_trials = [speed_trials[trial] for trial in valid_trial_idx if trial in low_speed_idx]
    high_init_speed_trials = [
        speed_trials[trial] for trial in valid_trial_idx
        if trial in high_init_speed_idx
    ]
    low_init_speed_trials = [
        speed_trials[trial] for trial in valid_trial_idx
        if trial in low_init_speed_idx
    ]
    high_accel_trials = [speed_trials[trial] for trial in valid_trial_idx if trial in high_accel_idx]
    low_accel_trials = [speed_trials[trial] for trial in valid_trial_idx if trial in low_accel_idx]

    if high_speed_trials and low_speed_trials:
        high_speed_lengths = [len(trial) for trial in high_speed_trials]
        padded_high_speed = np.full(
            (len(high_speed_trials), max(high_speed_lengths)), np.nan
        )
        for trial, speed in enumerate(high_speed_trials):
            padded_high_speed[trial, :len(speed)] = speed

        low_speed_lengths = [len(trial) for trial in low_speed_trials]
        padded_low_speed = np.full(
            (len(low_speed_trials), max(low_speed_lengths)), np.nan
        )
        for trial, speed in enumerate(low_speed_trials):
            padded_low_speed[trial, :len(speed)] = speed

        all_high_speed_speed.append(np.nanmean(padded_high_speed, axis=0))
        all_low_speed_speed.append(np.nanmean(padded_low_speed, axis=0))

    if high_init_speed_trials and low_init_speed_trials:
        high_init_speed_lengths = [len(trial) for trial in high_init_speed_trials]
        padded_high_init_speed = np.full(
            (len(high_init_speed_trials), max(high_init_speed_lengths)), np.nan
        )
        for trial, speed in enumerate(high_init_speed_trials):
            padded_high_init_speed[trial, :len(speed)] = speed

        low_init_speed_lengths = [len(trial) for trial in low_init_speed_trials]
        padded_low_init_speed = np.full(
            (len(low_init_speed_trials), max(low_init_speed_lengths)), np.nan
        )
        for trial, speed in enumerate(low_init_speed_trials):
            padded_low_init_speed[trial, :len(speed)] = speed

        all_high_init_speed_speed.append(np.nanmean(padded_high_init_speed, axis=0))
        all_low_init_speed_speed.append(np.nanmean(padded_low_init_speed, axis=0))

    if high_accel_trials and low_accel_trials:
        high_accel_lengths = [len(trial) for trial in high_accel_trials]
        padded_high_accel = np.full(
            (len(high_accel_trials), max(high_accel_lengths)), np.nan
        )
        for trial, speed in enumerate(high_accel_trials):
            padded_high_accel[trial, :len(speed)] = speed

        low_accel_lengths = [len(trial) for trial in low_accel_trials]
        padded_low_accel = np.full(
            (len(low_accel_trials), max(low_accel_lengths)), np.nan
        )
        for trial, speed in enumerate(low_accel_trials):
            padded_low_accel[trial, :len(speed)] = speed

        all_high_accel_speed.append(np.nanmean(padded_high_accel, axis=0))
        all_low_accel_speed.append(np.nanmean(padded_low_accel, axis=0))

    curr_high_speed_curve = []
    curr_low_speed_curve = []
    curr_high_init_speed_curve = []
    curr_low_init_speed_curve = []
    curr_high_accel_curve = []
    curr_low_accel_curve = []

    curr_high_speed_amp = []
    curr_low_speed_amp = []
    curr_high_init_speed_amp = []
    curr_low_init_speed_amp = []
    curr_high_accel_amp = []
    curr_low_accel_amp = []

    for clu in list(trains.keys()):
        if clu in ro_keys:
            curr_train = trains[clu]

            high_speed_cell_curve = [
                curr_train[trial] for trial in valid_trial_idx
                if trial in high_speed_idx
            ]
            low_speed_cell_curve = [
                curr_train[trial] for trial in valid_trial_idx
                if trial in low_speed_idx
            ]
            high_init_speed_cell_curve = [
                curr_train[trial] for trial in valid_trial_idx
                if trial in high_init_speed_idx
            ]
            low_init_speed_cell_curve = [
                curr_train[trial] for trial in valid_trial_idx
                if trial in low_init_speed_idx
            ]
            high_accel_cell_curve = [
                curr_train[trial] for trial in valid_trial_idx
                if trial in high_accel_idx
            ]
            low_accel_cell_curve = [
                curr_train[trial] for trial in valid_trial_idx
                if trial in low_accel_idx
            ]

            high_speed_cell_amp = [
                np.mean(curr_train[trial][RO_WINDOW[0]:RO_WINDOW[1]])
                for trial in valid_trial_idx if trial in high_speed_idx
            ]
            low_speed_cell_amp = [
                np.mean(curr_train[trial][RO_WINDOW[0]:RO_WINDOW[1]])
                for trial in valid_trial_idx if trial in low_speed_idx
            ]
            high_init_speed_cell_amp = [
                np.mean(curr_train[trial][RO_WINDOW[0]:RO_WINDOW[1]])
                for trial in valid_trial_idx if trial in high_init_speed_idx
            ]
            low_init_speed_cell_amp = [
                np.mean(curr_train[trial][RO_WINDOW[0]:RO_WINDOW[1]])
                for trial in valid_trial_idx if trial in low_init_speed_idx
            ]
            high_accel_cell_amp = [
                np.mean(curr_train[trial][RO_WINDOW[0]:RO_WINDOW[1]])
                for trial in valid_trial_idx if trial in high_accel_idx
            ]
            low_accel_cell_amp = [
                np.mean(curr_train[trial][RO_WINDOW[0]:RO_WINDOW[1]])
                for trial in valid_trial_idx if trial in low_accel_idx
            ]

            if high_speed_cell_curve and low_speed_cell_curve:
                curr_high_speed_curve.append(np.mean(high_speed_cell_curve, axis=0))
                curr_low_speed_curve.append(np.mean(low_speed_cell_curve, axis=0))
            if high_init_speed_cell_curve and low_init_speed_cell_curve:
                curr_high_init_speed_curve.append(
                    np.mean(high_init_speed_cell_curve, axis=0)
                )
                curr_low_init_speed_curve.append(
                    np.mean(low_init_speed_cell_curve, axis=0)
                )
            if high_accel_cell_curve and low_accel_cell_curve:
                curr_high_accel_curve.append(np.mean(high_accel_cell_curve, axis=0))
                curr_low_accel_curve.append(np.mean(low_accel_cell_curve, axis=0))

            if high_speed_cell_amp and low_speed_cell_amp:
                curr_high_speed_amp.append(np.mean(high_speed_cell_amp))
                curr_low_speed_amp.append(np.mean(low_speed_cell_amp))
            if high_init_speed_cell_amp and low_init_speed_cell_amp:
                curr_high_init_speed_amp.append(np.mean(high_init_speed_cell_amp))
                curr_low_init_speed_amp.append(np.mean(low_init_speed_cell_amp))
            if high_accel_cell_amp and low_accel_cell_amp:
                curr_high_accel_amp.append(np.mean(high_accel_cell_amp))
                curr_low_accel_amp.append(np.mean(low_accel_cell_amp))

    if curr_high_speed_curve and curr_low_speed_curve:
        all_high_speed_curve.append(np.mean(curr_high_speed_curve, axis=0))
        all_low_speed_curve.append(np.mean(curr_low_speed_curve, axis=0))

    mean_high_speed_amp = np.mean(curr_high_speed_amp)
    mean_low_speed_amp = np.mean(curr_low_speed_amp)
    if not np.isnan(mean_high_speed_amp) and not np.isnan(mean_low_speed_amp):
        all_high_speed_amp.append(np.mean(curr_high_speed_amp))
        all_low_speed_amp.append(np.mean(curr_low_speed_amp))

    if curr_high_init_speed_curve and curr_low_init_speed_curve:
        all_high_init_speed_curve.append(np.mean(curr_high_init_speed_curve, axis=0))
        all_low_init_speed_curve.append(np.mean(curr_low_init_speed_curve, axis=0))

    mean_high_init_speed_amp = np.mean(curr_high_init_speed_amp)
    mean_low_init_speed_amp = np.mean(curr_low_init_speed_amp)
    if not np.isnan(mean_high_init_speed_amp) and not np.isnan(mean_low_init_speed_amp):
        all_high_init_speed_amp.append(np.mean(curr_high_init_speed_amp))
        all_low_init_speed_amp.append(np.mean(curr_low_init_speed_amp))

    if curr_high_accel_curve and curr_low_accel_curve:
        all_high_accel_curve.append(np.mean(curr_high_accel_curve, axis=0))
        all_low_accel_curve.append(np.mean(curr_low_accel_curve, axis=0))

    mean_high_accel_amp = np.mean(curr_high_accel_amp)
    mean_low_accel_amp = np.mean(curr_low_accel_amp)
    if not np.isnan(mean_high_accel_amp) and not np.isnan(mean_low_accel_amp):
        all_high_accel_amp.append(np.mean(curr_high_accel_amp))
        all_low_accel_amp.append(np.mean(curr_low_accel_amp))

    session_mean_FR_per_trial = []
    session_init_speed_per_trial = []
    session_init_accel_per_trial = []

    for trial in valid_trial_idx:
        session_init_speed_per_trial.append(init_speed_trials[trial])
        session_init_accel_per_trial.append(acceleration_trials[trial])

        fr_list = []
        for clu in trains.keys():
            if clu in ro_keys:
                fr_list.append(
                    np.mean(trains[clu][trial][RO_WINDOW[0]:RO_WINDOW[1]])
                )
        if len(fr_list) > 0:
            session_mean_FR_per_trial.append(np.mean(fr_list))

    if len(session_init_speed_per_trial) > 3:
        all_init_speeds.extend(session_init_speed_per_trial)
        all_init_accel.extend(session_init_accel_per_trial)
        all_init_FR.extend(session_mean_FR_per_trial)

        session_init_speed_per_trial = np.array(session_init_speed_per_trial, float)
        session_init_accel_per_trial = np.array(session_init_accel_per_trial, float)
        session_mean_FR_per_trial = np.array(session_mean_FR_per_trial, float)

        slope, intercept, r, p, _ = linregress(session_init_speed_per_trial, session_mean_FR_per_trial)
        all_session_init_speed_FR_r.append(r)
        xfit = np.linspace(session_init_speed_per_trial.min(), session_init_speed_per_trial.max(), 2)
        yfit = intercept + slope * xfit

        fig, ax = plt.subplots(figsize=(1.6, 1.6))
        ax.scatter(session_init_speed_per_trial, session_mean_FR_per_trial,
                   s=12, color='orange', ec='none', alpha=0.7)
        ax.plot(xfit, yfit, color='black', lw=1)
        ax.text(0.05, 0.95,
                f'r = {r:.2f}\np = {p:.3g}',
                transform=ax.transAxes,
                ha='left', va='top', fontsize=7)
        ax.set(xlabel='Init. speed (cm/s)',
               ylabel='Run-onset FR (Hz)',
               title=recname)
        ax.spines[['top', 'right']].set_visible(False)

        for ext in ['.png', '.pdf']:
            fig.savefig(
                SINGLE_SESSION_INIT_SPEED_FR_STEM / f'{recname}_init_speed_FR_corr{ext}',
                dpi=300, bbox_inches='tight'
            )
        plt.close(fig)

        slope, intercept, r, p, _ = linregress(session_init_accel_per_trial, session_mean_FR_per_trial)
        all_session_init_accel_FR_r.append(r)
        xfit = np.linspace(session_init_accel_per_trial.min(), session_init_accel_per_trial.max(), 2)
        yfit = intercept + slope * xfit

        fig, ax = plt.subplots(figsize=(1.6, 1.6))
        ax.scatter(session_init_accel_per_trial, session_mean_FR_per_trial,
                   s=12, color='firebrick', ec='none', alpha=0.7)
        ax.plot(xfit, yfit, color='black', lw=1)
        ax.text(0.05, 0.95,
                f'r = {r:.2f}\np = {p:.3g}',
                transform=ax.transAxes,
                ha='left', va='top', fontsize=7)
        ax.set(xlabel='Init. accel. (cm/s²)',
               ylabel='Run-onset FR (Hz)',
               title=recname)
        ax.spines[['top', 'right']].set_visible(False)

        for ext in ['.png', '.pdf']:
            fig.savefig(
                SINGLE_SESSION_INIT_ACCEL_FR_STEM / f'{recname}_init_accel_FR_corr{ext}',
                dpi=300, bbox_inches='tight'
            )
        plt.close(fig)
    else:
        print_status('skipped', 'not enough trials')

#%% plotting - speed
padded_high_speed = np.full((len(all_high_speed_speed), 4000), np.nan)
for session, speed in enumerate(all_high_speed_speed):
    n_valid = min(len(speed), 4000)
    padded_high_speed[session, :n_valid] = speed[:n_valid]

padded_low_speed = np.full((len(all_low_speed_speed), 4000), np.nan)
for session, speed in enumerate(all_low_speed_speed):
    n_valid = min(len(speed), 4000)
    padded_low_speed[session, :n_valid] = speed[:n_valid]

high_speed_mean = np.mean(padded_high_speed, axis=0)
high_speed_sem = sem(padded_high_speed, axis=0)
low_speed_mean = np.mean(padded_low_speed, axis=0)
low_speed_sem = sem(padded_low_speed, axis=0)

fig, ax = plt.subplots(figsize=(1.8, 1.4))
high_speed_line, = ax.plot(XAXIS_SPEED_TIME, high_speed_mean, color='orange')
ax.fill_between(
    XAXIS_SPEED_TIME,
    high_speed_mean + high_speed_sem,
    high_speed_mean - high_speed_sem,
    alpha=.3,
    color='orange',
    edgecolor='none',
)
low_speed_line, = ax.plot(XAXIS_SPEED_TIME, low_speed_mean, color='navajowhite')
ax.fill_between(
    XAXIS_SPEED_TIME,
    low_speed_mean + low_speed_sem,
    low_speed_mean - low_speed_sem,
    alpha=.3,
    color='navajowhite',
    edgecolor='none',
)

ax.set(
    xlabel='Time from run onset (s)',
    xticks=[0, 2, 4],
    ylabel='Speed (cm/s)',
    ylim=(0, 75),
    title='High v low speed',
)
ax.legend([high_speed_line, low_speed_line], ['High', 'Low'], frameon=False)
for spine in ['top', 'right']:
    ax.spines[spine].set_visible(False)

for ext in ['.png', '.pdf']:
    fig.savefig(
        SPEED_CONTROL_FIGURE_STEM / f'high_low_speed_speed{ext}',
        dpi=300,
        bbox_inches='tight',
    )

padded_high_init_speed = np.full(
    (len(all_high_init_speed_speed), 4000), np.nan
)
for session, speed in enumerate(all_high_init_speed_speed):
    n_valid = min(len(speed), 4000)
    padded_high_init_speed[session, :n_valid] = speed[:n_valid]

padded_low_init_speed = np.full(
    (len(all_low_init_speed_speed), 4000), np.nan
)
for session, speed in enumerate(all_low_init_speed_speed):
    n_valid = min(len(speed), 4000)
    padded_low_init_speed[session, :n_valid] = speed[:n_valid]

high_init_speed_mean = np.mean(padded_high_init_speed, axis=0)
high_init_speed_sem = sem(padded_high_init_speed, axis=0)
low_init_speed_mean = np.mean(padded_low_init_speed, axis=0)
low_init_speed_sem = sem(padded_low_init_speed, axis=0)

fig, ax = plt.subplots(figsize=(1.8, 1.4))
high_init_speed_line, = ax.plot(
    XAXIS_SPEED_TIME, high_init_speed_mean, color='orange'
)
ax.fill_between(
    XAXIS_SPEED_TIME,
    high_init_speed_mean + high_init_speed_sem,
    high_init_speed_mean - high_init_speed_sem,
    alpha=.3,
    color='orange',
    edgecolor='none',
)
low_init_speed_line, = ax.plot(
    XAXIS_SPEED_TIME, low_init_speed_mean, color='navajowhite'
)
ax.fill_between(
    XAXIS_SPEED_TIME,
    low_init_speed_mean + low_init_speed_sem,
    low_init_speed_mean - low_init_speed_sem,
    alpha=.3,
    color='navajowhite',
    edgecolor='none',
)

ax.set(
    xlabel='Time from run onset (s)',
    xticks=[0, 2, 4],
    ylabel='Speed (cm/s)',
    ylim=(0, 75),
    title='High v low init. speed',
)
ax.legend([high_init_speed_line, low_init_speed_line], ['High', 'Low'], frameon=False)
for spine in ['top', 'right']:
    ax.spines[spine].set_visible(False)

for ext in ['.png', '.pdf']:
    fig.savefig(
        SPEED_CONTROL_FIGURE_STEM / f'high_low_init_speed_speed{ext}',
        dpi=300,
        bbox_inches='tight',
    )

padded_high_accel = np.full((len(all_high_accel_speed), 4000), np.nan)
for session, speed in enumerate(all_high_accel_speed):
    n_valid = min(len(speed), 4000)
    padded_high_accel[session, :n_valid] = speed[:n_valid]

padded_low_accel = np.full((len(all_low_accel_speed), 4000), np.nan)
for session, speed in enumerate(all_low_accel_speed):
    n_valid = min(len(speed), 4000)
    padded_low_accel[session, :n_valid] = speed[:n_valid]

high_accel_mean = np.mean(padded_high_accel, axis=0)
high_accel_sem = sem(padded_high_accel, axis=0)
low_accel_mean = np.mean(padded_low_accel, axis=0)
low_accel_sem = sem(padded_low_accel, axis=0)

fig, ax = plt.subplots(figsize=(1.8, 1.4))
high_accel_line, = ax.plot(XAXIS_SPEED_TIME, high_accel_mean, color='firebrick')
ax.fill_between(
    XAXIS_SPEED_TIME,
    high_accel_mean + high_accel_sem,
    high_accel_mean - high_accel_sem,
    alpha=.25,
    color='firebrick',
    edgecolor='none',
)
low_accel_line, = ax.plot(XAXIS_SPEED_TIME, low_accel_mean, color='lightcoral')
ax.fill_between(
    XAXIS_SPEED_TIME,
    low_accel_mean + low_accel_sem,
    low_accel_mean - low_accel_sem,
    alpha=.25,
    color='lightcoral',
    edgecolor='none',
)

ax.set(
    xlabel='Time from run onset (s)',
    xticks=[0, 2, 4],
    ylabel='Speed (cm/s)',
    ylim=(0, 75),
    title='High v low accel',
)
ax.legend([high_accel_line, low_accel_line], ['High', 'Low'], frameon=False)
for spine in ['top', 'right']:
    ax.spines[spine].set_visible(False)

for ext in ['.png', '.pdf']:
    fig.savefig(
        SPEED_CONTROL_FIGURE_STEM / f'high_low_accel_speed{ext}',
        dpi=300,
        bbox_inches='tight',
    )

#%% plotting - spike curves
high_speed_mean = np.mean(all_high_speed_curve, axis=0)[1250 * 2:1250 * 7]
high_speed_sem = sem(all_high_speed_curve, axis=0)[1250 * 2:1250 * 7]
low_speed_mean = np.mean(all_low_speed_curve, axis=0)[1250 * 2:1250 * 7]
low_speed_sem = sem(all_low_speed_curve, axis=0)[1250 * 2:1250 * 7]

fig, ax = plt.subplots(figsize=(1.8, 1.4))
high_speed_line, = ax.plot(XAXIS_SPIKE_TIME, high_speed_mean, color='orange')
ax.fill_between(
    XAXIS_SPIKE_TIME,
    high_speed_mean + high_speed_sem,
    high_speed_mean - high_speed_sem,
    alpha=.25,
    color='orange',
    edgecolor='none',
)
low_speed_line, = ax.plot(XAXIS_SPIKE_TIME, low_speed_mean, color='navajowhite')
ax.fill_between(
    XAXIS_SPIKE_TIME,
    low_speed_mean + low_speed_sem,
    low_speed_mean - low_speed_sem,
    alpha=.25,
    color='navajowhite',
    edgecolor='none',
)

ax.set(
    xlabel='Time from run onset (s)',
    ylabel='Firing rate (Hz)',
    ylim=(1.5, 5.7),
    title='High v low speed',
)
ax.legend([high_speed_line, low_speed_line], ['High', 'Low'], frameon=False)
for spine in ['top', 'right']:
    ax.spines[spine].set_visible(False)

for ext in ['.png', '.pdf']:
    fig.savefig(
        SPEED_CONTROL_FIGURE_STEM / f'high_low_speed_spike_curves{ext}',
        dpi=300,
        bbox_inches='tight',
    )

high_init_speed_mean = np.mean(all_high_init_speed_curve, axis=0)[1250 * 2:1250 * 7]
high_init_speed_sem = sem(all_high_init_speed_curve, axis=0)[1250 * 2:1250 * 7]
low_init_speed_mean = np.mean(all_low_init_speed_curve, axis=0)[1250 * 2:1250 * 7]
low_init_speed_sem = sem(all_low_init_speed_curve, axis=0)[1250 * 2:1250 * 7]

fig, ax = plt.subplots(figsize=(1.8, 1.4))
high_init_speed_line, = ax.plot(
    XAXIS_SPIKE_TIME, high_init_speed_mean, color='orange'
)
ax.fill_between(
    XAXIS_SPIKE_TIME,
    high_init_speed_mean + high_init_speed_sem,
    high_init_speed_mean - high_init_speed_sem,
    alpha=.25,
    color='orange',
    edgecolor='none',
)
low_init_speed_line, = ax.plot(
    XAXIS_SPIKE_TIME, low_init_speed_mean, color='navajowhite'
)
ax.fill_between(
    XAXIS_SPIKE_TIME,
    low_init_speed_mean + low_init_speed_sem,
    low_init_speed_mean - low_init_speed_sem,
    alpha=.25,
    color='navajowhite',
    edgecolor='none',
)

ax.set(
    xlabel='Time from run onset (s)',
    ylabel='Firing rate (Hz)',
    ylim=(1.5, 5.7),
    title='High v low init. speed',
)
ax.legend([high_init_speed_line, low_init_speed_line], ['High', 'Low'], frameon=False)
for spine in ['top', 'right']:
    ax.spines[spine].set_visible(False)

for ext in ['.png', '.pdf']:
    fig.savefig(
        SPEED_CONTROL_FIGURE_STEM / f'high_low_init_speed_spike_curves{ext}',
        dpi=300,
        bbox_inches='tight',
    )

high_accel_mean = np.mean(all_high_accel_curve, axis=0)[1250 * 2:1250 * 7]
high_accel_sem = sem(all_high_accel_curve, axis=0)[1250 * 2:1250 * 7]
low_accel_mean = np.mean(all_low_accel_curve, axis=0)[1250 * 2:1250 * 7]
low_accel_sem = sem(all_low_accel_curve, axis=0)[1250 * 2:1250 * 7]

fig, ax = plt.subplots(figsize=(1.8, 1.4))
high_accel_line, = ax.plot(XAXIS_SPIKE_TIME, high_accel_mean, color='firebrick')
ax.fill_between(
    XAXIS_SPIKE_TIME,
    high_accel_mean + high_accel_sem,
    high_accel_mean - high_accel_sem,
    alpha=.25,
    color='firebrick',
    edgecolor='none',
)
low_accel_line, = ax.plot(XAXIS_SPIKE_TIME, low_accel_mean, color='lightcoral')
ax.fill_between(
    XAXIS_SPIKE_TIME,
    low_accel_mean + low_accel_sem,
    low_accel_mean - low_accel_sem,
    alpha=.25,
    color='lightcoral',
    edgecolor='none',
)

ax.set(
    xlabel='Time from run onset (s)',
    ylabel='Firing rate (Hz)',
    ylim=(1.5, 5.7),
    title='High v low accel',
)
ax.legend([high_accel_line, low_accel_line], ['High', 'Low'], frameon=False)
for spine in ['top', 'right']:
    ax.spines[spine].set_visible(False)

for ext in ['.png', '.pdf']:
    fig.savefig(
        SPEED_CONTROL_FIGURE_STEM / f'high_low_accel_spike_curves{ext}',
        dpi=300,
        bbox_inches='tight',
    )

#%% plotting - amp
pf.plot_violin_with_scatter(all_low_speed_amp, all_high_speed_amp,
                            'navajowhite', 'orange',
                            xticklabels=['Low speed', 'High speed'],
                            ylabel='Firing rate (Hz)',
                            ylim=(0, 10),
                            print_statistics=True,
                            save=True,
                            savepath=SPEED_CONTROL_FIGURE_STEM / 'speed_35_45_55',
                            show=False,
                            close=True
                            )

pf.plot_violin_with_scatter(all_low_init_speed_amp, all_high_init_speed_amp,
                            'navajowhite', 'orange',
                            xticklabels=['Low init. speed', 'High init. speed'],
                            ylabel='Firing rate (Hz)',
                            ylim=(0, 10),
                            save=True,
                            savepath=SPEED_CONTROL_FIGURE_STEM / 'init_speed_25_35_45',
                            show=False,
                            close=True
                            )

pf.plot_violin_with_scatter(all_low_accel_amp, all_high_accel_amp,
                            'lightcoral', 'firebrick',
                            xticklabels=['Low init. accel.', 'High init. accel.'],
                            ylabel='Firing rate (Hz)',
                            ylim=(0, 10),
                            save=True,
                            savepath=SPEED_CONTROL_FIGURE_STEM / 'init_accel_60_80',
                            show=False,
                            close=True
                            )

#%% corr between init speed and FR
speed_rvals = np.array(all_session_init_speed_FR_r, float)
speed_rvals = speed_rvals[~np.isnan(speed_rvals)]
n_speed_sess = len(speed_rvals)

median_speed_r = np.median(speed_rvals)
mean_speed_r = np.mean(speed_rvals)
sem_speed_r = np.std(speed_rvals, ddof=1) / np.sqrt(n_speed_sess)
q25_speed_r, q75_speed_r = np.percentile(speed_rvals, [25, 75])
iqr_speed_r = q75_speed_r - q25_speed_r

w_speed_stat, p_speed_w = wilcoxon(speed_rvals, alternative='two-sided')
t_speed_stat, p_speed_t = ttest_1samp(speed_rvals, popmean=0)

print_statistics_section()
print(f'n_sessions = {n_speed_sess}')
print(f'median r = {median_speed_r:.3f}')
print(f'mean r +/- SEM = {mean_speed_r:.3f} +/- {sem_speed_r:.3f}')
print(f'IQR = [{q25_speed_r:.3f}, {q75_speed_r:.3f}] (IQR = {iqr_speed_r:.3f})')
print(f'wilcoxon vs 0: W = {w_speed_stat:.3f}, p = {p_speed_w:.3g}')
print(f't-test vs 0: T = {t_speed_stat:.3f}, p = {p_speed_t:.3g}')

fig, ax = plt.subplots(figsize=(2.0, 2.2))
parts = ax.violinplot(
    speed_rvals, positions=[1], showmeans=False, showmedians=True, showextrema=False
)

for pc in parts['bodies']:
    pc.set_facecolor('orange')
    pc.set_edgecolor('none')
    pc.set_alpha(0.35)

parts['cmedians'].set_color('k')
parts['cmedians'].set_linewidth(1.2)
ax.scatter(
    np.ones(n_speed_sess), speed_rvals,
    s=12, color='orange', ec='none', alpha=0.55, zorder=3
)
ax.axhline(0, color='gray', lw=1, ls='--')

ax.text(
    1.35,
    np.max(speed_rvals),
    f'Median = {median_speed_r:.2f}\n'
    f'IQR = [{q25_speed_r:.2f}, {q75_speed_r:.2f}]\n'
    f'{mean_speed_r:.2f} +/- {sem_speed_r:.2f}\n'
    f'Wilc {p_speed_w:.2e}\n'
    f'Ttest {p_speed_t:.2e}',
    ha='left',
    va='top',
    fontsize=7,
    color='forestgreen',
)

ax.set(
    xlim=(0.5, 1.5),
    xticks=[1],
    xticklabels=['corr(init. speed, RO FR)'],
    ylim=(-0.5, 0.5),
    ylabel='Correlation (r)',
    title='Across-session corr.',
)
ax.spines[['top', 'right', 'bottom']].set_visible(False)
plt.tight_layout()

#%% corr between init accel and FR
accel_rvals = np.array(all_session_init_accel_FR_r, float)
accel_rvals = accel_rvals[~np.isnan(accel_rvals)]
n_accel_sess = len(accel_rvals)

median_accel_r = np.median(accel_rvals)
mean_accel_r = np.mean(accel_rvals)
sem_accel_r = np.std(accel_rvals, ddof=1) / np.sqrt(n_accel_sess)
q25_accel_r, q75_accel_r = np.percentile(accel_rvals, [25, 75])
iqr_accel_r = q75_accel_r - q25_accel_r

w_accel_stat, p_accel_w = wilcoxon(accel_rvals, alternative='two-sided')
t_accel_stat, p_accel_t = ttest_1samp(accel_rvals, popmean=0)

print_statistics_section()
print(f'n_sessions = {n_accel_sess}')
print(f'median r = {median_accel_r:.3f}')
print(f'mean r +/- SEM = {mean_accel_r:.3f} +/- {sem_accel_r:.3f}')
print(f'IQR = [{q25_accel_r:.3f}, {q75_accel_r:.3f}] (IQR = {iqr_accel_r:.3f})')
print(f'wilcoxon vs 0: W = {w_accel_stat:.3f}, p = {p_accel_w:.3g}')
print(f't-test vs 0: T = {t_accel_stat:.3f}, p = {p_accel_t:.3g}')

fig, ax = plt.subplots(figsize=(2.0, 2.2))
parts = ax.violinplot(
    accel_rvals, positions=[1], showmeans=False, showmedians=True, showextrema=False
)

for pc in parts['bodies']:
    pc.set_facecolor('darkred')
    pc.set_edgecolor('none')
    pc.set_alpha(0.35)

parts['cmedians'].set_color('k')
parts['cmedians'].set_linewidth(1.2)
ax.scatter(
    np.ones(n_accel_sess), accel_rvals,
    s=12, color='darkred', ec='none', alpha=0.55, zorder=3
)
ax.axhline(0, color='gray', lw=1, ls='--')

ax.text(
    1.35,
    np.max(accel_rvals),
    f'Median = {median_accel_r:.2f}\n'
    f'IQR = [{q25_accel_r:.2f}, {q75_accel_r:.2f}]\n'
    f'{mean_accel_r:.2f} +/- {sem_accel_r:.2f}\n'
    f'Wilc {p_accel_w:.2e}\n'
    f'Ttest {p_accel_t:.2e}',
    ha='left',
    va='top',
    fontsize=7,
    color='forestgreen',
)

ax.set(
    xlim=(0.5, 1.5),
    xticks=[1],
    xticklabels=['corr(init. accel, RO FR)'],
    ylim=(-0.5, 0.5),
    ylabel='Correlation (r)',
    title='Across-session corr.',
)
ax.spines[['top', 'right', 'bottom']].set_visible(False)
plt.tight_layout()
