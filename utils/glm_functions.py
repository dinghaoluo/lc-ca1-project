# -*- coding: utf-8 -*-
'''
Created on Thu Sep 25 17:14:41 2025

behaviour predictors and GLM fitting for the LC analyses

@author: Dinghao Luo
'''

#%% imports
import numpy as np
import statsmodels.api as sm


#%% behaviour predictors
def lick_rate_last5s(lick_times, onset_time, window=5.0):
    '''
    compute lick rate (licks/s) in the last pre-onset window.

    parameters:
    - lick_times: 1d array of lick times (s)
    - onset_time: run onset time (s)
    - window: lookback window (s)

    returns:
    - rate: float (licks/sec)
    '''
    mask = (lick_times >= onset_time - window) & (lick_times < onset_time)
    n_licks = np.sum(mask)
    return n_licks / window

def first_lick_to_reward_last_trial(lick_times_trials, reward_times, ti):
    if ti == 0:
        return np.nan
    try:
        last_first_lick = lick_times_trials[ti-1][0]
    except IndexError:  # if no licks
        return np.nan
    last_rew = reward_times[ti-1]
    if np.isnan(last_rew): return np.nan
    return (last_rew - last_first_lick) / 1000.0  # convert ms to s

def time_since_last_reward(reward_times, onset_time, trial_index):
    last_reward_time = reward_times[trial_index - 1]
    if np.isnan(last_reward_time) or np.isnan(onset_time):
        return np.nan
    else:
        last_reward_time /= 1000.0
    return (onset_time - last_reward_time)

def mean_speed_prev_trial(timestamps_s, speeds_cm_s, run_onsets_s, ti):
    if ti == 0:
        return np.nan
    # trial boundaries defined by successive run onsets; last trial ends at this onset
    t_start = run_onsets_s[ti-1]
    t_end   = run_onsets_s[ti] if ti < len(run_onsets_s) else timestamps_s[-1]
    if not np.isfinite(t_start) or not np.isfinite(t_end) or t_end <= t_start:
        return np.nan
    mask = (timestamps_s >= t_start) & (timestamps_s < t_end)
    if not np.any(mask):
        return np.nan
    return float(np.nanmean(speeds_cm_s[mask]))

def preonset_rate(train, samp_freq=1250, onset_idx=3750, window=(2.5, 1.5)):
    '''
    compute the mean pre-onset firing rate in a configurable window.
    '''
    lo = int(onset_idx - window[0]*samp_freq)
    hi = int(onset_idx - window[1]*samp_freq)
    return float(np.nanmean(train[lo:hi]))


#%% target (run onset rates)
def run_onset_amplitude(spk_rate: np.ndarray, sr: float, onset_idx: int) -> float:
    '''
    compute the mean run-onset response in a symmetric peri-onset window.

    parameters:
    - spk_rate: spike rate vector (hz)
    - sr: sampling rate (hz)
    - onset_idx: sample index of run-onset

    returns:
    - amp: summed spike rate in window (float)
    '''
    half_win = int(0.5 * sr)
    lo = onset_idx - half_win
    hi = onset_idx + half_win
    return float(np.nanmean(spk_rate[lo:hi]))


#%% fit GLM
def fit_glm_log_gaussian(X: np.ndarray, y: np.ndarray, eps: float = 1e-6):
    '''
    fit a Gaussian GLM to a log-transformed target vector.

    parameters:
    - X: design matrix (n_trials, n_features)
    - y: target vector (n_trials,)
    - eps: small constant to avoid log(0)

    returns:
    - result: fitted statsmodels glm result
    '''
    y_log = np.log(y.astype(float) + eps)
    Xc = sm.add_constant(X, has_constant='add')
    fam = sm.families.Gaussian()
    model = sm.GLM(y_log, Xc, family=fam)
    return model.fit()
