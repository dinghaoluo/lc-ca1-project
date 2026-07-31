'''
Created on Fri Nov 21 10:04:59 2025
Modified on 1 July 2026

estimate the first point of stimulation effects for CA1 cells
add a finite-grid Bayesian latency estimate to the older four-of-five rule

@author: Dinghao Luo
'''

#%% imports
from pathlib import Path
import sys
import pickle
import argparse

import numpy as np
import pandas as pd
import scipy.io as sio
from scipy.stats import sem, ranksums, ks_2samp
from scipy.special import logsumexp, log_ndtr
import matplotlib.pyplot as plt
from tqdm import tqdm

repo_root = Path(__file__).resolve().parents[3]
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))

from behaviour_functions import process_txt, detect_run_onsets_teensy
from common_functions import mpl_formatting
from console_formatting import print_files_saved, print_session, print_statistics_section
import project_paths as pp
mpl_formatting()

import rec_list
bad_behaviour_sessions = {
    Path(path).name for path in rec_list.pathHPCbadbeh
}
paths = [
    path for path in rec_list.pathHPCLCopt
    if Path(path).name not in bad_behaviour_sessions
]


#%% command-line options
parser = argparse.ArgumentParser(
    description='Estimate LC-stimulation CA1 effect latencies.')
parser.add_argument(
    '--single-cell-plots',
    action='store_true',
    help='export per-cell mean and trial heatmap figures during processing',
)
args = parser.parse_args()
SAVE_SINGLE_CELL_PLOTS = args.single_cell_plots


#%% path stems
all_sess_stem    = pp.HPC_EPHYS_STEM / 'all_sessions'
mice_exp_stem    = pp.MICEEXP_ROOT
stim_effect_stem = pp.LC_STIM_CA1_EFFECT_FIGURES_STEM
single_cell_stem = stim_effect_stem / 'single_cells'
data_output_stem = repo_root / 'data' / 'analysis' / 'hpc' / 'stim_ctrl'
data_output_path = data_output_stem / 'lc_stim_ca1_effect_latency.pkl'
event_time_path  = data_output_stem / 'lc_stim_ca1_event_time_effects.csv'
posterior_path   = data_output_stem / 'lc_stim_ca1_effect_latency_bayesian.csv'

data_output_stem.mkdir(parents=True, exist_ok=True)
if SAVE_SINGLE_CELL_PLOTS:
    single_cell_stem.mkdir(parents=True, exist_ok=True)


#%% parameters
SAMP_FREQ = 1250  # Hz

BEF        = 1  # in s
AFT        = 4
MAX_LENGTH = (BEF + AFT) * SAMP_FREQ

XAXIS = np.arange(MAX_LENGTH) / SAMP_FREQ - BEF

OF_CONSTANT = (2**32-1)/1000  # overflow constant; this can only be hard-coded

# for single cell stats
ALPHA            = 0.05
BIN_MS           = 50  # ms
BIN_SIZE         = int((BIN_MS/1000) * SAMP_FREQ)
N_BINS           = MAX_LENGTH // BIN_SIZE
BIN_CENTRES      = (np.arange(N_BINS) * BIN_SIZE + BIN_SIZE/2) / SAMP_FREQ - BEF
SIGNAL_START_BIN = int((BEF * SAMP_FREQ) / BIN_SIZE)

CTRL_COLOR       = 'grey'
STIM_COLOR       = 'royalblue'

# for the probabilistic latency estimate
NO_RESPONSE_PRIOR     = .5
AMPLITUDE_PRIOR_SD_HZ = 5
MIN_EFFECT_SE_HZ      = .001
MIN_RESPONSE_BINS     = 3  # 300 ms in the 100-ms event-time table

#%% latency plots and time correction
def _plot_single_cell_mean(cluname, run_mean, run_sem, stim_mean, stim_sem,
                           effect_sign):
    fig, ax = plt.subplots(figsize=(3.2, 2.8))

    ax.plot(XAXIS, run_mean, c=CTRL_COLOR)
    ax.fill_between(XAXIS, run_mean+run_sem,
                    run_mean-run_sem,
                    color=CTRL_COLOR, edgecolor='none', alpha=.5)
    ax.plot(XAXIS, stim_mean, c=STIM_COLOR)
    ax.fill_between(XAXIS, stim_mean+stim_sem,
                    stim_mean-stim_sem,
                    color=STIM_COLOR, edgecolor='none', alpha=.5)

    ax.set(xlabel='Time from run onset (s)',
           ylabel='Firing rate (Hz)',
           title=cluname)

    ymin, ymax = ax.get_ylim()
    marker_bottom = ymin + 0.90*(ymax - ymin)
    marker_top    = ymin + 0.98*(ymax - ymin)
    for bi in range(N_BINS):
        x = BIN_CENTRES[bi]
        if x < 0:
            continue
        if effect_sign[bi] == 1:
            ax.plot([x, x], [marker_bottom, marker_top], color='red', lw=0.8)
        elif effect_sign[bi] == -1:
            ax.plot([x, x], [marker_bottom, marker_top], color='blue', lw=0.8)

    fig.tight_layout()
    fig.savefig(single_cell_stem / f'{cluname}.png',
                dpi=300, bbox_inches='tight')
    plt.close(fig)

def _plot_single_cell_trials(cluname, run_aligned, stim_aligned,
                             run_onset_online_converted,
                             stim_times_converted):
    vmin = np.nanmin([run_aligned.min(),  stim_aligned.min()])
    vmax = np.nanmax([run_aligned.max(),  stim_aligned.max()])

    fig, axs = plt.subplots(2, 1, figsize=(3, 4), sharex=True)

    axs[0].imshow(run_aligned, aspect='auto', interpolation='none',
                  extent=[-1, 4, len(run_onset_online_converted)+1, 1],
                  vmin=vmin, vmax=vmax)

    axs[1].imshow(stim_aligned, aspect='auto', interpolation='none',
                  extent=[-1, 4, len(stim_times_converted)+1, 1],
                  vmin=vmin, vmax=vmax)

    axs[0].set(title='Run-aligned')
    axs[1].set(xlabel='Time from run/stim. onset (s)',
               title='Stim.-aligned')

    for i in range(2):
        axs[i].set(ylabel='Trial #')

    fig.suptitle(cluname)
    fig.tight_layout()

    fig.savefig(single_cell_stem / f'{cluname}_single_trial.png',
                dpi=300, bbox_inches='tight')
    plt.close(fig)

# correct a single overflow event in monotonic teensy-style timestamps
def correct_teensy_overflow(times, of_constant):
    times = np.array(times, dtype=float)
    overflow = np.flatnonzero(np.diff(times) < 0)
    correction_time = overflow[0] + 1 if len(overflow) else None
    if correction_time is not None:
        times[correction_time:] += of_constant
    return times, correction_time

#%% load data
print('loading dataframes...')
cell_profiles = pd.read_pickle(
    pp.HPC_EPHYS_STEM / 'hpc_all_profiles.pkl')
df_pyr = cell_profiles[cell_profiles['cell_identity'] == 'pyr']


#%% main (HPCLC)
# sig_dict is cluname: first_sustained_time, first_sustained_sign
# first_sustained_time is the first time bin (in seconds) where sustained
#   activation/inhibition is detected
# first_sustained_sign is the sign of the effect
sig_dict = {}
session_records = []
cell_records = []

for path in paths:
    recname = Path(path).name
    print_session(recname)

    txtpath = mice_exp_stem / f'ANMD{recname[1:4]}r' / recname[:-3] / recname / f'{recname}T.txt'
    beh = process_txt(txtpath)

    ## ---- conversion fit ---- ##
    # we need to determine a linear mapping between teensy time and spike time
    # using cue time to align due to least variability
    rec_stem = mice_exp_stem / f'ANMD{recname[1:5]}' / recname[:14] / recname
    aligned_cue_path = rec_stem / f'{recname}_DataStructure_mazeSection1_TrialType1_alignCue_msess1.mat'
    aligned_cue = sio.loadmat(aligned_cue_path)['trialsCue'][0][0]

    cue_spike_time  = aligned_cue['startLfpInd'][0]
    cue_teensy_time = [trial[0][0] if trial else np.nan for trial in beh['movie_times']]
    # correct teensy overflow before fitting the clock conversion
    cue_teensy_time, correction_time = correct_teensy_overflow(
        cue_teensy_time, OF_CONSTANT)  # we note the correction time

    # linear map
    a, b = np.polyfit(cue_teensy_time, cue_spike_time, 1)

    ## ---- conversion fit ends ---- ##

    # get stim and ctrl idx
    stim_conds = [t[15] for t in beh['trial_statements']][1:]
    stim_idx = [trial for trial, cond in enumerate(stim_conds) if cond != '0']
    ctrl_idx = [trial+2 for trial in stim_idx]

    # get online run onsets
    speed_times = beh['speed_times'][1:]
    run_onset_online = detect_run_onsets_teensy(speed_times)
    run_onset_online = [run for trial, run in enumerate(run_onset_online)
                        if trial in ctrl_idx]
    run_onset_online, _ = correct_teensy_overflow(
        run_onset_online, OF_CONSTANT)

    # get stim times and truncate pulses
    pulse_times = np.array(beh['pulse_times'])
    diffs = np.diff(pulse_times)
    split_idx = np.where(diffs >= 1000)[0] + 1
    pulse_trains = np.split(pulse_times, split_idx)
    stim_times = [pulse_train[0] for pulse_train in pulse_trains]  # actual stim times
    stim_times, _ = correct_teensy_overflow(stim_times, OF_CONSTANT)

    # get spike maps
    spike_map_path = all_sess_stem / recname / f'{recname}_smoothed_spike_map.npy'
    spike_maps = np.load(spike_map_path, allow_pickle=True)
    max_time = len(spike_maps[0])

    ## ---- conversion (teensy-time to spike-time) ---- ##
    # IMPORTANT: if cue times needed to be corrected for overflow, we use that
    #    correction_time to determine if the other variables need to be forced
    #    to correct
    if correction_time is not None:  # correction for cue times; otherwise move forwards as normal
        if correction_time <= stim_idx[0]:  # earlier than first stim
            print('force correction for stim. times.')
            stim_times = [t + OF_CONSTANT for t in stim_times]
        if correction_time <= ctrl_idx[0]:  # earlier than first run
            print('force correction for run onsets')
            run_onset_online = [t + OF_CONSTANT for t in run_onset_online]
    # we only need to do this if correction_time is BEFORE stim. or run onsets
    # if AFTER, then we don't need to worry

    # actual conversions
    stim_times_converted       = [int(a*t + b) for t in stim_times]
    run_onset_online_converted = [int(a*t + b) for t in run_onset_online]
    ## ---- conversion ends ---- ##

    # ignore if not enough BEF or AFT
    stim_times_converted       = [stim for stim in stim_times_converted
                                  if BEF*SAMP_FREQ <= stim <= max_time-AFT*SAMP_FREQ]
    run_onset_online_converted = [run for run in run_onset_online_converted
                                  if BEF*SAMP_FREQ <= run <= max_time-AFT*SAMP_FREQ]

    # curr session df
    curr_df_pyr = df_pyr[df_pyr['recname'] == recname]

    session_records.append({
        'recname': recname,
        'source_path': str(path),
        'teensy_to_spike_slope': a,
        'teensy_to_spike_intercept': b,
        'cue_correction_time': correction_time,
        'stim_idx': stim_idx,
        'ctrl_idx': ctrl_idx,
        'stim_times_converted': stim_times_converted,
        'run_onset_online_converted': run_onset_online_converted,
        'n_stim_events': len(stim_times_converted),
        'n_run_events': len(run_onset_online_converted),
        'n_pyr_cells': len(curr_df_pyr),
    })

    ## ---- processing ---- ##
    for idx in tqdm(curr_df_pyr.index, desc='processing', total=len(curr_df_pyr)):
        cluname = idx

        clu_idx = int(cluname.split('clu')[1].split(' ')[0]) - 2  # actual index for retrieval
        spike_map = spike_maps[clu_idx, :]

        ## ---- alignment ---- ##
        run_aligned = np.zeros((len(run_onset_online_converted), MAX_LENGTH))
        for trial, run in enumerate(run_onset_online_converted):
            run_aligned[trial, :] = spike_map[run - BEF*SAMP_FREQ : run + AFT*SAMP_FREQ]

        stim_aligned = np.zeros((len(stim_times_converted), MAX_LENGTH))
        for trial, stim in enumerate(stim_times_converted):
            stim_aligned[trial, :] = spike_map[stim - BEF*SAMP_FREQ : stim + AFT*SAMP_FREQ]
        ## ---- alignment ends ---- ##

        run_mean  = np.mean(run_aligned, axis=0)
        run_sem   = sem(run_aligned, axis=0)
        stim_mean = np.mean(stim_aligned, axis=0)
        stim_sem  = sem(stim_aligned, axis=0)

        ## ---- effect test ---- ##
        pvals = np.ones(N_BINS)
        effect_sign = np.zeros(N_BINS, dtype=int)   # -1, 0, +1

        for bi in range(N_BINS):
            start = bi * BIN_SIZE
            end   = start + BIN_SIZE

            ctrl_bin = np.mean(run_aligned[:,  start:end], axis=1)
            stim_bin = np.mean(stim_aligned[:, start:end], axis=1)

            _, p_val = ranksums(ctrl_bin, stim_bin)
            pvals[bi] = p_val

            if p_val < ALPHA:
                if np.mean(stim_bin) > np.mean(ctrl_bin):
                    effect_sign[bi] = 1    # stim > ctrl
                else:
                    effect_sign[bi] = -1   # stim < ctrl

        first_sustained_time = None
        first_sustained_sign = 0   # +1 or -1

        for bi in range(SIGNAL_START_BIN, N_BINS - 4):
            window = effect_sign[bi : bi+5]

            # activation (stim > ctrl)
            if np.sum(window == 1) >= 4:
                first_sustained_time = float(BIN_CENTRES[bi])
                first_sustained_sign = 1
                break

            # suppression (stim < ctrl)
            if np.sum(window == -1) >= 4:
                first_sustained_time = float(BIN_CENTRES[bi])
                first_sustained_sign = -1
                break

        sig_dict[cluname] = [first_sustained_time, first_sustained_sign]
        cell_records.append({
            'recname': recname,
            'cluname': cluname,
            'clu_idx': clu_idx,
            'run_mean_profile': run_mean,
            'run_sem_profile': run_sem,
            'stim_mean_profile': stim_mean,
            'stim_sem_profile': stim_sem,
            'bin_pvals': pvals,
            'effect_sign': effect_sign,
            'first_sustained_time': first_sustained_time,
            'first_sustained_sign': first_sustained_sign,
            'n_run_trials': run_aligned.shape[0],
            'n_stim_trials': stim_aligned.shape[0],
        })
        ## ---- effect test ends ---- ##

        if SAVE_SINGLE_CELL_PLOTS:
            _plot_single_cell_mean(cluname, run_mean, run_sem, stim_mean,
                                   stim_sem, effect_sign)
            _plot_single_cell_trials(cluname, run_aligned, stim_aligned,
                                     run_onset_online_converted,
                                     stim_times_converted)
    ## ---- processing ends ---- ##


#%% summary
onset_times_act = np.asarray([clu[0] for clu in sig_dict.values() if clu[1] == 1])
onset_times_inh = np.asarray([clu[0] for clu in sig_dict.values() if clu[1] == -1])

onset_times_act = onset_times_act[~np.isnan(onset_times_act)]
onset_times_inh = onset_times_inh[~np.isnan(onset_times_inh)]

n_act = onset_times_act.size
n_inh = onset_times_inh.size

# medians
act_median = np.median(onset_times_act)
inh_median = np.median(onset_times_inh)

# MAD / sqrt(n)
act_mad = np.median(np.abs(onset_times_act - act_median))
inh_mad = np.median(np.abs(onset_times_inh - inh_median))

act_mad_sem = act_mad / np.sqrt(n_act)
inh_mad_sem = inh_mad / np.sqrt(n_inh)

# IQR
act_q25, act_q75 = np.percentile(onset_times_act, [25, 75])
inh_q25, inh_q75 = np.percentile(onset_times_inh, [25, 75])

# printout
print_statistics_section()
print('descriptive four-of-five latency summary:')

print(
    f'activation: n={n_act}, '
    f'median={act_median:.4g}s, '
    f'MAD/sqrt(n)={act_mad_sem:.4g}s, '
    f'IQR=[{act_q25:.4g}, {act_q75:.4g}]s'
)

print(
    f'inhibition: n={n_inh}, '
    f'median={inh_median:.4g}s, '
    f'MAD/sqrt(n)={inh_mad_sem:.4g}s, '
    f'IQR=[{inh_q25:.4g}, {inh_q75:.4g}]s'
)

#%% summary stats
ks_stat, ks_p = ks_2samp(onset_times_act, onset_times_inh)
print(f'ks-test: D={ks_stat:.4g}, p={ks_p:.2e}')

intermediate_payload = {
    'metadata': {
        'source_script': Path(__file__).name,
        'sample_frequency_hz': SAMP_FREQ,
        'before_seconds': BEF,
        'after_seconds': AFT,
        'max_length_samples': MAX_LENGTH,
        'bin_ms': BIN_MS,
        'bin_size_samples': BIN_SIZE,
        'alpha': ALPHA,
        'sustained_window_bins': 5,
        'sustained_min_significant_bins': 4,
        'xaxis_seconds': XAXIS,
        'bin_centres_seconds': BIN_CENTRES,
        'paths_attr': 'pathHPCLCopt',
    },
    'session_records': session_records,
    'cell_records': cell_records,
    'sig_dict': sig_dict,
    'latency_summary': {
        'method': 'rank-sum tests followed by the first four-of-five window',
        'activation_times': onset_times_act,
        'inhibition_times': onset_times_inh,
        'n_activation': n_act,
        'n_inhibition': n_inh,
        'activation_median': act_median,
        'inhibition_median': inh_median,
        'activation_mad_sem': act_mad_sem,
        'inhibition_mad_sem': inh_mad_sem,
        'activation_iqr': (act_q25, act_q75),
        'inhibition_iqr': (inh_q25, inh_q75),
        'ks_stat': ks_stat,
        'ks_p': ks_p,
    },
}
#%% probabilistic latency
event_time = pd.read_csv(event_time_path)
posterior_records = []

for cluname, cell_data in event_time.groupby('cluname', sort=False):
    cell_data = cell_data.sort_values('time_centre_s')
    pre = cell_data['time_centre_s'] < 0
    post = ~pre

    pre_effect = cell_data.loc[pre, 'stim_effect_hz'].to_numpy()
    pre_se = cell_data.loc[pre, 'stim_effect_se_hz'].to_numpy()
    # silent cells can have a numerically exact HAC error of zero
    pre_se = np.maximum(pre_se, MIN_EFFECT_SE_HZ)
    pre_weight = 1 / pre_se**2
    baseline = np.sum(pre_effect * pre_weight) / np.sum(pre_weight)
    baseline_var = 1 / np.sum(pre_weight)

    effect = cell_data.loc[post, 'stim_effect_hz'].to_numpy() - baseline
    effect_se = cell_data.loc[post, 'stim_effect_se_hz'].to_numpy()
    effect_se = np.maximum(effect_se, MIN_EFFECT_SE_HZ)
    weight = 1 / effect_se**2

    # every post-stimulation bin shares the same estimated pre-stimulation offset
    precision = (
        np.diag(weight)
        - np.outer(weight, weight) / (1 / baseline_var + np.sum(weight))
    )
    time_start = cell_data.loc[post, 'time_start_s'].to_numpy()

    onsets = []
    activation_log_bf = []
    inhibition_log_bf = []
    activation_amplitude = []
    inhibition_amplitude = []

    for start in range(len(effect) - MIN_RESPONSE_BINS + 1):
        for stop in range(start + MIN_RESPONSE_BINS, len(effect) + 1):
            interval_precision = np.sum(
                precision[start:stop, start:stop]
            )
            score = np.sum(precision[start:stop], axis=0) @ effect
            posterior_precision = (
                1 / AMPLITUDE_PRIOR_SD_HZ**2 + interval_precision
            )
            posterior_mean = score / posterior_precision
            posterior_sd = 1 / np.sqrt(posterior_precision)
            z = posterior_mean / posterior_sd

            log_bf = (
                -.5 * np.log(
                    AMPLITUDE_PRIOR_SD_HZ**2 * posterior_precision
                )
                + .5 * score**2 / posterior_precision
            )
            activation_log_bf.append(
                np.log(2) + log_bf + log_ndtr(z)
            )
            inhibition_log_bf.append(
                np.log(2) + log_bf + log_ndtr(-z)
            )

            log_density = -.5*z**2 - .5*np.log(2*np.pi)
            activation_amplitude.append(
                posterior_mean
                + posterior_sd * np.exp(log_density - log_ndtr(z))
            )
            inhibition_amplitude.append(
                posterior_mean
                - posterior_sd * np.exp(log_density - log_ndtr(-z))
            )
            onsets.append(time_start[start])

    onsets = np.asarray(onsets)
    activation_log_bf = np.asarray(activation_log_bf)
    inhibition_log_bf = np.asarray(inhibition_log_bf)
    activation_amplitude = np.asarray(activation_amplitude)
    inhibition_amplitude = np.asarray(inhibition_amplitude)

    interval_log_prior = np.log(
        (1 - NO_RESPONSE_PRIOR) / (2 * len(onsets))
    )
    log_weights = np.concatenate([
        [np.log(NO_RESPONSE_PRIOR)],
        interval_log_prior + activation_log_bf,
        interval_log_prior + inhibition_log_bf,
    ])
    weights = np.exp(log_weights - logsumexp(log_weights))
    p_no_response = weights[0]
    activation_weights = weights[1:1 + len(onsets)]
    inhibition_weights = weights[1 + len(onsets):]
    p_activation = np.sum(activation_weights)
    p_inhibition = np.sum(inhibition_weights)

    response_weights = activation_weights + inhibition_weights
    response_weights = response_weights / np.sum(response_weights)
    order = np.argsort(onsets)
    cumulative = np.cumsum(response_weights[order])
    onset_interval = onsets[order][
        np.searchsorted(cumulative, [.025, .5, .975])
    ]

    amplitude_mean = np.sum(
        activation_weights * activation_amplitude
        + inhibition_weights * inhibition_amplitude
    )
    probabilities = [p_no_response, p_activation, p_inhibition]
    posterior_effect = [
        'no response', 'activation', 'inhibition'
    ][int(np.argmax(probabilities))]

    four_of_five_onset, four_of_five_sign = sig_dict[cluname]
    four_of_five_effect = {
        -1: 'inhibition',
        0: 'no response',
        1: 'activation',
    }[four_of_five_sign]

    posterior_records.append({
        'animal': cell_data['animal'].iloc[0],
        'session': cell_data['session'].iloc[0],
        'cluname': cluname,
        'baseline_class': cell_data['baseline_class'].iloc[0],
        'n_cycles': cell_data['n_cycles'].iloc[0],
        'baseline_effect_hz': baseline,
        'baseline_effect_se_hz': np.sqrt(baseline_var),
        'p_no_response': p_no_response,
        'p_activation': p_activation,
        'p_inhibition': p_inhibition,
        'posterior_effect': posterior_effect,
        'onset_ci_low_s': onset_interval[0],
        'onset_median_s': onset_interval[1],
        'onset_ci_high_s': onset_interval[2],
        'amplitude_mean_hz': amplitude_mean,
        'four_of_five_effect': four_of_five_effect,
        'four_of_five_onset_s': four_of_five_onset,
    })

posterior_df = pd.DataFrame(posterior_records)
intermediate_payload['bayesian_latency'] = posterior_df

print('\nprobabilistic latency summary:')
print(posterior_df['posterior_effect'].value_counts().to_string())

with open(data_output_path, 'wb') as f:
    pickle.dump(intermediate_payload, f)
posterior_df.to_csv(posterior_path, index=False)
print_files_saved([
    ('intermediate data', data_output_path),
    ('probabilistic latency', posterior_path),
])
