'''
Created on 29 July 2026
Updated on 30 July 2026: fixed behavioural exclusion and analysis path

survival analysis of first licks

@author: Dinghao Luo
'''


#%% imports
import argparse
from pathlib import Path
import sys

from lifelines import AalenJohansenFitter, CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test
from lifelines.utils import restricted_mean_survival_time
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.io as sio

repo_root = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(repo_root))
sys.path.insert(0, str(repo_root / 'utils'))

from common_functions import mpl_formatting
import project_paths as pp
import rec_list

mpl_formatting()


#%% constants
SAMP_FREQ = 1250
ENTRY_S   = 1


#%% endpoint
def first_lick_endpoint(start, end, licks, pumps):
    '''
    define the first post-landmark lick and its censoring time
    '''
    lick_s = (np.asarray(licks, dtype=float).reshape(-1) - start) / SAMP_FREQ
    pump_s = (np.asarray(pumps, dtype=float).reshape(-1) - start) / SAMP_FREQ

    lick_s = lick_s[lick_s > 0]
    pump_s = pump_s[pump_s > 0]

    trial_end_s = (end - start) / SAMP_FREQ
    reward_s = pump_s[0] if len(pump_s) else np.inf
    censor_s = min(reward_s, trial_end_s)

    qualifying_licks = lick_s[(lick_s > ENTRY_S) & (lick_s < censor_s)]
    event = int(len(qualifying_licks) > 0)
    duration_s = qualifying_licks[0] if event else censor_s

    first_raw_lick_s = lick_s[0] if len(lick_s) else np.nan
    early_lick = int(np.any(lick_s <= ENTRY_S))
    reward_first = int(reward_s < trial_end_s and not event)

    return duration_s, event, reward_first, reward_s, trial_end_s, early_lick, first_raw_lick_s


#%% command line
parser = argparse.ArgumentParser()
parser.add_argument('--cohort', choices=['soma', 'terminal'], default='soma')
parser.add_argument('--rmst-horizon-s', type=float, default=5)
parser.add_argument('--n-bootstrap', type=int, default=2000)
parser.add_argument('--seed', type=int, default=42)
parser.add_argument(
    '--output-dir', type=Path,
    default=pp.BEHAVIOUR_STEM / 'lc_stim_firstlick_survival'
    )
parser.add_argument(
    '--figure-dir', type=Path,
    default=pp.BEHAVIOUR_FIGURES_STEM / 'lc_stim_firstlick_survival'
    )
args = parser.parse_args()

args.output_dir.mkdir(parents=True, exist_ok=True)
args.figure_dir.mkdir(parents=True, exist_ok=True)


#%% build matched trial table
paths = rec_list.pathHPCLCopt if args.cohort == 'soma' else rec_list.pathHPCLCtermopt
# bad-behaviour sessions are a pre-specified exclusion, not an analysis option
paths = [path for path in paths if path not in rec_list.pathHPCbadbeh]

records = []

for path in paths:
    sessname = Path(path).name
    animal = sessname[:5]
    rec_stem = pp.MICEEXP_ROOT / f'ANMD{sessname[1:5]}' / sessname[:14] / sessname

    info = sio.loadmat(
        rec_stem / f'{sessname}_DataStructure_mazeSection1_TrialType1_Info.mat'
        )
    align_run = sio.loadmat(
        rec_stem / f'{sessname}_DataStructure_mazeSection1_TrialType1_alignRun_msess1.mat'
        )

    pulse_method = np.asarray(info['beh'][0][0]['pulseMethod'][0]).reshape(-1)
    trials_run = align_run['trialsRun']
    starts = trials_run['startLfpInd'][0][0][0]
    ends   = trials_run['endLfpInd'][0][0][0]
    licks  = trials_run['lickLfpInd'][0][0][0]
    pumps  = trials_run['pumpLfpInd'][0][0][0]

    # MATLAB trial 0 is an empty placeholder, so pulseMethod keeps its physical indexing here
    stim_trials = np.where(pulse_method != 0)[0]
    stim_trials = stim_trials[stim_trials >= 1]

    for stim_trial in stim_trials:
        control_trial = stim_trial + 2

        # stim + 2 is the control used by the existing first-lick comparison
        if stim_trial - 1 < 1 or control_trial >= len(starts):
            continue
        if not (
                pulse_method[stim_trial - 1] == 0
                and pulse_method[stim_trial + 1] == 0
                and pulse_method[control_trial] == 0
                ):
            continue

        stim_endpoint = first_lick_endpoint(
            starts[stim_trial], ends[stim_trial],
            licks[stim_trial], pumps[stim_trial]
            )
        control_endpoint = first_lick_endpoint(
            starts[control_trial], ends[control_trial],
            licks[control_trial], pumps[control_trial]
            )

        # an early lick alone does not remove a trial; the matched cycle leaves only if
        # reward or trial end has already occurred by the 1 s landmark in either member
        if stim_endpoint[0] <= ENTRY_S or control_endpoint[0] <= ENTRY_S:
            continue

        for trial, is_stim, endpoint in [
                (stim_trial, 1, stim_endpoint),
                (control_trial, 0, control_endpoint)
                ]:
            records.append({
                'cohort': args.cohort,
                'animal': animal,
                'session': sessname,
                'cycle': stim_trial,
                'trial': trial,
                'condition': 'stim' if is_stim else 'control',
                'is_stim': is_stim,
                'entry_s': ENTRY_S,
                'duration_s': endpoint[0],
                'event': endpoint[1],
                'event_type': 1 if endpoint[1] else 2 if endpoint[2] else 0,
                'reward_first': endpoint[2],
                'reward_s': endpoint[3],
                'trial_end_s': endpoint[4],
                'lick_before_1s': endpoint[5],
                'first_raw_lick_s': endpoint[6]
                })

trials = pd.DataFrame(records)
trials['time_from_landmark_s'] = trials['duration_s'] - ENTRY_S

suffix = args.cohort
trials.to_csv(args.output_dir / f'{suffix}_firstlick_survival_trials.csv', index=False)


#%% Kaplan-Meier curves and log-rank test
stim = trials['is_stim'] == 1
control = ~stim

# Greenwood variance, with lifelines' log-log transform keeping intervals inside 0 and 1
km_stim = KaplanMeierFitter(label='stimulation')
km_control = KaplanMeierFitter(label='control')

km_stim.fit(
    trials.loc[stim, 'duration_s'], trials.loc[stim, 'event'],
    entry=trials.loc[stim, 'entry_s']
    )
km_control.fit(
    trials.loc[control, 'duration_s'], trials.loc[control, 'event'],
    entry=trials.loc[control, 'entry_s']
    )

logrank = logrank_test(
    # lifelines has no entry argument here; with one common landmark, shifting both
    # clocks by 1 s gives the same risk sets as the delayed-entry KM fits above
    trials.loc[stim, 'time_from_landmark_s'],
    trials.loc[control, 'time_from_landmark_s'],
    event_observed_A=trials.loc[stim, 'event'],
    event_observed_B=trials.loc[control, 'event']
    )


#%% Cox models
# all trials enter at the same landmark, so shifting the Cox clock to zero is the same delayed-entry analysis
cox_data = trials[['time_from_landmark_s', 'event', 'is_stim', 'session']].copy()

cph = CoxPHFitter()
cph.fit(
    cox_data[['time_from_landmark_s', 'event', 'is_stim']],
    duration_col='time_from_landmark_s', event_col='event'
    )

cph_clustered = CoxPHFitter()
cph_clustered.fit(
    cox_data,
    duration_col='time_from_landmark_s', event_col='event',
    cluster_col='session', robust=True
    )

print('\nproportional-hazards check')
cph.check_assumptions(
    cox_data[['time_from_landmark_s', 'event', 'is_stim']],
    p_value_threshold=0.05, show_plots=False
    )


#%% RMST and session bootstrap
# tau is only the reporting horizon; neither the KM curves nor trial follow-up are cut at this time
tau = args.rmst_horizon_s
# lifelines integrates from the first entry time, so this is already RMST from 1 s to tau
rmst_stim = restricted_mean_survival_time(km_stim, t=tau)
rmst_control = restricted_mean_survival_time(km_control, t=tau)
rmst_difference = rmst_stim - rmst_control

rng = np.random.default_rng(args.seed)
sessions = trials['session'].unique()
plot_times = np.linspace(ENTRY_S, trials['duration_s'].max(), 300)
bootstrap_rmst_differences = []
bootstrap_curve_differences = []

for _ in range(args.n_bootstrap):
    sampled_sessions = rng.choice(sessions, len(sessions), replace=True)
    sampled = pd.concat(
        [trials.loc[trials['session'] == session] for session in sampled_sessions],
        ignore_index=True
        )

    sampled_stim = sampled['is_stim'] == 1
    sampled_control = ~sampled_stim

    km_stim_boot = KaplanMeierFitter().fit(
        sampled.loc[sampled_stim, 'duration_s'],
        sampled.loc[sampled_stim, 'event'],
        entry=sampled.loc[sampled_stim, 'entry_s']
        )
    km_control_boot = KaplanMeierFitter().fit(
        sampled.loc[sampled_control, 'duration_s'],
        sampled.loc[sampled_control, 'event'],
        entry=sampled.loc[sampled_control, 'entry_s']
        )

    bootstrap_rmst_differences.append(
        restricted_mean_survival_time(km_stim_boot, t=tau)
        - restricted_mean_survival_time(km_control_boot, t=tau)
        )
    bootstrap_curve_differences.append(
        np.asarray(km_stim_boot.predict(plot_times))
        - np.asarray(km_control_boot.predict(plot_times))
        )

rmst_ci_low, rmst_ci_high = np.percentile(
    bootstrap_rmst_differences, [2.5, 97.5]
    )
curve_difference = (
    np.asarray(km_stim.predict(plot_times))
    - np.asarray(km_control.predict(plot_times))
    )
curve_ci_low, curve_ci_high = np.percentile(
    bootstrap_curve_differences, [2.5, 97.5], axis=0
    )


#%% competing-risk CIF
# here reward is an observed competing event rather than a censoring assumption
aj_stim = AalenJohansenFitter(seed=args.seed, label='stimulation')
aj_control = AalenJohansenFitter(seed=args.seed, label='control')
aj_stim.fit(
    trials.loc[stim, 'duration_s'], trials.loc[stim, 'event_type'],
    event_of_interest=1, entry=trials.loc[stim, 'entry_s']
    )
aj_control.fit(
    trials.loc[control, 'duration_s'], trials.loc[control, 'event_type'],
    event_of_interest=1, entry=trials.loc[control, 'entry_s']
    )
cif_stim_tau = float(aj_stim.predict(tau))
cif_control_tau = float(aj_control.predict(tau))
cif_difference_tau = cif_stim_tau - cif_control_tau


#%% save results
cox_row = cph.summary.loc['is_stim']
clustered_row = cph_clustered.summary.loc['is_stim']

results = pd.Series({
    'cohort': args.cohort,
    'n_animals': trials['animal'].nunique(),
    'n_sessions': trials['session'].nunique(),
    'n_cycles': trials['cycle'].count() // 2,
    'stim_first_lick_events': trials.loc[stim, 'event'].sum(),
    'control_first_lick_events': trials.loc[control, 'event'].sum(),
    'stim_reward_first': trials.loc[stim, 'reward_first'].sum(),
    'control_reward_first': trials.loc[control, 'reward_first'].sum(),
    'stim_censored': (~trials.loc[stim, 'event'].astype(bool)).sum(),
    'control_censored': (~trials.loc[control, 'event'].astype(bool)).sum(),
    'logrank_statistic': logrank.test_statistic,
    'logrank_p': logrank.p_value,
    'cox_hazard_ratio': cox_row['exp(coef)'],
    'cox_p': cox_row['p'],
    'clustered_cox_hazard_ratio': clustered_row['exp(coef)'],
    'clustered_cox_p': clustered_row['p'],
    'rmst_horizon_s': tau,
    'stim_rmst_1_to_tau_s': rmst_stim,
    'control_rmst_1_to_tau_s': rmst_control,
    'rmst_difference_s': rmst_difference,
    'rmst_difference_ci_low_s': rmst_ci_low,
    'rmst_difference_ci_high_s': rmst_ci_high,
    'stim_firstlick_cif_at_tau': cif_stim_tau,
    'control_firstlick_cif_at_tau': cif_control_tau,
    'firstlick_cif_difference_at_tau': cif_difference_tau
    })
results.to_csv(args.output_dir / f'{suffix}_firstlick_survival_summary.csv', header=False)

print('\n' + results.to_string())


#%% plot
fig, (ax, ax_difference) = plt.subplots(
    2, 1, figsize=(4.2, 5.2), sharex=True,
    gridspec_kw={'height_ratios': [3, 1.35]}
    )
km_control.plot_survival_function(ax=ax, ci_show=True, color='0.35')
km_stim.plot_survival_function(ax=ax, ci_show=True, color='#3182bd')
ax.axvline(tau, color='0.7', linestyle=':', linewidth=1, label='RMST horizon')
ax.set(
    ylabel='probability without a qualifying first lick',
    xlim=(ENTRY_S, None),
    ylim=(0, 1.02),
    title=f'LC-{args.cohort} stimulation'
    )
ax.legend(frameon=False)

# the lower panel shows the contrast itself; its bands are session-bootstrap bands
ax_difference.fill_between(
    plot_times, curve_ci_low, curve_ci_high,
    color='#3182bd', alpha=0.2, linewidth=0
    )
ax_difference.plot(
    plot_times, curve_difference,
    color='#3182bd', linewidth=1.5
    )
ax_difference.axhline(0, color='0.35', linewidth=0.8)
ax_difference.axvline(tau, color='0.7', linestyle=':', linewidth=1)
ax_difference.set(
    xlabel='time from run onset (s)',
    ylabel='stim - control',
    xlim=(ENTRY_S, None)
    )
fig.tight_layout()
fig.savefig(
    args.figure_dir / f'{suffix}_firstlick_survival.png',
    dpi=300, bbox_inches='tight'
    )

# the endpoint starts at 1 s, so this is the observable part of the first five seconds
ax_difference.set_xlim(ENTRY_S, 5)
fig.savefig(
    args.figure_dir / f'{suffix}_firstlick_survival_first5s.png',
    dpi=300, bbox_inches='tight'
    )
plt.close(fig)


#%% competing-risk CIF plot
fig_cif, ax_cif = plt.subplots(figsize=(4.2, 3.4))
aj_control.plot_cumulative_density(ax=ax_cif, ci_show=True, color='0.35')
aj_stim.plot_cumulative_density(ax=ax_cif, ci_show=True, color='#3182bd')
ax_cif.axvline(tau, color='0.7', linestyle=':', linewidth=1, label='RMST horizon')
ax_cif.set(
    xlabel='time from run onset (s)',
    ylabel='probability of first lick before reward',
    xlim=(ENTRY_S, None),
    ylim=(0, 1.02),
    title=f'LC-{args.cohort} first-lick CIF'
    )
ax_cif.legend(frameon=False)
fig_cif.tight_layout()
fig_cif.savefig(
    args.figure_dir / f'{suffix}_firstlick_cif.png',
    dpi=300, bbox_inches='tight'
    )
ax_cif.set_xlim(ENTRY_S, 5)
fig_cif.savefig(
    args.figure_dir / f'{suffix}_firstlick_cif_first5s.png',
    dpi=300, bbox_inches='tight'
    )
plt.close(fig_cif)
