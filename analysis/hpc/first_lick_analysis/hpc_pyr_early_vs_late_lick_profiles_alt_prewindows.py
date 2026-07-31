'''
Created on Mon Jun 22 2026

compare early- and late-lick HPC pyramidal-cell profiles after reclassifying
PyrUp/PyrDown cells with immediate pre-run windows.

@author: Dinghao Luo
'''

#%% imports
import argparse
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
import matplotlib.pyplot as plt
from scipy.stats import sem, ttest_ind

from common_functions import mpl_formatting
from console_formatting import print_files_saved, print_session, print_status, print_statistics_section
from plotting_functions import plot_violin_with_scatter
import first_lick_analysis_functions as flaf
import project_paths as pp
mpl_formatting()

import rec_list
bad_behs = rec_list.pathHPCbadbeh
paths    = rec_list.pathHPCLCopt + rec_list.pathHPCLCtermopt + rec_list.pathHPC_Raphi
paths    = [p for p in paths if p not in bad_behs]
recnames = [Path(path).name for path in paths]


#%% parameters
SAMP_FREQ     = 1250
RUN_ONSET_BIN = 3750
BEF = 1
AFT = 4

X_SEC_PLOT = np.arange(4000) / 1000.0
PYRUP_CLASS = 'run-onset ON'
PYRDOWN_CLASS = 'run-onset OFF'
EARLY_PYRUP = 'early PyrUp'
LATE_PYRUP = 'late PyrUp'
EARLY_PYRDOWN = 'early PyrDown'
LATE_PYRDOWN = 'late PyrDown'

WINDOW_SPECS = {
    'pre_1s': {
        'label': '[-1, 0] s',
        'folder': 'pre_1s',
        'class_col': 'class_1s',
        'ratio_col': 'pre_post_1s',
        'pre_window': (-1.0, 0.0),
    },
    'pre_half_s': {
        'label': '[-0.5, 0] s',
        'folder': 'pre_half_s',
        'class_col': 'class_half_s',
        'ratio_col': 'pre_post_half_s',
        'pre_window': (-0.5, 0.0),
    },
}
ORIGINAL_WINDOW_SPEC = {
    'key': 'original',
    'label': '[-1.5, -0.5] s',
    'folder': 'original_pre_1p5_to_0p5',
    'class_col': 'class',
    'ratio_col': 'pre_post',
    'pre_window': (-1.5, -0.5),
}

early_c = (168 / 255, 155 / 255, 202 / 255)
late_c  = (102 / 255, 83 / 255, 162 / 255)


#%% path stems
ORIGINAL_HPC_STEM = Path('Z:/Dinghao/code_dinghao/HPC_ephys')
ORIGINAL_BEHAVIOUR_ROOT = Path(
    'Z:/Dinghao/code_dinghao/behaviour/all_experiments'
)

parser = argparse.ArgumentParser(
    description='compare early- and late-lick HPC profiles after alternative PyrUp/PyrDown classifications.',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter,
)
parser.add_argument('--hpc-stem', type=Path, default=ORIGINAL_HPC_STEM,
                    help='folder containing HPC profile pickles and per-session train arrays')
parser.add_argument('--behaviour-root', type=Path, default=ORIGINAL_BEHAVIOUR_ROOT,
                    help='folder containing HPCLC, HPCLCterm, and HPCRaphi behaviour pickles')
parser.add_argument('--repo-inputs', action='store_true',
                    help='use generated in-repo profile, train, and behaviour files')
parser.add_argument('--windows', nargs='+', choices=sorted(WINDOW_SPECS),
                    default=['pre_1s', 'pre_half_s'],
                    help='alternative pre-run windows to plot')
parser.add_argument('--include-original', action='store_true',
                    help='also write the original [-1.5, -0.5] s classification as a reference folder')
args = parser.parse_args()
first_lick_stem = pp.HPC_EPHYS_FIGURES_STEM / 'first_lick_analysis_alt_prewindows'

if args.repo_inputs:
    HPC_stem            = pp.HPC_EPHYS_STEM
    all_exp_stem       = pp.behaviour_experiment_stem('HPCLC')
    all_exp_term_stem  = pp.behaviour_experiment_stem('HPCLCterm')
    all_exp_raphi_stem = pp.behaviour_experiment_stem('HPCRaphi')
else:
    HPC_stem            = args.hpc_stem
    all_exp_stem       = args.behaviour_root / 'HPCLC'
    all_exp_term_stem  = args.behaviour_root / 'HPCLCterm'
    all_exp_raphi_stem = args.behaviour_root / 'HPCRaphi'

WINDOWS = []
if args.include_original:
    WINDOWS.append(ORIGINAL_WINDOW_SPEC.copy())
for key in args.windows:
    spec = WINDOW_SPECS[key].copy()
    spec['key'] = key
    WINDOWS.append(spec)

def _classified_cell_sets(cell_profiles, cell_profiles_raphi, spec):
    class_col = spec['class_col']
    cols = ['recname', class_col]

    df_pyr = cell_profiles[cell_profiles['cell_identity'] == 'pyr']
    df_pyr_raphi = cell_profiles_raphi[cell_profiles_raphi['cell_identity'] == 'pyr']
    df_all = pd.concat([df_pyr[cols], df_pyr_raphi[cols]], axis=0)
    df_all = df_all[df_all['recname'].isin(set(recnames))]

    pyrup = set(df_all[df_all[class_col] == PYRUP_CLASS].index)
    pyrdown = set(df_all[df_all[class_col] == PYRDOWN_CLASS].index)
    return {
        PYRUP_CLASS: pyrup,
        PYRDOWN_CLASS: pyrdown,
    }

def _plot_speed_pre_matched(window_stem,
                            sess_early_speed_means_raw,
                            sess_late_speed_means_raw):
    if not sess_early_speed_means_raw or not sess_late_speed_means_raw:
        print('pre-matching speed plot skipped: no session speed traces')
        return

    E_raw = np.vstack(sess_early_speed_means_raw)
    L_raw = np.vstack(sess_late_speed_means_raw)

    E_raw_mean = np.mean(E_raw, axis=0)
    E_raw_sem  = sem(E_raw, axis=0)
    L_raw_mean = np.mean(L_raw, axis=0)
    L_raw_sem  = sem(L_raw, axis=0)

    fig, ax = plt.subplots(figsize=(2.1, 2.0))
    ax.plot(X_SEC_PLOT, E_raw_mean, c=early_c, label='early (<2.5 s)')
    ax.fill_between(X_SEC_PLOT, E_raw_mean + E_raw_sem, E_raw_mean - E_raw_sem,
                    color=early_c, edgecolor='none', alpha=.25)

    ax.plot(X_SEC_PLOT, L_raw_mean, c=late_c, label='late (2.5-3.5 s)')
    ax.fill_between(X_SEC_PLOT, L_raw_mean + L_raw_sem, L_raw_mean - L_raw_sem,
                    color=late_c, edgecolor='none', alpha=.25)

    ax.set(xlabel='Time from run onset (s)', xlim=(0, 4),
           ylabel='Speed (cm/s)', ylim=(0, 70),
           title='pre-matching speed')
    ax.legend(frameon=False, fontsize=7)
    for s in ['top', 'right']:
        ax.spines[s].set_visible(False)
    fig.tight_layout()

    for ext in ['.png', '.pdf']:
        fig.savefig(window_stem / f'speed_pre_matched{ext}',
                    dpi=300,
                    bbox_inches='tight')
    plt.close(fig)

def _plot_behaviour_summaries(window_stem,
                              sess_early_speed_scalars_raw,
                              sess_late_speed_scalars_raw,
                              sess_early_accel_means_raw,
                              sess_late_accel_means_raw):
    print_statistics_section()

    raw_E = np.asarray(sess_early_speed_scalars_raw, dtype=float)
    raw_L = np.asarray(sess_late_speed_scalars_raw, dtype=float)
    raw_mask = np.isfinite(raw_E) & np.isfinite(raw_L)
    raw_E = raw_E[raw_mask]
    raw_L = raw_L[raw_mask]

    if raw_E.size and raw_L.size:
        plot_violin_with_scatter(raw_E, raw_L,
                                 colour0=early_c, colour1=late_c,
                                 showscatter=False,
                                 xticklabels=['early', 'late'],
                                 stats_labels=[
                                     'pre-matching speed, first lick <2.5 s',
                                     'pre-matching speed, first lick 2.5-3.5 s',
                                 ],
                                 print_statistics=True,
                                 save=True,
                                 savepath=window_stem / 'all_run_onset_speed_violin')

    raw_ea = np.asarray(sess_early_accel_means_raw, dtype=float)
    raw_la = np.asarray(sess_late_accel_means_raw, dtype=float)
    raw_acc_mask = np.isfinite(raw_ea) & np.isfinite(raw_la)
    raw_ea = raw_ea[raw_acc_mask]
    raw_la = raw_la[raw_acc_mask]

    if raw_ea.size and raw_la.size:
        plot_violin_with_scatter(raw_ea, raw_la,
                                 colour0=early_c, colour1=late_c,
                                 showscatter=False,
                                 xticklabels=['early', 'late'],
                                 stats_labels=[
                                     'pre-matching initial acceleration, first lick <2.5 s',
                                     'pre-matching initial acceleration, first lick 2.5-3.5 s',
                                 ],
                                 print_statistics=True,
                                 save=True,
                                 savepath=window_stem / 'all_run_onset_acc_boxplot')

def _plot_speed_post_matched(window_stem,
                             sess_early_speed_means,
                             sess_late_speed_means):
    if not sess_early_speed_means or not sess_late_speed_means:
        print('post-matching speed plot skipped: no session speed traces')
        return

    E = np.vstack(sess_early_speed_means)
    L = np.vstack(sess_late_speed_means)

    E_trace_mean = np.mean(E, axis=0)
    E_trace_sem  = sem(E, axis=0)
    L_trace_mean = np.mean(L, axis=0)
    L_trace_sem  = sem(L, axis=0)

    fig, ax = plt.subplots(figsize=(2.1, 2.0))
    ax.plot(X_SEC_PLOT, E_trace_mean, c=early_c, label='early (<2.5 s)')
    ax.fill_between(X_SEC_PLOT, E_trace_mean + E_trace_sem, E_trace_mean - E_trace_sem,
                    color=early_c, edgecolor='none', alpha=.25)

    ax.plot(X_SEC_PLOT, L_trace_mean, c=late_c, label='late (2.5-3.5 s)')
    ax.fill_between(X_SEC_PLOT, L_trace_mean + L_trace_sem, L_trace_mean - L_trace_sem,
                    color=late_c, edgecolor='none', alpha=.25)

    bin_size = 500
    n_bins = 7
    E_bins = np.vstack([E[:, i * bin_size:(i + 1) * bin_size].mean(axis=1)
                        for i in range(n_bins)]).T
    L_bins = np.vstack([L[:, i * bin_size:(i + 1) * bin_size].mean(axis=1)
                        for i in range(n_bins)]).T

    pvals, tvals = [], []
    dz = []
    for i in range(n_bins):
        t, p = ttest_ind(E_bins[:, i], L_bins[:, i], equal_var=False, nan_policy='omit')
        tvals.append(t)
        pvals.append(p)
        d = E_bins[:, i] - L_bins[:, i]
        d = d[~np.isnan(d)]
        dz.append(np.mean(d) / np.std(d, ddof=1))
    pvals = np.array(pvals)

    ymax = max((E_trace_mean + E_trace_sem).max(), (L_trace_mean + L_trace_sem).max())
    ymin = min((E_trace_mean - E_trace_sem).min(), (L_trace_mean - L_trace_sem).min())
    yr = ymax - ymin if ymax > ymin else 1.0
    bar_y  = ymax + 0.06 * yr
    text_y = ymax + 0.08 * yr

    for i in range(n_bins):
        x_left  = i * 0.5 + 0.1
        x_right = (i + 1) * 0.5 - 0.1
        ax.hlines(bar_y, x_left, x_right, color='k', lw=1)
        ax.text((x_left + x_right) / 2, text_y, f'{pvals[i]:.3e}',
                ha='center', va='bottom', fontsize=3)

    ax.set(xlabel='Time from run onset (s)', xlim=(0, 4),
           ylabel='Speed (cm/s)', ylim=(0, 70),
           title='post-matching speed')
    ax.legend(frameon=False, fontsize=7)
    for s in ['top', 'right']:
        ax.spines[s].set_visible(False)
    fig.tight_layout()

    edges_ms = [(i * 500, (i + 1) * 500) for i in range(n_bins)]
    print('\npost-matching speed binwise t-tests:')
    print('bin | t | p | dz')
    for i, (lo, hi) in enumerate(edges_ms):
        print(f'{lo:>4d} to {hi:<4d} ms | {tvals[i]:.3f} | {pvals[i]:.4f} | {dz[i]:.3f}')

    for ext in ['.png', '.pdf']:
        fig.savefig(window_stem / f'speed_post_matched{ext}',
                    dpi=300,
                    bbox_inches='tight')
    plt.close(fig)

    mean_E = np.mean(E, axis=1)
    mean_L = np.mean(L, axis=1)
    plot_violin_with_scatter(mean_E, mean_L,
                             early_c, late_c,
                             ylabel='Speed (cm/s)',
                             xticklabels=['early', 'late'],
                             stats_labels=[
                                 'post-matching speed, first lick <2.5 s',
                                 'post-matching speed, first lick 2.5-3.5 s',
                             ],
                             print_statistics=True,
                             save=True,
                             savepath=window_stem / 'speed_post_matched_violin')

def _plot_acceleration_post_matched(window_stem,
                                    sess_early_accel_means,
                                    sess_late_accel_means):

    ea = np.asarray(sess_early_accel_means, dtype=float)
    la = np.asarray(sess_late_accel_means, dtype=float)

    mask = np.isfinite(ea) & np.isfinite(la)
    ea = ea[mask]
    la = la[mask]

    if ea.size and la.size:
        plot_violin_with_scatter(ea, la,
                                 early_c, late_c,
                                 ylabel='Acceleration (cm/s^2)',
                                 xticklabels=['early', 'late'],
                                 stats_labels=[
                                     'post-matching initial acceleration, first lick <2.5 s',
                                     'post-matching initial acceleration, first lick 2.5-3.5 s',
                                 ],
                                 print_statistics=True,
                                 save=True,
                                 savepath=window_stem / 'accel_post_matched_violin')

def _plot_mean_profiles(window_stem, spec, profiles):
    window_label = spec['label']
    valid_profiles = {
        key: value
        for key, value in profiles.items()
        if any(len(profile_list) for profile_list in value.values())
        }
    if not valid_profiles:
        print(f'profile plots skipped for {window_label}: no cells survived matching')
        return

    reclist = list(valid_profiles.keys())
    anmlist = [s.split('-')[0] for s in reclist]
    anmlist = np.unique(anmlist)
    print(f'\nprofile summary for {window_label}:')
    print(f'n_animals = {len(anmlist)}')
    print(f'n_sessions = {len(reclist)}')

    early_pyrup   = [arr for session in valid_profiles.values() for arr in session[EARLY_PYRUP]]
    late_pyrup    = [arr for session in valid_profiles.values() for arr in session[LATE_PYRUP]]
    early_pyrdown = [arr for session in valid_profiles.values() for arr in session[EARLY_PYRDOWN]]
    late_pyrdown  = [arr for session in valid_profiles.values() for arr in session[LATE_PYRDOWN]]

    XAXIS = np.arange(5 * SAMP_FREQ) / SAMP_FREQ - 1

    _plot_one_cell_class(window_stem, spec, 'PyrUp', early_pyrup, late_pyrup,
                         'lightcoral', 'firebrick', XAXIS)
    _plot_one_cell_class(window_stem, spec, 'PyrDown', early_pyrdown, late_pyrdown,
                         'violet', 'purple', XAXIS)

def _plot_one_cell_class(window_stem, spec, label, early_profiles, late_profiles,
                         colour_early, colour_late, xaxis):
    window_label = spec['label']
    if not early_profiles or not late_profiles:
        print(f'{label} plot skipped for {window_label}: no cells')
        return

    early_arr = np.asarray(early_profiles)
    late_arr = np.asarray(late_profiles)
    early_mean = np.mean(early_arr, axis=0)
    early_sem = sem(early_arr, axis=0)
    late_mean = np.mean(late_arr, axis=0)
    late_sem = sem(late_arr, axis=0)

    fig, ax = plt.subplots(figsize=(2.3, 2.0))
    ax.plot(xaxis, early_mean, c=colour_early, label='<2.5')
    ax.fill_between(xaxis, early_mean + early_sem, early_mean - early_sem,
                    color=colour_early, edgecolor='none', alpha=.25)
    ax.plot(xaxis, late_mean, c=colour_late, label='2.5~3.5')
    ax.fill_between(xaxis, late_mean + late_sem, late_mean - late_sem,
                    color=colour_late, edgecolor='none', alpha=.25)

    test_label = f'{label} cells, {window_label}'
    p_ind, p_rs, p_rel, p_wil = flaf.compute_binwise_test_suite(
        early_arr,
        late_arr,
        samp_freq=SAMP_FREQ,
        bef=BEF,
        start=-0.5,
        end=3.5,
        bin_size=1,
        label=test_label,
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
        label_fontsize=2.5,
    )

    ax.legend(fontsize=5, frameon=False)
    ax.set(xlabel='Time from run onset (s)', xlim=(-1, 4),
           ylabel='Firing rate (Hz)', ylim=(0.75, 3.15), yticks=[1, 2, 3])
    for s in ['top', 'right']:
        ax.spines[s].set_visible(False)

    fig.tight_layout()

    for ext in ['.png', '.pdf']:
        fig.savefig(
            window_stem / f'all_run_onset_{label}_mean_profiles{ext}',
            dpi=300,
            bbox_inches='tight')
    plt.close(fig)


#%% load data
print('loading dataframes...')

cell_profiles_path = HPC_stem / 'hpc_all_profiles.pkl'
cell_profiles = pd.read_pickle(cell_profiles_path)

cell_profiles_raphi_path = HPC_stem / 'hpc_all_profiles_raphi.pkl'
cell_profiles_raphi = pd.read_pickle(cell_profiles_raphi_path)

if HPC_stem == pp.HPC_EPHYS_STEM:
    expected_raphi_recnames = {
        Path(path).name for path in rec_list.pathHPC_Raphi
        if path not in bad_behs
    }
    found_raphi_recnames = set(cell_profiles_raphi['recname'].dropna().astype(str))
    missing_raphi_recnames = sorted(expected_raphi_recnames - found_raphi_recnames)
    if missing_raphi_recnames:
        raise RuntimeError(
            'hpc_all_profiles_raphi.pkl is incomplete: '
            f'missing {len(missing_raphi_recnames)} Raphi recording(s), '
            f'starting with {missing_raphi_recnames[0]}. '
            'Finish analysis\\hpc\\hpc_all_profiles_raphi.py before running this script.'
        )

class_members = {
    spec['key']: _classified_cell_sets(cell_profiles, cell_profiles_raphi, spec)
    for spec in WINDOWS
}

all_valid_clunames = set().union(*[
    members[PYRUP_CLASS] | members[PYRDOWN_CLASS]
    for members in class_members.values()
])
all_valid_clunames = [
    cluname for cluname in all_valid_clunames
    if cluname.split(' ')[0] in recnames
]
all_valid_clunames = sorted(all_valid_clunames, key=lambda x: x.split(' ')[0])

if not all_valid_clunames:
    raise RuntimeError(
        'No PyrUp or PyrDown cells found under the selected alternative windows. '
        'Check that the profile pickles were regenerated after adding the new columns.'
    )

#%% main
profiles = {spec['key']: {} for spec in WINDOWS}

sess_early_speed_means = []
sess_late_speed_means = []
sess_early_accel_means = []
sess_late_accel_means  = []
sess_early_speed_means_raw = []
sess_late_speed_means_raw = []
sess_early_speed_scalars_raw = []
sess_late_speed_scalars_raw = []
sess_early_accel_means_raw = []
sess_late_accel_means_raw = []

recname   = ''
skip_flag = False
current_session_profiles = {
    spec['key']: {
        EARLY_PYRUP: [],
        LATE_PYRUP: [],
        EARLY_PYRDOWN: [],
        LATE_PYRDOWN: [],
    }
    for spec in WINDOWS
}

for cluname in all_valid_clunames:
    temp_recname = cluname.split(' ')[0]

    if temp_recname != recname:
        if recname:
            for spec in WINDOWS:
                profiles[spec['key']][recname] = current_session_profiles[spec['key']]
        if recname and any(
                any(len(profile_list) for profile_list in session.values())
                for session in current_session_profiles.values()
                ):
            print_status('done')
        current_session_profiles = {
            spec['key']: {
                EARLY_PYRUP: [],
                LATE_PYRUP: [],
                EARLY_PYRDOWN: [],
                LATE_PYRDOWN: [],
            }
            for spec in WINDOWS
        }

        recname = temp_recname
        print_session(recname)
        skip_flag = False

        alignRun_path = None
        behPar_path = None
        for maze_section in range(6):
            raphi_stem = (
                pp.RAPHAEL_ROOT / f'ANM{recname[1:4]}' / recname[:-3]
                / recname
            )
            curr_alignRun_path = (
                raphi_stem / f'{recname}_DataStructure_mazeSection1_'
                f'TrialType1_alignRun_msess{maze_section}.mat'
            )
            curr_behPar_path = (
                raphi_stem / f'{recname}_DataStructure_mazeSection1_'
                f'TrialType1_behPar_msess{maze_section}.mat'
            )
            if curr_alignRun_path.exists() and curr_behPar_path.exists():
                alignRun_path = curr_alignRun_path
                behPar_path = curr_behPar_path
                break

        if alignRun_path is None:
            mice_stem = (
                pp.MICEEXP_ROOT / f'ANMD{recname[1:5]}' / recname[:-3]
                / recname
            )
            alignRun_path = pp.resolve_matlab_pipeline_file(
                mice_stem / f'{recname}_DataStructure_mazeSection1_'
                'TrialType1_alignRun_msess1.mat',
                recname,
            )
            behPar_path = pp.resolve_matlab_pipeline_file(
                mice_stem / f'{recname}_DataStructure_mazeSection1_'
                'TrialType1_behPar_msess1.mat',
                recname,
            )
        alignRun = sio.loadmat(alignRun_path)

        licks = alignRun['trialsRun']['lickLfpInd'][0][0][0][1:]
        starts = alignRun['trialsRun']['startLfpInd'][0][0][0][1:]
        tot_trial = licks.shape[0]

        behPar = sio.loadmat(behPar_path)
        bad_idx  = np.where(behPar['behPar'][0]['indTrBadBeh'][0] == 1)[1] - 1
        stim_idx = np.where(behPar['behPar'][0]['stimOn'][0] != 0)[1] - 1

        beh_path = all_exp_stem / f'{recname}.pkl'
        if not beh_path.exists():
            beh_path = all_exp_term_stem / f'{recname}.pkl'
        if not beh_path.exists():
            beh_path = all_exp_raphi_stem / f'{recname}.pkl'
        with open(beh_path, 'rb') as f:
            beh = pickle.load(f)
        if not skip_flag:
            first_licks = []
            for trial in range(tot_trial):
                lk = [l for l in licks[trial]
                      if l - starts[trial] > .5 * SAMP_FREQ]

                if len(lk) == 0:
                    first_licks.append(np.nan)
                else:
                    first_licks.extend(lk[0] - starts[trial])

            first_licks_sec = np.array(first_licks) / SAMP_FREQ

            early_trials = []
            late_trials = []
            for trial, t in enumerate(first_licks_sec):
                if trial in bad_idx or trial in stim_idx or trial - 1 in stim_idx or np.isnan(t):
                    continue
                if t < 2.5:
                    early_trials.append(trial)
                elif 2.5 < t < 3.5:
                    late_trials.append(trial)

            if len(early_trials) < 10 or len(late_trials) < 10:
                print_status('skipped', 'not enough trials')
                skip_flag = True
                continue
            train_path = HPC_stem / 'all_sessions' / recname / f'{recname}_all_trains_run.npy'
            if not train_path.exists():
                train_path = HPC_stem / 'all_sessions_raphi' / recname / f'{recname}_all_trains.npy'
            all_trains = np.load(train_path, allow_pickle=True).item()
            speed_times = beh['speed_times_aligned'][1:]

            early_speed_scalars = []
            late_speed_scalars = []
            early_accel_trials_raw = []
            late_accel_trials_raw = []

            for t in early_trials:
                if t >= len(speed_times):
                    continue
                sp = np.asarray([pt[1] for pt in speed_times[t]], dtype=float)
                if len(sp) < 1000 or not np.any(np.isfinite(sp)):
                    continue
                early_speed_scalars.append(np.nanmean(sp))
                early_accel_trials_raw.append(sp[999] - sp[0])

            for t in late_trials:
                if t >= len(speed_times):
                    continue
                sp = np.asarray([pt[1] for pt in speed_times[t]], dtype=float)
                if len(sp) < 1000 or not np.any(np.isfinite(sp)):
                    continue
                late_speed_scalars.append(np.nanmean(sp))
                late_accel_trials_raw.append(sp[999] - sp[0])

            if len(early_speed_scalars) > 10 and len(late_speed_scalars) > 10:
                e_speed = np.asarray(early_speed_scalars, dtype=float)
                l_speed = np.asarray(late_speed_scalars, dtype=float)
                e_accel = np.asarray(early_accel_trials_raw, dtype=float)
                l_accel = np.asarray(late_accel_trials_raw, dtype=float)
                if (
                    np.any(np.isfinite(e_speed))
                    and np.any(np.isfinite(l_speed))
                    and np.any(np.isfinite(e_accel))
                    and np.any(np.isfinite(l_accel))
                ):
                    sess_early_speed_scalars_raw.append(np.nanmean(e_speed))
                    sess_late_speed_scalars_raw.append(np.nanmean(l_speed))
                    sess_early_accel_means_raw.append(np.nanmean(e_accel))
                    sess_late_accel_means_raw.append(np.nanmean(l_accel))

            e_mean_sp_raw = flaf.compute_session_mean_speed(early_trials, speed_times, n=4000)
            l_mean_sp_raw = flaf.compute_session_mean_speed(late_trials, speed_times, n=4000)
            if e_mean_sp_raw is not None and l_mean_sp_raw is not None:
                sess_early_speed_means_raw.append(e_mean_sp_raw)
                sess_late_speed_means_raw.append(l_mean_sp_raw)

            E_bins, e_valid = flaf.compute_binned_speed_matrix(
                early_trials, speed_times, n_bins=7, bin_size=500
            )
            L_bins, l_valid = flaf.compute_binned_speed_matrix(
                late_trials, speed_times, n_bins=7, bin_size=500
            )

            matched_early, matched_late = [], []
            if len(E_bins) and len(L_bins):
                k = 1.5
                e_mu = E_bins.mean(axis=0)
                e_sd = E_bins.std(axis=0, ddof=0)
                e_low, e_high = e_mu - k * e_sd, e_mu + k * e_sd

                l_mu = L_bins.mean(axis=0)
                l_sd = L_bins.std(axis=0, ddof=0)
                l_low, l_high = l_mu - k * l_sd, l_mu + k * l_sd

                l_mask_in_early_bounds = np.all((L_bins >= e_low) & (L_bins <= e_high), axis=1)
                e_mask_in_late_bounds  = np.all((E_bins >= l_low) & (E_bins <= l_high), axis=1)

                matched_late  = [l_valid[i] for i in np.where(l_mask_in_early_bounds)[0]]
                matched_early = [e_valid[i] for i in np.where(e_mask_in_late_bounds)[0]]

            if len(matched_early) < 10 or len(matched_late) < 10:
                print_status('skipped', 'not enough speed-matched trials')
                skip_flag = True
                continue

            print_status('analysing', 'passed speed-matching filter')

            ACC_N = int(1 * SAMP_FREQ)
            early_acc_trials = []
            late_acc_trials  = []

            for t in matched_early:
                sp = [pt[1] for pt in speed_times[t]]
                if len(sp) >= ACC_N + 1:
                    s = np.asarray(sp[:ACC_N + 1], dtype=float)
                    a = np.diff(s) * SAMP_FREQ
                    early_acc_trials.append(np.nanmean(a))

            for t in matched_late:
                sp = [pt[1] for pt in speed_times[t]]
                if len(sp) >= ACC_N + 1:
                    s = np.asarray(sp[:ACC_N + 1], dtype=float)
                    a = np.diff(s) * SAMP_FREQ
                    late_acc_trials.append(np.nanmean(a))

            early_acc_trials = np.asarray(early_acc_trials)
            late_acc_trials  = np.asarray(late_acc_trials)

            if early_acc_trials.size and late_acc_trials.size:
                sess_early_accel_means.append(np.nanmean(early_acc_trials))
                sess_late_accel_means.append(np.nanmean(late_acc_trials))

            e_mean_sp = flaf.compute_session_mean_speed(matched_early, speed_times, n=4000)
            l_mean_sp = flaf.compute_session_mean_speed(matched_late, speed_times, n=4000)

            if e_mean_sp is not None and l_mean_sp is not None:
                sess_early_speed_means.append(e_mean_sp)
                sess_late_speed_means.append(l_mean_sp)

    if not skip_flag:
        trains = all_trains[cluname]
        early_profiles = np.nanmean(
            flaf.extract_run_onset_profiles(
                trains, matched_early, RUN_ONSET_BIN, SAMP_FREQ, BEF, AFT
            ),
            axis=0,
        )
        late_profiles = np.nanmean(
            flaf.extract_run_onset_profiles(
                trains, matched_late, RUN_ONSET_BIN, SAMP_FREQ, BEF, AFT
            ),
            axis=0,
        )
        if int(recname[1:4]) > 40 and 'r' not in recname:
            early_profiles = early_profiles * SAMP_FREQ
            late_profiles = late_profiles * SAMP_FREQ

        for spec in WINDOWS:
            members = class_members[spec['key']]
            if cluname in members[PYRUP_CLASS]:
                current_session_profiles[spec['key']][EARLY_PYRUP].append(early_profiles)
                current_session_profiles[spec['key']][LATE_PYRUP].append(late_profiles)
            if cluname in members[PYRDOWN_CLASS]:
                current_session_profiles[spec['key']][EARLY_PYRDOWN].append(early_profiles)
                current_session_profiles[spec['key']][LATE_PYRDOWN].append(late_profiles)

if recname:
    for spec in WINDOWS:
        profiles[spec['key']][recname] = current_session_profiles[spec['key']]
if recname and any(
        any(len(profile_list) for profile_list in session.values())
        for session in current_session_profiles.values()
        ):
    print_status('done')

first_lick_stem.mkdir(parents=True, exist_ok=True)

#%% figures
saved_stems = []
for spec in WINDOWS:
    window_stem = first_lick_stem / spec['folder']
    window_stem.mkdir(parents=True, exist_ok=True)
    label = spec['label']
    saved_stems.append((f'{label} figures', window_stem))

    print(f'\nwriting figures for classification window {label}')
    _plot_speed_pre_matched(window_stem,
                            sess_early_speed_means_raw,
                            sess_late_speed_means_raw)
    _plot_behaviour_summaries(window_stem,
                              sess_early_speed_scalars_raw,
                              sess_late_speed_scalars_raw,
                              sess_early_accel_means_raw,
                              sess_late_accel_means_raw)
    _plot_speed_post_matched(window_stem,
                             sess_early_speed_means,
                             sess_late_speed_means)
    _plot_acceleration_post_matched(window_stem,
                                    sess_early_accel_means,
                                    sess_late_accel_means)
    _plot_mean_profiles(window_stem, spec, profiles[spec['key']])

print_files_saved(saved_stems)
