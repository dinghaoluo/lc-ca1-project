# -*- coding: utf-8 -*-
'''
Created on 13 April 2026

compare first-lick distance or time between stimulation and matched
baseline/control trials in HPCLC and HPCLCterm sessions

@author: Dinghao Luo
'''

#%% imports
import argparse
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import scipy.io as sio
from scipy.stats import ranksums

repo_root = Path(__file__).resolve().parents[3]
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))

from common_functions import mpl_formatting
from console_formatting import print_session, print_status
from plotting_functions import plot_violin_with_scatter
import project_paths as pp
import rec_list
mpl_formatting()


#%% constants
SAMP_FREQ = 1250
MICEEXP_ROOT = pp.MICEEXP_ROOT

EXPERIMENT_CONFIGS = {
    'lc': {
        'dist_paths': rec_list.pathHPCLCopt[9:],
        'time_paths': rec_list.pathHPCLCopt,
        'dist_output': pp.LC_OPTO_EPHYS_FIGURES_STEM,
        'time_output': pp.BEHAVIOUR_FIGURES_STEM / 'HPC_opto',
        'dist_label': 'HPC_LCstim',
        'time_label': 'HPC_LC_stim',
    },
    'lcterm': {
        'dist_paths': rec_list.pathHPCLCtermopt_beh,
        'time_paths': rec_list.pathHPCLCtermopt,
        'dist_output': pp.HPC_OPTO_EPHYS_FIGURES_STEM,
        'time_output': pp.HPC_OPTO_EPHYS_FIGURES_STEM,
        'dist_label': 'HPC_LCterm_stim',
        'time_label': 'HPC_LCterm_stim',
    },
}


#%% plotting
def _plot_pair(values0, values1, title, xlabel, outdir, filename, xlim=None):
    fig, ax = plt.subplots(figsize=(4, 3))
    for spine in ['top', 'right', 'left']:
        ax.spines[spine].set_visible(False)
    ax.set(title=title, ylim=(0, 1.5), xlabel=xlabel)
    if xlim is not None:
        ax.set_xlim(*xlim)
    ax.set_yticks([.5, 1])
    ax.set_yticklabels(['baseline', 'stim'])
    ax.scatter(values0, [.5] * len(values0), color='grey')
    ax.scatter(values1, [1] * len(values1), color='darkblue')
    if len(values0) > 0 and len(values1) > 0:
        ax.plot([np.median(values0), np.median(values1)], [.5, 1],
                color='grey', alpha=.5)
    outdir.mkdir(parents=True, exist_ok=True)
    fig.savefig(outdir / filename, dpi=300, bbox_inches='tight')
    plt.close(fig)


#%% analyses
def run_distance(experiment='lc', n_bst=1000, comp_method='baseline'):
    config = EXPERIMENT_CONFIGS[experiment]
    sess_list = [Path(sess).name[-17:] for sess in config['dist_paths']]
    output_stem = config['dist_output']
    label = config['dist_label']

    print(f'\nbootstrap n = {n_bst}')
    print(f'comparison method = {comp_method}')

    all_licks_non_stim = []
    all_licks_stim = []

    for sessname in sess_list:
        print_session(sessname)
        rec_stem = MICEEXP_ROOT / f'ANMD{sessname[1:5]}' / sessname[:14] / sessname
        info_path = rec_stem / f'{sessname}_DataStructure_mazeSection1_TrialType1_Info.mat'
        align_path = rec_stem / f'{sessname}_DataStructure_mazeSection1_TrialType1_alignRun_msess1.mat'
        info = sio.loadmat(info_path)
        align_run = sio.loadmat(align_path)
        pulse_method = info['beh'][0][0]['pulseMethod'][0]
        trials_run = align_run['trialsRun']

        nonzero = pulse_method[np.where(pulse_method != 0)]
        stim_cond = int(np.asarray(nonzero[0]).squeeze().item())
        stim = [idx for idx, val in enumerate(pulse_method) if val == stim_cond]

        licks = trials_run['lickLfpInd'][0][0][0][1:]
        starts = trials_run['startLfpInd'][0][0][0][1:]
        dists = trials_run['xMM'][0][0][0][1:]
        first_licks = []
        for trial in range(licks.shape[0]):
            start = np.asarray(starts[trial]).squeeze().item()
            lick_indices = []
            for lick in licks[trial]:
                lick = np.asarray(lick).squeeze().item()
                if lick - start > SAMP_FREQ:
                    lick_indices.append(int(round(lick)))

            first_lick = 0.0
            dist_trace = np.asarray(dists[trial], dtype=float).squeeze()
            for lick_idx in lick_indices:
                rel_idx = int(round(lick_idx - start))
                if 0 <= rel_idx < len(dist_trace):
                    lick_dist = float(dist_trace[rel_idx]) / 10
                    if lick_dist > 30:
                        first_lick = lick_dist
                        break
            first_licks.append(first_lick)

        licks_stim = []
        for idx in stim:
            value_idx = idx - 1
            next_idx = idx + 1
            if not 0 <= value_idx < len(first_licks):
                continue
            if not 0 <= next_idx < len(first_licks) or first_licks[next_idx] == 0:
                continue
            if first_licks[value_idx] != 0:
                licks_stim.append(float(first_licks[value_idx]))
        if len(licks_stim) == 0:
            print_status('skipped', 'no valid stim licks')
            continue

        curr_licks_non_stim = []
        curr_licks_stim = []
        licks_non_stim = []
        for _ in range(n_bst):
            if comp_method == 'stim_cont':
                selected_non_stim = [idx + 2 for idx in stim]
            else:
                non_stim_trials = np.where(pulse_method == 0)[0]
                pool_size = min(max(stim[0] - 1, 1), len(non_stim_trials))
                selected_non_stim = non_stim_trials[
                    np.random.randint(0, pool_size, len(licks_stim))
                ]
            if comp_method == 'baseline':
                licks_non_stim = []
                for trial in selected_non_stim:
                    if first_licks[trial - 1] != 0:
                        licks_non_stim.append(float(first_licks[trial - 1]))
                    else:
                        licks_non_stim.append(float(first_licks[trial]))
            else:
                licks_non_stim = []
                for idx in selected_non_stim:
                    value_idx = int(idx) - 1
                    pre_stim_idx = int(idx) - 3
                    if (
                            0 <= value_idx < len(first_licks)
                            and 0 <= pre_stim_idx < len(first_licks)
                            and first_licks[value_idx] != 0
                            and first_licks[pre_stim_idx] != 0
                            ):
                        licks_non_stim.append(float(first_licks[value_idx]))
            if len(licks_non_stim) == 0:
                continue
            curr_licks_non_stim.append(licks_non_stim)
            curr_licks_stim.append(licks_stim)
            ranksums(licks_non_stim, licks_stim)

        if len(curr_licks_non_stim) == 0:
            print_status('skipped', 'no valid baseline licks')
            continue

        if stim_cond == 2:
            all_licks_non_stim.append(np.median(curr_licks_non_stim))
            all_licks_stim.append(np.median(curr_licks_stim))

        suffix = '' if comp_method == 'baseline' else '_stim_cont'
        outdir = output_stem / f'opto_lickdist_0{stim_cond}0{suffix}_{label}'
        _plot_pair(
            licks_non_stim,
            licks_stim,
            f'{sessname}, stim={stim_cond}',
            'dist. 1st lick (cm)',
            outdir,
            f'{sessname}.png',
            xlim=(30, 225),
        )

    if len(all_licks_non_stim) > 0:
        suffix = '' if comp_method == 'baseline' else '_stim_cont'
        summary_dir = output_stem / f'opto_lickdist_020{suffix}_{label}'
        summary_dir.mkdir(parents=True, exist_ok=True)
        plot_violin_with_scatter(
            all_licks_non_stim,
            all_licks_stim,
            'grey',
            'royalblue',
            xticklabels=('ctrl.', 'stim.'),
            ylabel='dist. 1st licks (cm)',
            print_statistics=True,
            save=True,
            savepath=str(summary_dir / 'summary_wilc'),
            show=False,
            close=True,
        )

def run_time(experiment='lc', n_bst=1000, comp_method='baseline'):
    config = EXPERIMENT_CONFIGS[experiment]
    sess_list = [Path(sess).name[-17:] for sess in config['time_paths']]
    output_stem = config['time_output']
    label = config['time_label']

    print(f'\nbootstrap n = {n_bst}')
    print(f'comparison method = {comp_method}')

    all_licks_non_stim = []
    all_licks_stim = []
    all_mspeeds_non_stim = []
    all_mspeeds_stim = []
    all_pspeeds_non_stim = []
    all_pspeeds_stim = []
    all_initacc_non_stim = []
    all_initacc_stim = []

    for sessname in sess_list:
        print_session(sessname)
        rec_stem = MICEEXP_ROOT / f'ANMD{sessname[1:5]}' / sessname[:14] / sessname
        info_path = rec_stem / f'{sessname}_DataStructure_mazeSection1_TrialType1_Info.mat'
        align_path = rec_stem / f'{sessname}_DataStructure_mazeSection1_TrialType1_alignRun_msess1.mat'
        info = sio.loadmat(info_path)
        align_run = sio.loadmat(align_path)
        pulse_method = info['beh'][0][0]['pulseMethod'][0]
        trials_run = align_run['trialsRun']

        nonzero = pulse_method[np.where(pulse_method != 0)]
        stim_cond = int(np.asarray(nonzero[0]).squeeze().item())
        stim = [idx for idx, val in enumerate(pulse_method) if val == stim_cond]

        licks = trials_run['lickLfpInd'][0][0][0][1:]
        starts = trials_run['startLfpInd'][0][0][0][1:]
        pumps = trials_run['pumpLfpInd'][0][0][0][1:]
        speeds = trials_run['speed_MMsec'][0][0][0][1:]

        gx_speed = np.arange(-50, 50, 1)
        sigma_speed = SAMP_FREQ / 100
        gaus_speed = [
            1 / (sigma_speed * np.sqrt(2 * np.pi)) * np.exp(-x**2 / (2 * sigma_speed**2))
            for x in gx_speed
        ]
        smoothed_speeds = [
            np.convolve(np.asarray(speed, dtype=float).squeeze(), gaus_speed, mode='same')
            for speed in speeds
        ]
        init_accels = [
            float(np.mean(np.gradient(speed)[:625]) * 10)
            for speed in smoothed_speeds
        ]

        first_licks = []
        licks_bef_rew = []
        mean_speeds = []
        peak_speeds = []
        for trial in range(licks.shape[0]):
            start = np.asarray(starts[trial]).squeeze().item()
            pump = (
                np.asarray(pumps[trial][0]).squeeze().item() - start
                if len(pumps[trial]) > 0 else 20000
            )
            lick_indices = []
            for lick in licks[trial]:
                lick = np.asarray(lick).squeeze().item()
                if lick - start > SAMP_FREQ:
                    lick_indices.append(lick)

            if len(lick_indices) > 0:
                first_licks.append(float((lick_indices[0] - start) / SAMP_FREQ))
            else:
                first_licks.append(0.0)
            licks_bef_rew.append(
                len([lick for lick in lick_indices if lick - start < pump])
            )

            speed_trace = np.asarray(speeds[trial], dtype=float).squeeze()
            mean_speeds.append(float(np.mean(speed_trace) / 10))
            peak_speeds.append(float(np.max(speed_trace) / 10))

        licks_stim = []
        for idx in stim:
            value_idx = idx - 1
            next_idx = idx + 1
            if not 0 <= value_idx < len(first_licks):
                continue
            if not 0 <= next_idx < len(first_licks) or first_licks[next_idx] == 0:
                continue
            if first_licks[value_idx] != 0:
                licks_stim.append(float(first_licks[value_idx]))
        licks_br_stim = [
            float(licks_bef_rew[int(idx) - 1])
            for idx in stim
            if 0 <= int(idx) - 1 < len(licks_bef_rew)
        ]
        if len(licks_stim) == 0:
            print_status('skipped', 'no valid stim licks')
            continue

        curr_licks_non_stim = []
        curr_licks_stim = []
        curr_mspeeds_non_stim = []
        curr_mspeeds_stim = []
        curr_pspeeds_non_stim = []
        curr_pspeeds_stim = []
        curr_initacc_non_stim = []
        curr_initacc_stim = []
        pval = []
        pval_mspeeds = []
        pval_pspeeds = []
        pval_bef_rew = []
        pval_initacc = []

        licks_non_stim = []
        mspeeds_non_stim = []
        mspeeds_stim = []
        pspeeds_non_stim = []
        pspeeds_stim = []
        initaccs_non_stim = []
        initaccs_stim = []

        for _ in range(n_bst):
            if comp_method == 'stim_cont':
                selected_non_stim = [idx + 2 for idx in stim]
            else:
                non_stim_trials = np.where(pulse_method == 0)[0]
                pool_size = min(max(stim[0] - 1, 1), len(non_stim_trials))
                selected_non_stim = non_stim_trials[
                    np.random.randint(0, pool_size, len(licks_stim))
                ]
            if comp_method == 'baseline':
                licks_non_stim = [
                    float(first_licks[int(idx) - 1])
                    for idx in selected_non_stim
                    if 0 <= int(idx) - 1 < len(first_licks)
                ]
                licks_br_non_stim = [
                    float(licks_bef_rew[int(idx) - 1])
                    for idx in selected_non_stim
                    if 0 <= int(idx) - 1 < len(licks_bef_rew)
                ]
            else:
                licks_non_stim = []
                for idx in selected_non_stim:
                    value_idx = int(idx) - 1
                    pre_stim_idx = int(idx) - 3
                    if (
                            0 <= value_idx < len(first_licks)
                            and 0 <= pre_stim_idx < len(first_licks)
                            and first_licks[value_idx] != 0
                            and first_licks[pre_stim_idx] != 0
                            ):
                        licks_non_stim.append(float(first_licks[value_idx]))
                licks_br_non_stim = [
                    float(licks_bef_rew[int(idx) - 1])
                    for idx in selected_non_stim
                    if 0 <= int(idx) - 1 < len(licks_bef_rew)
                ]

            if len(licks_non_stim) == 0:
                continue

            mspeeds_non_stim = [
                float(mean_speeds[int(idx) - 1])
                for idx in selected_non_stim
                if 0 <= int(idx) - 1 < len(mean_speeds)
            ]
            mspeeds_stim = [
                float(mean_speeds[int(idx) - 1])
                for idx in stim
                if 0 <= int(idx) - 1 < len(mean_speeds)
            ]
            pspeeds_non_stim = [
                float(peak_speeds[int(idx) - 1])
                for idx in selected_non_stim
                if 0 <= int(idx) - 1 < len(peak_speeds)
            ]
            pspeeds_stim = [
                float(peak_speeds[int(idx) - 1])
                for idx in stim
                if 0 <= int(idx) - 1 < len(peak_speeds)
            ]
            initaccs_non_stim = [
                float(init_accels[int(idx) - 1])
                for idx in selected_non_stim
                if 0 <= int(idx) - 1 < len(init_accels)
            ]
            initaccs_stim = [
                float(init_accels[int(idx) - 1])
                for idx in stim
                if 0 <= int(idx) - 1 < len(init_accels)
            ]

            curr_licks_non_stim.append(licks_non_stim)
            curr_licks_stim.append(licks_stim)
            curr_mspeeds_non_stim.append(mspeeds_non_stim)
            curr_mspeeds_stim.append(mspeeds_stim)
            curr_pspeeds_non_stim.append(pspeeds_non_stim)
            curr_pspeeds_stim.append(pspeeds_stim)
            curr_initacc_non_stim.append(initaccs_non_stim)
            curr_initacc_stim.append(initaccs_stim)

            pval.append(ranksums(licks_non_stim, licks_stim)[1])
            pval_mspeeds.append(ranksums(mspeeds_non_stim, mspeeds_stim)[1])
            pval_pspeeds.append(ranksums(pspeeds_non_stim, pspeeds_stim)[1])
            pval_initacc.append(ranksums(initaccs_non_stim, initaccs_stim)[1])
            pval_bef_rew.append(ranksums(licks_br_non_stim, licks_br_stim)[1])

        if len(curr_licks_non_stim) == 0:
            print_status('skipped', 'no valid baseline licks')
            continue

        if stim_cond == 2:
            all_licks_non_stim.append(np.median(curr_licks_non_stim))
            all_licks_stim.append(np.median(curr_licks_stim))
            all_mspeeds_non_stim.append(np.median(curr_mspeeds_non_stim))
            all_mspeeds_stim.append(np.median(curr_mspeeds_stim))
            all_pspeeds_non_stim.append(np.median(curr_pspeeds_non_stim))
            all_pspeeds_stim.append(np.median(curr_pspeeds_stim))
            all_initacc_non_stim.append(np.median(curr_initacc_non_stim))
            all_initacc_stim.append(np.median(curr_initacc_stim))

        suffix = '' if comp_method == 'baseline' else '_stim_cont'
        outdir = output_stem / f'opto_licktime_0{stim_cond}0{suffix}_{label}'
        title = f'{sessname}, stim={stim_cond}, p={np.mean(pval):.4f}, {np.mean(pval_bef_rew):.4f}'
        _plot_pair(licks_non_stim, licks_stim, title, 't 1st lick (s)', outdir, f'{sessname}.png')
        _plot_pair(
            mspeeds_non_stim,
            mspeeds_stim,
            f'{sessname}, stim={stim_cond}, p={np.mean(pval_mspeeds):.4f}',
            'mean velocity (cm/s)',
            outdir,
            f'{sessname}_control_velocity.png',
        )
        _plot_pair(
            pspeeds_non_stim,
            pspeeds_stim,
            f'{sessname}, stim={stim_cond}, p={np.mean(pval_pspeeds):.4f}',
            'peak velocity (cm/s)',
            outdir,
            f'{sessname}_control_peak_velocity.png',
        )
        _plot_pair(
            initaccs_non_stim,
            initaccs_stim,
            f'{sessname}, stim={stim_cond}, p={np.mean(pval_initacc):.4f}',
            'init. accel. (cm/s^2)',
            outdir,
            f'{sessname}_control_init_accel.png',
        )

    if len(all_licks_non_stim) > 0:
        suffix = '' if comp_method == 'baseline' else '_stim_cont'
        summary_dir = output_stem / f'opto_licktime_020{suffix}_{label}'
        summary_dir.mkdir(parents=True, exist_ok=True)
        for values0, values1, ylabel, filename in [
                (all_licks_non_stim, all_licks_stim,
                 '1st-lick time (s)', 'summary_wilc'),
                (all_mspeeds_non_stim, all_mspeeds_stim,
                 'Mean speed (m/s)', 'summary_wilc_control_velocity'),
                (all_pspeeds_non_stim, all_pspeeds_stim,
                 'Peak speed (m/s)', 'summary_wilc_control_peak_velocity'),
                (all_initacc_non_stim, all_initacc_stim,
                 'Init. accel. (m/s^2)', 'summary_wilc_control_init_accel'),
                ]:
            plot_violin_with_scatter(
                values0,
                values1,
                'grey',
                'royalblue',
                xticklabels=('ctrl.', 'stim.'),
                ylabel=ylabel,
                print_statistics=True,
                save=True,
                savepath=str(summary_dir / filename),
                show=False,
                close=True,
            )

def main(argv=None):
    parser = argparse.ArgumentParser(
        description='compare HPC first-lick distance/time for LC or LC-terminal stimulation.'
    )
    parser.add_argument('--metric', choices=['dist', 'time', 'all'], default='all')
    parser.add_argument('--experiment', choices=['lc', 'lcterm', 'all'], default='all')
    parser.add_argument('--n-bst', type=int, default=1000)
    parser.add_argument('--comp-method', choices=['baseline', 'stim_cont'], default='baseline')
    args = parser.parse_args(argv)

    metrics = ['dist', 'time'] if args.metric == 'all' else [args.metric]
    experiments = ['lc', 'lcterm'] if args.experiment == 'all' else [args.experiment]
    for curr_experiment in experiments:
        for curr_metric in metrics:
            print(f'{curr_experiment} {curr_metric}')
            if curr_metric == 'dist':
                run_distance(
                    curr_experiment,
                    n_bst=args.n_bst,
                    comp_method=args.comp_method,
                )
            else:
                run_time(
                    curr_experiment,
                    n_bst=args.n_bst,
                    comp_method=args.comp_method,
                )

if __name__ == '__main__':
    main()
