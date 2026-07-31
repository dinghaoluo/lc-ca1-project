# -*- coding: utf-8 -*-
'''
Created on Wed 13 Nov 16:53:32 2024
Modified on Thur 6 Mar 2025

summarise optogenetic experiments

@author: Dinghao Luo
'''

#%% imports
from pathlib import Path
import sys

repo_root = Path(__file__).resolve().parents[2]
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
utils_path = repo_root / 'utils'
if str(utils_path) not in sys.path:
    sys.path.insert(0, str(utils_path))

import numpy as np

from common_functions import mpl_formatting
from console_formatting import print_session
from behaviour_functions import process_behavioural_data
import plotting_functions as pf
import project_paths as pp
import rec_list

mpl_formatting()


#%% constants
OPTO_OUTPUT_ROOT = pp.BEHAVIOUR_LC_OPTO_FIGURES_STEM
OPTO_OUTPUT_ROOT.mkdir(parents=True, exist_ok=True)


#%% trial and session summaries

def _segment_opto_session(curr_txt, sessname):
    # find the opto bit
    optogenetic_protocol = [trial[15] for trial in curr_txt['trial_statements']]
    start_idx = end_idx = None

    # first stim trial
    for trial_idx, protocol in enumerate(optogenetic_protocol):
        if protocol != '0':
            start_idx = trial_idx
            break

    if start_idx is not None:
        # three zeros mark the end
        for trial_idx in range(start_idx, len(optogenetic_protocol) - 2):
            if (
                optogenetic_protocol[trial_idx] == '0'
                and optogenetic_protocol[trial_idx + 1] == '0'
                and optogenetic_protocol[trial_idx + 2] == '0'
            ):
                end_idx = trial_idx  # end before the three consecutive 0s
                break

    # opto carries on to the end
    if start_idx is not None and end_idx is None:
        end_idx = len(optogenetic_protocol)

    # no opto in this session
    if start_idx is None:
        print(f'{sessname}: no optogenetic stim. in this session')
        baseline_trials = list(range(len(optogenetic_protocol)))
        opto_stim_trials = []
        opto_ctrl_trials = []
    else:
        baseline_trials = list(range(0, start_idx))
        opto_stim_trials = list(range(start_idx, end_idx, 3))
        opto_ctrl_trials = [trial_idx + 2 for trial_idx in opto_stim_trials]

    return optogenetic_protocol, baseline_trials, opto_stim_trials, opto_ctrl_trials

def _build_condition_dicts(curr_txt, baseline_trials, opto_stim_trials, opto_ctrl_trials, exclude_pulse_times=False):
    baseline_dict = {
        key: [value[idx] for idx in baseline_trials if len(value) > 1]
        for key, value in curr_txt.items()
        if not (exclude_pulse_times and key == 'pulse_times')
    }
    opto_stim_dict = {
        key: [value[idx] for idx in opto_stim_trials if len(value) > 1]
        for key, value in curr_txt.items()
        if not (exclude_pulse_times and key == 'pulse_times')
    }
    opto_ctrl_dict = {
        key: [value[idx] for idx in opto_ctrl_trials if len(value) > 1]
        for key, value in curr_txt.items()
        if not (exclude_pulse_times and key == 'pulse_times')
    }
    return baseline_dict, opto_stim_dict, opto_ctrl_dict

def _extract_lick_times(lick_trials, run_onsets, scalar_licks=False):
    aligned_licks = []
    for trial_idx, licks in enumerate(lick_trials):
        start_time = float(np.asarray(run_onsets[trial_idx]).squeeze())
        lick_times = [
            float(np.asarray(lick if scalar_licks else lick[0]).squeeze())
            for lick in licks
        ]
        aligned_licks.append(
            [
                (lick_time - start_time) / 1000
                for lick_time in lick_times
                if lick_time > start_time + 1000
            ]
        )
    return aligned_licks

def _summarise_standard_session(curr_txt, sessname, exclude_pulse_times=False, include_accel=False):
    _, baseline_trials, opto_stim_trials, opto_ctrl_trials = _segment_opto_session(curr_txt, sessname)
    baseline_dict, opto_stim_dict, _ = _build_condition_dicts(
        curr_txt,
        baseline_trials,
        opto_stim_trials,
        opto_ctrl_trials,
        exclude_pulse_times=exclude_pulse_times,
    )

    ctrl_lick_times = _extract_lick_times(
        baseline_dict['lick_times'],
        baseline_dict['run_onsets'],
        scalar_licks=False,
    )
    stim_lick_times = _extract_lick_times(
        opto_stim_dict['lick_times'],
        opto_stim_dict['run_onsets'],
        scalar_licks=False,
    )
    ctrl_lick_distances = [
        [lick for lick in licks if lick > 30]
        for licks in baseline_dict['lick_distances_aligned']
    ]
    stim_lick_distances = [
        [lick for lick in licks if lick > 30]
        for licks in opto_stim_dict['lick_distances_aligned']
    ]

    ctrl_mean_speeds = []
    for trial_idx in baseline_trials:
        speed_times = curr_txt['speed_times_aligned'][trial_idx]
        if speed_times:
            ctrl_mean_speeds.append(np.mean([sample[1] for sample in speed_times]))
    stim_mean_speeds = []
    for trial_idx in opto_stim_trials:
        speed_times = curr_txt['speed_times_aligned'][trial_idx]
        if speed_times:
            stim_mean_speeds.append(np.mean([sample[1] for sample in speed_times]))

    ctrl_rewarded = [
        not np.isnan(curr_txt['reward_times'][trial_idx])
        for trial_idx in baseline_trials[1:]
    ]
    stim_rewarded = [
        not np.isnan(curr_txt['reward_times'][trial_idx])
        for trial_idx in opto_stim_trials
    ]

    summary = {
        'ctrl_first_lick_time': np.median([licks[0] for licks in ctrl_lick_times if licks]),
        'stim_first_lick_time': np.median([licks[0] for licks in stim_lick_times if licks]),
        'ctrl_first_lick_distance': np.median([licks[0] for licks in ctrl_lick_distances if licks]),
        'stim_first_lick_distance': np.median([licks[0] for licks in stim_lick_distances if licks]),
        'ctrl_mean_speed': np.mean(ctrl_mean_speeds),
        'stim_mean_speed': np.mean(stim_mean_speeds),
        'ctrl_reward_perc': sum(ctrl_rewarded) / len(ctrl_rewarded),
        'stim_reward_perc': sum(stim_rewarded) / len(stim_rewarded),
    }

    if include_accel:
        ctrl_mean_accels = []
        for trial_idx in baseline_trials:
            speed_times = curr_txt['speed_times_aligned'][trial_idx]
            if speed_times:
                speeds = [sample[1] for sample in speed_times]
                ctrl_mean_accels.append(np.mean(np.diff(speeds[:500])) * 1_000)
        stim_mean_accels = []
        for trial_idx in opto_stim_trials:
            speed_times = curr_txt['speed_times_aligned'][trial_idx]
            if speed_times:
                speeds = [sample[1] for sample in speed_times]
                stim_mean_accels.append(np.mean(np.diff(speeds[:500])) * 1_000)
        summary['ctrl_mean_accel'] = np.mean(ctrl_mean_accels)
        summary['stim_mean_accel'] = np.mean(stim_mean_accels)

    return summary

def _summarise_opto_condition(path_groups, stim_code, cond_signature=None, include_accel=False):
    ctrl_lick_times = []
    stim_lick_times = []
    ctrl_lick_distances = []
    stim_lick_distances = []
    ctrl_mean_speeds = []
    stim_mean_speeds = []
    ctrl_mean_accels = []
    stim_mean_accels = []
    ctrl_perc_rew = []
    stim_perc_rew = []

    animals = set()
    n_sess = 0

    for group_name, pathnames in path_groups:
        for idx, pathname in enumerate(pathnames):
            if group_name == 'pathLCBehopt':
                curr_cond = rec_list.condLCBehopt[idx]
                curr_sess = rec_list.sessLCBehopt[idx]
                if str(cond_signature)[1:-1] not in str(curr_cond)[1:-1]:
                    continue
                sessname = pathname[-13:]
                file_idx = curr_sess[curr_cond.index(stim_code)]
                text_file = rf'{pathname}\{sessname}-0{file_idx}\{sessname}-0{file_idx}T.txt'
                exclude_pulse_times = False
            else:
                sessname = pathname[-17:]
                text_file = rf'{pathname}\{sessname}T.txt'
                curr_txt = process_behavioural_data(text_file)
                optogenetic_protocol = [trial[15] for trial in curr_txt['trial_statements']]
                if str(stim_code) not in np.unique(optogenetic_protocol):
                    continue
                exclude_pulse_times = True
                summary = _summarise_standard_session(
                    curr_txt,
                    sessname,
                    exclude_pulse_times=exclude_pulse_times,
                    include_accel=include_accel,
                )
                print_session(sessname)
                animals.add(sessname.split('-')[0])
                n_sess -= -1

                if 2 < summary['ctrl_first_lick_time'] < 10 and 2 < summary['stim_first_lick_time'] < 10:
                    ctrl_lick_times.append(summary['ctrl_first_lick_time'])
                    stim_lick_times.append(summary['stim_first_lick_time'])

                if summary['ctrl_first_lick_distance'] > 100 and summary['stim_first_lick_distance'] > 100:
                    ctrl_lick_distances.append(summary['ctrl_first_lick_distance'])
                    stim_lick_distances.append(summary['stim_first_lick_distance'])

                ctrl_mean_speeds.append(summary['ctrl_mean_speed'])
                stim_mean_speeds.append(summary['stim_mean_speed'])
                ctrl_perc_rew.append(summary['ctrl_reward_perc'])
                stim_perc_rew.append(summary['stim_reward_perc'])

                if include_accel:
                    ctrl_mean_accels.append(summary['ctrl_mean_accel'])
                    stim_mean_accels.append(summary['stim_mean_accel'])
                continue

            print_session(sessname)
            animals.add(sessname.split('-')[0])
            n_sess -= -1
            curr_txt = process_behavioural_data(text_file)
            summary = _summarise_standard_session(
                curr_txt,
                sessname,
                exclude_pulse_times=exclude_pulse_times,
                include_accel=include_accel,
            )

            if 2 < summary['ctrl_first_lick_time'] < 10 and 2 < summary['stim_first_lick_time'] < 10:
                ctrl_lick_times.append(summary['ctrl_first_lick_time'])
                stim_lick_times.append(summary['stim_first_lick_time'])

            if summary['ctrl_first_lick_distance'] > 100 and summary['stim_first_lick_distance'] > 100:
                ctrl_lick_distances.append(summary['ctrl_first_lick_distance'])
                stim_lick_distances.append(summary['stim_first_lick_distance'])

            ctrl_mean_speeds.append(summary['ctrl_mean_speed'])
            stim_mean_speeds.append(summary['stim_mean_speed'])
            ctrl_perc_rew.append(summary['ctrl_reward_perc'])
            stim_perc_rew.append(summary['stim_reward_perc'])

            if include_accel:
                ctrl_mean_accels.append(summary['ctrl_mean_accel'])
                stim_mean_accels.append(summary['stim_mean_accel'])

    results = {
        'ctrl_lick_times': ctrl_lick_times,
        'stim_lick_times': stim_lick_times,
        'ctrl_lick_distances': ctrl_lick_distances,
        'stim_lick_distances': stim_lick_distances,
        'ctrl_mean_speeds': ctrl_mean_speeds,
        'stim_mean_speeds': stim_mean_speeds,
        'ctrl_perc_rew': ctrl_perc_rew,
        'stim_perc_rew': stim_perc_rew,
        'animals': animals,
        'n_sess': n_sess,
    }

    if include_accel:
        results['ctrl_mean_accels'] = ctrl_mean_accels
        results['stim_mean_accels'] = stim_mean_accels

    return results


#%% main

def main():

    #%% recording list
    pathLCBehopt = rec_list.pathLCBehopt
    pathLCopt = rec_list.pathLCopt
    pathHPCopt = rec_list.pathHPCLCopt
    pathHPCLCtermopt = rec_list.pathHPCLCtermopt_beh

    path_groups = [
        ('pathLCBehopt', pathLCBehopt),
        ('pathLCopt', pathLCopt),
        ('pathHPCopt', pathHPCopt),
    ]

    #%% 020
    results_020 = _summarise_opto_condition(
        path_groups,
        stim_code=2,
        cond_signature=[0, 2, 0],
        include_accel=True,
    )

    ctrl_lick_times_020 = results_020['ctrl_lick_times']
    stim_lick_times_020 = results_020['stim_lick_times']
    ctrl_lick_distances_020 = results_020['ctrl_lick_distances']
    stim_lick_distances_020 = results_020['stim_lick_distances']
    ctrl_mean_speeds_020 = results_020['ctrl_mean_speeds']
    stim_mean_speeds_020 = results_020['stim_mean_speeds']
    ctrl_mean_accels_020 = results_020['ctrl_mean_accels']
    stim_mean_accels_020 = results_020['stim_mean_accels']
    ctrl_perc_rew_020 = results_020['ctrl_perc_rew']
    stim_perc_rew_020 = results_020['stim_perc_rew']
    n_sess = results_020['n_sess']

    pf.plot_violin_with_scatter(
        ctrl_lick_times_020,
        stim_lick_times_020,
        'grey',
        'royalblue',
        xticklabels=['Ctrl.', 'Stim.'],
        ylabel='First-lick time (s)',
        print_statistics=True,
        save=True,
        savepath=str(OPTO_OUTPUT_ROOT / '020_lick_times'),
        dpi=300,
        show=False,
        close=True,
    )

    pf.plot_violin_with_scatter(
        ctrl_lick_distances_020,
        stim_lick_distances_020,
        'grey',
        'royalblue',
        xticklabels=['Ctrl.', 'Stim.'],
        ylabel='First-lick distance (cm)',
        print_statistics=True,
        save=True,
        savepath=str(OPTO_OUTPUT_ROOT / '020_lick_distances'),
        dpi=300,
        show=False,
        close=True,
    )

    pf.plot_violin_with_scatter(
        ctrl_mean_speeds_020,
        stim_mean_speeds_020,
        'grey',
        'royalblue',
        xticklabels=['Ctrl.', 'Stim.'],
        ylabel='Mean speed (cm/s)',
        print_statistics=True,
        save=True,
        savepath=str(OPTO_OUTPUT_ROOT / '020_mean_speed'),
        dpi=300,
        show=False,
        close=True,
    )

    pf.plot_violin_with_scatter(
        ctrl_mean_accels_020,
        stim_mean_accels_020,
        'grey',
        'royalblue',
        xticklabels=['Ctrl.', 'Stim.'],
        ylabel='Mean accel. (cm/s2)',
        print_statistics=True,
        save=True,
        savepath=str(OPTO_OUTPUT_ROOT / '020_mean_accel'),
        dpi=300,
        show=False,
        close=True,
    )

    rmv = [
        idx for idx in range(len(ctrl_perc_rew_020))
        if ctrl_perc_rew_020[idx] < 0.5 or stim_perc_rew_020[idx] < 0.5
    ]
    ctrl_perc_rew_020_cln = [
        value for idx, value in enumerate(ctrl_perc_rew_020) if idx not in rmv
    ]
    stim_perc_rew_020_cln = [
        value for idx, value in enumerate(stim_perc_rew_020) if idx not in rmv
    ]

    pf.plot_violin_with_scatter(
        ctrl_perc_rew_020_cln,
        stim_perc_rew_020_cln,
        'grey',
        'royalblue',
        xticklabels=['Ctrl.', 'Stim.'],
        ylabel='Reward perc.',
        ylim=(.5, 1.02),
        print_statistics=True,
        save=True,
        savepath=str(OPTO_OUTPUT_ROOT / '020_perc_rew'),
        dpi=300,
        show=False,
        close=True,
    )

    print(f'n_sess = {n_sess}')

    #%% 030
    results_030 = _summarise_opto_condition(
        path_groups,
        stim_code=3,
        cond_signature=[0, 3, 0],
        include_accel=False,
    )

    ctrl_lick_times_030 = results_030['ctrl_lick_times']
    stim_lick_times_030 = results_030['stim_lick_times']
    ctrl_lick_distances_030 = results_030['ctrl_lick_distances']
    stim_lick_distances_030 = results_030['stim_lick_distances']
    ctrl_mean_speeds_030 = results_030['ctrl_mean_speeds']
    stim_mean_speeds_030 = results_030['stim_mean_speeds']
    ctrl_perc_rew_030 = results_030['ctrl_perc_rew']
    stim_perc_rew_030 = results_030['stim_perc_rew']
    n_sess = results_030['n_sess']

    pf.plot_violin_with_scatter(
        ctrl_lick_times_030,
        stim_lick_times_030,
        'grey',
        'royalblue',
        xticklabels=['Ctrl.', 'Stim.'],
        ylabel='First-lick time (s)',
        print_statistics=True,
        save=True,
        savepath=str(OPTO_OUTPUT_ROOT / '030_lick_times'),
        dpi=300,
        show=False,
        close=True,
    )

    pf.plot_violin_with_scatter(
        ctrl_lick_distances_030,
        stim_lick_distances_030,
        'grey',
        'royalblue',
        xticklabels=['Ctrl.', 'Stim.'],
        ylabel='First-lick distance (cm)',
        print_statistics=True,
        save=True,
        savepath=str(OPTO_OUTPUT_ROOT / '030_lick_distances'),
        dpi=300,
        show=False,
        close=True,
    )

    pf.plot_violin_with_scatter(
        ctrl_mean_speeds_030,
        stim_mean_speeds_030,
        'grey',
        'royalblue',
        xticklabels=['Ctrl.', 'Stim.'],
        ylabel='Mean speed (cm/s)',
        print_statistics=True,
        save=True,
        savepath=str(OPTO_OUTPUT_ROOT / '030_mean_speed'),
        dpi=300,
        show=False,
        close=True,
    )

    rmv = [
        idx for idx in range(len(ctrl_perc_rew_030))
        if ctrl_perc_rew_030[idx] < 0.6 or stim_perc_rew_030[idx] < 0.6
    ]
    ctrl_perc_rew_030 = [
        value for idx, value in enumerate(ctrl_perc_rew_030) if idx not in rmv
    ]
    stim_perc_rew_030 = [
        value for idx, value in enumerate(stim_perc_rew_030) if idx not in rmv
    ]

    pf.plot_violin_with_scatter(
        ctrl_perc_rew_030,
        stim_perc_rew_030,
        'grey',
        'royalblue',
        xticklabels=['Ctrl.', 'Stim.'],
        ylabel='Reward perc.',
        ylim=(.5, 1.03),
        print_statistics=True,
        save=True,
        savepath=str(OPTO_OUTPUT_ROOT / '030_perc_rew'),
        dpi=300,
        show=False,
        close=True,
    )

    print(f'n_sess = {n_sess}')

    #%% 040
    all_opto_ctrl_num_licks_aft = []
    all_opto_stim_num_licks_aft = []
    n_sess = 0

    for idx, pathname in enumerate(pathLCBehopt):
        sessname = pathname[-13:]
        curr_cond = rec_list.condLCBehopt[idx]
        curr_sess = rec_list.sessLCBehopt[idx]

        if str([0, 4])[1:-1] not in str(curr_cond)[1:-1]:
            continue

        print_session(sessname)
        n_sess -= -1

        file_idx = curr_sess[curr_cond.index(4)]
        text_file = rf'{pathname}\{sessname}-0{file_idx}\{sessname}-0{file_idx}T.txt'
        curr_txt = process_behavioural_data(text_file)
        _, baseline_trials, opto_stim_trials, opto_ctrl_trials = _segment_opto_session(
            curr_txt,
            sessname
        )
        _, opto_stim_dict, opto_ctrl_dict = _build_condition_dicts(
            curr_txt,
            baseline_trials,
            opto_stim_trials,
            opto_ctrl_trials,
            exclude_pulse_times=False,
        )

        opto_stim_num_licks_aft = []
        for trial_idx, licks in enumerate(opto_stim_dict['lick_times']):
            pump_time = float(
                np.asarray(opto_stim_dict['reward_times'][trial_idx]).squeeze()
            )
            n_licks = 0
            for lick in licks:
                if float(np.asarray(lick[0]).squeeze()) > pump_time:
                    n_licks += 1
            opto_stim_num_licks_aft.append(n_licks)

        opto_ctrl_num_licks_aft = []
        for trial_idx, licks in enumerate(opto_ctrl_dict['lick_times']):
            pump_time = float(
                np.asarray(opto_ctrl_dict['reward_times'][trial_idx]).squeeze()
            )
            n_licks = 0
            for lick in licks:
                if float(np.asarray(lick[0]).squeeze()) > pump_time:
                    n_licks += 1
            opto_ctrl_num_licks_aft.append(n_licks)

        all_opto_ctrl_num_licks_aft.append(np.mean(opto_ctrl_num_licks_aft))
        all_opto_stim_num_licks_aft.append(np.mean(opto_stim_num_licks_aft))

    pf.plot_violin_with_scatter(
        all_opto_ctrl_num_licks_aft,
        all_opto_stim_num_licks_aft,
        'grey',
        'royalblue',
        xticklabels=['ctrl.', 'stim.'],
        ylabel='num. licks aft. reward',
        print_statistics=True,
        save=True,
        savepath=str(OPTO_OUTPUT_ROOT / '040_lick_aft_rew'),
        dpi=300,
        show=False,
        close=True,
    )

    print(f'n_sess = {n_sess}')

    #%% 020 terminal
    ctrl_lick_times_020_term = []
    stim_lick_times_020_term = []

    for pathname in pathHPCLCtermopt:
        sessname = pathname[-17:]
        text_file = rf'{pathname}\{sessname}T.txt'
        curr_txt = process_behavioural_data(text_file)
        optogenetic_protocol, baseline_trials, opto_stim_trials, opto_ctrl_trials = _segment_opto_session(
            curr_txt,
            sessname
        )
        if '2' not in np.unique(optogenetic_protocol):
            continue

        _, opto_stim_dict, opto_ctrl_dict = _build_condition_dicts(
            curr_txt,
            baseline_trials,
            opto_stim_trials,
            opto_ctrl_trials,
            exclude_pulse_times=False,
        )
        ctrl_lick_times = _extract_lick_times(
            opto_ctrl_dict['lick_times'],
            opto_ctrl_dict['run_onsets'],
            scalar_licks=True,
        )
        stim_lick_times = _extract_lick_times(
            opto_stim_dict['lick_times'],
            opto_stim_dict['run_onsets'],
            scalar_licks=True,
        )
        ctrl_first_lick_time = np.median([licks[0] for licks in ctrl_lick_times if licks])
        stim_first_lick_time = np.median([licks[0] for licks in stim_lick_times if licks])
        if ctrl_first_lick_time < 10 and stim_first_lick_time < 10:
            ctrl_lick_times_020_term.append(ctrl_first_lick_time)
            stim_lick_times_020_term.append(stim_first_lick_time)

    pf.plot_violin_with_scatter(
        ctrl_lick_times_020_term,
        stim_lick_times_020_term,
        'grey',
        'royalblue',
        xticklabels=['ctrl.', 'stim.'],
        ylabel='first-lick time (s)',
        save=False,
        savepath=str(OPTO_OUTPUT_ROOT / '020_lick_times'),
        dpi=300,
        show=False,
        close=True,
    )

if __name__ == '__main__':
    main()
