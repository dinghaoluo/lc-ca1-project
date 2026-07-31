# -*- coding: utf-8 -*-
'''
Created on 26 February 2026
Modified on 2 June 2026

summarise speed, lick, and reward changes after prazosin, propranolol, and
SCH23390

@author: Dinghao Luo
'''

#%% imports
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import sem, ttest_rel, wilcoxon

repo_root = Path(__file__).resolve().parents[2]
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))

from common_functions import mpl_formatting, replace_outlier, smooth_convolve
mpl_formatting()

from plotting_functions import plot_violin_with_scatter
from lick_time_utils import lick_time_map
import project_paths as pp
import behaviour_functions as bf

import rec_list


#%% global parameters
pharmacology_figures_stem = pp.PHARMACOLOGY_FIGURES_STEM
pharmacology_figures_stem.mkdir(parents=True, exist_ok=True)

tot_distance = 2200  # mm
XAXIS = np.arange(tot_distance) / 10  # cm
EARLY_LICK_TIME_WINDOW_S = (0.5, 2.0)
LATE_LICK_TIME_WINDOW_S = (2.0, 4.0)
FULL_LICK_TIME_WINDOW_S = (0.0, 4.0)
LICK_TIME_NORM_WINDOW_S = (0.0, 4.0)
LICK_DISTANCE_NORM_WINDOW_CM = (30.0, 219.0)
EARLY_LICK_DISTANCE_WINDOW_CM = (30.0, 100.0)
LICK_DISTANCE_STATS_WINDOW_CM = (120.0, 180.0)

def paired_stats(experiment_name, comparison_name, baseline_values, drug_values,
                 drug_label):
    baseline_values = np.asarray(baseline_values, dtype=float).ravel()
    drug_values = np.asarray(drug_values, dtype=float).ravel()
    finite = np.isfinite(baseline_values) & np.isfinite(drug_values)
    baseline_values = baseline_values[finite]
    drug_values = drug_values[finite]
    baseline_summary = {
        'n': len(baseline_values),
        'median': np.nanmedian(baseline_values),
        'q25': np.nanpercentile(baseline_values, 25),
        'q75': np.nanpercentile(baseline_values, 75),
        'mean': np.nanmean(baseline_values),
        'sem': sem(baseline_values, nan_policy='omit'),
    }
    drug_summary = {
        'n': len(drug_values),
        'median': np.nanmedian(drug_values),
        'q25': np.nanpercentile(drug_values, 25),
        'q75': np.nanpercentile(drug_values, 75),
        'mean': np.nanmean(drug_values),
        'sem': sem(drug_values, nan_policy='omit'),
    }
    wilc_stat, wilc_p = wilcoxon(baseline_values, drug_values)
    ttest_stat, ttest_p = ttest_rel(baseline_values, drug_values)
    wilc_p_text = f'{wilc_p:.2e}'
    ttest_p_text = f'{ttest_p:.2e}'

    print(f'\n[{experiment_name}] {comparison_name}: Baseline vs {drug_label}')
    print(
        '  Baseline: '
        f'n={baseline_summary["n"]}, '
        f'median={baseline_summary["median"]:.4g}, '
        f'IQR=[{baseline_summary["q25"]:.4g}, {baseline_summary["q75"]:.4g}], '
        f'mean+/-SEM={baseline_summary["mean"]:.4g}+/-{baseline_summary["sem"]:.4g}'
    )
    print(
        f'  {drug_label}: '
        f'n={drug_summary["n"]}, '
        f'median={drug_summary["median"]:.4g}, '
        f'IQR=[{drug_summary["q25"]:.4g}, {drug_summary["q75"]:.4g}], '
        f'mean+/-SEM={drug_summary["mean"]:.4g}+/-{drug_summary["sem"]:.4g}'
    )
    print(f'  Wilcoxon paired: W={wilc_stat:.4g}, p={wilc_p_text}')
    print(f'  Paired t-test:   t={ttest_stat:.4g}, p={ttest_p_text}')

    return {
        'wilcoxon_stat': wilc_stat,
        'wilcoxon_p': wilc_p,
        'ttest_stat': ttest_stat,
        'ttest_p': ttest_p,
        'baseline_summary': baseline_summary,
        'drug_summary': drug_summary,
    }

def paired_animal_mean_values(baseline_values, drug_values, animal_ids):
    baseline_values = np.asarray(baseline_values, dtype=float)
    drug_values = np.asarray(drug_values, dtype=float)
    animal_ids = np.asarray(animal_ids)

    animals = list(dict.fromkeys(animal_ids.tolist()))
    baseline_by_animal = []
    drug_by_animal = []

    for animal in animals:
        mask = animal_ids == animal
        baseline_by_animal.append(np.nanmean(baseline_values[mask]))
        drug_by_animal.append(np.nanmean(drug_values[mask]))

    return np.asarray(baseline_by_animal), np.asarray(drug_by_animal), animals

def plot_profile_pair(
        x_axis,
        baseline_profiles,
        drug_profiles,
        color,
        drug_name,
        save_dir,
        filename,
        xlabel,
        ylabel,
        xlim,
        ylim=None
        ):
    baseline_profiles = np.asarray(baseline_profiles, dtype=float)
    drug_profiles = np.asarray(drug_profiles, dtype=float)

    mean_baseline = np.nanmean(baseline_profiles, axis=0)
    mean_drug = np.nanmean(drug_profiles, axis=0)
    sem_baseline = sem(baseline_profiles, axis=0, nan_policy='omit')
    sem_drug = sem(drug_profiles, axis=0, nan_policy='omit')

    fig, ax = plt.subplots(figsize=(2, 1.7))

    lp, = ax.plot(x_axis, mean_baseline, color='grey')
    ax.fill_between(
        x_axis,
        mean_baseline + sem_baseline,
        mean_baseline - sem_baseline,
        color='grey',
        alpha=.15,
        edgecolor='none'
    )

    ld, = ax.plot(x_axis, mean_drug, color=color)
    ax.fill_between(
        x_axis,
        mean_drug + sem_drug,
        mean_drug - sem_drug,
        color=color,
        alpha=.15,
        edgecolor='none'
    )

    if ylim is None:
        ylim = (
            0,
            np.nanmax(np.concatenate([
                mean_baseline + sem_baseline,
                mean_drug + sem_drug,
            ])) * 1.1
        )

    ax.legend(
        [lp, ld],
        ['Baseline', drug_name],
        frameon=False,
        loc='upper left',
        fontsize=7,
        handlelength=1.6,
        borderaxespad=0.2
    )
    ax.set(xlim=xlim, ylim=ylim, xlabel=xlabel, ylabel=ylabel)

    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)

    save_dir = Path(save_dir)
    save_dir.mkdir(parents=True, exist_ok=True)
    for ext in ['.png', '.pdf']:
        fig.savefig(save_dir / f'{filename}{ext}', dpi=300, bbox_inches='tight')
    plt.close(fig)


#%% core function
def summarise_drug(
        drug_name,
        paths,
        sessions,
        color,
        save_folder,
        sch_mode=False
    ):
    '''
    sch_mode=False -> baseline=session[0], drug=session[1]
    sch_mode=True  -> baseline=session==1, drug=mean(session 2 & 3)
    '''

    mean_speeds_baseline = []
    mean_speeds_drug     = []
    mean_licks_baseline  = []
    mean_licks_drug      = []
    mean_lick_times_baseline = []
    mean_lick_times_drug     = []
    norm_licks_baseline  = []
    norm_licks_drug      = []
    norm_lick_times_baseline = []
    norm_lick_times_drug     = []
    frac_licks_baseline  = []
    frac_licks_drug      = []
    reward_baseline      = []
    reward_drug          = []
    valid_run_baseline   = []
    valid_run_drug       = []
    animal_ids           = []

    for pathname, sesslist in zip(paths, sessions):

        subsessions = []
        speed_drugs = []
        lick_drugs  = []
        lick_time_drugs = []
        norm_lick_drugs = []
        norm_lick_time_drugs = []
        frac_lick_drugs = []
        reward_drugs = []
        valid_run_drugs = []
        for part in Path(str(pathname).replace('\\', '/')).parts:
            if part.startswith('ANMD'):
                animal_id = part.replace('ANMD', 'A')
                break
        else:
            animal_id = Path(pathname).name.split('-')[0]

        for idx, sess in enumerate(sesslist):

            # -------- build txt path --------
            if sch_mode:
                sessname = Path(pathname).name
                recname = f'{sessname}-0{sess}'
                txtpath = rf'{pathname}\{recname}\{recname}T.txt'
            else:
                suffix = f'-0{sess}'
                txtpath = f'{pathname}{suffix}T.txt'

            file = bf.process_behavioural_data(txtpath)

            # -------- reward --------
            reward_times = file['reward_times'][1:-1]
            rewarded = [1 if not np.isnan(t) else 0 for t in reward_times]
            reward_val = np.mean(rewarded)

            # -------- speed --------
            speed_dist = np.array([
                replace_outlier(np.array(trial))
                for trial in file['speed_distances_aligned']
                if len(trial) > 0
            ])
            speed_mean = np.mean(speed_dist, axis=0)

            # prazosin/propranolol multiply by 1.8
            if not sch_mode:
                speed_mean *= 1.8

            # -------- licks --------
            lick_dist = np.array([
                smooth_convolve(np.array(trial), sigma=10) * 10
                for trial in file['lick_maps']
                if len(trial) > 0
            ])
            lick_mean = np.mean(lick_dist, axis=0)

            # -------- licks by time --------
            # no-lick trials still contribute zeros, but trials without a run onset do not
            lick_times = [
                lick_trial
                for lick_trial, run_onset in zip(
                    file['lick_times_aligned'],
                    file['run_onsets']
                )
                if np.isfinite(float(run_onset))
            ]
            lick_time_mean = np.nanmean(lick_time_map(lick_times), axis=0)
            run_onsets = np.asarray(file['run_onsets'], dtype=float)
            valid_run_frac = np.mean(np.isfinite(run_onsets))

            lick_density = lick_mean / 10
            distance_bin_width = np.nanmedian(np.diff(XAXIS))
            distance_norm_mask = (
                (XAXIS >= LICK_DISTANCE_NORM_WINDOW_CM[0])
                & (XAXIS <= LICK_DISTANCE_NORM_WINDOW_CM[1])
            )
            frac_lick_mean = lick_density / np.nansum(
                lick_density[distance_norm_mask] * distance_bin_width
            )

            subsessions.append({
                'idx': idx,
                'sess': sess,
                'speed_mean': speed_mean,
                'lick_mean': lick_mean,
                'lick_time_mean': lick_time_mean,
                'frac_lick_mean': frac_lick_mean,
                'reward_val': reward_val,
                'valid_run_frac': valid_run_frac,
            })

        lick_time_axis = np.arange(subsessions[0]['lick_time_mean'].size) / 1000
        distance_norm_mask = (
            (XAXIS >= LICK_DISTANCE_NORM_WINDOW_CM[0])
            & (XAXIS <= LICK_DISTANCE_NORM_WINDOW_CM[1])
        )
        shared_distance_peak = np.nanmax([
            np.nanmax((entry['lick_mean'] / 10)[distance_norm_mask])
            for entry in subsessions
        ])
        time_norm_mask = (
            (lick_time_axis >= LICK_TIME_NORM_WINDOW_S[0])
            & (lick_time_axis <= LICK_TIME_NORM_WINDOW_S[1])
        )
        shared_time_peak = np.nanmax([
            np.nanmax(entry['lick_time_mean'][time_norm_mask])
            for entry in subsessions
        ])

        for entry in subsessions:
            sess = entry['sess']
            idx = entry['idx']
            speed_mean = entry['speed_mean']
            lick_mean = entry['lick_mean']
            lick_time_mean = entry['lick_time_mean']
            norm_lick_mean = (lick_mean / 10) / shared_distance_peak
            norm_lick_time_mean = lick_time_mean / shared_time_peak
            frac_lick_mean = entry['frac_lick_mean']
            reward_val = entry['reward_val']
            valid_run_frac = entry['valid_run_frac']

            # -------- grouping logic --------
            if sch_mode:
                if sess == 1:
                    mean_speeds_baseline.append(speed_mean)
                    mean_licks_baseline.append(lick_mean)
                    mean_lick_times_baseline.append(lick_time_mean)
                    norm_licks_baseline.append(norm_lick_mean)
                    norm_lick_times_baseline.append(norm_lick_time_mean)
                    frac_licks_baseline.append(frac_lick_mean)
                    reward_baseline.append(reward_val)
                    valid_run_baseline.append(valid_run_frac)
                    animal_ids.append(animal_id)
                elif sess in [2, 3]:
                    speed_drugs.append(speed_mean)
                    lick_drugs.append(lick_mean)
                    lick_time_drugs.append(lick_time_mean)
                    norm_lick_drugs.append(norm_lick_mean)
                    norm_lick_time_drugs.append(norm_lick_time_mean)
                    frac_lick_drugs.append(frac_lick_mean)
                    reward_drugs.append(reward_val)
                    valid_run_drugs.append(valid_run_frac)
            else:
                if idx == 0:
                    mean_speeds_baseline.append(speed_mean)
                    mean_licks_baseline.append(lick_mean)
                    mean_lick_times_baseline.append(lick_time_mean)
                    norm_licks_baseline.append(norm_lick_mean)
                    norm_lick_times_baseline.append(norm_lick_time_mean)
                    frac_licks_baseline.append(frac_lick_mean)
                    reward_baseline.append(reward_val)
                    valid_run_baseline.append(valid_run_frac)
                    animal_ids.append(animal_id)
                elif idx == 1:
                    mean_speeds_drug.append(speed_mean)
                    mean_licks_drug.append(lick_mean)
                    mean_lick_times_drug.append(lick_time_mean)
                    norm_licks_drug.append(norm_lick_mean)
                    norm_lick_times_drug.append(norm_lick_time_mean)
                    frac_licks_drug.append(frac_lick_mean)
                    reward_drug.append(reward_val)
                    valid_run_drug.append(valid_run_frac)

        if sch_mode:
            mean_speeds_drug.append(np.mean(speed_drugs, axis=0))
            mean_licks_drug.append(np.mean(lick_drugs, axis=0))
            mean_lick_times_drug.append(np.mean(lick_time_drugs, axis=0))
            norm_licks_drug.append(np.nanmean(norm_lick_drugs, axis=0))
            norm_lick_times_drug.append(
                np.nanmean(norm_lick_time_drugs, axis=0)
            )
            frac_licks_drug.append(np.nanmean(frac_lick_drugs, axis=0))
            reward_drug.append(np.mean(reward_drugs))
            valid_run_drug.append(np.mean(valid_run_drugs))

    # convert to arrays
    mean_speeds_baseline = np.array(mean_speeds_baseline)
    mean_speeds_drug     = np.array(mean_speeds_drug)
    mean_licks_baseline  = np.array(mean_licks_baseline)
    mean_licks_drug      = np.array(mean_licks_drug)
    mean_lick_times_baseline = np.array(mean_lick_times_baseline)
    mean_lick_times_drug     = np.array(mean_lick_times_drug)
    norm_licks_baseline  = np.array(norm_licks_baseline)
    norm_licks_drug      = np.array(norm_licks_drug)
    norm_lick_times_baseline = np.array(norm_lick_times_baseline)
    norm_lick_times_drug     = np.array(norm_lick_times_drug)
    frac_licks_baseline  = np.array(frac_licks_baseline)
    frac_licks_drug      = np.array(frac_licks_drug)
    reward_baseline      = np.array(reward_baseline)
    reward_drug          = np.array(reward_drug)
    valid_run_baseline   = np.array(valid_run_baseline)
    valid_run_drug       = np.array(valid_run_drug)
    animal_ids           = np.array(animal_ids)

    save_dir = pharmacology_figures_stem / save_folder
    save_dir.mkdir(parents=True, exist_ok=True)

    early_time_slice = slice(
        int(EARLY_LICK_TIME_WINDOW_S[0] * 1000),
        int(EARLY_LICK_TIME_WINDOW_S[1] * 1000),
    )
    late_time_slice = slice(
        int(LATE_LICK_TIME_WINDOW_S[0] * 1000),
        int(LATE_LICK_TIME_WINDOW_S[1] * 1000),
    )
    full_time_slice = slice(
        int(FULL_LICK_TIME_WINDOW_S[0] * 1000),
        int(FULL_LICK_TIME_WINDOW_S[1] * 1000),
    )
    early_distance_mask = (
        (XAXIS >= EARLY_LICK_DISTANCE_WINDOW_CM[0]) &
        (XAXIS <= EARLY_LICK_DISTANCE_WINDOW_CM[1])
    )
    distance_stats_mask = (
        (XAXIS >= LICK_DISTANCE_STATS_WINDOW_CM[0]) &
        (XAXIS <= LICK_DISTANCE_STATS_WINDOW_CM[1])
    )
    distance_norm_mask = (
        (XAXIS >= LICK_DISTANCE_NORM_WINDOW_CM[0]) &
        (XAXIS <= LICK_DISTANCE_NORM_WINDOW_CM[1])
    )
    distance_bin_width = np.nanmedian(np.diff(XAXIS))

    early_lick_time_baseline = np.nanmean(
        mean_lick_times_baseline[:, early_time_slice],
        axis=1,
    )
    early_lick_time_drug = np.nanmean(
        mean_lick_times_drug[:, early_time_slice],
        axis=1,
    )
    early_lick_time_stats = paired_stats(
        drug_name,
        'early time-based lick rate, 0.5-2 s from run onset',
        early_lick_time_baseline,
        early_lick_time_drug,
        drug_name
    )
    late_lick_time_baseline = np.nanmean(
        mean_lick_times_baseline[:, late_time_slice],
        axis=1,
    )
    late_lick_time_drug = np.nanmean(
        mean_lick_times_drug[:, late_time_slice],
        axis=1,
    )
    paired_stats(
        drug_name,
        'late time-based lick rate, 2-4 s from run onset',
        late_lick_time_baseline,
        late_lick_time_drug,
        drug_name
    )
    full_lick_time_baseline = np.nanmean(
        mean_lick_times_baseline[:, full_time_slice],
        axis=1,
    )
    full_lick_time_drug = np.nanmean(
        mean_lick_times_drug[:, full_time_slice],
        axis=1,
    )
    paired_stats(
        drug_name,
        'time-based lick rate, 0-4 s from run onset',
        full_lick_time_baseline,
        full_lick_time_drug,
        drug_name
    )
    peak_lick_time_baseline = np.nanmax(
        mean_lick_times_baseline[:, full_time_slice],
        axis=1,
    )
    peak_lick_time_drug = np.nanmax(
        mean_lick_times_drug[:, full_time_slice],
        axis=1,
    )
    paired_stats(
        drug_name,
        'peak time-based lick rate, 0-4 s from run onset',
        peak_lick_time_baseline,
        peak_lick_time_drug,
        drug_name
    )
    paired_stats(
        drug_name,
        'finite run-onset trial fraction',
        valid_run_baseline,
        valid_run_drug,
        drug_name
    )
    lick_density_baseline = mean_licks_baseline / 10
    lick_density_drug = mean_licks_drug / 10
    distance_lick_baseline = np.nanmean(
        lick_density_baseline[:, distance_stats_mask],
        axis=1,
    )
    distance_lick_drug = np.nanmean(
        lick_density_drug[:, distance_stats_mask],
        axis=1,
    )
    paired_stats(
        drug_name,
        'distance-based lick density, 120-180 cm',
        distance_lick_baseline,
        distance_lick_drug,
        drug_name
    )
    peak_distance_lick_baseline = np.nanmax(
        lick_density_baseline[:, distance_norm_mask],
        axis=1,
    )
    peak_distance_lick_drug = np.nanmax(
        lick_density_drug[:, distance_norm_mask],
        axis=1,
    )
    paired_stats(
        drug_name,
        'peak distance-based lick density, 30-219 cm',
        peak_distance_lick_baseline,
        peak_distance_lick_drug,
        drug_name
    )

    # speed
    fig, ax = plt.subplots(figsize=(2, 1.7))

    ms_b = np.mean(mean_speeds_baseline, axis=0)
    ms_d = np.mean(mean_speeds_drug, axis=0)
    ss_b = sem(mean_speeds_baseline, axis=0)
    ss_d = sem(mean_speeds_drug, axis=0)

    ax.plot(XAXIS, ms_b, color='grey')
    ax.fill_between(XAXIS, ms_b+ss_b, ms_b-ss_b,
                    color='grey', alpha=.15, edgecolor='none')

    ax.plot(XAXIS, ms_d, color=color)
    ax.fill_between(XAXIS, ms_d+ss_d, ms_d-ss_d,
                    color=color, alpha=.15, edgecolor='none')

    ax.set(xlim=(0, 180), ylim=(0, 75),
           xlabel='Distance (cm)', ylabel='Speed (cm/s)')

    for s in ['top','right']:
        ax.spines[s].set_visible(False)

    fig.savefig(save_dir / 'speed_profile.png', dpi=300, bbox_inches='tight')
    fig.savefig(save_dir / 'speed_profile.pdf', dpi=300, bbox_inches='tight')

    # lick
    fig, ax = plt.subplots(figsize=(2, 1.7))

    ml_b = np.mean(mean_licks_baseline, axis=0) / 10
    ml_d = np.mean(mean_licks_drug, axis=0) / 10
    sl_b = sem(mean_licks_baseline, axis=0) / 10
    sl_d = sem(mean_licks_drug, axis=0) / 10

    ax.plot(XAXIS, ml_b, color='grey')
    ax.fill_between(XAXIS, ml_b+sl_b, ml_b-sl_b,
                    color='grey', alpha=.15, edgecolor='none')

    ax.plot(XAXIS, ml_d, color=color)
    ax.fill_between(XAXIS, ml_d+sl_d, ml_d-sl_d,
                    color=color, alpha=.15, edgecolor='none')

    ax.set(xlim=(30, 219), ylim=(0, 0.6),
           xlabel='Distance (cm)', ylabel='Lick density (count/cm)')

    for s in ['top','right']:
        ax.spines[s].set_visible(False)

    fig.savefig(save_dir / 'lick_profile.png', dpi=300, bbox_inches='tight')
    fig.savefig(save_dir / 'lick_profile.pdf', dpi=300, bbox_inches='tight')

    frac_early_distance_lick_baseline = np.nansum(
        frac_licks_baseline[:, early_distance_mask] * distance_bin_width,
        axis=1,
    )
    frac_early_distance_lick_drug = np.nansum(
        frac_licks_drug[:, early_distance_mask] * distance_bin_width,
        axis=1,
    )
    paired_stats(
        drug_name,
        'fraction of licks, 30-100 cm',
        frac_early_distance_lick_baseline,
        frac_early_distance_lick_drug,
        drug_name
    )
    frac_distance_lick_baseline = np.nansum(
        frac_licks_baseline[:, distance_stats_mask] * distance_bin_width,
        axis=1,
    )
    frac_distance_lick_drug = np.nansum(
        frac_licks_drug[:, distance_stats_mask] * distance_bin_width,
        axis=1,
    )
    paired_stats(
        drug_name,
        'fraction of licks, 120-180 cm',
        frac_distance_lick_baseline,
        frac_distance_lick_drug,
        drug_name
    )
    plot_profile_pair(
        XAXIS,
        frac_licks_baseline,
        frac_licks_drug,
        color,
        drug_name,
        save_dir,
        'lick_profile_fraction',
        xlabel='Distance (cm)',
        ylabel='Fraction of licks / cm',
        xlim=(30, 219)
    )

    plot_violin_with_scatter(
        frac_early_distance_lick_baseline,
        frac_early_distance_lick_drug,
        'grey',
        color,
        xticklabels=['Baseline', drug_name],
        ylabel='Fraction of licks',
        title='30-100 cm',
        figsize=(2.0, 2.0),
        showmainline=False,
        print_statistics=False,
        plot_statistics=True,
        save=True,
        savepath=save_dir / 'lick_distance_early_30_100cm_fraction'
    )

    plot_violin_with_scatter(
        frac_distance_lick_baseline,
        frac_distance_lick_drug,
        'grey',
        color,
        xticklabels=['Baseline', drug_name],
        ylabel='Fraction of licks',
        title='120-180 cm',
        figsize=(2.0, 2.0),
        showmainline=False,
        print_statistics=False,
        plot_statistics=True,
        save=True,
        savepath=save_dir / 'lick_distance_120_180cm_fraction'
    )

    norm_early_distance_lick_baseline = np.nanmean(
        norm_licks_baseline[:, early_distance_mask],
        axis=1,
    )
    norm_early_distance_lick_drug = np.nanmean(
        norm_licks_drug[:, early_distance_mask],
        axis=1,
    )
    paired_stats(
        drug_name,
        'normalised early distance lick density, 30-100 cm',
        norm_early_distance_lick_baseline,
        norm_early_distance_lick_drug,
        drug_name
    )
    norm_distance_lick_baseline = np.nanmean(
        norm_licks_baseline[:, distance_stats_mask],
        axis=1,
    )
    norm_distance_lick_drug = np.nanmean(
        norm_licks_drug[:, distance_stats_mask],
        axis=1,
    )
    paired_stats(
        drug_name,
        'normalised distance lick density, 120-180 cm',
        norm_distance_lick_baseline,
        norm_distance_lick_drug,
        drug_name
    )
    plot_profile_pair(
        XAXIS,
        norm_licks_baseline,
        norm_licks_drug,
        color,
        drug_name,
        save_dir,
        'lick_profile_norm',
        xlabel='Distance (cm)',
        ylabel='Norm. lick density',
        xlim=(30, 219),
        ylim=(0, 1.1)
    )

    plot_violin_with_scatter(
        norm_early_distance_lick_baseline,
        norm_early_distance_lick_drug,
        'grey',
        color,
        xticklabels=['Baseline', drug_name],
        ylabel='Norm. lick density',
        title='30-100 cm',
        figsize=(2.0, 2.0),
        showmainline=False,
        print_statistics=False,
        plot_statistics=True,
        save=True,
        savepath=save_dir / 'lick_distance_early_30_100cm_norm'
    )

    plot_violin_with_scatter(
        norm_distance_lick_baseline,
        norm_distance_lick_drug,
        'grey',
        color,
        xticklabels=['Baseline', drug_name],
        ylabel='Norm. lick density',
        title='120-180 cm',
        figsize=(2.0, 2.0),
        showmainline=False,
        print_statistics=False,
        plot_statistics=True,
        save=True,
        savepath=save_dir / 'lick_distance_120_180cm_norm'
    )

    # lick by time
    early_wilc_p = early_lick_time_stats['wilcoxon_p']
    early_ttest_p = early_lick_time_stats['ttest_p']
    early_wilc_p_text = f'{early_wilc_p:.2e}'
    early_ttest_p_text = f'{early_ttest_p:.2e}'
    early_lick_time_stats_text = (
        f'{EARLY_LICK_TIME_WINDOW_S[0]:g}-{EARLY_LICK_TIME_WINDOW_S[1]:g} s\n'
        f'wilc p={early_wilc_p_text}\n'
        f't-test p={early_ttest_p_text}'
    )
    time_axis = np.arange(mean_lick_times_baseline.shape[1]) / 1000
    mean_baseline = np.nanmean(mean_lick_times_baseline, axis=0)
    mean_drug = np.nanmean(mean_lick_times_drug, axis=0)
    sem_baseline = sem(mean_lick_times_baseline, axis=0, nan_policy='omit')
    sem_drug = sem(mean_lick_times_drug, axis=0, nan_policy='omit')

    fig, ax = plt.subplots(figsize=(2, 1.7))

    lp, = ax.plot(time_axis, mean_baseline, color='grey')
    ax.fill_between(
        time_axis,
        mean_baseline + sem_baseline,
        mean_baseline - sem_baseline,
        color='grey',
        alpha=.15,
        edgecolor='none'
    )

    ld, = ax.plot(time_axis, mean_drug, color=color)
    ax.fill_between(
        time_axis,
        mean_drug + sem_drug,
        mean_drug - sem_drug,
        color=color,
        alpha=.15,
        edgecolor='none'
    )

    ax.axvspan(
        EARLY_LICK_TIME_WINDOW_S[0],
        EARLY_LICK_TIME_WINDOW_S[1],
        facecolor='0.85',
        alpha=0.25,
        edgecolor='none',
        linewidth=0,
        zorder=0
    )
    ax.text(
        0.98,
        0.98,
        early_lick_time_stats_text,
        transform=ax.transAxes,
        ha='right',
        va='top',
        fontsize=6
    )
    ax.legend(
        [lp, ld],
        ['Baseline', drug_name],
        frameon=False,
        loc='upper left',
        fontsize=7,
        handlelength=1.6,
        borderaxespad=0.2
    )
    ax.set(
        xlim=(0, 4),
        ylim=(
            0,
            np.nanmax(np.concatenate([
                mean_baseline + sem_baseline,
                mean_drug + sem_drug,
            ])) * 1.1
        ),
        xlabel='Time from run-onset (s)',
        ylabel='Lick rate (Hz)'
    )

    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)

    fig.savefig(save_dir / 'lick_time_profile.png', dpi=300,
                bbox_inches='tight')
    fig.savefig(save_dir / 'lick_time_profile.pdf', dpi=300,
                bbox_inches='tight')
    plt.close(fig)

    norm_early_lick_time_baseline = np.nanmean(
        norm_lick_times_baseline[:, early_time_slice],
        axis=1,
    )
    norm_early_lick_time_drug = np.nanmean(
        norm_lick_times_drug[:, early_time_slice],
        axis=1,
    )
    paired_stats(
        drug_name,
        'normalised early lick rate, 0.5-2 s from run onset',
        norm_early_lick_time_baseline,
        norm_early_lick_time_drug,
        drug_name
    )
    norm_late_lick_time_baseline = np.nanmean(
        norm_lick_times_baseline[:, late_time_slice],
        axis=1,
    )
    norm_late_lick_time_drug = np.nanmean(
        norm_lick_times_drug[:, late_time_slice],
        axis=1,
    )
    paired_stats(
        drug_name,
        'normalised late lick rate, 2-4 s from run onset',
        norm_late_lick_time_baseline,
        norm_late_lick_time_drug,
        drug_name
    )
    norm_full_lick_time_baseline = np.nanmean(
        norm_lick_times_baseline[:, full_time_slice],
        axis=1,
    )
    norm_full_lick_time_drug = np.nanmean(
        norm_lick_times_drug[:, full_time_slice],
        axis=1,
    )
    paired_stats(
        drug_name,
        'normalised lick rate, 0-4 s from run onset',
        norm_full_lick_time_baseline,
        norm_full_lick_time_drug,
        drug_name
    )
    plot_profile_pair(
        time_axis,
        norm_lick_times_baseline,
        norm_lick_times_drug,
        color,
        drug_name,
        save_dir,
        'lick_time_profile_norm',
        xlabel='Time from run-onset (s)',
        ylabel='Norm. lick rate',
        xlim=(0, 4),
        ylim=(0, 1.1)
    )

    plot_violin_with_scatter(
        norm_early_lick_time_baseline,
        norm_early_lick_time_drug,
        'grey',
        color,
        xticklabels=['Baseline', drug_name],
        ylabel='Norm. lick rate',
        title='0.5-2 s',
        figsize=(2.0, 2.0),
        showmainline=False,
        print_statistics=False,
        plot_statistics=True,
        save=True,
        savepath=save_dir / 'lick_time_early_0p5_2s_norm'
    )

    plot_violin_with_scatter(
        norm_late_lick_time_baseline,
        norm_late_lick_time_drug,
        'grey',
        color,
        xticklabels=['Baseline', drug_name],
        ylabel='Norm. lick rate',
        title='2-4 s',
        figsize=(2.0, 2.0),
        showmainline=False,
        print_statistics=False,
        plot_statistics=True,
        save=True,
        savepath=save_dir / 'lick_time_late_2_4s_norm'
    )

    plot_violin_with_scatter(
        norm_full_lick_time_baseline,
        norm_full_lick_time_drug,
        'grey',
        color,
        xticklabels=['Baseline', drug_name],
        ylabel='Norm. lick rate',
        title='0-4 s',
        figsize=(2.0, 2.0),
        showmainline=False,
        print_statistics=False,
        plot_statistics=True,
        save=True,
        savepath=save_dir / 'lick_time_0_4s_norm'
    )

    plot_violin_with_scatter(
        early_lick_time_baseline,
        early_lick_time_drug,
        'grey',
        color,
        xticklabels=['Baseline', drug_name],
        ylabel='Lick rate (Hz)',
        title='Early licks, 0.5-2 s',
        figsize=(1.8, 1.8),
        print_statistics=False,
        plot_statistics=True,
        save=True,
        savepath=save_dir / 'lick_time_early_0p5_2s'
    )

    plot_violin_with_scatter(
        late_lick_time_baseline,
        late_lick_time_drug,
        'grey',
        color,
        xticklabels=['Baseline', drug_name],
        ylabel='Lick rate (Hz)',
        title='Late licks, 2-4 s',
        figsize=(1.8, 1.8),
        print_statistics=False,
        plot_statistics=True,
        save=True,
        savepath=save_dir / 'lick_time_late_2_4s'
    )

    plot_violin_with_scatter(
        full_lick_time_baseline,
        full_lick_time_drug,
        'grey',
        color,
        xticklabels=['Baseline', drug_name],
        ylabel='Lick rate (Hz)',
        title='Licks, 0-4 s',
        figsize=(1.8, 1.8),
        print_statistics=False,
        plot_statistics=True,
        save=True,
        savepath=save_dir / 'lick_time_0_4s'
    )

    if drug_name == 'SCH23390':
        animals = list(dict.fromkeys(animal_ids.tolist()))
        animal_early_baseline, animal_early_drug, _ = paired_animal_mean_values(
            early_lick_time_baseline,
            early_lick_time_drug,
            animal_ids
        )
        animal_late_baseline, animal_late_drug, _ = paired_animal_mean_values(
            late_lick_time_baseline,
            late_lick_time_drug,
            animal_ids
        )
        animal_reward_baseline, animal_reward_drug, _ = paired_animal_mean_values(
            reward_baseline,
            reward_drug,
            animal_ids
        )
        print('\n[SCH23390] animal-level time lick trace')
        for idx, animal in enumerate(animals):
            print(
                f'  {animal}: '
                f'n_dates={int(np.sum(animal_ids == animal))}, '
                f'early 0.5-2s '
                f'{animal_early_baseline[idx]:.4g}->'
                f'{animal_early_drug[idx]:.4g} Hz, '
                f'late 2-4s '
                f'{animal_late_baseline[idx]:.4g}->'
                f'{animal_late_drug[idx]:.4g} Hz, '
                f'reward '
                f'{animal_reward_baseline[idx]:.3g}->'
                f'{animal_reward_drug[idx]:.3g}'
            )

    # reward
    paired_stats(
        drug_name,
        'reward percentage',
        reward_baseline,
        reward_drug,
        drug_name
    )
    plot_violin_with_scatter(
        reward_baseline,
        reward_drug,
        'grey',
        color,
        xticklabels=['Baseline', drug_name],
        ylabel='Reward percentage',
        figsize=(1.8,1.8),
        print_statistics=False,
        save=True,
        savepath=save_dir / 'reward_percentage'
    )

    return reward_drug


#%% run all drugs
reward_prazosin = summarise_drug(
    'Prazosin',
    rec_list.pathAlphaBlocker,
    rec_list.sessAlphaBlocker,
    color='darkgreen',
    save_folder='prazosin',
    sch_mode=False
)

reward_propranolol = summarise_drug(
    'Propranolol',
    rec_list.pathBetaBlocker,
    rec_list.sessBetaBlocker,
    color='darkcyan',
    save_folder='propranolol',
    sch_mode=False
)

reward_SCH = summarise_drug(
    'SCH23390',
    rec_list.pathSCH,
    rec_list.sessSCH,
    color='#004D80',
    save_folder='sch23390',
    sch_mode=True
)


#%% compare rewards
fig, ax = plt.subplots(figsize=(2.2, 2.4))

drug_names = ['Prazosin', 'Propranolol', 'SCH23390']
colors = ['darkgreen', 'darkcyan', '#004D80']
data = [reward_prazosin, reward_propranolol, reward_SCH]

means = [np.mean(d) for d in data]
sems  = [sem(d) for d in data]

x = np.arange(3)

ax.bar(x, means, yerr=sems, capsize=2,
       color=colors, alpha=0.7, edgecolor='none')

# jittered scatter
for i, d in enumerate(data):
    jitter = np.random.uniform(-0.12, 0.12, size=len(d))
    ax.scatter(np.full(len(d), x[i]) + jitter,
               d,
               s=15,
               color='k',
               edgecolor='none',
               alpha=0.75)

ax.set(
    xticks=x,
    xticklabels=drug_names,
    ylabel='Reward percentage',
    ylim=(0, 1.035)
)

for s in ['top', 'right']:
    ax.spines[s].set_visible(False)

for ext in ['.png', '.pdf']:
    fig.savefig(pharmacology_figures_stem / f'reward_percentage_across_drugs{ext}',
                dpi=300, bbox_inches='tight')
