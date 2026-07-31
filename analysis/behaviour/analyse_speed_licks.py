# -*- coding: utf-8 -*-
'''
Created on Mon Apr 28 17:03:16 2025

analyse single-session speed and lick profiles

@author: Dinghao Luo
'''

#%% imports
from pathlib import Path
import sys

import numpy as np
import pickle
from scipy.stats import sem, ttest_1samp, wilcoxon

repo_root = Path(__file__).resolve().parents[2]
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))

from common_functions import mpl_formatting, replace_outlier, smooth_convolve
from console_formatting import print_files_saved, print_session, print_statistics_section, print_status
from lick_time_utils import lick_time_map
import project_paths as pp
mpl_formatting()

import rec_list


#%% paths and parameters
output_path = pp.BEHAVIOUR_SESSION_PROFILES_STEM / 'speed_lick_profiles_payload.npy'

N_SHUF  = 500
RAMP_T0 = 2.0
RAMP_T1 = 4.0
i0 = int(RAMP_T0 * 1000)
i1 = int(RAMP_T1 * 1000)


#%% bad behaviour exclusion
remove_recnames = [
    'A060r-20230602-01',
    'A062r-20230626-01',

    'A063r-20230708-01',
    'A063r-20230708-02',

    'A069r-20230905-01',
    'A069r-20230905-02',

    'A070r-20231109-01',
    'A070r-20231110-01',
    'A070r-20231115-01',
    'A070r-20231116-01',
    'A070r-20231117-01',

    'A078r-20240124-01',
    'A078r-20240125-01',
    'A078r-20240129-01',
    'A078r-20240130-01',
    'A078r-20240131-01',
    'A078r-20240201-01',
    'A078r-20240202-01',

    'A093i-20240620-01',
    'A093i-20240621-01',
    'A093i-20240625-01',
    'A093i-20240626-01',
    'A093i-20240708-01',
    'A093i-20240708-02',

    'A094i-20240701-01',
    'A094i-20240705-01',
    'A094i-20240705-02',
    'A094i-20240709-01',
    'A094i-20240711-01',
    'A094i-20240712-01',
    'A094i-20240717-01',
    'A094i-20240718-01',
    'A094i-20240718-02',
    'A094i-20240719-01',
    'A094i-20240719-02',
    'A094i-20240807-01',

    'A097i-20240826-01',
    'A097i-20240826-02',
    'A097i-20240827-01',
    'A097i-20240827-02',
    'A097i-20240829-01',
    'A097i-20240829-02',

    'A098i-20240923-01',
    'A098i-20240924-01',
    'A098i-20240926-01',
    'A098i-20240927-01',
    'A098i-20240927-02',
    'A098i-20241002-01',
    'A098i-20241004-01',
    'A098i-20241007-01',
    'A098i-20241023-01',

    'A101i-20241029-02',
    'A101i-20241031-01',
    'A101i-20241031-02',
    'A101i-20241101-01',
    'A101i-20241105-02',
    'A101i-20241106-01',
    'A101i-20241107-01',
    'A101i-20241107-02',
    'A101i-20241107-03',

    'A106i-20250128-01',
    'A106i-20250128-02'
    ]

def pad_trace(trace, length=5000):
    trace = np.asarray(trace, dtype=float)
    if len(trace) >= length:
        return trace[:length]

    padded = np.full(length, np.nan)
    padded[:len(trace)] = trace
    return padded


#%% analysis
all_speed_time_curves = []
all_lick_time_curves = []
session_payloads = {}
# predictive licking containers
all_real_slopes = []
all_shuf_means  = []
all_shuf_stds   = []

experiment_paths = {
    'HPCLC': rec_list.pathHPCLCopt,
    'HPCLCterm': rec_list.pathHPCLCtermopt,
    'LC': rec_list.pathLC,
    'HPCGRABNE': rec_list.pathHPCGRABNE,
    'LCHPCGCaMP': rec_list.pathLCHPCGCaMP,
    'HPCdLightLCOpto': rec_list.pathdLightLCOpto,
}
for exp_name in ['HPCLC',
                 'HPCLCterm',
                 'LC',
                 'HPCGRABNE',
                 'LCHPCGCaMP',
                 'HPCdLightLCOpto']:
    beh_stem = pp.behaviour_experiment_stem(exp_name)
    paths = experiment_paths[exp_name]

    for path in paths:
        recname = Path(path).name

        if recname in remove_recnames:
            print_session(recname)
            print_status('skipped', 'bad behaviour list')
            continue

        print_session(recname)

        beh_path = beh_stem / f'{recname}.pkl'
        with open(beh_path, 'rb') as f:
            beh = pickle.load(f)
        if not beh['run_onsets']:
            print_status('skipped', 'immobile session')
            continue

        # speed (temporal)
        speed_times_aligned = beh['speed_times_aligned']
        speed_aligned = [replace_outlier(np.array([s[1] for s in trial]))
                         for i, trial in enumerate(speed_times_aligned)
                         if trial]
        max_length_speeds = np.max([len(trial) for trial in speed_times_aligned])
        speed_arr = np.zeros((len(speed_aligned), max_length_speeds))
        for i in range(len(speed_aligned)):
            speed_arr[i, :len(speed_aligned[i])] = speed_aligned[i]
        mean_speed_times = np.nanmean(speed_arr, axis=0)
        sem_speed_times = sem(speed_arr, axis=0)

        # speed (spatial)
        speed_distances = np.array(
            [replace_outlier(np.array(trial))
            for i, trial in enumerate(beh['speed_distances_aligned'])
            if len(trial)>0]
            )

        # licks (spatial)
        lick_maps = np.array(
            [smooth_convolve(np.array(trial), sigma=10) * 10  # convert from mm to cm
            for i, trial in enumerate(beh['lick_maps'])
            if len(trial)>0]
            )

        mean_speeds = pad_trace(np.nanmean(speed_arr, axis=0))  # 5 s
        sem_speeds = pad_trace(sem(speed_arr, axis=0, nan_policy='omit'))
        speed_time_axis = np.arange(5 * 1000) / 1000  # 50 Hz

        mean_speeds_distances = np.nanmean(speed_distances, axis=0)
        sem_speeds_distances = sem(speed_distances, axis=0, nan_policy='omit')
        speed_distance_axis = np.arange(2200) / 10

        mean_lick_maps = np.nanmean(lick_maps, axis=0)
        sem_lick_maps = sem(lick_maps, axis=0, nan_policy='omit')
        lick_distance_axis = np.arange(2200) / 10

        # licks (temporal)
        lick_times_map = lick_time_map(beh['lick_times_aligned'])

        mean_lick_times_maps = np.nanmean(lick_times_map, axis=0)
        sem_lick_times_maps = sem(lick_times_map, axis=0, nan_policy='omit')
        lick_times_axis = np.arange(5000) / 1000

        # -------------------------------
        # predictive ramp quantification
        # -------------------------------
        x = lick_times_axis[i0:i1]
        y = mean_lick_times_maps[i0:i1]

        if np.all(np.isnan(y)):
            real_slope = np.nan
        else:
            real_slope = np.polyfit(x, y, 1)[0]

        shuf_slopes = []
        for _ in range(N_SHUF):

            # circularly shift EACH TRIAL'S temporal lick trace
            shuf = np.array([
                np.roll(tr, np.random.randint(tr.size))
                for tr in lick_times_map
            ])

            # average across trials
            shuf_mean = np.nanmean(shuf, axis=0)
            y_shuf = shuf_mean[i0:i1]

            if np.all(np.isnan(y_shuf)):
                slope_shuf = np.nan
            else:
                slope_shuf = np.polyfit(x, y_shuf, 1)[0]

            shuf_slopes.append(slope_shuf)

        shuf_slopes = np.array(shuf_slopes)
        shuf_slopes = shuf_slopes[~np.isnan(shuf_slopes)]

        if len(shuf_slopes) > 2:
            shuf_mean = np.mean(shuf_slopes)
            shuf_std  = np.std(shuf_slopes)
        else:
            shuf_mean = np.nan
            shuf_std  = np.nan

        # store
        all_real_slopes.append(real_slope)
        all_shuf_means.append(shuf_mean)
        all_shuf_stds.append(shuf_std)
        # ------------------------------------
        # predictive ramp quantification ends
        # ------------------------------------

        # store for global average (clip to 5000 samples = 5 s)
        all_speed_time_curves.append(pad_trace(mean_speed_times))
        all_lick_time_curves.append(mean_lick_times_maps[:5000])

        session_payloads[recname] = {
            'exp_name': exp_name,
            'speed_time_axis': speed_time_axis,
            'speed_distance_axis': speed_distance_axis,
            'lick_distance_axis': lick_distance_axis,
            'lick_times_axis': lick_times_axis,
            'speed_arr': speed_arr,
            'speed_distances': speed_distances,
            'mean_speeds': mean_speeds,
            'sem_speeds': sem_speeds,
            'mean_speed_times': mean_speed_times,
            'sem_speed_times': sem_speed_times,
            'mean_speeds_distances': mean_speeds_distances,
            'sem_speeds_distances': sem_speeds_distances,
            'mean_lick_maps': mean_lick_maps,
            'sem_lick_maps': sem_lick_maps,
            'mean_lick_times_maps': mean_lick_times_maps,
            'sem_lick_times_maps': sem_lick_times_maps
            }


#%% all
all_speed_time_curves = np.array(all_speed_time_curves)
all_lick_time_curves = np.array(all_lick_time_curves)

mean_speed_group = np.nanmean(all_speed_time_curves, axis=0)
sem_speed_group = sem(all_speed_time_curves, axis=0, nan_policy='omit')
mean_lick_group = np.nanmean(all_lick_time_curves, axis=0)
sem_lick_group = sem(all_lick_time_curves, axis=0, nan_policy='omit')

#%% predictive temporal lick ramp group statistics
tval, p_t = ttest_1samp(all_real_slopes, 0)
wstat, p_w = wilcoxon(all_real_slopes)

# shuffle + CI
shuf_mean = np.nanmean(all_shuf_means)
shuf_std  = np.nanmean(all_shuf_stds)

lower_95 = shuf_mean - 1.96 * shuf_std
upper_95 = shuf_mean + 1.96 * shuf_std

mean_r, sem_r = np.nanmean(all_real_slopes), sem(all_real_slopes)
ymax = np.max(all_real_slopes)

print_statistics_section()
print(f'n_sessions = {len(all_real_slopes)}')
print(f'mean +/- SEM = {mean_r:.4f} +/- {sem_r:.4f}')
print(f'ttest: t={tval:.4g}, p={p_t:.2e}')
print(f'wilcoxon: W={wstat:.4g}, p={p_w:.2e}')

output_path.parent.mkdir(parents=True, exist_ok=True)
np.save(
    output_path,
    {
        'sessions': session_payloads,
        'speed_time_axis': speed_time_axis,
        'lick_times_axis': lick_times_axis,
        'all_speed_time_curves': all_speed_time_curves,
        'all_lick_time_curves': all_lick_time_curves,
        'mean_speed_group': mean_speed_group,
        'sem_speed_group': sem_speed_group,
        'mean_lick_group': mean_lick_group,
        'sem_lick_group': sem_lick_group,
        'all_real_slopes': all_real_slopes,
        'all_shuf_means': all_shuf_means,
        'all_shuf_stds': all_shuf_stds,
        'tval': tval,
        'p_t': p_t,
        'wstat': wstat,
        'p_w': p_w,
        'shuf_mean': shuf_mean,
        'shuf_std': shuf_std,
        'lower_95': lower_95,
        'upper_95': upper_95,
        'mean_r': mean_r,
        'sem_r': sem_r,
        'ymax': ymax
        }
    )
print_files_saved([
    ('speed and lick data', output_path),
])
