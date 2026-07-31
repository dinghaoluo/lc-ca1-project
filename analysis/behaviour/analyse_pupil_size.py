# -*- coding: utf-8 -*-
'''
Created on 5 Dec 2025

Analyse pupil size changes aligned to run onset.

@author: Dinghao Luo
'''

#%% imports
from pathlib import Path
import sys

import numpy as np
from scipy.stats import ttest_1samp, wilcoxon, sem

repo_root = Path(__file__).resolve().parents[2]
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))

from common_functions import mpl_formatting, normalise, smooth_convolve
from console_formatting import print_files_saved, print_session, print_statistics_section
import project_paths as pp
mpl_formatting()

import behaviour_functions as bf


#%% paths and parameters
pupil_stem  = pp.PUPIL_TRACKING_INPUT_STEM
output_path = pp.PUPIL_TRACKING_OUTPUT_STEM / 'run_aligned_pupil_payload.npy'

SAMP_FREQ     = 30  # fps
SAMP_FREQ_BEH = 1000  # for Arduino

BEF = 1  # s
AFT = 14

N_SHUF = 500


#%% define sessions to loop over
recnames = [
    # 'A057-20230510-03',
    'A057-20230511-03',
    'A057-20230516-03',
    'A057-20230517-03',
    'A057-20230517-04',
    # 'A057-20230518-03',
    # 'A057-20230518-04',
    'A057-20230519-03',
    'A057-20230522-03',
    # 'A057-20230522-04',

    'A059-20230424-01',
    'A059-20230424-02',
    'A059-20230425-02',
    'A059-20230503-02',
    'A059-20230504-03',
    # 'A059-20230505-04',
    # 'A059-20230509-02',
    'A059-20230510-03',
    # 'A059-20230512-02',
    'A059-20230523-02',
    'A059-20230523-03',
    'A059-20230524-02',
    'A059-20230524-03',
    'A059-20230525-02',

    # 'A061-20230620-02',
    # 'A061-20230620-03',
    # 'A061-20230621-02',
    # 'A061-20230621-03',
    # 'A061-20230622-02'
]

#%% analysis
avg_start_traces = []
single_session_payload = {}

for recname in recnames:
    animal  = f'ANMD{recname[1:4]}'
    day     = recname[:-3]
    print_session(recname)

    face_path  = pupil_stem / animal / day / f'{recname}_proc.npy'
    ctime_path = pupil_stem / animal / day / f'{recname}_tsdict.npy'
    txt_path   = pp.MICEEXP_ROOT / animal / f'{recname}T.txt'

    face  = np.load(face_path, allow_pickle=True).item()
    ctime = np.load(ctime_path, allow_pickle=True).item()['ctime']
    ctime = [t*1000 for t in ctime]
    ctime = [t-ctime[0] for t in ctime]

    pupil_area = face['pupil'][0]['area']

    beh = bf.process_behavioural_data(txt_path)

    # get stim
    trial_statements = beh['trial_statements']
    opto_cds         = [t[15] for t in trial_statements]

    ctrl_idx = np.array([trial for trial, cond in enumerate(opto_cds)
                         if cond == '0' and trial > 1 and trial < len(opto_cds)-1])

    t_ST = np.array(beh['run_onsets'])

    logfile = open(txt_path, 'r')

    line = ['$']
    t_camsync = []

    while line[0].find('$') == 0:
        if line[0] == '$SY':
            t_camsync.append(float(line[1]))
        line = logfile.readline().rstrip('\n').split(',')
        if len(line) == 1:
            line = logfile.readline().rstrip('\n').split(',')

    t_ST = t_ST[ctrl_idx]

    tot_camsync = len(t_camsync)
    for sync in range(1, tot_camsync):
        dt_curr = t_camsync[sync] - t_camsync[sync - 1]
        if dt_curr > 0 and (dt_curr < 25 or dt_curr > 75):
            raise Exception(f'{recname}: sync #{sync} invalid dt')

    tot_frame = len(ctime)
    if tot_frame != pupil_area.shape[0]:
        raise Exception(f'{recname}: ctime length mismatch')

    ctime = np.array([t + t_camsync[0] for t in ctime])

    timebef = BEF * SAMP_FREQ_BEH
    timeaft = AFT * SAMP_FREQ_BEH

    frames_by_start, pupil_by_start = [], []
    for start in t_ST:
        if recname in ['A057-20230519-03',
                       'A057-20230517-04',
                       'A059-20230503-02',
                       'A059-20230523-03',
                       'A059-20230510-03']:
            start += 800
        if recname in ['A057-20230510-03']:
            start += 1000
        window = [start - timebef, start + timeaft]
        frame_curr = [f for f in range(tot_frame) if window[0] < ctime[f] < window[1]]
        if frame_curr:
            frames_by_start.append(frame_curr)
            pupil_by_start.append(smooth_convolve(pupil_area[frame_curr], sigma=SAMP_FREQ/20))  # smoothing

    min_len_start = min(len(p) for p in pupil_by_start)

    # organise data
    session_traces = [normalise(p[:min_len_start]) for p in pupil_by_start]

    avg_start = np.nanmean(session_traces, axis=0)
    sem_start = np.nanstd(session_traces, axis=0) / np.sqrt(len(pupil_by_start))
    x = np.arange(min_len_start) / SAMP_FREQ - BEF  # - 1 seccond

    avg_start_traces.append(avg_start)
    single_session_payload[recname] = {
        'x': x,
        'avg_start': avg_start,
        'sem_start': sem_start
        }


#%% quantify per-session pupil modulation
session_real_deltas = []
session_shuf_means  = []
session_shuf_stds   = []
session_shuf_bands  = []   # (low95, high95)

for trace in avg_start_traces:

    n = len(trace)
    t = np.arange(n) / SAMP_FREQ - BEF

    # analysis windows
    baseline_mask = (t >= -0.5) & (t <= 0)
    response_mask = (t >= 1.0) & (t <= 2.0)

    # real
    base = np.nanmean(trace[baseline_mask])
    resp = np.nanmean(trace[response_mask])
    real_delta = resp - base
    session_real_deltas.append(real_delta)

    # shuffle
    shuf_deltas = []
    for _ in range(N_SHUF):
        shift = np.random.randint(n)
        shuf_trace = np.roll(trace, shift)

        base_s = np.nanmean(shuf_trace[baseline_mask])
        resp_s = np.nanmean(shuf_trace[response_mask])

        shuf_deltas.append(resp_s - base_s)

    shuf_deltas = np.array(shuf_deltas)
    session_shuf_means.append(np.nanmean(shuf_deltas))
    session_shuf_stds.append(np.nanstd(shuf_deltas))

    low95  = np.percentile(shuf_deltas, 2.5)
    high95 = np.percentile(shuf_deltas, 97.5)
    session_shuf_bands.append((low95, high95))


#%% mean
trim_end = int(SAMP_FREQ * 1 + SAMP_FREQ * 4)

trimmed = np.array([t[:trim_end] for t in avg_start_traces if len(t) > 200])

# trimmed = np.array([
#     (trace[:trim_end] - np.nanmin(trace[:trim_end])) /
#     (np.nanmax(trace[:trim_end]) - np.nanmin(trace[:trim_end]))
#     for trace in avg_start_traces
# ])

grand_avg = np.nanmean(trimmed, axis=0)
grand_sem = np.nanstd(trimmed, axis=0) / np.sqrt(trimmed.shape[0])
x_axis = np.arange(trim_end) / 30 - 1  # same as before

#%% violin plot for deltas
# IQR
q1 = np.percentile(session_real_deltas, 25)
q3 = np.percentile(session_real_deltas, 75)
iqr = q3 - q1

tval, p_t = ttest_1samp(session_real_deltas, 0)
wstat, p_w = wilcoxon(session_real_deltas)

# shuffle + CI
shuf_mean = np.nanmean(session_shuf_means)
shuf_std  = np.nanmean(session_shuf_means)

lower_95 = shuf_mean - 1.96 * shuf_std
upper_95 = shuf_mean + 1.96 * shuf_std

mean_r, sem_r = np.nanmean(session_real_deltas), sem(session_real_deltas)
ymax = np.max(session_real_deltas)

print_statistics_section()
print(f'n_sessions = {len(session_real_deltas)}')
print(f'median = {np.nanmedian(session_real_deltas):.4f}')
print(f'IQR = [{q1:.4f}, {q3:.4f}]')
print(f'mean +/- SEM = {mean_r:.4f} +/- {sem_r:.4f}')
print(f'ttest: t={tval:.4g}, p={p_t:.2e}')
print(f'wilcoxon: W={wstat:.4g}, p={p_w:.2e}')

output_path.parent.mkdir(parents=True, exist_ok=True)
np.save(
    output_path,
    {
        'single_sessions': single_session_payload,
        'avg_start_traces': avg_start_traces,
        'x_axis': x_axis,
        'grand_avg': grand_avg,
        'grand_sem': grand_sem,
        'session_real_deltas': session_real_deltas,
        'session_shuf_means': session_shuf_means,
        'session_shuf_stds': session_shuf_stds,
        'session_shuf_bands': session_shuf_bands,
        'q1': q1,
        'q3': q3,
        'iqr': iqr,
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
    ('analysis data', output_path),
])
