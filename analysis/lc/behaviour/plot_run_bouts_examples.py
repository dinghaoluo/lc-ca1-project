# -*- coding: utf-8 -*-
'''
Created on Fri 27 June 17:16:12 2025

example plot: single session, targeted trial window

@author: Dinghao Luo
'''

#%% imports
from pathlib import Path
import sys

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy.io as sio
import mat73

repo_root = Path(__file__).resolve().parents[3]
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))

from common_functions import mpl_formatting, smooth_convolve, gaussian_kernel_unity
import project_paths as pp
mpl_formatting()


#%% paths parameters
mice_exp_stem = pp.MICEEXP_ROOT
run_bout_stem = pp.RUN_BOUTS_STEM
save_stem     = pp.LC_EPHYS_FIGURES_STEM / 'run_onset_v_run_bout' / 'single_cell_examples_rolling'
save_stem.mkdir(parents=True, exist_ok=True)

SAMP_FREQ = 1250

KERN_SPEED = gaussian_kernel_unity(sigma=SAMP_FREQ*0.03)


#%% recname
# note: t=30 will show trials 29-31

# recname = 'A067r-20230821-01'
# t = 30
# t = 157

# recname = 'A045r-20221207-02'
# t = 181

recname = 'A032r-20220802-02'
t = 75


#%% paths
base_path     = mice_exp_stem / f'ANMD{recname[1:5]}' / recname[:14] / recname
run_bout_path = run_bout_stem / f'{recname}_run_bouts_py.csv'

alignedRun_path = (mice_exp_stem / f'ANMD{recname[1:5]}' / recname[:14] / recname /
                   f'{recname}_DataStructure_mazeSection1_TrialType1_alignRun_msess1.mat')
alignedCue_path = (mice_exp_stem / f'ANMD{recname[1:5]}' / recname[:14] / recname /
                   f'{recname}_DataStructure_mazeSection1_TrialType1_alignCue_msess1.mat')
alignedRew_path = (mice_exp_stem / f'ANMD{recname[1:5]}' / recname[:14] / recname /
                   f'{recname}_DataStructure_mazeSection1_TrialType1_alignRew_msess1.mat')

behave_lfp_path = base_path / f'{recname}_BehavElectrDataLFP.mat'

clu_path = base_path / f'{recname}.clu.1'
res_path = base_path / f'{recname}.res.1'

cell_profiles = pd.read_pickle(pp.LC_EPHYS_STEM / 'LC_all_cell_profiles.pkl')


#%% load data
run_bout_table = pd.read_csv(run_bout_path)
beh_lfp        = mat73.loadmat(str(behave_lfp_path))
aligned        = sio.loadmat(str(alignedRun_path))['trialsRun'][0][0]
alignedCue     = sio.loadmat(str(alignedCue_path))['trialsCue'][0][0]  # for cue marking
alignedRew     = sio.loadmat(str(alignedRew_path))['trialsRew'][0][0]  # for rew marking

# read cues
cueLfpInd = alignedCue['startLfpInd'].flatten()

# read rewards
rewLfpInd = alignedRew['startLfpInd'].flatten()

# read other data
tracks = beh_lfp['Track']
laps = beh_lfp['Laps']

clusters = np.loadtxt(clu_path, dtype=int, skiprows=1)
spike_times = np.loadtxt(res_path, dtype=int) / (20_000 / 1_250)
spike_times = spike_times.astype(int)

unique_clus = [clu for clu in np.unique(clusters) if clu not in [0,1]]
clu_to_row = {clu: i for i, clu in enumerate(unique_clus)}

max_time = spike_times.max() + 1
spike_map = np.zeros((len(unique_clus), max_time), dtype=int)
for time, clu in zip(spike_times, clusters):
    if clu in [0,1]:
        continue
    spike_map[clu_to_row[clu], time] = 1

spike_array = np.array([smooth_convolve(spike_map[i], sigma=int(SAMP_FREQ*.05))
                        for i in range(len(unique_clus))])

lickLfp = laps['lickLfpInd']
lickLfp_flat = []
for trial in range(len(lickLfp)):
    if isinstance(lickLfp[trial][0], np.ndarray):
        for i in range(len(lickLfp[trial][0])):
            lickLfp_flat.append(int(lickLfp[trial][0][i]))
lickLfp_flat = np.array(lickLfp_flat)

speed_MMsec = tracks['speed_MMsecAll']
speed_MMsec[speed_MMsec < 0] = np.nan
speed_MMsec = pd.Series(speed_MMsec).interpolate().bfill().ffill().values
speed_MMsec = np.convolve(speed_MMsec, KERN_SPEED, mode='same')/10

startLfpInd = aligned['startLfpInd'][0]
endLfpInd   = aligned['endLfpInd'][0]


#%% select cells
selected_indices = []
for i, clu in enumerate(unique_clus):
    cluname = f'{recname} clu{clu}'
    profile = cell_profiles.loc[cluname]
    if profile['identity'] in ['tagged', 'putative'] and profile['run_onset_peak'] is True:
        selected_indices.append(i)

if not selected_indices:
    print(f'no tagged/putative run-onset peak cells found for {recname}')

print(f'{len(selected_indices)} tagged/putative run-onset peak cells selected.')


#%% plot trials
lfp_indices_t = np.arange(startLfpInd[t]-SAMP_FREQ, min(endLfpInd[t+2], len(speed_MMsec)))
lap_start = lfp_indices_t[0]
xaxis = np.arange(0, len(lfp_indices_t)) / SAMP_FREQ

fig, ax = plt.subplots(figsize=(len(lfp_indices_t)/4300, 1.55))
ax.set(xlabel='Time (s)', ylabel='Speed (cm/s)',
       ylim=(0, 1.5 * max(speed_MMsec[lfp_indices_t])),
       xlim=(0, len(lfp_indices_t) / SAMP_FREQ),
       title=f'{recname} trial {t}')
ax.plot(xaxis, speed_MMsec[lfp_indices_t], color='royalblue')

cueLfpInd_t = cueLfpInd[np.isin(cueLfpInd, lfp_indices_t)]
ax.vlines((cueLfpInd_t - lap_start)/SAMP_FREQ, 0, ax.get_ylim()[1], 'darkgrey', zorder=10)

rewLfpInd_t = rewLfpInd[np.isin(rewLfpInd, lfp_indices_t)]
ax.vlines((rewLfpInd_t - lap_start)/SAMP_FREQ, ax.get_ylim()[1], ax.get_ylim()[1]*.95, 'forestgreen', linewidth=1.5, zorder=10)

startLfpInd_t = startLfpInd[np.isin(startLfpInd, lfp_indices_t)]
ax.vlines((startLfpInd_t - lap_start)/SAMP_FREQ, 0, ax.get_ylim()[1], 'red', linestyle='dashed', zorder=10)

run_bout_t = run_bout_table.iloc[:,1][np.isin(run_bout_table.iloc[:,1], lfp_indices_t)]
ax.vlines((run_bout_t - lap_start)/SAMP_FREQ, 0, ax.get_ylim()[1], 'green', linestyle='dashed')

ax_spk = ax.twinx()
ax_spk.set_ylabel('FR (Hz)')

spike_subset = spike_array[selected_indices, :][:, lfp_indices_t] * SAMP_FREQ
mean_trace = np.mean(spike_subset, axis=0)
mean_trace = np.clip(mean_trace, 0, np.percentile(mean_trace, 99.5))
ax_spk.plot(xaxis, mean_trace, color='black', linewidth=1.2)
ax_spk.set_ylim(0, np.max(mean_trace) * 1.1)

licks_t = lickLfp_flat[np.isin(lickLfp_flat, lfp_indices_t)]
ax.vlines((licks_t - lap_start)/SAMP_FREQ, ax.get_ylim()[1]*.96, ax.get_ylim()[1] * 0.99, 'magenta')

ax.spines['top'].set_visible(False)
ax_spk.spines['top'].set_visible(False)

plt.tight_layout()
plt.show()

for ext in ['.pdf', '.png']:
    fig.savefig(save_stem / f'{recname}_t{t}{ext}',
                dpi=300,
                bbox_inches='tight')
