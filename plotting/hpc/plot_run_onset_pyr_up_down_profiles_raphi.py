# -*- coding: utf-8 -*-
'''
Created on Mon 10 Mar 15:04:01 2025
Modified on 21 Jan 2026

plot run-onset PyrUp and PyrDown cells for Raphi's data

@author: Dinghao Luo
'''

#%% imports
from pathlib import Path
import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import PowerNorm
from scipy.stats import sem
import pandas as pd

repo_root = Path(__file__).resolve().parents[2]
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))

from common_functions import mpl_formatting, normalise
import project_paths as pp
mpl_formatting()

import rec_list
paths = rec_list.pathHPC_Raphi


#%% paths and parameters
HPC_stem       = pp.HPC_EPHYS_STEM
run_onset_stem = pp.HPC_EPHYS_FIGURES_STEM / 'run_onset_response_raphi'
run_onset_stem.mkdir(parents=True, exist_ok=True)

RUN_ONSET_BIN = 3750  # in samples
SAMP_FREQ     = 1250  # in Hz

BEF = 1  # in seconds
AFT = 4  # in seconds

XAXIS = np.arange(-SAMP_FREQ * BEF, SAMP_FREQ * AFT) / SAMP_FREQ

PROF_WINDOW = [
    int(RUN_ONSET_BIN - SAMP_FREQ * BEF),
    int(RUN_ONSET_BIN + SAMP_FREQ * AFT)
    ]
PYRUP_CLASS = 'run-onset ON'
PYRDOWN_CLASS = 'run-onset OFF'


#%% load dataframe
print('loading dataframe...')

cell_profiles_path = HPC_stem / 'hpc_all_profiles_raphi.pkl'

cell_profiles = pd.read_pickle(cell_profiles_path)

# ------------
# scale (temp)
# ------------
profile_cols = [
    'prof_mean', 'prof_sem',
    'prof_stim_mean', 'prof_stim_sem',
    'prof_ctrl_mean', 'prof_ctrl_sem'
]

recnums = cell_profiles['recname'].str[1:4].astype(int)
scale_mask = recnums > 40

for col in profile_cols:
    cell_profiles.loc[scale_mask, col] = (
        cell_profiles.loc[scale_mask, col].apply(
            lambda x: x * SAMP_FREQ if isinstance(x, np.ndarray) else x
        )
    )

print(f'scaled {scale_mask.sum()} cells (recname > A040)')
# ------------
# scale (temp) ends
# ------------

df_pyr = cell_profiles[cell_profiles['cell_identity']=='pyr']  # pyramidal only

df_pyrup  = df_pyr[df_pyr['class']==PYRUP_CLASS]
df_pyrdown = df_pyr[df_pyr['class']==PYRDOWN_CLASS]


#%% statistics first
# unique sessions
sessions = df_pyr['recname'].unique()
n_sessions = len(sessions)

# unique animals (Axxx from recname[:4])
animals = df_pyr['recname'].str[:4].unique()
n_animals = len(animals)

# PyrUp/PyrDown
n_pyrup  = len(df_pyrup)
n_pyrdown = len(df_pyrdown)
n_tot = len(df_pyr)


#%% sorting
df_sorted = df_pyr.sort_values(by='pre_post')

pop_mat = df_sorted['prof_mean'].to_numpy()

pop_mat = np.asarray([normalise(cell[PROF_WINDOW[0] : PROF_WINDOW[1]])
    for cell in pop_mat
    if not np.isnan(normalise(cell[PROF_WINDOW[0] : PROF_WINDOW[1]])[0])])


#%% overall plot
fig, ax = plt.subplots(figsize=(2.4,1.9))

gim = ax.imshow(pop_mat, aspect='auto', cmap='Greys', norm=PowerNorm(gamma=0.8),
                extent=(-1, 4, 0, pop_mat.shape[0]))
ax.set(title=f'{n_pyrup}/{n_tot} PyrUp ({round(n_pyrup/n_tot, 4)})\n{n_pyrdown}/{n_tot} PyrDown ({round(n_pyrdown/n_tot, 4)})\nn_sess={n_sessions}, n_anm={n_animals}',
       xlabel='Time from run onset (s)',
       ylabel='Cell #', yticks=[4000, 8000, 12000])

plt.colorbar(gim, shrink=.5, ticks=[0, 1])

for ext in ['.png', '.pdf']:
    fig.savefig(
        run_onset_stem / f'all_run_onset_Greys{ext}',
        dpi=300,
        bbox_inches='tight'
        )


#%% mean profiles
pyrup_all = [cell.prof_mean[PROF_WINDOW[0] : PROF_WINDOW[1]] for cell in
          df_pyrup.itertuples(index=False)]
pyrup_all_mean = np.mean(pyrup_all, axis=0)
pyrup_all_sem  = sem(pyrup_all, axis=0)

pyrdown_all = [cell.prof_mean[PROF_WINDOW[0] : PROF_WINDOW[1]] for cell in
           df_pyrdown.itertuples(index=False)]
pyrdown_all_mean = np.mean(pyrdown_all, axis=0)
pyrdown_all_sem  = sem(pyrdown_all, axis=0)

# plotting
fig, ax = plt.subplots(figsize=(2.6,2))

pyrup_ln, = ax.plot(XAXIS, pyrup_all_mean, lw=1, c='firebrick')
ax.fill_between(XAXIS,
                pyrup_all_mean + pyrup_all_sem,
                pyrup_all_mean - pyrup_all_sem,
                color='firebrick', edgecolor='none', alpha=.3)
for p in ['top', 'right']:
    ax.spines[p].set_visible(False)

ax.set(xlabel='Time from run onset (s)', xticks=[0,2,4], xlim=(-1,4),
       ylabel='Firing rate (Hz)', yticks=[1, 1.5, 2, 2.5], ylim=(0.8, 2.6))

for ext in ['.png', '.pdf']:
    fig.savefig(run_onset_stem / f'all_PyrUp_curve_raphi{ext}',
        dpi=300,
        bbox_inches='tight'
        )

fig, ax = plt.subplots(figsize=(2.6,2))

pyrdown_ln, = ax.plot(XAXIS, pyrdown_all_mean, lw=1, c='purple')
ax.fill_between(XAXIS,
                pyrdown_all_mean + pyrdown_all_sem,
                pyrdown_all_mean - pyrdown_all_sem,
                color='purple', edgecolor='none', alpha=.3)
for p in ['top', 'right']:
    ax.spines[p].set_visible(False)

ax.set(xlabel='Time from run onset (s)', xticks=[0,2,4], xlim=(-1,4),
       ylabel='Firing rate (Hz)', yticks=[1, 1.5, 2], ylim=(0.8, 2.6))

for ext in ['.png', '.pdf']:
    fig.savefig(run_onset_stem / f'all_PyrDown_curve_raphi{ext}',
        dpi=300,
        bbox_inches='tight'
        )
