# -*- coding: utf-8 -*-
'''
Created on Fri Jun 13 17:09:31 2025

plot profiles of run-onset PyrUp/PyrDown cells based on raw spike rate change

@author: Dinghao Luo

'''

#%% imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
from pathlib import Path


#%% project imports and paths
repo_root = Path(__file__).resolve().parents[2]
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))
from common_functions import normalise, mpl_formatting
import project_paths as pp
mpl_formatting()

HPC_STEM = pp.HPC_EPHYS_STEM
RUN_ONSET_STEM = pp.HPC_EPHYS_FIGURES_STEM / 'run_onset_response'


#%% load dataframe
print('loading dataframe...')
df = pd.read_pickle(HPC_STEM / 'hpc_all_profiles.pkl')
df_pyr = df[(df['cell_identity']=='pyr') & (df['rectype']=='HPCLC')]


#%% parameters
xaxis = np.arange(-1, 4, 1/1250)


#%% calculate raw change
df_pyr['raw_change'] = df_pyr['prof_mean_MATLAB'].apply(
    lambda profile: np.mean(profile[4375:5625]) - np.mean(profile[1875:3125])
    )
df_pyr['ctrl_raw_change'] = df_pyr['prof_ctrl_mean_MATLAB'].apply(
    lambda profile: np.mean(profile[4375:5625]) - np.mean(profile[1875:3125])
    )
df_pyr['stim_raw_change'] = df_pyr['prof_stim_mean_MATLAB'].apply(
    lambda profile: np.mean(profile[4375:5625]) - np.mean(profile[1875:3125])
    )


#%% sort dataframe by baseline pre-post ratios
df_sorted = df_pyr.sort_values(by='raw_change', ascending=False)
df_sorted_ctrl = df_pyr.sort_values(by='ctrl_raw_change', ascending=False)
df_sorted_stim = df_pyr.sort_values(by='stim_raw_change', ascending=False)


#%% matrices
pop_mat_ctrl = df_sorted_ctrl['prof_ctrl_mean'].to_numpy()
pop_mat_ctrl = np.asarray([normalise(cell[2500:2500+5*1250]) for cell in pop_mat_ctrl])

pop_mat_stim = df_sorted_stim['prof_stim_mean'].to_numpy()
pop_mat_stim = np.asarray([normalise(cell[2500:2500+5*1250]) for cell in pop_mat_stim])


#%% overall plot
fig, ax = plt.subplots(figsize=(2.4,1.9))

cim = ax.imshow(pop_mat_ctrl, aspect='auto', cmap='Greys', interpolation='sinc',
                extent=(-1, 4, 0, pop_mat_ctrl.shape[0]))

ax.set(title='rate based -- ctrl.',
       xlabel='time from run-onset (s)',
       ylabel='cell #')

plt.colorbar(cim, shrink=.5, ticks=[0, 1])

RUN_ONSET_STEM.mkdir(parents=True, exist_ok=True)
for ext in ['.png', '.pdf']:
    fig.savefig(
        RUN_ONSET_STEM / f'HPCLC_run_onset_ctrl_ratebased_Greys{ext}',
        dpi=300,
        bbox_inches='tight'
        )

fig, ax = plt.subplots(figsize=(2.4,1.9))

sim = ax.imshow(pop_mat_stim, aspect='auto', cmap='Greys', interpolation='sinc',
                extent=(-1, 4, 0, pop_mat_stim.shape[0]))
ax.set(title='rate based -- stim.',
       xlabel='time from run-onset (s)',
       ylabel='cell #')

plt.colorbar(sim, shrink=.5, ticks=[0, 1])

for ext in ['.png', '.pdf']:
    fig.savefig(
        RUN_ONSET_STEM / f'HPCLC_run_onset_stim_ratebased_Greys{ext}',
        dpi=300,
        bbox_inches='tight'
        )
