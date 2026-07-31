# -*- coding: utf-8 -*-
'''
Created on Mon 10 Mar 15:04:01 2025

plot run-onset PyrUp and PyrDown cells

@author: Dinghao Luo
'''

#%% imports
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import sem
import pandas as pd
import sys
from pathlib import Path

repo_root = Path(__file__).resolve().parents[2]
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))
from common_functions import mpl_formatting, normalise_to_all
from plotting_functions import plot_violin_with_scatter
import project_paths as pp
mpl_formatting()

import rec_list
pathHPCLC = rec_list.pathHPCLCopt
pathHPCLCterm = rec_list.pathHPCLCtermopt


#%% parameters
run_onset_bin = 3750  # in samples
samp_freq = 1250  # in Hz
time_bef = 1  # in seconds
time_aft = 4  # in seconds
xaxis = np.arange(-samp_freq*time_bef, samp_freq*time_aft) / samp_freq
prof_window = (run_onset_bin-samp_freq*time_bef, run_onset_bin+samp_freq*time_aft)
PYRUP_CLASS = 'run-onset ON'
PYRDOWN_CLASS = 'run-onset OFF'
UNRESPONSIVE_CLASS = 'run-onset unresponsive'
HPC_STEM = pp.HPC_EPHYS_STEM
RUN_ONSET_STEM = pp.HPC_EPHYS_FIGURES_STEM / 'run_onset_response'

def save_figure(fig, filepath, **kwargs):
    filepath.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(filepath, **kwargs)

def plot_transition_matrix(matrix, title):
    fig, ax = plt.subplots(figsize=(4,4))

    values = matrix.to_numpy(dtype=float)
    image = ax.imshow(values, cmap='viridis', aspect='auto')

    for row in range(values.shape[0]):
        for col in range(values.shape[1]):
            ax.text(col, row, f'{values[row, col]:.2f}',
                    ha='center', va='center', color='w', fontsize=9)

    ax.set_xticks(np.arange(len(matrix.columns)))
    ax.set_yticks(np.arange(len(matrix.index)))
    ax.set_xticklabels(matrix.columns)
    ax.set_yticklabels(matrix.index)
    ax.set_xlabel('stim. class', fontsize=10)
    ax.set_ylabel('ctrl. class', fontsize=10)
    ax.set_title(title, fontsize=10, fontweight='bold')
    ax.tick_params(axis='x', labelsize=9, rotation=45)
    ax.tick_params(axis='y', labelsize=9)
    image.set_clim(0, 1)

    return fig, ax


#%% load dataframe
print('loading dataframe...')
cell_profiles = pd.read_pickle(HPC_STEM / 'hpc_all_profiles.pkl')
df_pyr = cell_profiles[cell_profiles['cell_identity']=='pyr']  # pyramidal only

df_pyrup = df_pyr[df_pyr['class']==PYRUP_CLASS]
df_pyrdown = df_pyr[df_pyr['class']==PYRDOWN_CLASS]


#%% plain and simple first--just the mean profiles
pyrup_all = [cell.prof_ctrl_mean[2500:3750+4*1250] for cell in
          df_pyr[df_pyr['class_ctrl']==PYRUP_CLASS]
          .itertuples(index=False)]
pyrup_all_mean = np.mean(pyrup_all, axis=0)
pyrup_all_sem = sem(pyrup_all, axis=0)

pyrdown_all = [cell.prof_ctrl_mean[2500:3750+4*1250] for cell in
           df_pyr[df_pyr['class_ctrl']==PYRDOWN_CLASS]
           .itertuples(index=False)]
pyrdown_all_mean = np.mean(pyrdown_all, axis=0)
pyrdown_all_sem = sem(pyrdown_all, axis=0)

fig, ax = plt.subplots(figsize=(2.3,2))

pyrup_ln, = ax.plot(xaxis, pyrup_all_mean, lw=1, c='firebrick')
ax.fill_between(xaxis,
                pyrup_all_mean + pyrup_all_sem,
                pyrup_all_mean - pyrup_all_sem,
                color='firebrick', edgecolor='none', alpha=.25)
for p in ['top', 'right']:
    ax.spines[p].set_visible(False)

ax.set(xlabel='time from run-onset (s)', xticks=[0,2,4], xlim=(-1,4),
       ylabel='spike rate (Hz)', yticks=[1,2,3,4], ylim=(1.2, 3.5))

for ext in ['.png', '.pdf']:
    save_figure(
        fig,
        RUN_ONSET_STEM / f'all_PyrUp_curve{ext}',
        dpi=300,
        bbox_inches='tight'
        )

fig, ax = plt.subplots(figsize=(2.3,2))

pyrdown_ln, = ax.plot(xaxis, pyrdown_all_mean, lw=1, c='purple')
ax.fill_between(xaxis,
                pyrdown_all_mean + pyrdown_all_sem,
                pyrdown_all_mean - pyrdown_all_sem,
                color='purple', edgecolor='none', alpha=.25)
for p in ['top', 'right']:
    ax.spines[p].set_visible(False)

ax.set(xlabel='time from run-onset (s)', xticks=[0,2,4], xlim=(-1,4),
       ylabel='spike rate (Hz)', yticks=[1,2,3,4], ylim=(1.2, 3.5))

for ext in ['.png', '.pdf']:
    save_figure(
        fig,
        RUN_ONSET_STEM / f'all_PyrDown_curve{ext}',
        dpi=300,
        bbox_inches='tight'
        )


#%% good v bad trials
pyrup_all_good = [cell.prof_good_mean for cell in
               df_pyr[df_pyr['class_ctrl']==PYRUP_CLASS]
               .itertuples(index=False) if len(cell.prof_good_mean)==12500]
pyrup_all_bad = [cell.prof_bad_mean for cell in
              df_pyr[df_pyr['class_ctrl']==PYRUP_CLASS]
              .itertuples(index=False) if len(cell.prof_bad_mean)==12500]

pyrup_all_good_mean = np.mean(pyrup_all_good, axis=0)[2500:3750+4*1250]
pyrup_all_good_sem = sem(pyrup_all_good, axis=0)[2500:3750+4*1250]
pyrup_all_bad_mean = np.mean(pyrup_all_bad, axis=0)[2500:3750+4*1250]
pyrup_all_bad_sem = sem(pyrup_all_bad, axis=0)[2500:3750+4*1250]

fig, ax = plt.subplots(figsize=(2.2, 1.6))
ax.plot(xaxis, pyrup_all_good_mean, c='firebrick')
ax.fill_between(xaxis, pyrup_all_good_mean+pyrup_all_good_sem,
                       pyrup_all_good_mean-pyrup_all_good_sem,
                       color='firebrick', edgecolor='none', alpha=.2)
ax.plot(xaxis, pyrup_all_bad_mean, c='grey')
ax.fill_between(xaxis, pyrup_all_bad_mean+pyrup_all_bad_sem,
                       pyrup_all_bad_mean-pyrup_all_bad_sem,
                       color='grey', edgecolor='none', alpha=.2)


#%% plot ctrl v stim profiles
pyrup_all_ctrl = [cell.prof_mean for cell in
               df_pyr[df_pyr['class_ctrl']==PYRUP_CLASS]
               .itertuples(index=False)]
pyrup_all_stim = [cell.prof_stim_mean for cell in
               df_pyr[df_pyr['class_stim']==PYRUP_CLASS]
               .itertuples(index=False)]
pyrdown_all_ctrl = [cell.prof_mean for cell in
               df_pyr[df_pyr['class_ctrl']==PYRDOWN_CLASS]
               .itertuples(index=False)]
pyrdown_all_stim = [cell.prof_stim_mean for cell in
               df_pyr[df_pyr['class_stim']==PYRDOWN_CLASS]
               .itertuples(index=False)]

pyrup_all_ctrl_mean = np.mean(pyrup_all_ctrl, axis=0)
pyrup_all_stim_mean = np.mean(pyrup_all_stim, axis=0)
pyrdown_all_ctrl_mean = np.mean(pyrdown_all_ctrl, axis=0)
pyrdown_all_stim_mean = np.mean(pyrdown_all_stim, axis=0)

pyrup_all_ctrl_sem = sem(pyrup_all_ctrl, axis=0)
pyrup_all_stim_sem = sem(pyrup_all_stim, axis=0)
pyrdown_all_ctrl_sem = sem(pyrdown_all_ctrl, axis=0)
pyrdown_all_stim_sem = sem(pyrdown_all_stim, axis=0)

# PyrUp cells
fig, ax = plt.subplots(figsize=(2.1,1.5))
pyrup_stim_ln, = ax.plot(
    xaxis,
    pyrup_all_stim_mean[prof_window[0]:prof_window[1]],
    color='firebrick', linewidth=1, zorder=10)
ax.fill_between(
    xaxis,
    pyrup_all_stim_mean[prof_window[0]:prof_window[1]]+pyrup_all_stim_sem[prof_window[0]:prof_window[1]],
    pyrup_all_stim_mean[prof_window[0]:prof_window[1]]-pyrup_all_stim_sem[prof_window[0]:prof_window[1]],
    alpha=.25, color='firebrick', edgecolor='none', zorder=10)
pyrup_ctrl_ln, = ax.plot(
    xaxis,
    pyrup_all_ctrl_mean[prof_window[0]:prof_window[1]],
    color='grey', linewidth=1)
ax.fill_between(
    xaxis,
    pyrup_all_ctrl_mean[prof_window[0]:prof_window[1]]+pyrup_all_ctrl_sem[prof_window[0]:prof_window[1]],
    pyrup_all_ctrl_mean[prof_window[0]:prof_window[1]]-pyrup_all_ctrl_sem[prof_window[0]:prof_window[1]],
    alpha=.25, color='grey', edgecolor='none')

ax.legend(
    [pyrup_stim_ln, pyrup_ctrl_ln], ['stim.', 'ctrl.'],
    frameon=False, fontsize=6)

ax.set(title='PyrUp',
       xlabel='time from run-onset (s)', ylabel='spike rate (Hz)',
       xlim=(-time_bef, time_aft), xticks=(0,2,4))
ax.title.set_fontsize(10)
for p in ['top', 'right']:
    ax.spines[p].set_visible(False)

for ext in ('.png', '.pdf'):
    save_figure(
        fig,
        RUN_ONSET_STEM / 'ctrl_stim' / f'PyrUp_ctrl_stim{ext}',
        dpi=300,
        bbox_inches='tight'
        )

# PyrDown cells
fig, ax = plt.subplots(figsize=(2.2,1.5))
pyrdown_stim_ln, = ax.plot(
    xaxis,
    pyrdown_all_stim_mean[prof_window[0]:prof_window[1]],
    color='purple', linewidth=1, zorder=10)
ax.fill_between(
    xaxis,
    pyrdown_all_stim_mean[prof_window[0]:prof_window[1]]+pyrdown_all_stim_sem[prof_window[0]:prof_window[1]],
    pyrdown_all_stim_mean[prof_window[0]:prof_window[1]]-pyrdown_all_stim_sem[prof_window[0]:prof_window[1]],
    alpha=.25, color='purple', edgecolor='none', zorder=10)
pyrdown_ctrl_ln, = ax.plot(
    xaxis,
    pyrdown_all_ctrl_mean[prof_window[0]:prof_window[1]],
    color='grey', linewidth=1)
ax.fill_between(
    xaxis,
    pyrdown_all_ctrl_mean[prof_window[0]:prof_window[1]]+pyrdown_all_ctrl_sem[prof_window[0]:prof_window[1]],
    pyrdown_all_ctrl_mean[prof_window[0]:prof_window[1]]-pyrdown_all_ctrl_sem[prof_window[0]:prof_window[1]],
    alpha=.25, color='grey', edgecolor='none')

ax.legend(
    [pyrdown_stim_ln, pyrdown_ctrl_ln], ['stim.', 'ctrl.'],
    frameon=False, fontsize=6)

ax.set(title='PyrDown',
       xlabel='time from run-onset (s)', ylabel='spike rate (Hz)',
       xlim=(-time_bef, time_aft), xticks=(0,2,4))

for p in ['top', 'right']:
    ax.spines[p].set_visible(False)

for ext in ('.png', '.pdf'):
    save_figure(
        fig,
        RUN_ONSET_STEM / 'ctrl_stim' / f'PyrDown_ctrl_stim{ext}',
        dpi=300,
        bbox_inches='tight'
        )


#%% transition matrices
class_labels = {
    PYRDOWN_CLASS: 'PyrDown',
    PYRUP_CLASS: 'PyrUp',
    UNRESPONSIVE_CLASS: 'unresponsive'}

transition_matrix = pd.crosstab(
    df_pyr['class_ctrl'], df_pyr['class_stim'],
    normalize='index'
    ).rename(index=class_labels, columns=class_labels)

fig, ax = plot_transition_matrix(
    transition_matrix,
    'cell class transition (ctrl. to stim.)',
    )

for ext in ('.png', '.pdf'):
    save_figure(fig,
                RUN_ONSET_STEM / 'ctrl_stim' / f'cell_class_transition_matrix_ctrl_stim{ext}',
                dpi=300, bbox_inches='tight')

# LC-opt
df_pyr_LCopt = df_pyr[df_pyr['rectype']=='HPCLC']
transition_matrix_LCopt = pd.crosstab(
    df_pyr_LCopt['class_ctrl'], df_pyr_LCopt['class_stim'],
    normalize='index'
    ).rename(index=class_labels, columns=class_labels)

fig, ax = plot_transition_matrix(
    transition_matrix_LCopt,
    'cell class transition LCopt (ctrl. to stim.)',
    )

for ext in ('.png', '.pdf'):
    save_figure(fig,
                RUN_ONSET_STEM / 'ctrl_stim' / f'cell_class_transition_matrix_ctrl_stim_LCopt{ext}',
                dpi=300, bbox_inches='tight')

# LC-opt
df_pyr_LCtermopt = df_pyr[df_pyr['rectype']=='HPCLCterm']
transition_matrix_LCtermopt = pd.crosstab(
    df_pyr_LCtermopt['class_ctrl'], df_pyr_LCtermopt['class_stim'],
    normalize='index'
    ).rename(index=class_labels, columns=class_labels)

fig, ax = plot_transition_matrix(
    transition_matrix_LCtermopt,
    'cell class transition LCtermopt (ctrl. to stim.)',
    )

for ext in ('.png', '.pdf'):
    save_figure(fig,
                RUN_ONSET_STEM / 'ctrl_stim' / f'cell_class_transition_matrix_ctrl_stim_LCtermopt{ext}',
                dpi=300, bbox_inches='tight')


#%% pre-post ratio for remaining PyrUp and new PyrUp cells
# step 1: identify PyrUp-PyrUp and extract pre-post
df_pyrup_pyrup = df_pyr[(df_pyr['class']==PYRUP_CLASS) &
                  (df_pyr['class_stim']==PYRUP_CLASS)]
pre_post_pyrup_pyrup = pd.to_numeric(df_pyrup_pyrup['pre_post_stim']).to_numpy()

# step 2: identify other/PyrDown-PyrUp and extract pre-post
df_other_pyrup = df_pyr[(df_pyr['class']!=PYRUP_CLASS) &
                     (df_pyr['class_stim']==PYRUP_CLASS)]
pre_post_other_pyrup = pd.to_numeric(df_other_pyrup['pre_post_stim']).to_numpy()

(RUN_ONSET_STEM / 'ctrl_stim' / 'new_PyrUp_cells_pre_post').parent.mkdir(parents=True, exist_ok=True)
plot_violin_with_scatter(
    pre_post_pyrup_pyrup, pre_post_other_pyrup,
    'firebrick', 'darkorange',
    showscatter=True,
    paired=False,
    ylabel='pre-post ratio',
    xticklabels=('cons.\nPyrUp', 'new\nPyrUp'),
    save=True,
    savepath=str(RUN_ONSET_STEM / 'ctrl_stim' / 'new_PyrUp_cells_pre_post')
    )


#%% persistent vs newly-induced
df_pyr_sorted = df_pyr.sort_values(by='pre_post_ctrl')

pyrup_pers_ctrl = [cell.prof_ctrl_mean for cell in
                df_pyr_sorted[df_pyr_sorted['class_ctrl']==PYRUP_CLASS]
                .itertuples(index=False)]
pyrup_pers_stim = [cell.prof_stim_mean for cell in
                df_pyr_sorted[df_pyr_sorted['class_ctrl']==PYRUP_CLASS]
                .itertuples(index=False)]
pyrup_pers_concat = [np.concatenate((pyrup_pers_ctrl[i], pyrup_pers_stim[i]))
                  for i in range(len(pyrup_pers_ctrl))]
pyrup_pers_ctrl_norm = [normalise_to_all(pyrup_pers_ctrl[i], pyrup_pers_concat[i])
                     for i in range(len(pyrup_pers_ctrl))]
pyrup_pers_stim_norm = [normalise_to_all(pyrup_pers_stim[i], pyrup_pers_concat[i])
                     for i in range(len(pyrup_pers_ctrl))]

pyrup_new_ctrl = [cell.prof_ctrl_mean for cell in
               df_pyr_sorted[(df_pyr_sorted['class_ctrl']!=PYRUP_CLASS) &
                             (df_pyr_sorted['class_stim']==PYRUP_CLASS)]
               .itertuples(index=False)]
pyrup_new_stim = [cell.prof_stim_mean for cell in
               df_pyr_sorted[(df_pyr_sorted['class_ctrl']!=PYRUP_CLASS) &
                             (df_pyr_sorted['class_stim']==PYRUP_CLASS)]
               .itertuples(index=False)]
pyrup_new_concat = [np.concatenate((pyrup_new_ctrl[i], pyrup_new_stim[i]))
                 for i in range(len(pyrup_new_ctrl))]
pyrup_new_ctrl_norm = [normalise_to_all(pyrup_new_ctrl[i], pyrup_new_concat[i])
                    for i in range(len(pyrup_new_ctrl))]
pyrup_new_stim_norm = [normalise_to_all(pyrup_new_stim[i], pyrup_new_concat[i])
                    for i in range(len(pyrup_new_ctrl))]

pyrdown_all_ctrl = [cell.prof_mean for cell in
               df_pyr[df_pyr['class_ctrl']==PYRDOWN_CLASS]
               .itertuples(index=False)]
pyrdown_all_stim = [cell.prof_stim_mean for cell in
               df_pyr[df_pyr['class_stim']==PYRDOWN_CLASS]
               .itertuples(index=False)]


#%% plotting
fig, axs = plt.subplots(2,1, figsize=(2,3))

axs[0].imshow(pyrup_pers_ctrl_norm, cmap='viridis', interpolation='none',
              extent=(-1, 4, 0, len(pyrup_pers_ctrl)))
axs[0].set_aspect(.003)
axs[1].imshow(pyrup_pers_stim_norm, cmap='viridis', interpolation='none',
              extent=(-1, 4, 0, len(pyrup_pers_ctrl)))
axs[1].set_aspect(.003)

save_figure(fig, RUN_ONSET_STEM / 'persistent_matrix.png',
            dpi=300, bbox_inches='tight')

fig, axs = plt.subplots(2,1, figsize=(2,3))

axs[0].imshow(pyrup_new_ctrl_norm, cmap='viridis', interpolation='none',
              extent=(-1, 4, 0, len(pyrup_new_ctrl)))
axs[0].set_aspect(.005)
axs[1].imshow(pyrup_new_stim_norm, cmap='viridis', interpolation='none',
              extent=(-1, 4, 0, len(pyrup_new_ctrl)))
axs[1].set_aspect(.005)

save_figure(fig, RUN_ONSET_STEM / 'new_matrix.png',
            dpi=300, bbox_inches='tight')
