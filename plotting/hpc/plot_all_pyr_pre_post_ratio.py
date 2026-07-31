# -*- coding: utf-8 -*-
'''
Created on Fri 13 Dec 08:59:31 2024

plot profiles of run-onset PyrUp/PyrDown cells
plot profiles in good trials versus bad trials
plot profiles in ctrl trials versus stim trials, 26 Dec 2024

@author: Dinghao Luo

'''

#%% imports
import numpy as np
from scipy.stats import sem
import pandas as pd
import matplotlib.pyplot as plt
import sys
from pathlib import Path


#%% project imports and paths
repo_root = Path(__file__).resolve().parents[2]
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))
from plotting_functions import plot_violin_with_scatter, plot_bar_with_paired_scatter
from common_functions import normalise, normalise_to_all, mpl_formatting
import project_paths as pp
mpl_formatting()

HPC_STEM = pp.HPC_EPHYS_STEM
RUN_ONSET_STEM = pp.HPC_EPHYS_FIGURES_STEM / 'run_onset_response'

def save_figure(fig, filepath, **kwargs):
    filepath.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(filepath, **kwargs)


#%% load dataframe
print('loading dataframe...')
df = pd.read_pickle(HPC_STEM / 'hpc_all_profiles.pkl')
df_pyr = df[df['cell_identity']=='pyr']


#%% parameters
xaxis = np.arange(-1, 4, 1/1250)
PYRUP_CLASS = 'run-onset ON'
PYRDOWN_CLASS = 'run-onset OFF'


#%% sort dataframe by baseline pre-post ratios
df_sorted = df_pyr.sort_values(by='pre_post')


#%% matrices
pop_mat = df_sorted['prof_mean'].to_numpy()
pop_mat = np.asarray([normalise(cell[2500:2500+5*1250]) for cell in pop_mat])

pyrup = df_sorted[df_sorted['class']==PYRUP_CLASS]
pyrdown = df_sorted[df_sorted['class']==PYRDOWN_CLASS]


#%% overall plot
fig, ax = plt.subplots(figsize=(2.4,1.9))

gim = ax.imshow(pop_mat, aspect='auto', cmap='Greys', interpolation='sinc',
                extent=(-1, 4, 0, pop_mat.shape[0]))
ax.set(title=f'{pyrup.shape[0]} PyrUp, {pyrdown.shape[0]} PyrDown',
       xlabel='time from run-onset (s)',
       ylabel='cell #')

plt.colorbar(gim, shrink=.5, ticks=[0, 1])

for ext in ['.png', '.pdf']:
    save_figure(
        fig,
        RUN_ONSET_STEM / f'all_run_onset_Greys{ext}',
        dpi=300,
        bbox_inches='tight'
        )

fig, ax = plt.subplots(figsize=(2.4,1.9))

cw = ax.imshow(pop_mat, aspect='auto', cmap='coolwarm', interpolation='sinc',
          extent=(-1, 4, 0, pop_mat.shape[0]))
ax.set(title=f'{pyrup.shape[0]} PyrUp, {pyrdown.shape[0]} PyrDown',
       xlabel='time from run-onset (s)',
       ylabel='cell #')

plt.colorbar(cw, shrink=.5, ticks=[0, 1])

for ext in ['.png', '.pdf']:
    save_figure(
        fig,
        RUN_ONSET_STEM / f'all_run_onset_coolwarm{ext}',
        dpi=300,
        bbox_inches='tight'
        )


#%% extract good and bad mean profiles
pyrdown_mat = pyrdown['prof_mean'].to_numpy()
pyrdown_mat = np.asarray([normalise(cell[2500:2500+5*1250]) for cell in pyrdown_mat])

pyrup_mat = pyrup['prof_mean'].to_numpy()
pyrup_mat = np.asarray([normalise(cell[2500:2500+5*1250]) for cell in pyrup_mat])


#%% plotting  (pyrup and pyrdown)
fig, axs = plt.subplots(1, 2, figsize=(5,2))
plt.subplots_adjust(wspace=.5)

axs[0].imshow(pyrup_mat, aspect='auto', cmap='Greys', interpolation='sinc',
              extent=(-1, 4, 0, pyrup_mat.shape[0]))
axs[1].imshow(pyrdown_mat, aspect='auto', cmap='Greys', interpolation='sinc',
              extent=(-1, 4, 0, pyrdown_mat.shape[0]))

axs[0].set(title='PyrUp')
axs[1].set(title='PyrDown')
for ax in axs:
    ax.set(xlabel='time from run-onset (s)',
           ylabel='cell #')

for ext in ['.png', '.pdf']:
    save_figure(
        fig,
        RUN_ONSET_STEM / f'all_run_onset_PyrUp_PyrDown_Greys{ext}',
        dpi=300,
        bbox_inches='tight'
        )

fig, axs = plt.subplots(1, 2, figsize=(5,2))
plt.subplots_adjust(wspace=.5)

axs[0].imshow(pyrup_mat, aspect='auto', cmap='coolwarm', interpolation='sinc',
              extent=(-1, 4, 0, pyrup_mat.shape[0]))
axs[1].imshow(pyrdown_mat, aspect='auto', cmap='coolwarm', interpolation='sinc',
              extent=(-1, 4, 0, pyrdown_mat.shape[0]))

axs[0].set(title='PyrUp')
axs[1].set(title='PyrDown')
for ax in axs:
    ax.set(xlabel='time from run-onset (s)',
           ylabel='cell #')

for ext in ['.png', '.pdf']:
    save_figure(
        fig,
        RUN_ONSET_STEM / f'all_run_onset_PyrUp_PyrDown_coolwarm{ext}',
        dpi=300,
        bbox_inches='tight'
        )


#%% good v bad
# normalise each cell's good and bad trial profiles by the pooled profile;
# this assumes each cell's output in the CA1 network is tuned to its spike rate
pyrup_good_index = [
    np.mean(profile[4375:5625]) / np.mean(profile[1875:3125])
    for profile in pyrup['prof_good_mean']
    ]
pyrup_bad_index = [
    np.mean(profile[4375:5625]) / np.mean(profile[1875:3125])
    for profile in pyrup['prof_bad_mean']
    ]
pyrdown_good_index = [
    np.mean(profile[1875:3125]) / np.mean(profile[4375:5625])
    for profile in pyrdown['prof_good_mean']
    ]
pyrdown_bad_index = [
    np.mean(profile[1875:3125]) / np.mean(profile[4375:5625])
    for profile in pyrdown['prof_bad_mean']
    ]

pyrup_empty_mask = [good.size!=0 and bad.size!=0 for good, bad
                 in zip(pyrup['prof_good_mean'].to_numpy(), pyrup['prof_bad_mean'].to_numpy())]

pyrup_good_mat = [cell[2500:2500+1250*5] for cell, valid
               in zip(pyrup['prof_good_mean'].to_numpy(), pyrup_empty_mask)
               if valid]
pyrup_bad_mat = [cell[2500:2500+1250*5] for cell, valid
              in zip(pyrup['prof_bad_mean'].to_numpy(), pyrup_empty_mask)
              if valid]

for i in range(len(pyrup_good_mat)):
    temp_pool = np.concatenate((pyrup_good_mat[i], pyrup_bad_mat[i]))
    pyrup_good_mat[i] = normalise_to_all(pyrup_good_mat[i], temp_pool)
    pyrup_bad_mat[i] = normalise_to_all(pyrup_bad_mat[i], temp_pool)

pyrup_good_mean, pyrup_good_error = np.nanmean(pyrup_good_mat, axis=0), sem(pyrup_good_mat, axis=0)
pyrup_bad_mean, pyrup_bad_error = np.nanmean(pyrup_bad_mat, axis=0), sem(pyrup_bad_mat, axis=0)

pyrdown_empty_mask = [good.size!=0 and bad.size!=0 for good, bad
                  in zip(pyrdown['prof_good_mean'].to_numpy(), pyrdown['prof_bad_mean'].to_numpy())]

pyrdown_good_mat = [cell[2500:2500+1250*5] for cell, valid
                in zip(pyrdown['prof_good_mean'].to_numpy(), pyrdown_empty_mask)
                if valid]
pyrdown_bad_mat = [cell[2500:2500+1250*5] for cell, valid
               in zip(pyrdown['prof_bad_mean'].to_numpy(), pyrdown_empty_mask)
               if valid]

for i in range(len(pyrdown_good_mat)):
    temp_pool = np.concatenate((pyrdown_good_mat[i], pyrdown_bad_mat[i]))
    pyrdown_good_mat[i] = normalise_to_all(pyrdown_good_mat[i], temp_pool)
    pyrdown_bad_mat[i] = normalise_to_all(pyrdown_bad_mat[i], temp_pool)

pyrdown_good_mean, pyrdown_good_error = np.nanmean(pyrdown_good_mat, axis=0), sem(pyrdown_good_mat, axis=0)
pyrdown_bad_mean, pyrdown_bad_error = np.nanmean(pyrdown_bad_mat, axis=0), sem(pyrdown_bad_mat, axis=0)


#%% plot good v bad profiles
fig, axs = plt.subplots(1, 2, figsize=(3.6,1.2))
plt.subplots_adjust(wspace=.5)

lg, = axs[0].plot(xaxis, pyrup_good_mean, c='firebrick', linewidth=1, zorder=10)
axs[0].fill_between(xaxis, pyrup_good_mean+pyrup_good_error,
                           pyrup_good_mean-pyrup_good_error,
                    color='firebrick', edgecolor='none', alpha=.35, zorder=10)

lb, = axs[0].plot(xaxis, pyrup_bad_mean, c='grey', linewidth=1)
axs[0].fill_between(xaxis, pyrup_bad_mean+pyrup_bad_error,
                           pyrup_bad_mean-pyrup_bad_error,
                    color='grey', edgecolor='none', alpha=.35)

axs[0].legend([lg, lb], ['good trials', 'bad trials'], frameon=False, fontsize=5)

lg, = axs[1].plot(xaxis, pyrdown_good_mean, c='purple', linewidth=1, zorder=10)
axs[1].fill_between(xaxis, pyrdown_good_mean+pyrdown_good_error,
                           pyrdown_good_mean-pyrdown_good_error,
                    color='purple', edgecolor='none', alpha=.35, zorder=10)

lb, = axs[1].plot(xaxis, pyrdown_bad_mean, c='grey', linewidth=1)
axs[1].fill_between(xaxis, pyrdown_bad_mean+pyrdown_bad_error,
                           pyrdown_bad_mean-pyrdown_bad_error,
                    color='grey', edgecolor='none', alpha=.35)

axs[1].legend([lg, lb], ['good trials', 'bad trials'], frameon=False, fontsize=5)

axs[0].set(ylim=(.2,.5))
axs[1].set(ylim=(.24,.39))
for ax in axs:
    ax.set(xlabel='time from run-onset (s)', xticks=(0,2,4),
           ylabel='norm. spike rate')
    for s in ['top', 'right']:
        ax.spines[s].set_visible(False)

for ext in ['.png', '.pdf']:
    save_figure(
        fig,
        RUN_ONSET_STEM / f'all_PyrUp_PyrDown_good_bad{ext}',
        dpi=300,
        bbox_inches='tight'
        )


#%% filter and plot indices
outlier_mask = [1 if
                not (0 < good < 10) or
                not (0 < bad < 10)
                else 0
                for good, bad in zip(pyrup_good_index, pyrup_bad_index)]
pyrup_good_index_filt, pyrup_bad_index_filt = zip(*[(good, bad) for i, (good, bad)
                                            in enumerate(zip(pyrup_good_index, pyrup_bad_index))
                                            if outlier_mask[i]==0])
(RUN_ONSET_STEM / 'all_PyrUp_good_bad_index_violinplot').parent.mkdir(parents=True, exist_ok=True)
plot_violin_with_scatter(pyrup_good_index_filt, pyrup_bad_index_filt, 'firebrick', 'grey',
                         xticklabels=('good\ntrials', 'bad\ntrials'),
                         ylabel='response index',
                         showscatter=False, plot_statistics=True,
                         ylim=(0, 5),
                         figsize=(1.4, 2),
                         save=True,
                         savepath=str(RUN_ONSET_STEM / 'all_PyrUp_good_bad_index_violinplot'))

outlier_mask = [1 if
                not (0 < good < 10) or
                not (0 < bad < 10)
                else 0
                for good, bad in zip(pyrdown_good_index, pyrdown_bad_index)]
pyrdown_good_index_filt, pyrdown_bad_index_filt = zip(*[(good, bad) for i, (good, bad)
                                            in enumerate(zip(pyrdown_good_index, pyrdown_bad_index))
                                            if outlier_mask[i]==0])
(RUN_ONSET_STEM / 'all_PyrDown_good_bad_index_violinplot').parent.mkdir(parents=True, exist_ok=True)
plot_violin_with_scatter(pyrdown_good_index_filt, pyrdown_bad_index_filt, 'purple', 'grey',
                         xticklabels=('good\ntrials', 'bad\ntrials'),
                         ylabel='response index',
                         showscatter=False, plot_statistics=True,
                         ylim=(0, 5),
                         figsize=(1.4, 2),
                         save=True,
                         savepath=str(RUN_ONSET_STEM / 'all_PyrDown_good_bad_index_violinplot'))


#%% extract separate experiments
df_HPCLC = df_pyr[df_pyr['rectype']=='HPCLC']
df_HPCLCterm = df_pyr[df_pyr['rectype']=='HPCLCterm']


#%% plot HPCLC
count_ctrl_pyrup = (df_HPCLC['pre_post_ctrl_MATLAB']<.8).sum()
count_stim_pyrup = (df_HPCLC['pre_post_stim_MATLAB']<.8).sum()

df_sorted_ctrl = df_HPCLC.sort_values(by='pre_post_ctrl')
df_sorted_stim = df_HPCLC.sort_values(by='pre_post_stim')

ctrl_profiles = df_sorted_ctrl['prof_ctrl_mean'].to_numpy()
ctrl_profiles = [normalise(cell[2500:2500+5*1250]) for cell in ctrl_profiles]

stim_profiles = df_sorted_stim['prof_stim_mean'].to_numpy()
stim_profiles = [normalise(cell[2500:2500+5*1250]) for cell in stim_profiles]

fig, axs = plt.subplots(1, 2, figsize=(4.4,2.5))
axs[0].imshow(ctrl_profiles, extent=(-1, 4, 0, len(ctrl_profiles)),
              aspect='auto', cmap='Greys', interpolation='none')
axs[0].set(title='ctrl.')
axs[1].imshow(stim_profiles, extent=(-1, 4, 0, len(ctrl_profiles)),
              aspect='auto', cmap='Greys', interpolation='none')
axs[1].set(title='stim.')

for i in [0,1]:
    axs[i].set(xlabel='time from run-onset (s)',
               ylabel='cell #')

fig.suptitle(f'ctrl. PyrUp: {count_ctrl_pyrup}, stim. PyrUp: {count_stim_pyrup}')

fig.tight_layout()

for ext in ['.png', '.pdf']:
    save_figure(fig, RUN_ONSET_STEM / f'HPCLC_run_onset_ctrl_stim{ext}',
                dpi=300)


#%% percentage of PyrUp cells within sessions
perc_pyrup_ctrl_per_session = df_HPCLC.groupby('recname')['pre_post_ctrl_MATLAB'].apply(
    lambda x: (x < 0.8).mean() * 100  # mean of booleans = percentage of True values
).to_numpy()
perc_pyrup_stim_per_session = df_HPCLC.groupby('recname')['pre_post_stim_MATLAB'].apply(
    lambda x: (x < 0.8).mean() * 100
).to_numpy()

# pyrup plot
fig, ax = plt.subplots(figsize=(1.6, 2.0), dpi=300)
plot_bar_with_paired_scatter(ax,
                             perc_pyrup_ctrl_per_session,
                             perc_pyrup_stim_per_session,
                             colors=('gray','firebrick'),
                             ylim=(15, 75),
                             title='PyrUp', ylabel='% PyrUp cells', xticklabels=('ctrl.','stim.'))
fig.tight_layout()
for ext in ['.png', '.pdf']:
    save_figure(fig, RUN_ONSET_STEM / f'HPCLC_PyrUp_perc_sess{ext}',
                dpi=300, bbox_inches='tight')

perc_pyrdown_ctrl_per_session = df_HPCLC.groupby('recname')['pre_post_ctrl_MATLAB'].apply(
    lambda x: (x > 1.25).mean() * 100  # mean of booleans = percentage of True values
).to_numpy()
perc_pyrdown_stim_per_session = df_HPCLC.groupby('recname')['pre_post_stim_MATLAB'].apply(
    lambda x: (x > 1.25).mean() * 100
).to_numpy()

# pyrdown plot
fig, ax = plt.subplots(figsize=(1.6, 2.0), dpi=300)
plot_bar_with_paired_scatter(ax,
                             perc_pyrdown_ctrl_per_session,
                             perc_pyrdown_stim_per_session,
                             colors=('gray','purple'),
                             ylim=(0, 50),
                             title='PyrDown', ylabel='% PyrDown cells', xticklabels=('ctrl.','stim.'))
fig.tight_layout()
for ext in ['.png', '.pdf']:
    save_figure(fig, RUN_ONSET_STEM / f'HPCLC_PyrDown_perc_sess{ext}',
                dpi=300, bbox_inches='tight')


#%% plot HPCLCterm
count_ctrl_term_pyrup = (df_HPCLCterm['pre_post_ctrl_MATLAB']<.8).sum()
count_stim_term_pyrup = (df_HPCLCterm['pre_post_stim_MATLAB']<.8).sum()

df_sorted_ctrl_term = df_HPCLCterm.sort_values(by='pre_post_ctrl')
df_sorted_stim_term = df_HPCLCterm.sort_values(by='pre_post_stim')

ctrl_term_profiles = df_sorted_ctrl_term['prof_ctrl_mean'].to_numpy()
ctrl_term_profiles = [normalise(cell[2500:2500+5*1250]) for cell in ctrl_term_profiles]

stim_term_profiles = df_sorted_stim_term['prof_stim_mean'].to_numpy()
stim_term_profiles = [normalise(cell[2500:2500+5*1250]) for cell in stim_term_profiles]

fig, axs = plt.subplots(1, 2, figsize=(4.4,2.5))
axs[0].imshow(ctrl_term_profiles, extent=(-1, 4, 0, len(ctrl_term_profiles)),
              aspect='auto', cmap='Greys', interpolation='none')
axs[0].set(title='ctrl.')
axs[1].imshow(stim_term_profiles, extent=(-1, 4, 0, len(ctrl_term_profiles)),
              aspect='auto', cmap='Greys', interpolation='none')
axs[1].set(title='stim.')

for i in [0,1]:
    axs[i].set(xlabel='time from run-onset (s)',
               ylabel='cell #')

fig.suptitle(f'ctrl. PyrUp: {count_ctrl_term_pyrup}, stim. PyrUp: {count_stim_term_pyrup}')

fig.tight_layout()

for ext in ['.png', '.pdf']:
    save_figure(fig, RUN_ONSET_STEM / f'HPCLCterm_run_onset_ctrl_stim{ext}',
                dpi=300)


#%% percentage of PyrUp cells within sessions
perc_pyrup_ctrl_term_per_session = df_HPCLCterm.groupby('recname')['pre_post_ctrl_MATLAB'].apply(
    lambda x: (x < 0.8).mean() * 100  # mean of booleans = percentage of True values
).to_numpy()
perc_pyrup_stim_term_per_session = df_HPCLCterm.groupby('recname')['pre_post_stim_MATLAB'].apply(
    lambda x: (x < 0.8).mean() * 100
).to_numpy()

# pyrup plot
fig, ax = plt.subplots(figsize=(1.6, 2.0), dpi=300)
plot_bar_with_paired_scatter(ax,
                             perc_pyrup_ctrl_term_per_session,
                             perc_pyrup_stim_term_per_session,
                             colors=('gray','firebrick'),
                             ylim=(4, 62),
                             title='PyrUp', ylabel='% PyrUp cells', xticklabels=('ctrl.','stim.'))
fig.tight_layout()
for ext in ['.png', '.pdf']:
    save_figure(fig, RUN_ONSET_STEM / f'HPCLCterm_PyrUp_perc_sess{ext}',
                dpi=300, bbox_inches='tight')

perc_pyrdown_ctrl_term_per_session = df_HPCLCterm.groupby('recname')['pre_post_ctrl_MATLAB'].apply(
    lambda x: (x > 1.25).mean() * 100  # mean of booleans = percentage of True values
).to_numpy()
perc_pyrdown_stim_term_per_session = df_HPCLCterm.groupby('recname')['pre_post_stim_MATLAB'].apply(
    lambda x: (x > 1.25).mean() * 100
).to_numpy()

# pyrdown plot
fig, ax = plt.subplots(figsize=(1.6, 2.0), dpi=300)
plot_bar_with_paired_scatter(ax,
                             perc_pyrdown_ctrl_term_per_session,
                             perc_pyrdown_stim_term_per_session,
                             colors=('gray','purple'),
                             ylim=(0, 35),
                             title='PyrDown', ylabel='% PyrDown cells', xticklabels=('ctrl.','stim.'))
fig.tight_layout()
for ext in ['.png', '.pdf']:
    save_figure(fig, RUN_ONSET_STEM / f'HPCLCterm_PyrDown_perc_sess{ext}',
                dpi=300, bbox_inches='tight')
