# -*- coding: utf-8 -*-
'''
Created on 20 May 2026
Modified on 24 June 2026

plot session and pooled speed and lick profiles

@author: Dinghao Luo
'''

#%% imports
from pathlib import Path
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

repo_root = Path(__file__).resolve().parents[2]
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))

from common_functions import mpl_formatting
import project_paths as pp
mpl_formatting()


#%% paths
payload_path = pp.BEHAVIOUR_SESSION_PROFILES_STEM / 'speed_lick_profiles_payload.npy'
all_sess_stem = pp.BEHAVIOUR_SESSION_PROFILES_FIGURES_STEM

#%% load data
payload = np.load(payload_path, allow_pickle=True).item()

#%% session plots
for recname, session_payload in payload['sessions'].items():
    exp_stem = all_sess_stem / session_payload['exp_name'].lower()

    fig, axs = plt.subplots(1, 3, figsize=(7.4, 2.4))
    speed_arr_plot = session_payload['speed_arr'][:, :len(session_payload['speed_time_axis'])]
    speed_arr_axis = session_payload['speed_time_axis'][:speed_arr_plot.shape[1]]

    axs[0].plot(session_payload['speed_time_axis'], session_payload['mean_speeds'],
                c='navy', zorder=10)
    axs[0].fill_between(
        session_payload['speed_time_axis'],
        session_payload['mean_speeds'] + session_payload['sem_speeds'],
        session_payload['mean_speeds'] - session_payload['sem_speeds'],
        color='navy', alpha=.2, edgecolor='none', zorder=10
        )
    axs[0].plot(
        speed_arr_axis,
        speed_arr_plot.T,
        c='grey', lw=.5, alpha=.05, zorder=1
        )
    axs[0].set(xlabel='time from run-onset (s)',
               xticks=[0, 2, 4],
               ylabel='speed (cm/s)',
               ylim=(0, np.nanmax(session_payload['speed_arr'][:, :5*1000]) + 1))

    axs[1].plot(session_payload['speed_distance_axis'],
                session_payload['mean_speeds_distances'],
                c='navy', zorder=10)
    axs[1].fill_between(
        session_payload['speed_distance_axis'],
        session_payload['mean_speeds_distances'] + session_payload['sem_speeds_distances'],
        session_payload['mean_speeds_distances'] - session_payload['sem_speeds_distances'],
        color='navy', alpha=.2, edgecolor='none', zorder=10
        )
    axs[1].plot(
        session_payload['speed_distance_axis'],
        session_payload['speed_distances'].T,
        c='grey', lw=.5, alpha=.05, zorder=1
        )
    axs[1].set(xlabel='distance (cm)',
               xticks=[0, 90, 180],
               ylabel='speed (cm/s)',
               ylim=(0, np.nanmax(session_payload['speed_distances'][:, :5*1000]) + 1))

    axs[2].plot(session_payload['lick_distance_axis'],
                session_payload['mean_lick_maps'],
                c='orchid', lw=1)
    axs[2].fill_between(
        session_payload['lick_distance_axis'],
        session_payload['mean_lick_maps'] + session_payload['sem_lick_maps'],
        session_payload['mean_lick_maps'] - session_payload['sem_lick_maps'],
        color='orchid', alpha=.2, edgecolor='none'
        )
    axs[2].axvspan(
        179.5, 220,
        edgecolor='none', facecolor='darkgreen', alpha=.15, zorder=10
        )
    ymax = np.nanmax(session_payload['mean_lick_maps'] + session_payload['sem_lick_maps']) * 1.01
    axs[2].set(xlabel='Distance (cm)',
               xticks=[90, 180],
               ylabel='Licks',
               ylim=(0, ymax))
    for i in range(3):
        for s in ['top', 'right']:
            axs[i].spines[s].set_visible(False)

    fig.suptitle(recname)
    fig.tight_layout()
    output_stem = exp_stem / f'{recname}_dist'
    output_stem.parent.mkdir(parents=True, exist_ok=True)
    for ext in ['.png', '.pdf']:
        fig.savefig(str(output_stem) + ext, dpi=300, bbox_inches='tight')
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(3.4, 2.3))

    ax.plot(session_payload['speed_distance_axis'],
            session_payload['mean_speeds_distances'],
            c='navy', zorder=10)
    ax.fill_between(
        session_payload['speed_distance_axis'],
        session_payload['mean_speeds_distances'] + session_payload['sem_speeds_distances'],
        session_payload['mean_speeds_distances'] - session_payload['sem_speeds_distances'],
        color='navy', alpha=.2, edgecolor='none'
        )

    axt = ax.twinx()
    axt.plot(session_payload['lick_distance_axis'],
             session_payload['mean_lick_maps'],
             c='orchid')
    axt.fill_between(
        session_payload['lick_distance_axis'],
        session_payload['mean_lick_maps'] + session_payload['sem_lick_maps'],
        session_payload['mean_lick_maps'] - session_payload['sem_lick_maps'],
        color='orchid', alpha=.2, edgecolor='none'
        )

    ax.axvspan(179.5, 220,
               edgecolor='none', facecolor='darkgreen', alpha=.15)

    ax.spines['top'].set_visible(False)
    axt.spines['top'].set_visible(False)

    ax.set(xlabel='Distance (cm)',
           xlim=(0, 220),
           ylabel='Speed (cm/s)',
           ylim=(0, np.nanmax(session_payload['mean_speeds_distances']) * 1.1))
    axt.set(xlim=(0, 220),
            ylabel='Licks',
            ylim=(0, np.nanmax(session_payload['mean_lick_maps']) * 1.1))

    fig.tight_layout()
    output_stem = exp_stem / f'{recname}_overlay_dist'
    output_stem.parent.mkdir(parents=True, exist_ok=True)
    for ext in ['.png', '.pdf']:
        fig.savefig(str(output_stem) + ext, dpi=300, bbox_inches='tight')
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(3.4, 2.3))

    ax.plot(session_payload['speed_time_axis'], session_payload['mean_speeds'],
            c='navy', zorder=10)
    ax.fill_between(
        session_payload['speed_time_axis'],
        session_payload['mean_speeds'] + session_payload['sem_speeds'],
        session_payload['mean_speeds'] - session_payload['sem_speeds'],
        color='navy', alpha=.2, edgecolor='none'
        )

    axt = ax.twinx()
    axt.plot(session_payload['lick_times_axis'], session_payload['mean_lick_times_maps'],
             c='orchid')
    axt.fill_between(
        session_payload['lick_times_axis'],
        session_payload['mean_lick_times_maps'] + session_payload['sem_lick_times_maps'],
        session_payload['mean_lick_times_maps'] - session_payload['sem_lick_times_maps'],
        color='orchid', alpha=.2, edgecolor='none'
        )

    ax.spines['top'].set_visible(False)
    axt.spines['top'].set_visible(False)

    ax.set(xlabel='Time from run-onset (s)',
           xlim=(0, 5),
           ylabel='Speed (cm/s)',
           ylim=(0, np.nanmax(session_payload['mean_speed_times']) * 1.1))
    axt.set(xlim=(0, 5),
            ylabel='Lick rate (Hz)',
            ylim=(0, np.nanmax(session_payload['mean_lick_times_maps']) * 1.1))

    fig.tight_layout()
    output_stem = exp_stem / f'{recname}_overlay_time'
    output_stem.parent.mkdir(parents=True, exist_ok=True)
    for ext in ['.png', '.pdf']:
        fig.savefig(str(output_stem) + ext, dpi=300, bbox_inches='tight')
    plt.close(fig)

#%% pooled speed and lick trace
fig, ax = plt.subplots(figsize=(2.8, 2.3))

ax.plot(payload['speed_time_axis'], payload['mean_speed_group'],
        c='navy', zorder=10)
ax.fill_between(
    payload['speed_time_axis'],
    payload['mean_speed_group'] + payload['sem_speed_group'],
    payload['mean_speed_group'] - payload['sem_speed_group'],
    color='navy', alpha=.2, edgecolor='none'
    )

axt = ax.twinx()
axt.plot(payload['lick_times_axis'], payload['mean_lick_group'],
         c='orchid')
axt.fill_between(
    payload['lick_times_axis'],
    payload['mean_lick_group'] + payload['sem_lick_group'],
    payload['mean_lick_group'] - payload['sem_lick_group'],
    color='orchid', alpha=.2, edgecolor='none'
    )

ax.spines['top'].set_visible(False)
axt.spines['top'].set_visible(False)

ax.set(xlabel='Time from run-onset (s)',
       xlim=(0, 5),
       ylabel='speed (cm/s)',
       ylim=(0, np.nanmax(payload['mean_speed_group']) * 1.1))
axt.set(xlim=(0, 5),
        ylabel='lick rate (Hz)',
        ylim=(0, np.nanmax(payload['mean_lick_group']) * 1.1))

fig.suptitle('all sessions')
fig.tight_layout()
output_stem = all_sess_stem / 'mean_overlay_time'
output_stem.parent.mkdir(parents=True, exist_ok=True)
for ext in ['.png', '.pdf']:
    fig.savefig(str(output_stem) + ext, dpi=300, bbox_inches='tight')
plt.close(fig)

#%% predictive lick slope
all_real_slopes = payload['all_real_slopes']
ymax = payload['ymax']
mean_r = payload['mean_r']
sem_r = payload['sem_r']
tval = payload['tval']
p_t = payload['p_t']
wstat = payload['wstat']
p_w = payload['p_w']

fig, ax = plt.subplots(figsize=(1.6, 2.2))

parts = ax.violinplot(
    all_real_slopes,
    positions=[1],
    showmeans=False,
    showmedians=True,
    showextrema=False
    )
for pc in parts['bodies']:
    pc.set_facecolor('orchid')
    pc.set_edgecolor('none')
    pc.set_alpha(0.35)
parts['cmedians'].set_color('k')
parts['cmedians'].set_linewidth(1.2)

ax.scatter(np.ones(len(all_real_slopes)), all_real_slopes,
           color='orchid', ec='none', s=10, alpha=0.5, zorder=3)

ax.axhline(payload['shuf_mean'], color='grey', lw=1, ls='--')

ax.fill_between(
    [0, 2],
    payload['lower_95'],
    payload['upper_95'],
    color='grey', alpha=0.2, edgecolor='none', zorder=0
    )

ax.text(1, ymax + 0.05*(ymax - np.min(all_real_slopes)),
        f'{mean_r:.2f} ± {sem_r:.2f}',
        ha='center', va='bottom', fontsize=7, color='orchid')

ax.text(1, np.min(all_real_slopes) - 0.10*(ymax - np.min(all_real_slopes)),
        f't(1-samp)={tval:.2f}, p={p_t:.2e}\n'
        f'Wilcoxon={wstat:.2f}, p={p_w:.2e}',
        ha='center', va='top', fontsize=6.5, color='black')

ax.set(xlim=(0.5, 1.5), xticks=[],
       ylabel='Predictive lick slope (Hz/s)',
       title='Predictive lick slope')
ax.spines[['top', 'right', 'bottom']].set_visible(False)

plt.tight_layout()
output_stem = all_sess_stem / 'predictive_lick_slope_2to4_violinplot'
output_stem.parent.mkdir(parents=True, exist_ok=True)
for ext in ['.png', '.pdf']:
    fig.savefig(str(output_stem) + ext, dpi=300, bbox_inches='tight')
plt.close(fig)
