# -*- coding: utf-8 -*-
'''
Created on 20 May 2026

plot run-aligned pupil traces and session-wise pre/post changes

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
payload_path = pp.PUPIL_TRACKING_OUTPUT_STEM / 'run_aligned_pupil_payload.npy'
output_stem = pp.PUPIL_TRACKING_FIGURES_STEM

#%% load data
payload = np.load(payload_path, allow_pickle=True).item()
output_stem.mkdir(parents=True, exist_ok=True)
single_session_stem = output_stem / 'single_sessions'
single_session_stem.mkdir(parents=True, exist_ok=True)

#%% single sessions
for recname, rec_payload in payload['single_sessions'].items():
    fig, ax = plt.subplots(figsize=(2.9,2.3))
    ax.set(xlabel='time (s)', xlim=(-1,4), ylabel='pupil size')
    ax.set_xticks([0,2,4])

    ax.plot(rec_payload['x'], rec_payload['avg_start'], 'k', lw=2)
    ax.fill_between(
        rec_payload['x'],
        rec_payload['avg_start'] + rec_payload['sem_start'],
        rec_payload['avg_start'] - rec_payload['sem_start'],
        color='k',
        edgecolor='none',
        alpha=.25
        )
    ax.set(title=recname)
    for p in ['right', 'top']:
        ax.spines[p].set_visible(False)

    fig.tight_layout()
    for ext in ['.png', '.pdf']:
        fig.savefig(
            single_session_stem / f'{recname}_run_aligned{ext}',
            dpi=500,
            bbox_inches='tight'
            )
    plt.close(fig)

#%% pooled trace
fig, ax = plt.subplots(figsize=(2.8, 2.3))
ax.plot(payload['x_axis'], payload['grand_avg'], 'k', lw=2)
ax.fill_between(
    payload['x_axis'],
    payload['grand_avg'] + payload['grand_sem'],
    payload['grand_avg'] - payload['grand_sem'],
    color='k',
    edgecolor='none',
    alpha=0.15
    )
ax.set(xlabel='Time from run onset (s)',
       ylabel='Norm. pupil size',
       title='Run-onset-aligned pupil trace',
       xlim=(-1, 4), xticks=[0,2,4])
for p in ['right', 'top']:
    ax.spines[p].set_visible(False)
fig.tight_layout()
for ext in ['.png', '.pdf']:
    fig.savefig(
        output_stem / f'mean_run_aligned_pupil{ext}',
        dpi=500,
        bbox_inches='tight'
        )
plt.close(fig)

#%% session-wise changes
session_real_deltas = payload['session_real_deltas']
q1 = payload['q1']
q3 = payload['q3']
iqr = payload['iqr']
ymax = payload['ymax']
mean_r = payload['mean_r']
sem_r = payload['sem_r']
tval = payload['tval']
p_t = payload['p_t']
wstat = payload['wstat']
p_w = payload['p_w']

fig, ax = plt.subplots(figsize=(1.6, 2.2))

parts = ax.violinplot(
    session_real_deltas,
    positions=[1],
    showmeans=False,
    showmedians=True,
    showextrema=False
    )
for pc in parts['bodies']:
    pc.set_facecolor('k')
    pc.set_edgecolor('none')
    pc.set_alpha(0.35)
parts['cmedians'].set_color('k')
parts['cmedians'].set_linewidth(1.2)

ax.scatter(np.ones(len(session_real_deltas)), session_real_deltas,
           color='k', ec='none', s=10, alpha=0.5, zorder=3)

ax.axhline(payload['shuf_mean'], color='grey', lw=1, ls='--')

ax.fill_between(
    [0, 2],
    payload['lower_95'],
    payload['upper_95'],
    color='grey', alpha=0.2, edgecolor='none', zorder=0
    )

ax.text(
    1,
    ymax + 0.05*(ymax - np.min(session_real_deltas)),
    f'{mean_r:.2f} ± {sem_r:.2f}',
    ha='center', va='bottom', fontsize=7, color='k'
    )

ax.text(
    1,
    q1 - 0.20 * (q3 - q1),
    f'IQR={iqr:.2f}',
    ha='center', va='top', fontsize=6.5, color='k'
    )

ax.text(
    1, np.min(session_real_deltas) - 0.10*(ymax - np.min(session_real_deltas)),
    f't(1-samp)={tval:.2f}, p={p_t:.2e}\n'
    f'Wilcoxon={wstat:.2f}, p={p_w:.2e}',
    ha='center', va='top', fontsize=6.5, color='k'
    )

ax.set(xlim=(0.5, 1.5), xticks=[],
       ylabel='Norm. pupil size change',
       title='Pupil size delta')
ax.spines[['top', 'right', 'bottom']].set_visible(False)

plt.tight_layout()
for ext in ['.png', '.pdf']:
    fig.savefig(
        output_stem / f'mean_run_aligned_deltas_violin{ext}',
        dpi=500,
        bbox_inches='tight'
        )
plt.close(fig)
