# -*- coding: utf-8 -*-
'''
Created on 20 May 2026

plot the pooled relationship between previous reward and first-lick timing

@author: Dinghao Luo
'''

#%% imports
from pathlib import Path
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import linregress

repo_root = Path(__file__).resolve().parents[2]
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))

from common_functions import mpl_formatting
import project_paths as pp
mpl_formatting()


#%% paths
payload_path = (
    pp.BEHAVIOUR_STEM
    / 'first_lick_since_last_reward'
    / 'time_since_reward_vs_first_lick_payload.npy'
    )
figure_stem = (
    pp.BEHAVIOUR_FIGURES_STEM
    / 'first_lick_since_last_reward'
    / 'time_since_reward_vs_first_lick'
    )

#%% load data
payload = np.load(payload_path, allow_pickle=True).item()
x = payload['x']
y = payload['y']

#%% plotting
fig, ax = plt.subplots(figsize=(4, 4))
ax.scatter(x, y, s=20, alpha=0.7, edgecolor='none')

fit = linregress(x, y)
xx = np.linspace(x.min(), x.max(), 200)
yy = fit.intercept + fit.slope * xx
ax.plot(xx, yy, linewidth=2)

stats_txt = (
    f'n = {len(x)}\n'
    f'y = {fit.slope:.3f}x + {fit.intercept:.3f}\n'
    f'r = {fit.rvalue:.3f}\n'
    f'p = {fit.pvalue:.3g}'
)
ax.text(
    0.02,
    0.98,
    stats_txt,
    transform=ax.transAxes,
    va='top',
    ha='left',
    fontsize=12,
    bbox=dict(boxstyle='round', facecolor='white', alpha=0.85, linewidth=0.5),
)

ax.set_xlabel('time since last reward (s)')
ax.set_ylabel('first lick time (s)')
ax.set_title('time since last reward vs first lick time')

plt.tight_layout()
figure_stem.parent.mkdir(parents=True, exist_ok=True)
for ext in ['.png', '.pdf']:
    fig.savefig(str(figure_stem) + ext, dpi=300, bbox_inches='tight')
plt.close(fig)
