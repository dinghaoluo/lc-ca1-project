# -*- coding: utf-8 -*-
'''
Created on Thu Feb 23 16:04:32 2023
Originally named egsess_lick.py

plot lick curves (beh example)

@author: Dinghao Luo
'''

#%% imports
import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio
from pathlib import Path
import sys

# for Illustrator
import matplotlib
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
plt.rcParams['font.family'] = 'Arial'

repo_root = Path(__file__).resolve().parents[3]
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))
import project_paths as pp
pp.LC_FIGURES_STEM.mkdir(parents=True, exist_ok=True)


#%% main
rec_used = 'A014-20211201'
print('recording used: {}'.format(rec_used))

recdata_path = pp.BEHAV_DATA_ANALYSIS_ROOT / f'ANMD{rec_used[1:4]}' / rec_used / f'{rec_used}_compSess.mat'

recdata = sio.loadmat(recdata_path)
all_licks = recdata['sessDataLick']

sess_used = 4
print('session used: {}'.format(sess_used))

# access the depths of lick data (how is this so deep jesus)
sess_licks_mean = all_licks['meanRun'][0][0][0,sess_used].reshape(-1)[:220]
sess_licks_sem = all_licks['SEMRun'][0][0][0,sess_used].reshape(-1)[:220]


#%% plotting
fig, ax = plt.subplots(figsize=(2.2,1.8))
ax.set(title='avg lick profile',
       xlabel='distance (cm)', ylabel='lick rate (Hz)',
       xlim=(30, 220), ylim=(0, 5.5))
for p in ['right', 'top']:
    ax.spines[p].set_visible(False)

xaxis = np.arange(220)
ax.plot(xaxis, sess_licks_mean, 'k')
ax.fill_between(xaxis,
                sess_licks_mean-sess_licks_sem,
                sess_licks_mean+sess_licks_sem,
                color='grey', alpha=0.25,
                edgecolor='none')

plt.show()
fig.savefig(pp.LC_FIGURES_STEM / 'egsess_lick.png',
            dpi=300,
            bbox_inches='tight')
fig.savefig(pp.LC_FIGURES_STEM / 'egsess_lick.pdf',
            bbox_inches='tight')

rew_ln, = ax.plot([180, 180], [0, 10], color='limegreen', alpha=.45)
ax.legend([rew_ln], ['reward'], loc='upper left', frameon=False, fontsize=8)

plt.show()
fig.savefig(pp.LC_FIGURES_STEM / 'egsess_lick_w_rew.png',
            dpi=300,
            bbox_inches='tight')
fig.savefig(pp.LC_FIGURES_STEM / 'egsess_lick_w_rew.pdf',
            bbox_inches='tight')
