# -*- coding: utf-8 -*-
'''
Created on Thu Feb 23 16:31:35 2023
Originally named egsess_speed_passive_raphi.py

plot speed curve (beh example)

@author: Dinghao Luo
'''

#%% imports
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['font.family'] = 'Arial'
import scipy.io as sio
from scipy.stats import sem
from pathlib import Path
import sys

repo_root = Path(__file__).resolve().parents[3]
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))
import project_paths as pp
pp.LC_FIGURES_STEM.mkdir(parents=True, exist_ok=True)


#%% main
rec_used = 'Raphi A009-20190128-01'
print('recording used: {}'.format(rec_used))

recdata_path = (
    pp.RAPHAEL_ROOT / 'ANM009' / 'A009-20190128' / 'A009-20190128-01'
    / 'A009-20190128-01_DataStructure_mazeSection1_TrialType1_runSpeedDist_Run1.mat'
)

recdata = sio.loadmat(recdata_path)
all_speeds = recdata['speedOverDist']

# access the depths of lick data (how is this so deep jesus)
sess_speeds_mean = np.mean(all_speeds, axis=0)/10
sess_speeds_sem = sem(all_speeds, axis=0)/10


#%% plotting
fig, ax = plt.subplots(figsize=(2.2,1.8))
ax.set(title='avg. velocity profile',
       xlabel='distance (cm)', ylabel='velocity (cm/s)',
       xlim=(0, 180), ylim=(0, 80),
       yticks=[0, 50])
for p in ['right', 'top']:
    ax.spines[p].set_visible(False)

xaxis = np.arange(1801)/10
ax.plot(xaxis, sess_speeds_mean, 'k')
ax.fill_between(xaxis,
                sess_speeds_mean-sess_speeds_sem,
                sess_speeds_mean+sess_speeds_sem,
                color='grey', alpha=.25)

fig.savefig(pp.LC_FIGURES_STEM / 'egsess_speed_passive.png',
            dpi=300,
            bbox_inches='tight')
