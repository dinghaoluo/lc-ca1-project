# -*- coding: utf-8 -*-
'''
Created on Thu Oct 24 18:22:01 2024

plot single cell waveform

@author: Dinghao Luo

'''

#%% imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
from pathlib import Path

repo_root = Path(__file__).resolve().parents[2]
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))
from common_functions import normalise
import project_paths as pp

# plotting parameters
import matplotlib
plt.rcParams['font.family'] = 'Arial'
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42


#%% cell properties
cell_prop = pd.read_pickle(pp.LC_ALL_STEM / 'LC_all_single_cell_properties.pkl')
tag_list = []; put_list = []
for clu in cell_prop.index:
    tg = cell_prop['tagged'][clu]
    pt = cell_prop['putative'][clu]

    if tg:
        tag_list.append(clu)
    if pt:
        put_list.append(clu)


#%% load waveforms
all_cells = np.load(
    pp.LC_ALL_STEM / 'LC_all_waveforms.npy',
    allow_pickle=True
    ).item()

#%% plotting
tot_clu = len(all_cells)
time_ax = np.arange(32)

all_cells_list = list(all_cells.items())
waveform_output_stem = pp.LC_ALL_FIGURES_STEM / 'single_cell_waveform'
waveform_output_stem.mkdir(parents=True, exist_ok=True)

for clu in all_cells_list:
    cluname = clu[0]

    spk = clu[1][0]
    spk_norm = normalise(spk)
    scaling_factor = spk_norm[0]/spk[0]
    sem = clu[1][1] * scaling_factor

    suffix = ''
    if cluname in tag_list: suffix=' tgd'
    if cluname in put_list: suffix=' put'

    fig, ax = plt.subplots(figsize=(1,1))
    ax.set(title=cluname+suffix)
    ax.plot(time_ax, spk_norm, color='k')
    ax.fill_between(time_ax, spk_norm+sem, spk_norm-sem, color='k', alpha=.25)
    ax.axis('off')

    plt.show()

    fig.savefig(waveform_output_stem / f'{cluname + suffix}.png',
                dpi=200, bbox_inches='tight')
    fig.savefig(waveform_output_stem / f'{cluname + suffix}.pdf',
                bbox_inches='tight')

    plt.close(fig)
