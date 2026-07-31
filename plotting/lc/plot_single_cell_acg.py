# -*- coding: utf-8 -*-
'''
Created on Thu Oct  3 16:36:21 2024

plot ACGs for single cells

@author: Dinghao Luo
'''

#%% imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import sys

# plotting parameters
import matplotlib
plt.rcParams['font.family'] = 'Arial'
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

repo_root = Path(__file__).resolve().parents[2]
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))
import project_paths as pp


#%% load cell properties
cell_prop = pd.read_pickle(pp.LC_ALL_STEM / 'LC_all_single_cell_properties.pkl')
tag_list = [clu for clu in cell_prop.index if cell_prop['tagged'][clu]]
put_list = [clu for clu in cell_prop.index if cell_prop['putative'][clu]]


#%% load ACGs
acgs = np.load(pp.LC_ALL_STEM / 'LC_all_acg_baseline.npy',
               allow_pickle=True).item()
output_stem = pp.LC_ALL_FIGURES_STEM / 'single_cell_ACG'
output_stem.mkdir(parents=True, exist_ok=True)


#%% main
xaxis = np.arange(-200, 200, 1)
for cluname, acg in acgs.items():
    suffix = ''
    if cluname in tag_list: suffix=' tgd'
    if cluname in put_list: suffix=' put'
    fig, ax = plt.subplots(figsize=(1,1))
    ax.plot(xaxis, acg[9800:10200], color='k')
    for s in ['left', 'top', 'right']: ax.spines[s].set_visible(False)
    ax.set(xlabel='lag (ms)', xticks=(-200, 0, 200),
           yticks=[], ylabel='', yticklabels='',
           title=cluname+suffix)
    for ext in ['png', 'pdf']:
        fig.savefig(output_stem / f'{cluname + suffix}.{ext}',
                    dpi=200, bbox_inches='tight')
    plt.close(fig)
