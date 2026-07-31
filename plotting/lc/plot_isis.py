# -*- coding: utf-8 -*-
'''
Created on Fri 7 Mar 16:05:21 2025

plot ISIs of single neurones

@author: Dinghao Luo
'''

#%% imports
import numpy as np
import matplotlib.pyplot as plt
import sys
from pathlib import Path

repo_root = Path(__file__).resolve().parents[2]
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))
from common_functions import mpl_formatting
import project_paths as pp
mpl_formatting()

import rec_list
paths = rec_list.pathLC

single_cell_isi_stem = pp.LC_EPHYS_FIGURES_STEM / 'single_cell_ISIs'
single_cell_isi_stem.mkdir(parents=True, exist_ok=True)


#%% plot single-cell ISIs
for path in paths:
    recname = Path(path).name
    print(f'plotting {recname}')

    sess_folder = pp.LC_EPHYS_STEM / 'all_sessions' / recname

    ISI_dict = np.load(
        sess_folder / f'{recname}_all_ISIs.npy',
        allow_pickle=True
        ).item()

    identity_dict = np.load(
        sess_folder / f'{recname}_all_identities.npy',
        allow_pickle=True
        ).item()

    for clu, ISIs in ISI_dict.items():

        fig, ax = plt.subplots(figsize=(3,2))

        cluname = f'{clu} tagged' if identity_dict[clu] else clu

        ax.hist([ISI/20000 for ISI in ISIs], bins=np.arange(0, 2, .01),
                color='k')

        ax.set(title=cluname,
               xlabel='ISI (s)',
               ylabel='frequency')

        fig.savefig(
            single_cell_isi_stem / f'{cluname}.png',
            dpi=300,
            bbox_inches='tight')

        plt.close(fig)
