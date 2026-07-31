# -*- coding: utf-8 -*-
'''
Created on Mon 10 July 14:57:34 2023
Modified on Mon 23 Dec 2024:
    - process all paths

summarise all cell waveforms from HPC recordings

@author: Dinghao Luo
'''

#%% imports
from pathlib import Path
import sys

repo_root = Path(__file__).resolve().parents[2]
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))

import numpy as np
import scipy.io as sio
from datetime import timedelta
from time import time
from console_formatting import print_files_saved, print_session
import project_paths as pp

import rec_list
pathHPCLC = rec_list.pathHPCLCopt
pathHPCLCterm = rec_list.pathHPCLCtermopt
paths = pathHPCLC + pathHPCLCterm
HPC_ALL_SESSIONS_STEM = pp.HPC_EPHYS_STEM / 'all_sessions'


#%% analysis
all_cells = {}

for pathname in paths:
    recname = Path(pathname).name
    print_session(recname)
    t0 = time()
    # load NeuronQuality file
    nq_file = sio.loadmat(Path(pathname) / f'{recname}.NeuronQuality.mat')

    # assume 6 shanks and load spikes into dict
    clu_count = 2
    for shank in range(1, 7):
        try:
            avspk = nq_file[f'nqShank{shank}']['AvSpk'][0][0]
        except KeyError:
            print(f'this recording does not have shank {shank}')
            continue
        tot_clu = avspk.shape[0]

        for clu in range(tot_clu):
            cluname = f'{recname} clu{clu_count} {shank} int{clu+2}'
            clu_count += 1
            all_cells[cluname] = avspk[clu, :]

    print('done; saving...')
    sess_folder = HPC_ALL_SESSIONS_STEM / recname
    sess_folder.mkdir(exist_ok=True)
    output_path = sess_folder / f'{recname}_all_waveforms.npy'
    np.save(output_path, all_cells)
    print_files_saved([
        ('waveforms', output_path),
    ])
    print(f'elapsed = {str(timedelta(seconds=int(time() - t0)))}')
