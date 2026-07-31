# -*- coding: utf-8 -*-
'''
Created on Fri Mar  7 15:45:04 2025

calculate ISIs of single neurones

@author: Dinghao Luo


'''

#%% imports
import numpy as np
import sys
import os
from pathlib import Path

repo_root = Path(__file__).resolve().parents[2]
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))
from common_functions import mpl_formatting
from spike_text_io import param2array
import project_paths as pp
mpl_formatting()

import rec_list
paths = rec_list.pathLC


#%% analysis
for pathname in paths:
    recname = pathname[-17:]
    print('\n\nProcessing {}'.format(recname))

    sess_folder = pp.LC_EPHYS_STEM / 'all_sessions' / recname
    os.makedirs(sess_folder, exist_ok=True)

    clu = param2array(rf'{pathname}/{recname}.clu.1')  # load .clu
    res = param2array(rf'{pathname}/{recname}.res.1')  # load .res

    all_clu = [int(c) for c in np.unique(clu)
               if c not in ('', '0', '1') and
               len(np.where(clu==c)[0]) != 1]
    tot_clu = len(all_clu)

    spike_dict = {
        f'{recname} clu{clu_idx}': [
            int(t)
            for i, t in enumerate(res[:-1])  # last element is ''
            if int(clu[i])==clu_idx]
        for clu_idx in np.arange(2, tot_clu+2)
        }

    ISI_dict = {
        key: np.diff(spikes).astype(np.int32)
        for key, spikes in spike_dict.items()
        }

    np.save(sess_folder / f'{recname}_all_spikes.npy',
            spike_dict)
    np.save(sess_folder / f'{recname}_all_ISIs.npy',
            ISI_dict)
