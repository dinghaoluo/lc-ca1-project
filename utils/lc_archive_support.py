'''
Created on Thu Jun 25 2026

load LC session dictionaries for archived raster plots

@author: Dinghao Luo
'''

#%% imports
from functools import lru_cache
import numpy as np

import project_paths as pp


#%% session dictionaries
@lru_cache(maxsize=4)
def load_lc_session_dict(recname, kind):
    path = (
        pp.LC_EPHYS_STEM
        / 'all_sessions'
        / recname
        / f'{recname}_{kind}.npy'
        )
    return np.load(path, allow_pickle=True).item()
