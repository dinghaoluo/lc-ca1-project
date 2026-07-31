# -*- coding: utf-8 -*-
'''
Created on Tue Jul  1 20:26:59 2025

plot example sess ref's and release map

@author: Dinghao Luo
'''

#%% imports
from pathlib import Path
import sys

import numpy as np
from PIL import Image

repo_root = Path(__file__).resolve().parents[2]
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))

import project_paths as pp


#%% main
recname = 'A126i-20250606-01'
processed_stem = pp.HPC_DLIGHT_LC_OPTO_STEM / 'all_sessions' / recname / 'processed_data'
output_stem = pp.FIGURES_ROOT / 'imaging' / 'figures_for_yingxue'
output_stem.mkdir(parents=True, exist_ok=True)

ref1 = np.load(processed_stem / 'ref_mat_ch1.npy', allow_pickle=True)
ref2 = np.load(processed_stem / 'ref_mat_ch2.npy', allow_pickle=True)
release_map = np.load(processed_stem / f'{recname}_release_map.npy', allow_pickle=True)

ref1_i = Image.fromarray(ref1)
ref1_i.save(output_stem / f'{recname}_ref1.tif')

ref2_i = Image.fromarray(ref2)
ref2_i.save(output_stem / f'{recname}_ref2.tif')

release_map_i = Image.fromarray(release_map)
release_map_i.save(output_stem / f'{recname}_release_map.tif')
