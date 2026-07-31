# -*- coding: utf-8 -*-
'''
Created on Sat Jul  5 09:18:43 2025
Originally named plot_16_bit_maps.py

save 16-bit maps of ref channels of example session

@author: Dinghao Luo
'''

#%% imports
from pathlib import Path
import sys

import numpy as np
import tifffile

repo_root = Path(__file__).resolve().parents[2]
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))

import project_paths as pp


#%% main
path = pp.TWO_PHOTON_ROOT / 'A170i' / 'A170i-20260116' / 'A170i-20260116-03'

binpath = path / 'suite2p' / 'plane0' / 'data.bin'
bin2path = path / 'suite2p' / 'plane0' / 'data_chan2.bin'
opspath = path / 'suite2p' / 'plane0' / 'ops.npy'

ops = np.load(opspath, allow_pickle=True).item()
tot_frames = ops['nframes']
shape = tot_frames, ops['Ly'], ops['Lx']

print('loading movies and saving references...')
mov = np.memmap(binpath, mode='r', dtype='int16', shape=(1500, shape[1], shape[2])).astype(np.float32)
mov2 = np.memmap(bin2path, mode='r', dtype='int16', shape=(1500, shape[1], shape[2])).astype(np.float32)

tot_frames = mov.shape[0]

ref1 = np.mean(mov, axis=0)
ref2 = np.mean(mov2, axis=0)

ref1_16 = ref1.astype(np.uint16)
ref2_16 = ref2.astype(np.uint16)

output_stem = pp.IMAGING_REF_CHANNEL_MAPS_FIGURES_STEM
output_stem.mkdir(parents=True, exist_ok=True)
tifffile.imwrite(output_stem / 'nLight_ref1.tiff', ref1_16)
tifffile.imwrite(output_stem / 'nLight_ref2.tiff', ref2_16)
