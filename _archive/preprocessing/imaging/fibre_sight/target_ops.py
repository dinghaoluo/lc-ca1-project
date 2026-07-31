'''
Created on 7 April 2026
Archived on 23 July 2026
build foreground and support targets for training

@author: Dinghao Luo
'''

#%% imports
import numpy as np
from scipy import ndimage as ndi


#%% targets
def make_target(mask, mode='foreground', support_radius=3):
    mask = np.asarray(mask).astype(bool)

    if mode == 'foreground':
        return mask[None, ...].astype(np.float32)

    if mode == 'foreground_support':
        support = ndi.binary_dilation(mask, iterations=int(support_radius))
        target = np.stack([mask, support], axis=0)
        return target.astype(np.float32)

    raise ValueError(f'unknown target mode: {mode}')
