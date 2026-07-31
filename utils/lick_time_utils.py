# -*- coding: utf-8 -*-
'''
Created on 1 June 2026

lick-rate maps for the pharmacology summaries

@author: Dinghao Luo
'''

import numpy as np

from common_functions import smooth_convolve

LICK_TIME_MS = 5000

def lick_time_map(lick_times_aligned, duration_ms=LICK_TIME_MS, sigma=10):
    '''
    Convert trial-wise lick times, aligned to run onset in ms, into a lick-rate map.
    '''
    lick_map = np.zeros((len(lick_times_aligned), duration_ms), dtype=float)
    for trial_idx, trial in enumerate(lick_times_aligned):
        for lick in trial:
            lick_ms = float(lick)
            if np.isnan(lick_ms):
                continue
            if 0 <= lick_ms < duration_ms:
                lick_map[trial_idx, int(lick_ms)] += 1000.0
        lick_map[trial_idx] = smooth_convolve(
            lick_map[trial_idx],
            sigma=sigma
        )

    return lick_map
