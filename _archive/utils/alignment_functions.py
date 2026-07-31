# -*- coding: utf-8 -*-
'''
Created on Wed Dec 18 17:19:22 2024

functions for alignment testing

@author: Dinghao Luo
'''


#%% definitions
def bootstrap_ratio(spike_arr,
                    bootstrap=500,
                    samp_freq=1250,
                    length=6*1250,
                    GPU_AVAILABLE=False,
                    VERBOSE=True):
    '''
    bootstrap post/pre ratios after independent circular trial shifts

    return percentile thresholds ordered [99.9, 99, 95, 50, 5, 1, 0.1]
    '''
    from tqdm import tqdm

    if GPU_AVAILABLE:
        import cupy as xp
        device = 'GPU'
    else:
        import numpy as xp
        device = 'CPU'

    tot_trials = spike_arr.shape[0]
    shuf_ratio = xp.zeros(bootstrap)

    indices = xp.arange(length)

    iterator = tqdm(
        range(bootstrap),
        desc=f'lick sensitivity ({device})',
        disable=not VERBOSE,
        )
    for shuf in iterator:
        rand_shifts = xp.random.randint(1, length, tot_trials)
        shifted_indices = (indices[None, :] - rand_shifts[:, None]) % length
        shuf_arr = spike_arr[xp.arange(tot_trials)[:, None], shifted_indices]
        shuf_result = xp.mean(shuf_arr, axis=0)

        shuf_ratio[shuf] = (
            xp.sum(shuf_result[length // 2:length // 2 + samp_freq]) /
            xp.sum(shuf_result[length // 2 - samp_freq:length // 2])
            )

    return xp.percentile(shuf_ratio, [99.9, 99, 95, 50, 5, 1, .1], axis=0)
