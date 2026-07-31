# -*- coding: utf-8 -*-
'''
Created on Sat Apr 11 2026

functions for HPC Bayesian decoding analyses.

@author: Dinghao Luo
'''

#%% imports
import numpy as np


#%% functions
def skipping_average(values, skip_step=250):
    '''
    average a 1d vector in consecutive fixed-width bins.

    the default 250 samples cover 200 ms at 1,250 Hz.
    '''
    n_int = np.floor_divide(len(values), skip_step)
    reduced = np.zeros(n_int)
    for i in range(n_int):
        reduced[i] = np.nanmean(values[i * skip_step:(i + 1) * skip_step])
    return reduced


def shuffle_mean(train, n_shuf=100, roll_bins=25):
    '''
    compute the mean of circularly shifted downsampled trains.
    '''
    shift = np.random.randint(1, roll_bins, n_shuf)
    shuf_array = np.zeros((n_shuf, roll_bins))

    for i in range(n_shuf):
        shuf_array[i, :] = np.roll(train, -shift[i])

    return np.mean(shuf_array, axis=0)
