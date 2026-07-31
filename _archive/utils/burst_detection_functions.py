# -*- coding: utf-8 -*-
'''
Created on 2 June 2026
functions extracted from LC_code/run_onset_burst_analysis/burst_detection.py,
    created on 13 March 2025
burst-detection functions for spike-train analyses

@author: Dinghao Luo
'''

#%% imports
import numpy as np


#%% burst-detection functions
def detect_bursts_from_isi(isis,
                           burst_thresh,
                           end_thresh,
                           min_spikes=3):
    '''
    detect bursts from inter-spike intervals.

    parameters:
    - isis: 1d array-like of inter-spike intervals in samples or seconds
    - burst_thresh: threshold below which a burst is initiated
    - end_thresh: threshold above which an active burst is terminated
    - min_spikes: minimum number of spikes required to keep a burst

    returns:
    - list of `(start_idx, end_idx)` tuples in isi-index coordinates
    '''
    bursts = []
    burst_start = None
    spike_count = 0

    for i, isi in enumerate(isis):
        if isi < burst_thresh:
            if burst_start is None:
                burst_start = i
                spike_count = 2
            else:
                spike_count += 1
        elif isi > end_thresh and burst_start is not None:
            if spike_count >= min_spikes:
                bursts.append((burst_start, i))
            burst_start = None
            spike_count = 0

    if burst_start is not None and spike_count >= min_spikes:
        bursts.append((burst_start, len(isis)))

    return bursts


def percentile_based_thresholds(isis,
                                burst_percentile=10,
                                end_factor=2.5,
                                samp_freq=20000,
                                max_thresh_seconds=0.2):
    '''
    compute burst onset and termination thresholds from the isi distribution.

    parameters:
    - isis: 1d array-like of inter-spike intervals in samples or seconds
    - burst_percentile: percentile used for the onset threshold
    - end_factor: multiplier applied to the onset threshold for burst termination
    - samp_freq: sampling frequency used when `isis` are expressed in samples
    - max_thresh_seconds: cap on the onset threshold, in seconds

    returns:
    - `(burst_thresh, end_thresh)`
    '''
    max_thresh = max_thresh_seconds * samp_freq
    burst_thresh = np.percentile(isis, burst_percentile)
    burst_thresh = min(burst_thresh, max_thresh)
    end_thresh = end_factor * burst_thresh
    return burst_thresh, end_thresh
