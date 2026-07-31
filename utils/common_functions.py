# -*- coding: utf-8 -*-
'''
Created on Thu Aug  4 14:54:04 2022

general functions for LC-CA1 scripts.

@author: Dinghao Luo
'''
from __future__ import annotations

#%% imports

import numpy as np


#%% colour constants
colour_tagged = (70 / 255, 101 / 255, 175 / 255)
colour_putative = (101 / 255, 82 / 255, 163 / 255)
colour_other = (169 / 255, 169 / 255, 169 / 255)


#%% gpu helpers
def get_GPU_availability():
    '''
    try to import CuPy and detect whether a cuda device is available.

    returns
    -------
    tuple
        `(cp_module_or_none, gpu_available, device_name_or_none)`
    '''
    try:
        import cupy as cp
    except ImportError:
        return None, False, None

    try:
        if cp.cuda.runtime.getDeviceCount() <= 0:
            print('GPU status: not available')
            return None, False, None

        cp.cuda.set_allocator(cp.cuda.MemoryPool().malloc)
        cp.cuda.set_pinned_memory_allocator(cp.cuda.PinnedMemoryPool().malloc)
        dev_id = cp.cuda.runtime.getDevice()
        props = cp.cuda.runtime.getDeviceProperties(dev_id)
        device_name = props['name'].decode()
        print(f'GPU status: available ({device_name})')
        return cp, True, device_name
    except RuntimeError:
        print('GPU status: not available')
        return None, False, None


#%% plotting helpers
def mpl_formatting():
    '''
    apply a minimal matplotlib style used across the project.
    '''
    import matplotlib

    matplotlib.rcParams.update(
        {
            'font.family': 'Arial',
            'pdf.fonttype': 42,
            'ps.fonttype': 42,
            'axes.labelsize': 10,
            'xtick.labelsize': 9,
            'ytick.labelsize': 9,
        }
    )


#%% normalisation helpers
def normalise(data, axis=1):
    '''
    normalise an array along a specified axis while ignoring NaN values.
    '''
    array = np.asarray(data, dtype=float)
    if array.ndim > 1:
        return np.apply_along_axis(
            lambda x: (x - np.nanmin(x)) / (np.nanmax(x) - np.nanmin(x)),
            axis=axis,
            arr=array,
        )
    return (array - np.nanmin(array)) / (np.nanmax(array) - np.nanmin(array))

def normalise_to_all(data, alldata):
    '''
    normalise a 1-d vector relative to the global min and max of another array.
    '''
    data = np.asarray(data, dtype=float)
    alldata = np.asarray(alldata, dtype=float)
    amin = np.nanmin(alldata)
    amax = np.nanmax(alldata)
    return (data - amin) / (amax - amin)

def get_trialtype_indices_from_stim_conds(stim_conds, control_offset=2):
    '''
    derive baseline, stimulation, and matched-control trial indices from stimulation labels.

    parameters:
    - stim_conds: iterable of stimulation-condition labels
    - control_offset: trial offset used to define matched control trials

    returns:
    - baseline_idx: list of baseline trial indices before the first stimulation trial
    - stim_idx: list of stimulation trial indices
    - ctrl_idx: list of valid control trial indices
    '''
    stim_idx = [trial for trial, cond in enumerate(stim_conds) if cond != '0']

    if not stim_idx:
        return list(range(len(stim_conds))), [], []

    ctrl_idx = [idx + control_offset for idx in stim_idx]
    return list(range(stim_idx[0])), stim_idx, ctrl_idx


#%% statistics helpers

def smooth_convolve(data, sigma=3, axis=1):
    '''
    apply gaussian smoothing to a 1-d or n-d array using convolution.
    '''
    data = np.asarray(data, dtype=float)
    kernel = gaussian_kernel_unity(sigma)
    pad_width = len(kernel) // 2

    if data.ndim == 1:
        data_padded = np.pad(data, pad_width, mode='reflect')
        return np.convolve(data_padded, kernel, mode='same')[pad_width:-pad_width]

    pad_config = [(0, 0)] * data.ndim
    pad_config[axis] = (pad_width, pad_width)
    data_padded = np.pad(data, pad_config, mode='reflect')
    smoothed = np.apply_along_axis(
        lambda x: np.convolve(x, kernel, mode='same'),
        axis=axis,
        arr=data_padded,
    )
    slice_config = [slice(None)] * data.ndim
    slice_config[axis] = slice(pad_width, -pad_width)
    return smoothed[tuple(slice_config)]

def gaussian_kernel_unity(sigma, GPU_AVAILABLE=False):
    '''
    generate a gaussian kernel normalised to sum to one.
    '''
    kernel_size = int(6 * sigma + 1)
    x = np.arange(kernel_size) - (kernel_size // 2)
    kernel = np.exp(-(x**2 / (2 * sigma**2)))
    kernel /= kernel.sum()

    if GPU_AVAILABLE:
        import cupy as cp

        return cp.asarray(kernel)
    return kernel

def replace_outlier(arr, method='std', k=5):
    '''
    replace outliers with linearly interpolated values.
    '''
    if arr.ndim == 2:
        arr = arr.ravel()

    if method == 'std':
        mean = np.mean(arr)
        std = np.std(arr)
        outliers = np.abs(arr - mean) > k * std
    elif method == 'mad':
        median = np.median(arr)
        mad = np.median(np.abs(arr - median))
        outliers = np.abs(arr - median) > k * mad
    else:
        raise ValueError('invalid method; choose \'std\' or \'mad\'')

    if not np.any(outliers):
        return arr

    from scipy.interpolate import interp1d

    indices = np.arange(len(arr))
    valid_indices = indices[~outliers]
    valid_values = arr[~outliers]
    interp_func = interp1d(valid_indices, valid_values, kind='linear', fill_value='extrapolate')
    arr[outliers] = interp_func(indices[outliers])
    return arr
