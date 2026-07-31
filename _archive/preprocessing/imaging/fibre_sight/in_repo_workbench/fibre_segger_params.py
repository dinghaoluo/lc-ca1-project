'''
Created on 12 May 2026
starting values for the MSER controls carried over from FibreSegger

@author: Dinghao Luo
'''

#%% imports
import numpy as np


#%% defaults
DEFAULT_FIBRE_SEGGER_PARAMS = {
    'solidity min': 0.1,
    'eccentricity min': 0.75,
    'thinness max': 0.8,
    'tophat kernel': 11,
    'clahe clip': 2.0,
    'MSER max variation': 1.2,
    'MSER delta': 5,
    'MSER min area': 30,
    'MSER max area': 15000,
    'aspect ratio min': 1.2,
    'clip-percentile': 99.0,
    'area min': 30,
    'MSER threshold': 85.0,
    }


#%% recommendations
def recommend_fibre_segger_params(image, roi_dict=None):
    # I still tune these controls against each image; this only avoids starting
    # from the same values when the contrast or hand-labelled ROI sizes differ
    image = np.asarray(image)
    params = dict(DEFAULT_FIBRE_SEGGER_PARAMS)
    contrast = robust_contrast(image)
    bright_fraction = get_bright_fraction(image)

    if contrast < 0.12:
        params['clahe clip'] = 3.0
        params['clip-percentile'] = 99.5
        params['MSER threshold'] = 80.0
    elif bright_fraction > 0.12:
        params['MSER threshold'] = 90.0
        params['MSER max variation'] = 0.9

    if roi_dict:
        areas = [len(roi['xpix']) for roi in roi_dict.values() if 'xpix' in roi]
        if areas:
            median_area = float(np.median(areas))
            params['MSER min area'] = int(max(10, round(median_area * 0.2)))
            params['area min'] = int(max(10, round(median_area * 0.2)))
            params['MSER max area'] = int(max(1000, round(np.percentile(areas, 95) * 3)))

    return params


def robust_contrast(image):
    low, high = np.nanpercentile(image, [5, 99])
    denom = max(abs(high), 1)
    return float((high - low) / denom)


def get_bright_fraction(image, percentile=95):
    threshold = np.nanpercentile(image, percentile)
    return float(np.mean(image >= threshold))
