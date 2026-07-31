'''
Created on 5 April 2026
convert between xpix/ypix ROI dictionaries and labelled images

@author: Dinghao Luo
'''

#%% imports
from pathlib import Path

import numpy as np

from ._repo import add_repo_paths

add_repo_paths()

from imaging_utility_functions import build_roi_mask


#%% loading
def load_roi_dict(path):
    roi_dict = np.load(path, allow_pickle=True).item()
    if not isinstance(roi_dict, dict):
        raise ValueError(f'not a ROI dict: {path}')
    return roi_dict


def save_roi_dict(roi_dict, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    np.save(path, roi_dict)


#%% conversion
# xpix/ypix is the format used by the older imaging scripts, so labels made
# here can return to the rest of the pipeline without a private model format
def clean_roi_dict(roi_dict, shape):
    cleaned = {}
    invalid_entries = 0

    for roi_id, roi in roi_dict.items():
        try:
            xpix = np.asarray(roi['xpix'], dtype=np.int64).ravel()
            ypix = np.asarray(roi['ypix'], dtype=np.int64).ravel()
        except Exception:
            invalid_entries += 1
            continue

        same_length = len(xpix) == len(ypix)
        if not same_length:
            invalid_entries += 1
            n_pix = min(len(xpix), len(ypix))
            xpix = xpix[:n_pix]
            ypix = ypix[:n_pix]

        in_bounds = (
            (xpix >= 0) &
            (xpix < shape[1]) &
            (ypix >= 0) &
            (ypix < shape[0])
            )
        if not np.all(in_bounds):
            invalid_entries += 1

        xpix = xpix[in_bounds]
        ypix = ypix[in_bounds]
        if len(xpix) == 0:
            continue

        cleaned[roi_id] = {'xpix': xpix, 'ypix': ypix}

    return cleaned, invalid_entries


def roi_dict_to_mask(roi_dict, shape):
    cleaned, invalid_entries = clean_roi_dict(roi_dict, shape)
    mask = build_roi_mask(cleaned, shape=shape)
    return mask.astype(bool), invalid_entries


def roi_dict_to_label(roi_dict, shape):
    cleaned, invalid_entries = clean_roi_dict(roi_dict, shape)
    labelled = np.zeros(shape, dtype=np.int32)
    roi_areas = []

    for new_id, roi in enumerate(cleaned.values(), start=1):
        labelled[roi['ypix'], roi['xpix']] = new_id
        roi_areas.append(len(roi['xpix']))

    return labelled, roi_areas, invalid_entries


def labels_to_roi_dict(labelled):
    labelled = np.asarray(labelled)
    roi_dict = {}

    next_id = 1
    for label_id in sorted(np.unique(labelled)):
        if label_id == 0:
            continue
        ypix, xpix = np.where(labelled == label_id)
        if len(xpix) == 0:
            continue
        roi_dict[next_id] = {
            'xpix': xpix.astype(np.int64),
            'ypix': ypix.astype(np.int64),
            }
        next_id += 1

    return roi_dict


def roi_summary(roi_dict, shape):
    labelled, roi_areas, invalid_entries = roi_dict_to_label(roi_dict, shape)
    mask = labelled > 0

    return {
        'roi_count': len(roi_areas),
        'positive_pixels': int(np.sum(mask)),
        'positive_fraction': float(np.mean(mask)),
        'median_roi_area_px': float(np.median(roi_areas)) if roi_areas else 0.0,
        'min_roi_area_px': int(np.min(roi_areas)) if roi_areas else 0,
        'max_roi_area_px': int(np.max(roi_areas)) if roi_areas else 0,
        'invalid_roi_entries': int(invalid_entries),
        }
