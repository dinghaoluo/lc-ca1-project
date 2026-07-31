'''
Created on 13 May 2026
Modified on 23 June 2026
MSER segmentation adapted from the hand-curated FibreSegger workflow of 11 April 2025

@author: Dinghao Luo
'''


#%% imports
import cv2
import numpy as np
from scipy.ndimage import median_filter
from skimage.measure import regionprops


#%% parameters
PARAMETER_SPECS = [
    {
        'name': 'solidity min',
        'kind': 'float',
        'default': 0.1,
        'minimum': 0.0,
        'maximum': 1.0,
        'step': 0.01,
        'decimals': 2,
    },
    {
        'name': 'eccentricity min',
        'kind': 'float',
        'default': 0.75,
        'minimum': 0.0,
        'maximum': 1.0,
        'step': 0.01,
        'decimals': 2,
    },
    {
        'name': 'thinness max',
        'kind': 'float',
        'default': 0.8,
        'minimum': 0.0,
        'maximum': 1.0,
        'step': 0.01,
        'decimals': 2,
    },
    {
        'name': 'tophat kernel',
        'kind': 'int',
        'default': 11,
        'minimum': 1,
        'maximum': 99,
        'step': 2,
    },
    {
        'name': 'clahe clip',
        'kind': 'float',
        'default': 2.0,
        'minimum': 0.1,
        'maximum': 10.0,
        'step': 0.1,
        'decimals': 2,
    },
    {
        'name': 'MSER max variation',
        'kind': 'float',
        'default': 1.2,
        'minimum': 0.0,
        'maximum': 10.0,
        'step': 0.05,
        'decimals': 2,
    },
    {
        'name': 'MSER delta',
        'kind': 'int',
        'default': 5,
        'minimum': 1,
        'maximum': 50,
        'step': 1,
    },
    {
        'name': 'MSER min area',
        'kind': 'int',
        'default': 30,
        'minimum': 1,
        'maximum': 5000,
        'step': 1,
    },
    {
        'name': 'MSER max area',
        'kind': 'int',
        'default': 15000,
        'minimum': 10,
        'maximum': 100000,
        'step': 100,
    },
    {
        'name': 'aspect ratio min',
        'kind': 'float',
        'default': 1.2,
        'minimum': 1.0,
        'maximum': 10.0,
        'step': 0.05,
        'decimals': 2,
    },
    {
        'name': 'clip-percentile',
        'kind': 'float',
        'default': 99.0,
        'minimum': 0.0,
        'maximum': 100.0,
        'step': 0.1,
        'decimals': 2,
    },
    {
        'name': 'area min',
        'kind': 'int',
        'default': 30,
        'minimum': 1,
        'maximum': 5000,
        'step': 1,
    },
    {
        'name': 'MSER threshold',
        'kind': 'float',
        'default': 85.0,
        'minimum': 0.0,
        'maximum': 100.0,
        'step': 0.1,
        'decimals': 2,
    },
]

PARAMETER_TOOLTIPS = {
    'solidity min': 'minimum filled-area ratio accepted for candidate ROIs',
    'eccentricity min': 'minimum elongation; higher values favour fibre-like shapes',
    'thinness max': 'upper compactness filter; lower values favour thin objects',
    'tophat kernel': 'background-removal kernel size before MSER detection',
    'clahe clip': 'local contrast equalisation strength',
    'MSER max variation': 'MSER stability filter; lower values are stricter',
    'MSER delta': 'MSER intensity step size',
    'MSER min area': 'smallest MSER component area kept, in pixels',
    'MSER max area': 'largest MSER component area kept, in pixels',
    'aspect ratio min': 'minimum long-axis to short-axis ratio',
    'clip-percentile': 'upper intensity percentile used before normalising',
    'area min': 'smallest final ROI area kept, in pixels',
    'MSER threshold': 'brightness cutoff used before MSER candidate selection',
}


def default_segment_params():
    return {spec['name']: spec['default'] for spec in PARAMETER_SPECS}


#%% image preparation
def enhance_contrast_u8(image, tophat_kernel=11, clahe_clip=2.0):
    image = np.asarray(image, dtype=np.float32)
    low, high = np.percentile(image, [1, 99])
    if high <= low:
        high = low + 1.0
    image01 = np.clip((image - low) / (high - low), 0, 1)

    kernel_size = int(tophat_kernel)
    if kernel_size % 2 == 0:
        kernel_size += 1

    kernel = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (kernel_size, kernel_size))
    tophat = cv2.morphologyEx((image01 * 255).astype(np.uint8), cv2.MORPH_TOPHAT, kernel)

    clahe = cv2.createCLAHE(clipLimit=float(clahe_clip), tileGridSize=(8, 8))
    return clahe.apply(tophat)


#%% segmentation
def run_mser_segmentation(image, params, fixed_rois=None):
    image = np.asarray(image)
    if image.ndim != 2:
        raise ValueError(f'expected a 2D image, got shape {image.shape}')

    fixed_rois = _normalise_fixed_rois(fixed_rois, image.shape)
    filtered_image = median_filter(image, size=(3, 3))
    image_u8 = enhance_contrast_u8(
        filtered_image,
        tophat_kernel=params['tophat kernel'],
        clahe_clip=params['clahe clip'],
    )

    # keep the cutoff relative to each channel-2 image
    threshold = np.percentile(image_u8, params['MSER threshold'])
    thresholded = image_u8.copy()
    thresholded[thresholded < threshold] = 0

    mser = _make_mser(params)
    regions, _ = mser.detectRegions(thresholded)
    candidate_label = _regions_to_label(regions, image_u8.shape)
    candidate_rois = _filter_candidates(candidate_label, params)

    roi_dict, labelled, fixed_ids = _combine_rois(
        image.shape,
        fixed_rois=fixed_rois,
        candidate_rois=candidate_rois,
    )
    stats = {
        'MSER regions': len(regions),
        'candidate ROIs': len(candidate_rois),
        'kept ROIs': len(roi_dict),
        'fixed ROIs': len(fixed_ids),
    }
    return roi_dict, labelled, fixed_ids, stats


def _make_mser(params):
    delta = int(max(1, round(params['MSER delta'])))
    min_area = int(max(5, round(params['MSER min area'])))
    max_area = int(max(min_area + 1, round(params['MSER max area'])))

    # lab workstations still span both OpenCV constructor signatures
    try:
        mser = cv2.MSER_create(delta, min_area, max_area)
    except TypeError:
        mser = cv2.MSER_create(_delta=delta, _min_area=min_area, _max_area=max_area)
    mser.setMaxVariation(float(params['MSER max variation']))
    return mser


def _regions_to_label(regions, shape):
    labelled = np.zeros(shape, dtype=np.int32)
    for roi_id, region in enumerate(regions, start=1):
        labelled[region[:, 1], region[:, 0]] = roi_id
    return labelled


def _filter_candidates(labelled, params):
    roi_list = []
    for region in regionprops(labelled):
        if region.area < params['area min']:
            continue

        eccentricity = region.eccentricity if np.isfinite(region.eccentricity) else 0.0
        solidity = region.solidity if np.isfinite(region.solidity) else 0.0
        if hasattr(region, 'axis_minor_length'):
            minor_length = region.axis_minor_length
            major_length = region.axis_major_length
        else:
            minor_length = region.minor_axis_length
            major_length = region.major_axis_length
        minor_axis = minor_length if minor_length > 1e-6 else 1e-6
        aspect_ratio = major_length / minor_axis
        perimeter = region.perimeter if region.perimeter > 1e-6 else 1e-6
        thinness = 4 * np.pi * region.area / (perimeter ** 2)

        if solidity < params['solidity min']:
            continue
        if eccentricity < params['eccentricity min']:
            continue
        if aspect_ratio < params['aspect ratio min']:
            continue
        if thinness > params['thinness max']:
            continue

        ypix, xpix = region.coords[:, 0], region.coords[:, 1]
        roi_list.append({'xpix': xpix, 'ypix': ypix})
    return roi_list


def _normalise_fixed_rois(fixed_rois, shape):
    if not fixed_rois:
        return []

    if isinstance(fixed_rois, dict):
        roi_iterable = fixed_rois.values()
    else:
        roi_iterable = fixed_rois

    cleaned = []
    for roi in roi_iterable:
        xpix = np.asarray(roi.get('xpix', []), dtype=np.int64).ravel()
        ypix = np.asarray(roi.get('ypix', []), dtype=np.int64).ravel()
        n_pix = min(len(xpix), len(ypix))
        xpix = xpix[:n_pix]
        ypix = ypix[:n_pix]
        in_bounds = (
            (xpix >= 0) &
            (xpix < shape[1]) &
            (ypix >= 0) &
            (ypix < shape[0])
        )
        xpix = xpix[in_bounds]
        ypix = ypix[in_bounds]
        if len(xpix) > 0:
            cleaned.append({'xpix': xpix, 'ypix': ypix})
    return cleaned


def _combine_rois(shape, fixed_rois, candidate_rois):
    labelled = np.zeros(shape, dtype=np.int32)
    roi_dict = {}
    fixed_ids = set()
    next_id = 1

    # fixed ROIs go in first so hand curation wins over a fresh MSER pass
    for roi in fixed_rois:
        xpix, ypix = roi['xpix'], roi['ypix']
        roi_dict[next_id] = {'xpix': xpix, 'ypix': ypix}
        labelled[ypix, xpix] = next_id
        fixed_ids.add(next_id)
        next_id += 1

    for roi in candidate_rois:
        xpix, ypix = roi['xpix'], roi['ypix']
        if np.any(labelled[ypix, xpix] > 0):
            continue
        roi_dict[next_id] = {'xpix': xpix, 'ypix': ypix}
        labelled[ypix, xpix] = next_id
        next_id += 1

    return roi_dict, labelled, fixed_ids
