'''
Created on 6 April 2026
Modified on 29 June 2026
scan labelled sessions and keep the train, validation, and test split

@author: Dinghao Luo
'''

#%% imports
from collections import defaultdict
from pathlib import Path
import csv

import numpy as np

from .roi_io import load_roi_dict, roi_summary


#%% columns
MANIFEST_COLUMNS = [
    'session',
    'animal',
    'processed_path',
    'image_path',
    'roi_path',
    'height',
    'width',
    'roi_count',
    'positive_pixels',
    'positive_fraction',
    'median_roi_area_px',
    'min_roi_area_px',
    'max_roi_area_px',
    'invalid_roi_entries',
    'included',
    'exclusion_reason',
    'split',
    ]


#%% scanning
def find_ref_paths(processed_path):
    # newer sessions use the canonical name; older exports kept the recording prefix
    canonical_paths = sorted(processed_path.glob('ref_mat_ch2.npy'))
    archive_paths = sorted(processed_path.glob('*_ref_mat_ch2.npy'))
    return canonical_paths + [
        path for path in archive_paths
        if path not in canonical_paths
        ]


def find_ref_roi_pair(processed_path):
    ref_paths = find_ref_paths(processed_path)
    roi_paths = sorted(processed_path.glob('*_ROI_dict.npy'))

    if not ref_paths or not roi_paths:
        return None, None

    ref_path = ref_paths[0]
    if ref_path.name == 'ref_mat_ch2.npy':
        prefix = processed_path.parent.name
    else:
        prefix = ref_path.name.replace('_ref_mat_ch2.npy', '')
    matched_rois = [path for path in roi_paths if path.name.startswith(prefix)]
    roi_path = matched_rois[0] if matched_rois else roi_paths[0]
    return ref_path, roi_path


def get_exclusion_reason(session_path):
    processed_path = session_path / 'processed_data'
    if not processed_path.exists():
        return 'no processed_data'

    ref_paths = find_ref_paths(processed_path)
    roi_paths = sorted(processed_path.glob('*_ROI_dict.npy'))
    has_1100 = bool(list(processed_path.glob('*1100*.npy')))

    if not ref_paths:
        return 'no channel-2 reference'
    if not roi_paths:
        return 'no ROI dict'
    if has_1100:
        return 'contains 1100 reference'
    return ''


def scan_session(session_path):
    session_path = Path(session_path)
    processed_path = session_path / 'processed_data'
    exclusion_reason = get_exclusion_reason(session_path)

    row = {
        'session': session_path.name,
        'animal': session_path.name.split('-')[0],
        'processed_path': str(processed_path),
        'image_path': '',
        'roi_path': '',
        'height': '',
        'width': '',
        'roi_count': 0,
        'positive_pixels': 0,
        'positive_fraction': 0,
        'median_roi_area_px': 0,
        'min_roi_area_px': 0,
        'max_roi_area_px': 0,
        'invalid_roi_entries': 0,
        'included': exclusion_reason == '',
        'exclusion_reason': exclusion_reason,
        'split': '',
        }

    if exclusion_reason:
        return row

    image_path, roi_path = find_ref_roi_pair(processed_path)
    image = np.load(image_path, mmap_mode='r')
    roi_dict = load_roi_dict(roi_path)
    summary = roi_summary(roi_dict, image.shape)

    row.update(summary)
    row.update({
        'image_path': str(image_path),
        'roi_path': str(roi_path),
        'height': int(image.shape[0]),
        'width': int(image.shape[1]),
        })
    return row


def scan_source_root(source_root):
    source_root = Path(source_root)
    rows = []

    for session_path in sorted(path for path in source_root.iterdir() if path.is_dir()):
        rows.append(scan_session(session_path))

    return rows


#%% splitting
def assign_session_splits(rows, val_fraction=0.15, test_fraction=0.15, seed=7):
    # Split whole sessions, never patches, so one recording cannot enter training and held-out data.
    rng = np.random.default_rng(seed)
    included_rows = [row for row in rows if as_bool(row['included'])]

    # make each animal contribute to the held-out splits when it has enough sessions
    grouped = defaultdict(list)
    for row in included_rows:
        grouped[row['animal']].append(row)

    for animal_rows in grouped.values():
        order = np.arange(len(animal_rows))
        rng.shuffle(order)

        n_rows = len(animal_rows)
        n_test = _split_count(n_rows, test_fraction)
        n_val = _split_count(n_rows - n_test, val_fraction)

        for rank, row_idx in enumerate(order):
            row = animal_rows[int(row_idx)]
            if rank < n_test:
                row['split'] = 'test'
            elif rank < n_test + n_val:
                row['split'] = 'val'
            else:
                row['split'] = 'train'

    for row in rows:
        if not as_bool(row['included']):
            row['split'] = ''

    return rows


def _split_count(n_rows, fraction):
    if fraction <= 0 or n_rows <= 1:
        return 0
    return max(1, int(round(n_rows * fraction)))


#%% io
def read_manifest(path, included_only=False, split=None):
    rows = []
    with open(path, 'r', newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            row = coerce_row(row)
            if included_only and not as_bool(row['included']):
                continue
            if split is not None and row.get('split', '') != split:
                continue
            rows.append(row)
    return rows


def write_manifest(rows, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)

    with open(path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=MANIFEST_COLUMNS)
        writer.writeheader()
        for row in rows:
            writer.writerow({col: row.get(col, '') for col in MANIFEST_COLUMNS})


def coerce_row(row):
    row = dict(row)

    for col in ['height', 'width', 'roi_count', 'positive_pixels',
                'min_roi_area_px', 'max_roi_area_px', 'invalid_roi_entries']:
        row[col] = int(float(row[col])) if row.get(col, '') != '' else 0

    for col in ['positive_fraction', 'median_roi_area_px']:
        row[col] = float(row[col]) if row.get(col, '') != '' else 0.0

    row['included'] = as_bool(row.get('included', False))
    return row


def as_bool(value):
    if isinstance(value, bool):
        return value
    return str(value).lower() in {'true', '1', 'yes', 'y'}
