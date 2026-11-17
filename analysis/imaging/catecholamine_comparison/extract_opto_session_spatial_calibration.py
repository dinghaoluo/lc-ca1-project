# -*- coding: utf-8 -*-
'''
Created on Sun Aug 30 01:33:41 2026

extract pixel geometry from the optogenetic ScanImage headers

@author: Dinghao Luo
'''

#%% imports
import sys
from pathlib import Path

repo_root = Path(__file__).resolve().parents[3]
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))

import project_paths as pp
import rec_list

import numpy as np
import pandas as pd
import tifffile


#%% paths
repo_root       = Path(__file__).resolve().parents[3]
COMPARISON_ROOT = repo_root / 'data' / 'catecholamine_comparison'
RAW_ROOT        = pp.TWO_PHOTON_ROOT
SPATIAL_INPUT   = COMPARISON_ROOT / 'spatial_tau_metrics.csv'
OUTPUT_PATH     = COMPARISON_ROOT / 'opto_session_spatial_calibration.csv'

EXPECTED_SESSIONS        = 47
EXPECTED_COMMON_GEOMETRY = 45
EXPECTED_EXCLUSIONS      = {
    'A152i-20251113-01',
    'A152i-20251114-01',
    }


#%% header extraction
def raw_recording_path(session):
    animal, date, _ = session.split('-')
    return RAW_ROOT / animal / f'{animal}-{date}' / session


def session_geometry(sensor, animal, session):
    recording_path = raw_recording_path(session)
    tiffs = sorted(recording_path.glob('*.tif'))
    if not tiffs:
        raise FileNotFoundError(f'no raw TIFF found for {session}')

    raw_tiff = tiffs[0]
    with tifffile.TiffFile(raw_tiff) as stack:
        metadata = stack.scanimage_metadata
        if metadata is None:
            raise ValueError(f'no ScanImage metadata found in {raw_tiff}')
        frame_data = metadata['FrameData']
        image_height_px, image_width_px = stack.pages[0].shape
        zoom = float(frame_data['SI.hRoiManager.scanZoomFactor'])

    common_zoom3 = (
        np.isclose(zoom, 3)
        and image_height_px == 512
        and image_width_px == 512
        )
    return {
        'sensor': sensor,
        'animal': animal,
        'session': session,
        'raw_tiff': str(raw_tiff),
        'raw_tiff_size_bytes': raw_tiff.stat().st_size,
        'image_height_px': int(image_height_px),
        'image_width_px': int(image_width_px),
        'scan_zoom_factor': zoom,
        'metadata_key': 'SI.hRoiManager.scanZoomFactor',
        'common_zoom3_geometry': common_zoom3,
        'chapter3_status': (
            'pixel-domain calibration diagnostic; omitted from Chapter 3'
            ),
        }


#%% run
def main():
    source = pd.read_csv(SPATIAL_INPUT)
    sessions = source[['sensor', 'animal', 'recname']].drop_duplicates()
    if sessions['recname'].duplicated().any():
        raise ValueError('a session appears under more than one sensor or animal')

    rows = [
        session_geometry(row.sensor, row.animal, row.recname)
        for row in sessions.itertuples(index=False)
        ]
    output = pd.DataFrame(rows).sort_values(
        ['sensor', 'animal', 'session'],
        ).reset_index(drop=True)
    if len(output) != source['recname'].nunique():
        raise AssertionError('not every spatial-analysis session received a header row')
    if len(output) != EXPECTED_SESSIONS:
        raise AssertionError(
            f'expected {EXPECTED_SESSIONS} session fields, found {len(output)}'
            )
    common = output['common_zoom3_geometry'].astype(bool)
    if common.sum() != EXPECTED_COMMON_GEOMETRY:
        raise AssertionError(
            f'expected {EXPECTED_COMMON_GEOMETRY} common-geometry fields, '
            f'found {common.sum()}'
            )
    exclusions = set(output.loc[~common, 'session'])
    if exclusions != EXPECTED_EXCLUSIONS:
        raise AssertionError(
            f'unexpected common-geometry exclusions: {sorted(exclusions)}'
            )

    output.to_csv(OUTPUT_PATH, index=False)
    print(f'wrote {len(output)} session geometries to {OUTPUT_PATH}')
    print(output.groupby(
        ['scan_zoom_factor', 'image_height_px', 'image_width_px'],
        ).size().to_string())
    print(
        '45/47 fields have common 512-by-512, zoom-three geometry; '
        'A152i-20251113-01 and A152i-20251114-01 are excluded. '
        'No physical conversion is made.'
        )


if __name__ == '__main__':
    main()
