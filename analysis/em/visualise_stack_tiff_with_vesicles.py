# -*- coding: utf-8 -*-
'''
Created on Mon May 25 17:43:15 2026

3D visualisation of EM axon profiles with vesicles and mitochondria

@author: Dinghao Luo
'''

#%% imports
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile
import napari
from scipy.ndimage import gaussian_filter, binary_dilation
from skimage.measure import marching_cubes


#%% data paths
DATA_ROOT_ORIG = Path(r'Z:\EM-DATA\Annotations of profiles'
                      r'\new exports from Dragonfly 6-16-2026')
DATA_ROOT_REVISED = Path(r'Z:\EM-DATA\Annotations of profiles'
                         r'\new revised profile 3 and 4')

PROFILES = {
    'profile 1': {
        'dir': DATA_ROOT_ORIG / 'profile 1',
        'axon_pattern': 'Profile 1 Axon Segmentation*.tiff',
        'mito_pattern': 'Profile 1 Mitochondria Segmentation*.tiff',
        'csv': 'Profile 1 Vesicle Positions.csv',
        'csv_sep': ',',
        'coords_in_nm': False,
    },
    'profile 2': {
        'dir': DATA_ROOT_ORIG / 'profile 2',
        'axon_pattern': 'Profile 2 Axon Segmentation*.tiff',
        'mito_pattern': 'Profile 2 Mitochondria Segmentation*.tiff',
        'csv': 'Profile 2 Vesicle Positions.csv',
        'csv_sep': ',',
        'coords_in_nm': False,
    },
    'profile 3': {
        'dir': DATA_ROOT_REVISED / 'new profile 3',
        'axon_pattern': 'Profile 3 Axon segmentation*.tiff',
        'mito_pattern': 'Profile 3 Mito segmentation*.tiff',
        'csv': 'Profile 3 Vesicle Positions.csv',
        'csv_sep': ';',
        'coords_in_nm': True,
    },
    'profile 4': {
        'dir': DATA_ROOT_REVISED / 'new profile 4',
        'axon_pattern': 'Profile 4 Axon Segmentation*.tiff',
        'mito_pattern': 'Profile 4 Mitochondria Segmentation*.tiff',
        'csv': 'Profile 4 Vesicle Positions.csv',
        'csv_sep': ';',
        'coords_in_nm': True,
    },
}

SCALE = (50, 5, 5)  # Z, Y, X in nm
SMOOTH = True
SMOOTH_SIGMA = (2, 1, 1)  # heavier in Z to interpolate across slices


#%% select profile
print('Available profiles:')
for i, name in enumerate(PROFILES, 1):
    print(f'  {i}. {name}')
print(f'  {len(PROFILES)+1}. all')

choice = input('Select profile: ').strip()
if choice == str(len(PROFILES)+1) or choice.lower() == 'all':
    selected = list(PROFILES.keys())
else:
    selected = [list(PROFILES.keys())[int(choice) - 1]]


#%% load and visualise
for profile_name in selected:
    paths = PROFILES[profile_name]
    profile_dir = paths['dir']

    axon_files = sorted(profile_dir.glob(paths['axon_pattern']))
    axon_vol = np.stack([tifffile.imread(f) for f in axon_files])
    axon_vol = (axon_vol > 0).astype(np.uint8)

    mito_files = sorted(profile_dir.glob(paths['mito_pattern']))
    mito_vol = np.stack([tifffile.imread(f) for f in mito_files])
    mito_vol = (mito_vol > 0).astype(np.uint8)

    df = pd.read_csv(profile_dir / paths['csv'], sep=paths['csv_sep'])
    scale = 1.0 if paths['coords_in_nm'] else 1000.0
    points = np.column_stack([
        df['PosZ'].values * scale / SCALE[0],
        df['PosY'].values * scale / SCALE[1],
        df['PosX'].values * scale / SCALE[2],
    ])

    vol_for_mesh = axon_vol.astype(np.float32)
    if SMOOTH:
        vol_for_mesh = binary_dilation(axon_vol, iterations=1).astype(np.float32)
        vol_for_mesh = gaussian_filter(vol_for_mesh, sigma=SMOOTH_SIGMA)
    axon_verts, axon_faces, _, _ = marching_cubes(
        vol_for_mesh,
        level=0.5,
        spacing=SCALE,
    )

    viewer = napari.Viewer(title=profile_name)
    viewer.add_surface(
        (axon_verts, axon_faces), name='axon',
        colormap='cyan', opacity=0.3,
    )

    vol_for_mesh = mito_vol.astype(np.float32)
    if SMOOTH:
        vol_for_mesh = binary_dilation(mito_vol, iterations=1).astype(np.float32)
        vol_for_mesh = gaussian_filter(vol_for_mesh, sigma=SMOOTH_SIGMA)
    mito_verts, mito_faces, _, _ = marching_cubes(
        vol_for_mesh,
        level=0.5,
        spacing=SCALE,
    )
    viewer.add_surface(
        (mito_verts, mito_faces), name='mitochondria',
        colormap='magenta', opacity=0.6,
    )

    viewer.add_points(
        points, name='vesicles', scale=SCALE,
        size=4, face_color='yellow', border_color='red',
        blending='additive',
    )
    viewer.dims.ndisplay = 3

napari.run()
