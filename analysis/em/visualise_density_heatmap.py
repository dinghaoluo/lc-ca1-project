# -*- coding: utf-8 -*-
'''
Created on Wed Jun  3 11:30:00 2026

3D visualisation of local vesicle density as a heatmap on the axon volume

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
DATA_ROOT = Path(r'Z:\EM-DATA\Annotations of profiles')

PROFILES = {
    '1st profile': {
        'tiff_pattern': '1st profile segmentation*.tiff',
        'csv': '1st profile vesicles annotation.csv',
    },
    '2nd profile': {
        'tiff_pattern': 'Segmentation 2nd profile*.tiff',
        'csv': 'Annotation 2nd profile.csv',
    },
    '3rd profile': {
        'tiff_pattern': 'Segmentation 3rd profile*.tiff',
        'csv': 'Annotation 3rd profile.csv',
    },
}

SCALE = (50, 5, 5)  # Z, Y, X in nm
KDE_SIGMA = (2, 8, 8)  # smoothing sigma in voxels for density field
CMAP = 'hot'


#%% load and visualise
for profile_name, info in PROFILES.items():
    profile_dir = DATA_ROOT / profile_name

    tiff_files = sorted(profile_dir.glob(info['tiff_pattern']))
    volume = np.stack([tifffile.imread(f) for f in tiff_files])
    volume = (volume > 0).astype(np.uint8)

    df = pd.read_csv(profile_dir / info['csv'])

    density_vol = np.zeros_like(volume, dtype=np.float32)
    for _, row in df.iterrows():
        xi = int(round(row['PosX'] / SCALE[2]))
        yi = int(round(row['PosY'] / SCALE[1]))
        zi = int(round(row['PosZ'] / SCALE[0]))
        density_vol[zi, yi, xi] += 1

    density_vol = gaussian_filter(density_vol, sigma=KDE_SIGMA)

    density_vol[volume == 0] = 0

    dmax = density_vol.max()
    density_vol /= dmax

    vol_dilated = binary_dilation(volume, iterations=1).astype(np.float32)
    vol_smooth = gaussian_filter(vol_dilated, sigma=(2, 1, 1))
    verts, faces, _, _ = marching_cubes(vol_smooth, level=0.5, spacing=SCALE)

    vert_voxels = np.round(verts / np.array(SCALE)).astype(int)
    vert_voxels[:, 0] = np.clip(vert_voxels[:, 0], 0, volume.shape[0] - 1)
    vert_voxels[:, 1] = np.clip(vert_voxels[:, 1], 0, volume.shape[1] - 1)
    vert_voxels[:, 2] = np.clip(vert_voxels[:, 2], 0, volume.shape[2] - 1)
    vert_density = density_vol[
        vert_voxels[:, 0], vert_voxels[:, 1], vert_voxels[:, 2]
    ]

    viewer = napari.Viewer(title=f'{profile_name} density')
    viewer.add_surface(
        (verts, faces, vert_density), name='density',
        colormap=CMAP, opacity=0.8,
    )
    viewer.add_points(
        np.column_stack([
            df['PosZ'].values / SCALE[0],
            df['PosY'].values / SCALE[1],
            df['PosX'].values / SCALE[2],
        ]),
        name='vesicles', scale=SCALE,
        size=3, face_color='white', border_color='white',
        blending='additive',
    )
    viewer.dims.ndisplay = 3

napari.run()
