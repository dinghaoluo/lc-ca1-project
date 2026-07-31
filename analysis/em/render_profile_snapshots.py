# -*- coding: utf-8 -*-
'''
Created on Thu Jun 18 11:00:00 2026

render 3D profile figures (axon + mitochondria + vesicles) using matplotlib

@author: Dinghao Luo
'''

#%% imports
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import tifffile
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from scipy.ndimage import gaussian_filter, binary_dilation
from skimage.measure import marching_cubes

repo_root = Path(__file__).resolve().parents[2]
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))

import project_paths as pp


#%% paths and parameters
DATA_ROOT_ORIG = Path(r'Z:\EM-DATA\Annotations of profiles'
                      r'\new exports from Dragonfly 6-16-2026')
DATA_ROOT_REVISED = Path(r'Z:\EM-DATA\Annotations of profiles'
                         r'\new revised profile 3 and 4')
FIGURES_STEM = pp.FIGURES_ROOT / 'em'
FIGURES_STEM.mkdir(parents=True, exist_ok=True)

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
SMOOTH_SIGMA = (2, 1, 1)
plt.rcParams['font.family'] = 'Arial'


#%% render each profile
for profile_name, paths in PROFILES.items():
    profile_dir = paths['dir']

    axon_files = sorted(profile_dir.glob(paths['axon_pattern']))
    axon_vol = np.stack([tifffile.imread(f) for f in axon_files])
    axon_vol = (axon_vol > 0).astype(np.uint8)

    mito_files = sorted(profile_dir.glob(paths['mito_pattern']))
    mito_vol = np.stack([tifffile.imread(f) for f in mito_files])
    mito_vol = (mito_vol > 0).astype(np.uint8)
    has_mito = mito_vol.max() > 0

    df = pd.read_csv(profile_dir / paths['csv'], sep=paths['csv_sep'])
    scale = 1.0 if paths['coords_in_nm'] else 1000.0
    ves_z = df['PosZ'].values * scale
    ves_y = df['PosY'].values * scale
    ves_x = df['PosX'].values * scale

    vol_for_mesh = binary_dilation(axon_vol, iterations=1).astype(np.float32)
    vol_for_mesh = gaussian_filter(vol_for_mesh, sigma=SMOOTH_SIGMA)
    axon_verts, axon_faces, _, _ = marching_cubes(
        vol_for_mesh,
        level=0.5,
        spacing=SCALE,
    )

    fig = plt.figure(figsize=(6, 5))
    ax = fig.add_subplot(111, projection='3d')

    axon_coll = Poly3DCollection(
        axon_verts[axon_faces], alpha=0.15, linewidths=0)
    axon_coll.set_facecolor((0, 0.8, 0.8))
    axon_coll.set_edgecolor('none')
    ax.add_collection3d(axon_coll)

    if has_mito:
        vol_for_mesh = binary_dilation(mito_vol, iterations=1).astype(np.float32)
        vol_for_mesh = gaussian_filter(vol_for_mesh, sigma=SMOOTH_SIGMA)
        mito_verts, mito_faces, _, _ = marching_cubes(
            vol_for_mesh,
            level=0.5,
            spacing=SCALE,
        )
        mito_coll = Poly3DCollection(
            mito_verts[mito_faces], alpha=0.6, linewidths=0)
        mito_coll.set_facecolor((0.8, 0.2, 0.8))
        mito_coll.set_edgecolor('none')
        ax.add_collection3d(mito_coll)

    ax.scatter(ves_z, ves_y, ves_x, s=12, c='gold', edgecolors='red',
               linewidths=0.3, alpha=0.9, zorder=5)

    ax.set_xlim(axon_verts[:, 0].min(), axon_verts[:, 0].max())
    ax.set_ylim(axon_verts[:, 1].min(), axon_verts[:, 1].max())
    ax.set_zlim(axon_verts[:, 2].min(), axon_verts[:, 2].max())
    ax.set_xlabel('Z (nm)', fontsize=7)
    ax.set_ylabel('Y (nm)', fontsize=7)
    ax.set_zlabel('X (nm)', fontsize=7)
    ax.set_title(profile_name, fontsize=9)
    ax.tick_params(labelsize=6)
    ax.view_init(elev=25, azim=45)

    for ext in ['.pdf', '.png']:
        fig.savefig(
            FIGURES_STEM / f'{profile_name.replace(" ", "_")}_3d{ext}',
            dpi=300, bbox_inches='tight')
    plt.close(fig)
    print(f'saved: {profile_name}')

print('done')
