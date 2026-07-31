# -*- coding: utf-8 -*-
'''
Created on Wed Jun  3 10:45:00 2026

vesicle clustering analysis: nearest-neighbour distances and Ripley's K
    against complete spatial randomness (CSR) within the axon volume

@author: Dinghao Luo
'''

#%% imports
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import tifffile
import matplotlib.pyplot as plt
from scipy.ndimage import distance_transform_edt
from scipy.spatial import KDTree

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

SCALE_XY = 5      # nm/pixel
SCALE_Z = 50      # nm/slice
N_SIMULATIONS = 500

plt.rcParams['font.family'] = 'Arial'


#%% spatial statistics
def sample_random_points_in_volume(volume, n_points, rng):
    axon_coords = np.argwhere(volume > 0)
    chosen = rng.choice(len(axon_coords), size=n_points, replace=True)
    voxels = axon_coords[chosen].astype(np.float64)
    voxels += rng.uniform(-0.5, 0.5, size=voxels.shape)

    pts_nm = np.column_stack([
        voxels[:, 2] * SCALE_XY,
        voxels[:, 1] * SCALE_XY,
        voxels[:, 0] * SCALE_Z,
    ])
    return pts_nm

def ripleys_k(coords_nm, radii, volume_um3):
    n = len(coords_nm)
    tree = KDTree(coords_nm)
    k_vals = np.zeros_like(radii, dtype=float)

    for i, r in enumerate(radii):
        counts = tree.query_ball_point(coords_nm, r, return_length=True)
        k_vals[i] = np.sum(counts - 1)

    k_vals = (volume_um3 * 1e9) / (n ** 2) * k_vals
    return k_vals


#%% load data
results = {}
for name in PROFILES:
    info = PROFILES[name]
    profile_dir = info['dir']
    axon_files = sorted(profile_dir.glob(info['axon_pattern']))
    axon_vol = (np.stack([tifffile.imread(f) for f in axon_files]) > 0).astype(np.uint8)
    mito_files = sorted(profile_dir.glob(info['mito_pattern']))
    mito_vol = (np.stack([tifffile.imread(f) for f in mito_files]) > 0).astype(np.uint8)
    cyto_vol = axon_vol & (~mito_vol)
    coords_nm = pd.read_csv(
        profile_dir / info['csv'],
        sep=info['csv_sep']
    )[['PosX', 'PosY', 'PosZ']].values
    if not info['coords_in_nm']:
        coords_nm = coords_nm * 1000
    cyto_voxels = np.sum(cyto_vol > 0)
    cyto_vol_um3 = cyto_voxels * SCALE_XY**2 * SCALE_Z / 1e9
    results[name] = {
        'axon_vol': axon_vol,
        'volume': cyto_vol,
        'coords_nm': coords_nm,
        'cyto_vol_um3': cyto_vol_um3,
    }

n_profiles = len(PROFILES)


#%% nearest-neighbour distance analysis
fig, axes = plt.subplots(1, n_profiles, figsize=(2.2 * n_profiles, 2.2),
                         constrained_layout=True)
rng = np.random.default_rng(42)

for idx, name in enumerate(PROFILES):
    coords = results[name]['coords_nm']
    volume = results[name]['volume']
    n_pts = len(coords)

    nn_obs = KDTree(coords).query(coords, k=2)[0][:, 1]

    nn_sims = np.zeros((N_SIMULATIONS, n_pts))
    for sim in range(N_SIMULATIONS):
        rand_pts = sample_random_points_in_volume(volume, n_pts, rng=rng)
        nn_sims[sim] = KDTree(rand_pts).query(rand_pts, k=2)[0][:, 1]

    sim_medians = np.median(nn_sims, axis=1)
    csr_lo = np.percentile(sim_medians, 2.5)
    csr_hi = np.percentile(sim_medians, 97.5)

    ax = axes[idx]
    ax.hist(nn_obs, bins=20, density=True, color='grey', edgecolor='white',
            linewidth=0.5, alpha=0.7, label='observed')
    all_csr = nn_sims.ravel()
    ax.hist(all_csr, bins=40, density=True, color='royalblue', alpha=0.3,
            edgecolor='none', label='CSR')
    ax.axvline(np.median(nn_obs), color='black', linestyle='--', linewidth=1)
    ax.set(xlabel='NN distance (nm)', ylabel='density')
    ax.set_title(name, fontsize=8)
    ax.spines[['top', 'right']].set_visible(False)
    ax.legend(fontsize=6, frameon=False)

    print(f'{name}: median NN = {np.median(nn_obs):.1f} nm, '
          f'CSR 95% CI = [{csr_lo:.1f}, {csr_hi:.1f}] nm')

for ext in ['.pdf', '.png']:
    fig.savefig(FIGURES_STEM / f'vesicle_nn_distances{ext}',
                dpi=300, bbox_inches='tight')
plt.show()


#%% Ripley's L function
fig, axes = plt.subplots(1, n_profiles, figsize=(2.2 * n_profiles, 2.2),
                         constrained_layout=True)

radii = np.linspace(20, 300, 30)
rng = np.random.default_rng(123)

for idx, name in enumerate(PROFILES):
    coords = results[name]['coords_nm']
    volume = results[name]['volume']
    vol_um3 = results[name]['cyto_vol_um3']
    n_pts = len(coords)

    k_obs = ripleys_k(coords, radii, vol_um3)
    l_obs = np.cbrt(3 * k_obs / (4 * np.pi)) - radii

    l_sims = np.zeros((N_SIMULATIONS, len(radii)))
    for sim in range(N_SIMULATIONS):
        rand_pts = sample_random_points_in_volume(volume, n_pts, rng=rng)
        k_sim = ripleys_k(rand_pts, radii, vol_um3)
        l_sims[sim] = np.cbrt(3 * k_sim / (4 * np.pi)) - radii

    l_lo = np.percentile(l_sims, 2.5, axis=0)
    l_hi = np.percentile(l_sims, 97.5, axis=0)

    ax = axes[idx]
    ax.fill_between(radii, l_lo, l_hi, color='royalblue', alpha=0.2,
                    label='CSR 95% envelope')
    ax.plot(radii, l_obs, color='black', linewidth=1.2, label='observed')
    ax.axhline(0, color='grey', linestyle=':', linewidth=0.7)
    ax.set(xlabel='r (nm)', ylabel='L(r) - r')
    ax.set_title(name, fontsize=8)
    ax.spines[['top', 'right']].set_visible(False)
    ax.legend(fontsize=6, frameon=False)

    # significance: at each radius, proportion of sims with L >= observed
    p_vals = np.mean(l_sims >= l_obs[None, :], axis=0)
    n_sig = np.sum(p_vals < 0.05)
    print(f'{name}: L(r) significant at {n_sig}/{len(radii)} radii '
          f'(p<0.05 at r={radii[p_vals < 0.05][0]:.0f}-'
          f'{radii[p_vals < 0.05][-1]:.0f} nm)'
          if n_sig > 0 else
          f'{name}: L(r) not significant at any radius')

for ext in ['.pdf', '.png']:
    fig.savefig(FIGURES_STEM / f'vesicle_ripleys_L{ext}',
                dpi=300, bbox_inches='tight')
plt.show()


#%% pairwise distance distribution
fig, axes = plt.subplots(1, n_profiles, figsize=(2.2 * n_profiles, 2.2),
                         constrained_layout=True)

for idx, name in enumerate(PROFILES):
    coords = results[name]['coords_nm']
    tree = KDTree(coords)
    pw_dists = tree.sparse_distance_matrix(tree, max_distance=500).toarray()
    pw_dists = pw_dists[np.triu_indices_from(pw_dists, k=1)]
    pw_dists = pw_dists[pw_dists > 0]

    ax = axes[idx]
    ax.hist(pw_dists, bins=30, color='grey', edgecolor='white',
            linewidth=0.5, density=True)
    ax.set(xlabel='pairwise distance (nm)', ylabel='density')
    ax.set_title(name, fontsize=8)
    ax.spines[['top', 'right']].set_visible(False)

for ext in ['.pdf', '.png']:
    fig.savefig(FIGURES_STEM / f'vesicle_pairwise_distances{ext}',
                dpi=300, bbox_inches='tight')
plt.show()


#%% cluster-to-membrane: are clustered vesicles near the membrane?
fig, axes = plt.subplots(1, n_profiles, figsize=(2.5 * n_profiles, 2.5),
                         constrained_layout=True)

for idx, name in enumerate(PROFILES):
    coords = results[name]['coords_nm']
    axon_vol = results[name]['axon_vol']

    membrane_dist = distance_transform_edt(
        axon_vol, sampling=(SCALE_Z, SCALE_XY, SCALE_XY))

    ves_mem_d = np.zeros(len(coords))
    for i, (px, py, pz) in enumerate(coords):
        xi = int(round(px / SCALE_XY))
        yi = int(round(py / SCALE_XY))
        zi = int(round(pz / SCALE_Z))
        zi = np.clip(zi, 0, axon_vol.shape[0] - 1)
        yi = np.clip(yi, 0, axon_vol.shape[1] - 1)
        xi = np.clip(xi, 0, axon_vol.shape[2] - 1)
        ves_mem_d[i] = membrane_dist[zi, yi, xi]

    # split by NN distance: 'clustered' = NN < median
    nn_dists = KDTree(coords).query(coords, k=2)[0][:, 1]
    median_nn = np.median(nn_dists)
    clustered = nn_dists < median_nn
    dispersed = ~clustered

    ax = axes[idx]
    bins = np.linspace(0, np.percentile(ves_mem_d, 99), 20)
    ax.hist(ves_mem_d[clustered], bins=bins, density=True, alpha=0.7,
            color='tomato', edgecolor='white', linewidth=0.5,
            label=f'clustered (NN<{median_nn:.0f}nm)')
    ax.hist(ves_mem_d[dispersed], bins=bins, density=True, alpha=0.4,
            color='steelblue', edgecolor='none',
            label=f'dispersed (NN>{median_nn:.0f}nm)')
    ax.axvline(np.median(ves_mem_d[clustered]), color='tomato',
               linestyle='--', linewidth=1)
    ax.axvline(np.median(ves_mem_d[dispersed]), color='steelblue',
               linestyle='--', linewidth=1)
    ax.set(xlabel='distance to membrane (nm)', ylabel='density')
    ax.set_title(name, fontsize=8)
    ax.spines[['top', 'right']].set_visible(False)
    ax.legend(fontsize=5.5, frameon=False)

    print(f'{name}: clustered median-to-membrane = '
          f'{np.median(ves_mem_d[clustered]):.1f} nm, '
          f'dispersed = {np.median(ves_mem_d[dispersed]):.1f} nm')

for ext in ['.pdf', '.png']:
    fig.savefig(FIGURES_STEM / f'cluster_to_membrane{ext}',
                dpi=300, bbox_inches='tight')
plt.show()
