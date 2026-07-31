# -*- coding: utf-8 -*-
'''
Created on Wed Jun  3 10:15:00 2026

vesicle density and distance-to-membrane analysis for EM axon profiles
    with mitochondria-aware cytoplasmic density

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
VOXEL_VOL_NM3 = SCALE_XY * SCALE_XY * SCALE_Z

plt.rcParams['font.family'] = 'Arial'


#%% load data
results = {}
for name in PROFILES:
    info = PROFILES[name]
    profile_dir = info['dir']
    axon_files = sorted(profile_dir.glob(info['axon_pattern']))
    axon_vol = (np.stack([tifffile.imread(f) for f in axon_files]) > 0).astype(np.uint8)
    mito_files = sorted(profile_dir.glob(info['mito_pattern']))
    mito_vol = (np.stack([tifffile.imread(f) for f in mito_files]) > 0).astype(np.uint8)
    coords_nm = pd.read_csv(
        profile_dir / info['csv'],
        sep=info['csv_sep']
    )[['PosX', 'PosY', 'PosZ']].values
    if not info['coords_in_nm']:
        coords_nm = coords_nm * 1000

    axon_voxels = np.sum(axon_vol > 0)
    mito_voxels = np.sum(mito_vol > 0)
    cyto_voxels = axon_voxels - mito_voxels

    n_ves = len(coords_nm)
    axon_vol_um3 = axon_voxels * VOXEL_VOL_NM3 / 1e9
    mito_vol_um3 = mito_voxels * VOXEL_VOL_NM3 / 1e9
    cyto_vol_um3 = cyto_voxels * VOXEL_VOL_NM3 / 1e9
    density_total = n_ves / axon_vol_um3
    density_cyto = n_ves / cyto_vol_um3

    # distance to axon membrane
    dist_field = distance_transform_edt(
        axon_vol, sampling=(SCALE_Z, SCALE_XY, SCALE_XY))
    ves_dists = np.zeros(len(coords_nm))
    for i, (px, py, pz) in enumerate(coords_nm):
        xi = int(round(px / SCALE_XY))
        yi = int(round(py / SCALE_XY))
        zi = int(round(pz / SCALE_Z))
        zi = np.clip(zi, 0, axon_vol.shape[0] - 1)
        yi = np.clip(yi, 0, axon_vol.shape[1] - 1)
        xi = np.clip(xi, 0, axon_vol.shape[2] - 1)
        ves_dists[i] = dist_field[zi, yi, xi]

    interior_dists = dist_field[axon_vol > 0].ravel()

    results[name] = {
        'axon_vol': axon_vol,
        'mito_vol': mito_vol,
        'coords_nm': coords_nm,
        'axon_vol_um3': axon_vol_um3,
        'mito_vol_um3': mito_vol_um3,
        'cyto_vol_um3': cyto_vol_um3,
        'density_total': density_total,
        'density_cyto': density_cyto,
        'mito_fraction': mito_voxels / axon_voxels,
        'ves_membrane_dists': ves_dists,
        'interior_dists': interior_dists,
    }
    print(f'{name}: {n_ves} vesicles, '
          f'axon={axon_vol_um3:.4f} um^3, '
          f'mito={100*mito_voxels/axon_voxels:.1f}%, '
          f'density(total)={density_total:.0f}, '
          f'density(cyto)={density_cyto:.0f} ves/um^3')


#%% density bar chart: total vs cytoplasmic
n_profiles = len(PROFILES)
names = list(PROFILES.keys())

fig, ax = plt.subplots(figsize=(4, 2.8), constrained_layout=True)
x = np.arange(n_profiles)
w = 0.35

d_total = [results[n]['density_total'] for n in names]
d_cyto = [results[n]['density_cyto'] for n in names]

ax.bar(x - w/2, d_total, w, color='grey', edgecolor='white',
       linewidth=0.5, label='total axon')
ax.bar(x + w/2, d_cyto, w, color='forestgreen', edgecolor='white',
       linewidth=0.5, label='cytoplasm only')

ax.set_xticks(x)
ax.set_xticklabels([f'P{i+1}' for i in range(n_profiles)], fontsize=7)
ax.set_ylabel('vesicles / μm³', fontsize=8)
ax.spines[['top', 'right']].set_visible(False)
ax.legend(fontsize=6.5, frameon=False)

for i in range(n_profiles):
    pct_increase = (d_cyto[i] - d_total[i]) / d_total[i] * 100
    ax.text(x[i] + w/2, d_cyto[i] + 20,
            f'+{pct_increase:.0f}%', ha='center', va='bottom', fontsize=5.5,
            color='forestgreen')

for ext in ['.pdf', '.png']:
    fig.savefig(FIGURES_STEM / f'vesicle_density_bar{ext}',
                dpi=300, bbox_inches='tight')
plt.show()


#%% mitochondria volume fraction
fig, ax = plt.subplots(figsize=(2.5, 2.5), constrained_layout=True)
mito_fracs = [results[n]['mito_fraction'] * 100 for n in names]
ax.bar(range(n_profiles), mito_fracs, color='orchid', edgecolor='white',
       linewidth=0.5, width=0.6)
ax.set_xticks(range(n_profiles))
ax.set_xticklabels([f'P{i+1}' for i in range(n_profiles)], fontsize=7)
ax.set_ylabel('mitochondria (% axon vol)', fontsize=8)
ax.spines[['top', 'right']].set_visible(False)

for i, f in enumerate(mito_fracs):
    ax.text(i, f + 0.5, f'{f:.1f}%', ha='center', va='bottom', fontsize=6.5)

for ext in ['.pdf', '.png']:
    fig.savefig(FIGURES_STEM / f'mitochondria_volume_fraction{ext}',
                dpi=300, bbox_inches='tight')
plt.show()


#%% distance-to-membrane histograms with CSR comparison
fig, axes = plt.subplots(1, n_profiles, figsize=(2.5 * n_profiles, 2.5),
                         constrained_layout=True)

for idx, name in enumerate(names):
    ves_d = results[name]['ves_membrane_dists']
    int_d = results[name]['interior_dists']

    ax = axes[idx]

    bins = np.linspace(0, max(ves_d.max(), np.percentile(int_d, 99)), 25)
    ax.hist(ves_d, bins=bins, density=True, color='grey', edgecolor='white',
            linewidth=0.5, alpha=0.7, label='observed')

    rng = np.random.default_rng(42)
    n_sub = min(50000, len(int_d))
    int_sub = rng.choice(int_d, size=n_sub, replace=False)
    ax.hist(int_sub, bins=bins, density=True, color='royalblue', alpha=0.25,
            edgecolor='none', label='CSR (uniform)')

    ax.axvline(np.median(ves_d), color='black', linestyle='--', linewidth=0.8)
    ax.set(xlabel='distance to membrane (nm)', ylabel='density')
    ax.set_title(name, fontsize=8)
    ax.spines[['top', 'right']].set_visible(False)
    ax.legend(fontsize=6, frameon=False)

    print(f'{name}: median dist to membrane = {np.median(ves_d):.1f} nm, '
          f'mean = {np.mean(ves_d):.1f} nm')

for ext in ['.pdf', '.png']:
    fig.savefig(FIGURES_STEM / f'vesicle_distance_to_membrane{ext}',
                dpi=300, bbox_inches='tight')
plt.show()


#%% distance to nearest mitochondrion
from scipy.stats import ks_2samp

fig, axes = plt.subplots(1, n_profiles, figsize=(2.5 * n_profiles, 2.5),
                         constrained_layout=True)

for idx, name in enumerate(names):
    mito_vol = results[name]['mito_vol']
    axon_vol = results[name]['axon_vol']
    coords_nm = results[name]['coords_nm']

    if mito_vol.max() == 0:
        axes[idx].set_title(f'{name}\n(no mito)', fontsize=8)
        continue

    # distance from mito surface outward into cytoplasm
    mito_dist_field = distance_transform_edt(
        ~(mito_vol > 0), sampling=(SCALE_Z, SCALE_XY, SCALE_XY))

    # vesicle distances to nearest mito
    ves_mito_d = np.zeros(len(coords_nm))
    for i, (px, py, pz) in enumerate(coords_nm):
        xi = int(round(px / SCALE_XY))
        yi = int(round(py / SCALE_XY))
        zi = int(round(pz / SCALE_Z))
        zi = np.clip(zi, 0, mito_vol.shape[0] - 1)
        yi = np.clip(yi, 0, mito_vol.shape[1] - 1)
        xi = np.clip(xi, 0, mito_vol.shape[2] - 1)
        ves_mito_d[i] = mito_dist_field[zi, yi, xi]

    # CSR: distance-to-mito for all cytoplasmic voxels
    cyto_mask = (axon_vol > 0) & (~(mito_vol > 0))
    cyto_mito_d = mito_dist_field[cyto_mask].ravel()

    ax = axes[idx]
    bins = np.linspace(
        0, max(ves_mito_d.max(), np.percentile(cyto_mito_d, 99)), 25)
    ax.hist(ves_mito_d, bins=bins, density=True, color='grey',
            edgecolor='white', linewidth=0.5, alpha=0.7, label='observed')

    rng = np.random.default_rng(99)
    n_sub = min(50000, len(cyto_mito_d))
    cyto_sub = rng.choice(cyto_mito_d, size=n_sub, replace=False)
    ax.hist(cyto_sub, bins=bins, density=True, color='orchid', alpha=0.25,
            edgecolor='none', label='CSR (cytoplasm)')

    ax.axvline(np.median(ves_mito_d), color='black', linestyle='--',
               linewidth=0.8)
    ax.set(xlabel='distance to mito (nm)', ylabel='density')
    ax.set_title(name, fontsize=8)
    ax.spines[['top', 'right']].set_visible(False)
    ax.legend(fontsize=6, frameon=False)

    print(f'{name}: median dist to mito = {np.median(ves_mito_d):.1f} nm')

for ext in ['.pdf', '.png']:
    fig.savefig(FIGURES_STEM / f'vesicle_distance_to_mito{ext}',
                dpi=300, bbox_inches='tight')
plt.show()


#%% mito-shell geometry: what fraction of each shell is cytoplasm?
SHELL_EDGES = np.arange(0, 351, 50)  # 0-50, 50-100, ..., 300-350 nm

fig, axes = plt.subplots(1, n_profiles, figsize=(2.5 * n_profiles, 2.8),
                         constrained_layout=True)

for idx, name in enumerate(names):
    mito_vol = results[name]['mito_vol']
    axon_vol = results[name]['axon_vol']

    if mito_vol.max() == 0:
        axes[idx].set_title(f'{name}\n(no mito)', fontsize=8)
        continue

    mito_dist_field = distance_transform_edt(
        ~(mito_vol > 0), sampling=(SCALE_Z, SCALE_XY, SCALE_XY))

    axon_mask = axon_vol > 0
    mito_mask = mito_vol > 0

    shell_cyto_frac = []
    shell_ves_frac = []
    coords_nm = results[name]['coords_nm']

    # vesicle distances to mito (recompute quickly from field)
    ves_mito_d = np.zeros(len(coords_nm))
    for i, (px, py, pz) in enumerate(coords_nm):
        xi = int(round(px / SCALE_XY))
        yi = int(round(py / SCALE_XY))
        zi = int(round(pz / SCALE_Z))
        zi = np.clip(zi, 0, mito_vol.shape[0] - 1)
        yi = np.clip(yi, 0, mito_vol.shape[1] - 1)
        xi = np.clip(xi, 0, mito_vol.shape[2] - 1)
        ves_mito_d[i] = mito_dist_field[zi, yi, xi]

    total_cyto = np.sum(axon_mask & ~mito_mask)
    n_ves = len(coords_nm)

    for j in range(len(SHELL_EDGES) - 1):
        lo, hi = SHELL_EDGES[j], SHELL_EDGES[j + 1]
        shell_mask = (mito_dist_field >= lo) & (mito_dist_field < hi)
        cyto_in_shell = np.sum(shell_mask & axon_mask & ~mito_mask)
        shell_cyto_frac.append(cyto_in_shell / total_cyto * 100)

        ves_in_shell = np.sum((ves_mito_d >= lo) & (ves_mito_d < hi))
        shell_ves_frac.append(ves_in_shell / n_ves * 100)

    ax = axes[idx]
    x = np.arange(len(shell_cyto_frac))
    w = 0.35
    ax.bar(x - w/2, shell_cyto_frac, w, color='royalblue', alpha=0.6,
           label='% cytoplasm')
    ax.bar(x + w/2, shell_ves_frac, w, color='gold', edgecolor='red',
           linewidth=0.4, label='% vesicles')
    labels = [f'{SHELL_EDGES[j]}-{SHELL_EDGES[j+1]}'
              for j in range(len(SHELL_EDGES) - 1)]
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=5.5, rotation=45, ha='right')
    ax.set_ylabel('% of total', fontsize=7)
    ax.set_title(name, fontsize=8)
    ax.spines[['top', 'right']].set_visible(False)
    ax.legend(fontsize=6, frameon=False)

    print(f'{name}: shell cytoplasm fracs = '
          + ', '.join(f'{v:.1f}%' for v in shell_cyto_frac))
    print(f'{name}: shell vesicle fracs   = '
          + ', '.join(f'{v:.1f}%' for v in shell_ves_frac))

for ext in ['.pdf', '.png']:
    fig.savefig(FIGURES_STEM / f'mito_shell_geometry{ext}',
                dpi=300, bbox_inches='tight')
plt.show()


#%% CDF of membrane distances with KS test
fig, axes = plt.subplots(1, n_profiles, figsize=(2.5 * n_profiles, 2.5),
                         constrained_layout=True)

for idx, name in enumerate(names):
    ves_d = results[name]['ves_membrane_dists']
    int_d = results[name]['interior_dists']

    rng = np.random.default_rng(42)
    n_sub = min(50000, len(int_d))
    int_sub = rng.choice(int_d, size=n_sub, replace=False)

    ves_sorted = np.sort(ves_d)
    csr_sorted = np.sort(int_sub)
    ves_cdf = np.arange(1, len(ves_sorted)+1) / len(ves_sorted)
    csr_cdf = np.arange(1, len(csr_sorted)+1) / len(csr_sorted)

    ax = axes[idx]
    ax.plot(ves_sorted, ves_cdf, color='black', linewidth=1.2,
            label='observed')
    ax.plot(csr_sorted, csr_cdf, color='royalblue', linewidth=1, alpha=0.7,
            label='CSR (uniform)')

    stat, pval = ks_2samp(ves_d, int_sub)
    pstr = 'p < 0.001' if pval < 0.001 else f'p = {pval:.3f}'
    ax.text(0.95, 0.15, f'KS = {stat:.3f}\n{pstr}',
            transform=ax.transAxes, ha='right', fontsize=6,
            family='monospace')

    ax.set(xlabel='distance to membrane (nm)', ylabel='CDF')
    ax.set_title(name, fontsize=8)
    ax.spines[['top', 'right']].set_visible(False)
    ax.legend(fontsize=6, frameon=False, loc='lower right')

for ext in ['.pdf', '.png']:
    fig.savefig(FIGURES_STEM / f'vesicle_distance_to_membrane_cdf{ext}',
                dpi=300, bbox_inches='tight')
plt.show()
