# -*- coding: utf-8 -*-
'''
Created on Tue Aug 19 12:00:00 2026

vesicle polarity analysis: testing whether vesicles inside the LC axon
    are positioned to face nearby neuropil synapses (directionality test)
    + literature-standard pool binning by membrane distance

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
from scipy.stats import mannwhitneyu, ks_2samp

repo_root = Path(__file__).resolve().parents[2]
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))

import project_paths as pp


#%% paths and parameters
DATA_ROOT_ORIG = Path(r'Z:\EM-DATA\Annotations of profiles'
                      r'\new exports from Dragonfly 6-16-2026')
DATA_ROOT_REVISED = Path(r'Z:\EM-DATA\Annotations of profiles'
                         r'\new revised profile 3 and 4')
SYNAPSE_ROOT = Path(r'Z:\EM-DATA\Annotations of profiles'
                    r'\Synapses center coordinates')
FIGURES_STEM = pp.FIGURES_ROOT / 'em'
FIGURES_STEM.mkdir(parents=True, exist_ok=True)

PROFILES = {
    'profile 1': {
        'dir': DATA_ROOT_ORIG / 'profile 1',
        'axon_pattern': 'Profile 1 Axon Segmentation*.tiff',
        'mito_pattern': 'Profile 1 Mitochondria Segmentation*.tiff',
        'ves_csv': 'Profile 1 Vesicle Positions.csv',
        'ves_sep': ',',
        'ves_in_nm': False,
        'syn_csv': SYNAPSE_ROOT / 'profile 1 synapse center coordinates.csv',
        'syn_sep': ',',
    },
    'profile 2': {
        'dir': DATA_ROOT_ORIG / 'profile 2',
        'axon_pattern': 'Profile 2 Axon Segmentation*.tiff',
        'mito_pattern': 'Profile 2 Mitochondria Segmentation*.tiff',
        'ves_csv': 'Profile 2 Vesicle Positions.csv',
        'ves_sep': ',',
        'ves_in_nm': False,
        'syn_csv': SYNAPSE_ROOT / 'Profile 2 synapse center coordinates.csv',
        'syn_sep': ';',
    },
    'profile 3': {
        'dir': DATA_ROOT_REVISED / 'new profile 3',
        'axon_pattern': 'Profile 3 Axon segmentation*.tiff',
        'mito_pattern': 'Profile 3 Mito segmentation*.tiff',
        'ves_csv': 'Profile 3 Vesicle Positions.csv',
        'ves_sep': ';',
        'ves_in_nm': True,
        'syn_csv': SYNAPSE_ROOT / 'Profile 3 Synapse center coordinates.csv',
        'syn_sep': ';',
    },
    'profile 4': {
        'dir': DATA_ROOT_REVISED / 'new profile 4',
        'axon_pattern': 'Profile 4 Axon Segmentation*.tiff',
        'mito_pattern': 'Profile 4 Mitochondria Segmentation*.tiff',
        'ves_csv': 'Profile 4 Vesicle Positions.csv',
        'ves_sep': ';',
        'ves_in_nm': True,
        'syn_csv': SYNAPSE_ROOT / 'Profile 4 Synapse center coordinates.csv',
        'syn_sep': ';',
    },
}

SCALE_XY = 5      # nm/pixel
SCALE_Z = 50      # nm/slice

# literature-informed pool boundaries (Bhatt et al. 2020; SynapsEM)
POOL_EDGES = [0, 30, 100, 200, np.inf]
POOL_LABELS = ['docked\n(0-30 nm)', 'proximal\n(30-100 nm)',
               'intermediate\n(100-200 nm)', 'reserve\n(>200 nm)']

plt.rcParams['font.family'] = 'Arial'


#%% load data
results = {}
for name, info in PROFILES.items():
    profile_dir = info['dir']

    axon_files = sorted(profile_dir.glob(info['axon_pattern']))
    axon_vol = (np.stack([tifffile.imread(f) for f in axon_files]) > 0
                ).astype(np.uint8)
    mito_files = sorted(profile_dir.glob(info['mito_pattern']))
    mito_vol = (np.stack([tifffile.imread(f) for f in mito_files]) > 0
                ).astype(np.uint8)

    # distance field + nearest-boundary indices
    dist_field, nearest_idx = distance_transform_edt(
        axon_vol, sampling=(SCALE_Z, SCALE_XY, SCALE_XY),
        return_indices=True)

    # vesicle coordinates -> nm
    ves_df = pd.read_csv(profile_dir / info['ves_csv'], sep=info['ves_sep'])
    ves_nm = ves_df[['PosX', 'PosY', 'PosZ']].values
    if not info['ves_in_nm']:
        ves_nm = ves_nm * 1000

    # synapse coordinates -> nm (all in micrometres in the CSVs)
    syn_df = pd.read_csv(info['syn_csv'], sep=info['syn_sep'])
    syn_nm = syn_df[['PosX', 'PosY', 'PosZ']].values * 1000

    results[name] = {
        'axon_vol': axon_vol,
        'dist_field': dist_field,
        'nearest_idx': nearest_idx,
        'ves_nm': ves_nm,
        'syn_nm': syn_nm,
    }
    print(f'{name}: {len(ves_nm)} vesicles, {len(syn_nm)} synapses loaded')


#%% compute polarity angles
# for each vesicle: angle between (vesicle → nearest membrane) and
# (vesicle → nearest synapse); if vesicles face synapses, angles cluster near 0
polarity_results = {}

for name, res in results.items():
    ves_nm = res['ves_nm']
    syn_nm = res['syn_nm']
    axon_vol = res['axon_vol']
    dist_field = res['dist_field']
    nearest_idx = res['nearest_idx']

    # vesicle voxel indices
    ves_vox = np.column_stack([
        np.clip((ves_nm[:, 2] / SCALE_Z).round().astype(int),
                0, axon_vol.shape[0] - 1),
        np.clip((ves_nm[:, 1] / SCALE_XY).round().astype(int),
                0, axon_vol.shape[1] - 1),
        np.clip((ves_nm[:, 0] / SCALE_XY).round().astype(int),
                0, axon_vol.shape[2] - 1),
    ])

    # nearest boundary coordinates in nm for each vesicle
    membrane_nm = np.zeros_like(ves_nm)
    ves_mem_dist = np.zeros(len(ves_nm))
    for i in range(len(ves_nm)):
        zi, yi, xi = ves_vox[i]
        bz = nearest_idx[0, zi, yi, xi]
        by = nearest_idx[1, zi, yi, xi]
        bx = nearest_idx[2, zi, yi, xi]
        membrane_nm[i] = [bx * SCALE_XY, by * SCALE_XY, bz * SCALE_Z]
        ves_mem_dist[i] = dist_field[zi, yi, xi]

    # direction vectors: vesicle → membrane and vesicle → nearest synapse
    vec_to_mem = membrane_nm - ves_nm
    norms_mem = np.linalg.norm(vec_to_mem, axis=1, keepdims=True)
    # avoid division by zero for vesicles at the membrane
    norms_mem[norms_mem == 0] = 1
    vec_to_mem_hat = vec_to_mem / norms_mem

    syn_tree = KDTree(syn_nm)
    ves_to_syn_dist, syn_idx = syn_tree.query(ves_nm)
    nearest_syn = syn_nm[syn_idx]
    vec_to_syn = nearest_syn - ves_nm
    norms_syn = np.linalg.norm(vec_to_syn, axis=1, keepdims=True)
    norms_syn[norms_syn == 0] = 1
    vec_to_syn_hat = vec_to_syn / norms_syn

    # angle between the two direction vectors
    dot = np.sum(vec_to_mem_hat * vec_to_syn_hat, axis=1)
    dot = np.clip(dot, -1, 1)
    angles_rad = np.arccos(dot)
    angles_deg = np.degrees(angles_rad)

    polarity_results[name] = {
        'angles_deg': angles_deg,
        'ves_mem_dist': ves_mem_dist,
        'ves_to_syn_dist': ves_to_syn_dist,
    }
    print(f'{name}: median polarity angle = {np.median(angles_deg):.1f}°, '
          f'mean = {np.mean(angles_deg):.1f}° (90° = random)')


#%% polarity angle distributions
fig, axes = plt.subplots(2, 2, figsize=(6, 5), constrained_layout=True)
axes = axes.ravel()

for idx, (name, res) in enumerate(results.items()):
    angles = polarity_results[name]['angles_deg']
    ax = axes[idx]

    bins = np.linspace(0, 180, 19)
    ax.hist(angles, bins=bins, density=True, color='grey', edgecolor='white',
            linewidth=0.5, alpha=0.8, label='observed')

    # uniform expectation on a sphere: P(θ) ∝ sin(θ)
    theta_mid = (bins[:-1] + bins[1:]) / 2
    sin_density = np.sin(np.radians(theta_mid))
    sin_density /= (sin_density * np.diff(bins)[0]).sum()
    ax.plot(theta_mid, sin_density, color='royalblue', linewidth=1.2,
            alpha=0.7, label='isotropic (sin θ)')

    ax.axvline(90, color='grey', linestyle=':', linewidth=0.7)
    ax.axvline(np.median(angles), color='black', linestyle='--', linewidth=0.8)

    # test: are angles smaller than expected under isotropy?
    # generate null: sample angles from sin(θ) distribution
    rng = np.random.default_rng(42 + idx)
    null_angles = np.degrees(np.arccos(1 - 2 * rng.random(10000)))
    _, p = mannwhitneyu(angles, null_angles, alternative='less')
    pstr = 'p < 0.001' if p < 0.001 else f'p = {p:.3f}'

    ax.text(0.95, 0.85, f'median: {np.median(angles):.0f}°\n'
            f'MWU vs isotropic\n{pstr}',
            transform=ax.transAxes, ha='right', va='top', fontsize=5.5,
            family='monospace')

    ax.set(xlabel='angle: membrane↔synapse (°)', ylabel='density')
    ax.set_title(f'{name} ({len(res["syn_nm"])} syn)', fontsize=8)
    ax.spines[['top', 'right']].set_visible(False)
    ax.legend(fontsize=6, frameon=False)

for ext in ['.pdf', '.png']:
    fig.savefig(FIGURES_STEM / f'vesicle_polarity_angle{ext}',
                dpi=300, bbox_inches='tight')
plt.show()


#%% pool-stratified synapse proximity
# bin vesicles by membrane distance (literature pools), then test
# whether each pool is closer to synapses than CSR
fig, axes = plt.subplots(2, 2, figsize=(6.5, 5.5), constrained_layout=True)
axes = axes.ravel()

for idx, (name, res) in enumerate(results.items()):
    ves_mem_dist = polarity_results[name]['ves_mem_dist']
    ves_syn_dist = polarity_results[name]['ves_to_syn_dist']

    pool_medians = []
    pool_counts = []
    for j in range(len(POOL_EDGES) - 1):
        lo, hi = POOL_EDGES[j], POOL_EDGES[j + 1]
        mask = (ves_mem_dist >= lo) & (ves_mem_dist < hi)
        n = mask.sum()
        pool_counts.append(n)
        if n > 0:
            pool_medians.append(np.median(ves_syn_dist[mask]))
        else:
            pool_medians.append(np.nan)

    ax = axes[idx]
    x = np.arange(len(POOL_LABELS))
    bars = ax.bar(x, pool_medians, color='steelblue', edgecolor='white',
                  linewidth=0.5, alpha=0.8)
    for i, (med, n) in enumerate(zip(pool_medians, pool_counts)):
        if not np.isnan(med):
            ax.text(i, med + 10, f'n={n}', ha='center', va='bottom',
                    fontsize=5.5)

    ax.set_xticks(x)
    ax.set_xticklabels(POOL_LABELS, fontsize=6)
    ax.set_ylabel('median dist to synapse (nm)', fontsize=7)
    ax.set_title(f'{name}', fontsize=8)
    ax.spines[['top', 'right']].set_visible(False)

    print(f'{name} pool counts: {pool_counts}, '
          f'median ves-to-syn per pool: '
          + ', '.join(f'{m:.0f}' if not np.isnan(m) else 'n/a'
                      for m in pool_medians))

for ext in ['.pdf', '.png']:
    fig.savefig(FIGURES_STEM / f'pool_stratified_synapse_distance{ext}',
                dpi=300, bbox_inches='tight')
plt.show()


#%% pool-stratified polarity: are membrane-proximal vesicles more polarized?
fig, axes = plt.subplots(2, 2, figsize=(6.5, 5.5), constrained_layout=True)
axes = axes.ravel()

for idx, (name, res) in enumerate(results.items()):
    ves_mem_dist = polarity_results[name]['ves_mem_dist']
    angles = polarity_results[name]['angles_deg']

    ax = axes[idx]
    pool_angle_medians = []
    for j in range(len(POOL_EDGES) - 1):
        lo, hi = POOL_EDGES[j], POOL_EDGES[j + 1]
        mask = (ves_mem_dist >= lo) & (ves_mem_dist < hi)
        if mask.sum() > 0:
            pool_angle_medians.append(np.median(angles[mask]))
        else:
            pool_angle_medians.append(np.nan)

    x = np.arange(len(POOL_LABELS))
    ax.bar(x, pool_angle_medians, color='coral', edgecolor='white',
           linewidth=0.5, alpha=0.8)
    ax.axhline(90, color='grey', linestyle=':', linewidth=0.7,
               label='isotropic expectation')
    ax.set_xticks(x)
    ax.set_xticklabels(POOL_LABELS, fontsize=6)
    ax.set_ylabel('median polarity angle (°)', fontsize=7)
    ax.set_title(f'{name}', fontsize=8)
    ax.spines[['top', 'right']].set_visible(False)
    ax.legend(fontsize=6, frameon=False)
    ax.set_ylim(0, 180)

    print(f'{name} pool polarity angles: '
          + ', '.join(f'{m:.0f}°' if not np.isnan(m) else 'n/a'
                      for m in pool_angle_medians))

for ext in ['.pdf', '.png']:
    fig.savefig(FIGURES_STEM / f'pool_stratified_polarity{ext}',
                dpi=300, bbox_inches='tight')
plt.show()


#%% axon surface synapse proximity map
# what fraction of the axon surface is within given distances of a synapse?
SURFACE_RADII = [200, 500, 1000, 2000]  # nm

fig, ax = plt.subplots(figsize=(5, 3), constrained_layout=True)

for name, res in results.items():
    axon_vol = res['axon_vol']
    syn_nm = res['syn_nm']

    # surface voxels: axon voxels adjacent to background
    from scipy.ndimage import binary_erosion
    interior = binary_erosion(axon_vol)
    surface_mask = axon_vol & (~interior)
    surface_vox = np.argwhere(surface_mask > 0)

    surface_nm = np.column_stack([
        surface_vox[:, 2] * SCALE_XY,
        surface_vox[:, 1] * SCALE_XY,
        surface_vox[:, 0] * SCALE_Z,
    ])

    syn_tree = KDTree(syn_nm)
    surf_to_syn, _ = syn_tree.query(surface_nm)

    fracs = [np.mean(surf_to_syn <= r) * 100 for r in SURFACE_RADII]
    ax.plot(SURFACE_RADII, fracs, 'o-', markersize=4, label=name)

    print(f'{name}: surface fraction within '
          + ', '.join(f'{r}nm={f:.1f}%' for r, f in zip(SURFACE_RADII, fracs)))

ax.set_xlabel('distance threshold (nm)', fontsize=8)
ax.set_ylabel('% axon surface within threshold of synapse', fontsize=8)
ax.spines[['top', 'right']].set_visible(False)
ax.legend(fontsize=6.5, frameon=False)

for ext in ['.pdf', '.png']:
    fig.savefig(FIGURES_STEM / f'axon_surface_synapse_proximity{ext}',
                dpi=300, bbox_inches='tight')
plt.show()


#%% vesicle enrichment at synapse-facing surface segments
# split axon surface into "synapse-facing" (<500nm from a synapse) and
# "non-facing" (>500nm); compare vesicle density near each
FACING_THRESHOLD = 500  # nm
VESICLE_SHELL = 100     # count vesicles within this distance of the surface

fig, axes = plt.subplots(1, 4, figsize=(9, 2.5), constrained_layout=True)

for idx, (name, res) in enumerate(results.items()):
    axon_vol = res['axon_vol']
    ves_nm = res['ves_nm']
    syn_nm = res['syn_nm']
    dist_field = res['dist_field']
    nearest_idx = res['nearest_idx']

    # get each vesicle's nearest membrane point
    ves_vox = np.column_stack([
        np.clip((ves_nm[:, 2] / SCALE_Z).round().astype(int),
                0, axon_vol.shape[0] - 1),
        np.clip((ves_nm[:, 1] / SCALE_XY).round().astype(int),
                0, axon_vol.shape[1] - 1),
        np.clip((ves_nm[:, 0] / SCALE_XY).round().astype(int),
                0, axon_vol.shape[2] - 1),
    ])

    membrane_pts = np.zeros((len(ves_nm), 3))
    for i in range(len(ves_nm)):
        zi, yi, xi = ves_vox[i]
        bz = nearest_idx[0, zi, yi, xi]
        by = nearest_idx[1, zi, yi, xi]
        bx = nearest_idx[2, zi, yi, xi]
        membrane_pts[i] = [bx * SCALE_XY, by * SCALE_XY, bz * SCALE_Z]

    # classify each vesicle's membrane projection point by synapse proximity
    syn_tree = KDTree(syn_nm)
    proj_to_syn, _ = syn_tree.query(membrane_pts)

    facing = proj_to_syn <= FACING_THRESHOLD
    non_facing = ~facing

    ax = axes[idx]
    counts = [facing.sum(), non_facing.sum()]
    labels = [f'facing\n(<{FACING_THRESHOLD}nm)', f'non-facing\n(>{FACING_THRESHOLD}nm)']
    bars = ax.bar(labels, counts, color=['forestgreen', 'lightgrey'],
                  edgecolor='white', linewidth=0.5)
    total = len(ves_nm)
    for i, c in enumerate(counts):
        ax.text(i, c + 0.5, f'{c/total*100:.0f}%', ha='center',
                va='bottom', fontsize=6.5)

    ax.set_ylabel('vesicle count', fontsize=7)
    ax.set_title(f'{name}', fontsize=8)
    ax.spines[['top', 'right']].set_visible(False)

    print(f'{name}: {facing.sum()}/{total} vesicles '
          f'({facing.sum()/total*100:.0f}%) project to synapse-facing surface')

for ext in ['.pdf', '.png']:
    fig.savefig(FIGURES_STEM / f'vesicle_facing_synapse{ext}',
                dpi=300, bbox_inches='tight')
plt.show()


#%% compare facing fraction to surface fraction (enrichment ratio)
print('\n--- enrichment test ---')
print('if vesicles are randomly placed, the fraction projecting to '
      f'synapse-facing surface should match the surface fraction within '
      f'{FACING_THRESHOLD}nm')

for name, res in results.items():
    axon_vol = res['axon_vol']
    syn_nm = res['syn_nm']
    ves_nm = res['ves_nm']

    from scipy.ndimage import binary_erosion
    interior = binary_erosion(axon_vol)
    surface_mask = axon_vol & (~interior)
    surface_vox = np.argwhere(surface_mask > 0)
    surface_nm = np.column_stack([
        surface_vox[:, 2] * SCALE_XY,
        surface_vox[:, 1] * SCALE_XY,
        surface_vox[:, 0] * SCALE_Z,
    ])
    syn_tree = KDTree(syn_nm)
    surf_to_syn, _ = syn_tree.query(surface_nm)
    surface_facing_frac = np.mean(surf_to_syn <= FACING_THRESHOLD)

    # vesicle facing fraction (from previous section)
    ves_vox = np.column_stack([
        np.clip((ves_nm[:, 2] / SCALE_Z).round().astype(int),
                0, axon_vol.shape[0] - 1),
        np.clip((ves_nm[:, 1] / SCALE_XY).round().astype(int),
                0, axon_vol.shape[1] - 1),
        np.clip((ves_nm[:, 0] / SCALE_XY).round().astype(int),
                0, axon_vol.shape[2] - 1),
    ])
    nearest_idx = res['nearest_idx']
    membrane_pts = np.zeros((len(ves_nm), 3))
    for i in range(len(ves_nm)):
        zi, yi, xi = ves_vox[i]
        membrane_pts[i] = [nearest_idx[2, zi, yi, xi] * SCALE_XY,
                           nearest_idx[1, zi, yi, xi] * SCALE_XY,
                           nearest_idx[0, zi, yi, xi] * SCALE_Z]
    proj_to_syn, _ = syn_tree.query(membrane_pts)
    vesicle_facing_frac = np.mean(proj_to_syn <= FACING_THRESHOLD)

    enrichment = vesicle_facing_frac / surface_facing_frac if surface_facing_frac > 0 else np.nan
    print(f'{name}: surface facing={surface_facing_frac*100:.1f}%, '
          f'vesicle facing={vesicle_facing_frac*100:.1f}%, '
          f'enrichment={enrichment:.2f}x')

print('\n=== polarity analysis complete ===')
