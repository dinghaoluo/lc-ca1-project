# -*- coding: utf-8 -*-
'''
Created on Tue Aug 12 2026

vesicle-to-synapse spatial analysis: two complementary tests
    1. proximity test — are vesicles closer to nearby synapses than CSR?
    2. loading variance test — do vesicles preferentially target specific
       synapses (uneven loading) vs uniform distribution?

only synapses within 1000nm of the axon surface are included; distant
    neuropil synapses annotated in the full EM field are excluded

@author: Dinghao Luo
'''

#%% imports
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import tifffile
import matplotlib.pyplot as plt
from scipy.ndimage import binary_erosion
from scipy.spatial import KDTree
from scipy.stats import ks_2samp, mannwhitneyu

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
N_SIMULATIONS = 200
SYNAPSE_THRESHOLD = 1000  # nm — only include synapses within this of axon surface

plt.rcParams['font.family'] = 'Arial'


#%% helper
def sample_random_points(cyto_voxels, n_points, rng):
    chosen = rng.choice(len(cyto_voxels), size=n_points, replace=True)
    voxels = cyto_voxels[chosen].astype(np.float64)
    voxels += rng.uniform(-0.5, 0.5, size=voxels.shape)
    return np.column_stack([
        voxels[:, 2] * SCALE_XY,
        voxels[:, 1] * SCALE_XY,
        voxels[:, 0] * SCALE_Z,
    ])


#%% load and filter
rng = np.random.default_rng(42)
results = {}

for name, info in PROFILES.items():
    profile_dir = info['dir']

    axon_files = sorted(profile_dir.glob(info['axon_pattern']))
    axon_vol = (np.stack([tifffile.imread(f) for f in axon_files]) > 0
                ).astype(np.uint8)
    mito_files = sorted(profile_dir.glob(info['mito_pattern']))
    mito_vol = (np.stack([tifffile.imread(f) for f in mito_files]) > 0
                ).astype(np.uint8)
    cyto_vol = axon_vol & (~mito_vol)
    cyto_voxels = np.argwhere(cyto_vol > 0)

    # axon surface via erosion
    interior = binary_erosion(axon_vol)
    surface_mask = axon_vol & (~interior)
    surface_vox = np.argwhere(surface_mask > 0)
    surface_nm = np.column_stack([
        surface_vox[:, 2] * SCALE_XY,
        surface_vox[:, 1] * SCALE_XY,
        surface_vox[:, 0] * SCALE_Z,
    ])
    surface_tree = KDTree(surface_nm)

    # vesicles
    ves_df = pd.read_csv(profile_dir / info['ves_csv'], sep=info['ves_sep'])
    ves_nm = ves_df[['PosX', 'PosY', 'PosZ']].values
    if not info['ves_in_nm']:
        ves_nm = ves_nm * 1000

    # synapses — filter to those within threshold of axon surface
    syn_df = pd.read_csv(info['syn_csv'], sep=info['syn_sep'])
    syn_all = syn_df[['PosX', 'PosY', 'PosZ']].values * 1000
    syn_to_surface, _ = surface_tree.query(syn_all)
    keep = syn_to_surface <= SYNAPSE_THRESHOLD
    syn_nm = syn_all[keep]

    results[name] = {
        'cyto_voxels': cyto_voxels,
        'ves_nm': ves_nm,
        'syn_nm': syn_nm,
        'n_total': len(syn_all),
        'n_kept': int(keep.sum()),
        'syn_dists': syn_to_surface,
    }
    print(f'{name}: kept {keep.sum()}/{len(syn_all)} synapses '
          f'within {SYNAPSE_THRESHOLD} nm of axon surface')


#%% compute vesicle-to-synapse distances and CSR
proximity_results = {}

for name, res in results.items():
    ves_nm = res['ves_nm']
    syn_nm = res['syn_nm']
    cyto_voxels = res['cyto_voxels']

    if len(syn_nm) < 2:
        print(f'{name}: too few nearby synapses ({len(syn_nm)}), skipping')
        proximity_results[name] = None
        continue

    syn_tree = KDTree(syn_nm)
    ves_to_syn, _ = syn_tree.query(ves_nm)

    csr_dists = np.zeros(N_SIMULATIONS * len(ves_nm))
    for sim in range(N_SIMULATIONS):
        rand_pts = sample_random_points(cyto_voxels, len(ves_nm), rng)
        d, _ = syn_tree.query(rand_pts)
        csr_dists[sim * len(ves_nm):(sim + 1) * len(ves_nm)] = d

    _, mwu_p = mannwhitneyu(ves_to_syn, csr_dists, alternative='less')
    ks_stat, ks_p = ks_2samp(ves_to_syn, csr_dists)

    proximity_results[name] = {
        'ves_to_syn': ves_to_syn,
        'csr_dists': csr_dists,
        'mwu_p': mwu_p,
        'ks_stat': ks_stat,
        'ks_p': ks_p,
    }

    pstr = 'p < 0.001' if mwu_p < 0.001 else f'p = {mwu_p:.3f}'
    print(f'{name} ({len(syn_nm)} syn): '
          f'median obs = {np.median(ves_to_syn):.0f} nm, '
          f'CSR = {np.median(csr_dists):.0f} nm, '
          f'MWU {pstr}, KS = {ks_stat:.3f}')


#%% histogram: vesicle-to-synapse distance vs CSR
fig, axes = plt.subplots(2, 2, figsize=(6, 5), constrained_layout=True)
axes = axes.ravel()

for idx, (name, res) in enumerate(results.items()):
    ax = axes[idx]
    pr = proximity_results[name]

    if pr is None:
        ax.text(0.5, 0.5, 'too few synapses', transform=ax.transAxes,
                ha='center', va='center', fontsize=8)
        ax.set_title(f'{name}', fontsize=8)
        continue

    ves_to_syn = pr['ves_to_syn']
    csr_dists = pr['csr_dists']
    mwu_p = pr['mwu_p']

    upper = max(np.percentile(ves_to_syn, 99), np.percentile(csr_dists, 95))
    bins = np.linspace(0, upper, 30)
    ax.hist(ves_to_syn, bins=bins, density=True, color='grey',
            edgecolor='white', linewidth=0.5, alpha=0.8, label='observed')
    ax.hist(csr_dists, bins=bins, density=True, color='royalblue',
            alpha=0.25, edgecolor='none', label='CSR')
    ax.axvline(np.median(ves_to_syn), color='black', linestyle='--',
               linewidth=0.8)
    ax.axvline(np.median(csr_dists), color='royalblue', linestyle='--',
               linewidth=0.8, alpha=0.6)

    pstr = 'p < 0.001' if mwu_p < 0.001 else f'p = {mwu_p:.3f}'
    n_syn = res['n_kept']
    ax.text(0.95, 0.85, f'n_syn = {n_syn}/{res["n_total"]}\n'
            f'median obs: {np.median(ves_to_syn):.0f} nm\n'
            f'median CSR: {np.median(csr_dists):.0f} nm\n'
            f'MWU {pstr}',
            transform=ax.transAxes, ha='right', fontsize=5.5,
            family='monospace', va='top')

    ax.set(xlabel='distance to nearest synapse (nm)', ylabel='density')
    ax.set_title(f'{name} (filtered)', fontsize=8)
    ax.spines[['top', 'right']].set_visible(False)
    ax.legend(fontsize=6, frameon=False)

for ext in ['.pdf', '.png']:
    fig.savefig(FIGURES_STEM / f'vesicle_to_synapse_distance_filtered{ext}',
                dpi=300, bbox_inches='tight')
plt.close(fig)


#%% CDF comparison
fig, axes = plt.subplots(2, 2, figsize=(6, 5), constrained_layout=True)
axes = axes.ravel()

for idx, (name, res) in enumerate(results.items()):
    ax = axes[idx]
    pr = proximity_results[name]

    if pr is None:
        ax.text(0.5, 0.5, 'too few synapses', transform=ax.transAxes,
                ha='center', va='center', fontsize=8)
        ax.set_title(f'{name}', fontsize=8)
        continue

    ves_to_syn = pr['ves_to_syn']
    csr_dists = pr['csr_dists']

    ves_sorted = np.sort(ves_to_syn)
    csr_sorted = np.sort(csr_dists)
    ax.plot(ves_sorted, np.linspace(0, 1, len(ves_sorted)),
            color='black', linewidth=1.2, label='observed')
    ax.plot(csr_sorted, np.linspace(0, 1, len(csr_sorted)),
            color='royalblue', linewidth=1, alpha=0.7, label='CSR')

    ks_stat = pr['ks_stat']
    ks_p = pr['ks_p']
    pstr = 'p < 0.001' if ks_p < 0.001 else f'p = {ks_p:.3f}'
    ax.text(0.95, 0.15, f'KS = {ks_stat:.3f}\n{pstr}',
            transform=ax.transAxes, ha='right', fontsize=6,
            family='monospace')

    n_syn = res['n_kept']
    ax.set(xlabel='distance to nearest synapse (nm)', ylabel='CDF')
    ax.set_title(f'{name} ({n_syn} syn, filtered)', fontsize=8)
    ax.spines[['top', 'right']].set_visible(False)
    ax.legend(fontsize=6, frameon=False, loc='lower right')

for ext in ['.pdf', '.png']:
    fig.savefig(FIGURES_STEM / f'vesicle_to_synapse_cdf_filtered{ext}',
                dpi=300, bbox_inches='tight')
plt.close(fig)


#%% shell enrichment with filtered synapses
SHELL_EDGES = np.arange(0, 1001, 100)

fig, axes = plt.subplots(2, 2, figsize=(7, 5.5), constrained_layout=True)
axes = axes.ravel()

for idx, (name, res) in enumerate(results.items()):
    ax = axes[idx]
    pr = proximity_results[name]

    if pr is None:
        ax.text(0.5, 0.5, 'too few synapses', transform=ax.transAxes,
                ha='center', va='center', fontsize=8)
        ax.set_title(f'{name}', fontsize=8)
        continue

    ves_nm = res['ves_nm']
    syn_nm = res['syn_nm']
    cyto_voxels = res['cyto_voxels']
    ves_to_syn = pr['ves_to_syn']

    syn_tree = KDTree(syn_nm)

    # subsample cytoplasm for volume estimation
    rng_sub = np.random.default_rng(77)
    n_sub = min(100000, len(cyto_voxels))
    sub_idx = rng_sub.choice(len(cyto_voxels), n_sub, replace=False)
    cyto_sub_nm = np.column_stack([
        cyto_voxels[sub_idx, 2] * SCALE_XY,
        cyto_voxels[sub_idx, 1] * SCALE_XY,
        cyto_voxels[sub_idx, 0] * SCALE_Z,
    ])
    cyto_to_syn, _ = syn_tree.query(cyto_sub_nm)

    n_ves = len(ves_nm)
    shell_ves_frac = []
    shell_vol_frac = []
    for j in range(len(SHELL_EDGES) - 1):
        lo, hi = SHELL_EDGES[j], SHELL_EDGES[j + 1]
        shell_ves_frac.append(
            np.sum((ves_to_syn >= lo) & (ves_to_syn < hi)) / n_ves * 100)
        shell_vol_frac.append(
            np.sum((cyto_to_syn >= lo) & (cyto_to_syn < hi)) / n_sub * 100)

    x = np.arange(len(shell_ves_frac))
    w = 0.35
    ax.bar(x - w/2, shell_vol_frac, w, color='royalblue', alpha=0.5,
           label='% cytoplasm vol')
    ax.bar(x + w/2, shell_ves_frac, w, color='gold', edgecolor='darkorange',
           linewidth=0.4, label='% vesicles')
    labels = [f'{SHELL_EDGES[j]}-{SHELL_EDGES[j+1]}'
              for j in range(len(SHELL_EDGES) - 1)]
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=5, rotation=45, ha='right')
    ax.set_ylabel('% of total', fontsize=7)
    ax.set_title(f'{name} ({res["n_kept"]} syn, filtered)', fontsize=8)
    ax.spines[['top', 'right']].set_visible(False)
    ax.legend(fontsize=6, frameon=False)

for ext in ['.pdf', '.png']:
    fig.savefig(FIGURES_STEM / f'synapse_shell_enrichment_filtered{ext}',
                dpi=300, bbox_inches='tight')
plt.close(fig)


#%% per-synapse loading variance test
# the proximity test saturates in dense-synapse geometries (P3/P4) where any
# random cytoplasm point is already near a synapse; the complementary question:
# do vesicles preferentially load *specific* synapses (uneven distribution)?
N_SIMS_LOADING = 500

fig, axes = plt.subplots(2, 2, figsize=(6, 5), constrained_layout=True)
axes = axes.ravel()

loading_results = {}
for idx, (name, res) in enumerate(results.items()):
    ax = axes[idx]
    syn_nm = res['syn_nm']
    ves_nm = res['ves_nm']
    cyto_voxels = res['cyto_voxels']

    if len(syn_nm) < 3:
        ax.text(0.5, 0.5, 'too few synapses', transform=ax.transAxes,
                ha='center', va='center', fontsize=8)
        ax.set_title(f'{name}', fontsize=8)
        loading_results[name] = None
        continue

    n_syn = len(syn_nm)
    n_ves = len(ves_nm)
    syn_tree = KDTree(syn_nm)

    _, ves_assignments = syn_tree.query(ves_nm)
    obs_counts = np.bincount(ves_assignments, minlength=n_syn)
    obs_cv = np.std(obs_counts) / np.mean(obs_counts)

    csr_cvs = np.zeros(N_SIMS_LOADING)
    for sim in range(N_SIMS_LOADING):
        rand_pts = sample_random_points(cyto_voxels, n_ves, rng)
        _, rand_assign = syn_tree.query(rand_pts)
        rand_counts = np.bincount(rand_assign, minlength=n_syn)
        csr_cvs[sim] = np.std(rand_counts) / np.mean(rand_counts)

    p_val = np.mean(csr_cvs >= obs_cv)
    loading_results[name] = {
        'obs_cv': obs_cv, 'csr_cvs': csr_cvs, 'p': p_val,
        'obs_counts': obs_counts,
    }

    pstr = 'p < 0.002' if p_val < 1/N_SIMS_LOADING else f'p = {p_val:.3f}'
    print(f'{name} ({n_syn} syn, {n_ves} ves): '
          f'CV = {obs_cv:.3f} vs CSR median = {np.median(csr_cvs):.3f}, '
          f'{pstr}')

    ax.hist(csr_cvs, bins=25, density=True, color='royalblue', alpha=0.4,
            edgecolor='none', label='CSR null')
    ax.axvline(obs_cv, color='red', linewidth=1.5, linestyle='-',
               label=f'observed (CV={obs_cv:.2f})')
    ax.axvline(np.median(csr_cvs), color='royalblue', linewidth=1,
               linestyle='--', alpha=0.7)
    ax.text(0.95, 0.85, f'n_syn={n_syn}, n_ves={n_ves}\n{pstr}',
            transform=ax.transAxes, ha='right', va='top', fontsize=6,
            family='monospace')
    ax.set(xlabel='CV of per-synapse vesicle count', ylabel='density')
    ax.set_title(f'{name}', fontsize=8)
    ax.spines[['top', 'right']].set_visible(False)
    ax.legend(fontsize=6, frameon=False)

for ext in ['.pdf', '.png']:
    fig.savefig(FIGURES_STEM / f'per_synapse_loading_variance{ext}',
                dpi=300, bbox_inches='tight')
plt.close(fig)


#%% summary
print('\n=== SUMMARY ===')
print('--- proximity test (vesicles closer to synapses than CSR?) ---')
for name in results:
    pr = proximity_results[name]
    res = results[name]
    if pr is None:
        print(f'{name}: SKIPPED (too few synapses within {SYNAPSE_THRESHOLD} nm)')
        continue
    pstr = 'p < 0.001' if pr['mwu_p'] < 0.001 else f'p = {pr["mwu_p"]:.4f}'
    print(f'  {name}: {res["n_kept"]}/{res["n_total"]} syn, '
          f'median obs={np.median(pr["ves_to_syn"]):.0f} vs '
          f'CSR={np.median(pr["csr_dists"]):.0f} nm, MWU {pstr}')

print('\n--- loading variance test (vesicles target specific synapses?) ---')
for name in results:
    lr = loading_results[name]
    if lr is None:
        print(f'  {name}: SKIPPED')
        continue
    pstr = 'p < 0.002' if lr['p'] < 1/N_SIMS_LOADING else f'p = {lr["p"]:.3f}'
    print(f'  {name}: CV={lr["obs_cv"]:.3f} vs CSR={np.median(lr["csr_cvs"]):.3f}, '
          f'{pstr}')

print('\n=== analysis complete ===')
