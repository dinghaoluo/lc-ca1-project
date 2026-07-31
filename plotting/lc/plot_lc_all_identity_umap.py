# -*- coding: utf-8 -*-
'''
Created on 20 May 2026
Modified on 30 June 2026

plot LC ACG UMAP and clustering diagnostics

@author: Dinghao Luo
'''

#%% imports
from pathlib import Path
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib as mpl

repo_root = Path(__file__).resolve().parents[2]
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))

from common_functions import normalise, mpl_formatting
from common_functions import colour_putative, colour_tagged, colour_other
import project_paths as pp
mpl_formatting()


#%% paths
UMAP_STEM = pp.LC_EPHYS_STEM / 'UMAP'
UMAP_FIGURES_STEM = pp.LC_EPHYS_FIGURES_STEM / 'umap'
payload_path = UMAP_STEM / 'LC_all_UMAP_plot_payload.npy'


#%% save figures
def save_figure_variants(fig, output_stem, dpi=300):
    output_stem = Path(output_stem)
    output_stem.parent.mkdir(parents=True, exist_ok=True)
    for ext in ('.png', '.pdf'):
        fig.savefig(
            str(output_stem) + ext,
            dpi=dpi,
            bbox_inches='tight'
            )


#%% load data
payload = np.load(payload_path, allow_pickle=True).item()

#%% distance to tagged cells
ACG_embedding = payload['ACG_embedding']
tagged_idx = payload['tagged_idx']
tagged_med = payload['tagged_med']
dist2mean = payload['dist2mean']

dist2mean_norm = normalise(
    np.concatenate((dist2mean, np.array([max(dist2mean)*1.1])))
    )
dist2mean_cmap = mpl.colormaps['ocean'](dist2mean_norm)[:len(ACG_embedding)]

fig, ax = plt.subplots(figsize=(3,2.6))
umap_scatter_tagged = ax.scatter(
    ACG_embedding[tagged_idx, 0],
    ACG_embedding[tagged_idx, 1],
    s=10, color=colour_tagged, ec='none'
    )
ax.scatter(
    ACG_embedding[:, 0],
    ACG_embedding[:, 1],
    s=10, c=dist2mean_cmap, alpha=.75, ec='none'
    )
umap_scatter_tgcom = ax.scatter(
    tagged_med[0],
    tagged_med[1],
    s=20, color='darkblue'
    )
ax.legend([umap_scatter_tagged, umap_scatter_tgcom],
          ['Tgd. Dbh+', 'Tgd. Dbh+ CoM'],
          frameon=False, loc='upper left', fontsize=8)

ax.set(title='UMAP embedding of ACGs',
       xlabel='First dim.', ylabel='Second dim.',
       xticks=[], yticks=[],
       xlim=payload['xlim'], ylim=payload['ylim'])
fig.tight_layout()
save_figure_variants(fig, UMAP_FIGURES_STEM / 'LC_all_UMAP_distance_to_tagged_med')
plt.close(fig)

#%% k-means binary
ACG_embedding = payload['ACG_embedding']
labels = payload['labels']

colours = []
for i in labels:
    if i==1:
        colours.append('k')
    else:
        colours.append('grey')

fig, ax = plt.subplots(figsize=(3,2.6))
ax.scatter(
    ACG_embedding[:,0],
    ACG_embedding[:,1],
    c=colours, s=10, ec='none', alpha=.5
    )

c1 = ax.scatter([], [], c='k', s=10, ec='none')
c2 = ax.scatter([], [], c='grey', s=10, ec='none')
ax.legend([c1, c2], ['cluster 1', 'cluster 2'], frameon=False, fontsize=8)

ax.set(title='UMAP embedding of ACGs',
       xlabel='First dim.', ylabel='Second dim.',
       xticks=[], yticks=[],
       xlim=payload['xlim'], ylim=payload['ylim'])
fig.tight_layout()
save_figure_variants(fig, UMAP_FIGURES_STEM / 'LC_all_UMAP_kmeans_binary')
plt.close(fig)

#%% k-means categorised
ACG_embedding = payload['ACG_embedding']
tagged_idx = payload['tagged_idx']
labels = payload['labels']

colours = []
for i in labels:
    if i==1:
        colours.append(colour_putative)
    else:
        colours.append(colour_other)

fig, ax = plt.subplots(figsize=(2.4,2.15))

ax.scatter(ACG_embedding[:,0], ACG_embedding[:,1],
           c=colours, s=8, ec='k', linewidth=.5, alpha=1)
ax.scatter(ACG_embedding[tagged_idx, 0], ACG_embedding[tagged_idx, 1],
           s=8, color=colour_tagged, ec='k', linewidth=.5, alpha=1)

ntgcolor = ax.scatter([], [],
                      s=10, ec='k', color=colour_other, linewidth=.5, alpha=1)
tgcolor = ax.scatter([], [],
                     s=10, ec='k', color=colour_tagged, linewidth=.5, alpha=1)
ptcolor = ax.scatter([], [],
                     s=10, ec='k', color=colour_putative, linewidth=.5, alpha=1)

ax.legend([tgcolor, ptcolor, ntgcolor],
          ['Tagged Dbh+', 'Putative Dbh+', 'Putative Dbh-'],
          frameon=False, fontsize=8)

ax.set(title='UMAP embedding of ACGs',
       xlabel='First dim.', ylabel='Second dim.',
       xticks=[], yticks=[],
       xlim=payload['xlim'], ylim=payload['ylim'])
fig.tight_layout()
save_figure_variants(fig, UMAP_FIGURES_STEM / 'LC_all_UMAP_kmeans_categorised')
plt.close(fig)

#%% DBSCAN
ACG_embedding = payload['ACG_embedding']
labels = payload['dbscan_labels']

unique_labels = sorted(set(int(label) for label in labels))
cmap = mpl.colormaps['tab10']
colours = []
for label in labels:
    if int(label) == -1:
        colours.append('lightgrey')
    else:
        colours.append(cmap(int(label) % 10))

fig, ax = plt.subplots(figsize=(3,2.6))
ax.scatter(
    ACG_embedding[:,0],
    ACG_embedding[:,1],
    c=colours, s=10, ec='none', alpha=.65
    )

handles = []
labels_text = []
for label in unique_labels:
    if label == -1:
        handles.append(ax.scatter([], [], c='lightgrey', s=10, ec='none'))
        labels_text.append('noise')
    else:
        handles.append(ax.scatter([], [], c=[cmap(label % 10)], s=10, ec='none'))
        labels_text.append(f'cluster {label}')
ax.legend(handles, labels_text, frameon=False, fontsize=8)

ax.set(title='UMAP embedding of ACGs',
       xlabel='First dim.', ylabel='Second dim.',
       xticks=[], yticks=[],
       xlim=payload['xlim'], ylim=payload['ylim'])
fig.tight_layout()
save_figure_variants(fig, UMAP_FIGURES_STEM / 'LC_all_UMAP_dbscan')
plt.close(fig)
