# -*- coding: utf-8 -*-
'''
Created on Fri Jun  9 18:50:11 2023
Modified on 26 Feb 2025

perform UMAP on ACG's from the general population to identify putative Dbh+
    neurones

@author: Dinghao Luo
'''

#%% imports
from pathlib import Path
import sys

import umap
import numpy as np
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import DBSCAN, KMeans

repo_root = Path(__file__).resolve().parents[2]
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))

from common_functions import normalise, mpl_formatting, gaussian_kernel_unity
mpl_formatting()

import project_paths as pp

import rec_list
paths = rec_list.pathLC


#%% parameters and path stems
ALL_SESS_STEM = pp.LC_EPHYS_STEM / 'all_sessions'
UMAP_STEM = pp.LC_EPHYS_STEM / 'UMAP'
KMEANS_PATH = UMAP_STEM / 'LC_all_UMAP_kmeans.npy'
KMEANS_FULL_PATH = UMAP_STEM / 'LC_all_UMAP_kmeans_full.npy'
DBSCAN_PATH = UMAP_STEM / 'LC_all_UMAP_dbscan.npy'
UMAP_PLOT_PAYLOAD_PATH = UMAP_STEM / 'LC_all_UMAP_plot_payload.npy'
DBSCAN_PARAMS = {
    'eps': 0.8,
    'min_samples': 8,
    }


#%% analysis
keys, identities, ACGs = [], [], []  # at this point there is only tagged and non-tagged in identities

for path in paths:
    recname = Path(path).name

    # identities
    curr_identities_path = ALL_SESS_STEM / recname / f'{recname}_all_identities.npy'
    curr_identities = np.load(curr_identities_path, allow_pickle=True).item()

    # ACGs
    curr_ACGs_path = ALL_SESS_STEM / recname / f'{recname}_all_ACGs.npy'
    curr_ACGs = np.load(curr_ACGs_path, allow_pickle=True).item()

    keys.extend(curr_identities.keys())
    identities.extend([int(i) for i in curr_identities.values()])
    ACGs.extend(curr_ACGs.values())

tagged_idx, non_tagged_idx = ([i for i, val in enumerate(identities) if val==1],
                              [i for i, val in enumerate(identities) if val==0])

# ----------------
# ACG UMAP starts
# ----------------
# create Gaussian kernel for smoothing
Gaussian = gaussian_kernel_unity(sigma=2)

# smooth the ACGs (only need one side)
ACG_arr = np.array(
    [normalise(np.convolve(ACG, Gaussian, mode='same')[9800:10000]) for ACG in ACGs]
    )

# initiate reducer
reducer = umap.UMAP(metric='cosine',
                    min_dist=0.0,
                    n_neighbors=30,
                    random_state=666)

# standard-scale the array first
scaled_ACG_arr = StandardScaler().fit_transform(ACG_arr)

# actual embedding
ACG_embedding = reducer.fit_transform(scaled_ACG_arr)
# --------------
# ACG UMAP ends
# --------------

# get coördinates of centres
tagged_med = np.median(ACG_embedding[tagged_idx,:], axis=0)

min_dist   = np.zeros(len(ACGs))
dist2mean  = np.zeros(len(ACGs))
for idx in non_tagged_idx:
    coörd = ACG_embedding[idx,:]

    mind = 100  # initialisation
    for j in tagged_idx:
        tg_coörd = ACG_embedding[j,:]
        new_mind = np.sqrt((coörd[0]-tg_coörd[0])**2 + (coörd[1]-tg_coörd[1])**2)
        if new_mind<mind:
            mind = new_mind

    min_dist[idx]  = mind
    dist2mean[idx] = np.sqrt((coörd[0]-tagged_med[0])**2 + (coörd[1]-tagged_med[1])**2)

xlower = min(ACG_embedding[:,0]) - 1.2
xupper = max(ACG_embedding[:,0]) + 1.2
ylower = min(ACG_embedding[:,1]) - .8
yupper = max(ACG_embedding[:,1]) + .8

# k-means
reducer_10 = umap.UMAP(metric='cosine',
                       min_dist=0.0,
                       n_neighbors=30,
                       n_components=2,  # identify 2 clusters
                       random_state=666)
ACG_embedding_10 = reducer_10.fit_transform(scaled_ACG_arr)

kmeans = KMeans(n_clusters=2, random_state=666, n_init='auto')
kmeans.fit(ACG_embedding_10)

labels = kmeans.labels_

# IMPORTANT: post-hoc relabelling based on tagged so the
#   cluster labels stay the same across runs (1 for tagged/putative and
#   0 for other)
if labels[tagged_idx[0]] == 0:  # if tagged (and thus putative) is 0
    labels = [1-x for x in labels]  # flip 1's and 0's
labels = np.asarray(labels, dtype=np.int32)

labels_dict = {keys[i]: labels[i] for i in range(len(keys))}

dbscan = DBSCAN(**DBSCAN_PARAMS)
dbscan_labels = dbscan.fit_predict(ACG_embedding).astype(np.int32)
dbscan_labels_dict = {
    keys[i]: dbscan_labels[i]
    for i in range(len(keys))
    }

UMAP_STEM.mkdir(parents=True, exist_ok=True)
np.save(KMEANS_PATH, labels_dict)
np.save(KMEANS_FULL_PATH, labels_dict)
np.save(DBSCAN_PATH, dbscan_labels_dict)
np.save(
    UMAP_PLOT_PAYLOAD_PATH,
    {
        'ACG_embedding': ACG_embedding,
        'tagged_idx': tagged_idx,
        'non_tagged_idx': non_tagged_idx,
        'tagged_med': tagged_med,
        'dist2mean': dist2mean,
        'labels': labels,
        'dbscan_labels': dbscan_labels,
        'dbscan_params': DBSCAN_PARAMS,
        'xlim': (xlower, xupper),
        'ylim': (ylower, yupper)
        }
    )
