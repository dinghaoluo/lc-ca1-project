# -*- coding: utf-8 -*-
'''
Created on Mon Jan 26 11:16:12 2026

Quantify the proportions of ROIs with significant release in each session

@author: Dinghao Luo
'''

#%% imports
from pathlib import Path
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import wilcoxon

repo_root = Path(__file__).resolve().parents[3]
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))

from common_functions import mpl_formatting
from console_formatting import print_files_saved, print_session, print_status
import project_paths as pp
mpl_formatting()

import rec_list

paths = rec_list.pathdLightLCOpto + \
        rec_list.pathdLightLCOptoDbhBlock


#%% parameters
ALPHA  = 0.05
MIN_RI = 0.1


#%% path stems
dLight_stem   = pp.HPC_DLIGHT_LC_OPTO_STEM
all_sess_stem = dLight_stem / 'all_sessions'


#%% main
all_release_proportions   = []  # proportions of axon ROIs that release in each session

for path in paths:
    recname = Path(path).name
    print_session(recname)

    pixel_RI_stim_path  = all_sess_stem / recname / 'processed_data' / f'{recname}_pixel_RI_stim.npy'
    roi_dict_path       = all_sess_stem / recname / 'processed_data' / f'{recname}_ROI_dict.npy'

    # load data
    print_status('loading data')
    pixel_RI_stim = np.load(pixel_RI_stim_path, allow_pickle=True)  # (512,512,40)
    roi_dict      = np.load(roi_dict_path, allow_pickle=True).item()

    # ---- identify releasing ROIs ----
    releasing_rois = {}

    for rid, roi in roi_dict.items():
        vals  = pixel_RI_stim[roi['ypix'], roi['xpix'], :]
        means = np.nanmean(vals, axis=0)  # mean over pixels
        means = [mean for mean in means if np.isfinite(mean)]  # filtering first
        if len(means) > 2:
            _, p = wilcoxon(means, alternative='greater')
            if p < ALPHA and np.mean(means) > MIN_RI:
                releasing_rois[rid] = roi

    if len(releasing_rois) == 0:
        print_status('skipped', 'no releasing ROI')
        continue
    # ---- identification ends ----

    # get proportion of releasing ROIs
    curr_prop = len(releasing_rois) / len(roi_dict)
    all_release_proportions.append(curr_prop)


#%% plotting
props = np.array(all_release_proportions)

# bins in proportion space [0, 1]
bins = np.linspace(0, 1, 21)

fig, ax = plt.subplots(figsize=(3, 2.4))

ax.hist(
    props,
    bins=bins,
    color='darkgreen',
    alpha=0.6,
    edgecolor='none'
)

# median + IQR
q1, med, q3 = np.percentile(props, [25, 50, 75])

ax.axvline(med, color='darkgreen', linestyle='--', lw=1)

ax.text(
    0.02, 0.95,
    f'median = {med:.4f}\nIQR = [{q1:.4f}, {q3:.4f}]',
    transform=ax.transAxes,
    va='top', ha='left',
    fontsize=7, color='darkgreen'
)

ax.set(
    xlabel='proportion of releasing ROIs per session',
    ylabel='session count',
    title='session-wise release proportion'
)

for s in ['top', 'right']:
    ax.spines[s].set_visible(False)

fig.tight_layout()

saved_paths = []
for ext in ['.png', '.pdf']:
    curr_path = pp.HPC_DLIGHT_LC_OPTO_FIGURES_STEM / f'session_release_proportion_hist{ext}'
    fig.savefig(
        curr_path,
        dpi=300,
        bbox_inches='tight'
    )
    saved_paths.append((ext[1:], curr_path))
print_files_saved(saved_paths)
