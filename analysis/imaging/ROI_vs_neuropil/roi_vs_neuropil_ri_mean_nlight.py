'''
Created on Mon Jun 22 2026

compare nLight ROI and neuropil response indices

@author: Dinghao Luo
'''

#%% imports
from pathlib import Path
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.ndimage import binary_dilation

repo_root = Path(__file__).resolve().parents[3]
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))

from common_functions import mpl_formatting
from console_formatting import print_session, print_status
import imaging_utility_functions as iuf
import project_paths as pp
from plotting_functions import plot_violin_with_scatter
mpl_formatting()

import rec_list
paths = rec_list.pathnLightLCOpto


#%% path stems
all_sess_stem = pp.HPC_NLIGHT_LC_OPTO_STEM / 'all_sessions'
save_stem     = pp.HPC_NLIGHT_LC_OPTO_FIGURES_STEM / 'ROI_vs_neuropil'


#%% parameters
# how far away from ROI to count as neuropil
DISTANCE_FROM_ROI = 9  # 9 pixels ~ 5 um

ALPHA  = 0.05
MIN_RI = 0.1


#%% analysis
# initialise containers
all_ROI_RIs      = []
all_ROI_RI2s     = []
all_neuropil_RIs = []
save_stem.mkdir(parents=True, exist_ok=True)

# loop
for path in paths:
    recname = Path(path).name
    print_session(recname)

    pixel_RI_stim_path  = all_sess_stem / recname / f'processed_data/{recname}_pixel_RI_stim.npy'
    pixel_RI2_stim_path = all_sess_stem / recname / f'processed_data/{recname}_pixel_RI_ch2_stim.npy'
    roi_dict_path       = all_sess_stem / recname / f'processed_data/{recname}_ROI_dict.npy'

    # load data
    print_status('loading data')
    pixel_RI_stim  = np.load(pixel_RI_stim_path, allow_pickle=True)
    pixel_RI2_stim = np.load(pixel_RI2_stim_path, allow_pickle=True)
    roi_dict       = np.load(roi_dict_path, allow_pickle=True).item()

    # ---- identify releasing ROIs ----
    releasing = iuf.identify_releasing_rois(pixel_RI_stim, roi_dict, alpha=ALPHA, min_ri=MIN_RI)

    if not releasing:
        print_status('skipped', 'no releasing ROI')
        continue
    # ---- identification ends ----

    # build mask of all ROIs
    releasing_mask = iuf.build_roi_mask(releasing)

    # build dilated mask and then the anti-mask
    ROI_mask = iuf.build_roi_mask(roi_dict)
    ROI_dilated = binary_dilation(ROI_mask, iterations=DISTANCE_FROM_ROI)
    anti_ROI_mask = ~ROI_dilated

    # get med of pixel_RI_stim
    pixel_RI_stim_med  = np.nanmedian(pixel_RI_stim, axis=2)
    pixel_RI2_stim_med = np.nanmedian(pixel_RI2_stim, axis=2)

    # get ROI and neuropil RI medians
    ROI_RI_med      = np.nanmedian(pixel_RI_stim_med[releasing_mask])
    ROI_RI2_med     = np.nanmedian(pixel_RI2_stim_med[releasing_mask])
    neuropil_RI_med = np.nanmedian(pixel_RI_stim_med[anti_ROI_mask])

    # append
    all_ROI_RIs.append(ROI_RI_med)
    all_ROI_RI2s.append(ROI_RI2_med)
    all_neuropil_RIs.append(neuropil_RI_med)

    # flat grey background + neuropil + ROI overlay
    fig, ax = plt.subplots(figsize=(5, 5))

    # neuropil overlay (slightly darker grey)
    neuropil_overlay = np.zeros((*anti_ROI_mask.shape, 4))
    neuropil_overlay[..., :3] = 0.6
    neuropil_overlay[..., 3] = anti_ROI_mask.astype(float)

    ax.imshow(neuropil_overlay, interpolation='nearest')

    # ROI overlay (dark green)
    roi_overlay = np.zeros((*releasing_mask.shape, 4))
    roi_overlay[..., 1] = 0.35
    roi_overlay[..., 3] = releasing_mask.astype(float)

    ax.imshow(roi_overlay, interpolation='nearest')

    ax.axis('off')
    ax.set_title(f'{recname}\nROI (dark green) on neuropil (grey)')
    plt.tight_layout()

    for ext in ['.png', '.pdf']:
        fig.savefig(save_stem / f'{recname}{ext}',
                    dpi=300,
                    bbox_inches='tight')

    plt.close(fig)


#%% statistics
plot_violin_with_scatter(all_neuropil_RIs, all_ROI_RIs,
                         'grey', 'darkgreen',
                         xticklabels=['Neuropil', 'ROI'],
                         ylabel='nLight RI',
                         save=True,
                         print_statistics=True,
                         savepath=save_stem / 'neuropil_vs_ROI_violin')

plot_violin_with_scatter(all_ROI_RI2s, all_ROI_RIs,
                         'darkred', 'darkgreen',
                         xticklabels=['ROI (ctrl.)', 'ROI'],
                         ylabel='RI',
                         save=True,
                         print_statistics=True,
                         savepath=save_stem / 'ROIctrl_vs_ROI_violin')
