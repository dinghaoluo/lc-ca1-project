# -*- coding: utf-8 -*-
"""
Created on Fri Sep 12 18:07:16 2025
Modified on 22 Jan 2026
Modified on 23 Mar 2026
    - recorded that the adaptive membrane-mask import was added with Jingyu Cao on 2026-03-23

analyse the dispersion of dLight signal after stim.
    dependent on extraction with hpc_dlight_lc_opto_extract.py
Modified to work on ROI vs neuropil over bins
Modified to use union of masks with dLight-expression-thresholded mask to
    extract traces

@author: Dinghao Luo
"""

#%% imports
from pathlib import Path
import sys

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import wilcoxon, sem
from scipy.optimize import curve_fit
from scipy.ndimage import binary_dilation

repo_root = Path(__file__).resolve().parent
while not (repo_root / 'utils').exists() and repo_root != repo_root.parent:
    repo_root = repo_root.parent
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))

from common_functions import mpl_formatting
import imaging_utility_functions as iuf
import project_paths as pp
mpl_formatting()


import rec_list
paths = rec_list.pathdLightLCOpto


#%% helper
def _exp_decay_fixed(t, A, tau, t0, B):
    return A * np.exp(-(t - t0) / tau) + B


#%% paths and parameters
dLight_stem   = pp.HPC_DLIGHT_LC_OPTO_STEM
all_sess_stem = dLight_stem / 'all_sessions'
save_stem      = pp.HPC_DLIGHT_LC_OPTO_FIGURES_STEM / 'ROI_vs_neuropil'
data_save_stem = pp.HPC_DLIGHT_LC_OPTO_STEM / 'ROI_vs_neuropil'

# how far away from ROI to count as neuropil
DISTANCE_FROM_ROI = 9  # 9 pixels ~ 5 um

ALPHA    = 0.05
MIN_RI   = 0.1
R2_THRES = 0.7
N_BINS   = 40
XAXIS    = np.arange(N_BINS) / 10


#%% main
all_ROI_RI_bins      = []
all_ROI_RI2_bins     = []
all_neuropil_RI_bins = []

all_ROI_RI_taus      = []
all_ROI_RI2_taus     = []
all_neuropil_RI_taus = []


# debug containers
debug_rows = []

def _safe_mean_trace(arr, mask):
    if np.sum(mask) == 0:
        return np.full(arr.shape[-1], np.nan)
    return np.nanmean(arr[mask, :], axis=0)

def _safe_std_trace(arr, mask):
    if np.sum(mask) == 0:
        return np.full(arr.shape[-1], np.nan)
    return np.nanstd(arr[mask, :], axis=0)

def _safe_sem_trace(arr, mask):
    n = np.sum(mask)
    if n == 0:
        return np.full(arr.shape[-1], np.nan)
    return np.nanstd(arr[mask, :], axis=0) / np.sqrt(n)


# main loop
for path in paths:
    recname = Path(path).name
    print(f'\n{recname}')

    ref_ch1_path       = all_sess_stem / recname / f'processed_data/{recname}_ref_mat_ch1.npy'
    pixel_RI_path      = all_sess_stem / recname / f'processed_data/{recname}_pixel_RI_bins.npy'
    pixel_RI_ch2_path  = all_sess_stem / recname / f'processed_data/{recname}_pixel_RI2_bins.npy'
    pixel_RI_stim_path = all_sess_stem / recname / f'processed_data/{recname}_pixel_RI_stim.npy'
    roi_path           = all_sess_stem / recname / f'processed_data/{recname}_ROI_dict.npy'

    if not ref_ch1_path.exists():
        print('No mean image for ch1; skipped')
        continue
    if not pixel_RI_path.exists():
        print('No pixel_RI array; skipped')
        continue
    if not pixel_RI_ch2_path.exists():
        print('No pixel_RI_ch2 array; skipped')
        continue
    if not pixel_RI_stim_path.exists():
        print('No pixel_RI_stim array; skipped')
        continue
    if not roi_path.exists():
        print('No roi_dict; skipped')
        continue

    ref_ch1           = np.load(ref_ch1_path, allow_pickle=True)
    pixel_RI_bins     = np.load(pixel_RI_path, allow_pickle=True)
    pixel_RI_ch2_bins = np.load(pixel_RI_ch2_path, allow_pickle=True)
    pixel_RI_stim     = np.load(pixel_RI_stim_path, allow_pickle=True)
    roi_dict          = np.load(roi_path, allow_pickle=True).item()

    H, W, _ = pixel_RI_bins.shape

    # ---- identify releasing ROIs ----
    releasing = iuf.identify_releasing_rois(pixel_RI_stim, roi_dict, alpha=ALPHA, min_ri=MIN_RI)

    if len(releasing) == 0:
        print('No releasing ROI; skipped')
        continue
    # ---- identification ends ----

    # new, dLight-expression-thresholded masks first, 23 March 2026
    thres_mask = iuf.generate_adaptive_membrane_mask(ref_ch1, visualize=0)

    # build masks without threshold first
    releasing_mask_raw = iuf.build_roi_mask(releasing)

    ROI_mask = iuf.build_roi_mask(roi_dict)
    ROI_dilated = binary_dilation(ROI_mask, iterations=DISTANCE_FROM_ROI)
    anti_ROI_mask_raw = ~ROI_dilated

    # thresholded versions
    releasing_mask_thres = releasing_mask_raw & thres_mask
    anti_ROI_mask_thres = anti_ROI_mask_raw & thres_mask

    # ------------------
    # DEBUG: compare masks
    # ------------------
    n_rel_raw = int(np.sum(releasing_mask_raw))
    n_rel_thres = int(np.sum(releasing_mask_thres))
    n_neu_raw = int(np.sum(anti_ROI_mask_raw))
    n_neu_thres = int(np.sum(anti_ROI_mask_thres))

    frac_rel_kept = n_rel_thres / n_rel_raw if n_rel_raw > 0 else np.nan
    frac_neu_kept = n_neu_thres / n_neu_raw if n_neu_raw > 0 else np.nan

    print(
        f'{recname} | ROI px raw={n_rel_raw}, thres={n_rel_thres}, kept={frac_rel_kept:.3f} | '
        f'neuropil px raw={n_neu_raw}, thres={n_neu_thres}, kept={frac_neu_kept:.3f}'
    )

    # ------------------
    # DEBUG: extract both versions
    # ------------------
    ROI_RI_bins_raw = _safe_mean_trace(pixel_RI_bins, releasing_mask_raw)
    ROI_RI2_bins_raw = _safe_mean_trace(pixel_RI_ch2_bins, releasing_mask_raw)
    neuropil_RI_bins_raw = _safe_mean_trace(pixel_RI_bins, anti_ROI_mask_raw)

    ROI_RI_bins = _safe_mean_trace(pixel_RI_bins, releasing_mask_thres)
    ROI_RI2_bins = _safe_mean_trace(pixel_RI_ch2_bins, releasing_mask_thres)
    neuropil_RI_bins = _safe_mean_trace(pixel_RI_bins, anti_ROI_mask_thres)

    # within-bin spread in ch2 across included pixels
    ROI_RI2_std_raw = _safe_std_trace(pixel_RI_ch2_bins, releasing_mask_raw)
    ROI_RI2_std_thres = _safe_std_trace(pixel_RI_ch2_bins, releasing_mask_thres)

    ROI_RI2_sem_raw = _safe_sem_trace(pixel_RI_ch2_bins, releasing_mask_raw)
    ROI_RI2_sem_thres = _safe_sem_trace(pixel_RI_ch2_bins, releasing_mask_thres)

    # how much the ch2 mean trace changed after thresholding
    ch2_diff = ROI_RI2_bins - ROI_RI2_bins_raw

    print(
        f'{recname} | ch2 mean abs diff={np.nanmean(np.abs(ch2_diff)):.5f} | '
        f'ch2 bin-std raw={np.nanmean(ROI_RI2_std_raw):.5f}, thres={np.nanmean(ROI_RI2_std_thres):.5f} | '
        f'ch2 bin-sem raw={np.nanmean(ROI_RI2_sem_raw):.5f}, thres={np.nanmean(ROI_RI2_sem_thres):.5f}'
    )

    debug_rows.append({
        'recname': recname,
        'n_rel_raw': n_rel_raw,
        'n_rel_thres': n_rel_thres,
        'frac_rel_kept': frac_rel_kept,
        'n_neu_raw': n_neu_raw,
        'n_neu_thres': n_neu_thres,
        'frac_neu_kept': frac_neu_kept,
        'ch2_mean_abs_diff': float(np.nanmean(np.abs(ch2_diff))),
        'ch2_bin_std_raw_mean': float(np.nanmean(ROI_RI2_std_raw)),
        'ch2_bin_std_thres_mean': float(np.nanmean(ROI_RI2_std_thres)),
        'ch2_bin_sem_raw_mean': float(np.nanmean(ROI_RI2_sem_raw)),
        'ch2_bin_sem_thres_mean': float(np.nanmean(ROI_RI2_sem_thres)),
        'ch2_trace_var_raw': float(np.nanvar(ROI_RI2_bins_raw)),
        'ch2_trace_var_thres': float(np.nanvar(ROI_RI2_bins)),
    })

    # build dilated mask and then the anti-mask
    ROI_mask      = iuf.build_roi_mask(roi_dict)
    ROI_dilated   = binary_dilation(ROI_mask, iterations=DISTANCE_FROM_ROI)
    anti_ROI_mask = ~ROI_dilated
    anti_ROI_mask = anti_ROI_mask & thres_mask  # new, intersection of anti_ROI_mask and thresholded, 23 Mar 2026

    # session plot
    fig_dbg, ax_dbg = plt.subplots(figsize=(4, 3))
    ax_dbg.plot(XAXIS, ROI_RI2_bins_raw, color='lightcoral', label='ch2 raw')
    ax_dbg.plot(XAXIS, ROI_RI2_bins, color='darkred', label='ch2 thresholded')
    ax_dbg.set_title(recname)
    ax_dbg.set_xlabel('Time from stim.-offset (s)')
    ax_dbg.set_ylabel('RI')
    ax_dbg.legend(frameon=False, fontsize=7)
    for s in ['top', 'right']:
        ax_dbg.spines[s].set_visible(False)
    fig_dbg.tight_layout()
    fig_dbg.savefig(save_stem / 'debug' / f'{recname}_debug_ch2_raw_vs_thres.png', dpi=200, bbox_inches='tight')
    plt.close(fig_dbg)


#%% debug summary
import pandas as pd

debug_df = pd.DataFrame(debug_rows)
data_save_stem.mkdir(parents=True, exist_ok=True)
debug_path = data_save_stem / 'debug_threshold_vs_raw.csv'
debug_df.to_csv(debug_path, index=False)

print('\n-----------------------\nthreshold debug summary\n-----------------------')
print(debug_df.describe(include='all'))

print('\nrecordings with biggest increase in ch2 trace variance after thresholding:')
tmp = debug_df.copy()
tmp['delta_ch2_trace_var'] = tmp['ch2_trace_var_thres'] - tmp['ch2_trace_var_raw']
print(tmp.sort_values('delta_ch2_trace_var', ascending=False)[
    ['recname', 'frac_rel_kept', 'n_rel_raw', 'n_rel_thres',
     'ch2_trace_var_raw', 'ch2_trace_var_thres', 'delta_ch2_trace_var',
     'ch2_bin_sem_raw_mean', 'ch2_bin_sem_thres_mean']
].head(10))
