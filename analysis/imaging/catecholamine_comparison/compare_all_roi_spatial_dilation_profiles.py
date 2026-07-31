# -*- coding: utf-8 -*-
'''
Created on Thu Jun 25 2026

compare dopamine and norepinephrine spatial spread around all detected axon ROIs

@author: Dinghao Luo
'''

#%% imports
from pathlib import Path
import csv
import sys

import numpy as np
from scipy.ndimage import binary_dilation
from scipy.optimize import curve_fit
from scipy.stats import ranksums, sem
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

repo_root = Path(__file__).resolve().parents[3]
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))

from common_functions import mpl_formatting
from console_formatting import print_files_saved, print_session, print_statistics_section, print_status
import imaging_utility_functions as iuf
import project_paths as pp
import rec_list
mpl_formatting()


#%% parameters
STRUCT      = np.ones((3, 3), dtype=bool)
DILATE_STEP = 1
MAX_DILATE  = 20

PIXEL_REDUCTION_THRESH = 0.5

WINDOW_BINS = {
    '0-1 s': range(0, 10),
    '1-2 s': range(10, 20),
    '2-3 s': range(20, 30),
    '3-4 s': range(30, 40),
    }

DATASETS = [
    {
        'transmitter': 'dopamine',
        'sensor': 'dLight',
        'paths': rec_list.pathdLightLCOpto,
        'stem': pp.HPC_DLIGHT_LC_OPTO_STEM,
        'colour': 'darkgreen',
        },
    {
        'transmitter': 'norepinephrine',
        'sensor': 'nLight',
        'paths': rec_list.pathnLightLCOpto,
        'stem': pp.HPC_NLIGHT_LC_OPTO_STEM,
        'colour': 'royalblue',
        },
    ]

save_stem = pp.FIGURES_ROOT / 'catecholamine_comparison'
data_stem = pp.DATA_ROOT / 'catecholamine_comparison'
save_stem.mkdir(parents=True, exist_ok=True)
data_stem.mkdir(parents=True, exist_ok=True)


#%% radial decay fits and plots
def _exp_decay(d, A, tau):
    return A * np.exp(-d / tau)

def plot_radial_summary(rows, normalise=False):
    fig, axs = plt.subplots(1, 4, figsize=(8.8, 2.2), sharex=True)

    for ax, wname in zip(axs, WINDOW_BINS.keys()):
        for dataset in DATASETS:
            curr = [
                row for row in rows
                if row['transmitter'] == dataset['transmitter'] and row['window'] == wname
                ]
            recnames = sorted(set(row['recname'] for row in curr))
            dilations = list(range(0, MAX_DILATE + 1, DILATE_STEP))
            mat = np.full((len(recnames), len(dilations)), np.nan)

            for i, recname in enumerate(recnames):
                rec_rows = [row for row in curr if row['recname'] == recname]
                by_dilation = {int(row['dilation_px']): float(row['mean_RI']) for row in rec_rows}
                for j, dilation in enumerate(dilations):
                    mat[i, j] = by_dilation[dilation]

                if normalise:
                    mat[i, :] = mat[i, :] / mat[i, 0]

            x = np.asarray(dilations, dtype=float)
            mean = np.nanmean(mat, axis=0)
            err = sem(mat, axis=0, nan_policy='omit')
            n_sess = int(np.sum(np.any(np.isfinite(mat), axis=1)))

            ax.plot(x, mean, color=dataset['colour'], linewidth=1.2,
                    label=f'{dataset["transmitter"]} (n={n_sess})')
            ax.fill_between(x, mean - err, mean + err,
                            color=dataset['colour'], alpha=0.18, edgecolor='none')

        ax.set(title=wname, xlabel='dilation (px)', xlim=(-1, 21))
        ax.spines[['top', 'right']].set_visible(False)

    ylabel = 'normalised RI' if normalise else 'mean RI'
    axs[0].set_ylabel(ylabel)
    axs[-1].legend(frameon=False, fontsize=6, loc='upper right')
    fig.tight_layout(w_pad=0.7)
    return fig


#%% main
curve_rows = []
tau_rows = []

for dataset in DATASETS:
    all_sess_stem = dataset['stem'] / 'all_sessions'

    for path in dataset['paths']:
        recname = Path(path).name
        print_session(f'{dataset["sensor"]}: {recname}')
        proc_path = all_sess_stem / recname / 'processed_data'

        ref_ch1_path = proc_path / 'ref_mat_ch1.npy'
        pixel_RI_path = proc_path / f'{recname}_pixel_RI_bins.npy'
        roi_dict_path = proc_path / f'{recname}_ROI_dict.npy'

        ref_ch1 = np.load(ref_ch1_path, allow_pickle=True)
        pixel_RI_bins = np.load(pixel_RI_path, allow_pickle=True)
        roi_dict = np.load(roi_dict_path, allow_pickle=True).item()

        thres_mask = iuf.generate_adaptive_membrane_mask(ref_ch1, visualize=0)
        ROI_mask_raw = iuf.build_roi_mask(roi_dict)
        ROI_mask = ROI_mask_raw & thres_mask

        n_raw = int(np.sum(ROI_mask_raw))
        n_thres = int(np.sum(ROI_mask))
        if n_raw == 0 or n_thres == 0:
            print_status('skipped', 'no ROI pixels after mask')
            continue

        pixel_reduction = 1 - (n_thres / n_raw)
        if pixel_reduction > PIXEL_REDUCTION_THRESH:
            print_status('skipped', f'pixel_reduction={pixel_reduction:.3f}')
            continue

        dilated_masks = {}
        ring_masks = {}
        session_curve_rows = []

        for dilation in range(0, MAX_DILATE + 1, DILATE_STEP):
            if dilation == 0:
                dilated_masks[0] = ROI_mask_raw
                ring_masks[0] = ROI_mask
            else:
                dm = binary_dilation(ROI_mask_raw, structure=STRUCT, iterations=dilation)
                dilated_masks[dilation] = dm
                ring_masks[dilation] = dm & ~dilated_masks[dilation - DILATE_STEP]
                ring_masks[dilation] = ring_masks[dilation] & ~ROI_mask_raw & thres_mask

            for wname, bins in WINDOW_BINS.items():
                vals = pixel_RI_bins[:, :, bins][ring_masks[dilation]]
                vals = vals[np.isfinite(vals)]
                mean_RI = float(np.mean(vals))
                n_pixels = int(vals.size)
                session_curve_rows.append({
                    'transmitter': dataset['transmitter'],
                    'sensor': dataset['sensor'],
                    'recname': recname,
                    'animal': recname.split('-')[0],
                    'window': wname,
                    'dilation_px': dilation,
                    'mean_RI': mean_RI,
                    'n_pixels': n_pixels,
                    'n_roi_total': len(roi_dict),
                    'pixel_reduction': float(pixel_reduction),
                    })

        session_tau_rows = []
        for wname in WINDOW_BINS.keys():
            curr = [row for row in session_curve_rows if row['window'] == wname]
            dilations = [row['dilation_px'] for row in curr]
            vals = [row['mean_RI'] for row in curr]
            dilations = np.asarray(dilations, dtype=float)
            vals = np.asarray(vals, dtype=float)

            if np.sum(np.isfinite(vals)) < 6:
                spatial_tau_px = np.nan
            else:
                centre = np.nanmin(vals)
                shifted = vals - centre
                valid = np.isfinite(shifted) & (shifted >= 0)
                if np.sum(valid) < 6 or np.nanmax(shifted[valid]) <= 0:
                    spatial_tau_px = np.nan
                else:
                    popt, _ = curve_fit(
                        _exp_decay,
                        dilations[valid],
                        shifted[valid],
                        bounds=([0, 0.5], [np.inf, 50]),
                        )
                    spatial_tau_px = float(popt[1])

            session_tau_rows.append({
                'transmitter': dataset['transmitter'],
                'sensor': dataset['sensor'],
                'recname': recname,
                'animal': recname.split('-')[0],
                'window': wname,
                'spatial_tau_px': spatial_tau_px,
                'n_roi_total': len(roi_dict),
                'pixel_reduction': float(pixel_reduction),
                })

        print_status('done', f'{len(roi_dict)} ROIs')
        curve_rows.extend(session_curve_rows)
        tau_rows.extend(session_tau_rows)


#%% save tables
curve_csv = data_stem / 'spatial_dilation_curves_all_rois.csv'
tau_csv = data_stem / 'spatial_tau_metrics_all_rois.csv'
with open(curve_csv, 'w', newline='', encoding='utf-8') as f:
    writer = csv.DictWriter(f, fieldnames=list(curve_rows[0].keys()))
    writer.writeheader()
    writer.writerows(curve_rows)

with open(tau_csv, 'w', newline='', encoding='utf-8') as f:
    writer = csv.DictWriter(f, fieldnames=list(tau_rows[0].keys()))
    writer.writeheader()
    writer.writerows(tau_rows)


#%% print summary
print_statistics_section()
for dataset in DATASETS:
    curr = [row for row in tau_rows if row['transmitter'] == dataset['transmitter']]
    print(f'{dataset["sensor"]}: {len(set(row["recname"] for row in curr))} sessions')
    for wname in WINDOW_BINS.keys():
        vals = np.asarray([
            float(row['spatial_tau_px']) for row in curr
            if row['window'] == wname
            ], dtype=float)
        vals = vals[np.isfinite(vals)]
        print(f'  {wname}: median tau={np.nanmedian(vals):.3f}, mean={np.nanmean(vals):.3f}')


#%% plot
saved = [('csv', curve_csv), ('csv', tau_csv)]

fig = plot_radial_summary(curve_rows, normalise=False)
for ext in ['.png', '.pdf']:
    out = save_stem / f'dopamine_norepinephrine_spatial_dilation_curves_all_rois{ext}'
    fig.savefig(out, dpi=300, bbox_inches='tight')
    saved.append((ext[1:], out))
plt.close(fig)

fig = plot_radial_summary(curve_rows, normalise=True)
for ext in ['.png', '.pdf']:
    out = save_stem / f'dopamine_norepinephrine_spatial_dilation_curves_all_rois_normalised{ext}'
    fig.savefig(out, dpi=300, bbox_inches='tight')
    saved.append((ext[1:], out))
plt.close(fig)

fig, axs = plt.subplots(1, 4, figsize=(8.8, 2.2), sharey=True)
rng = np.random.default_rng(241)

for ax, wname in zip(axs, WINDOW_BINS.keys()):
    groups = []
    labels = []
    colours = []

    for dataset in DATASETS:
        vals = [
            float(row['spatial_tau_px']) for row in tau_rows
            if row['transmitter'] == dataset['transmitter'] and row['window'] == wname
            ]
        vals = np.asarray(vals, dtype=float)
        vals = vals[np.isfinite(vals)]
        groups.append(vals)
        labels.append(dataset['transmitter'])
        colours.append(dataset['colour'])

    parts = ax.violinplot(groups, positions=[1, 2], showmedians=True, showextrema=False)
    for body, colour in zip(parts['bodies'], colours):
        body.set_facecolor(colour)
        body.set_edgecolor('none')
        body.set_alpha(0.6)
    parts['cmedians'].set_color('k')

    stat, p_val = ranksums(groups[0], groups[1])
    ax.text(0.05, 0.96, f'p={p_val:.2e}', transform=ax.transAxes,
            ha='left', va='top', fontsize=7)
    print(f'{wname}: tau ranksums z={stat:.3f}, p={p_val:.3e}')

    for xpos, vals, colour in zip([1, 2], groups, colours):
        jitter = rng.normal(0, 0.035, size=vals.size)
        ax.scatter(np.full(vals.size, xpos) + jitter, vals, s=12,
                   color=colour, edgecolor='none', alpha=0.65, zorder=3)

    ax.set(
        title=wname,
        xticks=[1, 2],
        xticklabels=labels,
        xlim=(0.5, 2.5),
        ylim=(0, 52),
        )
    ax.spines[['top', 'right']].set_visible(False)

axs[0].set_ylabel(r'spatial $\tau$ (px)')
fig.tight_layout(w_pad=0.7)
for ext in ['.png', '.pdf']:
    out = save_stem / f'dopamine_norepinephrine_spatial_tau_all_rois{ext}'
    fig.savefig(out, dpi=300, bbox_inches='tight')
    saved.append((ext[1:], out))
plt.close(fig)

print_files_saved(saved)
