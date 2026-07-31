# -*- coding: utf-8 -*-
'''
Created on Thu Jun 25 2026

compare dopamine and norepinephrine roi-versus-neuropil time courses

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
DISTANCE_FROM_ROI = 9

ALPHA  = 0.05
MIN_RI = 0.1
R2_THRES = 0.7
PIXEL_REDUCTION_THRESH = 0.5

N_BINS = 40
XAXIS = np.arange(N_BINS) / 10

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


#%% decay fits and plots
def _exp_decay_fixed(t, A, tau, t0, B):
    return A * np.exp(-(t - t0) / tau) + B

def fit_tau(trace):
    first_valid = np.where(np.isfinite(trace))[0][0]
    t0 = XAXIS[first_valid]
    B = np.nanmin(trace[first_valid:])
    t_fit = XAXIS[first_valid:]
    y_fit = trace[first_valid:]

    popt, _ = curve_fit(
        lambda t, A, tau: _exp_decay_fixed(t, A, tau, t0, B),
        t_fit,
        y_fit,
        p0=[y_fit[0] - B, 1.0],
        bounds=([0, 0], [np.inf, np.inf]),
        )

    A_fit, tau_fit = popt
    y_pred = _exp_decay_fixed(t_fit, A_fit, tau_fit, t0, B)
    ss_res = np.sum((y_fit - y_pred) ** 2)
    ss_tot = np.sum((y_fit - np.mean(y_fit)) ** 2)
    r2 = 1 - ss_res / ss_tot

    if r2 < R2_THRES:
        return np.nan, float(r2)

    return float(tau_fit), float(r2)

def plot_trace_panel(ax, rows, dataset):
    traces_roi = np.asarray([row['roi_trace'] for row in rows if row['transmitter'] == dataset['transmitter']], dtype=float)
    traces_roi2 = np.asarray([row['roi_ctrl_trace'] for row in rows if row['transmitter'] == dataset['transmitter']], dtype=float)
    traces_neu = np.asarray([row['neuropil_trace'] for row in rows if row['transmitter'] == dataset['transmitter']], dtype=float)

    mean_roi = np.nanmean(traces_roi, axis=0)
    err_roi = sem(traces_roi, axis=0, nan_policy='omit')
    mean_roi2 = np.nanmean(traces_roi2, axis=0)
    err_roi2 = sem(traces_roi2, axis=0, nan_policy='omit')
    mean_neu = np.nanmean(traces_neu, axis=0)
    err_neu = sem(traces_neu, axis=0, nan_policy='omit')

    ax.plot(XAXIS, mean_roi, color='darkgreen', linewidth=1.4, label='ROI')
    ax.fill_between(XAXIS, mean_roi - err_roi, mean_roi + err_roi,
                    color='darkgreen', alpha=0.20, edgecolor='none')

    ax.plot(XAXIS, mean_roi2, color='darkred', linewidth=1.0, label='ROI (ctrl.)')
    ax.fill_between(XAXIS, mean_roi2 - err_roi2, mean_roi2 + err_roi2,
                    color='darkred', alpha=0.16, edgecolor='none')

    ax.plot(XAXIS, mean_neu, color='grey', linewidth=1.2, label='Neuropil')
    ax.fill_between(XAXIS, mean_neu - err_neu, mean_neu + err_neu,
                    color='grey', alpha=0.18, edgecolor='none')

    ax.set(title=f'{dataset["sensor"]}  n={traces_roi.shape[0]}',
           xlabel='time from stim.-offset (s)',
           ylabel='RI')
    ax.legend(frameon=False, fontsize=6)
    ax.spines[['top', 'right']].set_visible(False)

def plot_tau_violin(ax, rows, key):
    groups = []
    labels = []
    colours = []
    for dataset in DATASETS:
        vals = [
            float(row[key]) for row in rows
            if row['transmitter'] == dataset['transmitter'] and np.isfinite(float(row[key]))
            ]
        vals = np.asarray(vals, dtype=float)
        groups.append(vals)
        labels.append(f'{dataset["transmitter"]}\nn={vals.size}')
        colours.append(dataset['colour'])

    parts = ax.violinplot(groups, positions=[1, 2], showmedians=True, showextrema=False)
    for body, colour in zip(parts['bodies'], colours):
        body.set_facecolor(colour)
        body.set_edgecolor('none')
        body.set_alpha(0.65)
    parts['cmedians'].set_color('k')

    rng = np.random.default_rng(137)
    for xpos, vals, colour in zip([1, 2], groups, colours):
        jitter = rng.normal(0, 0.035, size=vals.size)
        ax.scatter(np.full(vals.size, xpos) + jitter, vals, s=12, color=colour,
                   edgecolor='none', alpha=0.65, zorder=3)

    stat, p_val = ranksums(groups[0], groups[1])
    ax.text(0.05, 0.96, f'p={p_val:.2e}', transform=ax.transAxes,
            ha='left', va='top', fontsize=7)
    print(f'{key}: ranksums z={stat:.3f}, p={p_val:.3e}')

    ax.set(xticks=[1, 2], xticklabels=labels, xlim=(0.5, 2.5))
    ax.spines[['top', 'right']].set_visible(False)


#%% main
all_rows = []

for dataset in DATASETS:
    all_sess_stem = dataset['stem'] / 'all_sessions'

    for path in dataset['paths']:
        recname = Path(path).name
        print_session(f'{dataset["sensor"]}: {recname}')
        proc_path = all_sess_stem / recname / 'processed_data'

        ref_ch1_path = proc_path / 'ref_mat_ch1.npy'
        pixel_RI_path = proc_path / f'{recname}_pixel_RI_bins.npy'
        pixel_RI2_path = proc_path / f'{recname}_pixel_RI2_bins.npy'
        pixel_RI_stim_path = proc_path / f'{recname}_pixel_RI_stim.npy'
        roi_path = proc_path / f'{recname}_ROI_dict.npy'

        pixel_RI_bins = np.load(pixel_RI_path, allow_pickle=True)
        pixel_RI2_bins = np.load(pixel_RI2_path, allow_pickle=True)
        pixel_RI_stim = np.load(pixel_RI_stim_path, allow_pickle=True)
        roi_dict = np.load(roi_path, allow_pickle=True).item()
        ref_ch1 = np.load(ref_ch1_path, allow_pickle=True)

        releasing = iuf.identify_releasing_rois(pixel_RI_stim, roi_dict, alpha=ALPHA, min_ri=MIN_RI)
        if len(releasing) == 0:
            print_status('skipped', 'no releasing ROI')
            continue

        thres_mask = iuf.generate_adaptive_membrane_mask(ref_ch1, visualize=0)
        releasing_mask_raw = iuf.build_roi_mask(releasing)
        releasing_mask = releasing_mask_raw & thres_mask

        ROI_mask = iuf.build_roi_mask(roi_dict)
        anti_ROI_mask = ~binary_dilation(ROI_mask, iterations=DISTANCE_FROM_ROI)
        anti_ROI_mask = anti_ROI_mask & thres_mask

        n_raw = int(np.sum(releasing_mask_raw))
        n_thres = int(np.sum(releasing_mask))
        if n_raw == 0 or n_thres == 0:
            print_status('skipped', 'no ROI pixels')
            continue

        pixel_reduction = 1 - (n_thres / n_raw)
        if pixel_reduction > PIXEL_REDUCTION_THRESH:
            print_status('skipped', f'pixel_reduction={pixel_reduction:.3f}')
            continue

        ROI_RI_bins = np.nanmean(pixel_RI_bins[releasing_mask, :], axis=0)
        ROI_RI2_bins = np.nanmean(pixel_RI2_bins[releasing_mask, :], axis=0)
        neuropil_RI_bins = np.nanmean(pixel_RI_bins[anti_ROI_mask, :], axis=0)

        roi_tau, roi_r2 = fit_tau(ROI_RI_bins)
        roi2_tau, roi2_r2 = fit_tau(ROI_RI2_bins)
        neu_tau, neu_r2 = fit_tau(neuropil_RI_bins)

        row = {
            'transmitter': dataset['transmitter'],
            'sensor': dataset['sensor'],
            'recname': recname,
            'animal': recname.split('-')[0],
            'n_roi': len(releasing),
            'pixel_reduction': float(pixel_reduction),
            'roi_tau_s': roi_tau,
            'roi_r2': roi_r2,
            'roi_ctrl_tau_s': roi2_tau,
            'roi_ctrl_r2': roi2_r2,
            'neuropil_tau_s': neu_tau,
            'neuropil_r2': neu_r2,
            'roi_peak': float(np.nanmax(ROI_RI_bins)),
            'roi_ctrl_peak': float(np.nanmax(ROI_RI2_bins)),
            'neuropil_peak': float(np.nanmax(neuropil_RI_bins)),
            'roi_trace': ROI_RI_bins,
            'roi_ctrl_trace': ROI_RI2_bins,
            'neuropil_trace': neuropil_RI_bins,
            }

        all_rows.append(row)
        print_status('done', f'ROI tau={row["roi_tau_s"]:.3f}, neu tau={row["neuropil_tau_s"]:.3f}')


#%% save tables
csv_path = data_stem / 'roi_neuropil_timecourse_metrics.csv'
serialised = []
for row in all_rows:
    serialised.append({
        k: v for k, v in row.items()
        if not isinstance(v, np.ndarray)
        })

with open(csv_path, 'w', newline='', encoding='utf-8') as f:
    writer = csv.DictWriter(f, fieldnames=list(serialised[0].keys()))
    writer.writeheader()
    writer.writerows(serialised)


#%% print summary
print_statistics_section()
for dataset in DATASETS:
    curr = [row for row in all_rows if row['transmitter'] == dataset['transmitter']]
    print(f'{dataset["sensor"]}: {len(curr)} sessions')
    for key in ['roi_tau_s', 'roi_ctrl_tau_s', 'neuropil_tau_s', 'roi_peak', 'neuropil_peak']:
        vals = np.asarray([float(row[key]) for row in curr], dtype=float)
        vals = vals[np.isfinite(vals)]
        print(f'  {key}: n={vals.size}, median={np.nanmedian(vals):.4g}, mean={np.nanmean(vals):.4g}')


#%% plot
saved = [('csv', csv_path)]

fig, axs = plt.subplots(2, 2, figsize=(7.6, 4.8))
for ax, dataset in zip(axs[0], DATASETS):
    plot_trace_panel(ax, all_rows, dataset)
axs[0, 0].set_ylabel('RI')
axs[0, 1].set_ylabel('RI')

plot_tau_violin(axs[1, 0], all_rows, 'roi_tau_s')
axs[1, 0].set(title='ROI tau', ylabel='tau (s)')

plot_tau_violin(axs[1, 1], all_rows, 'neuropil_tau_s')
axs[1, 1].set(title='Neuropil tau', ylabel='tau (s)')

fig.tight_layout(w_pad=1.0, h_pad=1.0)

for ext in ['.png', '.pdf']:
    out = save_stem / f'dopamine_norepinephrine_roi_neuropil_timecourses{ext}'
    fig.savefig(out, dpi=300, bbox_inches='tight')
    saved.append((ext[1:], out))
plt.close(fig)

print_files_saved(saved)
