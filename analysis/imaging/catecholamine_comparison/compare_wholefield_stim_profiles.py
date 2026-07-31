# -*- coding: utf-8 -*-
'''
Created on Thu Jun 25 2026

compare dopamine and norepinephrine whole-field stimulation profiles

@author: Dinghao Luo
'''

#%% imports
from pathlib import Path
import csv
import sys

import numpy as np
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
from console_formatting import print_files_saved, print_session, print_statistics_section
import project_paths as pp
import rec_list
mpl_formatting()


#%% parameters
SAMP_FREQ = 30
BEF       = 2
AFT       = 10
XAXIS     = np.arange((BEF + AFT) * SAMP_FREQ) / SAMP_FREQ - BEF

BASELINE_IDX = (XAXIS >= -1.0) & (XAXIS <= -0.15)
EARLY_IDX    = (XAXIS >= 1.0) & (XAXIS <= 2.0)
LATE_IDX     = (XAXIS >= 4.0) & (XAXIS <= 8.0)
POST_IDX     = (XAXIS >= 0.0) & (XAXIS <= 10.0)

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


#%% comparison plot
def plot_unpaired_metric(ax, rows, key, ylabel):
    groups = []
    colours = []
    labels = []
    for dataset in DATASETS:
        vals = [float(row[key]) for row in rows
                if row['transmitter'] == dataset['transmitter'] and np.isfinite(float(row[key]))]
        groups.append(np.asarray(vals, dtype=float))
        colours.append(dataset['colour'])
        labels.append(dataset['transmitter'])

    parts = ax.violinplot(groups, positions=[1, 2], showmedians=True, showextrema=False)
    for body, colour in zip(parts['bodies'], colours):
        body.set_facecolor(colour)
        body.set_edgecolor('none')
        body.set_alpha(0.65)
    parts['cmedians'].set_color('k')

    rng = np.random.default_rng(17)
    for xpos, vals, colour in zip([1, 2], groups, colours):
        jitter = rng.normal(0, 0.035, size=vals.size)
        ax.scatter(np.full(vals.size, xpos) + jitter, vals, s=12, color=colour,
                   edgecolor='none', alpha=0.65, zorder=3)

    stat, p_val = ranksums(groups[0], groups[1])
    ax.text(0.05, 0.96, f'p={p_val:.2e}', transform=ax.transAxes,
            ha='left', va='top', fontsize=7)
    print(f'{key}: ranksums z={stat:.3f}, p={p_val:.3e}')

    ax.set(xticks=[1, 2], xticklabels=labels, ylabel=ylabel, xlim=(0.5, 2.5))
    ax.spines[['top', 'right']].set_visible(False)


#%% main
all_rows = []
trace_mats = {}

for dataset in DATASETS:
    traces = []
    all_sess_stem = dataset['stem'] / 'all_sessions'

    for path in dataset['paths']:
        recname = Path(path).name
        print_session(f'{dataset["sensor"]}: {recname}')
        proc_path = all_sess_stem / recname / 'processed_data'

        trace_path = proc_path / f'{recname}_wholefield_dFF_stim.npy'
        trace = np.load(trace_path, allow_pickle=True).astype(float)
        baseline = np.nanmean(trace[BASELINE_IDX])
        trace = (trace - baseline) * 100

        post = np.where(POST_IDX)[0]
        peak_idx = post[np.nanargmax(trace[post])]
        peak_val = float(trace[peak_idx])
        early = float(np.nanmean(trace[EARLY_IDX]))
        late = float(np.nanmean(trace[LATE_IDX]))

        auc = np.trapezoid(
            np.clip(trace[post], 0, None),
            XAXIS[post],
            )

        ratio = late / early

        half_val = peak_val * 0.5
        after = np.arange(peak_idx, trace.shape[0])
        below = after[trace[after] <= half_val]
        half_recovery = float(XAXIS[below[0]] - XAXIS[peak_idx])

        row = {
            'transmitter': dataset['transmitter'],
            'sensor': dataset['sensor'],
            'recname': recname,
            'animal': recname.split('-')[0],
            'peak_dff_percent': peak_val,
            'time_to_peak_s': float(XAXIS[peak_idx]),
            'early_mean_percent': early,
            'late_mean_percent': late,
            'late_early_ratio': float(ratio),
            'positive_auc_percent_s': float(auc),
            'half_recovery_s': half_recovery,
            }
        all_rows.append(row)
        traces.append(trace)

    trace_mats[dataset['transmitter']] = np.asarray(traces, dtype=float)

csv_path = data_stem / 'wholefield_stim_profile_metrics.csv'
with open(csv_path, 'w', newline='', encoding='utf-8') as f:
    writer = csv.DictWriter(f, fieldnames=list(all_rows[0].keys()))
    writer.writeheader()
    writer.writerows(all_rows)


#%% print summary
print_statistics_section()
for dataset in DATASETS:
    curr = [row for row in all_rows if row['transmitter'] == dataset['transmitter']]
    print(f'{dataset["sensor"]}: {len(curr)} sessions')
    for key in ['peak_dff_percent', 'time_to_peak_s', 'late_early_ratio', 'positive_auc_percent_s']:
        vals = np.asarray([float(row[key]) for row in curr], dtype=float)
        vals = vals[np.isfinite(vals)]
        print(f'  {key}: median={np.nanmedian(vals):.4g}, mean={np.nanmean(vals):.4g}')


#%% plot traces
fig, ax = plt.subplots(figsize=(3.6, 2.4))
for dataset in DATASETS:
    mat = trace_mats[dataset['transmitter']]
    mean = np.nanmean(mat, axis=0)
    err = sem(mat, axis=0, nan_policy='omit')
    ax.plot(XAXIS, mean, color=dataset['colour'], linewidth=1.5,
            label=f'{dataset["transmitter"]} ({dataset["sensor"]})')
    ax.fill_between(XAXIS, mean - err, mean + err,
                    color=dataset['colour'], alpha=0.2, edgecolor='none')

ax.axvline(0, color='k', linestyle='--', linewidth=0.8)
ax.set(xlabel='time from stim. (s)', ylabel='baseline-centred dF/F (%)', xticks=(0, 5, 10))
ax.legend(frameon=False, fontsize=7)
ax.spines[['top', 'right']].set_visible(False)
fig.tight_layout()

saved = []
for ext in ['.png', '.pdf']:
    out = save_stem / f'dopamine_norepinephrine_wholefield_traces{ext}'
    fig.savefig(out, dpi=300, bbox_inches='tight')
    saved.append((ext[1:], out))
plt.close(fig)


#%% plot scalar metrics
metrics = [
    ('peak_dff_percent', 'peak dF/F (%)'),
    ('time_to_peak_s', 'time to peak (s)'),
    ('late_early_ratio', 'late/early response'),
    ('positive_auc_percent_s', 'positive AUC (% s)'),
    ]

fig, axs = plt.subplots(1, 4, figsize=(8.0, 2.2))
for ax, (key, ylabel) in zip(axs, metrics):
    plot_unpaired_metric(ax, all_rows, key, ylabel)
fig.tight_layout(w_pad=0.8)

for ext in ['.png', '.pdf']:
    out = save_stem / f'dopamine_norepinephrine_wholefield_metrics{ext}'
    fig.savefig(out, dpi=300, bbox_inches='tight')
    saved.append((ext[1:], out))
plt.close(fig)

print_files_saved([('csv', csv_path), *saved])
