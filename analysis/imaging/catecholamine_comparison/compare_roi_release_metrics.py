# -*- coding: utf-8 -*-
'''
Created on Thu Jun 25 2026

compare dopamine and norepinephrine axon-ROI release metrics

@author: Dinghao Luo
'''

#%% imports
from pathlib import Path
import csv
import sys

import numpy as np
from scipy.stats import ranksums, sem, spearmanr, ttest_1samp
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
import project_paths as pp
import rec_list
mpl_formatting()


#%% parameters
ALPHA  = 0.05
MIN_RI = 0.1

MIN_ROI    = 5
MIN_TRIALS = 5
N_SHUF     = 500
RNG_SEED   = 42

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


#%% release summaries and plots
def top_fraction_contribution(x, frac=0.1):
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    x = np.clip(x, 0, None)
    if x.size == 0 or np.nansum(x) == 0:
        return np.nan
    n_top = max(1, int(np.ceil(x.size * frac)))
    return float(np.sum(np.sort(x)[-n_top:]) / np.sum(x))

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

    rng = np.random.default_rng(91)
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

for dataset in DATASETS:
    all_sess_stem = dataset['stem'] / 'all_sessions'

    for path in dataset['paths']:
        recname = Path(path).name
        print_session(f'{dataset["sensor"]}: {recname}')
        proc_path = all_sess_stem / recname / 'processed_data'
        pixel_path = proc_path / f'{recname}_pixel_RI_stim.npy'
        roi_path = proc_path / f'{recname}_ROI_dict.npy'

        pixel_RI_stim = np.load(pixel_path, allow_pickle=True)
        roi_dict = np.load(roi_path, allow_pickle=True).item()

        n_trials = pixel_RI_stim.shape[2]
        roi_ids = list(roi_dict.keys())
        roi_trial_RI = np.full((len(roi_ids), n_trials), np.nan)

        for roi_row, roi_id in enumerate(roi_ids):
            roi = roi_dict[roi_id]
            vals = pixel_RI_stim[roi['ypix'], roi['xpix'], :]
            roi_trial_RI[roi_row, :] = np.nanmean(vals, axis=0)

        good_roi = np.sum(np.isfinite(roi_trial_RI), axis=1) >= MIN_TRIALS
        roi_trial_RI = roi_trial_RI[good_roi, :]
        if roi_trial_RI.shape[0] < MIN_ROI:
            print_status('skipped', 'too few ROIs')
            continue

        n_trials = int(np.nanmax(np.sum(np.isfinite(roi_trial_RI), axis=1)))
        if n_trials < MIN_TRIALS:
            print_status('skipped', 'too few trials')
            continue

        roi_mean_RI = np.nanmean(roi_trial_RI, axis=1)
        valid = np.isfinite(roi_trial_RI)
        n_valid = np.sum(valid, axis=1)
        positive = np.sum((roi_trial_RI > 0) & valid, axis=1)
        roi_pos_prob = positive / n_valid
        releasing = np.zeros(roi_trial_RI.shape[0], dtype=bool)

        for roi_row in range(roi_trial_RI.shape[0]):
            vals = roi_trial_RI[roi_row, :]
            vals = vals[np.isfinite(vals)]
            if vals.size > 2:
                _, p_val = ttest_1samp(vals, 0, alternative='greater')
                if p_val < ALPHA and np.mean(vals) > MIN_RI:
                    releasing[roi_row] = True

        rng = np.random.default_rng(RNG_SEED)
        top10_shuf = []
        for _ in range(N_SHUF):
            shuf = np.array(roi_trial_RI, copy=True)
            for trial in range(shuf.shape[1]):
                shuf[:, trial] = rng.permutation(shuf[:, trial])
            top10_shuf.append(top_fraction_contribution(np.nanmean(shuf, axis=1)))
        top10_shuf = np.asarray(top10_shuf)
        split_shuf = []
        for _ in range(N_SHUF):
            odd = np.nanmean(roi_trial_RI[:, ::2], axis=1)
            even = np.nanmean(roi_trial_RI[:, 1::2], axis=1)
            even = rng.permutation(even)
            valid = np.isfinite(odd) & np.isfinite(even)
            split_shuf.append(float(spearmanr(odd[valid], even[valid])[0]))
        split_shuf = np.asarray(split_shuf)

        top10 = top_fraction_contribution(roi_mean_RI)
        gini_vals = roi_mean_RI[np.isfinite(roi_mean_RI)]
        gini_vals = gini_vals[gini_vals >= 0]
        if gini_vals.size == 0 or np.nansum(gini_vals) == 0:
            gini_positive_RI = np.nan
        else:
            gini_vals = np.sort(gini_vals)
            n = gini_vals.size
            gini_positive_RI = (2 * np.sum(np.arange(1, n + 1) * gini_vals) / (n * np.sum(gini_vals))) - ((n + 1) / n)

        odd = np.nanmean(roi_trial_RI[:, ::2], axis=1)
        even = np.nanmean(roi_trial_RI[:, 1::2], axis=1)
        valid = np.isfinite(odd) & np.isfinite(even)
        split_r = float(spearmanr(odd[valid], even[valid])[0])

        row = {
            'transmitter': dataset['transmitter'],
            'sensor': dataset['sensor'],
            'recname': recname,
            'animal': recname.split('-')[0],
            'n_roi': int(roi_trial_RI.shape[0]),
            'n_trials': n_trials,
            'prop_releasing': float(np.nanmean(releasing)),
            'mean_RI': float(np.nanmean(roi_mean_RI)),
            'median_RI': float(np.nanmedian(roi_mean_RI)),
            'mean_positive_prob': float(np.nanmean(roi_pos_prob)),
            'gini_positive_RI': float(gini_positive_RI),
            'top10_contribution': float(top10),
            'top10_shuf_median': float(np.nanmedian(top10_shuf)),
            'top10_minus_shuffle': float(top10 - np.nanmedian(top10_shuf)),
            'split_half_r': float(split_r),
            'split_half_shuf_median': float(np.nanmedian(split_shuf)),
            'split_half_minus_shuffle': float(split_r - np.nanmedian(split_shuf)),
            }

        all_rows.append(row)
        print_status('done', f"prop rel.={row['prop_releasing']:.3f}, mean RI={row['mean_RI']:.3f}")

csv_path = data_stem / 'roi_release_metrics.csv'
with open(csv_path, 'w', newline='', encoding='utf-8') as f:
    writer = csv.DictWriter(f, fieldnames=list(all_rows[0].keys()))
    writer.writeheader()
    writer.writerows(all_rows)


#%% print summary
print_statistics_section()
for dataset in DATASETS:
    curr = [row for row in all_rows if row['transmitter'] == dataset['transmitter']]
    print(f'{dataset["sensor"]}: {len(curr)} sessions')
    for key in ['prop_releasing', 'mean_RI', 'mean_positive_prob', 'top10_contribution', 'split_half_r']:
        vals = np.asarray([float(row[key]) for row in curr], dtype=float)
        vals = vals[np.isfinite(vals)]
        print(f'  {key}: median={np.nanmedian(vals):.4g}, mean={np.nanmean(vals):.4g}, sem={sem(vals):.4g}')


#%% plotting
metrics = [
    ('prop_releasing', 'proportion releasing'),
    ('mean_RI', 'mean ROI RI'),
    ('mean_positive_prob', 'positive trial probability'),
    ('gini_positive_RI', 'Gini of positive RI'),
    ('top10_contribution', 'top 10% contribution'),
    ('split_half_r', 'split-half Spearman r'),
    ]

fig, axs = plt.subplots(2, 3, figsize=(7.2, 4.5))
for ax, (key, ylabel) in zip(axs.ravel(), metrics):
    plot_unpaired_metric(ax, all_rows, key, ylabel)
fig.tight_layout(w_pad=0.8, h_pad=1.0)

saved = []
for ext in ['.png', '.pdf']:
    out = save_stem / f'dopamine_norepinephrine_roi_release_metrics{ext}'
    fig.savefig(out, dpi=300, bbox_inches='tight')
    saved.append((ext[1:], out))
plt.close(fig)

print_files_saved([('csv', csv_path), *saved])
