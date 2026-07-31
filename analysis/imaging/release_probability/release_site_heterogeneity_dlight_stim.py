# -*- coding: utf-8 -*-
'''
Created on Tue May 12 2026

quantify response concentration and odd-even stability across LC axonal ROIs,
using shuffled ROI identities and a nepicastat comparison

@author: Dinghao Luo
'''

#%% imports
from pathlib import Path
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.stats import sem, spearmanr, wilcoxon
from tqdm import tqdm

repo_root = Path(__file__).resolve().parents[3]
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))

import project_paths as pp
from common_functions import mpl_formatting
from console_formatting import (
    print_files_saved,
    print_session,
    print_statistics_section,
    print_status,
)
import rec_list

mpl_formatting()


#%% parameters
ALPHA  = 0.05
MIN_RI = 0.1

N_SHUF = 500
RNG_SEED = 42

MIN_ROI    = 5
MIN_TRIALS = 5

PRINT_SESSION_SUMMARY = True
SAVE_OUTPUT = False
SAVE_FIGURES = True


#%% paths
paths = rec_list.pathdLightLCOpto + \
        rec_list.pathdLightLCOptoDbhBlock

nepicastat_recs = {Path(path).name for path in rec_list.pathdLightLCOptoDbhBlock}

dLight_stem   = pp.HPC_DLIGHT_LC_OPTO_STEM
all_sess_stem = dLight_stem / 'all_sessions'
fig_stem      = pp.HPC_DLIGHT_LC_OPTO_RELEASE_SITE_HETEROGENEITY_FIGURES_STEM


#%% heterogeneity summaries and plots
def top_fraction_contribution(x, frac=0.1):
    '''Fraction of total positive response carried by the strongest ROIs.'''
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    x = np.clip(x, 0, None)
    if x.size == 0 or np.nansum(x) == 0:
        return np.nan
    n_top = max(1, int(np.ceil(x.size * frac)))
    return np.sum(np.sort(x)[-n_top:]) / np.sum(x)

def plot_real_vs_shuffle(results, key, shuffle_key, ylabel, title, savename):
    curr = [res for res in results if res['condition'] == 'dLight']
    real = np.asarray([res[key] for res in curr], dtype=float)
    shuf = np.asarray([res[shuffle_key] for res in curr], dtype=float)
    valid = np.isfinite(real) & np.isfinite(shuf)
    real = real[valid]
    shuf = shuf[valid]

    diff = real - shuf
    _, p_val = wilcoxon(diff, alternative='greater')

    fig, ax = plt.subplots(figsize=(2.2, 2.6))

    for r, s in zip(real, shuf):
        ax.plot([1, 2], [s, r],
                color='grey', alpha=.25, linewidth=.6, zorder=1)

    parts = ax.violinplot(
        [shuf, real],
        positions=[1, 2],
        showmeans=False,
        showmedians=True,
        showextrema=False
        )
    for pc, c in zip(parts['bodies'], ['grey', 'darkgreen']):
        pc.set_facecolor(c)
        pc.set_edgecolor('none')
        pc.set_alpha(.35)
    parts['cmedians'].set_color('k')
    parts['cmedians'].set_linewidth(1)

    ax.scatter(np.ones_like(shuf),
               shuf,
               color='grey', ec='none',
               s=10, alpha=.5, zorder=3)
    ax.scatter(np.ones_like(real) * 2,
               real,
               color='darkgreen', ec='none',
               s=10, alpha=.6, zorder=3)

    ymax = np.nanmax(np.concatenate([real, shuf]))
    ymin = np.nanmin(np.concatenate([real, shuf]))
    yrange = ymax - ymin
    ax.text(
        1.5, ymax + yrange * .08,
        f'p={p_val:.1e}',
        ha='center', va='bottom',
        fontsize=7
        )

    ax.set(
        xticks=[1, 2],
        xticklabels=['shuffle', 'real'],
        ylabel=ylabel,
        title=title
        )
    for s in ['top', 'right']:
        ax.spines[s].set_visible(False)

    fig.tight_layout()
    fig_stem.mkdir(parents=True, exist_ok=True)
    for ext in ['.png', '.pdf']:
        fig.savefig(fig_stem / f'{savename}{ext}', dpi=300, bbox_inches='tight')
    plt.close(fig)

def plot_heterogeneity_summary(results):
    curr = [res for res in results if res['condition'] == 'dLight']

    top10 = np.asarray([res['top10_contribution'] for res in curr], dtype=float)
    gini = np.asarray([res['gini_positive_RI'] for res in curr], dtype=float)
    split = np.asarray([res['split_half_r'] for res in curr], dtype=float)
    prop = np.asarray([res['prop_releasing'] for res in curr], dtype=float)

    fig, axs = plt.subplots(2, 2, figsize=(5.2, 4.2))
    axs = axs.ravel()

    axs[0].hist(
        top10,
        bins=np.linspace(0, .5, 16),
        color='darkgreen',
        alpha=.65,
        edgecolor='none'
        )
    axs[0].axvline(np.nanmedian(top10),
                   color='firebrick',
                   linestyle='--',
                   linewidth=1)
    axs[0].set(
        xlabel='top 10% contribution',
        ylabel='session count',
        title='response concentration'
        )

    axs[1].hist(
        split,
        bins=np.linspace(0, 1, 16),
        color='darkgreen',
        alpha=.65,
        edgecolor='none'
        )
    axs[1].axvline(np.nanmedian(split),
                   color='firebrick',
                   linestyle='--',
                   linewidth=1)
    axs[1].set(
        xlabel='odd-even Spearman r',
        ylabel='session count',
        title='hotspot stability'
        )

    axs[2].scatter(
        prop,
        top10,
        color='darkgreen',
        ec='none',
        s=16,
        alpha=.65
        )
    r_prop_top, p_prop_top = spearmanr(prop, top10, nan_policy='omit')
    axs[2].text(
        .04, .96,
        f'r={r_prop_top:.2f}\np={p_prop_top:.1e}',
        transform=axs[2].transAxes,
        ha='left', va='top',
        fontsize=7
        )
    axs[2].set(
        xlabel='prop. releasing ROIs',
        ylabel='top 10% contribution',
        title='recruitment vs concentration'
        )

    axs[3].scatter(
        gini,
        split,
        color='darkgreen',
        ec='none',
        s=16,
        alpha=.65
        )
    r_gini_split, p_gini_split = spearmanr(gini, split, nan_policy='omit')
    axs[3].text(
        .04, .96,
        f'r={r_gini_split:.2f}\np={p_gini_split:.1e}',
        transform=axs[3].transAxes,
        ha='left', va='top',
        fontsize=7
        )
    axs[3].set(
        xlabel='Gini coefficient',
        ylabel='odd-even Spearman r',
        title='inequality vs stability'
        )

    for ax in axs:
        for s in ['top', 'right']:
            ax.spines[s].set_visible(False)

    fig.tight_layout()
    fig_stem.mkdir(parents=True, exist_ok=True)
    for ext in ['.png', '.pdf']:
        fig.savefig(fig_stem / f'dLight_release_site_heterogeneity_summary{ext}', dpi=300, bbox_inches='tight')
    plt.close(fig)

def plot_condition_comparison(results):
    fig, axs = plt.subplots(1, 3, figsize=(5.6, 2.3))

    keys = [
        ('top10_contribution', 'top 10% contribution'),
        ('split_half_r', 'odd-even Spearman r'),
        ('gini_positive_RI', 'Gini coefficient'),
        ]

    for ax, (key, ylabel) in zip(axs, keys):
        vals = []
        for condition in ['dLight', 'nepicastat']:
            curr = [res for res in results if res['condition'] == condition]
            vals.append(np.asarray([res[key] for res in curr], dtype=float))

        parts = ax.violinplot(
            vals,
            positions=[1, 2],
            showmeans=False,
            showmedians=True,
            showextrema=False
            )
        for pc, c in zip(parts['bodies'], ['darkgreen', 'firebrick']):
            pc.set_facecolor(c)
            pc.set_edgecolor('none')
            pc.set_alpha(.35)
        parts['cmedians'].set_color('k')
        parts['cmedians'].set_linewidth(1)

        ax.scatter(np.ones_like(vals[0]),
                   vals[0],
                   color='darkgreen', ec='none',
                   s=10, alpha=.55)
        ax.scatter(np.ones_like(vals[1]) * 2,
                   vals[1],
                   color='firebrick', ec='none',
                   s=10, alpha=.55)

        ax.set(
            xticks=[1, 2],
            xticklabels=['dLight', 'nepi'],
            ylabel=ylabel
            )
        for s in ['top', 'right']:
            ax.spines[s].set_visible(False)

    fig.tight_layout()
    fig_stem.mkdir(parents=True, exist_ok=True)
    for ext in ['.png', '.pdf']:
        fig.savefig(fig_stem / f'dLight_release_site_heterogeneity_nepicastat_comparison{ext}', dpi=300, bbox_inches='tight')
    plt.close(fig)


#%% main
rng = np.random.default_rng(RNG_SEED)

all_results = []

for path in tqdm(paths, desc='sessions'):
    recname = Path(path).name
    condition = 'nepicastat' if recname in nepicastat_recs else 'dLight'

    pixel_RI_stim_path = all_sess_stem / recname / 'processed_data' / f'{recname}_pixel_RI_stim.npy'
    roi_dict_path      = all_sess_stem / recname / 'processed_data' / f'{recname}_ROI_dict.npy'

    pixel_RI_stim = np.load(pixel_RI_stim_path, allow_pickle=True)
    roi_dict      = np.load(roi_dict_path, allow_pickle=True).item()

    n_trials = pixel_RI_stim.shape[2]
    roi_ids = list(roi_dict.keys())
    # this matrix stores mean RI per axon ROI per trial
    roi_trial_RI = np.full((len(roi_ids), n_trials), np.nan)

    for roi_row, roi_id in enumerate(roi_ids):
        roi = roi_dict[roi_id]
        vals = pixel_RI_stim[roi['ypix'], roi['xpix'], :]
        roi_trial_RI[roi_row, :] = np.nanmean(vals, axis=0)

    good_roi = np.sum(np.isfinite(roi_trial_RI), axis=1) >= MIN_TRIALS
    roi_trial_RI = roi_trial_RI[good_roi, :]
    if roi_trial_RI.shape[0] < MIN_ROI or roi_trial_RI.shape[1] < MIN_TRIALS:
        if PRINT_SESSION_SUMMARY:
            print_session(recname)
            print_status('skipped', 'too few ROIs/trials')
        continue

    # use the same release criterion as the release-probability scripts
    releasing = np.zeros(roi_trial_RI.shape[0], dtype=bool)

    for roi_row in range(roi_trial_RI.shape[0]):
        means = roi_trial_RI[roi_row, :]
        means = means[np.isfinite(means)]
        if means.size > 2:
            _, p = wilcoxon(means, alternative='greater')
            if p < ALPHA and np.mean(means) > MIN_RI:
                releasing[roi_row] = True

    roi_mean_RI = np.nanmean(roi_trial_RI, axis=1)
    roi_positive_prob = np.nanmean(roi_trial_RI > 0, axis=1)

    # stable site identity: do high-response ROIs remain high in odd/even trials?
    odd_mean  = np.nanmean(roi_trial_RI[:, ::2], axis=1)
    even_mean = np.nanmean(roi_trial_RI[:, 1::2], axis=1)
    valid_split = np.isfinite(odd_mean) & np.isfinite(even_mean)
    if np.sum(valid_split) < 5:
        split_r = np.nan
    else:
        split_r, _ = spearmanr(odd_mean[valid_split], even_mean[valid_split])

    split_r_shuf = []
    top10_shuf = []
    for _ in range(N_SHUF):
        even_mean_shuf = rng.permutation(even_mean)
        valid_split = np.isfinite(odd_mean) & np.isfinite(even_mean_shuf)
        if np.sum(valid_split) < 5:
            split_r_shuf.append(np.nan)
        else:
            split_r_shuf.append(spearmanr(
                odd_mean[valid_split], even_mean_shuf[valid_split]
                )[0])
        # shuffle ROI identity independently within each trial
        shuf = np.array(roi_trial_RI, copy=True)
        for trial in range(shuf.shape[1]):
            shuf[:, trial] = rng.permutation(shuf[:, trial])
        top10_shuf.append(top_fraction_contribution(np.nanmean(shuf, axis=1)))

    split_r_shuf = np.asarray(split_r_shuf)
    top10_shuf   = np.asarray(top10_shuf)

    top10 = top_fraction_contribution(roi_mean_RI)
    top10_p = np.nanmean(top10_shuf >= top10)
    split_r_p = np.nanmean(split_r_shuf >= split_r)

    # Gini uses non-negative mean RI values only
    gini_vals = roi_mean_RI[np.isfinite(roi_mean_RI)]
    gini_vals = gini_vals[gini_vals >= 0]
    if gini_vals.size == 0 or np.nansum(gini_vals) == 0:
        gini_positive_RI = np.nan
    else:
        gini_vals = np.sort(gini_vals)
        n = gini_vals.size
        gini_positive_RI = (2 * np.sum((np.arange(1, n + 1)) * gini_vals) / (n * np.sum(gini_vals))) - ((n + 1) / n)

    result = {
        'recname': recname,
        'condition': condition,
        'n_roi': roi_trial_RI.shape[0],
        'n_trials': roi_trial_RI.shape[1],
        'prop_releasing': np.nanmean(releasing),
        'mean_RI': np.nanmean(roi_mean_RI),
        'median_RI': np.nanmedian(roi_mean_RI),
        'mean_positive_prob': np.nanmean(roi_positive_prob),
        'gini_positive_RI': gini_positive_RI,
        'top10_contribution': top10,
        'top10_shuf_median': np.nanmedian(top10_shuf),
        'top10_p': top10_p,
        'split_half_r': split_r,
        'split_half_shuf_median': np.nanmedian(split_r_shuf),
        'split_half_p': split_r_p,
    }
    all_results.append(result)

    if PRINT_SESSION_SUMMARY:
        print_session(recname)
        print_status(
            'done',
            f"{condition}: ROIs={result['n_roi']}, trials={result['n_trials']}, "
            f"prop rel.={result['prop_releasing']:.3f}, "
            f"top10={result['top10_contribution']:.3f} "
            f"(shuf={result['top10_shuf_median']:.3f}, p={result['top10_p']:.3f}), "
            f"split r={result['split_half_r']:.3f} "
            f"(shuf={result['split_half_shuf_median']:.3f}, p={result['split_half_p']:.3f})"
        )


#%% summary
print_statistics_section()
print('release-site heterogeneity summary:')
print(f'valid sessions: {len(all_results)}')

for condition in ['dLight', 'nepicastat']:
    curr = [res for res in all_results if res['condition'] == condition]
    print(f'\n{condition}:')
    for key in [
            'n_roi',
            'n_trials',
            'prop_releasing',
            'mean_RI',
            'mean_positive_prob',
            'gini_positive_RI',
            'top10_contribution',
            'top10_shuf_median',
            'split_half_r',
            'split_half_shuf_median',
            ]:
        vals = np.asarray([res[key] for res in curr], dtype=float)
        print(
            f'{key}: '
            f'{np.nanmean(vals):.3f} +/- {sem(vals, nan_policy="omit"):.3f}; '
            f'median={np.nanmedian(vals):.3f}'
        )

    vals = np.asarray([res['top10_contribution'] - res['top10_shuf_median'] for res in curr])
    _, p_top = wilcoxon(vals, alternative='greater')
    vals = np.asarray([res['split_half_r'] - res['split_half_shuf_median'] for res in curr])
    _, p_split = wilcoxon(vals, alternative='greater')

    print(f'top10 > shuffle, Wilcoxon p={p_top:.3e}')
    print(f'split-half r > shuffle, Wilcoxon p={p_split:.3e}')


#%% plotting
if SAVE_FIGURES:
    plot_real_vs_shuffle(
        all_results,
        key='top10_contribution',
        shuffle_key='top10_shuf_median',
        ylabel='top 10% contribution',
        title='response concentration',
        savename='dLight_top10_contribution_real_vs_shuffle'
        )
    plot_real_vs_shuffle(
        all_results,
        key='split_half_r',
        shuffle_key='split_half_shuf_median',
        ylabel='odd-even Spearman r',
        title='hotspot stability',
        savename='dLight_split_half_stability_real_vs_shuffle'
        )
    plot_heterogeneity_summary(all_results)
    plot_condition_comparison(all_results)
    print_files_saved([
        ('figures', fig_stem),
    ])


#%% optional save
if SAVE_OUTPUT:
    out_path = pp.HPC_DLIGHT_LC_OPTO_STEM / 'release_site_heterogeneity_dlight_stim.csv'
    keys = list(all_results[0].keys())
    with open(out_path, 'w', encoding='utf-8') as f:
        f.write(','.join(keys) + '\n')
        for res in all_results:
            f.write(','.join(str(res[key]) for key in keys) + '\n')
    print_files_saved([
        ('summary csv', out_path),
    ])
