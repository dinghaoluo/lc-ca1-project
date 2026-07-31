# -*- coding: utf-8 -*-
'''
Created on Fri Apr 10 2026

block shuffles and statistical summaries for the LC GLM analyses

@author: Dinghao Luo
'''

#%% imports
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib.pyplot as plt
import seaborn as sns


#%% block shuffles and collinearity summaries
def block_permute(arr, rng, block):
    '''permute within blocks; when block=None, permute the full array'''
    if block is None or block <= 1:
        return rng.permutation(arr)
    n = len(arr)
    n_blocks = int(np.ceil(n / block))
    idx = np.arange(n)
    blocks = [idx[i * block:(i + 1) * block] for i in range(n_blocks)]
    order = rng.permutation(n_blocks)
    new_idx = np.concatenate([blocks[i] for i in order])
    return arr[new_idx[:n]]

def save_figure(fig, output_stem):
    '''save a figure to png and pdf with the same stem.'''
    output_stem.parent.mkdir(parents=True, exist_ok=True)
    for ext in ['.png', '.pdf']:
        fig.savefig(f'{output_stem}{ext}', dpi=300, bbox_inches='tight')
    plt.close(fig)

def save_collinearity_heatmap(corr_df, cluname, output_stem):
    '''save one predictor-correlation heatmap.'''
    fig, ax = plt.subplots(figsize=(4.2, 3.6))
    sns.heatmap(
        corr_df, annot=True, fmt='.2f', cmap='coolwarm', vmin=-1, vmax=1,
        square=False, cbar_kws=None, ax=ax,
    )

    ax.set_title(f'predictor correlation: {cluname}')
    plt.tight_layout()
    save_figure(fig, output_stem)

def summarise_collinearity(pair_r_abs_cells,
                           pair_r_abs_by_sess,
                           vif_by_predictor,
                           glm_stem,
                           r_abs_thr,
                           output_prefix='model_comparison'):
    '''write shared collinearity summary tables and figures.'''
    pair_rows = []
    for (p1, p2), vals in pair_r_abs_cells.items():
        vals = np.array(vals, dtype=float)
        if len(vals) == 0:
            continue
        sem_val = stats.sem(vals, nan_policy='omit') if len(vals) > 1 else np.nan
        pair_rows.append({
            'pair': f'{p1} ↔ {p2}',
            'p1': p1,
            'p2': p2,
            'n_cells': len(vals),
            'median_|r|': np.median(vals),
            'mean_|r|': np.mean(vals),
            'sem_|r|': sem_val,
            f'prop_cells_|r|>={r_abs_thr}': np.mean(vals >= r_abs_thr)
        })

    if len(pair_rows) == 0:
        return

    pair_df_cells = pd.DataFrame(pair_rows).sort_values(
        [f'prop_cells_|r|>={r_abs_thr}', 'median_|r|'],
        ascending=[False, False]
    )

    sess_rows = []
    for (p1, p2), sess_dict in pair_r_abs_by_sess.items():
        per_sess_median = []
        for _, vals in sess_dict.items():
            if len(vals):
                per_sess_median.append(np.median(vals))
        if len(per_sess_median) == 0:
            continue
        per_sess_median = np.array(per_sess_median, dtype=float)
        sess_rows.append({
            'pair': f'{p1} ↔ {p2}',
            'p1': p1,
            'p2': p2,
            'n_sessions': len(per_sess_median),
            'median_of_session_medians_|r|': np.median(per_sess_median),
            f'prop_sessions_median_|r|>={r_abs_thr}': np.mean(per_sess_median >= r_abs_thr)
        })

    pair_df_sessions = pd.DataFrame(sess_rows).sort_values(
        [f'prop_sessions_median_|r|>={r_abs_thr}', 'median_of_session_medians_|r|'],
        ascending=[False, False]
    )

    summary_pairs = pd.merge(
        pair_df_cells,
        pair_df_sessions[['p1', 'p2', 'pair', 'n_sessions',
                          'median_of_session_medians_|r|',
                          f'prop_sessions_median_|r|>={r_abs_thr}']],
        on=['p1', 'p2', 'pair'],
        how='outer'
    ).sort_values(
        [f'prop_sessions_median_|r|>={r_abs_thr}', f'prop_cells_|r|>={r_abs_thr}', 'median_|r|'],
        ascending=[False, False, False]
    )

    summary_root = glm_stem / output_prefix / 'collinearity_summary'
    summary_root.mkdir(parents=True, exist_ok=True)
    pair_df_cells.to_csv(summary_root / 'pairwise_cells.csv', index=False)
    pair_df_sessions.to_csv(summary_root / 'pairwise_sessions.csv', index=False)
    summary_pairs.to_csv(summary_root / 'pairwise_summary_ranked.csv', index=False)

    all_preds = sorted({p for pair in pair_r_abs_cells for p in pair})
    mat = pd.DataFrame(np.eye(len(all_preds)), index=all_preds, columns=all_preds, dtype=float)
    for (p1, p2), vals in pair_r_abs_cells.items():
        if len(vals):
            med = float(np.median(vals))
            mat.loc[p1, p2] = med
            mat.loc[p2, p1] = med

    fig, ax = plt.subplots(figsize=(7, 6))
    sns.heatmap(
        mat, annot=True, fmt='.2f', cmap='coolwarm', vmin=0, vmax=1,
        square=True, cbar_kws={'label': 'median |r| across cells'}, ax=ax,
    )

    ax.set_title('pairwise median |r| across all cells')
    plt.tight_layout()
    save_figure(fig, summary_root / 'median_abs_r_heatmap')

    top_n = 10
    rank_sess = summary_pairs.dropna(
        subset=[f'prop_sessions_median_|r|>={r_abs_thr}']
    ).head(top_n)
    fig, ax = plt.subplots(figsize=(6.5, 3.8))
    ax.barh(rank_sess['pair'], rank_sess[f'prop_sessions_median_|r|>={r_abs_thr}'])
    ax.invert_yaxis()
    ax.set_xlabel(f'proportion of sessions with median |r| ≥ {r_abs_thr}')
    ax.set_title('top collinear predictor pairs (session-level)')
    ax.spines[['top', 'right']].set_visible(False)
    plt.tight_layout()
    save_figure(fig, summary_root / 'top_pairs_sessions')

    rank_mean = pair_df_cells.nlargest(top_n, 'mean_|r|')
    fig, ax = plt.subplots(figsize=(6.5, 3.8))
    ax.barh(rank_mean['pair'], rank_mean['mean_|r|'],
            xerr=rank_mean['sem_|r|'], capsize=3, color='#bdbdbd', edgecolor='k')
    ax.invert_yaxis()
    ax.set_xlabel('mean ± SEM of |r| across cells')
    ax.set_title('top collinear predictor pairs (cell-level)')
    ax.spines[['top', 'right']].set_visible(False)
    plt.tight_layout()
    save_figure(fig, summary_root / 'top_pairs_mean_sem')

    vif_rows = []
    for name, vals in vif_by_predictor.items():
        if len(vals):
            vals = np.array(vals, dtype=float)
            vif_rows.append({
                'predictor': name,
                'n_cells': len(vals),
                'median_VIF': np.median(vals),
                'mean_VIF': np.mean(vals),
                'prop_cells_VIF>=5': np.mean(vals >= 5),
                'prop_cells_VIF>=10': np.mean(vals >= 10)
            })
    if len(vif_rows):
        vif_df_summary = pd.DataFrame(vif_rows).sort_values(
            ['median_VIF', 'mean_VIF'],
            ascending=False
        )
        vif_df_summary.to_csv(summary_root / 'vif_summary.csv', index=False)
