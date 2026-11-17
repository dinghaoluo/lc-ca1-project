# -*- coding: utf-8 -*-
'''
Created on Thu Nov 12 10:28:09 2026

late-trial split-half amplitude distributions and early/late comparison

@author: Dinghao Luo
'''

#%% imports
import sys
from pathlib import Path

repo_root = Path(__file__).resolve().parents[3]
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))

import project_paths as pp
import rec_list

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from common_functions import mpl_formatting
mpl_formatting()
import numpy as np
import pandas as pd
from scipy.stats import ks_2samp
from threadpoolctl import threadpool_limits

from analyse_early_late_ecdf_and_prediction import condition
from grid_response_test import align_trials, keep_aligned_trials, trial_amplitudes

repo_root   = Path(__file__).resolve().parents[3]
OUTPUT_ROOT = pp.DLIGHT_HOTSPOTS_OUTPUT_ROOT
FIGURE_ROOT = pp.DLIGHT_HOTSPOTS_FIGURE_ROOT
SEED        = 42


#%% split and amplitudes
def main():
    saved = pd.read_csv(OUTPUT_ROOT / 'followup_early_late_amplitudes.csv',
                        float_precision='round_trip')
    saved = saved.loc[saved.variant == 'unmatched']
    membership = pd.read_parquet(OUTPUT_ROOT / 'aligned_early_late_per_recording.parquet')
    recordings = membership.loc[membership.unmatched_included, 'rec_id'].tolist()
    finite = saved.loc[np.isfinite(saved.amplitude), 'amplitude'].to_numpy()
    low, high = np.percentile(finite, [1, 99])
    low, high = np.floor(low / 0.05) * 0.05, np.ceil(high / 0.05) * 0.05
    panels, rows, skipped = [], [], []
    for recname in recordings:
        behaviour_path = condition.RAW_ROOT / 'behaviour_profile' / f'{recname}.pkl'
        behaviour = condition.load_behaviour(behaviour_path)
        _, _, late, _ = condition.valid_trials(behaviour)
        if len(late) < 10:
            skipped.append((recname, len(late)))
            print(f'Skipped {recname}: {len(late)} late trials', flush=True)
            continue
        behaviour, profile, traces = condition.load_recording(recname)
        run = np.asarray(behaviour['run_onset_frames'], dtype=int)
        permuted = np.random.RandomState(SEED).permutation(late)
        n = len(late) // 2
        halves = (permuted[:n], permuted[n:2 * n])
        amplitudes, missing = [], []
        for trials in halves:
            cleaned, keep = keep_aligned_trials(align_trials(traces, run[trials]))
            trial_amp, _, _ = trial_amplitudes(cleaned)
            values = np.nanmean(trial_amp, axis=1)
            missing.append(int(np.sum(~np.isfinite(values))))
            amplitudes.append(values[np.isfinite(values)])
        original = saved.loc[saved.rec_id == recname]
        original_values = []
        for group in ('early', 'late'):
            values = original.loc[original.group == group, 'amplitude'].to_numpy()
            original_values.append(values[np.isfinite(values)])
        split_d = float(ks_2samp(*amplitudes).statistic)
        early_late_d = float(ks_2samp(*original_values).statistic)
        dropped = str(int(permuted[-1]) + 1) if len(late) % 2 else 'none'
        rows.append(dict(rec_id=recname, animal=recname.split('-')[0], n_late=len(late),
                         n_half=n, dropped=dropped, missing_first=missing[0],
                         missing_second=missing[1], split_d=split_d, early_late_d=early_late_d,
                         n_first=len(amplitudes[0]), n_second=len(amplitudes[1]),
                         n_early=len(original_values[0]), n_late_acr=len(original_values[1])))
        panels.append((recname, amplitudes))
        print(f'{recname}: {n} trials per half, split KS D={split_d:.6g}, early/late D={early_late_d:.6g}', flush=True)
    table = pd.DataFrame(rows)

    # preserve the original unmatched figure dimensions and recording-list order
    ncol = 6
    nrow = int(np.ceil(len(recordings) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(ncol * 2.65, nrow * 2.15),
                             sharex=True, sharey=True, squeeze=False)
    for i, (recname, amplitudes) in enumerate(panels):
        ax = axes.flat[i]
        for values, label, colour in zip(amplitudes, ('half 1', 'half 2'), ('0.5', '#19a819')):
            x = np.sort(values)
            ax.step(np.r_[min(low, x[0]), x, max(high, x[-1])],
                    np.r_[0, np.arange(1, len(x) + 1) / len(x), 1],
                    where='post', color=colour, linewidth=1.2, label=label)
        ax.set_title(recname, fontsize=9)
        ax.set_xlim(low, high)
        ax.set_ylim(0, 1)
        ax.spines[['top', 'right']].set_visible(False)
        ax.tick_params(labelsize=8)
        if i + ncol >= len(panels):
            ax.set_xlabel(chr(916) + 'dF/F')
            ax.tick_params(labelbottom=True)
    for i in range(len(panels), nrow * ncol):
        axes.flat[i].set_visible(False)
    fig.supylabel('Cumulative fraction of ACRs', fontsize=11)
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', ncol=2, frameon=False)
    fig.tight_layout(rect=(0.025, 0, 1, 0.985))
    path = FIGURE_ROOT / 'split_half_ecdf_late.png'
    if path.exists():
        raise FileExistsError(f'Existing figure left unchanged: {path}')
    fig.savefig(path, dpi=200, facecolor='white')
    plt.close(fig)

    median_split, median_original = table.split_d.median(), table.early_late_d.median()
    split_low, split_high = np.percentile(table.split_d, [25, 75])
    original_low, original_high = np.percentile(table.early_late_d, [25, 75])
    print(f'Complete: n={len(table)}; split D median={median_split:.9g}, IQR [{split_low:.9g}, {split_high:.9g}]; early/late median={median_original:.9g}, IQR [{original_low:.9g}, {original_high:.9g}]', flush=True)


#%% run
if __name__ == '__main__':
    with threadpool_limits(limits=1):
        main()
