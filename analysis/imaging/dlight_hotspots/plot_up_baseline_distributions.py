'''
Created on Mon Nov 24 2026

plot Up and baseline single-trial dLight distributions

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
from threadpoolctl import threadpool_limits

from analyse_acr_consistency_aligned import load_recording, event_frames
from grid_response_test import align_trials, keep_aligned_trials

OUTPUT_ROOT = pp.DLIGHT_HOTSPOTS_OUTPUT_ROOT
FIGURE_ROOT = pp.DLIGHT_HOTSPOTS_FIGURE_ROOT
GROUPS = (
    ('Up response', 'up_response', '#228b22'),
    ('Up baseline', 'up_baseline', '#4477aa'),
    ('Other response', 'other_response', '0.5'),
)


#%% recording values
def recording_values(recname):
    behaviour, profile, traces = load_recording(recname)
    frames = event_frames(behaviour)
    aligned, keep = keep_aligned_trials(align_trials(traces, frames))
    response = (
        np.nanmean(aligned[:, :, 60:90], axis=2)
        - np.nanmean(aligned[:, :, 30:60], axis=2)
    )
    baseline = (
        np.nanmean(aligned[:, :, 30:60], axis=2)
        - np.nanmean(aligned[:, :, 0:30], axis=2)
    )
    up = profile['Up'].to_numpy(dtype=bool)
    values = {
        'up_response': response[up[:, None] & keep],
        'up_baseline': baseline[up[:, None] & keep],
        'other_response': response[~up[:, None] & keep],
    }
    return tuple(np.asarray(values[key], dtype=float).ravel()
                 for _, key, _ in GROUPS)


#%% figure
def ecdf(values, low, high):
    values = np.sort(values[np.isfinite(values)])
    return (
        np.r_[min(low, values[0]), values, max(high, values[-1])],
        np.r_[0, np.arange(1, len(values) + 1) / len(values), 1],
    )


def main():
    path = FIGURE_ROOT / 'up_baseline_distributions.png'
    if path.exists():
        raise FileExistsError(path)

    panels = []
    pooled = [[] for _ in GROUPS]
    for recname in rec_list.recDLightDbh:
        groups = recording_values(recname)
        panels.append((recname, groups))
        for index, values in enumerate(groups):
            pooled[index].extend(values[np.isfinite(values)].tolist())
        print(f'{recname}: Up response, Up baseline and Other response', flush=True)

    low, high = np.percentile(np.concatenate(pooled), [1, 99])
    ncol = 6
    nrow = int(np.ceil(len(panels) / ncol))
    fig, axes = plt.subplots(
        nrow, ncol, figsize=(ncol * 2.65, nrow * 2.15),
        sharex=True, sharey=True, squeeze=False,
    )
    for index, (recname, groups) in enumerate(panels):
        ax = axes.flat[index]
        for values, (label, _, colour) in zip(groups, GROUPS):
            finite = values[np.isfinite(values)]
            if not len(finite):
                continue
            x, y = ecdf(finite, low, high)
            ax.step(x, y, where='post', color=colour, linewidth=1.2, label=label)
        ax.set_title(recname, fontsize=9)
        ax.spines[['top', 'right']].set_visible(False)
        ax.tick_params(labelsize=8)
        if index + ncol >= len(panels):
            ax.set_xlabel('dF/F change')
            ax.tick_params(labelbottom=True)
    for ax in axes.flat[len(panels):]:
        ax.set_visible(False)
    axes[0, 0].set(xlim=(low, high), ylim=(0, 1), yticks=[0, 0.5, 1])
    fig.supylabel('Cumulative fraction of ACR-trials', fontsize=11)
    handles = [
        plt.Line2D([], [], color=colour, lw=1.2, label=label)
        for label, _, colour in GROUPS
    ]
    fig.legend(handles=handles, loc='upper center', ncol=3, frameon=False)
    fig.tight_layout(rect=(0.025, 0, 1, 0.977))
    FIGURE_ROOT.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=200, facecolor='white')
    plt.close(fig)

    n_animals = len({recname.split('-')[0] for recname, _ in panels})
    lines = [
        '',
        '## 24 November 2026: Up and baseline single-trial distributions',
        '',
        f'I used all valid trials from {len(panels)} recordings in {n_animals} animals '
        f'with the per-ACR full-window NaN rule. The pooled medians are '
        f'Up response={np.median(pooled[0]):.9g} (n={len(pooled[0])}), '
        f'Up baseline={np.median(pooled[1]):.9g} (n={len(pooled[1])}), and '
        f'Other response={np.median(pooled[2]):.9g} (n={len(pooled[2])}).',
        '',
        'Figure: figures/dlight_hotspots/up_baseline_distributions.png.',
        '',
    ]
    with (OUTPUT_ROOT / 'PROGRESS.md').open('a', encoding='utf-8') as stream:
        stream.write('\n'.join(lines))


if __name__ == '__main__':
    with threadpool_limits(limits=1):
        main()
