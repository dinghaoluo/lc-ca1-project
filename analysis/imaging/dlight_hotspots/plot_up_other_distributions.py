'''
Created on Wed Nov 18 22:17:17 2026

plot single-trial response distributions for saved Up and other ACRs

@author: Dinghao Luo
'''

#%% imports
import numpy as np
from scipy.ndimage import gaussian_filter1d
from threadpoolctl import threadpool_limits

from check_smoothed_template_and_co_recruitment import (
    load_recording, event_frames, recnames, align_trials, keep_aligned_trials,
    trial_amplitudes, quantify_grid_response, classify_up, OUTPUT_ROOT, N_FOLDS,
    N_SHUFFLE, SEED,
)
import project_paths as pp

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from common_functions import mpl_formatting
mpl_formatting()

FIGURE_ROOT = pp.DLIGHT_HOTSPOTS_FIGURE_ROOT
GROUPS      = ('Up', 'Other')
COLOURS     = ('#228b22', '0.5')


#%% held-out template coefficients
def template_amplitudes(recname, frames, profile, traces, aligned, keep):
    available = np.flatnonzero(np.any(keep, axis=0))
    time = (np.arange(180) - 60) / 30
    template_coef = np.full(keep.shape, np.nan)
    used_folds = 0
    # re-call training Up ACRs in each block, as in the smoothed template gain analysis
    for fold, test_trials in enumerate(np.array_split(available, N_FOLDS), 1):
        train_trials = np.setdiff1d(available, test_trials, assume_unique=True)
        if not len(train_trials) or not len(test_trials):
            print(f'{recname}, fold {fold}: empty training or test block', flush=True)
            continue
        training = quantify_grid_response(traces, frames[train_trials],
                                          n_shuffle=N_SHUFFLE, seed=SEED)
        up, _ = classify_up(training['amplitude'], training['effect_size'],
                            training['shuffle_amplitude'], profile['edge'].to_numpy())
        profiles = training['mean_profile'][up].copy()
        del training
        profiles -= profiles[:, 30:60].mean(axis=1)[:, None]
        norm = np.linalg.norm(profiles[:, 60:180], axis=1)
        usable = np.isfinite(profiles).all(axis=1) & np.isfinite(norm) & (norm > 0)
        if not np.any(usable):
            print(f'{recname}, fold {fold}: no usable training Up profile', flush=True)
            continue
        template = np.median(profiles[usable] / norm[usable, None], axis=0)
        template = gaussian_filter1d(template, sigma=3, mode='nearest')
        design = np.column_stack((template, np.ones(180), time))
        if np.linalg.matrix_rank(design) < 3:
            print(f'{recname}, fold {fold}: rank-deficient template fit', flush=True)
            continue
        held_keep = keep[:, test_trials]
        acr, trial = np.nonzero(held_keep)
        values = aligned[:, test_trials, :][held_keep]
        fitted = np.linalg.lstsq(design, values.T, rcond=None)[0][0]
        template_coef[acr, test_trials[trial]] = fitted
        used_folds += 1
    return template_coef, used_folds


#%% recording panels
def plot_distributions(panels, x_label, path):
    low, high = np.percentile(np.concatenate([x for _, groups in panels for x in groups]),
                              [1, 99])
    ncol = 6
    nrow = (len(panels) + ncol - 1) // ncol
    fig, axes = plt.subplots(nrow, ncol, figsize=(ncol * 2.65, nrow * 2.15),
                             sharex=True, sharey=True, squeeze=False)
    for i, (recname, groups) in enumerate(panels):
        ax = axes.flat[i]
        for values, label, colour in zip(groups, GROUPS, COLOURS):
            if not len(values):
                continue
            x = np.sort(values)
            ax.step(np.r_[min(low, x[0]), x, max(high, x[-1])],
                    np.r_[0, np.arange(1, len(x) + 1) / len(x), 1],
                    where='post', color=colour, linewidth=1.2, label=label)
        ax.set_title(recname, fontsize=9)
        ax.spines[['top', 'right']].set_visible(False)
        ax.tick_params(labelsize=8)
        if i + ncol >= len(panels):
            ax.set_xlabel(x_label)
            ax.tick_params(labelbottom=True)
    for ax in axes.flat[len(panels):]:
        ax.set_visible(False)
    axes[0, 0].set(xlim=(low, high), ylim=(0, 1), yticks=[0, 0.5, 1])
    fig.supylabel('Cumulative fraction of ACR-trials', fontsize=11)
    handles = [plt.Line2D([], [], color=colour, lw=1.2, label=label)
               for label, colour in zip(GROUPS, COLOURS)]
    fig.legend(handles=handles, loc='upper center', ncol=2, frameon=False)
    fig.tight_layout(rect=(0.025, 0, 1, 0.977))
    fig.savefig(path, dpi=200, facecolor='white')
    plt.close(fig)


def medians(groups):
    up, other = [np.median(values) if len(values) else np.nan for values in groups]
    return up, other, up - other


#%% run and report
def main():
    paths = (FIGURE_ROOT / 'up_other_dff_change_distributions.png',
             FIGURE_ROOT / 'up_other_template_coef_distributions.png')
    for path in paths:
        if path.exists():
            raise FileExistsError(path)

    panels = ([], [])
    rows = []
    for recname in recnames:
        behaviour, profile, traces = load_recording(recname)
        frames = event_frames(behaviour)
        aligned, keep = keep_aligned_trials(align_trials(traces, frames))
        contrast, _, _ = trial_amplitudes(aligned)
        coefficient, used_folds = template_amplitudes(recname, frames, profile, traces,
                                                      aligned, keep)
        saved_up = profile['Up'].to_numpy(dtype=bool)
        summaries = []
        for values, measure_panels in zip((contrast, coefficient), panels):
            groups = tuple(values[group & keep & np.isfinite(values)]
                           for group in (saved_up[:, None], ~saved_up[:, None]))
            measure_panels.append((recname, groups))
            summaries.append((*map(len, groups), *medians(groups)))
        rows.append((recname, used_folds, *summaries))
        print(f'{recname}: {used_folds}/{N_FOLDS} folds, '
              f'Up minus Other dF/F change={summaries[0][-1]:.6g}, '
              f'A={summaries[1][-1]:.6g}', flush=True)

    FIGURE_ROOT.mkdir(parents=True, exist_ok=True)
    for measure_panels, label, path in zip(panels, ('dF/F change', 'A'), paths):
        plot_distributions(measure_panels, label, path)

    n_animals = len({recname.split('-')[0] for recname, *_ in rows})
    lines = [
        '', '## 18 November 2026: single-trial distributions in saved Up and Other ACRs', '',
        f'I used all selected trials from {len(rows)} recordings in {n_animals} animals, '
        f'with the saved Up flags and the full-window NaN rule. The dF/F contrast is '
        f'mean frames 60:90 minus mean frames 30:60. The coefficient A comes from the '
        f'five-block held-out fit with the training Up template smoothed at sigma 3 frames, '
        f'an intercept and a linear time term. The curves use all finite values for '
        f'their measure and the axes show the pooled 1st to 99th percentile.', '',
    ]
    incomplete = [(name, folds) for name, folds, *_ in rows if folds < N_FOLDS]
    for name, folds in incomplete:
        lines.append(f'{name} contributed {folds}/{N_FOLDS} folds for A. Its dF/F contrast uses all kept trials.')
    if incomplete:
        lines.append('')
    for index, (measure_panels, label) in enumerate(zip(panels, ('dF/F change', 'A'))):
        groups = tuple(np.concatenate([groups[group] for _, groups in measure_panels])
                       for group in range(2))
        up, other, difference = medians(groups)
        by_recording = np.asarray([row[index + 2][2:] for row in rows])
        low, middle, high = np.nanpercentile(by_recording, [25, 50, 75], axis=0)
        lines.append(f'{label}, pooled over ACR-trials: Up median {up:.6g} (n={len(groups[0])}), '
                     f'Other median {other:.6g} (n={len(groups[1])}), Up minus Other {difference:.6g}. '
                     f'Giving each recording equal weight, the medians are {middle[0]:.6g} for Up '
                     f'and {middle[1]:.6g} for Other, with a median paired recording difference '
                     f'of {middle[2]:.6g} (IQR [{low[2]:.6g}, {high[2]:.6g}]).')

    lines.extend([
        '', '| Recording | Folds for A | dF/F Up n | dF/F Other n | dF/F Up median | dF/F Other median | dF/F Up minus Other | A Up n | A Other n | A Up median | A Other median | A Up minus Other |',
        '| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |',
    ])
    for name, folds, contrast, coefficient in rows:
        cells = ' | '.join(f'{value:.6g}' for value in (*contrast, *coefficient))
        lines.append(f'| {name} | {folds} | {cells} |')
    lines.extend(['', *(f'`figures/dlight_hotspots/{path.name}`' for path in paths), ''])
    with (OUTPUT_ROOT / 'PROGRESS.md').open('a', encoding='utf-8') as stream:
        stream.write('\n'.join(lines))


if __name__ == '__main__':
    with threadpool_limits(limits=1):
        main()
