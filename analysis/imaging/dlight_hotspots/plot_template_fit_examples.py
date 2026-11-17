# -*- coding: utf-8 -*-
'''
Created on Fri Nov 13 00:06:30 2026

held-out template fits for the existing single-trial examples

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
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import numpy as np
from threadpoolctl import threadpool_limits

from estimate_template_fit_gain import (load_recording, event_frames,
                                        align_trials, keep_aligned_trials,
                                        trial_amplitudes, quantify_grid_response,
                                        classify_up, N_FOLDS, N_SHUFFLE, SEED)
from plot_single_trial_examples import choose_recordings

repo_root   = Path(__file__).resolve().parents[3]
FIGURE_ROOT = pp.DLIGHT_HOTSPOTS_FIGURE_ROOT


#%% training-only templates
def fit_examples(recname):
    behaviour, profile, traces = load_recording(recname)
    frames = event_frames(behaviour)
    run = np.asarray(behaviour['run_onset_frames'], dtype=float)
    reward = np.asarray(behaviour['reward_times'], dtype=float)
    non_stop = np.asarray(behaviour['non_stop_trials'])
    non_fullstop = np.asarray(behaviour['non_fullstop_trials'])
    n = min(len(run), len(reward), len(non_stop), len(non_fullstop))
    valid = (np.isfinite(run[:n]) & (run[:n] != -1)
             & np.isfinite(reward[:n]) & (non_stop[:n] == 0)
             & (non_fullstop[:n] == 0))
    trial_numbers = np.flatnonzero(valid) + 1
    assert np.array_equal(run[:n][valid].astype(int), frames)
    first = (len(frames) - 5) // 2
    examples = np.arange(first, first + 5)
    aligned, keep = keep_aligned_trials(align_trials(traces, frames))
    contrast, _, _ = trial_amplitudes(aligned)
    average_amplitude = np.nanmean(contrast, axis=1)
    saved_up = np.flatnonzero(profile['Up'].to_numpy(dtype=bool))
    order = np.argsort(-average_amplitude[saved_up], kind='stable')
    selected = saved_up[order[:3]]
    if len(selected) != 3:
        raise ValueError(f'{recname}: fewer than three saved Up ACRs')
    coords = np.asarray(profile['roi_id'].tolist(), dtype=int)
    available = np.flatnonzero(np.any(keep, axis=0))
    blocks = np.array_split(available, N_FOLDS)
    time = (np.arange(180) - 60) / 30
    coefficients = np.full((len(profile), len(frames), 3), np.nan)
    fitted_curves = np.full((3, 5, 180), np.nan)
    templates, notes = {}, []
    for fold, test_trials in enumerate(blocks, 1):
        columns = np.flatnonzero(np.isin(examples, test_trials))
        if not len(columns):
            continue
        train_trials = np.setdiff1d(available, test_trials, assume_unique=True)
        assert not np.intersect1d(train_trials, test_trials).size
        training = quantify_grid_response(traces, frames[train_trials],
                                           n_shuffle=N_SHUFFLE, seed=SEED)
        up, _ = classify_up(training['amplitude'], training['effect_size'],
                             training['shuffle_amplitude'], profile['edge'].to_numpy())
        if not np.any(up):
            notes.append(f'{recname}: fold {fold} skipped, no training Up ACR.')
            del training
            continue
        profiles = training['mean_profile'][up].copy()
        profiles -= profiles[:, 30:60].mean(axis=1)[:, None]
        norm = np.linalg.norm(profiles[:, 60:180], axis=1)
        usable = np.isfinite(profiles).all(axis=1) & np.isfinite(norm) & (norm > 0)
        if not np.any(usable):
            notes.append(f'{recname}: fold {fold} skipped, no finite non-zero training Up template.')
            del training
            continue
        template = np.median(profiles[usable] / norm[usable, None], axis=0)
        design = np.column_stack((template, np.ones(180), time))
        if np.linalg.matrix_rank(design) < 3:
            notes.append(f'{recname}: fold {fold} skipped, template collinear with intercept and time.')
            del training
            continue
        held_keep = keep[:, test_trials]
        acr, trial = np.nonzero(held_keep)
        values = aligned[:, test_trials, :][held_keep]
        fitted = np.linalg.lstsq(design, values.T, rcond=None)[0]
        coefficients[acr, test_trials[trial], :] = fitted.T
        for column in columns:
            fitted_curves[:, column, :] = coefficients[selected, examples[column], :] @ design.T
        templates[fold] = template
        trials = ', '.join(str(x) for x in trial_numbers[examples[columns]])
        notes.append(f'{recname}: fold {fold}, {len(train_trials)} training trials, '
                     f'{len(test_trials)} held-out trials, {int(up.sum())} training Up ACRs, '
                     f'{int(usable.sum())} usable template ACRs; example behaviour trials {trials}.')
        del training
    missing_trials = examples[~np.isin(examples, available)]
    if len(missing_trials):
        notes.append(f'{recname}: example trials with no kept ACR, '
                     f'{trial_numbers[missing_trials].tolist()}; no fold or fit assigned.')

    # A row keeps its full signal range across all five trials and their fits.
    fig = plt.figure(figsize=(19, 8.3), layout='constrained')
    grid = fig.add_gridspec(3, 6, width_ratios=[1.2, 1, 1, 1, 1, 1])
    template_ax = fig.add_subplot(grid[:, 0])
    for fold, template in templates.items():
        template_ax.plot(time, template, linewidth=1.5, label=f'Fold {fold}')
    template_ax.axvline(0, color='black', linewidth=0.8, linestyle=':')
    template_ax.set(xlim=(-2, 4), xlabel='Time from run onset (s)',
                    ylabel='Fold template (normalised dF/F)', xticks=[-2, 0, 2, 4])
    if templates:
        template_ax.legend(frameon=False, loc='upper left')
    for row, acr in enumerate(selected):
        row_values = np.concatenate((aligned[acr, examples, :].ravel(),
                                     fitted_curves[row].ravel()))
        finite = row_values[np.isfinite(row_values)]
        limits = None
        if len(finite):
            low, high = finite.min(), finite.max()
            margin = 0.06 * (high - low) if high > low else 0.1
            limits = (low - margin, high + margin)
        for column, trial in enumerate(examples):
            ax = fig.add_subplot(grid[row, column + 1])
            ax.axvspan(-1, 0, facecolor='#6baed6', alpha=0.2, linewidth=0)
            ax.axvspan(0, 1, facecolor='#fdae6b', alpha=0.2, linewidth=0)
            ax.axvline(0, color='black', linewidth=0.8, linestyle=':')
            ax.plot(time, aligned[acr, trial], color='0.55', linewidth=1.0)
            ax.plot(time, fitted_curves[row, column], color='#228b22', linewidth=1.5)
            amplitude = coefficients[acr, trial, 0]
            change = contrast[acr, trial]
            ax.set_title(f'Trial {trial_numbers[trial]}\n'
                         f'A = {amplitude:.3g}, ΔdF/F = {change:.3g}', fontsize=10)
            ax.set(xlim=(-2, 4), xticks=[-2, 0, 2, 4])
            if limits is not None:
                ax.set_ylim(limits)
            if column == 0:
                y, x = coords[acr]
                ax.set_ylabel(f'ACR ({y}, {x})\ndF/F')
            if row == 2:
                ax.set_xlabel('Time from run onset (s)')
            ax.spines[['top', 'right']].set_visible(False)
    legend = [Line2D([], [], color='0.55', label='Single-trial trace'),
              Line2D([], [], color='#228b22', label='A × template + C + D × t'),
              Patch(facecolor='#6baed6', alpha=0.2, label='Baseline [-1, 0) s'),
              Patch(facecolor='#fdae6b', alpha=0.2, label='Response [0, 1) s'),
              Line2D([], [], color='black', linestyle=':', label='Run onset')]
    fig.legend(handles=legend, loc='outside lower center', ncol=5, frameon=False)
    path = FIGURE_ROOT / f'template_fit_examples_{recname}.png'
    if path.exists():
        raise FileExistsError(f'Existing figure left unchanged: {path}')
    fig.savefig(path, dpi=200, facecolor='white')
    plt.close(fig)
    trials = ', '.join(str(x) for x in trial_numbers[examples])
    positions = ', '.join(str(x + 1) for x in examples)
    notes.insert(0, f'{recname}: behaviour trial numbers (one-based) {trials}, '
                 f'valid-list positions {positions}; row ACR coordinates (zero-based grid row, column), '
                 + ', '.join(f'({coords[acr, 0]}, {coords[acr, 1]})' for acr in selected)
                 + '; saved profile indices '
                 + ', '.join(str(profile.index[acr]) for acr in selected)
                 + '; trial-averaged amplitudes (dF/F) '
                 + ', '.join(f'{average_amplitude[acr]:.9g}' for acr in selected) + '.')
    missing = [(tuple(coords[acr]), int(trial_numbers[trial]))
               for acr in selected for trial in examples if not keep[acr, trial]]
    notes.append(f'{recname}: missing kept windows among plotted ACR-trials: {missing or "none"}. '
                 f'Figure: figures/dlight_hotspots/{path.name}.')
    return notes


#%% run
def main():
    _, _, chosen = choose_recordings()
    for _, row in chosen:
        path = FIGURE_ROOT / f'template_fit_examples_{row.rec_id}.png'
        if path.exists():
            raise FileExistsError(f'Existing figure left unchanged: {path}')
    for label, row in chosen:
        notes = fit_examples(row.rec_id)
        for note in notes:
            print(f'{label}: {note}', flush=True)


if __name__ == '__main__':
    with threadpool_limits(limits=1):
        main()
