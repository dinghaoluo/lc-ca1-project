# -*- coding: utf-8 -*-
'''
Created on Sun Nov 15 11:06:25 2026

compare held-out A-d correlations and templates from other recordings of each animal

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

import numpy as np
import pandas as pd
from scipy.stats import pearsonr
from threadpoolctl import threadpool_limits

from estimate_template_fit_gain import (load_recording, event_frames,
                                        recnames, align_trials,
                                        keep_aligned_trials, trial_amplitudes,
                                        quantify_grid_response, classify_up, snr,
                                        OUTPUT_ROOT, N_FOLDS, N_SHUFFLE, SEED)
from check_smoothed_template_and_co_recruitment import gaussian_filter1d, quantiles

TIME = (np.arange(180) - 60) / 30


#%% profiles and fits
def normalised_profiles(profiles):
    profiles = profiles.copy()
    profiles -= profiles[:, 30:60].mean(axis=1)[:, None]
    norm = np.linalg.norm(profiles[:, 60:180], axis=1)
    usable = np.isfinite(profiles).all(axis=1) & np.isfinite(norm) & (norm > 0)
    return profiles[usable] / norm[usable, None]


def design_matrix(profiles):
    if len(profiles) == 0:
        return None
    template = np.median(profiles, axis=0)
    template = gaussian_filter1d(template, sigma=3, mode='nearest')
    design = np.column_stack((template, np.ones(180), TIME))
    return design if np.linalg.matrix_rank(design) == 3 else None


def blocked_coefficients(recname, profile, traces, frames, aligned, keep, contrast):
    available = np.flatnonzero(np.any(keep, axis=0))
    blocks = np.array_split(available, N_FOLDS)
    template_coef = np.full(keep.shape, np.nan)
    held_contrast = np.full(keep.shape, np.nan)
    notes = []
    for fold, test_trials in enumerate(blocks, 1):
        train_trials = np.setdiff1d(available, test_trials, assume_unique=True)
        reason, n_up, n_template = '', 0, 0
        design = None
        if not len(train_trials) or not len(test_trials):
            reason = 'empty training or held-out block'
        else:
            training = quantify_grid_response(traces, frames[train_trials],
                                               n_shuffle=N_SHUFFLE, seed=SEED)
            up, _ = classify_up(training['amplitude'], training['effect_size'],
                                 training['shuffle_amplitude'], profile['edge'].to_numpy())
            n_up = int(up.sum())
            if not n_up:
                reason = 'no training Up ACR'
            else:
                profiles = normalised_profiles(training['mean_profile'][up])
                n_template = len(profiles)
                if not n_template:
                    reason = 'no finite non-zero training Up template'
                else:
                    design = design_matrix(profiles)
                    if design is None:
                        reason = 'template collinear with intercept and time'
            del training
        n_pairs = 0
        if design is not None:
            held_keep = keep[:, test_trials]
            acr, trial = np.nonzero(held_keep)
            values = aligned[:, test_trials, :][held_keep]
            fitted = np.linalg.lstsq(design, values.T, rcond=None)[0][0]
            template_coef[acr, test_trials[trial]] = fitted
            held_contrast[acr, test_trials[trial]] = contrast[acr, test_trials[trial]]
            n_pairs = len(fitted)
        notes.append(dict(rec_id=recname, fold=fold, n_train=len(train_trials),
                          n_test=len(test_trials), n_up=n_up, n_template=n_template,
                          n_pairs=n_pairs, reason=reason))
    return template_coef, held_contrast, notes


#%% correlations and descriptive summaries
def acr_correlations(template_coef, contrast):
    finite = np.isfinite(template_coef) & np.isfinite(contrast)
    n = finite.sum(axis=1)
    template_coef = np.where(finite, template_coef, 0.0)
    dff_change = np.where(finite, contrast, 0.0)
    mean_a = np.divide(template_coef.sum(axis=1), n, out=np.zeros(len(n)), where=n > 0)
    mean_d = np.divide(dff_change.sum(axis=1), n, out=np.zeros(len(n)), where=n > 0)
    template_coef = np.where(finite, template_coef - mean_a[:, None], 0.0)
    dff_change = np.where(finite, contrast - mean_d[:, None], 0.0)
    ss_a, ss_d = np.sum(template_coef * template_coef, axis=1), np.sum(dff_change * dff_change, axis=1)
    denominator = np.sqrt(ss_a * ss_d)
    r = np.divide(np.sum(template_coef * dff_change, axis=1), denominator,
                   out=np.full(len(n), np.nan), where=(n >= 2) & (denominator > 0))
    r = np.clip(r, -1, 1)
    return r, n, ss_a, ss_d, finite


def summarise(recname, profile, template_coef, contrast, variant, pooled, pair_values,
              template_counts=None):
    r, n, ss_a, ss_d, finite = acr_correlations(template_coef, contrast)
    pair_values[variant][0].append(template_coef[finite])
    pair_values[variant][1].append(contrast[finite])
    gain = np.full(len(profile), np.nan)
    snr_d = np.full(len(profile), np.nan)
    if variant == 'other_recordings_same_animal':
        snr_a, n_a, _ = snr(template_coef)
        snr_d, n_d, _ = snr(contrast)
        assert np.array_equal(n_a, n_d)
        gain = np.divide(snr_a, snr_d, out=gain,
                         where=np.isfinite(snr_a) & np.isfinite(snr_d) & (snr_d > 0))
    rows = []
    saved_up = profile['Up'].to_numpy(dtype=bool)
    for label, group in (('Up', saved_up), ('other', ~saved_up)):
        valid_r = group & np.isfinite(r)
        valid_gain = group & np.isfinite(gain)
        pooled[variant][label]['r'].extend(r[valid_r].tolist())
        pooled[variant][label]['gain'].extend(gain[valid_gain].tolist())
        low, median, high = quantiles(r[valid_r])
        glow, gmedian, ghigh = quantiles(gain[valid_gain])
        row = dict(scope='recording', rec_id=recname, animal=recname.split('-')[0],
                   variant=variant, group=label, n_recordings=1, n_animals=1,
                   n_valid_acr=int(group.sum()), n_r=int(valid_r.sum()),
                   median_r=median, r_q25=low, r_q75=high,
                   n_pairs=int(finite[group].sum()),
                   n_fewer_two_pairs=int(np.sum(group & (n < 2))),
                   n_constant_a=int(np.sum(group & (n >= 2) & (ss_a == 0))),
                   n_constant_d=int(np.sum(group & (n >= 2) & (ss_d == 0))),
                   n_excluded_r=int(np.sum(group & ~np.isfinite(r))),
                   n_gain=int(valid_gain.sum()) if variant != 'blocked_smoothed' else np.nan,
                   median_gain=gmedian, gain_q25=glow, gain_q75=ghigh,
                   n_excluded_snr_d_nonpositive=int(np.sum(group & np.isfinite(snr_d) & (snr_d <= 0))) if variant != 'blocked_smoothed' else np.nan,
                   n_excluded_gain=int(np.sum(group & ~np.isfinite(gain))) if variant != 'blocked_smoothed' else np.nan,
                   pooled_pair_r=np.nan)
        if template_counts:
            row.update(template_counts)
        rows.append(row)
    return rows


#%% run both template variants
def main():
    output = OUTPUT_ROOT / 'template_variants.csv'
    if output.exists():
        raise FileExistsError(f'Existing output left unchanged: {output}')
    recordings = list(recnames)
    variants = ('blocked_smoothed', 'other_recordings_same_animal')
    pooled = {variant: {label: {'r': [], 'gain': []} for label in ('Up', 'other')}
              for variant in variants}
    pair_values = {variant: [[], []] for variant in variants}
    rows, contributors = [], {}
    for recname in recordings:
        behaviour, profile, traces = load_recording(recname)
        frames = event_frames(behaviour)
        aligned, keep = keep_aligned_trials(align_trials(traces, frames))
        contrast, _, _ = trial_amplitudes(aligned)
        template_coef, held_contrast, notes = blocked_coefficients(
            recname, profile, traces, frames, aligned, keep, contrast)
        result = summarise(recname, profile, template_coef, held_contrast,
                           'blocked_smoothed', pooled, pair_values)
        rows.extend(result)
        saved_up = profile['Up'].to_numpy(dtype=bool)
        up_profiles = np.nanmean(aligned[saved_up], axis=1)
        normalised = normalised_profiles(up_profiles)
        contributors[recname] = dict(profiles=normalised, n_saved=int(saved_up.sum()))
        print(f'{recname}: blocked median r Up={result[0]["median_r"]:.6g}, '
              f'other={result[1]["median_r"]:.6g}; '
              f'{len(normalised)} usable saved Up profiles', flush=True)
    for recname in recordings:
        animal = recname.split('-')[0]
        sources = [other for other in recordings
                   if other != recname and other.split('-')[0] == animal]
        source_profiles = np.concatenate([contributors[other]['profiles'] for other in sources], axis=0)
        counts = dict(n_template_recordings=len(sources),
                      n_template_acr=len(source_profiles),
                      n_template_saved_up=sum(contributors[other]['n_saved'] for other in sources))
        counts['n_template_excluded'] = counts['n_template_saved_up'] - counts['n_template_acr']
        design = design_matrix(source_profiles)
        behaviour, profile, traces = load_recording(recname)
        frames = event_frames(behaviour)
        aligned, keep = keep_aligned_trials(align_trials(traces, frames))
        contrast, _, _ = trial_amplitudes(aligned)
        template_coef = np.full(keep.shape, np.nan)
        held_contrast = np.full(keep.shape, np.nan)
        if design is not None:
            template_coef[keep] = np.linalg.lstsq(design, aligned[keep].T, rcond=None)[0][0]
            held_contrast[keep] = contrast[keep]
        rows.extend(summarise(recname, profile, template_coef, held_contrast,
                              'other_recordings_same_animal', pooled, pair_values, counts))
        print(f'{recname}: animal template from {counts["n_template_acr"]} ACRs, '
              f'gain Up={rows[-2]["median_gain"]:.6g}, other={rows[-1]["median_gain"]:.6g}', flush=True)
    recording_table = pd.DataFrame(rows)
    for variant in variants:
        for label in ('Up', 'other'):
            subset = recording_table.loc[(recording_table.variant == variant)
                                          & (recording_table.group == label)]
            low, median, high = quantiles(pooled[variant][label]['r'])
            glow, gmedian, ghigh = quantiles(pooled[variant][label]['gain'])
            row = dict(scope='pooled_acrs', rec_id='all', animal='all',
                       variant=variant, group=label,
                       n_recordings=int(subset.loc[subset.n_r > 0, 'rec_id'].nunique()),
                       n_animals=int(subset.loc[subset.n_r > 0, 'animal'].nunique()),
                       median_r=median, r_q25=low, r_q75=high,
                       median_gain=gmedian, gain_q25=glow, gain_q75=ghigh,
                       pooled_pair_r=np.nan)
            for column in ('n_valid_acr', 'n_r', 'n_pairs', 'n_fewer_two_pairs',
                           'n_constant_a', 'n_constant_d', 'n_excluded_r'):
                row[column] = int(subset[column].sum())
            for column in ('n_gain', 'n_excluded_snr_d_nonpositive', 'n_excluded_gain'):
                row[column] = int(subset[column].sum()) if variant != 'blocked_smoothed' else np.nan
            rows.append(row)
        template_coef, dff_change = [np.concatenate(parts) for parts in pair_values[variant]]
        r = float(pearsonr(template_coef, dff_change)[0])
        rows.append(dict(scope='pooled_pairs', rec_id='all', animal='all',
                         variant=variant, group='all', n_recordings=len(recordings),
                         n_animals=len({recname.split('-')[0] for recname in recordings}),
                         n_pairs=len(template_coef), pooled_pair_r=r))
    pd.DataFrame(rows).to_csv(output, index=False)
    print(pd.DataFrame(rows).loc[lambda table: table.scope != 'recording',
          ['scope', 'variant', 'group', 'median_r', 'r_q25', 'r_q75',
           'median_gain', 'gain_q25', 'gain_q75', 'pooled_pair_r', 'n_pairs']].to_string(index=False), flush=True)


#%% run
if __name__ == '__main__':
    with threadpool_limits(limits=1):
        main()
