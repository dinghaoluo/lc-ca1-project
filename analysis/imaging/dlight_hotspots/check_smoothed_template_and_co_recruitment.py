# -*- coding: utf-8 -*-
'''
Created on Fri Nov 13 17:53:44 2026

check smoothed template gain and trialwise Up versus other co-recruitment

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
from scipy.ndimage import gaussian_filter1d
from scipy.stats import pearsonr
from threadpoolctl import threadpool_limits

from estimate_template_fit_gain import (load_recording, event_frames,
                                        recnames, align_trials,
                                        keep_aligned_trials, trial_amplitudes,
                                        quantify_grid_response, classify_up, snr,
                                        OUTPUT_ROOT, N_FOLDS, N_SHUFFLE, SEED)


#%% correlations across kept trials
def group_mean(values):
    count = np.isfinite(values).sum(axis=0)
    return np.divide(np.nansum(values, axis=0), count,
                     out=np.full(values.shape[1], np.nan), where=count > 0)


def correlation(left, right):
    finite = np.isfinite(left) & np.isfinite(right)
    n = int(finite.sum())
    if n < 3:
        return np.nan, n, 'fewer than three paired trials'
    if np.ptp(left[finite]) == 0 or np.ptp(right[finite]) == 0:
        return np.nan, n, 'constant trial means'
    return float(pearsonr(left[finite], right[finite])[0]), n, ''


def co_recruitment(recname, profile, aligned, contrast):
    up = profile['Up'].to_numpy(dtype=bool)
    row = dict(rec_id=recname, animal=recname.split('-')[0], n_valid_trials=aligned.shape[1],
               n_up=int(up.sum()), n_other=int((~up).sum()),
               r_resp=np.nan, r_base=np.nan, n_resp_trials=0, n_base_trials=0,
               z_resp=np.nan, z_base=np.nan, z_resp_minus_base=np.nan,
               resp_reason='', base_reason='', status='ok')
    if not np.any(up) or np.all(up):
        row['status'] = 'no saved Up ACR' if not np.any(up) else 'no other valid ACR'
        return row
    baseline = aligned[:, :, 30:60].mean(axis=2) - aligned[:, :, 0:30].mean(axis=2)
    for label, values in (('resp', contrast), ('base', baseline)):
        r, n, reason = correlation(group_mean(values[up]), group_mean(values[~up]))
        row[f'r_{label}'] = r
        row[f'n_{label}_trials'] = n
        row[f'{label}_reason'] = reason
        if np.isfinite(r) and abs(r) < 1:
            row[f'z_{label}'] = float(np.arctanh(r))
        elif np.isfinite(r):
            row[f'{label}_reason'] = 'perfect correlation; Fisher z not finite'
    row['z_resp_minus_base'] = row['z_resp'] - row['z_base']
    if not np.isfinite(row['z_resp_minus_base']):
        row['status'] = 'unavailable paired Fisher z'
    return row


def quantiles(values):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    return np.percentile(values, [25, 50, 75]) if len(values) else (np.nan, np.nan, np.nan)


#%% blocked cross-validation
def main():
    output = OUTPUT_ROOT / 'smoothed_template_gain.csv'
    if output.exists():
        raise FileExistsError(f'Existing output left unchanged: {output}')
    co_output = OUTPUT_ROOT / 'co_recruitment.csv'
    if co_output.exists():
        raise FileExistsError(f'Existing output left unchanged: {co_output}')
    co_rows = []
    time = (np.arange(180) - 60) / 30
    rows = []
    for recname in list(recnames):
        behaviour, profile, traces = load_recording(recname)
        frames = event_frames(behaviour)
        aligned, keep = keep_aligned_trials(align_trials(traces, frames))
        contrast, _, _ = trial_amplitudes(aligned)
        co_rows.append(co_recruitment(recname, profile, aligned, contrast))
        available = np.flatnonzero(np.any(keep, axis=0))
        blocks = np.array_split(available, N_FOLDS)
        template_coef = np.full(keep.shape, np.nan)
        held_contrast = np.full(keep.shape, np.nan)
        used_folds = 0
        for fold, test_trials in enumerate(blocks, 1):
            train_trials = np.setdiff1d(available, test_trials, assume_unique=True)
            if not len(train_trials) or not len(test_trials):
                continue
            training = quantify_grid_response(traces, frames[train_trials],
                                               n_shuffle=N_SHUFFLE, seed=SEED)
            up, _ = classify_up(training['amplitude'], training['effect_size'],
                                 training['shuffle_amplitude'], profile['edge'].to_numpy())
            if not np.any(up):
                print(f'{recname}, fold {fold}: no training Up ACR; skipped', flush=True)
                del training
                continue
            profiles = training['mean_profile'][up].copy()
            profiles -= profiles[:, 30:60].mean(axis=1)[:, None]
            norm = np.linalg.norm(profiles[:, 60:180], axis=1)
            usable = np.isfinite(profiles).all(axis=1) & np.isfinite(norm) & (norm > 0)
            if not np.any(usable):
                del training
                continue
            template = np.median(profiles[usable] / norm[usable, None], axis=0)
            template = gaussian_filter1d(template, sigma=3, mode='nearest')
            design = np.column_stack((template, np.ones(180), time))
            if np.linalg.matrix_rank(design) < 3:
                del training
                continue
            held_keep = keep[:, test_trials]
            acr, trial = np.nonzero(held_keep)
            values = aligned[:, test_trials, :][held_keep]
            fitted = np.linalg.lstsq(design, values.T, rcond=None)[0][0]
            template_coef[acr, test_trials[trial]] = fitted
            held_contrast[acr, test_trials[trial]] = contrast[acr, test_trials[trial]]
            used_folds += 1
            del training
        snr_a, n_held, sd_a = snr(template_coef)
        snr_d, n_d, sd_d = snr(held_contrast)
        assert np.array_equal(n_held, n_d)
        gain = np.divide(snr_a, snr_d, out=np.full(len(profile), np.nan),
                         where=np.isfinite(snr_a) & np.isfinite(snr_d) & (snr_d > 0))
        saved_up = profile['Up'].to_numpy(dtype=bool)
        for label, group in (('Up', saved_up), ('other', ~saved_up)):
            finite = group & np.isfinite(gain)
            values = gain[finite]
            low, median, high = np.percentile(values, [25, 50, 75]) if len(values) else (np.nan, np.nan, np.nan)
            rows.append(dict(rec_id=recname, animal=recname.split('-')[0], group=label,
                             n_valid_trials=len(frames), n_trials_any_kept=len(available),
                             n_all_missing_trials=int(np.sum(~np.any(keep, axis=0))),
                             n_used_folds=used_folds, n_skipped_folds=N_FOLDS - used_folds,
                             n_valid_acr=int(group.sum()), n_gain=int(finite.sum()),
                             n_excluded_snr_d_nonpositive=int(np.sum(group & np.isfinite(snr_d) & (snr_d <= 0))),
                             n_fewer_two_held=int(np.sum(group & (n_held < 2))),
                             n_zero_sd_a=int(np.sum(group & (n_held >= 2) & (sd_a == 0))),
                             n_zero_sd_d=int(np.sum(group & (n_held >= 2) & (sd_d == 0))),
                             n_excluded_total=int(np.sum(group & ~np.isfinite(gain))),
                             median_gain=median, gain_q25=low, gain_q75=high))
        print(f'{recname}: {used_folds}/5 folds; Up gain={rows[-2]["median_gain"]:.6g} (n={rows[-2]["n_gain"]}), other={rows[-1]["median_gain"]:.6g} (n={rows[-1]["n_gain"]})', flush=True)
        co = co_rows[-1]
        print(f'{recname}: r_resp={co["r_resp"]:.6g}, r_base={co["r_base"]:.6g}, n_trials={co["n_resp_trials"]}; {co["status"]}', flush=True)
    table = pd.DataFrame(rows)
    table.to_csv(output, index=False)
    co_table = pd.DataFrame(co_rows)
    co_table.to_csv(co_output, index=False)


#%% run
if __name__ == '__main__':
    with threadpool_limits(limits=1):
        main()
