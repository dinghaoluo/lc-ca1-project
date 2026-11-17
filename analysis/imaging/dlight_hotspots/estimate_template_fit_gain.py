# -*- coding: utf-8 -*-
'''
Created on Thu Nov 12 23:31:03 2026

cross-validated template-fit SNR gain over the window contrast

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
from threadpoolctl import threadpool_limits

from analyse_acr_consistency_aligned import (load_recording, event_frames,
                                              recnames)
from grid_response_test import (align_trials, keep_aligned_trials, trial_amplitudes,
                                quantify_grid_response, classify_up)

repo_root   = Path(__file__).resolve().parents[3]
OUTPUT_ROOT = pp.DLIGHT_HOTSPOTS_OUTPUT_ROOT
N_FOLDS     = 5
N_SHUFFLE   = 1000
SEED        = 42


#%% held-out SNR
def snr(values):
    n = np.isfinite(values).sum(axis=1)
    mean = np.divide(np.nansum(values, axis=1), n,
                     out=np.full(len(n), np.nan), where=n > 0)
    squared = np.nansum((values - mean[:, None]) ** 2, axis=1)
    sd = np.sqrt(np.divide(squared, n - 1,
                           out=np.full(len(n), np.nan), where=n >= 2))
    ratio = np.divide(mean, sd, out=np.full(len(n), np.nan),
                      where=np.isfinite(sd) & (sd > 0))
    return ratio, n, sd


#%% blocked cross-validation
def main():
    output = OUTPUT_ROOT / 'template_fit_gain.csv'
    if output.exists():
        raise FileExistsError(f'Existing output left unchanged: {output}')
    time = (np.arange(180) - 60) / 30
    rows = []
    for recname in list(recnames):
        behaviour, profile, traces = load_recording(recname)
        frames = event_frames(behaviour)
        aligned, keep = keep_aligned_trials(align_trials(traces, frames))
        contrast, _, _ = trial_amplitudes(aligned)
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
    table = pd.DataFrame(rows)
    table.to_csv(output, index=False)


#%% run
if __name__ == '__main__':
    with threadpool_limits(limits=1):
        main()
