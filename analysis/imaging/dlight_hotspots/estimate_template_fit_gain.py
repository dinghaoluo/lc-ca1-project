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
    rows, folds, pooled = [], [], {'Up': [], 'other': []}
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
            fold_row = dict(rec_id=recname, fold=fold, n_train=len(train_trials),
                            n_test=len(test_trials), n_up=0, n_template=0,
                            n_held_pairs=0, reason='')
            if not len(train_trials) or not len(test_trials):
                fold_row['reason'] = 'empty training or held-out block'
                folds.append(fold_row)
                continue
            training = quantify_grid_response(traces, frames[train_trials],
                                               n_shuffle=N_SHUFFLE, seed=SEED)
            up, _ = classify_up(training['amplitude'], training['effect_size'],
                                 training['shuffle_amplitude'], profile['edge'].to_numpy())
            fold_row['n_up'] = int(up.sum())
            if not np.any(up):
                fold_row['reason'] = 'no training Up ACR'
                folds.append(fold_row)
                print(f'{recname}, fold {fold}: no training Up ACR; skipped', flush=True)
                del training
                continue
            profiles = training['mean_profile'][up].copy()
            profiles -= profiles[:, 30:60].mean(axis=1)[:, None]
            norm = np.linalg.norm(profiles[:, 60:180], axis=1)
            usable = np.isfinite(profiles).all(axis=1) & np.isfinite(norm) & (norm > 0)
            fold_row['n_template'] = int(usable.sum())
            if not np.any(usable):
                fold_row['reason'] = 'no finite non-zero training Up template'
                folds.append(fold_row)
                del training
                continue
            template = np.median(profiles[usable] / norm[usable, None], axis=0)
            design = np.column_stack((template, np.ones(180), time))
            if np.linalg.matrix_rank(design) < 3:
                fold_row['reason'] = 'template collinear with intercept and time'
                folds.append(fold_row)
                del training
                continue
            held_keep = keep[:, test_trials]
            acr, trial = np.nonzero(held_keep)
            values = aligned[:, test_trials, :][held_keep]
            fitted = np.linalg.lstsq(design, values.T, rcond=None)[0][0]
            template_coef[acr, test_trials[trial]] = fitted
            held_contrast[acr, test_trials[trial]] = contrast[acr, test_trials[trial]]
            fold_row['n_held_pairs'] = len(fitted)
            folds.append(fold_row)
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
            pooled[label].extend(values.tolist())
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
    fold_table = pd.DataFrame(folds)
    lines = ['', '## 12 November 2026: cross-validated template-fit SNR gain', '',
             'I used all 57 recordings, valid trials, run-onset frames and saved valid ACRs from analyse_acr_consistency_aligned.py. Its imported loader applies grid_response_test.postprocess_dff to the full grid of tdTomato-regressed dLight traces before selecting ACRs; no red term was added. Alignment calls align_trials and keep_aligned_trials over [-2, 4) s, 180 frames at 30 Hz, with frame 60 at run onset and the complete-window NaN rule.', '',
             'Kept trials vary between ACRs, so I defined shared recording-level folds by splitting the ordered trials with at least one kept ACR into five contiguous blocks using array_split. Every kept ACR-trial belongs to one held-out block, neighbouring trials do not interleave between folds, and each ACR retains its own NaN mask within those blocks. Selected trials with no kept ACR are omitted from fold construction. Fold sizes differ by at most one, with the earlier blocks receiving the remainder.', '',
             'Each fold calls quantify_grid_response on training frames only, with 1,000 shuffles and seed 42, and classify_up applies the validated 95th-percentile trial-averaged shuffle threshold and effect size > 0.05. Saved Up flags do not select template ACRs. Each training Up mean trace has its frames 30:60 mean subtracted, is divided by its L2 norm over frames 60:180, and contributes to the framewise median over all 180 frames; the median template is not renormalised. Non-finite or zero-norm training profiles cannot define a template and are counted in the fold table.', '',
             'For every held-out kept ACR-trial, ordinary least squares fits the full 180-frame trace against [template, 1, t], with t=(frame - 60)/30 seconds; A is the template template_coef. The same trial supplies dff_change from trial_amplitudes, mean frames 60:90 minus mean frames 30:60. Per ACR, all successful held-out estimates are combined before SNR_A and SNR_d are computed as mean / sample SD (ddof=1), then gain=SNR_A/SNR_d. ACRs with finite SNR_d <= 0 are excluded and counted; fewer than two held-out estimates, zero SDs or other non-finite gains are also excluded. A negative SNR_A is retained if SNR_d > 0. Exclusion columns for different causes can overlap; n_excluded_total is the unduplicated count.', '',
             'Saved Up and other valid ACRs define the descriptive reporting groups only. Pooled ACR medians give every included ACR equal weight; medians of recording medians give every recording equal weight, and both are reported to distinguish these two meanings of pooled. IQRs describe the corresponding distribution; no statistical tests were requested or run.', '']
    for label in ('Up', 'other'):
        group = table.loc[table.group == label]
        reported = group.loc[np.isfinite(group.median_gain)]
        values = np.asarray(pooled[label])
        q25, median, q75 = np.percentile(values, [25, 50, 75])
        r25, rmedian, r75 = np.percentile(reported.median_gain, [25, 50, 75])
        lines.append(f'{label}: pooled ACR gain median={median:.9g}, IQR [{q25:.9g}, {q75:.9g}], n={len(values)} ACRs from n={len(reported)} recordings, n={reported.animal.nunique()} animals. Median of recording medians={rmedian:.9g}, recording IQR [{r25:.9g}, {r75:.9g}]. Excluded finite SNR_d <= 0: {int(group.n_excluded_snr_d_nonpositive.sum())}; fewer than two held-out estimates: {int(group.n_fewer_two_held.sum())}; zero SD_A: {int(group.n_zero_sd_a.sum())}; zero SD_d: {int(group.n_zero_sd_d.sum())}; total excluded ACRs: {int(group.n_excluded_total.sum())} of {int(group.n_valid_acr.sum())}.')
    lines += ['', f'Skipped folds: {int((fold_table.reason != "").sum())} of {len(fold_table)}; no training Up: {int((fold_table.reason == "no training Up ACR").sum())}. All-missing selected trials omitted from fold construction: {int(table.loc[table.group == "Up", "n_all_missing_trials"].sum())}.', '',
              '| recording | group | ACRs with gain / valid | median gain | gain IQR | excluded SNR_d <= 0 | total excluded | folds used |',
              '| --- | --- | ---: | ---: | --- | ---: | ---: | ---: |']
    for row in table.itertuples():
        lines.append(f'| {row.rec_id} | {row.group} | {row.n_gain} / {row.n_valid_acr} | {row.median_gain:.9g} | [{row.gain_q25:.9g}, {row.gain_q75:.9g}] | {row.n_excluded_snr_d_nonpositive} | {row.n_excluded_total} | {row.n_used_folds} |')
    lines += ['', 'Training-only fold template counts:', '',
              '| recording | fold | training trials | held-out trials | training Up ACRs | template ACRs | held-out ACR-trials fitted | exclusion |',
              '| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |']
    for row in fold_table.itertuples():
        lines.append(f'| {row.rec_id} | {row.fold} | {row.n_train} | {row.n_test} | {row.n_up} | {row.n_template} | {row.n_held_pairs} | {row.reason or "none"} |')
    lines += ['', 'Table: outputs/da_hotspots/template_fit_gain.csv. Training Up selection and template estimation use no held-out trials; full-recording preprocessing is retained exactly as specified by the existing analysis, and the saved reporting groups were originally defined from these recordings. Gain can be large when positive SNR_d is close to zero; no threshold beyond the requested SNR_d > 0 was applied. Fold templates retain their specified individual L2 scaling before the A coefficients are combined across folds. No figure or other analysis was produced.']
    with (OUTPUT_ROOT / 'PROGRESS.md').open('a', encoding='utf-8') as stream:
        stream.write(chr(10).join(lines) + chr(10))
    print(chr(10).join(lines[14:18]), flush=True)


#%% run
if __name__ == '__main__':
    with threadpool_limits(limits=1):
        main()
