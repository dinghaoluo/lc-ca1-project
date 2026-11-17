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

import re

import numpy as np
import pandas as pd
from scipy.ndimage import gaussian_filter1d
from scipy.stats import pearsonr, wilcoxon
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
    previous_log = (OUTPUT_ROOT / 'PROGRESS.md').read_text(encoding='utf-8')
    co_rows = []
    time = (np.arange(180) - 60) / 30
    rows, folds, pooled = [], [], {'Up': [], 'other': []}
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
            template = gaussian_filter1d(template, sigma=3, mode='nearest')
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
        co = co_rows[-1]
        print(f'{recname}: r_resp={co["r_resp"]:.6g}, r_base={co["r_base"]:.6g}, n_trials={co["n_resp_trials"]}; {co["status"]}', flush=True)
    table = pd.DataFrame(rows)
    table.to_csv(output, index=False)
    co_table = pd.DataFrame(co_rows)
    co_table.to_csv(co_output, index=False)
    fold_table = pd.DataFrame(folds)
    lines = ['', '## 13 November 2026: smoothed templates and trial co-recruitment', '',
             'I repeated estimate_template_fit_gain.py with its imported inputs, full-grid postprocessing, saved valid ACR mask, valid trials, 180-frame [-2, 4) s alignment and complete-window NaN rule unchanged. The five shared contiguous blocks use array_split on the ordered trials with at least one kept ACR. Training-only quantify_grid_response uses 1,000 shuffles and seed 42, with the same classify_up rule; training Up profiles are baseline-subtracted over frames 30:60 and individually L2-normalised over frames 60:180 before the median. The sole fit change is gaussian_filter1d(template, sigma=3, mode=nearest) after that median and before the rank check and ordinary least-squares fit on [template, 1, t]. The smoothed template is not renormalised. Held-out estimates, sample SD (ddof=1), SNR_d > 0 cutoff and other undefined-gain exclusions are unchanged.', '',
             'The pooled gain below is SNR_A / SNR_d across all included ACRs, with each ACR receiving equal weight, as in the earlier pooled result. The unsmoothed values are read from the existing 12 November log entry, not recomputed. IQRs describe ACR dispersion; no smoothed-versus-unsmoothed statistical test was requested or run.', '',
             '| group | smoothed pooled median | smoothed IQR | unsmoothed pooled median | unsmoothed IQR | ACRs | recordings | animals |',
             '| --- | ---: | --- | ---: | --- | ---: | ---: | ---: |']
    exclusion_notes = []
    for label in ('Up', 'other'):
        group = table.loc[table.group == label]
        reported = group.loc[np.isfinite(group.median_gain)]
        low, median, high = quantiles(pooled[label])
        match = re.search(r'^' + label + r': pooled ACR gain median=([^,]+), IQR \[([^,]+), ([^\]]+)\]', previous_log, flags=re.MULTILINE)
        old = [float(value) for value in match.groups()] if match else [np.nan] * 3
        lines.append(f'| {label} | {median:.9g} | [{low:.9g}, {high:.9g}] | {old[0]:.9g} | [{old[1]:.9g}, {old[2]:.9g}] | {len(pooled[label])} | {len(reported)} | {reported.animal.nunique()} |')
        exclusion_notes.append(f'{label} exclusions: finite SNR_d <= 0: {int(group.n_excluded_snr_d_nonpositive.sum())}; fewer than two held-out estimates: {int(group.n_fewer_two_held.sum())}; zero SD_A: {int(group.n_zero_sd_a.sum())}; zero SD_d: {int(group.n_zero_sd_d.sum())}; total excluded: {int(group.n_excluded_total.sum())} of {int(group.n_valid_acr.sum())} ACRs. Cause counts may overlap; total is unduplicated.')
        if match is None:
            exclusion_notes.append(f'{label}: earlier pooled unsmoothed log entry missing; comparison left blank.')
    lines += [''] + exclusion_notes
    lines += ['', f'Skipped training folds: {int((fold_table.reason != "").sum())} of {len(fold_table)}. All-missing valid trial columns omitted from fold construction: {int(table.loc[table.group == "Up", "n_all_missing_trials"].sum())}.', '',
              'Co-recruitment uses all valid trials, irrespective of fold success: dff_change=mean frames 60:90 minus mean frames 30:60, and n=mean frames 30:60 minus mean frames 0:30, with the same per-ACR complete-window mask. Each trial mean gives equal weight to its kept ACRs within the saved Up or other group; a trial without a kept ACR in either group cannot enter that correlation. Pearson correlations use finite paired trial means separately for response and baseline. Recordings with fewer than three paired trials, constant group means or non-finite Fisher z are logged and excluded where required; an exact perfect correlation is not clipped. Saved Up labels were selected using these recordings, and the composition of kept ACRs can vary between trials.', '']
    for label in ('resp', 'base'):
        group = co_table.loc[np.isfinite(co_table[f'r_{label}'])]
        low, median, high = quantiles(group[f'r_{label}'])
        lines.append(f'r_{label}: pooled recording median={median:.9g}, IQR [{low:.9g}, {high:.9g}], n={len(group)} recordings, n={group.animal.nunique()} animals.')
    paired = co_table.loc[np.isfinite(co_table.z_resp_minus_base)]
    differences = paired.z_resp_minus_base.to_numpy()
    if len(differences) >= 2:
        if np.all(differences == 0):
            statistic, pvalue = 0.0, 1.0
        else:
            test = wilcoxon(differences, alternative='two-sided')
            statistic, pvalue = float(test.statistic), float(test.pvalue)
    else:
        statistic, pvalue = np.nan, np.nan
    low, median, high = quantiles(differences)
    lines += [f'Fisher-z(r_resp) minus Fisher-z(r_base): median difference={median:.9g}, IQR [{low:.9g}, {high:.9g}]; two-sided Wilcoxon signed-rank W={statistic:.9g}, p={pvalue:.12g}, n={len(paired)} recordings, n={paired.animal.nunique()} animals (zero_method=wilcox, method=auto).', '',
              '| recording | saved Up | other | paired response trials | r_resp | paired baseline trials | r_base | Fisher-z difference | status |',
              '| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |']
    for row in co_table.itertuples():
        reason = '; '.join(value for value in (row.status, row.resp_reason, row.base_reason) if value)
        lines.append(f'| {row.rec_id} | {row.n_up} | {row.n_other} | {row.n_resp_trials} | {row.r_resp:.9g} | {row.n_base_trials} | {row.r_base:.9g} | {row.z_resp_minus_base:.9g} | {reason} |')
    lines += ['', 'Per-recording smoothed gain:', '',
              '| recording | group | ACRs with gain / valid | median | IQR | excluded SNR_d <= 0 | total excluded | used folds |',
              '| --- | --- | ---: | ---: | --- | ---: | ---: | ---: |']
    for row in table.itertuples():
        lines.append(f'| {row.rec_id} | {row.group} | {row.n_gain} / {row.n_valid_acr} | {row.median_gain:.9g} | [{row.gain_q25:.9g}, {row.gain_q75:.9g}] | {row.n_excluded_snr_d_nonpositive} | {row.n_excluded_total} | {row.n_used_folds} |')
    lines += ['', 'Training-only fold counts:', '',
              '| recording | fold | training trials | held-out trials | training Up | usable template ACRs | fitted pairs | exclusion |',
              '| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |']
    for row in fold_table.itertuples():
        lines.append(f'| {row.rec_id} | {row.fold} | {row.n_train} | {row.n_test} | {row.n_up} | {row.n_template} | {row.n_held_pairs} | {row.reason or "none"} |')
    lines += ['', 'Outputs: outputs/da_hotspots/smoothed_template_gain.csv and outputs/da_hotspots/co_recruitment.csv. No figures or other analyses were produced; existing scripts and earlier outputs were left unchanged.']
    with (OUTPUT_ROOT / 'PROGRESS.md').open('a', encoding='utf-8') as stream:
        stream.write(chr(10).join(lines) + chr(10))
    print(chr(10).join(lines[:20]), flush=True)
    print(f'Co-recruitment: n={len(paired)}, W={statistic:.9g}, p={pvalue:.12g}', flush=True)


#%% run
if __name__ == '__main__':
    with threadpool_limits(limits=1):
        main()
