# -*- coding: utf-8 -*-
'''
Created on Wed Nov 11 15:14:00 2026

early and late amplitude distributions and single-trial predictions

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
from scipy.stats import wilcoxon
from threadpoolctl import threadpool_limits

import compare_early_late_aligned as condition
import analyse_acr_consistency_aligned as consistency
from grid_response_test import (align_trials, keep_aligned_trials,
                                trial_amplitudes, quantify_grid_response)

repo_root   = Path(__file__).resolve().parents[3]
OUTPUT_ROOT = pp.DLIGHT_HOTSPOTS_OUTPUT_ROOT
FIGURE_ROOT = pp.DLIGHT_HOTSPOTS_FIGURE_ROOT
N_SHUFFLE   = 1000
SEED        = 42
THRESHOLDS  = (90, 95, 99)


#%% early and late amplitudes
def early_late_figure(values, recordings, variant, low, high):
    n_rec = len(recordings)
    ncol = 3 if variant == 'speed_matched' else 6
    nrow = int(np.ceil(n_rec / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(ncol * 2.65, nrow * 2.15),
                             sharex=True, sharey=True, squeeze=False)
    for i, recname in enumerate(recordings.rec_id):
        ax = axes.flat[i]
        panel = values.loc[values.rec_id == recname]
        for group, colour in (('early', '0.5'), ('late', '#19a819')):
            x = np.sort(panel.loc[panel.group == group, 'amplitude'].dropna().to_numpy())
            ax.step(np.r_[min(low, x[0]), x, max(high, x[-1])], np.r_[0, np.arange(1, len(x) + 1) / len(x), 1],
                    where='post', color=colour, linewidth=1.2, label=group)
        ax.set_title(recname, fontsize=9)
        ax.set_xlim(low, high)
        ax.set_ylim(0, 1)
        ax.spines[['top', 'right']].set_visible(False)
        ax.tick_params(labelsize=8)
        if i + ncol >= n_rec:
            ax.set_xlabel(chr(916) + 'dF/F')
            ax.tick_params(labelbottom=True)
    for i in range(n_rec, nrow * ncol):
        axes.flat[i].set_visible(False)
    fig.supylabel('Cumulative fraction of ACRs', fontsize=11)
    handles, labels = axes.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper center', ncol=2, frameon=False)
    fig.tight_layout(rect=(0.025, 0, 1, 0.965 if nrow == 3 else 0.985))
    filename = f'early_late_ecdf_{variant}.png'
    fig.savefig(FIGURE_ROOT / filename, dpi=200, facecolor='white')
    plt.close(fig)
    return filename

def early_late():
    saved = pd.read_parquet(OUTPUT_ROOT / 'aligned_early_late_per_recording.parquet').set_index('rec_id')
    amplitudes, counts, selections = [], [], []
    notes = []
    for recname in list(condition.recnames):
        try:
            behaviour, profile, traces = condition.load_recording(recname)
        except (OSError, KeyError, ValueError) as error:
            notes.append(f'{recname}: early/late input unavailable: {error}')
            continue
        valid, early, late, _ = condition.valid_trials(behaviour)
        vectors = condition.speed_vectors(behaviour, len(valid))
        matched_early, matched_late, reason = condition.speed_match(early, late, vectors)
        run = np.asarray(behaviour['run_onset_frames'], dtype=int)
        for variant, first, second, exclusion in (
                ('speed_matched', matched_early, matched_late, reason),
                ('unmatched', early, late, '' if early and late else 'no early or late trials')):
            included = not exclusion
            saved_column = 'speed_included' if variant == 'speed_matched' else 'unmatched_included'
            if included != bool(saved.loc[recname, saved_column]):
                notes.append(f'{recname}, {variant}: imported selection and saved inclusion disagree; imported selection retained.')
            row = dict(rec_id=recname, animal=recname.split('-')[0], variant=variant,
                       n_early_available=len(first), n_late_available=len(second),
                       included=included, exclusion=exclusion, n_early_used=0, n_late_used=0)
            if included:
                n = min(len(first), len(second))
                rng = np.random.RandomState(SEED)
                first_draw = list(first) if len(first) <= len(second) else list(rng.choice(first, n, replace=False))
                second_draw = list(second) if len(second) <= len(first) else list(rng.choice(second, n, replace=False))
                row.update(n_early_used=n, n_late_used=n)
                for group, trials in (('early', first_draw), ('late', second_draw)):
                    cleaned, keep = keep_aligned_trials(align_trials(traces, run[trials]))
                    trial_amp, _, _ = trial_amplitudes(cleaned)
                    amplitude = np.nanmean(trial_amp, axis=1)
                    row[f'n_{group}_missing_amplitude'] = int(np.sum(~np.isfinite(amplitude)))
                    for roi, value, kept in zip(profile['roi_id'], amplitude, keep.sum(axis=1)):
                        amplitudes.append(dict(rec_id=recname, animal=recname.split('-')[0],
                                               variant=variant, group=group, roi_y=int(roi[0]),
                                               roi_x=int(roi[1]), amplitude=value, n_kept=int(kept)))
                    for trial in trials:
                        selections.append(dict(rec_id=recname, variant=variant, group=group,
                                               trial_index=int(trial), run_onset_frame=int(run[trial])))
            counts.append(row)
        print(f'ECDF {recname}: {len(early)}/{len(late)} unmatched, {len(matched_early)}/{len(matched_late)} speed, {reason}', flush=True)
    amplitudes = pd.DataFrame(amplitudes)
    counts = pd.DataFrame(counts)
    amplitudes.to_csv(OUTPUT_ROOT / 'followup_early_late_amplitudes.csv', index=False)
    counts.to_csv(OUTPUT_ROOT / 'followup_early_late_counts.csv', index=False)
    pd.DataFrame(selections).to_csv(OUTPUT_ROOT / 'followup_early_late_trials.csv', index=False)
    lines = ['', '## 11 November 2026: early and late amplitude ECDFs', '',
             'I used the imported early/late selection and speed-matching functions, then retained the smaller group in full and drew once without replacement from the larger group. A fresh RandomState(42) was used for each recording and variant; resetting between variants resolves the unspecified ordering of the two comparisons. Equal-sized groups were both retained in full, without consuming a draw.',
             '', 'Amplitudes are response [0, 1] s minus baseline [-1, 0] s around run onset, in dF/F units, averaged over complete finite [-2, 4] s windows separately for each saved valid ACR. The early/late selection retains its original treatment of run-onset -1; those short windows fail the imported NaN rule. No shuffles or statistical tests were used.', '']
    for variant in ('speed_matched', 'unmatched'):
        values = amplitudes.loc[amplitudes.variant == variant]
        finite = values.amplitude.to_numpy()
        finite = finite[np.isfinite(finite)]
        low, high = np.percentile(finite, [1, 99])
        low, high = np.floor(low / 0.05) * 0.05, np.ceil(high / 0.05) * 0.05
        recordings = counts.loc[(counts.variant == variant) & counts.included]
        filename = early_late_figure(values, recordings, variant, low, high)
        outside = 100 * np.mean((finite < low) | (finite > high))
        missing = int(values.amplitude.isna().sum())
        lines.append(f'{variant}: n={len(recordings)} recordings, n={recordings.animal.nunique()} animals; shared x-limits [{low:.2f}, {high:.2f}], rounded outward from the pooled 1st and 99th percentiles in increments of 0.05. Of {len(finite)} finite plotted amplitudes, {outside:.4f}% lie outside these limits; {missing} ACR-group amplitudes are missing because no trial survived the NaN rule. The ECDF denominator includes all finite values, including those outside the displayed range. Figure: `figures/dlight_hotspots/{filename}`.')
    lines += ['', 'Trial counts used in each group, before the per-ACR NaN rule:', '',
              '| recording | variant | early | late | missing early amplitudes | missing late amplitudes |',
              '| --- | --- | ---: | ---: | ---: | ---: |']
    for row in counts.loc[counts.included].itertuples():
        lines.append(f'| {row.rec_id} | {row.variant} | {row.n_early_used} | {row.n_late_used} | {int(row.n_early_missing_amplitude)} | {int(row.n_late_missing_amplitude)} |')
    lines += ['', 'Tables: `followup_early_late_counts.csv`, `followup_early_late_trials.csv` (zero-based behaviour indices), and `followup_early_late_amplitudes.csv`, all under `outputs/da_hotspots/`.', '',
              'Ambiguities or missing inputs: ' + ('; '.join(notes) if notes else 'no missing inputs or inclusion differences; the seed reset between variants is stated above.')]
    with (OUTPUT_ROOT / 'PROGRESS.md').open('a', encoding='utf-8') as stream:
        stream.write(chr(10).join(lines) + chr(10))
    print('ECDF figures and log complete', flush=True)


#%% injected single-trial responses
def predictions():
    check = pd.read_csv(OUTPUT_ROOT / 'aligned_check_overlap.csv', float_precision='round_trip').set_index(['rec_id', 'threshold'])
    rows, draws, notes = [], [], []
    for recname in list(consistency.recnames):
        try:
            behaviour, profile, traces = consistency.load_recording(recname)
        except (OSError, KeyError, ValueError) as error:
            notes.append(f'{recname}: prediction input unavailable: {error}')
            continue
        frames = consistency.event_frames(behaviour)
        result = quantify_grid_response(traces, frames, n_shuffle=N_SHUFFLE, seed=SEED, return_single_trial=True)
        keep = result['keep']
        kept = keep.sum(axis=1)
        up = profile['Up'].to_numpy(dtype=bool)
        injected = result['shuffle_single_trial_amplitude'] + result['amplitude'][:, None, None]
        injected_effect_ok = injected / result['pooled_std'][:, None, None] > 0.05
        for q in THRESHOLDS:
            threshold, observed, circular = consistency.response_matrices(result, q)
            del circular
            predicted = (injected > threshold[:, None, None]) & injected_effect_ok & keep[:, :, None]
            rate_draws = predicted.sum(axis=(0, 1)) / keep.sum()
            observed_rate = observed.sum() / keep.sum()
            probability = np.divide(observed.sum(axis=1), kept, out=np.full(len(kept), np.nan), where=kept > 0)
            predicted_probability = np.divide(predicted.sum(axis=1).mean(axis=1), kept,
                                               out=np.full(len(kept), np.nan), where=kept > 0)
            overlap, overlap_draws, n_pairs = consistency.jaccard_summary(observed, predicted, keep)
            row = dict(rec_id=recname, animal=recname.split('-')[0], threshold=q, n_valid_acr=len(profile),
                       n_trials=len(frames), n_kept_pairs=int(keep.sum()), n_zero_kept_acr=int(np.sum(kept == 0)),
                       observed_rate=observed_rate, predicted_rate=float(np.mean(rate_draws)),
                       predicted_rate_low=float(np.percentile(rate_draws, 2.5)),
                       predicted_rate_high=float(np.percentile(rate_draws, 97.5)),
                       observed_overlap=overlap, predicted_overlap=float(np.nanmedian(overlap_draws)),
                       predicted_overlap_low=float(np.nanpercentile(overlap_draws, 2.5)),
                       predicted_overlap_high=float(np.nanpercentile(overlap_draws, 97.5)),
                       n_observed_pairs=n_pairs, n_finite_predicted_overlap=int(np.isfinite(overlap_draws).sum()))
            for label, group in (('up', up), ('other', ~up)):
                valid_group = group & (kept > 0)
                row[f'n_{label}_acr'] = int(group.sum())
                row[f'n_{label}_acr_with_trials'] = int(valid_group.sum())
                row[f'observed_{label}_probability'] = float(np.mean(probability[valid_group])) if valid_group.any() else np.nan
                row[f'predicted_{label}_probability'] = float(np.mean(predicted_probability[valid_group])) if valid_group.any() else np.nan
            if (recname, q) in check.index:
                reference = check.loc[(recname, q)]
                for label in ('median', 'low', 'high'):
                    row[f'permutation_null_{label}'] = reference[f'permutation_null_{label}']
                row['aligned_check_observed_overlap_difference'] = overlap - reference.observed_overlap
                if not np.isclose(overlap, reference.observed_overlap, rtol=0, atol=1e-12):
                    notes.append(f'{recname}, {q}th threshold: observed Jaccard differs from check by {overlap - reference.observed_overlap:.9g}.')
            else:
                notes.append(f'{recname}, {q}th threshold: count-preserving permutation reference missing.')
                for label in ('median', 'low', 'high'):
                    row[f'permutation_null_{label}'] = np.nan
            rows.append(row)
            draws.append(pd.DataFrame(dict(rec_id=recname, threshold=q, shuffle=np.arange(N_SHUFFLE),
                                           predicted_rate=rate_draws, predicted_overlap=overlap_draws)))
            print(f'Prediction {recname}, {q}: rate {observed_rate:.6f} -> {row["predicted_rate"]:.6f}; overlap {overlap:.6f} -> {row["predicted_overlap"]:.6f}', flush=True)
            del predicted
        pd.DataFrame(rows).to_csv(OUTPUT_ROOT / 'followup_predicted_vs_observed.csv', index=False)
        del result, injected, injected_effect_ok
    table = pd.DataFrame(rows)
    pd.concat(draws, ignore_index=True).to_csv(OUTPUT_ROOT / 'followup_prediction_draws.csv', index=False)
    tests = []
    for q in THRESHOLDS:
        subset = table.loc[table.threshold == q]
        for measure in ('rate', 'overlap'):
            pairs = subset.loc[np.isfinite(subset[f'observed_{measure}']) & np.isfinite(subset[f'predicted_{measure}'])]
            difference = pairs[f'predicted_{measure}'].to_numpy() - pairs[f'observed_{measure}'].to_numpy()
            if np.all(difference == 0):
                statistic, p_value = 0.0, 1.0
            else:
                test = wilcoxon(difference, alternative='two-sided')
                statistic, p_value = float(test.statistic), float(test.pvalue)
            nonzero = pairs[f'observed_{measure}'] != 0
            relative = difference[nonzero] / pairs.loc[nonzero, f'observed_{measure}'].to_numpy()
            tests.append(dict(threshold=q, measure=measure, n_recordings=len(pairs), n_animals=pairs.animal.nunique(),
                              observed_median=pairs[f'observed_{measure}'].median(),
                              predicted_median=pairs[f'predicted_{measure}'].median(),
                              median_difference=np.median(difference), median_relative_difference=np.median(relative),
                              n_relative=len(relative), n_zero_observed=int((~nonzero).sum()), W=statistic, p=p_value))
    tests = pd.DataFrame(tests)
    tests.to_csv(OUTPUT_ROOT / 'followup_prediction_tests.csv', index=False)
    prediction_figure(table)
    lines = ['', '## 11 November 2026: single-trial predictions from mean response amplitude', '',
             f'I used the imported consistency trial set and CPU test for {table.rec_id.nunique()} recordings from {table.animal.nunique()} animals, retaining the complete-window NaN rule and 1,000 circular shuffles with seed 42 reset for each recording. Each shuffled single-trial amplitude received the observed per-ACR trial-averaged amplitude, in the same dF/F units. The observed response matrix, existing per-ACR percentile threshold and pooled SD were unchanged; an injected response had to exceed that threshold and have injected amplitude / pooled SD > 0.05.', '',
             'The overall responding rate is total responding ACR-trials / total kept ACR-trials, with prediction averaged over the 1,000 injected shuffle rates. Group probabilities first divide responding trials by kept trials for each ACR, then average equally across ACRs, separately for saved Up and other valid ACRs. ACRs with no kept trials are excluded from these group means. The saved Up result is descriptive only: its mean amplitudes were selected to be large, so its prediction is an upper bound.', '',
             'Jaccard overlap calls the existing function, restricting each trial pair to ACRs kept on both trials and selecting non-empty unions separately in each responding matrix. Prediction is the median of the 1,000 injected overlaps, with a 2.5th to 97.5th percentile band; the count-preserving permutation median and band are copied from aligned_check_overlap.csv alongside it.', '',
             'Paired two-sided Wilcoxon tests across recordings use predicted minus observed. Relative difference is (predicted - observed) / observed, summarised by its median; zero observed values are omitted only from this relative summary. Recording medians and median paired differences are listed below; rate differences are percentage points, and overlap differences are Jaccard units.', '']
    for row in tests.itertuples():
        scale = 100 if row.measure == 'rate' else 1
        lines.append(f'{row.threshold}th threshold, {row.measure}: observed recording median={scale * row.observed_median:.6g}, predicted recording median={scale * row.predicted_median:.6g}; median paired difference={scale * row.median_difference:.6g}, median relative difference={100 * row.median_relative_difference:.4f}%; W={row.W:g}, p={row.p:.9g}, n={row.n_recordings} recordings, n={row.n_animals} animals. Relative summary n={row.n_relative}; zero observed omitted={row.n_zero_observed}.')
    lines += ['', 'Descriptive group probabilities (mean across ACRs within each recording, then recording medians; all values in percent):', '',
              '| threshold | Up observed | Up predicted | other observed | other predicted |',
              '| ---: | ---: | ---: | ---: | ---: |']
    for q in THRESHOLDS:
        subset = table.loc[table.threshold == q]
        values = [100 * subset[f'{kind}_{group}_probability'].median() for group in ('up', 'other') for kind in ('observed', 'predicted')]
        lines.append(f'| {q} | ' + ' | '.join(f'{value:.5f}' for value in values) + ' |')
    main = table.loc[table.threshold == 95]
    pooled_obs = np.average(main.observed_rate, weights=main.n_kept_pairs)
    pooled_pred = np.average(main.predicted_rate, weights=main.n_kept_pairs)
    lines += ['', f'At the 95th threshold, the pooled observed responding rate is {100 * pooled_obs:.5f}% and the predicted mean rate is {100 * pooled_pred:.5f}% over {int(main.n_kept_pairs.sum())} kept ACR-trials and {int(main.n_valid_acr.sum())} valid ACRs. All per-recording rates, group probabilities, overlap bands and check permutation references are in followup_predicted_vs_observed.csv; the 1,000 injected values are retained in followup_prediction_draws.csv, and all six tests are in followup_prediction_tests.csv.', '',
              f'The largest absolute difference between recomputed observed overlap and the check is {table.aligned_check_observed_overlap_difference.abs().max():.9g}; {int(main.n_zero_kept_acr.sum())} ACRs have no kept trials. Finite predicted overlap counts range from {int(table.n_finite_predicted_overlap.min())} to {int(table.n_finite_predicted_overlap.max())} per recording and threshold.', '',
              'The prediction assumes the same response amplitude on every trial, so variation in response amplitude between trials is outside what this injection represents; the saved Up selection also favours large injected means, as stated above. The observed mean is estimated from the same trials being described, so these are descriptive predictions rather than held-out forecasts. No additional analyses or controls were run.', '',
              'Ambiguities or missing inputs: ' + ('; '.join(notes) if notes else 'none.'), '',
              'Figure: `figures/dlight_hotspots/single_trial_predicted_vs_observed.png`; the two panels show per-recording rate and overlap at the 95th threshold, with the identity line.']
    with (OUTPUT_ROOT / 'PROGRESS.md').open('a', encoding='utf-8') as stream:
        stream.write(chr(10).join(lines) + chr(10))
    print(tests.to_string(index=False), flush=True)


#%% observed and predicted figure
def prediction_figure(table):
    rows = table.loc[table.threshold == 95]
    fig, axes = plt.subplots(1, 2, figsize=(8.4, 4.0))
    colours = ('#19a819', '#4477aa', '#aa6677', '#aa8833')
    for ax, measure, label in zip(axes, ('rate', 'overlap'), ('responding rate (%)', 'mean Jaccard overlap')):
        scale = 100 if measure == 'rate' else 1
        for animal, colour in zip(rows.animal.unique(), colours):
            group = rows.loc[rows.animal == animal]
            ax.scatter(scale * group[f'observed_{measure}'], scale * group[f'predicted_{measure}'],
                       s=23, color=colour, alpha=0.8, label=animal)
        finite = scale * rows[[f'observed_{measure}', f'predicted_{measure}']].to_numpy()
        finite = finite[np.isfinite(finite)]
        upper = finite.max() * 1.06
        ax.plot([0, upper], [0, upper], color='0.5', linewidth=1, linestyle='--', label='identity')
        ax.set_xlim(0, upper)
        ax.set_ylim(0, upper)
        ax.set_aspect('equal', adjustable='box')
        ax.set_xlabel(f'Observed {label}')
        ax.set_ylabel(f'Predicted {label}')
        ax.spines[['top', 'right']].set_visible(False)
    axes[0].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(FIGURE_ROOT / 'single_trial_predicted_vs_observed.png', dpi=200, facecolor='white')
    plt.close(fig)


#%% run
if __name__ == '__main__':
    with threadpool_limits(limits=1):
        early_late()
        predictions()
