'''
Created on Sat Nov 21 2026

characterise single-trial dLight responses and noise-limited models

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
from scipy.stats import spearmanr, wilcoxon
from threadpoolctl import threadpool_limits

from analyse_acr_consistency_aligned import load_recording, event_frames
from grid_response_test import quantify_grid_response


#%% constants
OUTPUT_ROOT = pp.DLIGHT_HOTSPOTS_OUTPUT_ROOT
FIGURE_ROOT = pp.DLIGHT_HOTSPOTS_FIGURE_ROOT
N_SHUFFLE = 1000
SEED = 42
GROUPS = (('up', True), ('other', False))
MODELS = ('observed', 'constant', 'sparse_0p5', 'sparse_0p25', 'sparse_0p1')


#%% small summaries
def median_iqr(values):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if not len(values):
        return np.nan, np.nan, np.nan
    return tuple(np.percentile(values, [25, 50, 75]))


def safe_median(values):
    values = np.asarray(values, dtype=float)
    return float(np.nanmedian(values)) if np.any(np.isfinite(values)) else np.nan


def safe_spearman(left, right):
    finite = np.isfinite(left) & np.isfinite(right)
    if finite.sum() < 3:
        return np.nan, int(finite.sum()), 'fewer than three trials'
    if np.ptp(left[finite]) == 0 or np.ptp(right[finite]) == 0:
        return np.nan, int(finite.sum()), 'constant trial summary'
    return float(spearmanr(left[finite], right[finite]).statistic), int(finite.sum()), ''


def fisher_z(value):
    return float(np.arctanh(value)) if np.isfinite(value) and abs(value) < 1 else np.nan


def per_acr_statistics(values, shuffled, responding):
    values = np.asarray(values, dtype=float)
    finite = np.isfinite(values)
    observed = values[finite]
    mean = np.mean(observed) if len(observed) else np.nan
    sd = np.std(observed, ddof=1) if len(observed) >= 2 else np.nan
    p90 = np.percentile(observed, 90) if len(observed) else np.nan
    probability = np.mean(responding[finite]) if len(observed) else np.nan

    shuffled = np.asarray(shuffled, dtype=float)[finite]
    shuffle_sd = np.nanstd(shuffled, axis=0, ddof=1) if len(observed) >= 2 else np.full(N_SHUFFLE, np.nan)
    shuffle_p90 = np.nanpercentile(shuffled, 90, axis=0) if len(observed) else np.full(N_SHUFFLE, np.nan)
    sd_null = safe_median(shuffle_sd)
    p90_null = safe_median(shuffle_p90)
    sd_ratio = sd / sd_null if np.isfinite(sd) and np.isfinite(sd_null) and sd_null != 0 else np.nan
    p90_ratio = p90 / p90_null if np.isfinite(p90) and np.isfinite(p90_null) and p90_null != 0 else np.nan
    return dict(mean=mean, sd=sd, p90=p90, response_probability=probability,
                shuffle_sd_median=sd_null, shuffle_p90_median=p90_null,
                sd_ratio=sd_ratio, p90_ratio=p90_ratio, n_trials=int(len(observed)))


#%% trial summaries and noise models
def up_trial_summaries(change, responding, keep, up):
    valid = keep[up]
    values = np.where(valid, change[up], np.nan)
    responses = responding[up] & valid
    n_valid = valid.sum(axis=0)
    mean = np.divide(np.nansum(values, axis=0), n_valid,
                     out=np.full(values.shape[1], np.nan), where=n_valid > 0)
    median = np.nanmedian(values, axis=0)
    fraction = np.divide(responses.sum(axis=0), n_valid,
                         out=np.full(values.shape[1], np.nan), where=n_valid > 0)
    return mean, median, fraction, n_valid


def model_statistics(change, shuffle_single, keep, threshold, pooled_std, up, rng):
    rows = []
    acr_rows = []
    for acr in np.flatnonzero(up):
        valid = keep[acr]
        observed = change[acr, valid]
        noise = shuffle_single[acr, valid, :]
        if not len(observed):
            continue
        observed_stats = per_acr_statistics(
            change[acr], shuffle_single[acr],
            (change[acr] > threshold[acr])
            & (np.divide(change[acr], pooled_std[acr],
                         out=np.full(change.shape[1], np.nan),
                         where=np.isfinite(pooled_std[acr]) & (pooled_std[acr] != 0)) > 0.05)
        )
        acr_rows.append(dict(model='observed', acr=acr,
                             sd_ratio=observed_stats['sd_ratio'],
                             p90_ratio=observed_stats['p90_ratio'],
                             response_probability=observed_stats['response_probability']))
        signals = {
            'constant': np.full(noise.shape, observed.mean()),
            'sparse_0p5': (rng.random_sample(noise.shape) < 0.5) * (observed.mean() / 0.5),
            'sparse_0p25': (rng.random_sample(noise.shape) < 0.25) * (observed.mean() / 0.25),
            'sparse_0p1': (rng.random_sample(noise.shape) < 0.1) * (observed.mean() / 0.1),
        }
        for model, signal in signals.items():
            predicted = noise + signal
            sd = np.std(predicted, axis=0, ddof=1)
            p90 = np.percentile(predicted, 90, axis=0)
            ratio_sd = np.divide(sd, observed_stats['shuffle_sd_median'],
                                 out=np.full(N_SHUFFLE, np.nan),
                                 where=np.isfinite(observed_stats['shuffle_sd_median'])
                                 & (observed_stats['shuffle_sd_median'] != 0))
            ratio_p90 = np.divide(p90, observed_stats['shuffle_p90_median'],
                                  out=np.full(N_SHUFFLE, np.nan),
                                  where=np.isfinite(observed_stats['shuffle_p90_median'])
                                  & (observed_stats['shuffle_p90_median'] != 0))
            ratio = np.divide(predicted, pooled_std[acr],
                              out=np.full(predicted.shape, np.nan),
                              where=np.isfinite(pooled_std[acr]) & (pooled_std[acr] != 0))
            response_probability = np.mean(
                (predicted > threshold[acr]) & (ratio > 0.05), axis=0)
            acr_rows.append(dict(model=model, acr=acr,
                                 sd_ratio=safe_median(ratio_sd),
                                 p90_ratio=safe_median(ratio_p90),
                                 response_probability=safe_median(response_probability)))
    for model in MODELS:
        selected = [row for row in acr_rows if row['model'] == model]
        rows.append(dict(model=model, n_up_acr=len(selected),
                         median_sd_ratio=safe_median([row['sd_ratio'] for row in selected]),
                         median_p90_ratio=safe_median([row['p90_ratio'] for row in selected]),
                         median_response_probability=safe_median(
                             [row['response_probability'] for row in selected])))
    return rows, acr_rows


#%% recording analysis
def recording_summary(recname):
    behaviour, profile, traces = load_recording(recname)
    frames = event_frames(behaviour)
    result = quantify_grid_response(traces, frames, n_shuffle=N_SHUFFLE,
                                    seed=SEED, return_single_trial=True)
    keep = result['keep']
    change = result['trial_amplitude']
    shuffled = result['shuffle_single_trial_amplitude']
    pooled_std = result['pooled_std']
    threshold = np.nanpercentile(shuffled.reshape(len(profile), -1), 95, axis=1)
    effect = np.divide(change, pooled_std[:, None],
                       out=np.full(change.shape, np.nan),
                       where=np.isfinite(pooled_std[:, None]) & (pooled_std[:, None] != 0))
    responding = (change > threshold[:, None]) & (effect > 0.05) & keep
    up = profile['Up'].to_numpy(dtype=bool)

    point_rows = []
    row = dict(rec_id=recname, animal=recname.split('-')[0],
               n_trials=len(frames), n_valid_acr=len(profile),
               n_up_acr=int(up.sum()), n_other_acr=int((~up).sum()),
               n_kept_pairs=int(keep.sum()), n_all_missing_trials=int((~keep.any(axis=0)).sum()))
    acr_values = {}
    for label, group in GROUPS:
        selected = up == group
        stats = [per_acr_statistics(change[acr], shuffled[acr],
                                     responding[acr]) for acr in np.flatnonzero(selected)]
        for key in ('mean', 'sd', 'p90', 'response_probability', 'sd_ratio', 'p90_ratio'):
            values = [stat[key] for stat in stats]
            row[f'{key}_{label}_median'] = safe_median(values)
            row[f'n_{key}_{label}'] = int(np.sum(np.isfinite(values)))
        acr_values[label] = stats
        for acr in np.flatnonzero(selected):
            values = change[acr, keep[acr]]
            if len(values) >= 2:
                point_rows.append(dict(rec_id=recname, group=label,
                                       mean=np.mean(values),
                                       sd=np.std(values, ddof=1)))

    response_mean, response_median, response_fraction, response_n = up_trial_summaries(
        change, responding, keep, up)
    baseline_change = (np.nanmean(result['aligned'][:, :, 30:60], axis=2)
                       - np.nanmean(result['aligned'][:, :, 0:30], axis=2))
    baseline_response = ((baseline_change > threshold[:, None]) & (
        np.divide(baseline_change, pooled_std[:, None],
                  out=np.full(baseline_change.shape, np.nan),
                  where=np.isfinite(pooled_std[:, None]) & (pooled_std[:, None] != 0)) > 0.05)
                        & keep)
    base_mean, base_median, base_fraction, base_n = up_trial_summaries(
        baseline_change, baseline_response, keep, up)
    for label, means, medians, fractions, counts in (
            ('resp', response_mean, response_median, response_fraction, response_n),
            ('base', base_mean, base_median, base_fraction, base_n)):
        row[f'{label}_mean_median'] = safe_median(means)
        row[f'{label}_median_median'] = safe_median(medians)
        row[f'{label}_fraction_median'] = safe_median(fractions)
        rho, n, reason = safe_spearman(means, medians)
        row[f'rho_{label}_mean_median'] = rho
        row[f'n_{label}_mean_median'] = n
        row[f'{label}_mean_median_reason'] = reason
        rho, n, reason = safe_spearman(means, fractions)
        row[f'rho_{label}_mean_fraction'] = rho
        row[f'n_{label}_mean_fraction'] = n
        row[f'{label}_mean_fraction_reason'] = reason

    row['z_resp_mean_median'] = fisher_z(row['rho_resp_mean_median'])
    row['z_base_mean_median'] = fisher_z(row['rho_base_mean_median'])
    row['z_resp_mean_fraction'] = fisher_z(row['rho_resp_mean_fraction'])
    row['z_base_mean_fraction'] = fisher_z(row['rho_base_mean_fraction'])
    row['z_diff_mean_median'] = row['z_resp_mean_median'] - row['z_base_mean_median']
    row['z_diff_mean_fraction'] = row['z_resp_mean_fraction'] - row['z_base_mean_fraction']

    rng = np.random.RandomState(SEED)
    model_rows, model_acr_rows = model_statistics(
        change, shuffled, keep, threshold, pooled_std, up, rng)
    model_lookup = {entry['model']: entry for entry in model_rows}
    observed_sd = model_lookup['observed']['median_sd_ratio']
    observed_probability = model_lookup['observed']['median_response_probability']
    candidates = [model for model in MODELS if model != 'observed']
    row['closest_sd_model'] = min(
        candidates, key=lambda model: abs(model_lookup[model]['median_sd_ratio'] - observed_sd)
        if np.isfinite(model_lookup[model]['median_sd_ratio']) and np.isfinite(observed_sd) else np.inf)
    row['closest_response_probability_model'] = min(
        candidates, key=lambda model: abs(model_lookup[model]['median_response_probability'] - observed_probability)
        if np.isfinite(model_lookup[model]['median_response_probability']) and np.isfinite(observed_probability) else np.inf)
    for model in model_rows:
        model.update(rec_id=recname, animal=recname.split('-')[0],
                     closest_sd_model=row['closest_sd_model'],
                     closest_response_probability_model=row['closest_response_probability_model'])
    for entry in model_acr_rows:
        entry.update(rec_id=recname, animal=recname.split('-')[0])
    return row, model_rows, model_acr_rows, point_rows


#%% figure and report
def make_figure(points):
    path = FIGURE_ROOT / 'acr_mean_sd.png'
    if path.exists():
        raise FileExistsError(path)
    fig, ax = plt.subplots(figsize=(5.0, 4.0))
    for label, colour in (('up', '#228b22'), ('other', '0.5')):
        subset = points[points['group'] == label]
        ax.scatter(subset['mean'], subset['sd'], s=9, alpha=0.55,
                   color=colour, label=label.capitalize(), linewidths=0)
    ax.set_xlabel('single-trial dF/F change, mean')
    ax.set_ylabel('single-trial dF/F change, SD')
    ax.legend(frameon=False)
    ax.spines[['top', 'right']].set_visible(False)
    fig.tight_layout()
    fig.savefig(path, dpi=200, facecolor='white')
    plt.close(fig)


def append_log(table, model_table, model_acr, points):
    lines = [
        '', '## 21 November 2026: single-trial response characterisation', '',
        f'I used the saved Up flags, all {table.n_trials.sum()} valid trials from '
        f'{len(table)} recordings in {table.animal.nunique()} animals, with Jingyu\'s '
        f'full-window NaN rule, the 95th-percentile per-ACR shuffle threshold and '
        f'effect size above 0.05. The response window is mean frames 60:90 minus '
        f'baseline frames 30:60; the baseline comparison is mean frames 30:60 minus '
        f'frames 0:30.', '',
    ]
    for key, label in (('mean', 'mean'), ('sd', 'SD'), ('p90', '90th percentile'),
                       ('response_probability', 'response probability'),
                       ('sd_ratio', 'SD relative to shuffle median'),
                       ('p90_ratio', '90th percentile relative to shuffle median')):
        paired = table[['animal', f'{key}_up_median', f'{key}_other_median']].dropna()
        difference = paired.iloc[:, 1].to_numpy() - paired.iloc[:, 2].to_numpy()
        test = wilcoxon(difference, alternative='two-sided') if len(difference) else None
        low, median, high = median_iqr(difference)
        lines.append(f'{label}, Up minus Other: median difference={median:.9g}, IQR '
                     f'[{low:.9g}, {high:.9g}]; W={test.statistic if test else np.nan:.9g}, '
                     f'p={test.pvalue if test else np.nan:.12g}; n={len(paired)} recordings, '
                     f'n={paired.animal.nunique()} animals.')
    for label, text in (
        ('mean_median', 'mean against median'),
        ('mean_fraction', 'mean against response fraction')):
        for window, suffix in (('resp', 'response'), ('base', 'baseline')):
            values = table[f'rho_{window}_{label}'].to_numpy()
            finite = np.isfinite(values)
            low, median, high = median_iqr(values[finite])
            lines.append(f'{suffix} {text}: median rho={median:.9g}, IQR '
                         f'[{low:.9g}, {high:.9g}]; n={finite.sum()} recordings, '
                         f'n={table.loc[finite, "animal"].nunique()} animals.')
        paired = table[['animal', f'z_resp_{label}', f'z_base_{label}']].dropna()
        difference = paired.iloc[:, 1].to_numpy() - paired.iloc[:, 2].to_numpy()
        test = wilcoxon(difference, alternative='two-sided') if len(difference) else None
        low, median, high = median_iqr(difference)
        lines.append(f'Fisher-z response minus baseline, {text}: median={median:.9g}, '
                     f'IQR [{low:.9g}, {high:.9g}]; W={test.statistic if test else np.nan:.9g}, '
                     f'p={test.pvalue if test else np.nan:.12g}; n={len(paired)} recordings, '
                     f'n={paired.animal.nunique()} animals.')

    lines.append('')
    lines.append('Noise-limited and sparse models, pooled over Up ACRs:')
    for model, group in model_acr.groupby('model', sort=False):
        line = [f'{model}: n={len(group)}']
        for key, label in (('sd_ratio', 'SD ratio'), ('p90_ratio', '90th-percentile ratio'),
                           ('response_probability', 'response probability')):
            low, median, high = median_iqr(group[key])
            line.append(f'{label} median={median:.9g}, IQR [{low:.9g}, {high:.9g}]')
        lines.append('; '.join(line) + '.')
    for column, label in (('closest_sd_model', 'closest median SD ratio'),
                          ('closest_response_probability_model', 'closest median response probability')):
        counts = table[column].value_counts()
        lines.append(f'{label} by recording: ' + ', '.join(
            f'{model}={int(counts.get(model, 0))}' for model in MODELS[1:]) +
            f'; n={len(table)} recordings, n={table.animal.nunique()} animals.')
    lines += ['', 'Figure: figures/dlight_hotspots/acr_mean_sd.png.',
              'Per-recording summaries: outputs/da_hotspots/single_trial_characterisation.csv.',
              'Per-recording model summaries: outputs/da_hotspots/single_trial_models.csv.', '']
    with (OUTPUT_ROOT / 'PROGRESS.md').open('a', encoding='utf-8') as stream:
        stream.write('\n'.join(lines))


#%% run
def main():
    outputs = (OUTPUT_ROOT / 'single_trial_characterisation.csv',
               OUTPUT_ROOT / 'single_trial_models.csv')
    for path in outputs:
        if path.exists():
            raise FileExistsError(path)
    table_rows, model_rows, model_acr_rows, point_rows = [], [], [], []
    for recname in rec_list.recDLightDbh:
        row, models, model_acr, points = recording_summary(recname)
        table_rows.append(row)
        model_rows.extend(models)
        model_acr_rows.extend(model_acr)
        point_rows.extend(points)
        print(f'{recname}: {row["n_trials"]} trials, {row["n_valid_acr"]} ACRs', flush=True)
    table = pd.DataFrame(table_rows)
    model_table = pd.DataFrame(model_rows)
    model_acr = pd.DataFrame(model_acr_rows)
    points = pd.DataFrame(point_rows)
    table.to_csv(OUTPUT_ROOT / 'single_trial_characterisation.csv', index=False)
    model_table.to_csv(OUTPUT_ROOT / 'single_trial_models.csv', index=False)
    make_figure(points)
    append_log(table, model_table, model_acr, points)


if __name__ == '__main__':
    with threadpool_limits(limits=1):
        main()
