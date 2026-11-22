'''
Created on Sat Nov 21 2026

check whether Up-ACR fraction correlations follow a constant-response model

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
from scipy.stats import spearmanr, wilcoxon
from threadpoolctl import threadpool_limits

from analyse_acr_consistency_aligned import load_recording, event_frames
from grid_response_test import quantify_grid_response


#%% constants
OUTPUT_ROOT = pp.DLIGHT_HOTSPOTS_OUTPUT_ROOT
N_SHUFFLE = 1000
SEED = 42


#%% summaries
def safe_spearman(left, right):
    finite = np.isfinite(left) & np.isfinite(right)
    n = int(finite.sum())
    if n < 3 or np.ptp(left[finite]) == 0 or np.ptp(right[finite]) == 0:
        return np.nan, n
    return float(spearmanr(left[finite], right[finite]).statistic), n


def fisher_z(value):
    return float(np.arctanh(value)) if np.isfinite(value) and abs(value) < 1 else np.nan


def trial_summaries(change, responding, keep, up):
    valid = keep[up]
    values = np.where(valid, change[up], np.nan)
    responses = responding[up] & valid
    n_valid = valid.sum(axis=0)
    mean = np.divide(np.nansum(values, axis=0), n_valid,
                     out=np.full(values.shape[1], np.nan), where=n_valid > 0)
    fraction = np.divide(responses.sum(axis=0), n_valid,
                         out=np.full(values.shape[1], np.nan), where=n_valid > 0)
    return mean, fraction, n_valid


def constant_model_correlations(change, shuffled, keep, threshold, pooled_std, up):
    n_trials = change.shape[1]
    n_up = int(up.sum())
    model_values = np.full((n_up, n_trials, N_SHUFFLE), np.nan)
    model_responses = np.zeros((n_up, n_trials, N_SHUFFLE), dtype=bool)
    model_keep = np.zeros((n_up, n_trials), dtype=bool)
    for index, acr in enumerate(np.flatnonzero(up)):
        valid = keep[acr]
        observed_mean = np.mean(change[acr, valid]) if np.any(valid) else np.nan
        model_values[index, valid, :] = shuffled[acr, valid, :] + observed_mean
        model_keep[index, valid] = True
        effect = np.divide(model_values[index], pooled_std[acr],
                           out=np.full((n_trials, N_SHUFFLE), np.nan),
                           where=np.isfinite(pooled_std[acr]) & (pooled_std[acr] != 0))
        model_responses[index] = (
            (model_values[index] > threshold[acr])
            & (effect > 0.05)
            & model_keep[index, :, None]
        )

    correlations = np.full(N_SHUFFLE, np.nan)
    for shuffle in range(N_SHUFFLE):
        values = model_values[:, :, shuffle]
        valid = model_keep
        count = valid.sum(axis=0)
        mean = np.divide(np.nansum(values, axis=0), count,
                         out=np.full(n_trials, np.nan), where=count > 0)
        fraction = np.divide(model_responses[:, :, shuffle].sum(axis=0), count,
                             out=np.full(n_trials, np.nan), where=count > 0)
        correlations[shuffle] = safe_spearman(mean, fraction)[0]
    return correlations


#%% one recording
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

    response_mean, response_fraction, response_n = trial_summaries(
        change, responding, keep, up)
    baseline_change = (np.nanmean(result['aligned'][:, :, 30:60], axis=2)
                       - np.nanmean(result['aligned'][:, :, 0:30], axis=2))
    baseline_effect = np.divide(
        baseline_change, pooled_std[:, None],
        out=np.full(baseline_change.shape, np.nan),
        where=np.isfinite(pooled_std[:, None]) & (pooled_std[:, None] != 0))
    baseline_responding = (
        (baseline_change > threshold[:, None])
        & (baseline_effect > 0.05)
        & keep
    )
    baseline_mean, baseline_fraction, baseline_n = trial_summaries(
        baseline_change, baseline_responding, keep, up)

    observed_response, n_response = safe_spearman(response_mean, response_fraction)
    observed_baseline, n_baseline = safe_spearman(baseline_mean, baseline_fraction)
    constant = constant_model_correlations(
        change, shuffled, keep, threshold, pooled_std, up)
    constant = constant[np.isfinite(constant)]
    constant_q25, constant_median, constant_q75 = (
        np.percentile(constant, [25, 50, 75]) if len(constant)
        else (np.nan, np.nan, np.nan)
    )
    row = dict(
        rec_id=recname,
        animal=recname.split('-')[0],
        n_trials=len(frames),
        n_valid_acr=len(profile),
        n_up_acr=int(up.sum()),
        n_kept_pairs=int(keep.sum()),
        n_response_trials=n_response,
        n_baseline_trials=n_baseline,
        rho_observed_response=observed_response,
        rho_observed_baseline=observed_baseline,
        rho_constant_q25=constant_q25,
        rho_constant_median=constant_median,
        rho_constant_q75=constant_q75,
        n_constant_shuffles=len(constant),
    )
    row['z_observed_response'] = fisher_z(observed_response)
    row['z_constant_median'] = fisher_z(constant_median)
    row['z_observed_minus_constant'] = (
        row['z_observed_response'] - row['z_constant_median']
    )
    return row


#%% report
def append_log(table):
    finite_response = table.rho_observed_response.notna()
    finite_baseline = table.rho_observed_baseline.notna()
    finite_constant = table.rho_constant_median.notna()
    paired = table[['animal', 'z_observed_response',
                    'z_constant_median']].dropna()
    difference = paired.z_observed_response.to_numpy() - paired.z_constant_median.to_numpy()
    test = wilcoxon(difference, alternative='two-sided') if len(difference) else None

    lines = [
        '',
        '## 21 November 2026: constant-response mechanism for Up-ACR fractions',
        '',
        f'I used the saved Up flags, all {table.n_trials.sum()} valid trials from '
        f'{len(table)} recordings in {table.animal.nunique()} animals, with Jingyu\'s '
        f'full-window NaN rule, the 95th-percentile per-ACR shuffle threshold and '
        f'effect size above 0.05. For each Up ACR and each of the 1,000 shuffles, '
        f'the constant model adds that ACR\'s observed response mean to its shuffled '
        f'single-trial value, exactly as in the noise model comparison.',
        '',
    ]
    for column, label, finite in (
            ('rho_observed_response', 'observed response', finite_response),
            ('rho_observed_baseline', 'observed baseline', finite_baseline),
            ('rho_constant_median', 'constant model median', finite_constant)):
        values = table.loc[finite, column].to_numpy()
        low, median, high = np.percentile(values, [25, 50, 75]) if len(values) else (np.nan,) * 3
        lines.append(
            f'{label}: median rho={median:.9g}, IQR [{low:.9g}, {high:.9g}]; '
            f'n={len(values)} recordings, n={table.loc[finite, "animal"].nunique()} animals.'
        )
    low, median, high = np.percentile(difference, [25, 50, 75]) if len(difference) else (np.nan,) * 3
    lines += [
        f'Fisher-z(observed response) minus Fisher-z(constant model): median={median:.9g}, '
        f'IQR [{low:.9g}, {high:.9g}]; two-sided paired Wilcoxon '
        f'W={test.statistic if test else np.nan:.9g}, '
        f'p={test.pvalue if test else np.nan:.12g}; '
        f'n={len(paired)} recordings, n={paired.animal.nunique()} animals.',
        '',
        'Per-recording results: outputs/da_hotspots/fraction_mechanism.csv.',
        '',
    ]
    with (OUTPUT_ROOT / 'PROGRESS.md').open('a', encoding='utf-8') as stream:
        stream.write('\n'.join(lines))


#%% run
def main():
    output = OUTPUT_ROOT / 'fraction_mechanism.csv'
    if output.exists():
        raise FileExistsError(output)
    rows = []
    for recname in rec_list.recDLightDbh:
        row = recording_summary(recname)
        rows.append(row)
        print(f'{recname}: observed response rho={row["rho_observed_response"]:.6g}, '
              f'constant median={row["rho_constant_median"]:.6g}', flush=True)
    table = pd.DataFrame(rows)
    table.to_csv(output, index=False)
    append_log(table)


if __name__ == '__main__':
    with threadpool_limits(limits=1):
        main()
