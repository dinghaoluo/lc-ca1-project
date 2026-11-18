'''
Created on Wed Nov 18 00:41:25 2026

compare response and baseline co-recruitment of responding ACR counts

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

OUTPUT_ROOT = pp.DLIGHT_HOTSPOTS_OUTPUT_ROOT
N_SHUFFLE   = 1000
SEED        = 42


#%% responding counts across trials
def count_correlation(up_count, other_count):
    if len(up_count) < 3:
        return np.nan, 'fewer than three trials'
    if np.ptp(up_count) == 0 or np.ptp(other_count) == 0:
        return np.nan, 'constant responding count'
    return float(spearmanr(up_count, other_count)[0]), ''


def recording_summary(recname):
    behaviour, profile, traces = load_recording(recname)
    frames = event_frames(behaviour)
    result = quantify_grid_response(traces, frames, n_shuffle=N_SHUFFLE,
                                    seed=SEED, return_single_trial=True)
    shuffled = result['shuffle_single_trial_amplitude']
    threshold = np.nanpercentile(shuffled.reshape(len(profile), -1), 95, axis=1)
    del shuffled
    del result['shuffle_single_trial_amplitude']
    baseline_change = (result['aligned'][:, :, 30:60].mean(axis=2)
                       - result['aligned'][:, :, 0:30].mean(axis=2))
    saved_up = profile['Up'].to_numpy(dtype=bool)
    row = dict(rec_id=recname, animal=recname.split('-')[0],
               n_trials=len(frames), n_valid_acr=len(profile), n_up=int(saved_up.sum()),
               n_other=int((~saved_up).sum()), n_kept_pairs=int(result['keep'].sum()),
               n_all_missing_trials=int((~result['keep'].any(axis=0)).sum()),
               rho_resp=np.nan, rho_base=np.nan, z_resp=np.nan, z_base=np.nan,
               z_resp_minus_base=np.nan, resp_reason='', base_reason='', status='ok')
    if not saved_up.any() or saved_up.all():
        row['status'] = 'no saved Up ACR' if not saved_up.any() else 'no other valid ACR'
        return row
    for label, change in (('resp', result['trial_amplitude']), ('base', baseline_change)):
        responding = ((change > threshold[:, None])
                      & (change / result['pooled_std'][:, None] > 0.05)
                      & result['keep'])
        up_count = responding[saved_up].sum(axis=0)
        other_count = responding[~saved_up].sum(axis=0)
        rho, reason = count_correlation(up_count, other_count)
        row[f'rho_{label}'] = rho
        row[f'{label}_reason'] = reason
        if np.isfinite(rho) and abs(rho) < 1:
            row[f'z_{label}'] = float(np.arctanh(rho))
        elif np.isfinite(rho):
            row[f'{label}_reason'] = 'perfect correlation; Fisher z not finite'
    row['z_resp_minus_base'] = row['z_resp'] - row['z_base']
    if not np.isfinite(row['z_resp_minus_base']):
        row['status'] = 'unavailable paired Fisher z'
    return row


#%% run and report
def main():
    output = OUTPUT_ROOT / 'count_co_recruitment.csv'
    if output.exists():
        raise FileExistsError(output)
    rows, failures = [], []
    for recname in rec_list.recDLightDbh:
        try:
            row = recording_summary(recname)
            print(f'{recname}: rho response={row["rho_resp"]:.6g}, baseline={row["rho_base"]:.6g}; {row["status"]}', flush=True)
        except (OSError, ValueError, KeyError) as error:
            failures.append(f'{recname}: {type(error).__name__}: {error}')
            row = dict(rec_id=recname, animal=recname.split('-')[0], status=failures[-1])
            print(failures[-1], flush=True)
        rows.append(row)
    table = pd.DataFrame(rows)
    table.to_csv(output, index=False)
    lines = ['', '## 18 November 2026: co-recruitment of responding ACR counts', '',
             'I used the consistency loader, valid run-onset trials, saved valid ACRs and saved Up flags unchanged, with 1,000 circular shifts and seed 42. The 95th-percentile per-ACR threshold pools single-trial shuffle amplitudes exactly as in response_matrices; the same threshold and original response-window pooled SD are applied to response mean frames 60:90 minus 30:60 and baseline mean frames 30:60 minus 0:30, with effect size >0.05 and the existing kept-trial mask. Spearman correlations compare saved Up and other responding counts across all selected trial columns, including zero counts where no ACR was kept; counts are not normalised by availability.', '']
    for label in ('resp', 'base'):
        group = table.loc[np.isfinite(table[f'rho_{label}'])]
        values = group[f'rho_{label}'].to_numpy()
        low, median, high = np.percentile(values, [25, 50, 75]) if len(values) else [np.nan] * 3
        lines.append(f'rho_{label}: median={median:.9g}, IQR [{low:.9g}, {high:.9g}]; n={len(group)} recordings, n={group.animal.nunique()} animals.')
    paired = table.loc[np.isfinite(table.z_resp_minus_base)]
    difference = paired.z_resp_minus_base.to_numpy()
    if len(difference) == 0:
        statistic, pvalue = np.nan, np.nan
    elif np.all(difference == 0):
        statistic, pvalue = 0.0, 1.0
    else:
        test = wilcoxon(difference, alternative='two-sided')
        statistic, pvalue = float(test.statistic), float(test.pvalue)
    low, median, high = np.percentile(difference, [25, 50, 75]) if len(difference) else [np.nan] * 3
    lines.append(f'Fisher-z(rho_resp) minus Fisher-z(rho_base): median={median:.9g}, IQR [{low:.9g}, {high:.9g}]; two-sided paired Wilcoxon W={statistic:.9g}, p={pvalue:.12g}; n={len(paired)} recordings, n={paired.animal.nunique()} animals. IQRs describe recording dispersion; perfect correlations are not clipped.')
    lines += ['', f'Attempted n={len(table)} recordings, n={table.animal.nunique()} animals.']
    for row in table.loc[table.status.ne('ok')].itertuples():
        lines.append(f'{row.rec_id}: {row.status}; response: {getattr(row, "resp_reason", "")}; baseline: {getattr(row, "base_reason", "")}.')
    lines += failures
    lines += ['', 'Per-recording results: outputs/da_hotspots/count_co_recruitment.csv.', '']
    with (OUTPUT_ROOT / 'PROGRESS.md').open('a', encoding='utf-8') as stream:
        stream.write(chr(10).join(lines))
    print(chr(10).join(lines[4:8]), flush=True)


if __name__ == '__main__':
    with threadpool_limits(limits=1):
        main()
