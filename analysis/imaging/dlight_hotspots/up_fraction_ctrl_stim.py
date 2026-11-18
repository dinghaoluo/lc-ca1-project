'''
Created on Wed Nov 18 00:41:25 2026

compare ctrl and stim Up ACR fractions with equal trial counts

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

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon
from threadpoolctl import threadpool_limits

import compare_ctrl_stim_dlight as comparison
from grid_response_test import align_trials, pooled_std, quantify_grid_response, classify_up

OUTPUT_ROOT  = pp.DLIGHT_HOTSPOTS_OUTPUT_ROOT
SESSION_ROOT = OUTPUT_ROOT / 'ctrl_stim_sessions'
N_SHUFFLE    = 1000
N_DRAWS      = 20
SEED         = 42


#%% stimulation windows and circular shifts
def window_mean(values):
    count = np.isfinite(values).sum(axis=-1)
    return np.divide(np.nansum(values, axis=-1), count,
                     out=np.full(count.shape, np.nan), where=count > 0)


def stimulation_response(traces, event_frames, frame_mask, response_start):
    aligned = align_trials(traces, event_frames)
    indices = np.asarray(event_frames)[:, None] + np.arange(-60, 120)
    inside = (indices >= 0) & (indices < traces.shape[1])
    masked = np.zeros(indices.shape, dtype=bool)
    masked[inside] = frame_mask[indices[inside]]
    aligned = np.where(masked[None, :, :], np.nan, aligned)
    keep = np.all(np.isfinite(aligned) | masked[None, :, :], axis=2)
    aligned[~keep] = np.nan
    baseline = aligned[:, :, 30:60]
    response = aligned[:, :, response_start:response_start + 30]
    amplitude = np.nanmean(window_mean(response) - window_mean(baseline), axis=1)
    effect_size = amplitude / pooled_std(baseline, response)

    finite = np.isfinite(aligned)
    doubled = np.concatenate((np.where(finite, aligned, 0.0),
                              np.where(finite, aligned, 0.0)), axis=2)
    doubled_count = np.concatenate((finite, finite), axis=2)
    sums = np.concatenate((np.zeros((*keep.shape, 1)),
                           np.cumsum(doubled, axis=2)), axis=2)
    counts = np.concatenate((np.zeros((*keep.shape, 1), dtype=int),
                             np.cumsum(doubled_count, axis=2)), axis=2)
    rng = np.random.RandomState(SEED)
    shuffled = np.full((len(traces), N_SHUFFLE), np.nan)
    for start in range(0, N_SHUFFLE, 64):
        stop = min(start + 64, N_SHUFFLE)
        shifts = rng.randint(1, 180, size=(stop - start, len(event_frames)))
        means = []
        for first_frame in (30, response_start):
            first = (first_frame - shifts) % 180
            last = first + 30
            total = (np.take_along_axis(sums[:, None], last[None, :, :, None], axis=3)
                     - np.take_along_axis(sums[:, None], first[None, :, :, None], axis=3))[..., 0]
            count = (np.take_along_axis(counts[:, None], last[None, :, :, None], axis=3)
                     - np.take_along_axis(counts[:, None], first[None, :, :, None], axis=3))[..., 0]
            means.append(np.divide(total, count, out=np.full(total.shape, np.nan),
                                    where=count > 0))
        shuffled[:, start:stop] = np.nanmean(means[1] - means[0], axis=2)
    up, _ = classify_up(amplitude, effect_size, shuffled, np.zeros(len(traces), dtype=bool))
    return up, keep


#%% equal-count session summaries
def session_summary(row):
    recname = row.session
    folder = SESSION_ROOT / recname
    labels = pd.read_parquet(folder / f'{recname}_per_acr_trial.parquet')
    coords = labels[['grid_y', 'grid_x']].drop_duplicates().to_numpy(dtype=int)
    grids = coords[:, 0] * comparison.GRID_N + coords[:, 1]
    ctrl = labels.loc[labels.condition.eq('ctrl'),
                      ['trial_index', 'onset_frame']].drop_duplicates()
    ctrl_frames = ctrl.onset_frame.to_numpy(dtype=int)
    events = comparison.load_events(recname)
    stim_frames = events.onset_frame.to_numpy(dtype=int)
    if len(ctrl_frames) < len(stim_frames):
        raise ValueError('fewer ctrl trials than stimulation events')
    response_start = 60 + int(round(row.max_stim_duration_s * comparison.SAMP_FREQ)) + 1
    if response_start + 30 > 180:
        raise ValueError('stimulation response extends beyond the alignment window')
    rng = np.random.RandomState(SEED)
    draws = [rng.choice(ctrl_frames, len(stim_frames), replace=False) for _ in range(N_DRAWS)]
    result = dict(rec_id=recname, animal=row.animal, n_acr=len(grids),
                  n_ctrl_available=len(ctrl_frames), n_stim_events=len(stim_frames),
                  n_ctrl_per_draw=len(stim_frames), n_draws=N_DRAWS,
                  stim_response_start_s=(response_start - 60) / comparison.SAMP_FREQ,
                  status='ok')
    for channel, filename in (('dlight', 'dff_corrected_dlight.npy'), ('red', 'dff_red.npy')):
        traces = np.load(folder / 'single_trial_regression' / filename, mmap_mode='r')
        traces = np.asarray(traces).reshape(1024, -1)[grids]
        frame_mask = comparison.build_stim_mask(events, traces.shape[1])
        stim_up, stim_keep = stimulation_response(traces, stim_frames, frame_mask, response_start)
        ctrl_percent, ctrl_kept = [], []
        for frames in draws:
            tested = quantify_grid_response(traces, frames, n_shuffle=N_SHUFFLE, seed=SEED)
            up, _ = classify_up(tested['amplitude'], tested['effect_size'],
                                 tested['shuffle_amplitude'], np.zeros(len(grids), dtype=bool))
            ctrl_percent.append(100 * up.mean())
            ctrl_kept.append(tested['keep'].sum())
        result[f'{channel}_ctrl_up_pct'] = float(np.mean(ctrl_percent))
        result[f'{channel}_stim_up_pct'] = float(100 * stim_up.mean())
        result[f'{channel}_stim_minus_ctrl_pp'] = result[f'{channel}_stim_up_pct'] - result[f'{channel}_ctrl_up_pct']
        result[f'{channel}_ctrl_mean_kept_pairs'] = float(np.mean(ctrl_kept))
        result[f'{channel}_stim_kept_pairs'] = int(stim_keep.sum())
    return result


#%% run and report
def main():
    output = OUTPUT_ROOT / 'up_fraction_ctrl_stim.csv'
    if output.exists():
        raise FileExistsError(output)
    rows, failures = [], []
    for row in comparison.final_sessions().itertuples():
        try:
            result = session_summary(row)
            print(f'{row.session}: dLight ctrl/stim={result["dlight_ctrl_up_pct"]:.6g}/{result["dlight_stim_up_pct"]:.6g}%; '
                  f'red={result["red_ctrl_up_pct"]:.6g}/{result["red_stim_up_pct"]:.6g}%', flush=True)
        except (OSError, ValueError, KeyError) as error:
            failures.append(f'{row.session}: {type(error).__name__}: {error}')
            result = dict(rec_id=row.session, animal=row.animal, status=failures[-1])
            print(failures[-1], flush=True)
        rows.append(result)
    table = pd.DataFrame(rows)
    table.to_csv(output, index=False)
    lines = ['', '## 18 November 2026: Up ACR fractions in ctrl and stim', '',
             'Stimulation adaptation: detected-onset [-2, +4) s windows retain NaNs inside the inclusive stimulation mask and reject an ACR-event only for non-finite frames outside it; baseline is frames 30:60, response starts at 60 + round(max_stim_duration_s * 30) + 1 (the first frame after the maximum inclusive mask), and both observed and circular-shift windows use NaN-ignoring means and the existing NaN-ignoring pooled SD and Up rule.', '',
             'Ctrl uses quantify_grid_response unchanged on saved method-0 run-onset trials; the saved dLight ACR coordinates and already-prepared dff_corrected_dlight/dff_red arrays are used for both channels. Each session resets RandomState(42), draws the stimulation-event count of ctrl trials without replacement 20 times, and reuses the same draws for both channels; each response test uses 1,000 shifts and seed 42. Stim uses every detected event once. Percentages use the full saved valid ACR set, including ACRs with no retained trial as non-Up. Pooled summaries give each session equal weight; IQRs describe session dispersion.', '']
    for channel in ('dlight', 'red'):
        ctrl, stim = f'{channel}_ctrl_up_pct', f'{channel}_stim_up_pct'
        paired = table.loc[np.isfinite(table[ctrl]) & np.isfinite(table[stim])]
        difference = paired[stim].to_numpy() - paired[ctrl].to_numpy()
        if len(difference) == 0:
            statistic, pvalue = np.nan, np.nan
        elif np.all(difference == 0):
            statistic, pvalue = 0.0, 1.0
        else:
            test = wilcoxon(difference, alternative='two-sided')
            statistic, pvalue = float(test.statistic), float(test.pvalue)
        for condition, column in (('ctrl', ctrl), ('stim', stim)):
            values = paired[column].to_numpy()
            low, median, high = np.percentile(values, [25, 50, 75]) if len(values) else [np.nan] * 3
            lines.append(f'{channel}, {condition}: median Up={median:.9g}%, IQR [{low:.9g}, {high:.9g}]%; n={len(paired)} sessions, n={paired.animal.nunique()} animals.')
        median = float(np.median(difference)) if len(difference) else np.nan
        lines.append(f'{channel}, stim against ctrl: median paired difference={median:.9g} percentage points; two-sided paired Wilcoxon W={statistic:.9g}, p={pvalue:.12g}; n={len(paired)} sessions, n={paired.animal.nunique()} animals.')
    lines += ['', 'Missing inputs or exclusions: ' + ('; '.join(failures) if failures else 'none.'),
              '', 'Per-session results: outputs/da_hotspots/up_fraction_ctrl_stim.csv.', '']
    with (OUTPUT_ROOT / 'PROGRESS.md').open('a', encoding='utf-8') as stream:
        stream.write(chr(10).join(lines))
    print(chr(10).join(lines[-9:]), flush=True)


if __name__ == '__main__':
    with threadpool_limits(limits=1):
        main()
