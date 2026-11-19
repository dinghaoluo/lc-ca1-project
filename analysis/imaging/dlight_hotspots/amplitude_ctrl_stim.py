'''
Created on Wed Nov 18 22:17:17 2026

compare matched-trial ctrl and stim amplitudes by ctrl Up status

@author: Dinghao Luo
'''

#%% imports
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon
from threadpoolctl import threadpool_limits

from up_fraction_ctrl_stim import (comparison, window_mean, SESSION_ROOT,
                                   OUTPUT_ROOT, N_SHUFFLE, N_DRAWS, SEED)
from grid_response_test import align_trials, classify_up, quantify_grid_response


#%% stimulation amplitude
def stim_amplitude(traces, frames, frame_mask, response_start):
    aligned = align_trials(traces, frames)
    indices = frames[:, None] + np.arange(-60, 120)
    inside = (indices >= 0) & (indices < traces.shape[1])
    masked = np.zeros(indices.shape, dtype=bool)
    masked[inside] = frame_mask[indices[inside]]
    aligned = np.where(masked[None, :, :], np.nan, aligned)
    keep = np.all(np.isfinite(aligned) | masked[None, :, :], axis=2)
    aligned[~keep] = np.nan
    baseline = window_mean(aligned[:, :, 30:60])
    response = window_mean(aligned[:, :, response_start:response_start + 30])
    return window_mean(response - baseline)


#%% session amplitudes
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
    response_start = 60 + int(round(row.max_stim_duration_s * comparison.SAMP_FREQ)) + 1
    if response_start + 30 > 180:
        raise ValueError(f'{recname}: stimulation response extends beyond the alignment window')
    rng = np.random.RandomState(SEED)
    draws = np.array([rng.choice(len(ctrl_frames), len(stim_frames), replace=False)
                      for _ in range(N_DRAWS)])

    for channel, filename in (('dlight', 'dff_corrected_dlight.npy'), ('red', 'dff_red.npy')):
        traces = np.load(folder / 'single_trial_regression' / filename, mmap_mode='r')
        traces = traces.reshape(comparison.GRID_N ** 2, -1)[grids]
        tested = quantify_grid_response(traces, ctrl_frames, n_shuffle=N_SHUFFLE, seed=SEED)
        up, _ = classify_up(tested['amplitude'], tested['effect_size'],
                            tested['shuffle_amplitude'], np.zeros(len(grids), dtype=bool))
        # the NaN rule acts on each aligned trial before drawing, so the same amplitudes can be reused
        ctrl_amplitude = window_mean(window_mean(tested['trial_amplitude'][:, draws]))
        del tested
        frame_mask = comparison.build_stim_mask(events, traces.shape[1])
        stim = stim_amplitude(traces, stim_frames, frame_mask, response_start)

        for group, subset in (('ctrl_up', up), ('other', ~up)):
            paired = subset & np.isfinite(ctrl_amplitude) & np.isfinite(stim)
            ctrl_values, stim_values = ctrl_amplitude[paired], stim[paired]
            medians = (np.median(ctrl_values), np.median(stim_values),
                       np.median(stim_values - ctrl_values)) if np.any(paired) else (np.nan,) * 3
            yield dict(rec_id=recname, animal=row.animal, channel=channel, group=group,
                       n_ctrl_trials=len(ctrl_frames), n_stim_events=len(stim_frames),
                       n_valid_acr=int(subset.sum()), n_paired_acr=int(paired.sum()),
                       median_ctrl=medians[0], median_stim=medians[1],
                       median_stim_minus_ctrl=medians[2])


#%% run and report
def main():
    output = OUTPUT_ROOT / 'amplitude_ctrl_stim.csv'
    if output.exists():
        raise FileExistsError(output)

    rows = []
    for row in comparison.final_sessions().itertuples():
        rows.extend(session_summary(row))
        print(f'{row.session}: amplitudes computed for both channels', flush=True)
    table = pd.DataFrame(rows)
    lines = [
        '', '## 18 November 2026: ctrl and stim amplitudes by ctrl Up status', '',
        'I used the same saved ACRs, trials, windows, inclusive stimulation mask and 20 '
        f'equal-count ctrl draws (seed 42) as in the Up fraction comparison, across '
        f'{table.rec_id.nunique()} sessions from {table.animal.nunique()} animals. I called '
        f'ctrl Up separately in each channel using all saved ctrl trials. Each ACR contributes '
        f'its mean response minus baseline in stim and its mean across the ctrl draws that '
        f'gave it a finite amplitude. The session medians use the same finite ACR pairs in '
        f'both conditions, and the paired median difference is taken across ACRs within '
        f'each session.', '',
        'Pooled summaries give each session equal weight. The Wilcoxon test pairs the '
        'session ctrl and stim medians, and the difference column summarises the within-session '
        'medians of the individual ACR differences.', '',
        '| Channel | Group | n sessions | n animals | Ctrl median [IQR] | Stim median [IQR] | Paired ACR difference median [IQR] | Wilcoxon W | Two-sided p |',
        '| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |',
    ]
    for (channel, group), subset in table.groupby(['channel', 'group'], sort=False):
        paired = subset.dropna(subset=['median_ctrl', 'median_stim'])
        stats = []
        for column in ('median_ctrl', 'median_stim', 'median_stim_minus_ctrl'):
            low, median, high = paired[column].quantile([0.25, 0.5, 0.75])
            stats.append(f'{median:.6g} [{low:.6g}, {high:.6g}]')
        difference = paired.median_stim.to_numpy() - paired.median_ctrl.to_numpy()
        if not len(difference):
            statistic, pvalue = np.nan, np.nan
        elif np.all(difference == 0):
            statistic, pvalue = 0.0, 1.0
        else:
            statistic, pvalue = wilcoxon(difference, alternative='two-sided')
        lines.append(f'| {channel} | {group} | {len(paired)} | {paired.animal.nunique()} | '
                     f'{" | ".join(stats)} | {statistic:.6g} | {pvalue:.9g} |')
    lines += ['', 'Per-session results in `outputs/da_hotspots/amplitude_ctrl_stim.csv`.',
              '']
    table.to_csv(output, index=False)
    with (OUTPUT_ROOT / 'PROGRESS.md').open('a', encoding='utf-8') as stream:
        stream.write('\n'.join(lines))
    print('\n'.join(lines[9:]), flush=True)


if __name__ == '__main__':
    with threadpool_limits(limits=1):
        main()
