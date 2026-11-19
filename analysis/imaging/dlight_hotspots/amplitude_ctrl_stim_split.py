'''
Created on Wed Nov 18 22:42:59 2026

compare ctrl and stim amplitudes after splitting the ctrl trials

@author: Dinghao Luo
'''

#%% imports
import numpy as np
import pandas as pd
from scipy.stats import wilcoxon
from threadpoolctl import threadpool_limits

from amplitude_ctrl_stim import (stim_amplitude, comparison, window_mean, SESSION_ROOT,
                                 OUTPUT_ROOT, N_SHUFFLE, N_DRAWS, SEED)
from grid_response_test import classify_up, quantify_grid_response


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
    halves = np.array_split(rng.permutation(ctrl_frames), 2)
    draws = []
    short_halves = 0
    for half in halves:
        if len(half) < len(stim_frames):
            draws.append(np.arange(len(half))[None, :])
            short_halves += 1
        else:
            draws.append(np.array([rng.choice(len(half), len(stim_frames), replace=False)
                                   for _ in range(N_DRAWS)]))

    rows, disagreements = [], {}
    for channel, filename in (('dlight', 'dff_corrected_dlight.npy'), ('red', 'dff_red.npy')):
        traces = np.load(folder / 'single_trial_regression' / filename, mmap_mode='r')
        traces = traces.reshape(comparison.GRID_N ** 2, -1)[grids]
        up_by_half, amplitude_by_half = [], []
        for half, selected in zip(halves, draws):
            tested = quantify_grid_response(traces, half, n_shuffle=N_SHUFFLE, seed=SEED)
            up, _ = classify_up(tested['amplitude'], tested['effect_size'],
                                tested['shuffle_amplitude'], np.zeros(len(grids), dtype=bool))
            up_by_half.append(up)
            amplitude_by_half.append(window_mean(window_mean(tested['trial_amplitude'][:, selected])))
            del tested

        frame_mask = comparison.build_stim_mask(events, traces.shape[1])
        stim = stim_amplitude(traces, stim_frames, frame_mask, response_start)

        up = np.column_stack(up_by_half)
        # the group is called on the half opposite the amplitude in each direction
        held_ctrl = np.column_stack(amplitude_by_half[::-1])
        disagreements[channel] = int(np.count_nonzero(up[:, 0] != up[:, 1]))
        for group, membership in (('ctrl_up', up), ('other', ~up)):
            ctrl_amplitude = window_mean(np.where(membership, held_ctrl, np.nan))
            paired = np.isfinite(ctrl_amplitude) & np.isfinite(stim)
            ctrl_values, stim_values = ctrl_amplitude[paired], stim[paired]
            medians = (np.median(ctrl_values), np.median(stim_values),
                       np.median(stim_values - ctrl_values)) if np.any(paired) else (np.nan,) * 3
            rows.append(dict(rec_id=recname, animal=row.animal, channel=channel, group=group,
                             n_ctrl_trials=len(ctrl_frames), n_stim_events=len(stim_frames),
                             n_valid_acr=int(membership.any(axis=1).sum()),
                             n_paired_acr=int(paired.sum()),
                             median_ctrl=medians[0], median_stim=medians[1],
                             median_stim_minus_ctrl=medians[2]))
    return rows, short_halves, disagreements


#%% run and report
def main():
    output = OUTPUT_ROOT / 'amplitude_ctrl_stim_split.csv'
    if output.exists():
        raise FileExistsError(output)

    rows = []
    short_sessions = short_halves = 0
    disagreements = {'dlight': 0, 'red': 0}
    for row in comparison.final_sessions().itertuples():
        session_rows, short, discordant = session_summary(row)
        rows.extend(session_rows)
        short_sessions += int(short > 0)
        short_halves += short
        for channel in disagreements:
            disagreements[channel] += discordant[channel]
        print(f'{row.session}: {short} ctrl halves shorter than the stim event count', flush=True)
    table = pd.DataFrame(rows)
    lines = [
        '', '## 18 November 2026: ctrl and stim amplitudes after splitting ctrl trials', '',
        f'I split the saved ctrl trials at random (seed {SEED}) in {table.rec_id.nunique()} '
        f'sessions from {table.animal.nunique()} animals, then called Up separately in each '
        f'channel on one half and estimated ctrl amplitude on the other, with the roles '
        f'swapped for the second estimate. Each held-out half supplies {N_DRAWS} draws '
        f'without replacement, matched to the number of stim events where possible. I used '
        f'all trials in {short_halves} shorter halves across {short_sessions} sessions.',
        '',
        f'I averaged each ACR\'s finite held-out estimates within its group before taking '
        f'the session medians. When the Up call differed between halves, the ACR contributed '
        f'once to each group with the corresponding held-out estimate. This occurred for '
        f'{disagreements["dlight"]} dLight and {disagreements["red"]} red ACRs. The stim '
        f'amplitude uses the same events and mask as above and is the same in both split '
        f'directions.',
        '',
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
    lines += ['', 'Per-session results in `outputs/da_hotspots/amplitude_ctrl_stim_split.csv`.', '']

    # insert after the full-ctrl results and preserve the existing log bytes
    anchor = b'Per-session results in `outputs/da_hotspots/amplitude_ctrl_stim.csv`.'
    with (OUTPUT_ROOT / 'PROGRESS.md').open('r+b') as stream:
        existing = stream.read()
        insert_at = existing.index(b'\n', existing.index(anchor)) + 1
        newline = b'\r\n' if existing[insert_at - 2:insert_at] == b'\r\n' else b'\n'
        table.to_csv(output, index=False)
        stream.seek(insert_at)
        stream.write(newline.join(line.encode('utf-8') for line in lines) + existing[insert_at:])
    print('\n'.join(lines), flush=True)


if __name__ == '__main__':
    with threadpool_limits(limits=1):
        main()
