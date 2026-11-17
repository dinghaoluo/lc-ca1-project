# -*- coding: utf-8 -*-
'''
Created on Mon Aug 31 18:46:58 2026

extract complete representative whole-field traces from the first observed

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
from scipy.ndimage import percentile_filter

repo_root = Path(__file__).resolve().parents[3]
for module_path in (repo_root, repo_root / 'utils'):
    if str(module_path) not in sys.path:
        sys.path.insert(0, str(module_path))

from behaviour_functions import find_nearest
import imaging_pipeline_functions as ipf


#%% paths and parameters
output_path = (
    repo_root / 'data' / 'catecholamine_comparison' / 'representative_inputs'
    / 'representative_post_shutter_traces.npz'
)
recordings = {
    'dLight': pp.TWO_PHOTON_ROOT / 'A152i/A152i-20251020/A152i-20251020-02',
    'nLight': pp.TWO_PHOTON_ROOT / 'A171i/A171i-20260220/A171i-20260220-01',
}
log_paths = {
    'dLight': pp.MICEEXP_ANMD152_A152_20251020_02T_TXT_PATH,
    'nLight': pp.MICEEXP_ANMD171_A171_20260220_01T_TXT_PATH,
}
expected_trials = {'dLight': 32, 'nLight': 18}
samp_freq = 30


#%% whole-field trace extraction
def whole_field_trace(root, channel):
    plane = root / 'suite2p' / 'plane0'
    ops = np.load(plane / 'ops.npy', allow_pickle=True).item()
    shape = (ops['nframes'], ops['Ly'], ops['Lx'])
    filename = 'data.bin' if channel == 1 else 'data_chan2.bin'
    movie = np.memmap(plane / filename, mode='r', dtype='int16', shape=shape)
    return np.sum(movie, axis=(1, 2), dtype=np.float64)


def dlight_shutter_edges(red_dff, pulse_trains, frame_times, pulse_info):
    onsets_raw, offsets_raw = ipf.detect_step_pairs(
        red_dff, zthr=10,
        min_interval_frames=pulse_info['min_interval_frames'],
    )
    onsets = np.asarray(onsets_raw, dtype=int) - 3
    offsets = np.asarray(offsets_raw, dtype=int) + 3
    valid = (offsets - onsets) < pulse_info['min_interval_frames']
    onsets, offsets = onsets[valid], offsets[valid]

    if abs(len(onsets) - len(pulse_trains)) > 10:
        onsets_raw, offsets_raw = ipf.detect_step_pairs(
            red_dff, zthr=20,
            min_interval_frames=pulse_info['min_interval_frames'],
        )
        onsets = np.asarray(onsets_raw, dtype=int) - 3
        offsets = np.asarray(offsets_raw, dtype=int) + 3
        valid = (offsets - onsets) < pulse_info['min_interval_frames']
        onsets, offsets = onsets[valid], offsets[valid]

    if abs(len(onsets) - len(pulse_trains)) > 10:
        pulse_frames = [
            [find_nearest(pulse, frame_times) for pulse in train]
            for train in pulse_trains
        ]
        onsets = np.asarray([frames[0] - 3 for frames in pulse_frames])
        offsets = np.asarray([frames[-1] + 3 for frames in pulse_frames])
        valid = (offsets - onsets) < pulse_info['min_interval_frames']
        onsets, offsets = onsets[valid], offsets[valid]

    return onsets, offsets


def local_shutter_onset(trace, pulse_start, pulse_end):
    base_start = max(0, pulse_start - samp_freq)
    base_end = max(base_start + 1, pulse_start - 5)
    search_start = max(0, pulse_start - int(.5 * samp_freq))
    search_end = min(len(trace), pulse_end + int(.25 * samp_freq))

    baseline = np.nanmedian(trace[base_start:base_end])
    search_trace = trace[search_start:search_end]
    low_level = np.nanpercentile(search_trace, 5)
    high_level = np.nanpercentile(search_trace, 95)
    closure_sign = -1 if abs(baseline - low_level) >= abs(high_level - baseline) else 1
    d = np.diff(trace)

    for zthr in [10, 5]:
        local_diff = np.diff(trace[search_start:search_end])
        median_diff = np.nanmedian(local_diff)
        mad = np.nanmedian(np.abs(local_diff - median_diff))
        candidates = (
            np.where(np.abs((local_diff - median_diff) / mad) > zthr)[0]
            + search_start + 1
        )
        candidates = candidates[(candidates > 0) & (candidates <= len(d))]
        if len(candidates) == 0:
            continue

        matched = candidates[closure_sign * d[candidates - 1] > 0]
        if len(matched):
            return int(matched[np.argmax(closure_sign * d[matched - 1])]) - 1
        return int(candidates[np.argmin(abs(candidates - pulse_start))]) - 1

    raise ValueError('no local shutter onset')


def local_shutter_recovery(trace, onset, pulse_end):
    stable_frames = 3
    buffer_frames = 2
    base_start = max(0, onset - samp_freq)
    base_end = max(base_start + 1, onset - 5)
    search_start = max(onset + 1, pulse_end)
    search_end = min(len(trace), pulse_end + int(2.5 * samp_freq))

    baseline = np.nanmedian(trace[base_start:base_end])
    closure_trace = trace[onset:min(search_end, onset + int(1.5 * samp_freq))]
    low_level = np.nanpercentile(closure_trace, 5)
    high_level = np.nanpercentile(closure_trace, 95)
    closure_low = abs(baseline - low_level) >= abs(high_level - baseline)
    shutter_level = low_level if closure_low else high_level
    threshold = shutter_level + .8 * (baseline - shutter_level)

    for frame in range(search_start, search_end - stable_frames + 1):
        segment = trace[frame:frame + stable_frames]
        if closure_low and np.all(segment >= threshold):
            return frame + buffer_frames
        if not closure_low and np.all(segment <= threshold):
            return frame + buffer_frames

    raise ValueError('no local shutter recovery')


def nlight_shutter_edges(red_raw, pulse_trains, frame_times, pulse_info):
    onsets, offsets = [], []
    for train in pulse_trains:
        pulse_start = int(find_nearest(train[0], frame_times))
        pulse_end = min(
            len(red_raw) - 1,
            pulse_start + max(1, pulse_info['total_train_duration_frames']),
        )
        onset = local_shutter_onset(red_raw, pulse_start, pulse_end)
        offset = local_shutter_recovery(red_raw, onset, pulse_end)
        onsets.append(onset)
        offsets.append(offset)
    return np.asarray(onsets), np.asarray(offsets)


#%% run
def main():
    traces = {}
    for sensor, root in recordings.items():
        print(f'{sensor}: extracting whole-field traces', flush=True)
        green_raw = whole_field_trace(root, channel=1)
        red_raw = whole_field_trace(root, channel=2)

        print(f'{sensor}: calculating 15th-percentile baseline', flush=True)
        green_baseline = percentile_filter(
            green_raw, percentile=15, size=9000, mode='reflect')
        green_dff = (green_raw - green_baseline) / green_baseline

        log = ipf.process_txt_nobeh(log_paths[sensor])
        frame_times = log['frame_times']
        pulse_info = ipf.extract_opto_pulse_metadata(log, samp_freq)
        pulse_trains = pulse_info['pulse_trains']

        if sensor == 'dLight':
            red_baseline = percentile_filter(
                red_raw, percentile=15, size=9000, mode='reflect')
            red_dff = (red_raw - red_baseline) / red_baseline
            onsets, offsets = dlight_shutter_edges(
                red_dff, pulse_trains, frame_times, pulse_info)
            first_observed = offsets + 1
        else:
            onsets, offsets = nlight_shutter_edges(
                red_raw, pulse_trains, frame_times, pulse_info)
            first_observed = offsets

        valid = [
            (onset, first)
            for onset, first in zip(onsets, first_observed)
            if onset - samp_freq >= 0 and first + 301 <= len(green_dff)
        ]
        if len(valid) != expected_trials[sensor]:
            raise AssertionError(
                f'{sensor}: expected {expected_trials[sensor]} trials, '
                f'recovered {len(valid)}'
            )

        trial_matrix = np.empty((len(valid), 301), dtype=np.float32)
        valid_first = np.asarray([first for onset, first in valid], dtype=int)
        valid_onsets = np.asarray([onset for onset, first in valid], dtype=int)
        for row, (onset, first) in enumerate(valid):
            baseline = np.mean(green_dff[onset - samp_freq:onset])
            trial_matrix[row] = green_dff[first:first + 301] - baseline
            next_onset = valid_onsets[row + 1] if row + 1 < len(valid_onsets) else len(green_dff)
            overlap = max(0, min(301, next_onset - first))
            trial_matrix[row, overlap:] = np.nan

        traces[sensor] = {
            'session': root.name,
            'samp_freq': samp_freq,
            'time_s': np.arange(301, dtype=np.float32) / samp_freq,
            'trial_dff': trial_matrix,
            'shutter_onset_frames': np.asarray(
                [onset for onset, first in valid], dtype=int),
            'first_post_shutter_frames': np.asarray(
                [first for onset, first in valid], dtype=int),
        }

    np.savez_compressed(output_path, traces=traces)
    print(f'saved {output_path}')


if __name__ == '__main__':
    main()
