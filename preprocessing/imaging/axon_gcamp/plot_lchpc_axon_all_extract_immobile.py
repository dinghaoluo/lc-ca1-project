# -*- coding: utf-8 -*-
'''
Created on 22 May 2026
Modified on 29 June 2026

plot saved immobile LC-HPC axon-GCaMP ROI extractions

@author: Dinghao Luo
'''

#%% imports
import argparse
import pickle
import sys
from pathlib import Path

import numpy as np
from scipy.stats import sem
from tqdm import tqdm

repo_root = Path(__file__).resolve().parents[3]
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))

script_dir = Path(__file__).resolve().parent
if str(script_dir) not in sys.path:
    sys.path.insert(0, str(script_dir))

import lchpc_axon_support as support
import project_paths as pp
from common_functions import mpl_formatting, normalise, smooth_convolve
from plot_lchpc_axon_all_extract import save_roi_metadata_plots
mpl_formatting()

import rec_list


#%% paths and parameters
paths = rec_list.pathLCHPCGCaMPImmobile

AXON_GCAMP_IMMOBILE_STEM = pp.LCHPC_AXON_IMMOBILE_STEM
ALL_SESSION_STEM = AXON_GCAMP_IMMOBILE_STEM / 'all_sessions'
AXON_GCAMP_IMMOBILE_FIGURE_STEM = pp.LCHPC_AXON_IMMOBILE_FIGURES_STEM
ALL_SESSION_FIGURE_STEM = AXON_GCAMP_IMMOBILE_FIGURE_STEM / 'all_sessions'
BEH_STEM = pp.behaviour_experiment_stem('LCHPCGCaMPImmobile')

SAMP_FREQ = 30
BEF = 3
AFT = 7
XAXIS = np.arange((BEF + AFT) * SAMP_FREQ) / SAMP_FREQ - BEF

EVENT_SPECS = {
    'cue': {
        'frame_key': 'start_cue_frames',
        'data_stem': 'RO_aligned',
        'out_dir': 'cue_aligned_single_roi',
        'desc': 'cues',
        'event_title': 'cue-aligned',
        'xlabel': 'time from cue (s)',
        'dual_axis_ch2': False,
        },
    'rew': {
        'frame_key': 'reward_frames',
        'data_stem': 'rew_aligned',
        'out_dir': 'rew_aligned_single_roi',
        'desc': 'rewards',
        'event_title': 'rew-aligned',
        'xlabel': 'time from rew (s)',
        'dual_axis_ch2': True,
        },
    }


#%% plotting
def save_event_alignment_plots(
        recname, proc_data_path, figure_path,
        stat, ref_im, beh, event_key):
    event_spec = EVENT_SPECS[event_key]
    data_stem = event_spec['data_stem']

    aligned = np.load(proc_data_path / f'{data_stem}_dict.npy', allow_pickle=True).item()
    aligned_mean = np.load(proc_data_path / f'{data_stem}_mean_dict.npy', allow_pickle=True).item()
    aligned_ch2 = np.load(proc_data_path / f'{data_stem}_ch2_dict.npy', allow_pickle=True).item()
    aligned_mean_ch2 = np.load(
        proc_data_path / f'{data_stem}_mean_ch2_dict.npy',
        allow_pickle=True
        ).item()

    event_frames = np.asarray(beh[event_spec['frame_key']])
    event_frames = event_frames[event_frames != -1]
    total_events = 0
    last_frame = float('-inf')
    for frame in event_frames:
        if frame > last_frame:
            total_events += 1
            last_frame = frame
    out_dir = figure_path / event_spec['out_dir']
    out_dir.mkdir(parents=True, exist_ok=True)

    for roiname in tqdm(
            sorted(aligned, key=lambda name: int(name.split(' ')[1])),
                        desc=f"plotting {recname} {event_spec['desc']}"):
        roi = int(roiname.split(' ')[1])
        roi_aligned = aligned[roiname]
        roi_aligned_ch2 = aligned_ch2[roiname]

        support.save_roi_alignment_plots(
            ref_im=ref_im,
            stat=stat,
            roi=roi,
            aligned_im=np.asarray([
                normalise(smooth_convolve(trace))
                for trace in roi_aligned
                ]),
            aligned_mean=aligned_mean[roiname],
            aligned_sem=sem(roi_aligned, axis=0),
            aligned_mean_ch2=aligned_mean_ch2[roiname],
            aligned_sem_ch2=sem(roi_aligned_ch2, axis=0),
            xaxis=XAXIS,
            before=BEF,
            after=AFT,
            total_events=total_events,
            out_dir=out_dir,
            event_title=event_spec['event_title'],
            xlabel=event_spec['xlabel'],
            dual_axis_ch2=event_spec['dual_axis_ch2'],
            )

def plot_session(rec_path, event_keys, include_roi_metadata=True):
    recname = rec_path[-17:]
    print(f'\nplotting {recname}')

    proc_path = ALL_SESSION_STEM / recname
    proc_data_path = proc_path / 'processed_data'
    figure_path = ALL_SESSION_FIGURE_STEM / recname

    if include_roi_metadata:
        stat, ref_im = save_roi_metadata_plots(
            rec_path, proc_data_path, figure_path, recname
            )
    else:
        stat = np.load(
            pp.resolve_suite2p_session_stem(rec_path) / 'plane0' / 'stat.npy',
            allow_pickle=True
            )
        ref_im = np.load(proc_data_path / 'ref_mat_ch1.npy', allow_pickle=True)

    with open(BEH_STEM / f'{recname}.pkl', 'rb') as f:
        beh = pickle.load(f)

    for event_key in event_keys:
        save_event_alignment_plots(
            recname, proc_data_path, figure_path,
            stat, ref_im, beh, event_key
            )


#%% main
def main(argv=None):
    parser = argparse.ArgumentParser(
        description='plot immobile LC-HPC axon-GCaMP ROI figures from saved extraction outputs.'
        )
    parser.add_argument(
        '--recording',
        dest='recording_filter',
        help='optional substring filter for recname or recording path',
        )
    parser.add_argument(
        '--event',
        choices=['all', 'cue', 'rew'],
        default='all',
        help='which alignment figures to regenerate',
        )
    parser.add_argument(
        '--skip-roi-metadata',
        action='store_true',
        help='skip ROI/reference validation plots and regenerate alignment plots only',
        )
    args = parser.parse_args(argv)

    event_keys = list(EVENT_SPECS) if args.event == 'all' else [args.event]

    selected_paths = paths
    if args.recording_filter:
        selected_paths = [
            rec_path for rec_path in paths
            if (args.recording_filter in rec_path
                or args.recording_filter in rec_path[-17:])
            ]
        if not selected_paths:
            raise ValueError(f'no recordings matched filter: {args.recording_filter}')

    for rec_path in selected_paths:
        plot_session(
            rec_path,
            event_keys,
            include_roi_metadata=not args.skip_roi_metadata,
            )

if __name__ == '__main__':
    main()
