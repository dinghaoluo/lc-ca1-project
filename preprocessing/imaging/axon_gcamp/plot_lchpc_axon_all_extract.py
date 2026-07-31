# -*- coding: utf-8 -*-
'''
Created on 21 May 2026
Modified on 29 June 2026

plot saved LC-HPC axon-GCaMP ROI extractions without recalculating dF/F or
aligned traces

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

import matplotlib
matplotlib.use('Agg')

repo_root = Path(__file__).resolve().parents[3]
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))

import lchpc_axon_support as support
import project_paths as pp
import imaging_pipeline_functions as ipf
from common_functions import mpl_formatting, normalise, smooth_convolve
mpl_formatting()

import rec_list


#%% paths and parameters
paths = rec_list.pathLCHPCGCaMP[1:]

AXON_GCAMP_STEM = pp.LCHPC_AXON_STEM
ALL_SESSION_STEM = AXON_GCAMP_STEM / 'all_sessions'
AXON_GCAMP_FIGURE_STEM = pp.LCHPC_AXON_FIGURES_STEM
ALL_SESSION_FIGURE_STEM = AXON_GCAMP_FIGURE_STEM / 'all_sessions'
BEH_STEM = pp.behaviour_experiment_stem('LCHPCGCaMP')

SAMP_FREQ = 30
BEF = 3
AFT = 7
XAXIS = np.arange((BEF + AFT) * SAMP_FREQ) / SAMP_FREQ - BEF

EVENT_SPECS = {
    'run': {
        'frame_key': 'run_onset_frames',
        'drop_missing': True,
        'data_stem': 'all_run',
        'out_dir': 'ro_aligned_single_roi',
        'desc': 'run onset',
        'event_title': 'run-onset-aligned',
        'xlabel': 'Time from run-onset (s)',
        'overview_xlabel': 'Time from run onset (s)',
        'dual_axis_ch2': False,
        },
    'rew': {
        'frame_key': 'reward_frames',
        'drop_missing': True,
        'data_stem': 'all_rew',
        'out_dir': 'rew_aligned_single_roi',
        'desc': 'rewards',
        'event_title': 'rew-aligned',
        'xlabel': 'time from rew (s)',
        'overview_xlabel': None,
        'dual_axis_ch2': True,
        },
    'cue': {
        'frame_key': 'start_cue_frames',
        'drop_missing': False,
        'data_stem': 'all_cue',
        'out_dir': 'cue_aligned_single_roi',
        'desc': 'cues',
        'event_title': 'cue-aligned',
        'xlabel': 'time from cue (s)',
        'overview_xlabel': None,
        'dual_axis_ch2': True,
        },
    }


#%% plotting
def save_roi_metadata_plots(rec_path, proc_data_path, figure_path, recname):
    stat = np.load(pp.resolve_suite2p_session_stem(rec_path) / 'plane0' / 'stat.npy',
                   allow_pickle=True)
    ref_im = np.load(proc_data_path / 'ref_mat_ch1.npy', allow_pickle=True)
    ref_ch2_im = np.load(proc_data_path / 'ref_mat_ch2.npy', allow_pickle=True)
    ipf.save_reference_figure(
        ref_im, figure_path, recname, channel=1
        )
    ipf.save_reference_figure(
        ref_ch2_im, figure_path, recname, channel=2
        )

    valid_rois_dict = np.load(
        proc_data_path / 'valid_ROIs_dict.npy',
        allow_pickle=True
        ).item()
    valid_rois = set(valid_rois_dict)
    constituent_rois = {
        roi
        for roi_list in valid_rois_dict.values()
        for roi in roi_list
        }

    support.calculate_and_plot_overlap_indices(
        ref_im, ref_ch2_im,
        stat, valid_rois, recname, figure_path
        )
    support.get_roi_coord_dict(
        ref_im, ref_ch2_im,
        stat, constituent_rois, recname, figure_path
        )
    support.get_roi_coord_dict(
        ref_im, ref_ch2_im,
        stat, valid_rois, recname, figure_path
        )

    return stat, ref_im

def save_event_alignment_plots(
        recname, proc_path, proc_data_path, figure_path,
        stat, ref_im, beh, event_key):
    event_spec = EVENT_SPECS[event_key]
    data_stem = event_spec['data_stem']
    aligned = np.load(proc_path / f'{recname}_{data_stem}.npy', allow_pickle=True).item()
    aligned_mean = np.load(proc_path / f'{recname}_{data_stem}_mean.npy', allow_pickle=True).item()
    aligned_ch2 = np.load(proc_path / f'{recname}_{data_stem}_ch2.npy', allow_pickle=True).item()
    aligned_mean_ch2 = np.load(
        proc_path / f'{recname}_{data_stem}_ch2_mean.npy',
        allow_pickle=True
        ).item()

    event_frames = np.asarray(beh[event_spec['frame_key']])
    if event_spec['drop_missing']:
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
            overview_xlabel=event_spec['overview_xlabel'],
            dual_axis_ch2=event_spec['dual_axis_ch2'],
            )

def plot_session(rec_path, event_keys, include_roi_metadata=True):
    recname = Path(rec_path).name
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
            recname, proc_path, proc_data_path, figure_path,
            stat, ref_im, beh, event_key
            )


#%% main
def main(argv=None):
    parser = argparse.ArgumentParser(
        description='plot LC-HPC axon-GCaMP ROI figures from saved extraction outputs.'
        )
    parser.add_argument(
        '--recording',
        dest='recording_filter',
        help='optional substring filter for recname or recording path',
        )
    parser.add_argument(
        '--event',
        choices=['all', 'run', 'rew', 'cue'],
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
                or args.recording_filter in Path(rec_path).name)
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
