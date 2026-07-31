# -*- coding: utf-8 -*-
'''
Created on Mon Nov  4 14:26:31 2024
Modified drastically in Feb 2025:
    - added GPU support and mostly fixed memory-leak problems

process and align ROI traces
this script is specific to axon/dendrite GCaMP processing because activity
    merging produces the structures below

overview of the merged data structure:
    - stat.npy contains (per usual) info about the ROIs; after merging and
        appending the merges, each ROI has 'imerge' and 'inmerge':
        - 'imerge' tells one the constituent ROI(s) of a merged ROI and
        - 'inmerge' tells one of which ROI(s) the current ROI is a part
    - therefore, if one ROI conforms to (1) iscell==True and (2) inmerge<1,
        this ROI is considered valid (sorted)

@author: Dinghao Luo
'''

#%% imports
import argparse
import gc
from pathlib import Path

import numpy as np
import pickle
import sys
from tqdm import tqdm
from time import time
from datetime import timedelta
from scipy.stats import sem

repo_root = Path(__file__).resolve().parents[3]
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))
behaviour_root = repo_root / 'preprocessing' / 'behaviour'
if str(behaviour_root) not in sys.path:
    sys.path.insert(0, str(behaviour_root))
import imaging_pipeline_functions as ipf
import lchpc_axon_support as support
import project_paths as pp

from common_functions import (
    get_GPU_availability,
    mpl_formatting,
    )
from console_formatting import print_session
mpl_formatting()

import rec_list
paths = rec_list.pathLCHPCGCaMPImmobile

BEHAVIOUR_FRAME_KEYS = ('start_cue_frames', 'reward_frames')


#%% parameters
FIGURE_SESSION_STEM = pp.LCHPC_AXON_IMMOBILE_FIGURES_STEM / 'all_sessions'

SAMP_FREQ = 30
BEF = 3
AFT = 7  # in seconds
XAXIS = np.arange((BEF + AFT) * SAMP_FREQ) / SAMP_FREQ - BEF

# correction coefficient for neuropil subtraction
CORR_COEF = 0.7


#%% GPU acceleration
cp, GPU_AVAILABLE, device_name = get_GPU_availability()
print(f'GPU_AVAILABLE={GPU_AVAILABLE}')


#%% processing function
def process_session(path, make_plots=False):
    recname = path[-17:]
    print_session(recname)

    plane_stem = pp.resolve_suite2p_session_stem(path) / 'plane0'
    ops_path = plane_stem / 'ops.npy'
    bin_path = plane_stem / 'data.bin'
    bin2_path = plane_stem / 'data_chan2.bin'
    stat_path = plane_stem / 'stat.npy'
    F_path = plane_stem / 'F.npy'
    F2_path = plane_stem / 'F_chan2.npy'
    Fneu_path = plane_stem / 'Fneu.npy'
    F2neu_path = plane_stem / 'Fneu_chan2.npy'

    # folders to put processed data and optional single-session plots
    proc_path = pp.LCHPC_AXON_IMMOBILE_STEM / 'all_sessions' / recname
    proc_data_path = proc_path / 'processed_data'
    figure_path = FIGURE_SESSION_STEM / recname

    proc_path.mkdir(parents=True, exist_ok=True)
    proc_data_path.mkdir(parents=True, exist_ok=True)

    beh_path = pp.behaviour_experiment_stem('LCHPCGCaMPImmobile') / f'{recname}.pkl'
    with open(beh_path, 'rb') as f:
        beh = pickle.load(f)
    missing_keys = [key for key in BEHAVIOUR_FRAME_KEYS if key not in beh]
    if missing_keys:
        raise KeyError(
            f'immobile behaviour data for {recname} is missing {missing_keys}; '
            'run preprocessing/behaviour/process_behaviour.py with '
            '--dataset LCHPCGCaMPImmobile'
        )
    if make_plots:
        figure_path.mkdir(parents=True, exist_ok=True)

    # load files
    stat = np.load(stat_path, allow_pickle=True)
    if 'inmerge' not in stat[0]:
        print('no merging detected; skipped')
        return

    ops = np.load(ops_path, allow_pickle=True).item()

    tot_frames = ops['nframes']

    # behaviour file
    start_cue_frames = beh['start_cue_frames']
    pump_frames = beh['reward_frames']

    # filtering
    start_cue_frames = np.asarray(start_cue_frames)
    start_cue_frames = start_cue_frames[start_cue_frames != -1]  # filter out the no-clear-onset trials

    pump_frames = np.asarray(pump_frames)
    pump_frames = pump_frames[pump_frames != -1]  # filter out the no-rew trials

    # reference images
    ref_path = proc_data_path / 'ref_mat_ch1.npy'
    ref_ch2_path = proc_data_path / 'ref_mat_ch2.npy'

    if not ref_path.exists() or not ref_ch2_path.exists():
        shape = (tot_frames, ops['Ly'], ops['Lx'])

        mov = np.memmap(bin_path, mode='r', dtype='int16', shape=shape)
        ref_im = ipf.plot_reference(
            mov,
            recname=recname,
            channel=1,
            outpath=proc_path,
            GPU_AVAILABLE=GPU_AVAILABLE,
            save_figure=make_plots,
            figure_outpath=figure_path
        )
        mov._mmap.close()
        del mov
        gc.collect()

        mov2 = np.memmap(bin2_path, mode='r', dtype='int16', shape=shape)
        ref_ch2_im = ipf.plot_reference(
            mov2,
            recname=recname,
            channel=2,
            outpath=proc_path,
            GPU_AVAILABLE=GPU_AVAILABLE,
            save_figure=make_plots,
            figure_outpath=figure_path
        )
        mov2._mmap.close()
        del mov2
        gc.collect()

    else:
        ref_im = np.load(ref_path, allow_pickle=True)
        ref_ch2_im = np.load(ref_ch2_path, allow_pickle=True)

        if make_plots:
            ipf.save_reference_figure(
                ref_im, figure_path, recname, channel=1
            )
            ipf.save_reference_figure(
                ref_ch2_im, figure_path, recname, channel=2
            )

    # 19 Mar 2025: we also want to process the constituent ROIs
    (
        valid_rois_dict,
        valid_rois,
        _,
        all_rois,
        constituent_rois_coord_dict,
        valid_rois_coord_dict,
        ) = support.prepare_merged_roi_metadata(
            ref_im, ref_ch2_im, stat, recname, proc_path,
            plot=make_plots,
            figure_path=figure_path
        )

    # -----------------
    # dF/F calculation
    # -----------------
    # match the active-running axon extractor: neuropil subtraction followed
    # by a rolling percentile baseline, 4 Feb 2026.
    t0 = time()

    # load traces
    F     = np.load(F_path, allow_pickle=True)
    F2    = np.load(F2_path, allow_pickle=True)
    Fneu  = np.load(Fneu_path, allow_pickle=True)
    F2neu = np.load(F2neu_path, allow_pickle=True)

    # correction using neuropil
    F_corr  = F - Fneu * CORR_COEF
    F2_corr = F2 - F2neu * CORR_COEF

    # calculate dF/F
    F_dFF = ipf.calculate_dFF_percentile(F_corr,
                                         t_axis=1,
                                         window_size=9000,
                                         GPU_AVAILABLE=GPU_AVAILABLE,
                                         device_name=device_name,
                                         progress_desc='Ch1 dF/F baseline')
    F2_dFF = ipf.calculate_dFF_percentile(F2_corr,
                                          window_size=9000,
                                          t_axis=1,
                                          GPU_AVAILABLE=GPU_AVAILABLE,
                                          progress_desc='Ch2 dF/F baseline')

    # saving to disk
    np.save(rf'{proc_data_path}\F_dFF.npy', F_dFF)
    np.save(rf'{proc_data_path}\F2_dFF.npy', F2_dFF)

    print(f'dF/F computed and saved ({timedelta(seconds=int(time()-t0))})')
    # ----------------------
    # dF/F calculation ends
    # ----------------------

    # cue-aligned extraction
    if make_plots:
        cue_path = figure_path / 'cue_aligned_single_roi'
        cue_path.mkdir(parents=True, exist_ok=True)

    # valid ROIs RO dict
    all_cue_dict = {}
    all_mean_cue_dict = {}
    all_cue_ch2_dict = {}
    all_mean_cue_ch2_dict = {}

    # constituent ROIs RO dict
    all_cue_const_dict = {}
    all_mean_cue_const_dict = {}
    all_cue_const_ch2_dict = {}
    all_mean_cue_const_ch2_dict = {}

    for roi in tqdm(all_rois, desc='aligning to cue-'):
        ca = F_dFF[roi]
        ref_ca = F2_dFF[roi]

        cue_aligned, cue_aligned_im, cue_aligned_ch2, _, tot_cue = (
            support.align_dual_channel_traces(
                ca,
                ref_ca,
                start_cue_frames,
                BEF,
                AFT,
                SAMP_FREQ,
                tot_frames,
                )
            )

        if roi in valid_rois:  # if ROI resulted from a final merge
            all_cue_dict[f'ROI {roi}'] = cue_aligned
            all_cue_ch2_dict[f'ROI {roi}'] = cue_aligned_ch2

            cue_aligned_mean = np.mean(cue_aligned, axis=0)
            cue_aligned_sem = sem(cue_aligned, axis=0)
            all_mean_cue_dict[f'ROI {roi}'] = cue_aligned_mean

            cue_aligned_mean_ch2 = np.mean(cue_aligned_ch2, axis=0)
            cue_aligned_sem_ch2 = sem(cue_aligned_ch2, axis=0)
            all_mean_cue_ch2_dict[f'ROI {roi}'] = cue_aligned_mean_ch2

            if make_plots:
                support.save_roi_alignment_plots(
                    ref_im=ref_im,
                    stat=stat,
                    roi=roi,
                    aligned_im=cue_aligned_im,
                    aligned_mean=cue_aligned_mean,
                    aligned_sem=cue_aligned_sem,
                    aligned_mean_ch2=cue_aligned_mean_ch2,
                    aligned_sem_ch2=cue_aligned_sem_ch2,
                    xaxis=XAXIS,
                    before=BEF,
                    after=AFT,
                    total_events=tot_cue,
                    out_dir=cue_path,
                    event_title='cue-aligned',
                    xlabel='time from cue (s)',
                    dual_axis_ch2=False,
                    )

        else:  # if ROI is constituent to a valid ROI
            all_cue_const_dict[f'ROI {roi}'] = cue_aligned
            all_cue_const_ch2_dict[f'ROI {roi}'] = cue_aligned_ch2

    # save valid ROI data
    valid_rois_dict = {f'ROI {roi}': valid_rois_dict[roi]
                       for roi in valid_rois_dict}  # rename keys to align with other dicts
    np.save(rf'{proc_data_path}\valid_ROIs_dict.npy', valid_rois_dict)
    np.save(rf'{proc_data_path}\valid_ROIs_coord_dict.npy', valid_rois_coord_dict)
    np.save(rf'{proc_data_path}\constituent_ROIs_coord_dict.npy', constituent_rois_coord_dict)

    # save merged ROIs
    np.save(rf'{proc_data_path}\RO_aligned_dict.npy', all_cue_dict)
    np.save(rf'{proc_data_path}\RO_aligned_mean_dict.npy', all_mean_cue_dict)
    np.save(rf'{proc_data_path}\RO_aligned_ch2_dict.npy', all_cue_ch2_dict)
    np.save(rf'{proc_data_path}\RO_aligned_mean_ch2_dict.npy', all_mean_cue_ch2_dict)

    # save constituent ROIs
    np.save(rf'{proc_data_path}\RO_aligned_const_dict.npy', all_cue_const_dict)
    np.save(rf'{proc_data_path}\RO_aligned_const_mean_dict.npy', all_mean_cue_const_dict)
    np.save(rf'{proc_data_path}\RO_aligned_const_ch2_dict.npy', all_cue_const_ch2_dict)
    np.save(rf'{proc_data_path}\RO_aligned_const_mean_ch2_dict.npy', all_mean_cue_const_ch2_dict)

    # reward-aligned extraction
    if make_plots:
        rew_path = figure_path / 'rew_aligned_single_roi'
        rew_path.mkdir(parents=True, exist_ok=True)

    # valid ROIs RO dict
    all_rew_dict = {}
    all_mean_rew_dict = {}
    all_rew_ch2_dict = {}
    all_mean_rew_ch2_dict = {}

    # constituent ROIs RO dict
    all_rew_const_dict = {}
    all_mean_rew_const_dict = {}
    all_rew_const_ch2_dict = {}
    all_mean_rew_const_ch2_dict = {}

    for roi in tqdm(all_rois, desc='aligning to rewards'):
        ca = F_dFF[roi]
        ref_ca = F2_dFF[roi]

        rew_aligned, rew_aligned_im, rew_aligned_ch2, _, tot_pump = (
            support.align_dual_channel_traces(
                ca,
                ref_ca,
                pump_frames,
                BEF,
                AFT,
                SAMP_FREQ,
                tot_frames,
                )
            )

        if roi in valid_rois:  # if ROI resulted from a final merge
            all_rew_dict[f'ROI {roi}'] = rew_aligned
            all_rew_ch2_dict[f'ROI {roi}'] = rew_aligned_ch2

            rew_aligned_mean = np.mean(rew_aligned, axis=0)
            rew_aligned_sem = sem(rew_aligned, axis=0)
            all_mean_rew_dict[f'ROI {roi}'] = rew_aligned_mean

            rew_aligned_mean_ch2 = np.mean(rew_aligned_ch2, axis=0)
            rew_aligned_sem_ch2 = sem(rew_aligned_ch2, axis=0)
            all_mean_rew_ch2_dict[f'ROI {roi}'] = rew_aligned_mean_ch2

            if make_plots:
                support.save_roi_alignment_plots(
                    ref_im=ref_im,
                    stat=stat,
                    roi=roi,
                    aligned_im=rew_aligned_im,
                    aligned_mean=rew_aligned_mean,
                    aligned_sem=rew_aligned_sem,
                    aligned_mean_ch2=rew_aligned_mean_ch2,
                    aligned_sem_ch2=rew_aligned_sem_ch2,
                    xaxis=XAXIS,
                    before=BEF,
                    after=AFT,
                    total_events=tot_pump,
                    out_dir=rew_path,
                    event_title='rew-aligned',
                    xlabel='time from rew (s)',
                    dual_axis_ch2=True,
                    )

        else:  # if ROI is constituent to a valid ROI
            all_rew_const_dict[f'ROI {roi}'] = rew_aligned
            all_rew_const_ch2_dict[f'ROI {roi}'] = rew_aligned_ch2

    # save merged ROIs
    np.save(rf'{proc_data_path}\rew_aligned_dict.npy', all_rew_dict)
    np.save(rf'{proc_data_path}\rew_aligned_mean_dict.npy', all_mean_rew_dict)
    np.save(rf'{proc_data_path}\rew_aligned_ch2_dict.npy', all_rew_ch2_dict)
    np.save(rf'{proc_data_path}\rew_aligned_mean_ch2_dict.npy', all_mean_rew_ch2_dict)

    # save constituent ROIs
    np.save(rf'{proc_data_path}\rew_aligned_const_dict.npy', all_rew_const_dict)
    np.save(rf'{proc_data_path}\rew_aligned_const_mean_dict.npy', all_mean_rew_const_dict)
    np.save(rf'{proc_data_path}\rew_aligned_const_ch2_dict.npy', all_rew_const_ch2_dict)
    np.save(rf'{proc_data_path}\rew_aligned_const_mean_ch2_dict.npy', all_mean_rew_const_ch2_dict)

#%% main
def main(argv=None):
    parser = argparse.ArgumentParser(
        description='extract immobile LC-HPC axon-GCaMP aligned ROI data.'
    )
    parser.add_argument(
        '--recording',
        dest='recording_filter',
        help='optional substring filter for recname or recording path',
    )
    parser.add_argument(
        '--plots',
        dest='make_plots',
        action='store_true',
        help='also regenerate per-session figures during extraction',
    )
    parser.add_argument(
        '--no-plots',
        dest='make_plots',
        action='store_false',
        help=argparse.SUPPRESS,
    )
    parser.set_defaults(make_plots=False)
    args = parser.parse_args(argv)

    session_paths = paths
    if args.recording_filter:
        session_paths = [
            rec_path for rec_path in paths
            if args.recording_filter in rec_path
            or args.recording_filter in rec_path[-17:]
        ]
        if not session_paths:
            raise ValueError(f'no recordings matched filter: {args.recording_filter}')

    for path in session_paths:
        process_session(
            path,
            make_plots=args.make_plots,
        )
        recname = Path(path).name
        print(f'clearing CPU/GPU memory after {recname}...', flush=True)
        gc.collect()
        if GPU_AVAILABLE:
            cp.cuda.Stream.null.synchronize()
            cp.get_default_memory_pool().free_all_blocks()
            cp.get_default_pinned_memory_pool().free_all_blocks()

if __name__ == '__main__':
    main()
