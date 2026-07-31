# -*- coding: utf-8 -*-
'''
Created on Mon Nov  4 14:26:31 2024
Modified drastically in Feb 2025:
    - added GPU support and mostly fixed memory-leak problems
Modified on 19 Jan 2026

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
import sys

import numpy as np
import pickle
from tqdm import tqdm
from time import time
from datetime import timedelta
from scipy.stats import sem

repo_root = Path(__file__).resolve().parents[3]
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))

import imaging_pipeline_functions as ipf
import project_paths as pp
import lchpc_axon_support as support

from common_functions import mpl_formatting, get_GPU_availability
from console_formatting import print_session
mpl_formatting()

import rec_list
pathHPCLCGCaMP = rec_list.pathLCHPCGCaMP


#%% paths and parameters
AXON_GCAMP_STEM = pp.LCHPC_AXON_STEM
ALL_SESSION_STEM = AXON_GCAMP_STEM / 'all_sessions'
FIGURE_SESSION_STEM = pp.LCHPC_AXON_FIGURES_STEM / 'all_sessions'

BEH_STEM = pp.behaviour_experiment_stem('LCHPCGCaMP')

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
    recname = Path(path).name
    print_session(recname)

    plane_stem = pp.resolve_suite2p_session_stem(path) / 'plane0'
    ops_path   = plane_stem / 'ops.npy'
    bin_path   = plane_stem / 'data.bin'
    bin2_path  = plane_stem / 'data_chan2.bin'
    stat_path  = plane_stem / 'stat.npy'
    F_path     = plane_stem / 'F.npy'
    F2_path    = plane_stem / 'F_chan2.npy'
    Fneu_path  = plane_stem / 'Fneu.npy'
    F2neu_path = plane_stem / 'Fneu_chan2.npy'

    # folders to put processed data and optional single-session plots
    proc_path = ALL_SESSION_STEM / recname
    proc_path.mkdir(parents=True, exist_ok=True)

    proc_data_path = proc_path / 'processed_data'
    proc_data_path.mkdir(parents=True, exist_ok=True)
    figure_path = FIGURE_SESSION_STEM / recname

    if make_plots:
        figure_path.mkdir(parents=True, exist_ok=True)

    # load files
    stat = np.load(stat_path, allow_pickle=True)
    if 'inmerge' not in stat[0]:
        print('no merging detected; skipped')  # detect merging
        return

    ops = np.load(ops_path, allow_pickle=True).item()

    tot_frames = ops['nframes']

    # behaviour file
    # load beh file
    beh_path = BEH_STEM / f'{recname}.pkl'
    with open(beh_path, 'rb') as f:
        beh = pickle.load(f)

    run_frames = beh['run_onset_frames']
    rew_frames = beh['reward_frames']
    cue_frames = beh['start_cue_frames']

    # filtering
    run_frames = np.asarray(run_frames)
    run_frames = run_frames[run_frames != -1]  # filter out the no-clear-onset trials

    rew_frames = np.asarray(rew_frames)
    rew_frames = rew_frames[rew_frames != -1]  # filter out the no-rew trials

    cue_frames = np.asarray(cue_frames)

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
    # now we calculate using a rolling percentile, 4 Feb 2026
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
    np.save(proc_data_path / 'F_dFF.npy', F_dFF)
    np.save(proc_data_path / 'F2_dFF.npy', F2_dFF)

    print(f'dF/F computed and saved ({timedelta(seconds=int(time()-t0))})')
    # ----------------------
    # dF/F calculation ends
    # ----------------------

    # RO-aligned extraction
    if make_plots:
        run_path = figure_path / 'ro_aligned_single_roi'
        run_path.mkdir(parents=True, exist_ok=True)

    # valid ROIs RO dict
    all_run_dict = {}
    all_mean_run_dict = {}
    all_run_ch2_dict = {}
    all_mean_run_ch2_dict = {}

    # constituent ROIs RO dict
    all_run_const_dict = {}
    all_mean_run_const_dict = {}
    all_run_const_ch2_dict = {}
    all_mean_run_const_ch2_dict = {}

    for roi in tqdm(all_rois, desc='Aligning to run onset'):
        ca = F_dFF[roi]
        ref_ca = F2_dFF[roi]

        run_aligned, run_aligned_im, run_aligned_ch2, _, tot_run = (
            support.align_dual_channel_traces(
                ca,
                ref_ca,
                run_frames,
                BEF,
                AFT,
                SAMP_FREQ,
                tot_frames,
                )
            )

        if roi in valid_rois:  # if ROI resulted from a final merge
            all_run_dict[f'ROI {roi}'] = run_aligned
            all_run_ch2_dict[f'ROI {roi}'] = run_aligned_ch2

            run_aligned_mean = np.mean(run_aligned, axis=0)
            run_aligned_sem = sem(run_aligned, axis=0)
            all_mean_run_dict[f'ROI {roi}'] = run_aligned_mean

            run_aligned_mean_ch2 = np.mean(run_aligned_ch2, axis=0)
            run_aligned_sem_ch2 = sem(run_aligned_ch2, axis=0)
            all_mean_run_ch2_dict[f'ROI {roi}'] = run_aligned_mean_ch2

            if make_plots:
                support.save_roi_alignment_plots(
                    ref_im=ref_im,
                    stat=stat,
                    roi=roi,
                    aligned_im=run_aligned_im,
                    aligned_mean=run_aligned_mean,
                    aligned_sem=run_aligned_sem,
                    aligned_mean_ch2=run_aligned_mean_ch2,
                    aligned_sem_ch2=run_aligned_sem_ch2,
                    xaxis=XAXIS,
                    before=BEF,
                    after=AFT,
                    total_events=tot_run,
                    out_dir=run_path,
                    event_title='run-onset-aligned',
                    xlabel='Time from run-onset (s)',
                    overview_xlabel='Time from run onset (s)',
                    dual_axis_ch2=False,
                    )

        else:  # if ROI is constituent to a valid ROI
            all_run_const_dict[f'ROI {roi}'] = run_aligned
            all_run_const_ch2_dict[f'ROI {roi}'] = run_aligned_ch2

    # save valid ROI data
    valid_rois_dict = {f'ROI {roi}': valid_rois_dict[roi]
                       for roi in valid_rois_dict}  # rename keys to align with other dicts
    np.save(rf'{proc_data_path}\valid_ROIs_dict.npy', valid_rois_dict)
    np.save(rf'{proc_data_path}\valid_ROIs_coord_dict.npy', valid_rois_coord_dict)
    np.save(rf'{proc_data_path}\constituent_ROIs_coord_dict.npy', constituent_rois_coord_dict)

    # save merged ROIs
    np.save(proc_path / f'{recname}_all_run.npy', all_run_dict)
    np.save(proc_path / f'{recname}_all_run_mean.npy', all_mean_run_dict)
    np.save(proc_path / f'{recname}_all_run_ch2.npy', all_run_ch2_dict)
    np.save(proc_path / f'{recname}_all_run_ch2_mean.npy', all_mean_run_ch2_dict)

    # save constituent ROIs
    np.save(rf'{proc_data_path}\RO_aligned_const_dict.npy', all_run_const_dict)
    np.save(rf'{proc_data_path}\RO_aligned_const_mean_dict.npy', all_mean_run_const_dict)
    np.save(rf'{proc_data_path}\RO_aligned_const_ch2_dict.npy', all_run_const_ch2_dict)
    np.save(rf'{proc_data_path}\RO_aligned_const_mean_ch2_dict.npy', all_mean_run_const_ch2_dict)

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

    for roi in tqdm(all_rois, desc='Aligning to rewards'):
        ca = F_dFF[roi]
        ref_ca = F2_dFF[roi]

        rew_aligned, rew_aligned_im, rew_aligned_ch2, _, tot_rew = (
            support.align_dual_channel_traces(
                ca,
                ref_ca,
                rew_frames,
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
                    total_events=tot_rew,
                    out_dir=rew_path,
                    event_title='rew-aligned',
                    xlabel='time from rew (s)',
                    dual_axis_ch2=True,
                    )

        else:  # if ROI is constituent to a valid ROI
            all_rew_const_dict[f'ROI {roi}'] = rew_aligned
            all_rew_const_ch2_dict[f'ROI {roi}'] = rew_aligned_ch2

    # save merged ROIs
    np.save(proc_path / f'{recname}_all_rew.npy', all_rew_dict)
    np.save(proc_path / f'{recname}_all_rew_mean.npy', all_mean_rew_dict)
    np.save(proc_path / f'{recname}_all_rew_ch2.npy', all_rew_ch2_dict)
    np.save(proc_path / f'{recname}_all_rew_ch2_mean.npy', all_mean_rew_ch2_dict)

    # save constituent ROIs
    np.save(rf'{proc_data_path}\rew_aligned_const_dict.npy', all_rew_const_dict)
    np.save(rf'{proc_data_path}\rew_aligned_const_mean_dict.npy', all_mean_rew_const_dict)
    np.save(rf'{proc_data_path}\rew_aligned_const_ch2_dict.npy', all_rew_const_ch2_dict)
    np.save(rf'{proc_data_path}\rew_aligned_const_mean_ch2_dict.npy', all_mean_rew_const_ch2_dict)

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

    for roi in tqdm(all_rois, desc='aligning to cues'):
        ca = F_dFF[roi]
        ref_ca = F2_dFF[roi]

        cue_aligned, cue_aligned_im, cue_aligned_ch2, _, tot_cue = (
            support.align_dual_channel_traces(
                ca,
                ref_ca,
                cue_frames,
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
                    dual_axis_ch2=True,
                    )

        else:  # if ROI is constituent to a valid ROI
            all_cue_const_dict[f'ROI {roi}'] = cue_aligned
            all_cue_const_ch2_dict[f'ROI {roi}'] = cue_aligned_ch2

    # save merged ROIs
    np.save(proc_path / f'{recname}_all_cue.npy', all_cue_dict)
    np.save(proc_path / f'{recname}_all_cue_mean.npy', all_mean_cue_dict)
    np.save(proc_path / f'{recname}_all_cue_ch2.npy', all_cue_ch2_dict)
    np.save(proc_path / f'{recname}_all_cue_ch2_mean.npy', all_mean_cue_ch2_dict)

    # save constituent ROIs
    np.save(rf'{proc_data_path}\cue_aligned_const_dict.npy', all_cue_const_dict)
    np.save(rf'{proc_data_path}\cue_aligned_const_mean_dict.npy', all_mean_cue_const_dict)
    np.save(rf'{proc_data_path}\cue_aligned_const_ch2_dict.npy', all_cue_const_ch2_dict)
    np.save(rf'{proc_data_path}\cue_aligned_const_mean_ch2_dict.npy', all_mean_cue_const_ch2_dict)

#%% main
def main(argv=None):
    parser = argparse.ArgumentParser(
        description='extract LC-HPC axon-GCaMP aligned ROI data.'
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

    session_paths = pathHPCLCGCaMP
    if args.recording_filter:
        session_paths = [
            rec_path for rec_path in pathHPCLCGCaMP
            if args.recording_filter in rec_path
            or args.recording_filter in Path(rec_path).name
        ]
        if not session_paths:
            raise ValueError(f'no recordings matched filter: {args.recording_filter}')

    for rec_path in session_paths:
        process_session(
            rec_path,
            make_plots=args.make_plots,
        )
        recname = Path(rec_path).name
        print(f'clearing CPU/GPU memory after {recname}...', flush=True)
        gc.collect()
        if GPU_AVAILABLE:
            cp.cuda.Stream.null.synchronize()
            cp.get_default_memory_pool().free_all_blocks()
            cp.get_default_pinned_memory_pool().free_all_blocks()

if __name__ == '__main__':
    main()
