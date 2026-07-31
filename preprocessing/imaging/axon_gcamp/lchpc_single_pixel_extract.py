# -*- coding: utf-8 -*-
'''
Created on Mon Mar 24 17:46:39 2025

extract single-pixel fluorescence traces after spatial filtering

@author: Dinghao Luo
'''

#%% imports
import argparse
from pathlib import Path

import numpy as np
from tqdm import tqdm
from time import time
from datetime import timedelta
import sys

repo_root = Path(__file__).resolve().parents[3]
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))
import imaging_pipeline_functions as ipf
import lchpc_axon_support as support

from common_functions import get_GPU_availability, mpl_formatting
from console_formatting import print_session
import project_paths as pp
mpl_formatting()

import rec_list
paths = rec_list.pathLCHPCGCaMP

# match the active-running axon extractor's rolling-percentile baseline window
DFF_WINDOW_SIZE = 9000


#%% GPU acceleration
_, GPU_AVAILABLE, device_name = get_GPU_availability()
print(f'GPU_AVAILABLE={GPU_AVAILABLE}')


#%% main
def main(argv=None):
    parser = argparse.ArgumentParser(
        description='extract LC-HPC axon-GCaMP single-pixel fluorescence traces.'
    )
    parser.add_argument(
        '--recording',
        dest='recording_filter',
        help='optional substring filter for recname or recording path',
    )
    parser.add_argument(
        '--median-chunk-size',
        type=int,
        default=500,
        help='frames per spatial median-filter chunk (default: 500)',
    )
    args = parser.parse_args(argv)

    session_paths = paths
    if args.recording_filter:
        session_paths = [
            rec_path for rec_path in paths
            if args.recording_filter in rec_path
            or args.recording_filter in Path(rec_path).name
        ]
        if not session_paths:
            raise ValueError(f'no recordings matched filter: {args.recording_filter}')

    for path in session_paths:
        recname = path[-17:]
        print_session(recname)
        t0 = time()

        plane_stem = pp.resolve_suite2p_session_stem(path) / 'plane0'
        ops_path = plane_stem / 'ops.npy'
        bin_path = plane_stem / 'data.bin'
        stat_path = plane_stem / 'stat.npy'

        # folder to put processed data
        proc_path = pp.LCHPC_AXON_STEM / 'all_sessions' / recname
        proc_data_path = proc_path / 'processed_data'
        out_path = proc_data_path / 'roi_pixel_dFF.npy'

        proc_path.mkdir(parents=True, exist_ok=True)
        proc_data_path.mkdir(parents=True, exist_ok=True)

        # load files
        stat = np.load(stat_path, allow_pickle=True)
        if 'inmerge' not in stat[0]:
            print('no merging detected; skipped')
            continue

        # get roi idx
        valid_rois_dict = support.filter_valid_rois(stat)
        valid_rois = [*valid_rois_dict]  # using lists because of the later serial comprehension

        # get x and y pix idx and put in a dict
        roi_coords_dict = {
            f'ROI {roi}': [stat[roi]['xpix'], stat[roi]['ypix']]
            for roi in valid_rois
            }

        # read in bin files
        ops = np.load(ops_path, allow_pickle=True).item()
        tot_frames = ops['nframes']
        shape = (tot_frames, ops['Ly'], ops['Lx'])
        mov = np.memmap(bin_path, mode='r', dtype='int16', shape=shape)

        # Median-filter spatially in chunks and keep only ROI pixels. This
        # avoids materialising the full filtered movie for long recordings.
        print('channel 1')
        roi_pixel_traces = support.extract_spatial_median_roi_traces(
            mov,
            roi_coords_dict,
            size=5,
            GPU_AVAILABLE=GPU_AVAILABLE,
            chunk_size=args.median_chunk_size,
            )

        roi_pixels_dict = {}
        print('extracting pixel traces...')
        for roi in tqdm(roi_coords_dict, desc='single-pixel ROIs'):
            t1 = time()
            pixx = roi_coords_dict[roi][0]
            pixy = roi_coords_dict[roi][1]

            print(f'{roi}: calculating pixel dF/F...', flush=True)
            pixels_traces = roi_pixel_traces.pop(roi)
            pixels_dFF = ipf.calculate_dFF_percentile(
                pixels_traces,
                t_axis=0,
                window_size=DFF_WINDOW_SIZE,
                GPU_AVAILABLE=GPU_AVAILABLE,
                device_name=device_name,
                progress_desc=f'{roi} pixel dF/F baseline',
                )

            curr_roi_dict = {
                'coord': list(zip(pixx, pixy)),
                'dFF': pixels_dFF
                }
            roi_pixels_dict[roi] = curr_roi_dict

            print(f'{roi} done ({timedelta(seconds=int(time()-t1))})')

        np.save(out_path, roi_pixels_dict)
        print(f'{recname} done ({timedelta(seconds=int(time()-t0))})')
if __name__ == '__main__':
    main()
