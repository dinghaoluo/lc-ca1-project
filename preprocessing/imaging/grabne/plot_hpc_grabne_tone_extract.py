# -*- coding: utf-8 -*-
'''
Created on 19 May 2026

plot saved GRABNE tone-aligned data

@author: Dinghao Luo
'''

#%% imports
import argparse
from pathlib import Path
import sys

import numpy as np

repo_root = Path(__file__).resolve().parents[3]
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))

script_dir = Path(__file__).resolve().parent
if str(script_dir) not in sys.path:
    sys.path.insert(0, str(script_dir))

from common_functions import mpl_formatting
from hpc_grabne_tone_extract import (
    all_sess_stem,
    all_sess_fig_stem,
    save_tone_figures,
    )
mpl_formatting()

import rec_list


#%% paths
paths = rec_list.pathGRABNETone + rec_list.pathGRABNEToneDbhBlock


#%% main
def main(argv=None):
    parser = argparse.ArgumentParser(
        description='plot GRABNE tone-locked figures from saved extraction outputs.'
        )
    parser.add_argument(
        '--recording',
        dest='recording_filter',
        help='optional substring filter for recname or recording path',
        )
    args = parser.parse_args(argv)

    selected_paths = paths
    if args.recording_filter:
        selected_paths = [
            path for path in paths
            if (args.recording_filter in path
                or args.recording_filter in Path(path).name)
            ]
        if not selected_paths:
            raise ValueError(f'no recordings matched filter: {args.recording_filter}')

    for path in selected_paths:
        recname = Path(path).name
        savepath = all_sess_stem / recname
        figurepath = all_sess_fig_stem / recname
        proc_data_path = savepath / 'processed_data'
        print(f'\nplotting {recname}')

        payload = np.load(
            proc_data_path / f'{recname}_tone_plot_payload.npy',
            allow_pickle=True
            ).item()
        save_tone_figures(figurepath, recname, payload)

if __name__ == '__main__':
    main()
