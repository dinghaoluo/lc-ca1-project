# -*- coding: utf-8 -*-
'''
regenerate GRABNE significant-ROI figures from saved extraction outputs.

@author: Dinghao Luo
'''

#%% imports
import argparse
from pathlib import Path
import sys

import numpy as np


repo_root = Path(__file__).resolve().parent
while not (repo_root / 'utils').exists() and repo_root != repo_root.parent:
    repo_root = repo_root.parent
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))

script_dir = Path(__file__).resolve().parent
if str(script_dir) not in sys.path:
    sys.path.insert(0, str(script_dir))

from hpc_extract_significant_roi import (
    PLOT_PAYLOAD_STEM,
    get_recording_paths,
    save_sig_roi_plots,
    )


#%% command line
def build_arg_parser():
    parser = argparse.ArgumentParser(
        description='plot GRABNE significant-ROI figures from saved extraction outputs.'
        )
    parser.add_argument(
        '--recording',
        dest='recording_filter',
        help='optional substring filter for recname or recording path',
        )
    return parser


#%% main
def main(argv=None):
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    for rec_path in get_recording_paths(args.recording_filter):
        recname = rec_path[-17:]
        print(f'\nplotting {recname}')
        payload = np.load(
            PLOT_PAYLOAD_STEM / f'{recname}_sig_roi_plot_payload.npy',
            allow_pickle=True
            ).item()
        save_sig_roi_plots(recname, payload)


if __name__ == '__main__':
    main()
