# -*- coding: utf-8 -*-
"""
Created on Mon 13 May 17:13:42 2024
Modified on Tue 25 June 15:09:40 2024

This code combines grid_extract.py (Dinghao) and after_suite2p.py (Jingyu)

@authors: Dinghao Luo, Jingyu Cao
@modifiers: Dinghao Luo, Jingyu Cao
    - removed dead external path bootstraps and reused the shared gpu helper
"""

#%% imports
import argparse
import sys
from datetime import timedelta
from pathlib import Path
from time import time

import pandas as pd

repo_root = Path(__file__).resolve().parent
while not (repo_root / 'utils').exists() and repo_root != repo_root.parent:
    repo_root = repo_root.parent
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))
if str(repo_root / '_archive' / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / '_archive' / 'utils'))

import imaging_pipeline_functions as ipf
import imaging_pipeline_main_functions as ipmf
import project_paths as pp
from common_functions import get_GPU_availability

import rec_list


#%% constants
ROI_MODES = {
    '1': 'suite2p',
    '2': 'grid',
    'suite2p': 'suite2p',
    'grid': 'grid',
}

BEHAVIOUR_DF_PATH = pp.behaviour_session_pickle('all_HPCGRABNE_sessions.pkl')

DEFAULT_CONFIG = {
    'plot_ref': True,
    'align_run': 1,
    'align_rew': 1,
    'align_cue': 0,
    'smooth': 1,
    'bef': 1,
    'aft': 4,
    'dff': 1,
    'plot_heatmap': 1,
    'plot_trace': 1,
    'stride': 496,
    'border': 8,
    'save_grids': 1,
}

ROI_PROMPT = 'process with...\n1: suite2p ROIs\n2: grid ROIs\n'


#%% helper functions
def build_arg_parser():
    parser = argparse.ArgumentParser(
        description='run the hippocampal imaging extraction pipeline.'
    )
    parser.add_argument(
        '--roi-mode',
        dest='roi_mode',
        choices=['1', '2', 'suite2p', 'grid'],
        help='roi extraction mode; use 1/suite2p or 2/grid',
    )
    parser.add_argument(
        '--recording',
        dest='recording_filter',
        help='optional substring filter for recname or recording path',
    )
    return parser


def prompt_for_roi_mode():
    return input(ROI_PROMPT).strip()


def normalise_roi_mode(selection):
    if selection not in ROI_MODES:
        raise ValueError('invalid roi mode; only 1, 2, suite2p, and grid are supported')
    return ROI_MODES[selection]


def build_runtime_config(roi_mode):
    config = dict(DEFAULT_CONFIG)
    config['roi_mode'] = roi_mode
    if roi_mode == 'grid':
        config['fit'] = ipf.check_stride_border(config['stride'], config['border'])
    return config


def format_run_summary(config):
    if config['roi_mode'] == 'suite2p':
        return f'''
processing suite2p ROIs...
    align_run = {config['align_run']}
    align_rew = {config['align_rew']}
    align_cue = {config['align_cue']}
    bef = {config['bef']}
    aft = {config['aft']}
    smooth = {config['smooth']}
'''.strip()

    return f'''
processing grid ROIs...
    align_run = {config['align_run']}
    align_rew = {config['align_rew']}
    align_cue = {config['align_cue']}
    bef = {config['bef']}
    aft = {config['aft']}
    stride = {config['stride']}
    border = {config['border']}
    smooth = {config['smooth']}
    dFF = {config['dff']}
'''.strip()


def load_behaviour_dataframe():
    print('loading behaviour dataframe...')
    try:
        return pd.read_pickle(BEHAVIOUR_DF_PATH)
    except FileNotFoundError:
        print('loading failed: no behavioural dataframe found\n')
        return None


def resolve_recording_context(rec_path):
    if 'Dinghao' in rec_path:
        reg_path = rec_path + r'\processed\suite2p\plane0'
        recname = rec_path[-17:]
        txt_path = str(pp.MICEEXP_ROOT / f'ANMD{recname[1:4]}' / f'{recname[:4]}{recname[5:]}T.txt')
        return recname, reg_path, txt_path

    if 'Jingyu' in rec_path:
        reg_path = rec_path + r'\RegOnly\suite2p\plane0'
        recname = rec_path[-17:-3] + '-' + rec_path[-2:]
        txt_path = r'Z:\Jingyu\mice-expdata\{}\A{}T.txt'.format(
            rec_path[-23:-18],
            recname[2:],
        )
        return recname, reg_path, txt_path

    raise ValueError(f'unrecognised recording root for path: {rec_path}')


def get_recording_paths(recording_filter=None):
    paths = rec_list.pathHPCGRABNE
    if not recording_filter:
        return paths

    filtered_paths = [
        path for path in paths
        if recording_filter in path or recording_filter in path[-17:]
    ]
    if not filtered_paths:
        raise ValueError(f'no recordings matched filter: {recording_filter}')
    return filtered_paths


def run_recording_pipeline(rec_path, config, behaviour_df, gpu_available):
    recname, reg_path, txt_path = resolve_recording_context(rec_path)
    print(f'\nprocessing {recname}')
    start = time()

    beh = []
    if behaviour_df is not None:
        try:
            print('reading session behavioural data...')
            beh = behaviour_df.loc[recname]
        except KeyError:
            beh = []

    if config['roi_mode'] == 'suite2p':
        ipmf.run_suite2p_pipeline(
            rec_path,
            recname,
            reg_path,
            txt_path,
            config['plot_ref'],
            config['plot_heatmap'],
            config['plot_trace'],
            config['smooth'],
            config['dff'],
            config['bef'],
            config['aft'],
            config['align_run'],
            config['align_rew'],
            config['align_cue'],
        )
    else:
        ipmf.run_grid_pipeline(
            rec_path,
            recname,
            reg_path,
            txt_path,
            beh,
            config['stride'],
            config['border'],
            config['plot_ref'],
            config['smooth'],
            config['dff'],
            config['save_grids'],
            config['bef'],
            config['aft'],
            config['align_run'],
            config['align_rew'],
            config['align_cue'],
            gpu_available,
        )

    print(f'{recname} done ({timedelta(seconds=int(time() - start))})')


def main(argv=None):
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    roi_selection = args.roi_mode or prompt_for_roi_mode()
    roi_mode = normalise_roi_mode(roi_selection)
    config = build_runtime_config(roi_mode)

    cp, gpu_available, device_name = get_GPU_availability()
    _ = cp, device_name

    print(format_run_summary(config))
    behaviour_df = load_behaviour_dataframe()

    for rec_path in get_recording_paths(args.recording_filter):
        run_recording_pipeline(rec_path, config, behaviour_df, gpu_available)


if __name__ == '__main__':
    main()
