# -*- coding: utf-8 -*-
'''
Created on Tue Jun  2 2026

run the behaviour, ephys, imaging, and saved-output plotting stages

the imaging branch descends from the Suite2p/grid wrapper created on
13 May 2024 and modified on 25 June 2024 by Dinghao Luo and Jingyu Cao

@author: Dinghao Luo
'''


#%% imports
import argparse
from pathlib import Path
import sys

repo_root = Path(__file__).resolve().parents[1]

if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))

from pipeline_runner_functions import display_command, python_command, run_pipeline, select_stage_defs


#%% stages
STAGE_DEFS = [
    {
        'key': 'behaviour',
        'label': 'Behaviour preprocessing',
        'branch': 'behaviour',
        'script': Path('preprocessing') / 'behaviour' / 'process_behaviour.py',
        'args': ['--all'],
        },
    {
        'key': 'lc_ephys',
        'label': 'LC ephys extraction',
        'branch': 'ephys',
        'script': Path('preprocessing') / 'ephys' / 'lc' / 'lc_all_extract.py',
        },
    {
        'key': 'hpc_ephys',
        'label': 'HPC ephys extraction',
        'branch': 'ephys',
        'script': Path('preprocessing') / 'ephys' / 'hpc' / 'hpc_all_extract.py',
        },
    {
        'key': 'hpc_ephys_raphi',
        'label': 'Raphi HPC ephys extraction',
        'branch': 'ephys',
        'script': Path('preprocessing') / 'ephys' / 'hpc' / 'hpc_all_extract_raphi.py',
        },
    {
        'key': 'suite2p_registration',
        'label': 'Suite2p registration',
        'branch': 'imaging',
        'script': Path('preprocessing') / 'imaging' / 'suite2p_registration.py',
        },
    {
        'key': 'axon_gcamp',
        'label': 'Axon GCaMP extraction',
        'branch': 'imaging',
        'script': Path('preprocessing') / 'imaging' / 'axon_gcamp' / 'lchpc_axon_all_extract.py',
        },
    {
        'key': 'axon_gcamp_immobile',
        'label': 'Immobile axon GCaMP extraction',
        'branch': 'imaging',
        'script': Path('preprocessing') / 'imaging' / 'axon_gcamp' / 'lchpc_axon_all_extract_immobile.py',
        },
    {
        'key': 'axon_gcamp_pixels',
        'label': 'Axon GCaMP single-pixel extraction',
        'branch': 'imaging',
        'script': Path('preprocessing') / 'imaging' / 'axon_gcamp' / 'lchpc_single_pixel_extract.py',
        },
    {
        'key': 'dlight_lc_opto',
        'label': 'dLight LC-opto extraction',
        'branch': 'imaging',
        'script': Path('preprocessing') / 'imaging' / 'dlight' / 'hpc_dlight_lc_opto_extract.py',
        },
    {
        'key': 'grabne_lc_opto',
        'label': 'GRABNE LC-opto extraction',
        'branch': 'imaging',
        'script': Path('preprocessing') / 'imaging' / 'grabne' / 'hpc_grabne_lc_opto_extract.py',
        },
    {
        'key': 'grabne_tone',
        'label': 'GRABNE tone extraction',
        'branch': 'imaging',
        'script': Path('preprocessing') / 'imaging' / 'grabne' / 'hpc_grabne_tone_extract.py',
        },
    {
        'key': 'nlight_lc_opto',
        'label': 'nLight LC-opto extraction',
        'branch': 'imaging',
        'script': Path('preprocessing') / 'imaging' / 'nlight' / 'hpc_nlight_lc_opto_extract.py',
        },
    {
        'key': 'plot_axon_gcamp',
        'label': 'Plot axon GCaMP extraction',
        'branch': 'imaging_plots',
        'script': Path('preprocessing') / 'imaging' / 'axon_gcamp' / 'plot_lchpc_axon_all_extract.py',
        },
    {
        'key': 'plot_axon_gcamp_immobile',
        'label': 'Plot immobile axon GCaMP extraction',
        'branch': 'imaging_plots',
        'script': Path('preprocessing') / 'imaging' / 'axon_gcamp' / 'plot_lchpc_axon_all_extract_immobile.py',
        },
    {
        'key': 'plot_dlight_lc_opto',
        'label': 'Plot dLight LC-opto extraction',
        'branch': 'imaging_plots',
        'script': Path('preprocessing') / 'imaging' / 'dlight' / 'plot_hpc_dlight_lc_opto_extract.py',
        },
    {
        'key': 'plot_grabne_lc_opto',
        'label': 'Plot GRABNE LC-opto extraction',
        'branch': 'imaging_plots',
        'script': Path('preprocessing') / 'imaging' / 'grabne' / 'plot_hpc_grabne_lc_opto_extract.py',
        },
    {
        'key': 'plot_grabne_tone',
        'label': 'Plot GRABNE tone extraction',
        'branch': 'imaging_plots',
        'script': Path('preprocessing') / 'imaging' / 'grabne' / 'plot_hpc_grabne_tone_extract.py',
        },
    {
        'key': 'plot_nlight_lc_opto',
        'label': 'Plot nLight LC-opto extraction',
        'branch': 'imaging_plots',
        'script': Path('preprocessing') / 'imaging' / 'nlight' / 'plot_hpc_nlight_lc_opto_extract.py',
        },
    ]
STAGE_KEYS = [stage['key'] for stage in STAGE_DEFS]


#%% main
def main(argv=None):
    parser = argparse.ArgumentParser(
        description='run preprocessing pipeline'
    )
    parser.add_argument(
        '-v', '--verbose', '-verbose', action='store_true',
        help='print output from each stage',
    )
    parser.add_argument(
        '--dry-run', action='store_true',
        help='print commands without running them',
    )
    parser.add_argument(
        '--branch',
        choices=['all', 'behaviour', 'ephys', 'imaging', 'imaging_plots'],
        default='all',
        help='pipeline branch to run',
    )
    parser.add_argument(
        '--from-stage', choices=STAGE_KEYS,
        help='start at this stage',
    )
    parser.add_argument(
        '--only', choices=STAGE_KEYS, action='append', default=[],
        help='run only this stage; can be used more than once',
    )
    args = parser.parse_args(argv)
    selected_defs = select_stage_defs(STAGE_DEFS, args, parser)
    stages = [
        {
            'key': stage_def['key'],
            'label': stage_def['label'],
            'command': python_command(
                repo_root,
                stage_def['script'],
            ),
            'display': display_command(stage_def['script']),
        }
        for stage_def in selected_defs
        ]

    run_pipeline(
        stages,
        repo_root,
        'Preprocessing pipeline',
        verbose=args.verbose,
        dry_run=args.dry_run,
        )

if __name__ == '__main__':
    main(sys.argv[1:])
