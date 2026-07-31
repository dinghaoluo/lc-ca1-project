# -*- coding: utf-8 -*-
'''
Created on Tue Jun  2 2026

run electrophysiology and imaging optogenetic summaries

@author: Dinghao Luo
'''


#%% imports
import argparse
from pathlib import Path
import sys

repo_root = Path(__file__).resolve().parents[2]

if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))

from pipeline_runner_functions import display_command, python_command, run_pipeline, select_stage_defs


#%% stages
STAGE_DEFS = [
    {
        'key': 'ephys_summary',
        'label': 'Ephys optogenetics summary',
        'branch': 'ephys',
        'script': Path('analysis') / 'optogenetics' / 'summarise_opto.py',
        },
    {
        'key': 'imaging_summary',
        'label': 'Imaging-animal optogenetics summary',
        'branch': 'imaging',
        'script': Path('analysis') / 'imaging' / 'optogenetics' / 'summarise_opto_imaging.py',
        },
    ]
STAGE_KEYS = [stage['key'] for stage in STAGE_DEFS]


#%% main
def main(argv=None):
    parser = argparse.ArgumentParser(
        description='run optogenetics summary pipeline'
    )
    parser.add_argument(
        '-v',
        '--verbose',
        '-verbose',
        action='store_true',
        help='print output from each stage',
    )
    parser.add_argument(
        '--dry-run',
        action='store_true',
        help='print commands without running them',
    )
    parser.add_argument(
        '--branch',
        dest='branch',
        choices=['all', 'ephys', 'imaging'],
        default='all',
        help='pipeline branch to run',
    )
    parser.add_argument(
        '--from-stage',
        choices=STAGE_KEYS,
        help='start at this stage',
    )
    parser.add_argument(
        '--only',
        choices=STAGE_KEYS,
        action='append',
        default=[],
        help='run only this stage; can be used more than once',
    )
    args = parser.parse_args(argv)
    selected_defs = select_stage_defs(STAGE_DEFS, args, parser)
    stages = [
        {
            'key': stage_def['key'],
            'label': stage_def['label'],
            'command': python_command(repo_root, stage_def['script']),
            'display': display_command(stage_def['script']),
        }
        for stage_def in selected_defs
        ]

    run_pipeline(
        stages,
        repo_root,
        'Optogenetics summary pipeline',
        verbose=args.verbose,
        dry_run=args.dry_run,
        )

if __name__ == '__main__':
    main(sys.argv[1:])
