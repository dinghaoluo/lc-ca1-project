# -*- coding: utf-8 -*-
'''
Created on Tue Jun  2 2026

run behaviour summaries and LC/HPC optogenetic controls

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
        'key': 'speed_licks',
        'label': 'Speed and lick profiles',
        'branch': 'core',
        'script': Path('analysis') / 'behaviour' / 'analyse_speed_licks.py',
        },
    {
        'key': 'pupil_size',
        'label': 'Pupil size',
        'branch': 'core',
        'script': Path('analysis') / 'behaviour' / 'analyse_pupil_size.py',
        },
    {
        'key': 'first_lick',
        'label': 'First lick since last reward',
        'branch': 'core',
        'script': Path('analysis') / 'behaviour' / 'first_lick_since_last_reward.py',
        },
    {
        'key': 'off_target_bouts',
        'label': 'Off-target run bouts',
        'branch': 'core',
        'script': Path('analysis') / 'behaviour' / 'off_target_run_bouts.py',
        },
    {
        'key': 'hpc_opto_speed',
        'label': 'HPC opto speed controls',
        'branch': 'controls',
        'script': Path('analysis') / 'behaviour_control' / 'hpc_opto_speed_controls.py',
        },
    {
        'key': 'lc_controls',
        'label': 'LC controls',
        'branch': 'controls',
        'script': Path('analysis') / 'behaviour_control' / 'lc_controls.py',
        },
    {
        'key': 'lc_opto_speed',
        'label': 'LC opto speed controls',
        'branch': 'controls',
        'script': Path('analysis') / 'behaviour_control' / 'lc_opto_speed_controls.py',
        },
    ]
STAGE_KEYS = [stage['key'] for stage in STAGE_DEFS]


#%% main
def main(argv=None):
    parser = argparse.ArgumentParser(
        description='run behaviour/control pipeline'
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
        choices=['all', 'core', 'controls'],
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
        'Behaviour/control pipeline',
        verbose=args.verbose,
        dry_run=args.dry_run,
        )

if __name__ == '__main__':
    main(sys.argv[1:])
