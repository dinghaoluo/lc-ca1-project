# -*- coding: utf-8 -*-
'''
Created on Tue Jun  2 2026

run prazosin, propranolol, and SCH23390 summaries

@author: Dinghao Luo
'''


#%% imports
import argparse
from pathlib import Path
import sys

repo_root = Path(__file__).resolve().parents[2]

if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))

from pipeline_runner_functions import display_command, python_command, run_pipeline


#%% stages
STAGE_DEFS = [
    {
        'key': 'drug_summary',
        'label': 'Pharmacology summary',
        'script': Path('analysis') / 'pharmacology' / 'summarise_all_drugs.py',
        },
    ]

#%% main
def main(argv=None):
    parser = argparse.ArgumentParser(description='run pharmacology pipeline')
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
    args = parser.parse_args(argv)
    stages = [
        {
            'key': stage_def['key'],
            'label': stage_def['label'],
            'command': python_command(repo_root, stage_def['script']),
            'display': display_command(stage_def['script']),
        }
        for stage_def in STAGE_DEFS
        ]

    run_pipeline(
        stages,
        repo_root,
        'Pharmacology pipeline',
        verbose=args.verbose,
        dry_run=args.dry_run,
        )

if __name__ == '__main__':
    main(sys.argv[1:])
