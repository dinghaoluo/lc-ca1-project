# -*- coding: utf-8 -*-
'''
Created on 6 Jul 2026

run CA1 place-cell, time-cell, and sequence analyses

@author: Dinghao Luo
'''

from __future__ import annotations

#%% imports
import argparse
from pathlib import Path
import sys

repo_root = Path(__file__).resolve().parents[3]

if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))

from pipeline_runner_functions import display_command, python_command, run_pipeline


#%% stages
STAGE_DEFS = [
    {
        'key': 'metrics',
        'label': 'CA1 place/time sequence metrics',
        'script': (
            Path('analysis')
            / 'hpc'
            / 'sequence_analysis'
            / 'build_hpc_sequence_place_time_metrics.py'
        ),
    },
    {
        'key': 'plots',
        'label': 'CA1 place/time sequence plots',
        'script': Path('plotting') / 'hpc' / 'plot_hpc_sequence_place_time.py',
    },
]
STAGE_KEYS = [stage['key'] for stage in STAGE_DEFS]


#%% main
def main(argv=None):
    parser = argparse.ArgumentParser(
        description='run CA1 place/time sequence analysis pipeline'
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
        '--only',
        choices=STAGE_KEYS,
        action='append',
        default=[],
        help='run only this stage; can be used more than once',
    )
    parser.add_argument(
        '--from-stage',
        choices=STAGE_KEYS,
        help='start at this stage',
    )
    parser.add_argument(
        '--n-shuf',
        type=int,
        default=None,
        help='pass through to the metrics stage',
    )
    parser.add_argument(
        '--recording',
        help='pass through to the metrics stage',
    )
    parser.add_argument(
        '--max-recordings',
        type=int,
        help='pass through to the metrics stage',
    )
    args = parser.parse_args(argv)
    selected = list(STAGE_DEFS)
    selected_keys = [stage['key'] for stage in selected]
    if args.from_stage is not None:
        selected = selected[selected_keys.index(args.from_stage):]

    if args.only:
        wanted = set(args.only)
        selected = [
            stage
            for stage in selected
            if stage['key'] in wanted
        ]

    if len(selected) == 0:
        parser.error('no stages selected')
    stages = []
    for stage in selected:
        script_args = []
        if stage['key'] == 'metrics':
            if args.n_shuf is not None:
                script_args.extend(['--n-shuf', str(args.n_shuf)])
            if args.recording:
                script_args.extend(['--recording', args.recording])
            if args.max_recordings is not None:
                script_args.extend(['--max-recordings', str(args.max_recordings)])

        stages.append({
            'key': stage['key'],
            'label': stage['label'],
            'command': python_command(repo_root, stage['script'], script_args),
            'display': display_command(stage['script'], script_args),
        })
    run_pipeline(
        stages,
        repo_root,
        'CA1 sequence analysis pipeline',
        verbose=args.verbose,
        dry_run=args.dry_run,
    )

if __name__ == '__main__':
    main(sys.argv[1:])
