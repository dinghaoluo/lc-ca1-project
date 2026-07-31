# -*- coding: utf-8 -*-
'''
Created on Tue Jun  2 2026

run HPC profiles, waveforms, first-lick, stimulation, sequence, and latency analyses

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
        'key': 'profiles',
        'label': 'HPC profiles',
        'branch': 'core',
        'script': Path('analysis') / 'hpc' / 'hpc_all_profiles.py',
        'profile_stage': True,
        },
    {
        'key': 'profiles_raphi',
        'label': 'Raphi HPC profiles',
        'branch': 'core',
        'script': Path('analysis') / 'hpc' / 'hpc_all_profiles_raphi.py',
        'profile_stage': True,
        },
    {
        'key': 'waveforms',
        'label': 'HPC waveforms',
        'branch': 'core',
        'script': Path('analysis') / 'hpc' / 'hpc_all_waveforms.py',
        },
    {
        'key': 'lick_comp',
        'label': 'HPC opto lick comparison',
        'branch': 'analyses',
        'script': Path('analysis') / 'hpc' / 'behaviour' / 'hpc_lc_stim_lick_comp.py',
        },
    {
        'key': 'early_late_lick',
        'label': 'Early/late lick profiles',
        'branch': 'analyses',
        'script': Path('analysis') / 'hpc' / 'first_lick_analysis' / 'hpc_pyr_early_vs_late_lick_profiles_full.py',
        },
    {
        'key': 'stim_corr',
        'label': 'LC-stim CA1 behaviour correlation',
        'branch': 'analyses',
        'script': Path('analysis') / 'hpc' / 'stim_ctrl' / 'analyse_lc_stim_ca1_pyr_up_down_behaviour_corr.py',
        },
    {
        'key': 'plot_stim_corr',
        'label': 'Plot LC-stim CA1 behaviour correlation',
        'branch': 'analyses',
        'script': Path('analysis') / 'hpc' / 'stim_ctrl' / 'plot_lc_stim_ca1_pyr_up_down_behaviour_corr.py',
        },
    {
        'key': 'sequence_metrics',
        'label': 'CA1 place/time sequence metrics',
        'branch': 'analyses',
        'script': Path('analysis') / 'hpc' / 'sequence_analysis' / 'build_hpc_sequence_place_time_metrics.py',
        },
    {
        'key': 'plot_sequence_metrics',
        'label': 'Plot CA1 place/time sequence metrics',
        'branch': 'analyses',
        'script': Path('plotting') / 'hpc' / 'plot_hpc_sequence_place_time.py',
        },
    {
        'key': 'latency',
        'label': 'LC-stim CA1 effect latency',
        'branch': 'analyses',
        'script': Path('analysis') / 'hpc' / 'stim_ctrl' / 'estimate_lc_stim_ca1_effect_latency.py',
        'latency_stage': True,
        },
    {
        'key': 'plot_latency',
        'label': 'Plot LC-stim CA1 effect latency',
        'branch': 'analyses',
        'script': Path('analysis') / 'hpc' / 'stim_ctrl' / 'plot_lc_stim_ca1_effect_latency.py',
        },
    ]
STAGE_KEYS = [stage['key'] for stage in STAGE_DEFS]


#%% main
def main(argv=None):
    parser = argparse.ArgumentParser(
        description='run HPC processing pipeline'
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
        choices=['all', 'core', 'analyses'],
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
    parser.add_argument(
        '--n-workers',
        '--workers',
        dest='n_workers',
        type=int,
        default=None,
        help='number of workers for profile stages',
        )
    parser.add_argument(
        '--single-cell-latency-plots',
        action='store_true',
        help='export per-cell latency plots',
        )
    args = parser.parse_args(argv)
    selected_defs = select_stage_defs(STAGE_DEFS, args, parser)
    stages = []
    for stage_def in selected_defs:
        script_args = []
        if stage_def.get('profile_stage'):
            if args.n_workers is not None:
                script_args.extend(['--n-workers', str(args.n_workers)])
        if stage_def.get('latency_stage') and args.single_cell_latency_plots:
            script_args.append('--single-cell-plots')

        stages.append({
            'key': stage_def['key'],
            'label': stage_def['label'],
            'command': python_command(repo_root, stage_def['script'], script_args),
            'display': display_command(stage_def['script'], script_args),
        })

    run_pipeline(
        stages,
        repo_root,
        'HPC pipeline',
        verbose=args.verbose,
        dry_run=args.dry_run,
        )

if __name__ == '__main__':
    main(sys.argv[1:])
