# -*- coding: utf-8 -*-
'''
Created on Tue Jun  2 2026

run behaviour, LC/HPC cell-summary, and LC-HPC axon-profile analyses

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
        'key': 'speed_licks',
        'label': 'Speed and lick profiles',
        'branch': 'behaviour',
        'script': Path('analysis') / 'behaviour' / 'analyse_speed_licks.py',
        },
    {
        'key': 'pupil_size',
        'label': 'Pupil size',
        'branch': 'behaviour',
        'script': Path('analysis') / 'behaviour' / 'analyse_pupil_size.py',
        },
    {
        'key': 'first_lick',
        'label': 'First lick since last reward',
        'branch': 'behaviour',
        'script': Path('analysis') / 'behaviour' / 'first_lick_since_last_reward.py',
        },
    {
        'key': 'off_target_bouts',
        'label': 'Off-target run bouts',
        'branch': 'behaviour',
        'script': Path('analysis') / 'behaviour' / 'off_target_run_bouts.py',
        },
    {
        'key': 'lc_waveforms',
        'label': 'LC waveforms and ACGs',
        'branch': 'lc',
        'script': Path('analysis') / 'lc' / 'lc_all_waveforms_acgs.py',
        },
    {
        'key': 'lc_spikes',
        'label': 'LC spike times and ISIs',
        'branch': 'lc',
        'script': Path('analysis') / 'lc' / 'lc_all_spikes_isis.py',
        },
    {
        'key': 'lc_identity',
        'label': 'LC identity UMAP',
        'branch': 'lc',
        'script': Path('analysis') / 'lc' / 'lc_all_identity_umap.py',
        },
    {
        'key': 'lc_profiles',
        'label': 'LC cell profiles',
        'branch': 'lc',
        'script': Path('analysis') / 'lc' / 'lc_all_profiles.py',
        'profile_stage': True,
        },
    {
        'key': 'hpc_profiles',
        'label': 'HPC profiles',
        'branch': 'hpc',
        'script': Path('analysis') / 'hpc' / 'hpc_all_profiles.py',
        'profile_stage': True,
        },
    {
        'key': 'hpc_profiles_raphi',
        'label': 'Raphi HPC profiles',
        'branch': 'hpc',
        'script': Path('analysis') / 'hpc' / 'hpc_all_profiles_raphi.py',
        'profile_stage': True,
        },
    {
        'key': 'hpc_waveforms',
        'label': 'HPC waveforms',
        'branch': 'hpc',
        'script': Path('analysis') / 'hpc' / 'hpc_all_waveforms.py',
        },
    {
        'key': 'imaging_profiles',
        'label': 'Axon GCaMP profiles',
        'branch': 'imaging',
        'script': Path('analysis') / 'imaging' / 'lchpc_axon_all_profiles.py',
        },
    ]
STAGE_KEYS = [stage['key'] for stage in STAGE_DEFS]


#%% main
def main(argv=None):
    parser = argparse.ArgumentParser(description='run core analysis pipeline')
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
        choices=['all', 'behaviour', 'lc', 'hpc', 'imaging'],
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
    args = parser.parse_args(argv)
    selected_defs = select_stage_defs(STAGE_DEFS, args, parser)
    stages = []
    for stage_def in selected_defs:
        script_args = []
        if stage_def.get('profile_stage') and args.n_workers is not None:
            script_args.extend(['--n-workers', str(args.n_workers)])

        stages.append({
            'key': stage_def['key'],
            'label': stage_def['label'],
            'command': python_command(repo_root, stage_def['script'], script_args),
            'display': display_command(stage_def['script'], script_args),
        })

    run_pipeline(
        stages,
        repo_root,
        'Core analysis pipeline',
        verbose=args.verbose,
        dry_run=args.dry_run,
        )

if __name__ == '__main__':
    main(sys.argv[1:])
