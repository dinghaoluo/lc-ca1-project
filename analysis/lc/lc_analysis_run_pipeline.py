# -*- coding: utf-8 -*-
'''
Created on Tue Jun  2 2026

run LC alignment, behaviour, burst, GLM, tonic, optogenetic, run-bout, and tagging analyses

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
        'key': 'alignment',
        'label': 'LC alignment heatmap',
        'branch': 'alignment',
        'script': Path('analysis') / 'lc' / 'alignment_analysis' / 'analyse_alignment_with_heatmap_run_cue_rew_aligned.py',
        },
    {
        'key': 'first_lick',
        'label': 'LC first-lick profiles',
        'branch': 'alignment',
        'script': Path('analysis') / 'lc' / 'first_lick_analysis' / 'all_earlyvlate_ro_peak_fixed_threshold.py',
        },
    {
        'key': 'opto_first_lick',
        'label': 'LC opto first lick',
        'branch': 'behaviour',
        'script': Path('analysis') / 'lc' / 'behaviour' / 'lc_opto_first_lick_profile.py',
        },
    {
        'key': 'opto_lick_properties',
        'label': 'LC opto lick properties',
        'branch': 'behaviour',
        'script': Path('analysis') / 'lc' / 'behaviour' / 'lc_opto_ctrl_vs_stim_lick_properties.py',
        },
    {
        'key': 'example_licks',
        'label': 'LC example licks',
        'branch': 'behaviour',
        'script': Path('analysis') / 'lc' / 'behaviour' / 'lc_example_session_licks.py',
        },
    {
        'key': 'example_licks_raphi',
        'label': 'Raphi example licks',
        'branch': 'behaviour',
        'script': Path('analysis') / 'lc' / 'behaviour' / 'lc_example_session_licks_passive_raphi.py',
        },
    {
        'key': 'example_speed',
        'label': 'LC example speed',
        'branch': 'behaviour',
        'script': Path('analysis') / 'lc' / 'behaviour' / 'lc_example_session_speed.py',
        },
    {
        'key': 'example_speed_raphi',
        'label': 'Raphi example speed',
        'branch': 'behaviour',
        'script': Path('analysis') / 'lc' / 'behaviour' / 'lc_example_session_speed_passive_raphi.py',
        },
    {
        'key': 'good_trial_percentage',
        'label': 'LC opto good-trial percentage',
        'branch': 'behaviour',
        'script': Path('analysis') / 'lc' / 'behaviour' / 'lc_opto_good_trial_percentage.py',
        },
    {
        'key': 'lick_distance',
        'label': 'LC opto lick distance',
        'branch': 'behaviour',
        'script': Path('analysis') / 'lc' / 'behaviour' / 'lc_opto_lick_distance.py',
        },
    {
        'key': 'hpc_lc_lick_distance',
        'label': 'HPC-LC stim lick distance',
        'branch': 'behaviour',
        'script': Path('analysis') / 'lc' / 'behaviour' / 'hpc_lc_stim_lick_distance.py',
        },
    {
        'key': 'lick_history',
        'label': 'LC opto lick history',
        'branch': 'behaviour',
        'script': Path('analysis') / 'lc' / 'behaviour' / 'lc_opto_lick_history.py',
        },
    {
        'key': 'lick_history_comp',
        'label': 'LC opto lick-history comparison',
        'branch': 'behaviour',
        'script': Path('analysis') / 'lc' / 'behaviour' / 'lc_opto_lick_history_comparison.py',
        },
    {
        'key': 'lick_time',
        'label': 'LC opto lick time',
        'branch': 'behaviour',
        'script': Path('analysis') / 'lc' / 'behaviour' / 'lc_opto_lick_time.py',
        },
    {
        'key': 'cue_start',
        'label': 'Cue-start difference',
        'branch': 'behaviour',
        'script': Path('analysis') / 'lc' / 'behaviour' / 'plot_cue_start_difference.py',
        },
    {
        'key': 'run_bouts',
        'label': 'LC run bouts',
        'branch': 'behaviour',
        'script': Path('analysis') / 'lc' / 'behaviour' / 'plot_run_bouts.py',
        },
    {
        'key': 'run_bout_examples',
        'label': 'LC run-bout examples',
        'branch': 'behaviour',
        'script': Path('analysis') / 'lc' / 'behaviour' / 'plot_run_bouts_examples.py',
        },
    {
        'key': 'single_trial_example',
        'label': 'LC single-trial example',
        'branch': 'behaviour',
        'script': Path('analysis') / 'lc' / 'behaviour' / 'plot_single_trial_example.py',
        },
    {
        'key': 'burst_baseline_rate',
        'label': 'Burst amplitude vs baseline rate',
        'branch': 'burst',
        'script': Path('analysis') / 'lc' / 'GLM' / 'burst_amplitude_vs_baseline_rate.py',
        },
    {
        'key': 'burst_reward_interval',
        'label': 'Burst amplitude vs reward interval',
        'branch': 'burst',
        'script': Path('analysis') / 'lc' / 'GLM' / 'burst_amplitude_vs_reward_interval.py',
        },
    {
        'key': 'burst_reward_interval_binned',
        'label': 'Binned burst amplitude vs reward interval',
        'branch': 'burst',
        'script': Path('analysis') / 'lc' / 'GLM' / 'burst_amplitude_vs_reward_interval_binned.py',
        },
    {
        'key': 'glm_permutation',
        'label': 'LC GLM permutation',
        'branch': 'glm',
        'script': Path('analysis') / 'lc' / 'GLM' / 'glm_lc_beh_permutation.py',
        },
    {
        'key': 'glm_permutation_full',
        'label': 'LC GLM permutation full',
        'branch': 'glm',
        'script': Path('analysis') / 'lc' / 'GLM' / 'glm_lc_beh_permutation_full.py',
        },
    {
        'key': 'tonic_fft',
        'label': 'LC tonic FFT',
        'branch': 'tonic',
        'script': Path('analysis') / 'lc' / 'GLM' / 'tonic_fft_lc.py',
        },
    {
        'key': 'stim_response',
        'label': 'LC stim response',
        'branch': 'opto',
        'script': Path('analysis') / 'lc' / 'ephys_opto' / 'analyse_stim_response.py',
        },
    {
        'key': 'run_onset_bout',
        'label': 'LC run-onset vs bout',
        'branch': 'run_bouts',
        'script': Path('analysis') / 'lc' / 'run_onset_v_run_bout' / 'lc_run_onset_vs_bout_burst_amplitude.py',
        },
    {
        'key': 'tagging_latency',
        'label': 'LC tagging latency',
        'branch': 'tagging',
        'script': Path('analysis') / 'lc' / 'tagging_analysis' / 'tagging_latency.py',
        },
    ]
STAGE_KEYS = [stage['key'] for stage in STAGE_DEFS]


#%% main
def main(argv=None):
    parser = argparse.ArgumentParser(description='run LC analysis pipeline')
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
        choices=[
            'all',
            'alignment',
            'behaviour',
            'burst',
            'glm',
            'tonic',
            'opto',
            'run_bouts',
            'tagging',
        ],
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
        'LC analysis pipeline',
        verbose=args.verbose,
        dry_run=args.dry_run,
        )

if __name__ == '__main__':
    main(sys.argv[1:])
