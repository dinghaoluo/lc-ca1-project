# -*- coding: utf-8 -*-
'''
Created on Tue Jun  2 2026

regenerate figures from saved behaviour, HPC, imaging, and LC outputs

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
        'key': 'behaviour_speed_licks',
        'label': 'Behaviour speed and lick plots',
        'branch': 'behaviour',
        'script': Path('plotting') / 'behaviour' / 'plot_speed_licks.py',
        },
    {
        'key': 'behaviour_pupil',
        'label': 'Behaviour pupil plots',
        'branch': 'behaviour',
        'script': Path('plotting') / 'behaviour' / 'plot_pupil_size.py',
        },
    {
        'key': 'behaviour_first_lick',
        'label': 'Behaviour first-lick plots',
        'branch': 'behaviour',
        'script': Path('plotting') / 'behaviour' / 'plot_first_lick_since_last_reward.py',
        },
    {
        'key': 'behaviour_example_session',
        'label': 'Behaviour example-session plots',
        'branch': 'behaviour',
        'script': Path('plotting') / 'behaviour' / 'plot_example_session.py',
        },
    {
        'key': 'behaviour_example_trials',
        'label': 'Behaviour example-trial plots',
        'branch': 'behaviour',
        'script': Path('plotting') / 'behaviour' / 'plot_example_trials.py',
        },
    {
        'key': 'behaviour_immobile',
        'label': 'Immobile behaviour plots',
        'branch': 'behaviour',
        'script': Path('plotting') / 'behaviour' / 'plot_immobile.py',
        },
    {
        'key': 'behaviour_speeds',
        'label': 'Behaviour speed plots',
        'branch': 'behaviour',
        'script': Path('plotting') / 'behaviour' / 'plot_speeds.py',
        },
    {
        'key': 'behaviour_trial_by_trial',
        'label': 'Behaviour trial-by-trial plots',
        'branch': 'behaviour',
        'script': Path('plotting') / 'behaviour' / 'plot_trial_by_trial.py',
        },
    {
        'key': 'hpc_ctrl_stim_profiles',
        'label': 'HPC ctrl/stim profile plots',
        'branch': 'hpc',
        'script': Path('plotting') / 'hpc' / 'plot_all_ctrl_stim_profiles.py',
        },
    {
        'key': 'hpc_ctrl_stim_rasters',
        'label': 'HPC ctrl/stim raster plots',
        'branch': 'hpc',
        'script': Path('plotting') / 'hpc' / 'plot_all_ctrl_stim_rasters.py',
        },
    {
        'key': 'hpc_heatmap_dist',
        'label': 'HPC distance heatmaps',
        'branch': 'hpc',
        'script': Path('plotting') / 'hpc' / 'plot_all_pyr_heatmap_dist.py',
        },
    {
        'key': 'hpc_info_ctrl_stim',
        'label': 'HPC information ctrl/stim plots',
        'branch': 'hpc',
        'script': Path('plotting') / 'hpc' / 'plot_all_pyr_info_ctrl_stim.py',
        },
    {
        'key': 'hpc_pre_post_ratio',
        'label': 'HPC pre/post ratio plots',
        'branch': 'hpc',
        'script': Path('plotting') / 'hpc' / 'plot_all_pyr_pre_post_ratio.py',
        },
    {
        'key': 'hpc_pre_post_raw_change',
        'label': 'HPC pre/post raw-change plots',
        'branch': 'hpc',
        'script': Path('plotting') / 'hpc' / 'plot_all_pyr_pre_post_raw_change.py',
        },
    {
        'key': 'hpc_run_onset_up_down',
        'label': 'HPC run-onset up/down plots',
        'branch': 'hpc',
        'script': Path('plotting') / 'hpc' / 'plot_run_onset_pyr_up_down_profiles.py',
        },
    {
        'key': 'hpc_run_onset_up_down_raphi',
        'label': 'Raphi HPC run-onset up/down plots',
        'branch': 'hpc',
        'script': Path('plotting') / 'hpc' / 'plot_run_onset_pyr_up_down_profiles_raphi.py',
        },
    {
        'key': 'hpc_sequence_place_time',
        'label': 'HPC place/time sequence plots',
        'branch': 'hpc',
        'script': Path('plotting') / 'hpc' / 'plot_hpc_sequence_place_time.py',
        },
    {
        'key': 'imaging_example_refs',
        'label': 'Imaging reference TIFF examples',
        'branch': 'imaging',
        'script': Path('plotting') / 'imaging' / 'example_sess_refs_release_tiff.py',
        },
    {
        'key': 'imaging_ref_maps',
        'label': 'Imaging reference maps',
        'branch': 'imaging',
        'script': Path('plotting') / 'imaging' / 'plot_ref_channel_16bit_maps.py',
        },
    {
        'key': 'imaging_dlight_single_axon',
        'label': 'dLight single-axon stim plots',
        'branch': 'imaging',
        'script': Path('plotting') / 'imaging' / 'plot_dlight_lc_opto_single_axon_stim_profiles.py',
        },
    {
        'key': 'imaging_nlight_single_axon',
        'label': 'nLight single-axon stim plots',
        'branch': 'imaging',
        'script': Path('plotting') / 'imaging' / 'plot_nlight_lc_opto_single_axon_stim_profiles.py',
        },
    {
        'key': 'imaging_lick_profile',
        'label': 'Imaging lick profile plots',
        'branch': 'imaging',
        'script': Path('plotting') / 'imaging' / 'plot_lick_profile.py',
        },
    {
        'key': 'imaging_lick_pumps',
        'label': 'Imaging lick-to-pump plots',
        'branch': 'imaging',
        'script': Path('plotting') / 'imaging' / 'plot_lick_profile_to_pumps.py',
        },
    {
        'key': 'imaging_pooled_heatmap',
        'label': 'Axon GCaMP pooled heatmaps',
        'branch': 'imaging',
        'script': Path('plotting') / 'imaging' / 'plot_pooled_heatmap_axon_gcamp.py',
        },
    {
        'key': 'imaging_raw_traces',
        'label': 'Axon GCaMP raw-trace plots',
        'branch': 'imaging',
        'script': Path('plotting') / 'imaging' / 'plot_raw_traces_axon_gcamp.py',
        },
    {
        'key': 'imaging_raw_trace_examples',
        'label': 'Axon GCaMP example-trial raw traces',
        'branch': 'imaging',
        'script': Path('plotting') / 'imaging' / 'plot_raw_traces_axon_gcamp_example_trials.py',
        },
    {
        'key': 'imaging_sorted_grids',
        'label': 'Imaging sorted grid heatmaps',
        'branch': 'imaging',
        'script': Path('plotting') / 'imaging' / 'plot_sorted_heatmaps_grids.py',
        },
    {
        'key': 'imaging_sorted_rois',
        'label': 'Imaging sorted ROI heatmaps',
        'branch': 'imaging',
        'script': Path('plotting') / 'imaging' / 'plot_sorted_heatmaps_rois.py',
        },
    {
        'key': 'imaging_std_heatmap',
        'label': 'Imaging standard-deviation heatmap',
        'branch': 'imaging',
        'script': Path('plotting') / 'imaging' / 'plot_std_heatmap.py',
        },
    {
        'key': 'imaging_whole_field',
        'label': 'Imaging whole-field plots',
        'branch': 'imaging',
        'script': Path('plotting') / 'imaging' / 'plot_whole_field.py',
        },
    {
        'key': 'imaging_dlight_summary',
        'label': 'dLight LC-opto summary plots',
        'branch': 'imaging',
        'script': Path('plotting') / 'imaging' / 'summarise_dlight_lc_opto_all.py',
        },
    {
        'key': 'imaging_dlight_ctrl_inh',
        'label': 'dLight control/inhibition summary plots',
        'branch': 'imaging',
        'script': Path('plotting') / 'imaging' / 'summarise_dlight_lc_opto_ctrl_inh.py',
        },
    {
        'key': 'imaging_nlight_summary',
        'label': 'nLight LC-opto summary plots',
        'branch': 'imaging',
        'script': Path('plotting') / 'imaging' / 'summarise_nlight_lc_opto_all.py',
        },
    {
        'key': 'lc_acgs_heatmap',
        'label': 'LC ACG and heatmap plots',
        'branch': 'lc',
        'script': Path('plotting') / 'lc' / 'plot_acgs_and_heatmap.py',
        },
    {
        'key': 'lc_tagged_vs_putative',
        'label': 'LC tagged-vs-putative plots',
        'branch': 'lc',
        'script': Path('plotting') / 'lc' / 'plot_lc_tagged_vs_putative.py',
        },
    {
        'key': 'lc_ctrl_stim_profiles',
        'label': 'LC ctrl/stim profile plots',
        'branch': 'lc',
        'script': Path('plotting') / 'lc' / 'plot_ctrl_stim_profiles.py',
        },
    {
        'key': 'lc_isis',
        'label': 'LC ISI plots',
        'branch': 'lc',
        'script': Path('plotting') / 'lc' / 'plot_isis.py',
        },
    {
        'key': 'lc_first_lick_profiles',
        'label': 'LC first-lick-sensitive profile plots',
        'branch': 'lc',
        'script': Path('plotting') / 'lc' / 'plot_lc_first_lick_sensitive_profiles.py',
        },
    {
        'key': 'lc_first_lick_rasters',
        'label': 'LC early/late first-lick rasters',
        'branch': 'lc',
        'script': Path('plotting') / 'lc' / 'plot_rasters_1st_lick_ordered_early_late_only.py',
        },
    {
        'key': 'lc_run_cue_reward_rasters',
        'label': 'LC run/cue/reward rasters',
        'branch': 'lc',
        'script': Path('plotting') / 'lc' / 'plot_rasters_run_cue_rew_aligned.py',
        },
    {
        'key': 'lc_run_onset_profiles',
        'label': 'LC run-onset burst profile plots',
        'branch': 'lc',
        'script': Path('plotting') / 'lc' / 'plot_runonset_burst_and_non_burst_profiles.py',
        },
    {
        'key': 'lc_single_cell_acg',
        'label': 'LC single-cell ACG plots',
        'branch': 'lc',
        'script': Path('plotting') / 'lc' / 'plot_single_cell_acg.py',
        },
    {
        'key': 'lc_single_cell_waveform',
        'label': 'LC single-cell waveform plots',
        'branch': 'lc',
        'script': Path('plotting') / 'lc' / 'plot_single_cell_waveform.py',
        },
    {
        'key': 'lc_tagged_good_bad',
        'label': 'LC tagged good/bad raster plots',
        'branch': 'lc',
        'script': Path('plotting') / 'lc' / 'plot_tagged_example_good_bad_raster.py',
        },
    {
        'key': 'lc_tagging_responses',
        'label': 'LC tagging-response plots',
        'branch': 'lc',
        'script': Path('plotting') / 'lc' / 'plot_tagging_responses.py',
        },
    {
        'key': 'lc_trial_profiles',
        'label': 'LC trial-profile plots',
        'branch': 'lc',
        'script': Path('plotting') / 'lc' / 'plot_lc_trial_profiles.py',
        },
    {
        'key': 'lc_raster_last_reward_current',
        'label': 'LC previous-to-current reward rasters',
        'branch': 'lc',
        'script': Path('plotting') / 'lc' / 'rasters' / 'all_raster_last_reward_to_current_trial.py',
        },
    {
        'key': 'lc_raster_last_reward',
        'label': 'LC last-reward rasters',
        'branch': 'lc',
        'script': Path('plotting') / 'lc' / 'rasters' / 'all_raster_last_rew_ordered.py',
        },
    {
        'key': 'lc_raster_lick_ordered',
        'label': 'LC lick-ordered rasters',
        'branch': 'lc',
        'script': Path('plotting') / 'lc' / 'rasters' / 'all_raster_lick_ordered.py',
        },
    {
        'key': 'lc_raster_lick_only',
        'label': 'LC lick-ordered raster-only plots',
        'branch': 'lc',
        'script': Path('plotting') / 'lc' / 'rasters' / 'all_raster_lick_ordered_raster_only.py',
        },
    {
        'key': 'lc_raster_lick_reward_sensitivity',
        'label': 'LC lick/reward-sensitivity rasters',
        'branch': 'lc',
        'script': Path('plotting') / 'lc' / 'rasters' / 'all_raster_lick_reward_sensitivity.py',
        },
    {
        'key': 'lc_raster_reward_ordered',
        'label': 'LC reward-ordered rasters',
        'branch': 'lc',
        'script': Path('plotting') / 'lc' / 'rasters' / 'all_raster_rew_ordered.py',
        },
    {
        'key': 'lc_raster_reward_to_run',
        'label': 'LC reward-to-run-onset rasters',
        'branch': 'lc',
        'script': Path('plotting') / 'lc' / 'rasters' / 'all_raster_reward_to_run_onset_ordered.py',
        },
    ]
STAGE_KEYS = [stage['key'] for stage in STAGE_DEFS]


#%% main
def main(argv=None):
    parser = argparse.ArgumentParser(
        description='run plotting pipeline'
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
        choices=['all', 'behaviour', 'hpc', 'imaging', 'lc'],
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
        'Plotting pipeline',
        verbose=args.verbose,
        dry_run=args.dry_run,
        )

if __name__ == '__main__':
    main(sys.argv[1:])
