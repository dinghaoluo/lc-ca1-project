# -*- coding: utf-8 -*-
'''
Created on Tue Jun  2 2026

run dLight and nLight alignment, control, dispersion, release, and ROI-neuropil analyses

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
        'label': 'Imaging alignment heatmap',
        'branch': 'alignment',
        'script': Path('analysis') / 'imaging' / 'alignment_analysis' / 'analyse_alignment_with_heatmap_run_cue_rew_aligned.py',
        },
    {
        'key': 'dlight_expression',
        'label': 'dLight expression control',
        'branch': 'dlight',
        'script': Path('analysis') / 'imaging' / 'controls' / 'dlight_expression_control.py',
        },
    {
        'key': 'dlight_inhibition',
        'label': 'dLight inhibition control',
        'branch': 'dlight',
        'script': Path('analysis') / 'imaging' / 'dLight_inhibition' / 'hpc_dlight_lc_inh_stim_ctrl_run.py',
        },
    {
        'key': 'dlight_dispersion',
        'label': 'dLight spatial dispersion',
        'branch': 'dlight',
        'script': Path('analysis') / 'imaging' / 'dLight_stim_dispersion' / 'single_roi_binned_dilation_spatial_tau.py',
        },
    {
        'key': 'dlight_lc_opto_release',
        'label': 'dLight LC-opto release profiles',
        'branch': 'dlight',
        'script': Path('analysis') / 'imaging' / 'optogenetics' / 'dlight_lc_opto_release_stim_ctrl.py',
        },
    {
        'key': 'release_probability_prop',
        'label': 'dLight significant-release proportion',
        'branch': 'dlight',
        'script': Path('analysis') / 'imaging' / 'release_probability' / 'prop_signif_release_dlight_stim.py',
        },
    {
        'key': 'release_probability',
        'label': 'dLight release probability',
        'branch': 'dlight',
        'script': Path('analysis') / 'imaging' / 'release_probability' / 'release_probability_dlight_stim.py',
        },
    {
        'key': 'release_heterogeneity',
        'label': 'dLight release-site heterogeneity',
        'branch': 'dlight',
        'script': Path('analysis') / 'imaging' / 'release_probability' / 'release_site_heterogeneity_dlight_stim.py',
        },
    {
        'key': 'roi_neuropil_mean',
        'label': 'ROI-vs-neuropil mean',
        'branch': 'dlight',
        'script': Path('analysis') / 'imaging' / 'ROI_vs_neuropil' / 'roi_vs_neuropil_ri_mean.py',
        },
    {
        'key': 'roi_neuropil_time',
        'label': 'ROI-vs-neuropil over time',
        'branch': 'dlight',
        'script': Path('analysis') / 'imaging' / 'ROI_vs_neuropil' / 'roi_vs_neuropil_ri_over_time.py',
        },
    {
        'key': 'nlight_expression',
        'label': 'nLight expression control',
        'branch': 'nlight',
        'script': Path('analysis') / 'imaging' / 'controls' / 'nlight_expression_control.py',
        },
    {
        'key': 'nlight_dispersion',
        'label': 'nLight spatial dispersion',
        'branch': 'nlight',
        'script': Path('analysis') / 'imaging' / 'nLight_stim_dispersion' / 'single_roi_binned_dilation_spatial_tau.py',
        },
    {
        'key': 'nlight_lc_opto_release',
        'label': 'nLight LC-opto release profiles',
        'branch': 'nlight',
        'script': Path('analysis') / 'imaging' / 'optogenetics' / 'nlight_lc_opto_release_stim_ctrl.py',
        },
    {
        'key': 'nlight_release_probability_prop',
        'label': 'nLight significant-release proportion',
        'branch': 'nlight',
        'script': Path('analysis') / 'imaging' / 'release_probability' / 'prop_signif_release_nlight_stim.py',
        },
    {
        'key': 'nlight_release_probability',
        'label': 'nLight release probability',
        'branch': 'nlight',
        'script': Path('analysis') / 'imaging' / 'release_probability' / 'release_probability_nlight_stim.py',
        },
    {
        'key': 'nlight_release_heterogeneity',
        'label': 'nLight release-site heterogeneity',
        'branch': 'nlight',
        'script': Path('analysis') / 'imaging' / 'release_probability' / 'release_site_heterogeneity_nlight_stim.py',
        },
    {
        'key': 'nlight_roi_neuropil_mean',
        'label': 'nLight ROI-vs-neuropil mean',
        'branch': 'nlight',
        'script': Path('analysis') / 'imaging' / 'ROI_vs_neuropil' / 'roi_vs_neuropil_ri_mean_nlight.py',
        },
    {
        'key': 'nlight_roi_neuropil_time',
        'label': 'nLight ROI-vs-neuropil over time',
        'branch': 'nlight',
        'script': Path('analysis') / 'imaging' / 'ROI_vs_neuropil' / 'roi_vs_neuropil_ri_over_time_nlight.py',
        },
    ]
STAGE_KEYS = [stage['key'] for stage in STAGE_DEFS]


#%% main
def main(argv=None):
    parser = argparse.ArgumentParser(
        description='run imaging analysis pipeline'
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
        choices=['all', 'alignment', 'dlight', 'nlight'],
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
        'Imaging analysis pipeline',
        verbose=args.verbose,
        dry_run=args.dry_run,
        )

if __name__ == '__main__':
    main(sys.argv[1:])
