# -*- coding: utf-8 -*-
'''
Created on 26 May 2026
Modified on 22 June 2026

plot saved nLight LC-opto stimulation summaries, release maps, and
axon-reference figures

@author: Dinghao Luo
'''

#%% imports
import argparse
from pathlib import Path
import sys

repo_root = Path(__file__).resolve().parents[3]
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))

from common_functions import mpl_formatting
import opto_imaging_plotting as optoplot
import project_paths as pp
mpl_formatting()

import rec_list


#%% paths
default_paths = rec_list.pathnLightLCOpto

all_sess_stem = pp.HPC_NLIGHT_LC_OPTO_STEM / 'all_sessions'
all_sess_fig_stem = pp.HPC_NLIGHT_LC_OPTO_FIGURES_STEM / 'all_sessions'


#%% main
def main(argv=None):
    parser = argparse.ArgumentParser(
        description='plot nLight LC-opto figures from saved extraction outputs.'
        )
    parser.add_argument(
        '--recording',
        dest='recording_filter',
        help='optional substring filter for recname or recording path',
        )
    parser.add_argument(
        '--skip-stim-summary',
        action='store_true',
        help='skip aligned whole-field and violin figures',
        )
    parser.add_argument(
        '--skip-release-maps',
        action='store_true',
        help='skip pixel-wise release-map figures',
        )
    parser.add_argument(
        '--skip-axon-reference',
        action='store_true',
        help='skip optional 1100-nm axon-reference figure',
        )
    args = parser.parse_args(argv)

    selected_paths = default_paths
    if args.recording_filter:
        selected_paths = [
            path for path in default_paths
            if (args.recording_filter in path
                or args.recording_filter in Path(path).name)
            ]
        if not selected_paths:
            raise ValueError(f'no recordings matched filter: {args.recording_filter}')

    for path in selected_paths:
        recname = Path(path).name
        savepath = all_sess_stem / recname
        figurepath = all_sess_fig_stem / recname
        print(f'\nplotting {recname}')
        optoplot.plot_session_outputs(
            savepath,
            recname,
            release_mode='paired',
            plot_stim_summary=not args.skip_stim_summary,
            plot_release_maps=not args.skip_release_maps,
            plot_axon_reference=not args.skip_axon_reference,
            output_path=figurepath,
            )

if __name__ == '__main__':
    main()
