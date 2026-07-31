# -*- coding: utf-8 -*-
'''
Created on 4 June 16:30:51 2024
Modified on 8 July 2025
run batch registration with the Wang-lab Suite2p fork

fork: https://github.com/the-wang-lab/suite2p-wang-lab
accepts arbitrary directory lists for batch registration

@author: Dinghao Luo
'''

#%% imports
import argparse
import shutil
from pathlib import Path
import sys

repo_root = Path(__file__).resolve().parents[2]
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))
from console_formatting import print_session, print_status
import project_paths as pp
import imaging_utility_functions as iuf

import rec_list

paths = (
    rec_list.pathHPCGRABNE +
    rec_list.pathGRABNELCOpto +
    rec_list.pathGRABNELCOptoDbhBlock +
    rec_list.pathGRABNETone +
    rec_list.pathGRABNEToneDbhBlock +
    rec_list.pathdLightLCOpto +
    rec_list.pathdLightLCOptoDbhBlock +
    rec_list.pathdLightLCOptoCtrl +
    rec_list.pathdLightLCOptoInh +
    rec_list.pathnLightLCOpto
    )


#%% run all sessions
def main(argv=None):
    parser = argparse.ArgumentParser(
        description='run Wang-lab Suite2p registration over imaging recordings.'
        )
    parser.add_argument(
        '--recording',
        dest='recording_filter',
        help='optional substring filter for recname or recording path',
        )
    args = parser.parse_args(argv)

    selected_paths = paths
    if args.recording_filter:
        selected_paths = [
            path for path in paths
            if (args.recording_filter in path
                or args.recording_filter in Path(path).name)
            ]
        if not selected_paths:
            raise ValueError(f'no recordings matched filter: {args.recording_filter}')

    for path in selected_paths:
        recname = Path(path).name
        print_session(recname)

        raw_session_path = Path(path)
        if (
                (raw_session_path / 'suite2p' / 'plane0' / 'ops.npy').exists()
                or (
                    raw_session_path
                    / 'processed'
                    / 'suite2p'
                    / 'plane0'
                    / 'ops.npy'
                    ).exists()
                ):
            suite2p_stem = pp.resolve_suite2p_session_stem(path)
            mirror_stem = pp.SUITE2P_REGISTRATION_STEM / recname / 'suite2p'
            if (mirror_stem / 'plane0' / 'ops.npy').exists():
                print_status('skipped', f'session already registered at {suite2p_stem}')
                print_status('mirror present', mirror_stem)
                continue

            mirror_stem.parent.mkdir(parents=True, exist_ok=True)
            if mirror_stem.exists():
                shutil.rmtree(mirror_stem)
            shutil.copytree(suite2p_stem, mirror_stem)
            print_status('done', f'raw-side registration mirrored to {mirror_stem}')
            print_status('skipped', 'raw-side registration found')
            continue

        iuf.run_suite2p_registration(path)

if __name__ == '__main__':
    main()
