# -*- coding: utf-8 -*-
'''
Created on Fri Jul 12 17:35:35 2024
Modified on Wed 30 Apr 18:02:32 2025

process and save behaviour files as session dictionaries

run directly and enter a dataset code or name; processed pickles are written to
`lc-ca1-project/data/behaviour/all_experiments/<dataset_folder>`, where the
folder name is lowercase. The same selection can be passed with `--dataset`,
and `--all` processes datasets in the current figure-sweep order.

@author: Dinghao Luo
'''

#%% imports
import argparse
import pickle
import sys
from datetime import timedelta
from pathlib import Path
from time import time

repo_root = Path(__file__).resolve().parents[2]
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))

# import pre-processing functions
import behaviour_functions as bf
import project_paths as pp


#%% recording list
import rec_list


#%% constants
DATA_ROOT = repo_root / 'data'
DEFAULT_OUTPUT_ROOT = DATA_ROOT / 'behaviour' / 'all_experiments'

DATASET_SPECS = {
    '1': {
        'name': 'LC',
        'menu_label': 'LC tag',
        'paths': rec_list.pathLC,
        'txt_mode': 'standard',
        'processor': bf.process_behavioural_data,
    },
    '2': {
        'name': 'HPCRaphi',
        'menu_label': 'HPC Raphael Heldman',
        'paths': rec_list.pathHPC_Raphi,
        'txt_mode': 'raphi',
        'processor': bf.process_behavioural_data,
    },
    '3': {
        'name': 'HPCLC',
        'menu_label': 'HPC + LC act.',
        'paths': rec_list.pathHPCLCopt,
        'txt_mode': 'standard',
        'processor': bf.process_behavioural_data,
    },
    '4': {
        'name': 'HPCLCterm',
        'menu_label': 'HPC + LC-terminal act.',
        'paths': rec_list.pathHPCLCtermopt,
        'txt_mode': 'standard',
        'processor': bf.process_behavioural_data,
    },
    '5': {
        'name': 'HPCGRABNE',
        'menu_label': 'HPC GRABNE',
        'paths': rec_list.pathHPCGRABNE,
        'txt_mode': 'imaging',
        'processor': bf.process_behavioural_data_imaging,
    },
    '6': {
        'name': 'LCHPCGCaMP',
        'menu_label': 'LC-HPC axon-GCaMP',
        'paths': rec_list.pathLCHPCGCaMP,
        'txt_mode': 'imaging',
        'processor': bf.process_behavioural_data_imaging,
    },
    '7': {
        'name': 'HPCdLightLCOpto',
        'menu_label': 'HPC dLight + LC act.',
        'paths': rec_list.pathdLightLCOpto,
        'txt_mode': 'imaging',
        'processor': bf.process_behavioural_data_imaging,
    },
    '8': {
        'name': 'HPCdLightLCOptoInh',
        'menu_label': 'HPC dLight + LC inh.',
        'paths': rec_list.pathdLightLCOptoInh,
        'txt_mode': 'imaging',
        'processor': bf.process_behavioural_data_imaging,
    },
    '9': {
        'name': 'HPCnLightLCOpto',
        'menu_label': 'HPC nLight + LC act.',
        'paths': rec_list.pathnLightLCOpto,
        'txt_mode': 'imaging',
        'processor': bf.process_behavioural_data_imaging,
    },
    '10': {
        'name': 'LCHPCGCaMPImmobile',
        'menu_label': 'LC-HPC axon-GCaMP (immobile)',
        'paths': rec_list.pathLCHPCGCaMPImmobile,
        'txt_mode': 'immobile',
        'processor': bf.process_behavioural_data_immobile_imaging,
    },
}

DATASET_ALIASES = {}
for code, spec in DATASET_SPECS.items():
    DATASET_ALIASES[''.join(
        char for char in spec['name'].lower() if char.isalnum()
        )] = code
    DATASET_ALIASES[''.join(
        char for char in spec['menu_label'].lower() if char.isalnum()
        )] = code

DATASET_ALIASES.update({
    'all': 'all',
    '11': 'all',
})

DATASET_PROMPT = '\n'.join(
    ['process which list? (enter a dataset code/name, or "all")'] +
    [f'- {code}. {spec["name"]}: {spec["menu_label"]}'
     for code, spec in DATASET_SPECS.items()] +
    ['- 11. process all']
)


#%% main
def main(argv=None):
    parser = argparse.ArgumentParser(
        description='process and save behavioural session dictionaries.'
    )
    parser.add_argument(
        'selection',
        nargs='?',
        help='dataset code (1-10) or dataset name; use "all" or 11 to process every supported dataset',
    )
    parser.add_argument(
        '--dataset',
        dest='dataset',
        help='dataset code or dataset name; accepted as an explicit alternative to the positional selection',
    )
    parser.add_argument(
        '--all',
        action='store_true',
        help='process every supported dataset without prompting',
    )
    parser.add_argument(
        '--output-root',
        default=str(DEFAULT_OUTPUT_ROOT),
        help='directory that contains the per-dataset output folders',
    )
    args = parser.parse_args(argv)

    selection = 'all' if args.all else (args.dataset or args.selection)
    if selection is None:
        selection = input(f'{DATASET_PROMPT}\n> ').strip()

    code = selection.strip()
    if code not in DATASET_SPECS:
        code = DATASET_ALIASES.get(
            ''.join(char for char in code.lower() if char.isalnum())
        )
    if code is None:
        raise ValueError(
            'not a valid input; enter 1-10, 11 for all, \'all\', or a supported dataset name'
        )

    output_root = Path(args.output_root)
    specs = DATASET_SPECS.values() if code == 'all' else [DATASET_SPECS[code]]

    for spec in specs:
        if code == 'all':
            print(f'processing {spec["name"]}...\n')

        output_folder = output_root / pp.BEHAVIOUR_EXPERIMENT_FOLDER_NAMES[
            spec['name']
        ]
        for pathname in spec['paths']:
            recname = Path(pathname).name
            print(f'\nprocessing {recname}...')
            output_path = output_folder / f'{recname}.pkl'

            if spec['txt_mode'] == 'standard':
                txt_path = (
                    pp.MICEEXP_ROOT / f'ANMD{recname[1:5]}'
                    / recname[-17:-3]
                    / recname[-17:]
                    / f'{recname[-17:]}T.txt'
                )
            elif spec['txt_mode'] in {'imaging', 'immobile'}:
                txt_path = (
                    pp.MICEEXP_ROOT / f'ANMD{recname[1:4]}'
                    / f'{recname[:4]}{recname[5:]}T.txt'
                )
            else:
                txt_path = (
                    pp.RAPHAEL_ROOT / f'ANM{recname[1:4]}'
                    / recname[:13]
                    / recname
                    / f'{recname}T.txt'
                )

            session_start = time()
            try:
                behavioural_data = spec['processor'](str(txt_path))
            except (ValueError, FileNotFoundError) as e:
                print(f'  skipped ({e})')
                continue
            print(
                'session finished '
                f'({timedelta(seconds=int(time() - session_start))})'
            )

            save_start = time()
            output_folder.mkdir(parents=True, exist_ok=True)
            with open(output_path, 'wb') as f:
                pickle.dump(behavioural_data, f)
            print(
                'session saved '
                f'({timedelta(seconds=int(time() - save_start))})'
            )

if __name__ == '__main__':
    main()
