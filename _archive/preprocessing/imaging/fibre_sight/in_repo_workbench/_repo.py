'''
Created on 6 April 2026
Modified on 23 June 2026
paths used by the command-line tools and workbench

@author: Dinghao Luo
'''

#%% imports
from pathlib import Path
import sys


#%% repo paths
def get_repo_root():
    return Path(__file__).resolve().parents[3]


def add_repo_paths():
    repo_root = get_repo_root()
    utils_root = repo_root / 'utils'

    for path in [repo_root, utils_root]:
        path_str = str(path)
        if path_str not in sys.path:
            sys.path.insert(0, path_str)

    return repo_root


def default_source_root():
    add_repo_paths()
    import project_paths as pp

    return pp.HPC_DLIGHT_LC_OPTO_STEM / 'all_sessions'


def default_output_root():
    return get_repo_root() / 'data' / 'imaging' / 'fibre_sight'


def default_figure_root():
    return get_repo_root() / 'figures' / 'imaging' / 'fibre_sight'
