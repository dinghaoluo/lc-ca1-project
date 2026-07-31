# -*- coding: utf-8 -*-
'''
run the shared behaviour parser for LCHPCGCaMPImmobile recordings.

This archived wrapper selects experiment code 10 and writes to the same
LCHPCGCaMPImmobile behaviour output directory.

@author: Dinghao Luo
'''

#%% imports
import sys
from pathlib import Path

repo_root = Path(__file__).resolve().parent
while (
        not (repo_root / 'preprocessing' / 'behaviour' / 'process_behaviour.py').exists()
        and repo_root != repo_root.parent
        ):
    repo_root = repo_root.parent
behaviour_root = repo_root / 'preprocessing' / 'behaviour'
if str(behaviour_root) not in sys.path:
    sys.path.insert(0, str(behaviour_root))

from process_behaviour import main


#%% main
if __name__ == '__main__':
    main(['10'])
