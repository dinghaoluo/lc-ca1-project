# -*- coding: utf-8 -*-
"""
Created on 4 June 16:30:51 2024

**NOTE: THIS SCRIPT RUNS THE CUSTOMISED SUITE2P-WANG-LAB INSTEAD OF SUITE2P**
**SUITE2P-WANG-LAB CAN BE ACCESSED HERE: 
    https://github.com/the-wang-lab/suite2p-wang-lab**

@author: Dinghao Luo
    - merged registration into script as a function
    - completely revamped main block to read in rec_list
        - if needs be, can import a separate list of parameters fit for each 
          recording session
    - changed the main run function 

"""

#%% imports 
import sys
import os
from pathlib import Path

repo_root = Path(__file__).resolve().parent
while not (repo_root / 'utils').exists() and repo_root != repo_root.parent:
    repo_root = repo_root.parent
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))
import imaging_utility_functions as iuf

external_code_root = Path(r'Z:\Dinghao\code_dinghao')
if str(external_code_root) not in sys.path:
    sys.path.append(str(external_code_root))
import rec_list
pathGRABNE = rec_list.pathHPCGRABNE


#%% run all sessions
def main():
    for path in pathGRABNE:
        sessname = path[-17:]
        print('\n{}'.format(sessname))
    
        reg_path = path+r'\processed'
        if not os.path.exists(reg_path):  # if registration has not been performed
            iuf.run_suite2p_registration(path)
        else:
            print('session already registered')
        # if 'no tiffs' is raised, most likely it is due to typos in pathnames
    
        roi_path = reg_path+r'\suite2p\plane0\stat.npy'
        if not os.path.exists(roi_path):  # if roi extraction has not been performed 
            iuf.run_suite2p_roi_extraction(path)
        else:
            print('session ROIs already extracted')


if __name__ == '__main__':
    main()
