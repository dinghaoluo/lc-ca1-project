# -*- coding: utf-8 -*-
'''
Created on Wed 27 Sept 14:44:27 2023
Modified on Fri 10 Nov
Modified on Fri 20 Dec 2024:
    - merged everything together (HPCLC, HPCLCterm, pyr, int, you name it!)
    - use the hpc_all_profiles.pkl dataframe for info; discard everything else

compare all HPC cell's spiking profile between baseline ctrl and stim

@author: Dinghao Luo
'''

#%% imports
import numpy as np
import argparse
import gc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import sys
from pathlib import Path

repo_root = Path(__file__).resolve().parents[2]
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))
from common_functions import mpl_formatting
from console_formatting import print_session
import project_paths as pp
mpl_formatting()

import rec_list
pathHPCLC = rec_list.pathHPCLCopt
pathHPCLCterm = rec_list.pathHPCLCtermopt
paths = pathHPCLC + pathHPCLCterm


#%% parameters
run_onset_bin = 3750  # in samples
samp_freq = 1250  # in Hz
time_bef = 1  # in seconds
time_aft = 4  # in seconds
xaxis = np.arange(-samp_freq*time_bef, samp_freq*time_aft) / samp_freq
prof_window = (run_onset_bin-samp_freq*time_bef, run_onset_bin+samp_freq*time_aft)
PYRUP_CLASS = 'run-onset ON'
PYRDOWN_CLASS = 'run-onset OFF'
HPC_STEM = pp.HPC_EPHYS_STEM
HPC_FIGURES_STEM = pp.HPC_EPHYS_FIGURES_STEM

#%% arguments
parser = argparse.ArgumentParser()
parser.add_argument(
    '--from-rec',
    default=None,
    help='start plotting from this recording'
    )
args = parser.parse_args()

#%% load dataframe
print('loading dataframe...')
cell_profiles = pd.read_pickle(HPC_STEM / 'hpc_all_profiles.pkl')
started = args.from_rec is None
if args.from_rec is not None:
    recnames = {cluname[:17] for cluname in cell_profiles.index}
    if args.from_rec not in recnames:
        raise ValueError(f'{args.from_rec} not found')


#%% main
recname = cell_profiles.index[0]
for cluname in cell_profiles.index:
    curr_recname = cluname[:17]
    if not started:
        if curr_recname != args.from_rec:
            continue
        started = True

    if curr_recname!=recname:
        recname = curr_recname
        print_session(recname)

    clu_n = cluname.split(' ')[1]  # for plot titles

    cell = cell_profiles.loc[cluname]

    cell_identity = cell['cell_identity']  # 'pyr' or 'int'
    response_class = cell['class_ctrl']
    if response_class == PYRUP_CLASS:
        rt = 'RO-PyrUp'
    elif response_class == PYRDOWN_CLASS:
        rt = 'RO-PyrDown'
    else:
        rt = 'unresp.'

    ctrl_mean = cell['prof_ctrl_mean']
    ctrl_sem = cell['prof_ctrl_sem']
    stim_mean = cell['prof_stim_mean']
    stim_sem = cell['prof_stim_sem']

    # plotting
    fig, ax = plt.subplots(figsize=(1.6,1.2))

    ctrlln, = ax.plot(
        xaxis,
        ctrl_mean[prof_window[0]:prof_window[1]],
        color='grey')
    ax.fill_between(
        xaxis,
        ctrl_mean[prof_window[0]:prof_window[1]]+ctrl_sem[prof_window[0]:prof_window[1]],
        ctrl_mean[prof_window[0]:prof_window[1]]-ctrl_sem[prof_window[0]:prof_window[1]],
        alpha=.25, color='grey', edgecolor='none')
    stimln, = ax.plot(
        xaxis,
        stim_mean[prof_window[0]:prof_window[1]],
        color='royalblue')
    ax.fill_between(
        xaxis,
        stim_mean[prof_window[0]:prof_window[1]]+stim_sem[prof_window[0]:prof_window[1]],
        stim_mean[prof_window[0]:prof_window[1]]-stim_sem[prof_window[0]:prof_window[1]],
        alpha=.25, color='royalblue', edgecolor='none')

    ax.legend(
        [ctrlln, stimln], ['ctrl.', 'stim.'],
        frameon=False, fontsize=6)
    for p in ['top', 'right']:
        ax.spines[p].set_visible(False)
    ax.set(title=f'{recname}\n{clu_n} {rt}',
           xlabel='time from run-onset (s)', ylabel='spike rate (Hz)',
           xlim=(-time_bef, time_aft), xticks=(0,2,4))
    ax.title.set_fontsize(10)

    # pyr or int folder?
    pyr_dir = HPC_FIGURES_STEM / 'all_sessions' / recname / 'profiles_ctrl_stim_pyr'
    int_dir = HPC_FIGURES_STEM / 'all_sessions' / recname / 'profiles_ctrl_stim_int'

    for ext in ['.png', '.pdf']:
        if cell_identity == 'pyr':
            filepath = pyr_dir / f'{cluname}{ext}'
            filepath.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(filepath, dpi=300, bbox_inches='tight')
            if cell['rectype'] == 'HPCLC':
                filepath = HPC_FIGURES_STEM / 'single_cell_ctrl_stim_profiles' / 'HPC_LC_pyr' / f'{cluname}{ext}'
                filepath.parent.mkdir(parents=True, exist_ok=True)
                fig.savefig(filepath, dpi=300, bbox_inches='tight')
            if cell['rectype'] == 'HPCLCterm':
                filepath = HPC_FIGURES_STEM / 'single_cell_ctrl_stim_profiles' / 'HPC_LCterm_pyr' / f'{cluname}{ext}'
                filepath.parent.mkdir(parents=True, exist_ok=True)
                fig.savefig(filepath, dpi=300, bbox_inches='tight')
        else:
            filepath = int_dir / f'{cluname}{ext}'
            filepath.parent.mkdir(parents=True, exist_ok=True)
            fig.savefig(filepath, dpi=300, bbox_inches='tight')
            if cell['rectype'] == 'HPCLC':
                filepath = HPC_FIGURES_STEM / 'single_cell_ctrl_stim_profiles' / 'HPC_LC_int' / f'{cluname}{ext}'
                filepath.parent.mkdir(parents=True, exist_ok=True)
                fig.savefig(filepath, dpi=300, bbox_inches='tight')
            if cell['rectype'] == 'HPCLCterm':
                filepath = HPC_FIGURES_STEM / 'single_cell_ctrl_stim_profiles' / 'HPC_LCterm_int' / f'{cluname}{ext}'
                filepath.parent.mkdir(parents=True, exist_ok=True)
                fig.savefig(filepath, dpi=300, bbox_inches='tight')

    fig.clear()
    plt.close(fig)
    gc.collect()
