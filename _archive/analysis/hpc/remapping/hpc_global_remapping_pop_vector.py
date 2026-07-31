# -*- coding: utf-8 -*-
"""
Created on Fri May  9 19:43:16 2025

checking global remapping between baseline and ctrl/stim

@author: Dinghao Luo
"""

#%% imports
import numpy as np
import pickle
import matplotlib.pyplot as plt
import sys
import pandas as pd
from pathlib import Path
from scipy.stats import pearsonr

repo_root = Path(__file__).resolve().parent
while not (repo_root / 'utils').exists() and repo_root != repo_root.parent:
    repo_root = repo_root.parent
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))
from common_functions import mpl_formatting, normalise_to_all
from plotting_functions import plot_violin_with_scatter
import project_paths as pp
mpl_formatting()

import rec_list
pathHPCLC = rec_list.pathHPCLCopt
pathHPCLCterm = rec_list.pathHPCLCtermopt
PYRUP_CLASS   = 'run-onset ON'
PYRDOWN_CLASS = 'run-onset OFF'
HPC_EPHYS_STEM = pp.HPC_EPHYS_STEM
BEHAVIOUR_EXPERIMENTS_STEM = pp.BEHAVIOUR_EXPERIMENTS_STEM


#%% load dataframe
print('loading dataframe...')
cell_profiles = pd.read_pickle(HPC_EPHYS_STEM / 'hpc_all_profiles.pkl')
df_pyr = cell_profiles[cell_profiles['cell_identity']=='pyr']  # pyramidal only

df_pyrup   = df_pyr[df_pyr['class'] == PYRUP_CLASS]
df_pyrdown = df_pyr[df_pyr['class'] == PYRDOWN_CLASS]


#%% main
r_A_B = []
r_AB_ctrl = []
r_AB_stim = []
r_ctrl_stim = []

r_A_B_pyrup = []
r_AB_ctrl_pyrup = []
r_AB_stim_pyrup = []
r_ctrl_stim_pyrup = []

r_A_B_pyrdown = []
r_AB_ctrl_pyrdown = []
r_AB_stim_pyrdown = []
r_ctrl_stim_pyrdown = []

for path in pathHPCLCterm:
    recname = path[-17:]
    print(f'n{recname}')

    trains = np.load(
        HPC_EPHYS_STEM / 'all_sessions' / recname / f'{recname}_all_trains.npy',
        allow_pickle=True
        ).item()

    hpc_lc_beh_path = BEHAVIOUR_EXPERIMENTS_STEM / 'HPCLC' / f'{recname}.pkl'
    if hpc_lc_beh_path.exists():
        with open(hpc_lc_beh_path, 'rb') as f:
            beh = pickle.load(f)
    else:
        hpc_lcterm_beh_path = BEHAVIOUR_EXPERIMENTS_STEM / 'HPCLCterm' / f'{recname}.pkl'
        with open(hpc_lcterm_beh_path, 'rb') as f:
            beh = pickle.load(f)

    stim_conds = [t[15] for t in beh['trial_statements']][1:]
    stim_idx = [trial for trial, cond in enumerate(stim_conds)
                if cond!='0']
    ctrl_idx = [trial+2 for trial in stim_idx]
    baseline_idx = list(np.arange(stim_idx[0]))

    curr_df_pyr = df_pyr[df_pyr['recname']==recname]
    curr_df_pyrup = curr_df_pyr[curr_df_pyr['class'] == PYRUP_CLASS]
    curr_df_pyrdown = curr_df_pyr[curr_df_pyr['class'] == PYRDOWN_CLASS]
    pyr_list = curr_df_pyr.index.tolist()
    pyrup_list = curr_df_pyrup.index.tolist()
    pyrdown_list = curr_df_pyrdown.index.tolist()

    r_A_B_curr = []
    r_A_B_curr_pyrup = []
    r_A_B_curr_pyrdown = []
    for i in range(5):  # hard-coded for now; bootstrapping for baseline half-halves
        np.random.shuffle(baseline_idx)

        # split into two equal halves
        midpoint = len(baseline_idx) // 2
        A_idx = baseline_idx[:midpoint]
        B_idx = baseline_idx[midpoint:]

        pop_vector_A = []
        pop_vector_B = []
        pop_vector_A_pyrup = []
        pop_vector_B_pyrup = []
        pop_vector_A_pyrdown = []
        pop_vector_B_pyrdown = []
        for cluname in pyr_list:  # accumulate pop vectors
            mean_A = list(np.mean(trains[cluname][A_idx, :], axis=0))
            mean_B = list(np.mean(trains[cluname][B_idx, :], axis=0))
            pop_vector_A.extend(mean_A)
            pop_vector_B.extend(mean_B)

            if cluname in pyrup_list:
                pop_vector_A_pyrup.extend(mean_A)
                pop_vector_B_pyrup.extend(mean_B)
            if cluname in pyrdown_list:
                pop_vector_A_pyrdown.extend(mean_A)
                pop_vector_B_pyrdown.extend(mean_B)

        r, p = pearsonr(pop_vector_A, pop_vector_B)  # half-half corr
        r_A_B_curr.append(r)
        if len(pyrup_list)>0:
            r, p = pearsonr(pop_vector_A_pyrup, pop_vector_B_pyrup)
            r_A_B_curr_pyrup.append(r)
        if len(pyrdown_list)>0:
            r, p = pearsonr(pop_vector_A_pyrdown, pop_vector_B_pyrdown)
            r_A_B_curr_pyrdown.append(r)

    r_A_B.append(np.mean(r_A_B_curr))  # mean over bootstrapped r's
    if len(pyrup_list)>0:
        r_A_B_pyrup.append(np.mean(r_A_B_curr_pyrup))
    if len(pyrdown_list)>0:
        r_A_B_pyrdown.append(np.mean(r_A_B_curr_pyrdown))

    # next, calculate AB ctrl and AB stim
    pop_vector_AB = []
    pop_vector_ctrl = []
    pop_vector_stim = []
    pop_vector_AB_pyrup = []
    pop_vector_ctrl_pyrup = []
    pop_vector_stim_pyrup = []
    pop_vector_AB_pyrdown = []
    pop_vector_ctrl_pyrdown = []
    pop_vector_stim_pyrdown = []
    for cluname in pyr_list:
        mean_AB = list(np.mean(trains[cluname][baseline_idx, :], axis=0))
        mean_ctrl = list(np.mean(trains[cluname][ctrl_idx, :], axis=0))
        mean_stim = list(np.mean(trains[cluname][stim_idx, :], axis=0))
        pop_vector_AB.extend(mean_AB)
        pop_vector_ctrl.extend(mean_ctrl)
        pop_vector_stim.extend(mean_stim)

        if cluname in pyrup_list:
            pop_vector_AB_pyrup.extend(mean_AB)
            pop_vector_ctrl_pyrup.extend(mean_ctrl)
            pop_vector_stim_pyrup.extend(mean_stim)
        if cluname in pyrdown_list:
            pop_vector_AB_pyrdown.extend(mean_AB)
            pop_vector_ctrl_pyrdown.extend(mean_ctrl)
            pop_vector_stim_pyrdown.extend(mean_stim)

    r, p = pearsonr(pop_vector_AB, pop_vector_ctrl)
    r_AB_ctrl.append(r)
    r, p = pearsonr(pop_vector_AB, pop_vector_stim)
    r_AB_stim.append(r)
    r, p = pearsonr(pop_vector_ctrl, pop_vector_stim)
    r_ctrl_stim.append(r)

    if len(pyrup_list)>0:
        r, p = pearsonr(pop_vector_AB_pyrup, pop_vector_ctrl_pyrup)
        r_AB_ctrl_pyrup.append(r)
        r, p = pearsonr(pop_vector_AB_pyrup, pop_vector_stim_pyrup)
        r_AB_stim_pyrup.append(r)
        r, p = pearsonr(pop_vector_ctrl_pyrup, pop_vector_stim_pyrup)
        r_ctrl_stim_pyrup.append(r)

    if len(pyrdown_list)>0:
        r, p = pearsonr(pop_vector_AB_pyrdown, pop_vector_ctrl_pyrdown)
        r_AB_ctrl_pyrdown.append(r)
        r, p = pearsonr(pop_vector_AB_pyrdown, pop_vector_stim_pyrdown)
        r_AB_stim_pyrdown.append(r)
        r, p = pearsonr(pop_vector_ctrl_pyrdown, pop_vector_stim_pyrdown)
        r_ctrl_stim_pyrdown.append(r)


#%% statistics
from scipy.stats import ranksums
# create figure with 3 subplots
fig, axes = plt.subplots(1, 3, figsize=(15, 5))
conditions = ['overall', 'PyrUp', 'PyrDown']

# prepare data lists
data_sets = [
    [r_A_B, r_AB_ctrl, r_AB_stim, r_ctrl_stim],
    [r_A_B_pyrup, r_AB_ctrl_pyrup, r_AB_stim_pyrup, r_ctrl_stim_pyrup],
    [r_A_B_pyrdown, r_AB_ctrl_pyrdown, r_AB_stim_pyrdown, r_ctrl_stim_pyrdown]
]

group_names = ['A vs B', 'AB ctrl', 'AB stim', 'ctrl vs stim']

# plot each subplot
for ax, condition, data in zip(axes, conditions, data_sets):
    ax.boxplot(data, labels=group_names)
    ax.set_title(f'{condition} comparisons')
    ax.set_ylabel('value')

    # perform ranksum tests and annotate
    for i in range(len(data)):
        for j in range(i + 1, len(data)):
            stat, p = ranksums(data[i], data[j])
            y_max = max(max(data[i]), max(data[j]))
            # ax.plot(np.arange(1,5), data, ls='dashed')
            ax.text((i + j) / 2 + 1, y_max * 1.05, f'p={p:.3f}',
                    ha='center', va='bottom', fontsize=8)

plt.tight_layout()
plt.show()
