# -*- coding: utf-8 -*-
"""
Created on Fri Aug 30 16:25:24 2024

PyrUp cell single-cell profile

@author: Dinghao Luo

"""

#%% imports
import pandas as pd
import numpy as np
import scipy.io as sio
from math import log
from scipy.stats import sem, poisson
import matplotlib.pyplot as plt
import sys
from pathlib import Path

# plotting parameters
import matplotlib
plt.rcParams['font.family'] = 'Arial'
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

repo_root = Path(__file__).resolve().parent
while not (repo_root / 'utils').exists() and repo_root != repo_root.parent:
    repo_root = repo_root.parent
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))
import plotting_functions as pf
import project_paths as pp

HPC_ALL_STEM = pp.HPC_ALL_STEM
HPC_ALL_FIGURES_STEM = pp.HPC_ALL_FIGURES_STEM


def save_figure(fig, filepath, **kwargs):
    filepath.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(filepath, **kwargs)


#%% run HPC-LC or HPC-LCterm
HPC_LC = 1

# load paths to recordings
import rec_list
if HPC_LC:
    pathHPC = rec_list.pathHPCLCopt
elif not HPC_LC:
    pathHPC = rec_list.pathHPCLCtermopt


#%% pre-post ratio dataframe
df = pd.read_pickle(HPC_ALL_STEM / 'HPC_LC_stim_stimcont_diff_profiles_pyr_only.pkl')


#%% main
for pathname in pathHPC:
    recname = pathname[-17:]
    if recname=='A063r-20230708-02' or recname=='A063r-20230708-01':  # lick detection problems
        continue
    print(recname)

    trains = list(np.load(HPC_ALL_STEM / recname / f'HPC_all_info_{recname}.npy',
                          allow_pickle=True).item().values())
    tot_trial = len(trains[0])

    rasters = list(np.load(HPC_ALL_STEM / f'HPC_all_rasters_npy_simp{recname}.npy',
                           allow_pickle=True).item().values())

    # determine if each cell is pyramidal or intern
    info = sio.loadmat('{}{}_DataStructure_mazeSection1_TrialType1_Info.mat'.format(pathname, recname))
    rec_info = info['rec'][0][0]
    spike_rate = rec_info['firingRate'][0]
    intern_id = rec_info['isIntern'][0]
    pyr_id = [not(clu) for clu in intern_id]
    tot_clu = len(pyr_id)
    tot_pyr = sum(pyr_id)

    # classify PyrUp/PyrDown neurones
    pyrup = []; pyrdown = []
    for cluname, row in df.iterrows():
        if cluname.split(' ')[0]==recname:
            clu_ID = int(cluname.split(' ')[1][3:])
            if row['ctrl_pre_post']>=1.25:
                pyrdown.append(clu_ID-2)
            if row['ctrl_pre_post']<=.8:
                pyrup.append(clu_ID-2)  # HPCLC_all_train.py adds 2 to the ID, so we subtracts 2 here

    for pyr in pyrup:
        curr_train = trains[pyr]
        temp = np.zeros((tot_trial, 8*1250))
        for trial in range(tot_trial):
            trial_length = len(curr_train[trial])
            if trial_length<8*1250 and trial_length>0:
                temp[trial, :trial_length] = curr_train[trial][:8*1250]  # use [-3, 5]
            else:
                temp[trial, :] = curr_train[trial][:8*1250]
        mean_prof = np.mean(temp, axis=0)*1250
        sem_prof = sem(temp, axis=0)*1250
        mean_sr = spike_rate[pyr]

        fig, ax = plt.subplots(figsize=(2,1.9)); xaxis=np.arange(1250*8)/1250-3
        lnprof, = ax.plot(xaxis, mean_prof, c='royalblue', lw=1)
        ax.fill_between(xaxis, mean_prof+sem_prof,
                               mean_prof-sem_prof,
                        alpha=.1, edgecolor='none', color='royalblue')
        ax.axhline(y=mean_sr, color='grey', alpha=.5, label='mean', linestyle='dashed')
        ax.legend(fontsize=8, frameon=False, loc='upper right')
        scale_min, scale_max = pf.scale_min_max(mean_prof[2500:8750], sem_prof[2500:8750])
        ax.set(xlabel='time (s)', xlim=(-1,4), xticks=[0,2,4],
               ylabel='spike rate (Hz)', ylim=(scale_min, scale_max),
               title='{} clu{}'.format(recname, pyr))
        for s in ['top', 'right']: ax.spines[s].set_visible(False)
        fig.tight_layout()
        plt.show(fig)
        save_figure(fig, HPC_ALL_FIGURES_STEM / f'HPC_pyrup_profiles{recname}_clu{pyr}.png',
                    dpi=300, bbox_inches='tight')
        plt.close(fig)

        # example Poisson probability plot if needed
        xaxis = np.arange(mean_sr*5)
        y = []
        for x in xaxis:
            y.append(poisson.pmf(x, mean_sr))
        fig, ax = plt.subplots(figsize=(2.8,2))
        ax.stem(xaxis,y,
                linefmt='black', markerfmt='k.', basefmt=' ')
        ax.set(xlabel='spikes',
               ylabel='probability',
               title=cluname)
        for s in ['top', 'right']: ax.spines[s].set_visible(False)
        fig.tight_layout()
        plt.show(fig)
        save_figure(fig, HPC_ALL_FIGURES_STEM / f'HPC_pyrup_profiles{recname}_clu{pyr}_poisson.png',
                    dpi=300, bbox_inches='tight')
        plt.close(fig)
        # example ends

        deviation = np.zeros((tot_trial, 16))
        for trial in range(tot_trial):
            curr_raster = rasters[pyr][trial]
            for tbin, t in enumerate(np.linspace(0,8,16+1)[:-1]):  # 3 seconds before, 5 seconds after
                curr_bin = sum(curr_raster[int(t*1250):int((t+1)*1250)])
                coeff=1
                if curr_bin<mean_sr: coeff=-1
                deviation[trial, tbin] = -log(poisson.pmf(curr_bin, mean_sr))*coeff  # polarity achieved with coeff

        mean_dev = np.mean(deviation, axis=0)
        sem_dev = sem(deviation, axis=0)
        fig, ax = plt.subplots(figsize=(2,1.9)); xaxis=np.linspace(-2.5,5.5,16+1)[:-1]
        ax.plot(xaxis, mean_dev, c='royalblue', lw=1)
        ax.fill_between(xaxis, mean_dev+sem_dev,
                               mean_dev-sem_dev,
                        color='royalblue', edgecolor='none', alpha=.1)
        scale_min, scale_max = pf.scale_min_max(mean_dev[3:14], sem_dev[3:14])
        ax.set(xlabel='time (s)', xticks=[0,2,4], xlim=(-1,4),
               ylabel='Poisson dev.', ylim=(scale_min, scale_max),
               title='{} clu{}'.format(recname, pyr))
        for s in ['top', 'right']: ax.spines[s].set_visible(False)
        fig.tight_layout()
        plt.show(fig)
        save_figure(fig, HPC_ALL_FIGURES_STEM / f'HPC_pyrup_profiles{recname}_clu{pyr}_Poisson_deviation.png',
                    dpi=300, bbox_inches='tight')
        plt.close(fig)
