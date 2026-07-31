# -*- coding: utf-8 -*-
'''
Created on Thu Jun  1 16:29:45 2023
Modified on Tue 24 Sept 2024

plot first-lick-ordered LC rasters without mean profiles

@author: Dinghao Luo
'''

#%% imports
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.io as sio
from pathlib import Path
import sys

# plotting parameters
import matplotlib
plt.rcParams['font.family'] = 'Arial'
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

repo_root = Path(__file__).resolve().parents[3]
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))
from console_formatting import print_session
from lc_archive_support import load_lc_session_dict
import project_paths as pp


#%% load data
cell_prop = pd.read_pickle(pp.LC_ALL_STEM / 'LC_all_single_cell_properties.pkl')
output_stem = pp.LC_ALL_FIGURES_STEM / 'single_cell_raster_by_first_licks' / 'raster_only'
output_stem.mkdir(parents=True, exist_ok=True)


#%% specify RO peaking putative Dbh cells
clu_list = list(cell_prop.index)

tag_list = []; put_list = []
tag_rop_list = []; put_rop_list = []
for clu in cell_prop.index:
    tg = cell_prop['tagged'][clu]
    pt = cell_prop['putative'][clu]
    rop = cell_prop['peakness'][clu]

    if tg:
        tag_list.append(clu)
        if rop:
            tag_rop_list.append(clu)
    if pt:
        put_list.append(clu)
        if rop:
            put_rop_list.append(clu)


#%% main
noStim = 'N'

for cluname in clu_list:
    print_session(cluname)
    raster = load_lc_session_dict(cluname[:17], 'all_rasters')[cluname]
    train = load_lc_session_dict(cluname[:17], 'all_trains')[cluname]

    recname = cluname[:17]
    rec_stem = pp.MICEEXP_ROOT / f'ANMD{cluname[1:5]}' / cluname[:14] / recname
    filename = rec_stem / f'{recname}_DataStructure_mazeSection1_TrialType1_alignRun_msess1.mat'
    alignRun = sio.loadmat(filename)

    licks = alignRun['trialsRun']['lickLfpInd'][0][0][0][1:]
    starts = alignRun['trialsRun']['startLfpInd'][0][0][0][1:]
    tot_trial = licks.shape[0]

    behParf = rec_stem / f'{recname}_DataStructure_mazeSection1_TrialType1_behPar_msess1.mat'
    behPar = sio.loadmat(behParf)
    stimOn = behPar['behPar']['stimOn'][0][0][0][1:]
    stimOn_ind = np.where(stimOn!=0)[0]+1
    bad_beh_ind = np.where(behPar['behPar'][0]['indTrBadBeh'][0]==1)[1]-1

    first_licks = []
    for trial in range(tot_trial):
        lk = [l for l in licks[trial] if l-starts[trial] > 1250]  # only if the animal does not lick in the first second (carry-over licks)
        if len(lk)==0:
            first_licks.append(10000)
        else:
            first_licks.extend(lk[0]-starts[trial])

    temp = list(np.arange(tot_trial))
    licks_ordered, temp_ordered = zip(*sorted(zip(first_licks, temp)))

    if noStim=='Y' or noStim=='y':
        temp_ordered = [t for t in temp_ordered if t not in stimOn_ind]
        tot_trial = len(temp_ordered)  # reset tot_trial if noStim

    suffix = ' '
    if cluname in tag_list: suffix = ' tagged Dbh+'
    if cluname in put_list: suffix = ' putative Dbh+'
    clutitle = cluname + suffix

    # plotting
    fig, ax = plt.subplots(figsize=(3,2.2))
    ax.set(xticks=[0, 2, 4], xlim=(-1, 6), xlabel='time (s)',
           yticks=[1, 50, 100], ylabel='trial # by first licks',
           title=clutitle)
    for p in ['top', 'right']:
        ax.spines[p].set_visible(False)

    line_counter = 0  # for counting scatter plot lines
    for trial in range(tot_trial):
        if temp_ordered[trial] in bad_beh_ind:
            continue
        curr_raster = raster[temp_ordered[trial]]
        curr_train = train[temp_ordered[trial]]

        curr_trial = np.where(raster[temp_ordered[trial]]==1)[0]
        curr_trial = [(s-3750)/1250 for s in curr_trial if s>2500]  # starts from -1 s

        c = 'grey'
        calpha = 0.7
        dotsize = 1
        if (noStim=='N' or noStim=='n') and stimOn[temp_ordered[trial]]==1:
            c = 'royalblue'
            calpha = 1.0
            dotsize = 2

        ax.scatter(curr_trial, [line_counter+1]*len(curr_trial),
                   color=c, alpha=calpha, s=dotsize)
        ax.plot([licks_ordered[trial]/1250, licks_ordered[trial]/1250],
                [line_counter, line_counter+1],
                linewidth=2, color='orchid')

        line_counter+=1

    plt.show()

    if cluname in tag_list:
        fig.savefig(output_stem / f'{cluname}_tagged.png',
                    dpi=300,
                    bbox_inches='tight')
        fig.savefig(output_stem / f'{cluname}_tagged.pdf',
                    bbox_inches='tight')
    elif cluname in put_list:
        fig.savefig(output_stem / f'{cluname}_putative.png',
                    dpi=300,
                    bbox_inches='tight')
        fig.savefig(output_stem / f'{cluname}_putative.pdf',
                    bbox_inches='tight')
    else:
        fig.savefig(output_stem / f'{cluname}.png',
                    dpi=300,
                    bbox_inches='tight')
        fig.savefig(output_stem / f'{cluname}.pdf',
                    bbox_inches='tight')
