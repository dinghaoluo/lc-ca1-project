# -*- coding: utf-8 -*-
'''
Created on Thu Jun  1 16:29:45 2023

Does the RO-peak have anything to do with licking?

@author: Dinghao Luo
'''

#%% imports
import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio
import pandas as pd
from scipy.stats import ttest_rel, ranksums, wilcoxon
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
output_no_stim = pp.LC_ALL_FIGURES_STEM / 'single_cell_raster_by_first_licks_noStim'
output_stim = pp.LC_ALL_FIGURES_STEM / 'single_cell_raster_by_first_licks'
output_no_stim.mkdir(parents=True, exist_ok=True)
output_stim.mkdir(parents=True, exist_ok=True)


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


#%% shuffle function
def cir_shuf(train, length, num_shuf=1000):
    tot_t = len(train)
    shuf_array = np.zeros([num_shuf, length])
    for i in range(num_shuf):
        rand_shift = np.random.randint(1, tot_t/2)
        shuf_array[i,:] = np.roll(train, -rand_shift)[:length]

    return np.mean(shuf_array, axis=0)


#%% main
noStim = 'N'

lick_sensitive = []
lick_sensitive_type = []

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

    # plotting
    fig, axs = plt.subplot_mosaic('AAAABBCC',figsize=(6.4,2.2))
    axs['A'].set(xlabel='time (s)', ylabel='trial # by first licks',
                 xlim=(-1, 7))
    for p in ['top', 'right']:
        axs['A'].spines[p].set_visible(False)

    pre_rate = []; post_rate = []
    pre_rate_shuf = []; post_rate_shuf = []
    ratio = []; ratio_shuf = []

    line_counter = 0  # for counting scatter plot lines
    for trial in range(tot_trial):
        if temp_ordered[trial] in bad_beh_ind:
            continue
        curr_raster = raster[temp_ordered[trial]]
        curr_train = train[temp_ordered[trial]]
        window = [licks_ordered[trial]+3750-625, licks_ordered[trial]+3750, licks_ordered[trial]+3750+625]
        pre_rate.append(sum(curr_train[window[0]:window[1]])*2)  # times 2 because it's half a second
        post_rate.append(sum(curr_train[window[1]:window[2]])*2)
        if sum(curr_train[window[1]:window[2]])*2!=0 and sum(curr_train[window[0]:window[1]])*2!=0:
            ratio.append(sum(curr_train[window[0]:window[1]])*2/sum(curr_train[window[1]:window[2]])*2)

        # shuffle
        length = len(curr_train)
        shuf_train = cir_shuf(curr_train, length, num_shuf=100)
        pre_rate_shuf.append(sum(shuf_train[window[0]:window[1]])*2)
        post_rate_shuf.append(sum(shuf_train[window[1]:window[2]])*2)
        if sum(shuf_train[window[1]:window[2]])*2!=0 and sum(shuf_train[window[0]:window[1]])*2!=0:
            ratio_shuf.append(sum(shuf_train[window[0]:window[1]])*2/sum(shuf_train[window[1]:window[2]])*2)

        curr_trial = np.where(raster[temp_ordered[trial]]==1)[0]
        curr_trial = [(s-3750)/1250 for s in curr_trial if s>2500]  # starts from -1 s

        c = 'grey'
        calpha = 0.7
        dotsize = 1
        if (noStim=='N' or noStim=='n') and stimOn[temp_ordered[trial]]==1:
            c = 'royalblue'
            calpha = 1.0
            dotsize = 2

        axs['A'].scatter(curr_trial, [line_counter+1]*len(curr_trial),
                         color=c, alpha=calpha, s=dotsize)
        axs['A'].plot([licks_ordered[trial]/1250, licks_ordered[trial]/1250],
                      [line_counter, line_counter+1],
                      linewidth=2, color='orchid')

        line_counter+=1

    axs['A'].set(yticks=[1, 50, 100], xticks=[0, 2, 4],
                 xlim=(-1, 6))

    # t-test and pre-post comp.
    t, tpval = ttest_rel(pre_rate, post_rate)
    w, wpval = wilcoxon(pre_rate, post_rate)
    t_ratio_res = ranksums(ratio, ratio_shuf)
    pval_ratio = t_ratio_res[1]
    if tpval<.01 and wpval<.01:
        lick_sensitive.append(True)
        if np.median(pre_rate)>np.median(post_rate):
            lick_sensitive_type.append('inhibition')
            axs['A'].set(title='{} inhibition'.format(cluname))
        else:
            lick_sensitive_type.append('excitation')
            axs['A'].set(title='{} excitation'.format(cluname))
    else:
        lick_sensitive.append(False)
        lick_sensitive_type.append('none')
        axs['A'].set(title=cluname)

    for p in ['top', 'right', 'bottom']:
        axs['B'].spines[p].set_visible(False)
    axs['B'].set_xticklabels(['pre', 'post'], minor=False)
    axs['B'].set(ylabel='spike rate (Hz)',
                 title='tpval={}\nwpval={}'.format(round(tpval,4), round(wpval, 4)))

    bp = axs['B'].boxplot([pre_rate, post_rate],
                          positions=[.5, 1],
                          patch_artist=True,
                          notch='True')
    colors = ['coral', 'darkcyan']
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
    bp['fliers'][0].set(marker ='v',
                    color ='#e7298a',
                    markersize=2,
                    alpha=0.5)
    bp['fliers'][1].set(marker ='o',
                    color ='#e7298a',
                    markersize=2,
                    alpha=0.5)
    for median in bp['medians']:
        median.set(color='darkred',
                   linewidth=1)

    for p in ['top', 'right', 'bottom']:
        axs['C'].spines[p].set_visible(False)
    axs['C'].set_xticklabels(['pre-\npost', 'pre-\npost-\nshuf'], minor=False)
    axs['C'].set(ylabel='pre-post ratio',
                 title='p[ratio]={}'.format(round(pval_ratio,4)))

    bp = axs['C'].boxplot([ratio, ratio_shuf],
                          positions=[.5, 1],
                          patch_artist=True,
                          notch='True')
    colors = ['royalblue', 'grey']
    for patch, color in zip(bp['boxes'], colors):
        patch.set_facecolor(color)
    bp['fliers'][0].set(marker ='v',
                    color ='#e7298a',
                    markersize=2,
                    alpha=0.5)
    bp['fliers'][1].set(marker ='o',
                    color ='#e7298a',
                    markersize=2,
                    alpha=0.5)
    for median in bp['medians']:
        median.set(color='darkred',
                   linewidth=1)

    plt.subplots_adjust(right=0.8)
    plt.grid(False)
    plt.show()

    if noStim=='Y' or noStim=='y':
        if cluname in tag_list:
            fig.savefig(output_no_stim / f'{cluname}_tagged.png',
                        dpi=300,
                        bbox_inches='tight')
            fig.savefig(output_no_stim / f'{cluname}_tagged.pdf',
                        bbox_inches='tight')
        elif cluname in put_list:
            fig.savefig(output_no_stim / f'{cluname}_putative.png',
                        dpi=300,
                        bbox_inches='tight')
            fig.savefig(output_no_stim / f'{cluname}_putative.pdf',
                        bbox_inches='tight')
        else:
            fig.savefig(output_no_stim / f'{cluname}.png',
                        dpi=300,
                        bbox_inches='tight')
            fig.savefig(output_no_stim / f'{cluname}.pdf',
                        bbox_inches='tight')
    else:
        if cluname in tag_list:
            fig.savefig(output_stim / f'{cluname}_tagged.png',
                        dpi=300,
                        bbox_inches='tight')
            fig.savefig(output_stim / f'{cluname}_tagged.pdf',
                        bbox_inches='tight')
        elif cluname in put_list:
            fig.savefig(output_stim / f'{cluname}_putative.png',
                        dpi=300,
                        bbox_inches='tight')
            fig.savefig(output_stim / f'{cluname}_putative.pdf',
                        bbox_inches='tight')
        else:
            fig.savefig(output_stim / f'{cluname}.png',
                        dpi=300,
                        bbox_inches='tight')
            fig.savefig(output_stim / f'{cluname}.pdf',
                        bbox_inches='tight')

    plt.close(fig)


#%% figure code, to plot density of stim trials
density = []
for trial in temp_ordered:
    if stimOn[trial]==1:
        density.append(1)
    else:
        density.append(0)
density_ind = np.where(np.array(density)==1)
density_ind = [0-s for s in density_ind]

fig, ax = plt.subplots(figsize=(10,2))

for p in ['top','right','left','bottom']:
    ax.spines[p].set_visible(False)
ax.set(yticks=[]); ax.set(xticks=[])

ax.hist(density_ind, bins=24, edgecolor='k', color='royalblue', linewidth=3)

import seaborn as sns
sns.set_style('whitegrid')
sns.kdeplot(density_ind, bw=0.5, ax=ax)

fig.tight_layout()
plt.show()

pp.LC_FIGURES_STEM.mkdir(parents=True, exist_ok=True)
fig.savefig(pp.LC_FIGURES_STEM / 'eg_session_stimdensity.png',
            dpi=500,
            bbox_inches='tight')


#%% save to dataframe
cell_prop = cell_prop.assign(lick_sensitive=pd.Series(lick_sensitive).values)
cell_prop = cell_prop.assign(lick_sensitive_type=pd.Series(lick_sensitive_type).values)

cell_prop.to_pickle(pp.LC_ALL_STEM / 'LC_all_single_cell_properties.pkl')
