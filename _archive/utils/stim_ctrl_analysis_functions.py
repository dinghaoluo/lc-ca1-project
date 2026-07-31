# -*- coding: utf-8 -*-
'''
Created on Sat Apr 11 2026

functions for HPC stim-control analyses.

@author: Dinghao Luo
'''

#%% imports
import numpy as np
from scipy.stats import sem, ranksums, ttest_ind, wilcoxon, ttest_rel

from console_formatting import (
    print_binwise_header,
    print_binwise_row,
)


#%% functions
def correct_teensy_overflow(times, of_constant):
    '''
    correct a single overflow event in monotonic teensy-style timestamps.
    '''
    times = np.array(times, dtype=float)

    correction_flag = 0
    correction_time = None
    corrected = [times[0]]
    for i in range(1, len(times)):
        if correction_flag == 1:
            corrected.append(times[i] + of_constant)
        elif times[i] >= times[i - 1] and correction_flag == 0:
            corrected.append(times[i])
        else:
            correction_flag = 1
            correction_time = i
            corrected.append(times[i] + of_constant)

    return corrected, correction_time


def annotate_binwise_pvals(ax, pval_tuple, bin_edges, base_y=4.05, dy=0.12, star=False):
    '''
    annotate a four-test ctrl-vs-stim p-value table above a binwise profile plot.
    '''
    p_ranksums, p_ttest_ind, p_wilcoxon, p_ttest_rel = pval_tuple
    all_tests = [
        ('RS', p_ranksums),
        ('tt_ind', p_ttest_ind),
        ('Wil', p_wilcoxon),
        ('tt_rel', p_ttest_rel),
    ]

    mids = (bin_edges[:-1] + bin_edges[1:]) / 2
    mids = (mids - 3750) / 1250

    for row, (label, pvals) in enumerate(all_tests):
        y = base_y - row * dy
        for mid, p in zip(mids, pvals):
            if star:
                if p < 0.0001:
                    text = '****'
                elif p < 0.001:
                    text = '***'
                elif p < 0.01:
                    text = '**'
                elif p < 0.05:
                    text = '*'
                else:
                    text = 'n.s.'
            else:
                text = f'{p:.1e}'

            ax.text(mid, y, text, ha='center', va='bottom', fontsize=3, color='k')


def compute_binwise_ctrl_stim_tests(ctrl_traces, stim_traces, bin_edges, bin_labels, label):
    '''
    compute the four binwise ctrl-vs-stim tests used in the HPC stim-control plots.
    '''
    pvals_ranksums, pvals_ttest_ind = [], []
    pvals_wilcoxon, pvals_ttest_rel = [], []

    paired = len(ctrl_traces) == len(stim_traces)

    print_binwise_header(label, 'ctrl', 'stim', ['rs', 'tt_ind', 'wil', 'tt_rel'])

    for b in range(len(bin_edges) - 1):
        start, end = bin_edges[b], bin_edges[b + 1]

        ctrl_bin = np.array([np.mean(tr[start:end]) for tr in ctrl_traces], dtype=float)
        stim_bin = np.array([np.mean(tr[start:end]) for tr in stim_traces], dtype=float)

        c = ctrl_bin[~np.isnan(ctrl_bin)]
        s = stim_bin[~np.isnan(stim_bin)]

        c_mean = np.mean(c) if c.size else np.nan
        s_mean = np.mean(s) if s.size else np.nan
        c_sem = sem(c) if c.size > 1 else np.nan
        s_sem = sem(s) if s.size > 1 else np.nan

        _, p_rs = ranksums(c, s, nan_policy='omit')
        _, p_ti = ttest_ind(c, s, nan_policy='omit')

        if paired and len(c) == len(s):
            _, p_w = wilcoxon(c, s, nan_policy='omit')
            _, p_tr = ttest_rel(c, s, nan_policy='omit')
        else:
            p_w = np.nan
            p_tr = np.nan

        pvals_ranksums.append(p_rs)
        pvals_ttest_ind.append(p_ti)
        pvals_wilcoxon.append(p_w)
        pvals_ttest_rel.append(p_tr)

        bin_label = str(bin_labels[b]).strip()
        bin_label = bin_label.replace(' - ', ' to ')
        bin_label = bin_label.replace('\u2013', ' to ').replace('\u2014', ' to ')
        print_binwise_row(
            bin_label,
            c_mean,
            c_sem,
            s_mean,
            s_sem,
            [p_rs, p_ti, p_w, p_tr],
        )

    return (np.array(pvals_ranksums), np.array(pvals_ttest_ind), np.array(pvals_wilcoxon), np.array(pvals_ttest_rel))

