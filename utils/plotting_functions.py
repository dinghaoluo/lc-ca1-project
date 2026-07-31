# -*- coding: utf-8 -*-
'''
Created on Wed Aug 14 18:22:04 2024

Modified on Fri Apr 10 2026

plotting functions to save us from the chaos of having to code out the plotting
    section again and again without salvation which has been extremely painful
    and i do not understand why i did not do this earlier

@author: Dinghao Luo
'''

#%% imports
import matplotlib.pyplot as plt
import numpy as np
import sys
from scipy.stats import ranksums, sem, ttest_ind, ttest_rel, wilcoxon


#%% plotting helpers
def add_scale_bar(ax, x_start, y_start, x_len, y_len, color='k', lw=2):
    '''
    add horizontal and vertical scale bars to an existing axis.
    '''
    # horizontal scale bar (time)
    ax.plot([x_start, x_start + x_len], [y_start, y_start], color=color, lw=lw, solid_capstyle='butt')
    # vertical scale bar (dF/F)
    ax.plot([x_start, x_start], [y_start, y_start + y_len], color=color, lw=lw, solid_capstyle='butt')

def plot_bar_with_paired_scatter(
        ax, ctrl_vals, stim_vals, colors=('grey', 'firebrick'),
        title='', ylabel='% cells', xticklabels=('ctrl.', 'stim.'),
        ylim=None
        ):
    '''
    plot paired bar means with sem, individual paired scatter points,
    and paired statistical comparisons.

    parameters:
    - ax: matplotlib axis
        axis to draw the plot on
    - ctrl_vals: array-like
        control condition values (paired with stim_vals)
    - stim_vals: array-like
        stimulation condition values
    - colors: tuple
        colours for control and stimulation bars
    - title: str
        plot title
    - ylabel: str
        y-axis label
    - xticklabels: tuple
        x-axis category labels
    - ylim: tuple or none
        y-axis limits
    '''

    def annotate(axis, x1, x2, y, text, yrange):
        bump = 0.01 * yrange
        axis.plot([x1, x1, x2, x2], [y - 0.5 * bump, y, y, y - 0.5 * bump], lw=0.8, color='k')
        axis.text((x1 + x2) / 2, y + 0.5 * bump, text, ha='center', va='bottom', fontsize=6)

    ctrl = np.asarray(ctrl_vals)
    stim = np.asarray(stim_vals)

    # bars: mean +/- sem
    means = [np.nanmean(ctrl), np.nanmean(stim)]
    errs = [sem(ctrl, nan_policy='omit'), sem(stim, nan_policy='omit')]

    barc = ax.bar(
        [0, 1], means, yerr=errs, capsize=2, width=0.8, edgecolor='0.2', linewidth=1.0,
        alpha=.6, zorder=2,
        error_kw={'elinewidth': 0.6, 'capthick': 0.6, 'ecolor': 'k'}
    )
    barc.patches[0].set_facecolor(colors[0])
    barc.patches[1].set_facecolor(colors[1])

    # paired points + connecting lines (no jitter)
    for y0, y1 in zip(ctrl, stim):
        ax.plot([0, 1], [y0, y1], lw=0.8, color='k', alpha=.3, zorder=3)
    ax.scatter(np.zeros(len(ctrl)), ctrl, s=10, color=colors[0], edgecolor='none', alpha=.5, zorder=4)
    ax.scatter(np.ones(len(stim)), stim, s=10, color=colors[1], edgecolor='none', alpha=.5, zorder=4)

    # axes + limits
    ylims = (0, 100) if ylim is None else ylim
    ax.set(
        xticks=[0, 1],
        xticklabels=xticklabels,
        ylabel=ylabel,
        title=title,
        xlim=(-0.5, 1.5),
        ylim=ylims,
    )
    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)

    # stats (paired)
    w_stat, w_p = wilcoxon(ctrl, stim, alternative='two-sided', zero_method='wilcox', mode='auto')
    t_stat, t_p = ttest_rel(ctrl, stim)

    # annotations
    yrange = ylims[1] - ylims[0]
    top_data_val = max(np.nanmax(ctrl), np.nanmax(stim))
    top_bar_val = max(
        means[0] + errs[0],
        means[1] + errs[1],
    )
    top_y = max(top_data_val, top_bar_val)

    y1 = top_y
    y2 = y1 + 0.06 * yrange
    w_stars = 'ns' if w_p >= 0.05 else ('*' if w_p >= 0.01 else ('**' if w_p >= 0.001 else ('***' if w_p >= 1e-4 else '****')))
    t_stars = 'ns' if t_p >= 0.05 else ('*' if t_p >= 0.01 else ('**' if t_p >= 0.001 else ('***' if t_p >= 1e-4 else '****')))
    annotate(ax, 0, 1, y1, f'wilcoxon p={w_p:.4g} ({w_stars})', yrange)
    annotate(ax, 0, 1, y2, f't-test p={t_p:.4g} ({t_stars})', yrange)
    ax.set_ylim(ylims[0], max(ylims[1], y2 + 0.10 * yrange))

    print(f'wilcoxon: W={w_stat:.4g}, p={w_p:.2e}')
    print(f'ttest:    t={t_stat:.4g}, p={t_p:.2e}')

def plot_box_with_scatter(ctrl_data, stim_data, xlabel, savepath,
                          title='', show_scatter=True,
                          ctrl_color='grey', stim_color='royalblue'):
    '''
    plot horizontal boxplots with optional scatter overlay.

    parameters:
    - ctrl_data: array-like
        control condition values
    - stim_data: array-like
        stimulation condition values
    - xlabel: str
        x-axis label
    - savepath: str
        file path (without extension) for saving
    - title: str
        plot title
    - show_scatter: bool
        whether to overlay individual data points
    - ctrl_color: str or tuple
        colour for control box
    - stim_color: str or tuple
        colour for stimulation box
    '''
    ctrl_data = np.asarray(ctrl_data)
    stim_data = np.asarray(stim_data)

    fig, ax = plt.subplots(figsize=(2.6, 1.4))

    boxplots = ax.boxplot(
        [stim_data, ctrl_data],
        vert=False,
        positions=[2, 1],
        widths=0.25,
        patch_artist=True,
        medianprops={'color': 'black'},
        capprops={'color': 'black'},
        whiskerprops={'color': 'black'},
        flierprops={'marker': 'o', 'color': 'black', 'markersize': 3},
    )

    for patch, color in zip(boxplots['boxes'], (stim_color, ctrl_color)):
        patch.set_facecolor(color)
        patch.set_alpha(1)

    if show_scatter:
        ax.scatter(stim_data, np.full(stim_data.shape, 1.7), s=20, c=stim_color, alpha=0.75, ec='none')
        ax.scatter(ctrl_data, np.full(ctrl_data.shape, 1.3), s=20, c=ctrl_color, alpha=0.75, ec='none')

    for spine in ('top', 'right', 'left'):
        ax.spines[spine].set_visible(False)

    ax.set(
        title=title,
        xlabel=xlabel,
        yticks=(1, 2),
        yticklabels=('ctrl.', 'stim.'),
        ylim=(0.6, 2.4),
    )

    for ext in ('.png', '.pdf'):
        fig.savefig(f'{savepath}{ext}', dpi=300, bbox_inches='tight')

    plt.close(fig)

def plot_violin_with_scatter(data0, data1, colour0, colour1,
                             paired=True, alpha=.25,
                             xticklabels=('data0', 'data1'), ylabel=' ', ylim=None,
                             title=' ', stats_labels=None,
                             showscatter=False, showmainline=True,
                             print_statistics=False, plot_statistics=True,
                             save=False, savepath=' ', dpi=300,
                             figsize=(1.8, 2.2), show=True, close=False):
    '''
    plot half-violins with optional scatter and statistical comparisons.

    parameters
    ----------
    data0 : array-like
        values for the first dataset (plotted at x=1)
    data1 : array-like
        values for the second dataset (plotted at x=2)
    colour0 : str or tuple
        colour for the first dataset
    colour1 : str or tuple
        colour for the second dataset
    paired : bool, optional
        if true, use paired statistical tests; otherwise, use unpaired tests (default: True)
    alpha : float, optional
        transparency level for scatter points and lines (default: 0.25)
    xticklabels : tuple, optional
        labels for the x-axis ticks (default: ('data0', 'data1'))
    ylabel : str, optional
        label for the y-axis (default: ' ')
    ylim : tuple, optional
        limits for the y-axis (default: None)
    title : str, optional
        title for the plot (default: ' ')
    showscatter : bool, optional
        if true, scatter individual data points on the plot (default: False)
    showmainline : bool, optional
        if true, draw a line connecting median values (default: True)
    print_statistics : bool, optional
        if true, print statistical results in the console (default: False)
    plot_statistics : bool, optional
        if true, display statistical results on the plot (default: True)
    save : bool, optional
        if true, save the plot as a .png and .pdf file (default: False)
    savepath : str, optional
        path to save the plot (default: ' ')
    dpi : int, optional
        resolution for the saved image (default: 300)
    figsize : tuple, optional
        figure size (default: (1.8, 2.2))
    show : bool, optional
        if true, display the figure interactively (default: True)
    close : bool, optional
        if true, close the figure before returning (default: False)

    returns
    -------
    none
    '''
    data0 = np.asarray(data0)
    data1 = np.asarray(data1)

    stats_labels = stats_labels if stats_labels is not None else xticklabels
    label0 = str(stats_labels[0])
    label1 = str(stats_labels[1])

    fig, ax = plt.subplots(figsize=figsize)

    vp = ax.violinplot(
        [data0, data1],
        positions=[1.1, 1.9],
        showmedians=True,
        showextrema=False,
    )

    vp['bodies'][0].set_color(colour0)
    vp['bodies'][1].set_color(colour1)

    vp['cmedians'].set_color('k')
    vp['cmedians'].set_linewidth(2)
    ax.scatter(1.25, np.median(data0), s=30, color=colour0, ec='none', lw=.5, zorder=2)
    ax.scatter(1.75, np.median(data1), s=30, color=colour1, ec='none', lw=.5, zorder=2)
    if paired:
        for y0, y1 in zip(data0, data1):
            ax.plot([1.25, 1.75], [y0, y1], color='grey', alpha=alpha, linewidth=1, zorder=1)
    if showmainline:
        ax.plot([1.25, 1.75], [np.median(data0), np.median(data1)], color='k', linewidth=2, zorder=1)

    for i in [0, 1]:
        vp['bodies'][i].set_edgecolor('none')
        vp['bodies'][i].set_alpha(.75)
        body = vp['bodies'][i]
        mid_x = np.mean(body.get_paths()[0].vertices[:, 0])
        if i == 0:
            body.get_paths()[0].vertices[:, 0] = np.clip(body.get_paths()[0].vertices[:, 0], -np.inf, mid_x)
        else:
            body.get_paths()[0].vertices[:, 0] = np.clip(body.get_paths()[0].vertices[:, 0], mid_x, np.inf)

    if showscatter:
        ax.scatter([1.25] * len(data0), data0, s=10, color=colour0, ec='none', lw=.5, alpha=alpha)
        ax.scatter([1.75] * len(data1), data1, s=10, color=colour1, ec='none', lw=.5, alpha=alpha)

    # print mean +/- sem above violins
    mean0, mean1 = np.nanmean(data0), np.nanmean(data1)
    sem0, sem1 = sem(data0, nan_policy='omit'), sem(data1, nan_policy='omit')
    med0, med1 = np.nanmedian(data0), np.nanmedian(data1)
    q25_0, q75_0 = np.percentile(data0, [25, 75])
    q25_1, q75_1 = np.percentile(data1, [25, 75])

    y_max = max(np.nanmax(data0), np.nanmax(data1))
    y_min = min(np.nanmin(data0), np.nanmin(data1))
    y_span = y_max - y_min
    y_offset = (ylim[1] - ylim[0]) * 0.05 if ylim else 0.05 * y_span

    ax.text(
        1.1,
        y_max + y_offset,
        f'Med = {med0:.2f}\n'
        f'IQR = [{q25_0:.2f}, {q75_0:.2f}]\n'
        f'{mean0:.2f} +/- {sem0:.2f}',
        ha='center',
        va='bottom',
        fontsize=7,
        color=colour0,
    )
    ax.text(
        1.9,
        y_max + y_offset,
        f'Med = {med1:.2f}\n'
        f'IQR = [{q25_1:.2f}, {q75_1:.2f}]\n'
        f'{mean1:.2f} +/- {sem1:.2f}',
        ha='center',
        va='bottom',
        fontsize=7,
        color=colour1,
    )

    ax.set(xlim=(.5, 2.5))

    if ylim is not None:
        ax.set(ylim=ylim)
        y_top = ylim[1]
        y_bottom = ylim[0]
    else:
        y_top = max(np.max(data0), np.max(data1))
        y_bottom = min(np.min(data0), np.min(data1))
    y_range_tot = y_top - y_bottom

    if paired:
        wilc_stat, wilc_p = wilcoxon(data0, data1, alternative='two-sided', zero_method='wilcox', mode='auto')
        ttest_stat, ttest_p = ttest_rel(data0, data1)
        wilc_p_str = '{:.2e}'.format(wilc_p)
        ttest_p_str = '{:.2e}'.format(ttest_p)
        if print_statistics:
            print(f'\n{label0}:')
            print(f'  median = {med0:.4g}')
            print(f'  IQR    = [{q25_0:.4g}, {q75_0:.4g}]')
            print(f'  mean +/- SEM = {mean0:.4g} +/- {sem0:.4g}')
            print(f'{label1}:')
            print(f'  median = {med1:.4g}')
            print(f'  IQR    = [{q25_1:.4g}, {q75_1:.4g}]')
            print(f'  mean +/- SEM = {mean1:.4g} +/- {sem1:.4g}')
            print(f'wilcoxon: W={wilc_stat:.4g}, p={wilc_p_str}')
            print(f'ttest:    t={ttest_stat:.4g}, p={ttest_p_str}')
        if plot_statistics:
            stat_y = y_top + 0.05 * y_range_tot
            ax.plot([1.1, 1.9], [stat_y, stat_y], c='k', lw=.5)
            ax.text(
                1.5,
                stat_y,
                f'wilc_p={wilc_p_str}\nttest_p={ttest_p_str}',
                ha='center',
                va='bottom',
                color='k',
                fontsize=8,
            )
    else:
        wilc_stat, wilc_p = ranksums(data0, data1)
        ttest_stat, ttest_p = ttest_ind(data0, data1)
        wilc_p_str = '{:.2e}'.format(wilc_p)
        ttest_p_str = '{:.2e}'.format(ttest_p)
        if print_statistics:
            print(f'\n{label0}:')
            print(f'  median = {med0:.4g}')
            print(f'  IQR    = [{q25_0:.4g}, {q75_0:.4g}]')
            print(f'  mean +/- SEM = {mean0:.4g} +/- {sem0:.4g}')
            print(f'{label1}:')
            print(f'  median = {med1:.4g}')
            print(f'  IQR    = [{q25_1:.4g}, {q75_1:.4g}]')
            print(f'  mean +/- SEM = {mean1:.4g} +/- {sem1:.4g}')
            print(f'ranksums: z={wilc_stat:.4g}, p={wilc_p_str}')
            print(f'ttest:    t={ttest_stat:.4g}, p={ttest_p_str}')
        if plot_statistics:
            stat_y = y_top + 0.05 * y_range_tot
            ax.plot([1.1, 1.9], [stat_y, stat_y], c='k', lw=.5)
            ax.text(
                1.5,
                stat_y,
                f'ranksums_p={wilc_p_str}\nttest_p={ttest_p_str}',
                ha='center',
                va='bottom',
                color='k',
                fontsize=8,
            )

    ax.set(xticks=[1, 2], xticklabels=xticklabels, ylabel=ylabel, title=title)

    for spine in ['top', 'right', 'bottom']:
        ax.spines[spine].set_visible(False)

    fig.tight_layout()
    plt.grid(False)

    if save:
        for ext in ['.png', '.pdf']:
            fig.savefig(f'{savepath}{ext}', dpi=dpi, bbox_inches='tight')

    if show and not bool(getattr(sys.modules.get('__main__'), '__file__', None)):
        plt.show()

    if close:
        plt.close(fig)
