'''
Created on 28 April 2026
Modified on 24 June 2026

plot PyrUp and PyrDown stimulation-control profiles and their session-wise
first-lick correlations, together with the acute mixed-effects estimate

@author: Dinghao Luo
'''

#%% imports
from pathlib import Path
import argparse
import pickle
import sys

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import linregress, sem

repo_root = Path(__file__).resolve().parents[3]
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))

from common_functions import mpl_formatting
from console_formatting import print_files_saved
mpl_formatting()


#%% parameters
SAVE_EXTENSIONS = ('.pdf', '.png')

CTRL_PYRUP_COLOR   = 'firebrick'
CTRL_PYRDOWN_COLOR = 'purple'
STIM_PYRUP_COLOR   = (87/255, 90/255, 187/255)
STIM_PYRDOWN_COLOR = (78/255, 84/255, 206/255)

default_input_path = (
    repo_root / 'data' / 'analysis' / 'hpc' / 'stim_ctrl'
    / 'lc_stim_ca1_pyr_up_down_behaviour_corr.pkl')
default_output_stem = (
    repo_root / 'figures' / 'fig_2_lc_stim_ca1_resp_and_beh_corr')


#%% profile and correlation plots
def _plot_group_profiles(
        ctrl_mean, ctrl_sem, stim_mean, stim_sem, pvals,
        title, ctrl_label, stim_label, ctrl_color, stim_color, filename,
        xaxis, bin_edges, output_stem):
    fig, ax = plt.subplots(figsize=(2.6, 2))
    ax.plot(xaxis, ctrl_mean, label=ctrl_label, color=ctrl_color)
    ax.fill_between(xaxis, ctrl_mean + ctrl_sem, ctrl_mean - ctrl_sem,
                    color=ctrl_color, alpha=.15)
    ax.plot(xaxis, stim_mean, label=stim_label, color=stim_color)
    ax.fill_between(xaxis, stim_mean + stim_sem, stim_mean - stim_sem,
                    color=stim_color, alpha=.15)

    # four ctrl-v-stim tests across the 1 s bins
    p_ranksums, p_ttest_ind, p_wilcoxon, p_ttest_rel = pvals
    all_tests = [
        ('RS', p_ranksums),
        ('tt_ind', p_ttest_ind),
        ('Wil', p_wilcoxon),
        ('tt_rel', p_ttest_rel),
    ]
    mids = (bin_edges[:-1] + bin_edges[1:]) / 2
    mids = (mids - 3750) / 1250

    for row, (_, pvals) in enumerate(all_tests):
        y = 4.05 - row * .12
        for mid, p in zip(mids, pvals):
            text = f'{p:.1e}'
            ax.text(mid, y, text, ha='center', va='bottom', fontsize=3, color='k')

    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)
    ax.set(xlabel='Time from run onset (s)', xticks=[0, 2, 4],
           ylabel='Firing rate (Hz)', yticks=[1, 2, 3, 4], ylim=(1, 4.1))
    ax.set_title(title, fontsize=10)
    ax.legend(fontsize=5, frameon=False)
    for ext in SAVE_EXTENSIONS:
        fig.savefig((output_stem / filename).with_suffix(ext), dpi=300, bbox_inches='tight')
    plt.close(fig)

def _plot_time_regression(records, amp_key, amp_bounds, time_bounds,
                          x_label, color, filename, output_stem):
    plot_records = []
    for record in records:
        amplitude = record[amp_key]
        lick_time = record['lick_time_delta']
        if amplitude is None or lick_time is None:
            continue
        if np.isnan(amplitude) or np.isnan(lick_time):
            continue
        amp_low, amp_high = amp_bounds
        time_low, time_high = time_bounds
        if amp_low is not None and amplitude < amp_low:
            continue
        if amp_high is not None and amplitude > amp_high:
            continue
        if time_low is not None and lick_time < time_low:
            continue
        if time_high is not None and lick_time > time_high:
            continue
        plot_records.append(record)

    x_values = [record[amp_key] for record in plot_records]
    y_values = [record['lick_time_delta'] for record in plot_records]
    slope, intercept, r, p, _ = linregress(x_values, y_values)

    fig, ax = plt.subplots(figsize=(2.4, 2.2))
    ax.scatter(x_values, y_values, color=color, s=30, alpha=0.8)
    fit_x = np.linspace(min(x_values), max(x_values), 100)
    ax.plot(fit_x, intercept + slope*fit_x, color='k', lw=1)
    ax.text(0.05, 0.95, f'$R = {r:.2f}$\n$p = {p:.3g}$\n$n = {len(x_values)}$',
            transform=ax.transAxes, ha='left', va='top', fontsize=9)
    ax.set_xlabel(x_label)
    ax.set_ylabel('Delta lick time (stim-ctrl)')
    for spine in ['top', 'right']:
        ax.spines[spine].set_visible(False)
    fig.tight_layout()
    for ext in SAVE_EXTENSIONS:
        fig.savefig((output_stem / filename).with_suffix(ext), dpi=300, bbox_inches='tight')
    plt.close(fig)

def main():
    parser = argparse.ArgumentParser(
        description='plot HPCLC pyramidal PyrUp/PyrDown behaviour correlations.'
    )
    parser.add_argument(
        '--input', type=Path, default=default_input_path,
        help='analysis data produced by analyse_lc_stim_ca1_pyr_up_down_behaviour_corr.py',
    )
    parser.add_argument(
        '--output-dir', type=Path, default=default_output_stem,
        help='directory for exported figure files',
    )
    args = parser.parse_args()
    with open(args.input, 'rb') as f:
        payload = pickle.load(f)

    output_stem = args.output_dir
    output_stem.mkdir(parents=True, exist_ok=True)

    metadata = payload['metadata']
    sample_frequency = metadata['sample_frequency_hz']
    xaxis = np.arange(-sample_frequency, 4*sample_frequency) / sample_frequency
    plot_window = slice(2*sample_frequency, 7*sample_frequency)
    bin_edges = metadata['bin_edges']
    exp = metadata['experiment']

    ctrl_pyrup = payload['profile_groups']['pyrup']['ctrl_profiles']
    stim_pyrup = payload['profile_groups']['pyrup']['stim_profiles']
    ctrl_pyrdown = payload['profile_groups']['pyrdown']['ctrl_profiles']
    stim_pyrdown = payload['profile_groups']['pyrdown']['stim_profiles']

    mean_ctrl_pyrup = np.mean(ctrl_pyrup, axis=0)[plot_window]
    sem_ctrl_pyrup = sem(ctrl_pyrup, axis=0)[plot_window]
    mean_stim_pyrup = np.mean(stim_pyrup, axis=0)[plot_window]
    sem_stim_pyrup = sem(stim_pyrup, axis=0)[plot_window]
    mean_ctrl_pyrdown = np.mean(ctrl_pyrdown, axis=0)[plot_window]
    sem_ctrl_pyrdown = sem(ctrl_pyrdown, axis=0)[plot_window]
    mean_stim_pyrdown = np.mean(stim_pyrdown, axis=0)[plot_window]
    sem_stim_pyrdown = sem(stim_pyrdown, axis=0)[plot_window]

    _plot_group_profiles(
        mean_ctrl_pyrup, sem_ctrl_pyrup,
        mean_stim_pyrup, sem_stim_pyrup,
        payload['binwise_tests']['pyrup'],
        f'{exp}\nPyrUp remainers',
        'ctrl. PyrUp', 'stim. PyrUp',
        CTRL_PYRUP_COLOR, STIM_PYRUP_COLOR,
        f'{exp}_PyrUp',
        xaxis, bin_edges, output_stem,
    )

    _plot_group_profiles(
        mean_ctrl_pyrdown, sem_ctrl_pyrdown,
        mean_stim_pyrdown, sem_stim_pyrdown,
        payload['binwise_tests']['pyrdown'],
        f'{exp}\nPyrDown remainers',
        'ctrl. PyrDown', 'stim. PyrDown',
        CTRL_PYRDOWN_COLOR, STIM_PYRDOWN_COLOR,
        f'{exp}_PyrDown',
        xaxis, bin_edges, output_stem,
    )

    session_records = payload['session_delta_records']
    _plot_time_regression(
        session_records,
        'amp_pyrup_delta_mean',
        metadata['pyrup_amp_bounds_hz'],
        (None, None),
        'Delta PyrUp (Hz)',
        CTRL_PYRUP_COLOR,
        f'{exp}_PyrUp_delta_amp_time',
        output_stem,
    )

    _plot_time_regression(
        session_records,
        'amp_pyrdown_delta_mean',
        metadata['pyrdown_amp_bounds_hz'],
        metadata['pyrdown_delta_bounds_ms'],
        'Delta PyrDown (Hz)',
        CTRL_PYRDOWN_COLOR,
        f'{exp}_PyrDown_delta_amp_time',
        output_stem,
    )

    mixed_effect = next(
        record for record in payload['mixed_model_fixed_effects']
        if record['term'] == 'is_stim'
    )
    estimate = mixed_effect['estimate']
    ci_low = mixed_effect['ci_low']
    ci_high = mixed_effect['ci_high']
    p_value = mixed_effect['p']
    n_cells = mixed_effect['n_cells']
    n_sessions = mixed_effect['n_sessions']

    fig, ax = plt.subplots(figsize=(2.5, 1.7))
    ax.axvline(0, color='grey', lw=.8, ls='--')
    ax.errorbar(
        estimate, 0,
        xerr=np.array([[estimate-ci_low], [ci_high-estimate]]),
        fmt='o', color=STIM_PYRUP_COLOR, capsize=3,
    )
    ax.set(
        xlabel='Stim. - flanking ctrl. (Hz)',
        yticks=[0], yticklabels=['all Pyr'], ylim=(-.55, .55),
    )
    ax.text(
        .98, .95,
        f'p = {p_value:.3g}\n'
        f'{n_cells} cells, {n_sessions} sessions',
        transform=ax.transAxes, ha='right', va='top', fontsize=7,
    )
    for spine in ['top', 'right', 'left']:
        ax.spines[spine].set_visible(False)
    fig.tight_layout()
    for ext in SAVE_EXTENSIONS:
        fig.savefig(
            (output_stem / f'{exp}_acute_firing_mixedlm').with_suffix(ext),
            dpi=300, bbox_inches='tight',
        )
    plt.close(fig)

    print_files_saved([
        ('figures', output_stem),
    ], gap=1)

if __name__ == '__main__':
    main()
