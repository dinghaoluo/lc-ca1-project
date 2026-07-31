'''
Created on 27 April 2026
Modified on 24 June 2026

compare PyrUp and PyrDown firing changes with stimulation-related changes in
first-lick time and distance across HPCLC sessions

add a first mixed-effects estimate using the two flanking controls

@author: Dinghao Luo
'''

#%% imports
from pathlib import Path
import argparse
import pickle
import sys

import numpy as np
import pandas as pd
from scipy.stats import linregress, ranksums, sem, ttest_ind, ttest_rel, wilcoxon
import statsmodels.formula.api as smf

repo_root = Path(__file__).resolve().parents[3]
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))

import project_paths as pp
import rec_list
from console_formatting import (
    print_binwise_header, print_binwise_row, print_files_saved, print_session,
    print_statistics_section,
)


#%% parameters
SAMP_FREQ       = 1250  # Hz
RUN_ONSET       = 3 * SAMP_FREQ
DELTA_THRES     = 0.5  # Hz
PRE_WINDOW      = slice(RUN_ONSET - SAMP_FREQ, RUN_ONSET)
# the scalar summary starts after the immediate stimulation transient
AMP_WINDOW      = slice(RUN_ONSET + int(1.5*SAMP_FREQ),
                        RUN_ONSET + int(2.5*SAMP_FREQ))

EVENT_BIN_MS   = 100
EVENT_BIN_SIZE = int(EVENT_BIN_MS / 1000 * SAMP_FREQ)
EVENT_START_S  = -0.5
EVENT_END_S    = 3.5
EVENT_START    = RUN_ONSET + int(EVENT_START_S * SAMP_FREQ)
EVENT_END      = RUN_ONSET + int(EVENT_END_S * SAMP_FREQ)
N_EVENT_BINS   = (EVENT_END - EVENT_START) // EVENT_BIN_SIZE

PYRUP_CLASS   = 'run-onset ON'
PYRDOWN_CLASS = 'run-onset OFF'

PYRUP_AMP_BOUNDS     = (-0.5, 2.0)          # Hz
PYRDOWN_AMP_BOUNDS   = (-2.0, DELTA_THRES)  # Hz
PYRDOWN_DELTA_BOUNDS = (-400.0, 1000.0)     # ms

bin_edges = np.arange(RUN_ONSET + int(-0.5*SAMP_FREQ),
                      RUN_ONSET + int(3.5*SAMP_FREQ) + 1,
                      SAMP_FREQ)
bin_labels = ['-0.5 - 0.5 s',
              ' 0.5 - 1.5 s',
              ' 1.5 - 2.5 s',
              ' 2.5 - 3.5 s']

all_sess_stem = pp.HPC_EPHYS_STEM / 'all_sessions'
HPCLC_BEH_STEM = pp.behaviour_experiment_stem('HPCLC')
data_output_stem = repo_root / 'data' / 'analysis' / 'hpc' / 'stim_ctrl'
default_output_path = (
    data_output_stem / 'lc_stim_ca1_pyr_up_down_behaviour_corr.pkl')
default_session_csv_path = (
    data_output_stem / 'lc_stim_ca1_pyr_up_down_behaviour_corr_sessions.csv')
cell_output_path = data_output_stem / 'lc_stim_ca1_acute_firing_cells.csv'
event_time_output_path = data_output_stem / 'lc_stim_ca1_event_time_effects.csv'
model_output_path = data_output_stem / 'lc_stim_ca1_acute_firing_mixedlm.csv'


#%% correlation and binwise statistics
def _fit_session_records(records, amp_key, y_key, amp_bounds, y_bounds):
    filtered = []
    for record in records:
        amplitude = record[amp_key]
        value = record[y_key]
        if np.isnan(amplitude) or np.isnan(value):
            continue
        amp_low, amp_high = amp_bounds
        value_low, value_high = y_bounds
        if amp_low is not None and amplitude < amp_low:
            continue
        if amp_high is not None and amplitude > amp_high:
            continue
        if value_low is not None and value < value_low:
            continue
        if value_high is not None and value > value_high:
            continue
        filtered.append(record)

    x_values = np.array([record[amp_key] for record in filtered], dtype=float)
    y_values = np.array([record[y_key] for record in filtered], dtype=float)
    slope, intercept, r, p, _ = linregress(x_values, y_values)
    return {
        'n': len(filtered),
        'slope': slope,
        'intercept': intercept,
        'r': r,
        'p': p,
    }

def compute_binwise_ctrl_stim_tests(ctrl_traces, stim_traces, bin_edges, bin_labels, label):
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
        bin_label = bin_label.replace('?', ' to ').replace('?', ' to ')
        print_binwise_row(
            bin_label,
            c_mean,
            c_sem,
            s_mean,
            s_sem,
            [p_rs, p_ti, p_w, p_tr],
        )

    return (
        np.array(pvals_ranksums), np.array(pvals_ttest_ind),
        np.array(pvals_wilcoxon), np.array(pvals_ttest_rel),
    )

def main(argv=None):
    parser = argparse.ArgumentParser(
        description='analyse HPCLC pyramidal PyrUp/PyrDown behaviour correlations.'
    )
    parser.add_argument(
        '--output', type=Path, default=default_output_path,
        help='pickle output path for the analysis data',
    )
    parser.add_argument(
        '--session-csv', type=Path, default=default_session_csv_path,
        help='CSV output path for per-session amplitude/behaviour deltas',
    )
    args = parser.parse_args(argv)

    print('loading dataframes...')
    cell_profiles = pd.read_pickle(pp.HPC_EPHYS_STEM / 'hpc_all_profiles.pkl')
    df_pyr = cell_profiles[cell_profiles['cell_identity'] == 'pyr']

    mean_prof_ctrl_pyrup, mean_prof_stim_pyrup     = [], []
    mean_prof_ctrl_pyrdown, mean_prof_stim_pyrdown = [], []
    session_delta_records = []
    cell_records = []
    event_time_records = []

    bad_behaviour_sessions = {
        Path(path).name for path in rec_list.pathHPCbadbeh
    }
    for path in rec_list.pathHPCLCopt:
        recname = Path(path).name
        if recname in bad_behaviour_sessions:
            continue
        print_session(recname)

        train_path = all_sess_stem / recname / f'{recname}_all_trains_run.npy'
        trains = np.load(train_path, allow_pickle=True).item()
        n_trials = next(iter(trains.values())).shape[0]

        beh_path = HPCLC_BEH_STEM / f'{recname}.pkl'
        with open(beh_path, 'rb') as f:
            beh = pickle.load(f)

        first_lick_times = []
        for licks, onset in zip(
                beh['lick_times'][1:n_trials + 1],
                beh['run_onsets'][1:n_trials + 1],
                ):
            eligible_licks = [
                lick[0] for lick in licks
                if np.isfinite(onset) and lick[0] > onset + 1000
            ]
            # keep empty trials in place or every later trial index shifts
            first_lick_times.append(
                eligible_licks[0] - onset if eligible_licks else np.nan
            )
        first_lick_times = np.asarray(first_lick_times)
        stim_conds = [
            beh['trial_statements'][trial + 1][15]
            for trial in range(n_trials)
        ]
        stim_idx = [trial for trial, condition in enumerate(stim_conds) if condition != '0']

        # the first stimulation has no preceding stimulated cycle
        stim_cycles = [
            trial for trial in stim_idx
            if (
                trial >= 3
                and trial + 2 < n_trials
                and stim_conds[trial - 3] != '0'
                and all(stim_conds[idx] == '0'
                        for idx in [trial - 2, trial - 1, trial + 1, trial + 2])
            )
        ]
        stim_cycles = np.asarray(stim_cycles, dtype=int)
        ctrl_before = stim_cycles - 1
        ctrl_after = stim_cycles + 2

        lick_effect = (
            first_lick_times[stim_cycles]
            - 2/3 * first_lick_times[ctrl_before]
            - 1/3 * first_lick_times[ctrl_after]
        )
        ctrl_stim_lick_time_delta = np.nanmedian(lick_effect)

        curr_df_pyr = df_pyr[df_pyr['recname'] == recname]
        sess_ctrl_pyrup, sess_stim_pyrup     = [], []
        sess_ctrl_pyrdown, sess_stim_pyrdown = [], []
        sess_effect_pyrup, sess_effect_pyrdown = [], []

        for idx, session in curr_df_pyr.iterrows():
            cluname = idx
            baseline_class = session['class']
            cell_trains = trains[cluname][:n_trials]

            pre_rate = np.mean(cell_trains[:, PRE_WINDOW], axis=1)
            post_rate = np.mean(cell_trains[:, AMP_WINDOW], axis=1)
            response_delta = post_rate - pre_rate

            stim_response = np.mean(response_delta[stim_cycles])
            ctrl_response = np.mean(
                2/3 * response_delta[ctrl_before]
                + 1/3 * response_delta[ctrl_after]
            )
            for condition, is_stim, response in [
                    ('ctrl', 0, ctrl_response),
                    ('stim', 1, stim_response),
                    ]:
                cell_records.append({
                    'animal': recname.split('-')[0],
                    'session': recname,
                    'cluname': cluname,
                    'cell_id': f'{recname}:{cluname}',
                    'baseline_class': baseline_class,
                    'condition': condition,
                    'is_stim': is_stim,
                    'response_delta_hz': response,
                    'n_cycles': len(stim_cycles),
                })

            mean_prof_ctrl = np.mean(
                2/3 * cell_trains[ctrl_before]
                + 1/3 * cell_trains[ctrl_after],
                axis=0,
            )
            mean_prof_stim = np.mean(cell_trains[stim_cycles], axis=0)
            stim_effect = stim_response - ctrl_response

            if baseline_class == PYRUP_CLASS:
                sess_ctrl_pyrup.append(mean_prof_ctrl)
                sess_stim_pyrup.append(mean_prof_stim)
                mean_prof_ctrl_pyrup.append(mean_prof_ctrl)
                mean_prof_stim_pyrup.append(mean_prof_stim)
                sess_effect_pyrup.append(stim_effect)

            if baseline_class == PYRDOWN_CLASS:
                sess_ctrl_pyrdown.append(mean_prof_ctrl)
                sess_stim_pyrdown.append(mean_prof_stim)
                mean_prof_ctrl_pyrdown.append(mean_prof_ctrl)
                mean_prof_stim_pyrdown.append(mean_prof_stim)
                sess_effect_pyrdown.append(stim_effect)

            event_rate = np.mean(
                cell_trains[:, EVENT_START:EVENT_END].reshape(
                    n_trials, N_EVENT_BINS, EVENT_BIN_SIZE
                ),
                axis=2,
            )
            cycle_effects = (
                event_rate[stim_cycles]
                - 2/3 * event_rate[ctrl_before]
                - 1/3 * event_rate[ctrl_after]
            )

            for bin_idx in range(N_EVENT_BINS):
                effect = cycle_effects[:, bin_idx]
                effect_mean = np.mean(effect)

                # one control reappears in the next contrast, and its covariance
                # enters twice in the variance of the mean
                residual = effect - effect_mean
                gamma_0 = np.dot(residual, residual) / len(effect)
                gamma_1 = np.dot(residual[1:], residual[:-1]) / len(effect)
                effect_se = np.sqrt(
                    max((gamma_0 + 2*gamma_1) / len(effect), 0)
                )
                time_start = EVENT_START_S + bin_idx * EVENT_BIN_MS / 1000
                event_time_records.append({
                    'animal': recname.split('-')[0],
                    'session': recname,
                    'cluname': cluname,
                    'baseline_class': baseline_class,
                    'time_start_s': time_start,
                    'time_centre_s': time_start + EVENT_BIN_MS / 2000,
                    'stim_effect_hz': effect_mean,
                    'stim_effect_se_hz': effect_se,
                    'n_cycles': len(stim_cycles),
                })

        amp_pyrup_delta_mean = np.mean(sess_effect_pyrup)
        amp_pyrdown_delta_mean = np.mean(sess_effect_pyrdown)

        session_delta_records.append({
            'recname': recname,
            'amp_pyrup_delta_mean': amp_pyrup_delta_mean,
            'amp_pyrdown_delta_mean': amp_pyrdown_delta_mean,
            'lick_time_delta': ctrl_stim_lick_time_delta,
            'n_pyrup': len(sess_ctrl_pyrup),
            'n_pyrdown': len(sess_ctrl_pyrdown),
            'n_stim_cycles': len(stim_cycles),
        })

    # the sandwich correction happens before the model; this first model only
    # asks whether the paired stimulation and control means differ
    cell_df = pd.DataFrame(cell_records)
    mixed_model = smf.mixedlm(
        'response_delta_hz ~ is_stim',
        cell_df,
        groups=cell_df['session'],
        re_formula='1',
        vc_formula={'cell': '0 + C(cell_id)'},
    )
    mixed_result = mixed_model.fit(reml=True, method='lbfgs')
    mixed_ci = mixed_result.conf_int()
    model_records = []
    for term in mixed_result.fe_params.index:
        model_records.append({
            'term': term,
            'estimate': mixed_result.fe_params[term],
            'standard_error': mixed_result.bse_fe[term],
            'z': mixed_result.tvalues[term],
            'p': mixed_result.pvalues[term],
            'ci_low': mixed_ci.loc[term, 0],
            'ci_high': mixed_ci.loc[term, 1],
            'converged': mixed_result.converged,
            'n_sessions': cell_df['session'].nunique(),
            'n_cells': cell_df['cell_id'].nunique(),
        })
    model_df = pd.DataFrame(model_records)

    print_statistics_section()
    stim_result = model_df[model_df['term'] == 'is_stim'].iloc[0]
    print(
        'mixed model stimulation effect: '
        f'{stim_result.estimate:.4g} Hz, '
        f'95% CI [{stim_result.ci_low:.4g}, '
        f'{stim_result.ci_high:.4g}], '
        f'p={stim_result.p:.3g}'
    )
    pvals_pyrup = compute_binwise_ctrl_stim_tests(
        mean_prof_ctrl_pyrup, mean_prof_stim_pyrup,
        bin_edges, bin_labels, 'PyrUp')
    pvals_pyrdown = compute_binwise_ctrl_stim_tests(
        mean_prof_ctrl_pyrdown, mean_prof_stim_pyrdown,
        bin_edges, bin_labels, 'PyrDown')

    regression_fits = {
        'pyrup_delta_amp_vs_lick_time': _fit_session_records(
            session_delta_records,
            'amp_pyrup_delta_mean',
            'lick_time_delta',
            PYRUP_AMP_BOUNDS,
            (None, None),
        ),
        'pyrdown_delta_amp_vs_lick_time': _fit_session_records(
            session_delta_records,
            'amp_pyrdown_delta_mean',
            'lick_time_delta',
            PYRDOWN_AMP_BOUNDS,
            PYRDOWN_DELTA_BOUNDS,
        ),
    }
    for label, fit in regression_fits.items():
        print(
            f'{label}: n={fit["n"]}, r={fit["r"]:.6f}, '
            f'p={fit["p"]:.6g}')

    payload = {
        'metadata': {
            'experiment': 'HPCLC',
            'sample_frequency_hz': SAMP_FREQ,
            'run_onset_sample': RUN_ONSET,
            'pre_window_samples': (PRE_WINDOW.start, PRE_WINDOW.stop),
            'amp_window_samples': (AMP_WINDOW.start, AMP_WINDOW.stop),
            'amp_window_s': (
                (AMP_WINDOW.start - RUN_ONSET) / SAMP_FREQ,
                (AMP_WINDOW.stop - RUN_ONSET) / SAMP_FREQ,
            ),
            'event_bin_ms': EVENT_BIN_MS,
            'event_window_s': (EVENT_START_S, EVENT_END_S),
            'control_definition': '2/3 stim-1 + 1/3 stim+2',
            'excluded_bad_behaviour_sessions': sorted(bad_behaviour_sessions),
            'mixed_model_converged': mixed_result.converged,
            'pyrup_amp_bounds_hz': PYRUP_AMP_BOUNDS,
            'pyrdown_amp_bounds_hz': PYRDOWN_AMP_BOUNDS,
            'pyrdown_delta_bounds_ms': PYRDOWN_DELTA_BOUNDS,
            'bin_edges': bin_edges,
            'bin_labels': bin_labels,
        },
        'profile_groups': {
            'pyrup': {
                'ctrl_profiles': mean_prof_ctrl_pyrup,
                'stim_profiles': mean_prof_stim_pyrup,
            },
            'pyrdown': {
                'ctrl_profiles': mean_prof_ctrl_pyrdown,
                'stim_profiles': mean_prof_stim_pyrdown,
            },
        },
        'binwise_tests': {
            'pyrup': [np.asarray(pvals, dtype=float) for pvals in pvals_pyrup],
            'pyrdown': [np.asarray(pvals, dtype=float) for pvals in pvals_pyrdown],
        },
        'session_delta_records': session_delta_records,
        'regression_fits': regression_fits,
        'mixed_model_fixed_effects': model_records,
    }

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, 'wb') as f:
        pickle.dump(payload, f, protocol=pickle.HIGHEST_PROTOCOL)
    pd.DataFrame(payload['session_delta_records']).to_csv(
        args.session_csv, index=False)
    cell_df.to_csv(cell_output_path, index=False)
    pd.DataFrame(event_time_records).to_csv(event_time_output_path, index=False)
    model_df.to_csv(model_output_path, index=False)
    print_files_saved([
        ('analysis data', args.output),
        ('session summary', args.session_csv),
        ('acute firing cells', cell_output_path),
        ('event-time effects', event_time_output_path),
        ('mixed model', model_output_path),
    ])

if __name__ == '__main__':
    main()
