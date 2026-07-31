# -*- coding: utf-8 -*-
'''
Created on Mon Oct 20 12:20:52 2025

The full GLM for predictor sieving

@author: Dinghao Luo

'''

#%% imports
from pathlib import Path
import sys

import pickle
import numpy as np
import pandas as pd
from collections import defaultdict
from sklearn.preprocessing import StandardScaler
from scipy import stats
import scipy.io as sio
import statsmodels.api as sm
from statsmodels.stats.outliers_influence import variance_inflation_factor

import matplotlib.pyplot as plt

repo_root = Path(__file__).resolve().parents[3]
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))

from common_functions import mpl_formatting
from console_formatting import print_session
import glm_functions as gf
import glm_permutation_functions as gpf
import project_paths as pp
mpl_formatting()

import rec_list
paths = rec_list.pathLC


#%% paths & params
LC_stem       = pp.LC_EPHYS_STEM
all_sess_stem = LC_stem / 'all_sessions'
GLM_stem      = pp.LC_EPHYS_FIGURES_STEM / 'GLM'

LC_beh_stem   = pp.behaviour_experiment_stem('LC')

(GLM_stem / 'single_cell_pred').mkdir(parents=True, exist_ok=True)
(GLM_stem / 'model_comparison').mkdir(parents=True, exist_ok=True)

SAMP_FREQ     = 1250  # Hz
SAMP_FREQ_BEH = 1000  # Hz
RUN_ONSET_IDX = 3 * SAMP_FREQ

eps = 1e-6

scaler = StandardScaler()

# permutation controls
n_shuffles = 500
rng_seed   = 42
rng        = np.random.default_rng(rng_seed)

# for block shuffle
block_size = None  # e.g., 5 to preserve short-range structure


#%% load cell table
print('loading data...')
cell_prop_path = LC_stem / 'LC_all_cell_profiles.pkl'
cell_prop = pd.read_pickle(cell_prop_path)


#%% accumulation structures
# observed per-predictor across cells
obs_delta_r2 = {}   # predictor -> list of values across cells
obs_delta_aic = {}
obs_lr_stat  = {}

# shuffle accumulators for per-predictor means (sums & counts per shuffle)
shuf_sum_delta_r2 = {}
shuf_cnt_delta_r2 = {}
shuf_sum_delta_aic = {}
shuf_cnt_delta_aic = {}
shuf_sum_lr = {}
shuf_cnt_lr  = {}

# full-model r² across cells
obs_r2 = []  # list across cells
shuf_r2_sum  = np.zeros(n_shuffles)  # per-shuffle mean numerators
shuf_r2_count= np.zeros(n_shuffles)  # per-shuffle valid counts

# store coefficients for reference (no inferential test here)
all_results = []
all_results_lr = []  # per-cell, per-predictor observed Δ-r²/aic/lr (for saving)


#%% main
pair_r_abs_cells = defaultdict(list)  # (p1,p2) -> list of |r| across cells
pair_r_cells     = defaultdict(list)  # signed r (optional, for sign info)
pair_r_abs_by_sess = defaultdict(lambda: defaultdict(list))  # (p1,p2)[recname] -> list of |r| in that session

vif_by_predictor = defaultdict(list)  # predictor -> list of VIFs across cells

# threshold to call a pair 'collinear'
r_abs_thr = 0.7

for path in paths:
    recname = Path(path).name
    print_session(recname)

    # load behaviour
    beh_path = LC_beh_stem / f'{recname}.pkl'
    with open(beh_path, 'rb') as f:
        beh = pickle.load(f)

    timestamps_ms = beh['upsampled_timestamps_ms']
    timestamps_s  = timestamps_ms / 1000.0
    speeds_cm_s   = beh['upsampled_speed_cm_s']
    lick_times    = [[lick[0] for lick in trial] for trial in beh['lick_times'][1:]]
    lick_times_s  = [lick/SAMP_FREQ_BEH for trial in lick_times for lick in trial]
    reward_times  = beh['reward_times'][1:]
    run_onsets    = beh['run_onsets'][1:]
    run_onsets_s  = [ro/SAMP_FREQ_BEH for ro in run_onsets]
    trials_sts    = beh['trial_statements'][1:]

    # find opto trials
    opto_idx = [i for i, t in enumerate(trials_sts) if t[15] != '0']
    first_opto = opto_idx[0] if opto_idx else len(trials_sts) - 1

    # new--get run-onset aligned to full spike maps
    rec_stem = pp.MICEEXP_ROOT / f'ANMD{recname[1:5]}' / recname[:14] / recname

    aligned_run_path = rec_stem / f'{recname}_DataStructure_mazeSection1_TrialType1_alignRun_msess1.mat'
    aligned_run = sio.loadmat(aligned_run_path)['trialsRun'][0][0]
    run_onsets_spike = aligned_run['startLfpInd'][0][1:]

    # get spike maps
    sess_stem = all_sess_stem / recname
    spike_maps = np.load(sess_stem / f'{recname}_smoothed_spike_map.npy', allow_pickle=True)

    # design matrix (behaviour-only base)
    beh_rows, kept_trials = [], []
    # for ti, onset in enumerate(run_onsets[:-1]):
    for ti, onset in enumerate(run_onsets[:first_opto]):
        if ti in opto_idx or ti-1 in opto_idx or np.isnan(onset):
            continue
        onset_time = onset / SAMP_FREQ_BEH

        # last trial's speed
        prev_speed = gf.mean_speed_prev_trial(timestamps_s, speeds_cm_s, run_onsets_s, ti)

        # time since last reward
        time_since_rew = gf.time_since_last_reward(reward_times, onset_time, ti)

        # first-lick-to-reward (prediction error-ish)
        pred_err = gf.first_lick_to_reward_last_trial(lick_times, reward_times, ti)

        # licks in the last 5 seconds
        licks_prev_5s = gf.lick_rate_last5s(lick_times_s, onset_time)

        ## ---- FEATURE LIST HERE ---- ##
        feats = [
            time_since_rew,
            pred_err,
            # trial_idx,  # removed c. collinearity criteria
            prev_speed,
            licks_prev_5s  # removed c. collinearity criteria
            ]
        ## ---- FEATURE LIST ENDS ---- ##

        beh_rows.append(feats)
        kept_trials.append(ti)

    if len(kept_trials) < 5:
        print(f'{recname}: too few kept trials; skipping session')
        continue

    X_behav_base = np.vstack(beh_rows)

    ## ---- FEATURE NAMES HERE ---- ##
    beh_names_base = [
        'Time since last reward',
        'Prediction error',
        # 'Trial index',  # removed c. collinearity criteria
        'Prev. trial speed',
        'Prev. 5 s licks'  # removed c. collinearity criteria
        ]
    ## ---- FEATURE NAMES END ---- ##

    # single cell spiking data
    all_trains_path = all_sess_stem / recname / f'{recname}_all_trains.npy'
    all_trains = np.load(all_trains_path, allow_pickle=True).item()

    curr_cell_prop = cell_prop[cell_prop['sessname'] == recname]
    for cluname, row in curr_cell_prop.iterrows():
        if row['identity'] == 'other' or not row['run_onset_peak']:
            continue

        trains = all_trains[cluname]

        # amplitudes
        amp_rows = []
        for ti in kept_trials:
            amp = gf.run_onset_amplitude(trains[ti], SAMP_FREQ, RUN_ONSET_IDX)
            amp_rows.append(amp)
        amp_rows = np.array(amp_rows, dtype=float)

        # pre-onset FR
        preonset_rows = []
        for ti in kept_trials:
            preonset = gf.preonset_rate(trains[ti])
            preonset_rows.append(preonset)
        preonset_rows = np.array(preonset_rows, dtype=float)

        # new--get baseline rate for -10~-2.5 s
        clu_idx = int(cluname.split('clu')[-1]) - 2  # for retrieving spike map
        base_5_rows  = []
        base_10_rows = []
        base_15_rows = []
        base_20_rows = []
        base_30_rows = []
        for ti in kept_trials:
            base_5_start  = int(run_onsets_spike[ti] - SAMP_FREQ * 5)
            base_10_start = int(run_onsets_spike[ti] - SAMP_FREQ * 10)
            base_15_start = int(run_onsets_spike[ti] - SAMP_FREQ * 15)
            base_20_start = int(run_onsets_spike[ti] - SAMP_FREQ * 20)
            base_30_start = int(run_onsets_spike[ti] - SAMP_FREQ * 30)
            base_end   = int(run_onsets_spike[ti] - SAMP_FREQ * 2.5)
            base_5  = np.nanmean(spike_maps[clu_idx, base_5_start : base_end])
            base_10 = np.nanmean(spike_maps[clu_idx, base_10_start : base_end])
            base_15 = np.nanmean(spike_maps[clu_idx, base_15_start : base_end])
            base_20 = np.nanmean(spike_maps[clu_idx, base_20_start : base_end])
            base_30 = np.nanmean(spike_maps[clu_idx, base_30_start : base_end])
            base_5_rows.append(base_5)
            base_10_rows.append(base_10)
            base_15_rows.append(base_15)
            base_20_rows.append(base_20)
            base_30_rows.append(base_30)
        base_5_rows = np.array(base_5_rows, dtype=float)
        base_10_rows = np.array(base_10_rows, dtype=float)
        base_15_rows = np.array(base_15_rows, dtype=float)
        base_20_rows = np.array(base_20_rows, dtype=float)
        base_30_rows = np.array(base_30_rows, dtype=float)

        # previous amplitude
        prev_amp = np.concatenate([[np.nan], amp_rows[:-1]])

        ## ---- FEATURE LIST (PER CELL) ---- ##
        X = np.column_stack([
            X_behav_base,  # see above
            base_5_rows,
            base_10_rows,
            base_15_rows,
            base_20_rows,
            base_30_rows,
            prev_amp,  # removed c. collinearity criteria
            preonset_rows  # removed c. collinearity criteria
            ])
        ## ---- FEATURE LIST (PER CELL) ENDS ---- ##

        ## ---- FEATURE NAMES (PER CELL) HERE ---- ##
        beh_names = beh_names_base + [
            'Baseline 5 s',
            'Baseline 10 s',
            'Baseline 15 s',
            'Baseline 20 s',
            'Baseline 30 s',
            'Prev. run-onset amp.',  # removed c. collinearity criteria
            'Pre-onset FR'  # removed c. collinearity criteria
            ]
        ## ---- FEATURE NAMES (PER CELL) END ---- ##

        ## new--test for collinearity
        # correlations and VIF use only trials with complete predictors
        X_diag = np.asarray(X)
        X_diag = X_diag[np.isfinite(X_diag).all(axis=1)]
        if X_diag.shape[0] >= 3:
            # zero-variance predictors have no correlation or VIF
            col_std = X_diag.std(axis=0, ddof=1)
            keep_cols = col_std > 0
            Xc = X_diag[:, keep_cols]
            names_c = [n for n, k in zip(beh_names, keep_cols) if k]

            if Xc.shape[1] >= 2:
                corr = pd.DataFrame(Xc, columns=names_c).corr()

                gpf.save_collinearity_heatmap(
                    corr,
                    cluname,
                    GLM_stem / 'model_comparison' / 'collinearity_heatmaps' / f'{cluname}_corr_heatmap'
                )

                # accumulate pairwise |r| across cells + by session
                cols = corr.columns.tolist()
                for i in range(len(cols)):
                    for j in range(i+1, len(cols)):
                        p1, p2 = cols[i], cols[j]
                        key = tuple(sorted((p1, p2)))
                        r_ij = corr.loc[p1, p2]
                        if np.isfinite(r_ij):
                            pair_r_abs_cells[key].append(abs(r_ij))
                            pair_r_cells[key].append(r_ij)
                            pair_r_abs_by_sess[key][recname].append(abs(r_ij))

            # VIF uses the same predictors as the correlation matrix
            X_vif = Xc
            if X_vif.shape[1] >= 2:
                vif_vals = [
                    variance_inflation_factor(X_vif, i)
                    for i in range(X_vif.shape[1])
                    ]
                for name, v in zip(names_c, vif_vals):
                    if np.isfinite(v):
                        vif_by_predictor[name].append(float(v))
        ## new--end of collinearity tests

        y = np.clip(amp_rows, eps, None)

        row_mask = np.isfinite(y) & np.all(np.isfinite(X), axis=1)
        row_mask &= ~(np.all(X == 0, axis=1))
        X = X[row_mask]
        y = y[row_mask]

        if len(y) < X.shape[1] + 2:
            print(f'{cluname}: too few valid trials; skipping')
            continue

        # z-score predictors once
        Xz = scaler.fit_transform(X)

        # fit full model
        res_full = gf.fit_glm_log_gaussian(Xz, y, eps=eps)

        y_pred_log = res_full.predict(sm.add_constant(Xz))
        y_pred     = np.exp(y_pred_log) - eps
        r_val      = np.corrcoef(y, y_pred)[0, 1]
        r2_val     = r_val**2
        obs_r2.append(r2_val)

        # store coefficients
        for name, coef in zip(beh_names, res_full.params[1:]):  # skip intercept
            all_results.append({'rec': recname,
                                'cell_id': cluname,
                                'predictor': name,
                                'coef': coef})

        # observed drop-one stats for this cell
        percell_obs_delta_r2  = {}
        percell_obs_delta_aic = {}
        percell_obs_lr        = {}

        for drop_name in beh_names:
            keep_idx = [i for i, n in enumerate(beh_names) if n != drop_name]
            X_red = Xz[:, keep_idx]
            res_red = gf.fit_glm_log_gaussian(X_red, y, eps=eps)

            y_pred_red_log = res_red.predict(sm.add_constant(X_red))
            y_pred_red     = np.exp(y_pred_red_log) - eps
            r_val_red      = np.corrcoef(y, y_pred_red)[0, 1]
            r2_red         = r_val_red**2

            delta_r2  = r2_val - r2_red
            llf_full  = res_full.llf
            llf_red   = res_red.llf
            lr_stat   = 2 * (llf_full - llf_red)
            delta_aic = res_red.aic - res_full.aic

            percell_obs_delta_r2[drop_name]  = delta_r2
            percell_obs_delta_aic[drop_name] = delta_aic
            percell_obs_lr[drop_name]        = lr_stat

            all_results_lr.append({'rec': recname,
                                   'cell_id': cluname,
                                   'dropped': drop_name,
                                   'delta_r2': delta_r2,
                                   'aic_full': res_full.aic,
                                   'aic_reduced': res_red.aic,
                                   'delta_aic': delta_aic,
                                   'lr_stat': lr_stat})

        # accumulate observed into population containers
        for name in percell_obs_delta_r2:
            obs_delta_r2.setdefault(name, []).append(percell_obs_delta_r2[name])
            obs_delta_aic.setdefault(name, []).append(percell_obs_delta_aic[name])
            obs_lr_stat.setdefault(name, []).append(percell_obs_lr[name])

            # initialise shuffle accumulators if first time we see this predictor
            if name not in shuf_sum_delta_r2:
                shuf_sum_delta_r2[name]  = np.zeros(n_shuffles)
                shuf_cnt_delta_r2[name]  = np.zeros(n_shuffles)
                shuf_sum_delta_aic[name] = np.zeros(n_shuffles)
                shuf_cnt_delta_aic[name] = np.zeros(n_shuffles)
                shuf_sum_lr[name]        = np.zeros(n_shuffles)
                shuf_cnt_lr[name]        = np.zeros(n_shuffles)

        # build shuffle nulls for this cell
        for s in range(n_shuffles):
            # shuffle y (or block-shuffle)
            y_shuf = gpf.block_permute(y, rng, block_size)

            # full model on shuffled
            res_full_sh = gf.fit_glm_log_gaussian(Xz, y_shuf, eps=eps)

            y_pred_sh_log = res_full_sh.predict(sm.add_constant(Xz))
            y_pred_sh     = np.exp(y_pred_sh_log) - eps
            r_val_sh      = np.corrcoef(y_shuf, y_pred_sh)[0, 1]
            if np.isfinite(r_val_sh):
                r2_sh = r_val_sh**2
                shuf_r2_sum[s]   += r2_sh
                shuf_r2_count[s] += 1

            # drop-one on shuffled
            for drop_name in percell_obs_delta_r2.keys():
                keep_idx = [i for i, n in enumerate(beh_names) if n != drop_name]
                X_red = Xz[:, keep_idx]
                res_red_sh = gf.fit_glm_log_gaussian(X_red, y_shuf, eps=eps)

                # Δr² on shuffled
                y_pred_red_sh_log = res_red_sh.predict(sm.add_constant(X_red))
                y_pred_red_sh     = np.exp(y_pred_red_sh_log) - eps
                r_val_red_sh      = np.corrcoef(y_shuf, y_pred_red_sh)[0, 1]
                if not np.isfinite(r_val_red_sh) or not np.isfinite(r_val_sh):
                    continue
                delta_r2_sh = (r_val_sh**2) - (r_val_red_sh**2)

                # Δaic and lr on shuffled
                llf_full_sh = res_full_sh.llf
                llf_red_sh  = res_red_sh.llf
                lr_sh       = 2 * (llf_full_sh - llf_red_sh)
                delta_aic_sh= res_red_sh.aic - res_full_sh.aic

                # accumulate predictor-wise shuffle sums/counts for population mean
                shuf_sum_delta_r2[drop_name][s]  += delta_r2_sh
                shuf_cnt_delta_r2[drop_name][s]  += 1
                shuf_sum_delta_aic[drop_name][s] += delta_aic_sh
                shuf_cnt_delta_aic[drop_name][s] += 1
                shuf_sum_lr[drop_name][s]        += lr_sh
                shuf_cnt_lr[drop_name][s]        += 1


#%% summarise collinearity across all cells/sessions
gpf.summarise_collinearity(
    pair_r_abs_cells,
    pair_r_abs_by_sess,
    vif_by_predictor,
    GLM_stem,
    r_abs_thr,
    output_prefix='model_comparison'
)


#%% aggregate to population summaries and empirical p-values
results_df = pd.DataFrame(all_results)
lr_df      = pd.DataFrame(all_results_lr)

# full-model r² population mean ± sem
obs_r2 = np.array(obs_r2, dtype=float)
mean_r2_obs = np.nanmean(obs_r2) if len(obs_r2) else np.nan
sem_r2_obs  = stats.sem(obs_r2, nan_policy='omit') if len(obs_r2) > 1 else np.nan

# per-shuffle mean r² across cells
shuf_r2_mean = shuf_r2_sum / shuf_r2_count
shuf_r2_mean = shuf_r2_mean[np.isfinite(shuf_r2_mean)]

# empirical p for r² mean
if len(shuf_r2_mean):
    p_r2 = (1 + np.sum(shuf_r2_mean >= mean_r2_obs)) / (1 + len(shuf_r2_mean))
else:
    p_r2 = np.nan

# per-predictor: compute observed mean±sem and shuffle-mean distribution
pred_order = []
pred_stats = []

for name in sorted(obs_delta_r2.keys()):
    # observed
    vals_r2   = np.array(obs_delta_r2[name], dtype=float)
    vals_aic  = np.array(obs_delta_aic[name], dtype=float)
    vals_lr   = np.array(obs_lr_stat[name],  dtype=float)

    m_r2  = np.nanmean(vals_r2) if len(vals_r2) else np.nan
    s_r2  = stats.sem(vals_r2, nan_policy='omit') if len(vals_r2) > 1 else np.nan
    m_aic = np.nanmean(vals_aic) if len(vals_aic) else np.nan
    s_aic = stats.sem(vals_aic, nan_policy='omit') if len(vals_aic) > 1 else np.nan
    m_lr  = np.nanmean(vals_lr) if len(vals_lr) else np.nan
    s_lr  = stats.sem(vals_lr, nan_policy='omit') if len(vals_lr) > 1 else np.nan

    # shuffle-mean distributions
    sh_mean_r2  = shuf_sum_delta_r2[name]  / shuf_cnt_delta_r2[name]
    sh_mean_aic = shuf_sum_delta_aic[name] / shuf_cnt_delta_aic[name]
    sh_mean_lr  = shuf_sum_lr[name]        / shuf_cnt_lr[name]

    sh_mean_r2  = sh_mean_r2[np.isfinite(sh_mean_r2)]
    sh_mean_aic = sh_mean_aic[np.isfinite(sh_mean_aic)]
    sh_mean_lr  = sh_mean_lr[np.isfinite(sh_mean_lr)]

    p_emp_r2  = (1 + np.sum(sh_mean_r2  >= m_r2 )) / (1 + len(sh_mean_r2 )) if len(sh_mean_r2 ) else np.nan
    p_emp_aic = (1 + np.sum(sh_mean_aic >= m_aic)) / (1 + len(sh_mean_aic)) if len(sh_mean_aic) else np.nan
    p_emp_lr  = (1 + np.sum(sh_mean_lr  >= m_lr )) / (1 + len(sh_mean_lr )) if len(sh_mean_lr ) else np.nan

    pred_order.append(name)
    pred_stats.append({
        'predictor': name,
        'mean_delta_r2': m_r2, 'sem_delta_r2': s_r2, 'p_emp_delta_r2': p_emp_r2,
        'mean_delta_aic': m_aic, 'sem_delta_aic': s_aic, 'p_emp_delta_aic': p_emp_aic,
        'mean_lr': m_lr, 'sem_lr': s_lr, 'p_emp_lr': p_emp_lr
    })

pred_df = pd.DataFrame(pred_stats).set_index('predictor')

# order by descending
order_r2 = pred_df['mean_delta_r2'].sort_values(ascending=False).index.tolist()
order_aic= pred_df['mean_delta_aic'].sort_values(ascending=False).index.tolist()
order_lr = pred_df['mean_lr'].sort_values(ascending=False).index.tolist()


#%% plotting - ΔR² by predictor (population mean; p from shuffle means)
m_obs = pred_df.loc[order_r2, 'mean_delta_r2'].values
s_obs = pred_df.loc[order_r2, 'sem_delta_r2'].values

# compute shuffle percentile thresholds
null_thr = {}
for name in order_r2:
    sh_mean = shuf_sum_delta_r2[name] / shuf_cnt_delta_r2[name]  # sum / cell count
    sh_mean = sh_mean[np.isfinite(sh_mean)]
    null_thr[name] = {
        '95':   np.percentile(sh_mean, 95),
        '99':   np.percentile(sh_mean, 99),
        '99.9': np.percentile(sh_mean, 99.9),
        '99.99':np.percentile(sh_mean, 99.99),
        'mean': np.mean(sh_mean)
    }

# assign stars
stars = []
for name, mo in zip(order_r2, m_obs):
    t = null_thr[name]
    if mo >= t['99.99']:
        stars.append('****')
    elif mo >= t['99.9']:
        stars.append('***')
    elif mo >= t['99']:
        stars.append('**')
    elif mo >= t['95']:
        stars.append('*')
    else:
        stars.append('')

# plotting
fig, ax = plt.subplots(figsize=(4.2, 3))

ax.axvline(0, color='k', lw=0.75)

ax.barh(order_r2, m_obs, xerr=s_obs,
        color='#bdbdbd', edgecolor='k', capsize=2, height=0.6, label='Observed')

# overlay shuffle
means = [null_thr[n]['mean'] for n in order_r2]
ax.barh(order_r2, means,
        color='white', edgecolor='k', lw=0.8, height=0.6, alpha=0.5,
        label='Shuffle mean')

for i, (mo, so, st) in enumerate(zip(m_obs, s_obs, stars)):
    if st:
        ax.text(mo + np.sign(mo if mo != 0 else 1)*0.002, i+.1, st, va='center', ha='left',
                fontsize=8, color='k', fontweight='bold')
    ax.text(mo + np.sign(mo if mo != 0 else 1)*0.002, i + .25,
            f'{mo:.4f}±{so:.4f}', va='center', ha='left', fontsize=5)

# labels & style
ax.set(
    xlabel='Unique variance explained (ΔR²)',
    ylabel='',
    title='Drop-one-out ΔR² vs shuffle-null thresholds'
)
ax.spines[['top','right']].set_visible(False)
ax.legend(frameon=False, fontsize=8)

plt.tight_layout()
for ext in ['.pdf', '.png']:
    fig.savefig(GLM_stem / f'model_comparison/deltaR2_by_predictor_permutation{ext}',
                dpi=300,
                bbox_inches='tight')


#%% plotting - ΔAIC by predictor (percentile-based stars)
m_obs = pred_df.loc[order_aic, 'mean_delta_aic'].values
s_obs = pred_df.loc[order_aic, 'sem_delta_aic'].values

# compute shuffle percentile thresholds
null_thr = {}
for name in order_aic:
    sh_mean = shuf_sum_delta_aic[name] / shuf_cnt_delta_aic[name]
    sh_mean = sh_mean[np.isfinite(sh_mean)]
    null_thr[name] = {
        '95':   np.percentile(sh_mean, 95),
        '99':   np.percentile(sh_mean, 99),
        '99.9': np.percentile(sh_mean, 99.9),
        '99.99':np.percentile(sh_mean, 99.99),
        'mean': np.mean(sh_mean)
    }

# assign stars
stars = []
for name, mo in zip(order_aic, m_obs):
    t = null_thr[name]
    if mo >= t['99.99']:
        stars.append('****')
    elif mo >= t['99.9']:
        stars.append('***')
    elif mo >= t['99']:
        stars.append('**')
    elif mo >= t['95']:
        stars.append('*')
    else:
        stars.append('')

# plotting
fig, ax = plt.subplots(figsize=(4.2, 3))
ax.axvline(0, color='k', lw=0.75)

ax.barh(order_aic, m_obs, xerr=s_obs,
        color='#bdbdbd', edgecolor='k', capsize=2, height=0.6, label='Observed')

# overlay shuffle mean
means = [null_thr[n]['mean'] for n in order_aic]
ax.barh(order_aic, means,
        color='white', edgecolor='k', lw=0.8, height=0.6, alpha=0.5, label='Shuffle mean')

# annotate
for i, (mo, so, st) in enumerate(zip(m_obs, s_obs, stars)):
    if st:
        ax.text(mo + np.sign(mo if mo != 0 else 1)*0.002, i+.1, st, va='center', ha='left',
                fontsize=8, color='k', fontweight='bold')
    ax.text(mo + np.sign(mo if mo != 0 else 1)*0.002, i + .25,
            f'{mo:.4f}±{so:.4f}', va='center', ha='left', fontsize=5)

ax.set(
    xlabel='ΔAIC (Reduced - Full)',
    ylabel='',
    title='Model AIC comparison vs shuffle-null thresholds'
)
ax.spines[['top','right']].set_visible(False)
ax.legend(frameon=False, fontsize=8)

plt.tight_layout()
for ext in ['.pdf', '.png']:
    fig.savefig(GLM_stem / f'model_comparison/model_AIC_comparison_permutation{ext}',
                dpi=300, bbox_inches='tight')


#%% plotting - Likelihood-ratio statistic (percentile-based stars)
m_obs = pred_df.loc[order_lr, 'mean_lr'].values
s_obs = pred_df.loc[order_lr, 'sem_lr'].values

# compute shuffle percentile thresholds
null_thr = {}
for name in order_lr:
    sh_mean = shuf_sum_lr[name] / shuf_cnt_lr[name]
    sh_mean = sh_mean[np.isfinite(sh_mean)]
    null_thr[name] = {
        '95':   np.percentile(sh_mean, 95),
        '99':   np.percentile(sh_mean, 99),
        '99.9': np.percentile(sh_mean, 99.9),
        '99.99':np.percentile(sh_mean, 99.99),
        'mean': np.mean(sh_mean)
    }

# assign stars
stars = []
for name, mo in zip(order_lr, m_obs):
    t = null_thr[name]
    if mo >= t['99.99']:
        stars.append('****')
    elif mo >= t['99.9']:
        stars.append('***')
    elif mo >= t['99']:
        stars.append('**')
    elif mo >= t['95']:
        stars.append('*')
    else:
        stars.append('')

# plotting
fig, ax = plt.subplots(figsize=(4.2, 3))
ax.axvline(0, color='k', lw=0.75)

ax.barh(order_lr, m_obs, xerr=s_obs,
        color='#bdbdbd', edgecolor='k', capsize=2, height=0.6, label='Observed')

# overlay shuffle mean
means = [null_thr[n]['mean'] for n in order_lr]
ax.barh(order_lr, means,
        color='white', edgecolor='k', lw=0.8, height=0.6, alpha=0.5, label='Shuffle mean')

# annotate
for i, (mo, so, st) in enumerate(zip(m_obs, s_obs, stars)):
    if st:
        ax.text(mo + np.sign(mo if mo != 0 else 1)*0.002, i+.1, st, va='center', ha='left',
                fontsize=8, color='k', fontweight='bold')
    ax.text(mo + np.sign(mo if mo != 0 else 1)*0.002, i + .25,
            f'{mo:.4f}±{so:.4f}', va='center', ha='left', fontsize=5)

ax.set(
    xlabel='2×Δlog-likelihood (Λ)',
    ylabel='',
    title='Likelihood-ratio vs shuffle-null thresholds'
)
ax.spines[['top','right']].set_visible(False)
ax.legend(frameon=False, fontsize=8)

plt.tight_layout()
for ext in ['.pdf', '.png']:
    fig.savefig(GLM_stem / f'model_comparison/likelihood_ratio_stat_permutation{ext}',
                dpi=300, bbox_inches='tight')


#%% plotting - regression coefficients (population mean ± sem)
coef_df = pd.DataFrame(all_results)

# aggregate mean ± sem per predictor
coef_summary = (
    coef_df.groupby('predictor')['coef']
    .agg(['mean', 'sem'])
    .reset_index()
    .sort_values('mean', ascending=False)
)

# plotting
fig, ax = plt.subplots(figsize=(4.2, 3))
ax.axvline(0, color='k', lw=0.75)

ax.barh(coef_summary['predictor'], coef_summary['mean'],
        xerr=coef_summary['sem'], color='#bdbdbd',
        edgecolor='k', capsize=2, height=0.6)

# annotate mean ± sem
for i, (m, s) in enumerate(zip(coef_summary['mean'], coef_summary['sem'])):
    ax.text(m + np.sign(m if m != 0 else 1)*0.002, i + .1,
            f'{m:.4f}±{s:.4f}', va='center', ha='left', fontsize=5)

ax.set(
    xlabel='Regression coefficient (ß)',
    ylabel='',
    title='Mean ± SEM regression coefficients across LC cells'
)
ax.spines[['top','right']].set_visible(False)
plt.tight_layout()

for ext in ['.pdf', '.png']:
    fig.savefig(GLM_stem / f'model_comparison/coefficients_mean_sem{ext}',
                dpi=300, bbox_inches='tight')
