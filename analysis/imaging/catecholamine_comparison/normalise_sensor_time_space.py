# -*- coding: utf-8 -*-
'''
Created on Mon Jun 29 2026

normalise dLight and nLight comparisons against source-paper sensor scales

@author: Dinghao Luo
'''

#%% imports
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import csv
import re
import sys

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from scipy.stats import ranksums, sem
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

repo_root = Path(__file__).resolve().parents[3]
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))

from common_functions import mpl_formatting
from console_formatting import print_files_saved, print_statistics_section
import project_paths as pp
mpl_formatting()


#%% parameters
DATASETS = [
    {
        'transmitter': 'dopamine',
        'sensor': 'dLight',
        'source_sensor': 'dLight3.6',
        'colour': 'darkgreen',
        },
    {
        'transmitter': 'norepinephrine',
        'sensor': 'nLight',
        'source_sensor': 'nLightG',
        'colour': 'royalblue',
        },
    ]

WINDOW_ORDER = ['0-1 s', '1-2 s', '2-3 s', '3-4 s']
WINDOW_CENTRES = np.asarray([0.5, 1.5, 2.5, 3.5])

NUMBER_PATTERN = r'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?'

data_stem = pp.DATA_ROOT / 'catecholamine_comparison'
save_stem = pp.FIGURES_ROOT / 'catecholamine_comparison'
source_stem = (
    pp.FIGURES_ROOT
    / 'fig_6_lc_da_ca1_model'
    / 'parameter_sensitivity'
    / 'source_paper_metrics'
    )

data_stem.mkdir(parents=True, exist_ok=True)
save_stem.mkdir(parents=True, exist_ok=True)


#%% sensor scales
@dataclass
class SensorScale:
    transmitter: str
    sensor: str
    source_sensor: str
    ligand: str
    affinity_nm: float
    dff_max_fraction: float
    hill: float
    tau_off_s: float
    half_decay_s: float
    source_note: str

def _hill_response(conc, bottom, top, ec50, hill):
    conc = np.asarray(conc, dtype=float)
    return bottom + top * (conc ** hill) / ((ec50 ** hill) + (conc ** hill))

def _fit_titration_ec50(conc, resp):
    conc = np.asarray(conc, dtype=float)
    resp = np.asarray(resp, dtype=float)
    good = np.isfinite(conc) & np.isfinite(resp)
    conc = conc[good]
    resp = resp[good]

    if conc.size < 5 or np.nanmax(conc) <= 0:
        raise ValueError('nLight titration does not contain enough valid concentrations')

    popt, _ = curve_fit(
        _hill_response,
        conc,
        resp,
        p0=[float(np.nanmin(resp)), float(np.nanmax(resp)), 1e-6, 1.0],
        bounds=([-np.inf, 0, 1e-10, 0.1], [np.inf, np.inf, 1e-2, 4.0]),
        maxfev=20000,
        )

    bottom, top, ec50, hill = popt
    return float(ec50), float(hill), float(bottom + top)

def _extract_nlight_constants():
    path = source_stem / 'nlight_fig1.xlsx'
    # nLightG source: Nature Methods 2023 Fig. 1 source data, sheet
    # 'e) titrations nLightG, GRABNE'. The sheet stores replicate titration
    # traces, so I fit a Hill curve here instead of hard-coding an EC50.
    tit = pd.read_excel(path, sheet_name='e) titrations nLightG, GRABNE', header=None)
    conc_all = pd.to_numeric(tit.iloc[:, 0], errors='coerce')
    first_block = []
    for idx, val in conc_all.items():
        if np.isfinite(val):
            first_block.append(idx)
        elif len(first_block):
            break

    tit = tit.loc[first_block, :]
    conc = pd.to_numeric(tit.iloc[:, 0], errors='coerce').to_numpy(dtype=float)
    reps = tit.iloc[:, 1:7].apply(pd.to_numeric, errors='coerce')
    resp = reps.mean(axis=1, skipna=True).to_numpy(dtype=float)
    ec50_m, hill, dff_max = _fit_titration_ec50(conc, resp)

    # nLightG kinetics source: Nature Methods 2023 Fig. 1 source data, sheet
    # 'i) kinetic comparisons'. It reports tau_on and tau_off directly in s.
    kin = pd.read_excel(path, sheet_name='i) kinetic comparisons', header=None)
    lower = kin.astype(str).apply(lambda col: col.str.strip().str.lower())
    nlight_rows = np.where(lower.eq('nlightg').any(axis=1))[0]
    start = int(nlight_rows[-1])
    mean_rows = np.where(lower.iloc[start:, :].eq('mean').any(axis=1))[0]
    row = kin.iloc[start + int(mean_rows[0]), :]
    vals = pd.to_numeric(row, errors='coerce').dropna().to_numpy(dtype=float)
    tau_off_s = float(vals[1])

    return SensorScale(
        transmitter='norepinephrine',
        sensor='nLight',
        source_sensor='nLightG',
        ligand='NE',
        affinity_nm=float(ec50_m * 1e9),
        dff_max_fraction=float(dff_max),
        hill=float(hill),
        tau_off_s=float(tau_off_s),
        half_decay_s=float(tau_off_s * np.log(2)),
        source_note=(
            'Nature Methods 2023 nLight Fig. 1 source data: '
            'NE titration fit from sheet e; tau_off from sheet i'
            ),
        )

def _extract_dlight_constants():
    table2 = source_stem / 'dlight_rs7313638_supp_2.xlsx'
    table6 = source_stem / 'dlight_rs7313638_supp_6.xlsx'

    source_note = (
        'Research Square rs-7313638/v1 supplementary table 2 for '
        'dLight3.6 DA EC50 and dynamic range; supplementary table 6 for '
        'in vivo half-decay'
        )

    # dLight3.6 source: Research Square rs-7313638/v1 supplementary
    # table 2. This is the in vitro DA row, giving EC50/Kd and dynamic
    # range. The table does not give a Hill coefficient, so the
    # occupancy-to-C/EC50 conversion below uses Hill=1 for dLight.
    df = pd.read_excel(table2, sheet_name='Table2', header=None)
    rows = df.astype(str).apply(
        lambda row: row.str.contains('dLight 3.6', case=False, na=False).any(),
        axis=1,
        )
    row = df.loc[rows].iloc[0, :]
    strings = [str(x) for x in row.to_list()]
    nm_vals = []
    pct_vals = []
    for x in strings:
        if re.search(r'nm', x, flags=re.IGNORECASE):
            nm_vals.append(float(re.search(NUMBER_PATTERN, x).group(0)))
        if '%' in x:
            pct_vals.append(float(re.search(NUMBER_PATTERN, x).group(0)))
    affinity_nm = nm_vals[0]
    dff_max_fraction = pct_vals[0] / 100

    # dLight3.6 source: Research Square rs-7313638/v1 supplementary
    # table 6. It reports in vivo off kinetic half-decay, not tau_off.
    df = pd.read_excel(table6, sheet_name='Table6', header=None)
    rows = df.astype(str).apply(
        lambda row: row.str.contains('dLight3.6', case=False, na=False).any(),
        axis=1,
        )
    strings = [str(x) for x in df.loc[rows].iloc[0, :].to_list()]
    ms_vals = []
    for x in strings:
        if re.search(r'ms', x, flags=re.IGNORECASE):
            ms_vals.append(float(re.search(NUMBER_PATTERN, x).group(0)))
    half_decay_s = ms_vals[-1] / 1000

    return SensorScale(
        transmitter='dopamine',
        sensor='dLight',
        source_sensor='dLight3.6',
        ligand='DA',
        affinity_nm=float(affinity_nm),
        dff_max_fraction=float(dff_max_fraction),
        hill=1.0,
        tau_off_s=float(half_decay_s / np.log(2)),
        half_decay_s=float(half_decay_s),
        source_note=source_note,
        )


#%% normalisation and summaries
def add_sensor_columns(df, scales):
    out = df.copy()
    for col in [
            'sensor_affinity_nm', 'sensor_dff_max_fraction', 'sensor_hill',
            'sensor_tau_off_s', 'sensor_half_decay_s',
            ]:
        out[col] = np.nan
    for col in ['source_sensor', 'ligand', 'source_note']:
        out[col] = ''

    for transmitter, scale in scales.items():
        idx = out['transmitter'] == transmitter
        out.loc[idx, 'source_sensor'] = scale.source_sensor
        out.loc[idx, 'ligand'] = scale.ligand
        out.loc[idx, 'sensor_affinity_nm'] = scale.affinity_nm
        out.loc[idx, 'sensor_dff_max_fraction'] = scale.dff_max_fraction
        out.loc[idx, 'sensor_hill'] = scale.hill
        out.loc[idx, 'sensor_tau_off_s'] = scale.tau_off_s
        out.loc[idx, 'sensor_half_decay_s'] = scale.half_decay_s
        out.loc[idx, 'source_note'] = scale.source_note

    return out

def add_apparent_ligand_columns(df, dff_col, prefix):
    out = df.copy()
    dff = pd.to_numeric(out[dff_col], errors='coerce').to_numpy(dtype=float)
    dff_max = pd.to_numeric(out['sensor_dff_max_fraction'], errors='coerce').to_numpy(dtype=float)
    hill = pd.to_numeric(out['sensor_hill'], errors='coerce').to_numpy(dtype=float)
    affinity = pd.to_numeric(out['sensor_affinity_nm'], errors='coerce').to_numpy(dtype=float)

    occupancy = dff / dff_max
    c_over_ec50 = (occupancy / (1 - occupancy)) ** (1 / hill)

    out[f'{prefix}_sensor_occupancy'] = occupancy
    out[f'{prefix}_C_over_EC50'] = c_over_ec50
    out[f'{prefix}_apparent_ligand_nm'] = c_over_ec50 * affinity
    return out

def summarise_by_transmitter(df, key):
    rows = []
    for dataset in DATASETS:
        vals = pd.to_numeric(
            df.loc[df['transmitter'] == dataset['transmitter'], key],
            errors='coerce',
            ).dropna().to_numpy(dtype=float)
        rows.append({
            'transmitter': dataset['transmitter'],
            'sensor': dataset['sensor'],
            'metric': key,
            'n': int(vals.size),
            'median': float(np.nanmedian(vals)),
            'mean': float(np.nanmean(vals)),
            'sem': float(sem(vals, nan_policy='omit')),
            })
    return rows

def plot_violin(ax, df, key, ylabel, log_y=False):
    groups = []
    for dataset in DATASETS:
        vals = pd.to_numeric(
            df.loc[df['transmitter'] == dataset['transmitter'], key],
            errors='coerce',
            ).dropna().to_numpy(dtype=float)
        groups.append(vals)
    colours = [dataset['colour'] for dataset in DATASETS]
    labels = [
        f'{dataset["transmitter"]}\n{dataset["source_sensor"]}'
        for dataset in DATASETS
        ]

    clean_groups = [
        vals[vals > 0] if log_y else vals
        for vals in groups
        ]

    parts = ax.violinplot(
        clean_groups,
        positions=[1, 2],
        showmedians=True,
        showextrema=False,
        )
    for body, colour in zip(parts['bodies'], colours):
        body.set_facecolor(colour)
        body.set_edgecolor('none')
        body.set_alpha(0.62)
    parts['cmedians'].set_color('k')

    rng = np.random.default_rng(729)
    for xpos, vals, colour in zip([1, 2], clean_groups, colours):
        jitter = rng.normal(0, 0.035, size=vals.size)
        ax.scatter(
            np.full(vals.size, xpos) + jitter,
            vals,
            s=12,
            color=colour,
            edgecolor='none',
            alpha=0.65,
            zorder=3,
            )

    stat, p_val = ranksums(clean_groups[0], clean_groups[1])
    ax.text(
        0.05,
        0.96,
        f'p={p_val:.2e}',
        transform=ax.transAxes,
        ha='left',
        va='top',
        fontsize=7,
        )
    print(f'{key}: ranksums z={stat:.3f}, p={p_val:.3e}')

    if log_y:
        ax.set_yscale('log')
    ax.set(
        xticks=[1, 2],
        xticklabels=labels,
        xlim=(0.5, 2.5),
        ylabel=ylabel,
        )
    ax.spines[['top', 'right']].set_visible(False)

def spatial_matrix(df, transmitter, source, key='spatial_tau_px'):
    curr = df[
        (df['transmitter'] == transmitter)
        & (df['spatial_source'] == source)
        ].copy()
    recnames = sorted(curr['recname'].dropna().unique())
    mat = np.full((len(recnames), len(WINDOW_ORDER)), np.nan)

    for i, recname in enumerate(recnames):
        rec_rows = curr[curr['recname'] == recname]
        for j, window in enumerate(WINDOW_ORDER):
            vals = pd.to_numeric(
                rec_rows.loc[rec_rows['window'] == window, key],
                errors='coerce',
                ).dropna().to_numpy(dtype=float)
            mat[i, j] = vals[0]

    return mat

def plot_spatial_lines(ax, spatial, source, key, ylabel, ylim=None):
    for dataset in DATASETS:
        mat = spatial_matrix(spatial, dataset['transmitter'], source, key=key)
        mean = np.nanmean(mat, axis=0)
        err = sem(mat, axis=0, nan_policy='omit')
        n_sess = int(np.sum(np.any(np.isfinite(mat), axis=1)))
        ax.plot(
            WINDOW_CENTRES,
            mean,
            color=dataset['colour'],
            linewidth=1.3,
            marker='o',
            markersize=3.5,
            label=f'{dataset["transmitter"]} (n={n_sess})',
            )
        ax.fill_between(
            WINDOW_CENTRES,
            mean - err,
            mean + err,
            color=dataset['colour'],
            alpha=0.18,
            edgecolor='none',
            )

    ax.set(
        xticks=WINDOW_CENTRES,
        xticklabels=WINDOW_ORDER,
        xlabel='post-stim. window',
        ylabel=ylabel,
        )
    if ylim is not None:
        ax.set_ylim(*ylim)
    ax.spines[['top', 'right']].set_visible(False)


#%% main
def main():
    scales = {
        'dopamine': _extract_dlight_constants(),
        'norepinephrine': _extract_nlight_constants(),
        }
    rows = []
    for dataset in DATASETS:
        scale = scales[dataset['transmitter']]
        rows.append({
            'transmitter': scale.transmitter,
            'sensor': scale.sensor,
            'source_sensor': scale.source_sensor,
            'ligand': scale.ligand,
            'affinity_nm': scale.affinity_nm,
            'dff_max_fraction': scale.dff_max_fraction,
            'hill': scale.hill,
            'tau_off_s': scale.tau_off_s,
            'half_decay_s': scale.half_decay_s,
            'source_note': scale.source_note,
            })
    constants_path = data_stem / 'sensor_normalisation_constants.csv'
    with open(constants_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    release = pd.read_csv(data_stem / 'roi_release_metrics.csv')
    whole = pd.read_csv(data_stem / 'wholefield_stim_profile_metrics.csv')
    roi_time = pd.read_csv(data_stem / 'roi_neuropil_timecourse_metrics.csv')

    keys = ['transmitter', 'sensor', 'recname', 'animal']
    session = release.merge(whole, on=keys, how='outer')
    session = session.merge(roi_time, on=keys, how='outer', suffixes=('', '_roi_time'))
    session = add_sensor_columns(session, scales)

    # The old ROI scripts store RI as fractional dF/F, whereas the whole-field
    # profile table stores percent dF/F. Both are reduced to an apparent
    # fractional occupancy before the C/EC50 mapping.
    session = add_apparent_ligand_columns(session, 'mean_RI', 'roi_mean')
    session['wholefield_peak_dff_fraction'] = session['peak_dff_percent'] / 100
    session = add_apparent_ligand_columns(
        session,
        'wholefield_peak_dff_fraction',
        'wholefield_peak',
        )

    session['roi_tau_sensor_ratio'] = (
        session['roi_tau_s'] / session['sensor_tau_off_s']
        )
    session['wholefield_half_sensor_ratio'] = (
        session['half_recovery_s'] / session['sensor_half_decay_s']
        )
    rows = []
    for fname, source in [
            ('spatial_tau_metrics.csv', 'releasing_rois'),
            ('spatial_tau_metrics_all_rois.csv', 'all_rois'),
            ]:
        curr = pd.read_csv(data_stem / fname)
        curr['spatial_source'] = source
        rows.append(curr)

    spatial = pd.concat(rows, ignore_index=True)
    spatial = add_sensor_columns(spatial, scales)

    spatial['spatial_tau_rel_0_1'] = np.nan
    for (source, _), idx in spatial.groupby(['spatial_source', 'recname']).groups.items():
        curr = spatial.loc[idx]
        base = pd.to_numeric(
            curr.loc[curr['window'] == '0-1 s', 'spatial_tau_px'],
            errors='coerce',
            ).dropna()
        spatial.loc[idx, 'spatial_tau_rel_0_1'] = (
            spatial.loc[idx, 'spatial_tau_px'] / float(base.iloc[0])
            )

    session_path = data_stem / 'sensor_normalised_session_metrics.csv'
    spatial_path = data_stem / 'sensor_normalised_spatial_metrics.csv'
    session.to_csv(session_path, index=False)
    spatial.to_csv(spatial_path, index=False)

    saved = [
        ('csv', constants_path),
        ('csv', session_path),
        ('csv', spatial_path),
        ]

    fig, axs = plt.subplots(2, 3, figsize=(8.8, 5.2))

    plot_violin(
        axs[0, 0],
        session,
        'wholefield_peak_C_over_EC50',
        'whole-field peak C/EC50',
        log_y=True,
        )
    axs[0, 0].set_title('source-scale amplitude')

    plot_violin(
        axs[0, 1],
        session,
        'roi_mean_C_over_EC50',
        'mean ROI C/EC50',
        log_y=True,
        )
    axs[0, 1].set_title('ROI release amplitude')

    plot_violin(
        axs[0, 2],
        session,
        'roi_tau_sensor_ratio',
        r'ROI decay $\tau$ / sensor $\tau$',
        log_y=True,
        )
    axs[0, 2].set_title('ROI time scale')

    plot_violin(
        axs[1, 0],
        session,
        'wholefield_half_sensor_ratio',
        'whole-field half-recovery / sensor half-time',
        log_y=True,
        )
    axs[1, 0].set_title('whole-field time scale')

    plot_spatial_lines(
        axs[1, 1],
        spatial,
        'releasing_rois',
        'spatial_tau_px',
        r'spatial $\tau$ (px)',
        )
    axs[1, 1].set_title('spread around releasing ROIs')

    plot_spatial_lines(
        axs[1, 2],
        spatial,
        'releasing_rois',
        'spatial_tau_rel_0_1',
        r'spatial $\tau$ / 0-1 s',
        )
    axs[1, 2].axhline(1, color='0.3', linestyle='--', linewidth=0.8)
    axs[1, 2].set_title('spatial retention to 4 s')
    axs[1, 2].legend(frameon=False, fontsize=6, loc='upper left')

    fig.tight_layout(w_pad=0.9, h_pad=1.1)
    for ext in ['.png', '.pdf']:
        out = save_stem / f'dopamine_norepinephrine_sensor_normalised_time_space{ext}'
        fig.savefig(out, dpi=300, bbox_inches='tight')
        saved.append((ext[1:], out))
    plt.close(fig)

    fig, axs = plt.subplots(2, 2, figsize=(7.4, 4.8), sharex=True)

    plot_spatial_lines(
        axs[0, 0],
        spatial,
        'releasing_rois',
        'spatial_tau_px',
        r'spatial $\tau$ (px)',
        )
    axs[0, 0].set_title('releasing ROI anchor')

    plot_spatial_lines(
        axs[0, 1],
        spatial,
        'all_rois',
        'spatial_tau_px',
        r'spatial $\tau$ (px)',
        )
    axs[0, 1].set_title('all ROI anchor')
    axs[0, 1].legend(frameon=False, fontsize=6, loc='upper left')

    plot_spatial_lines(
        axs[1, 0],
        spatial,
        'releasing_rois',
        'spatial_tau_rel_0_1',
        r'spatial $\tau$ / 0-1 s',
        )
    axs[1, 0].axhline(1, color='0.3', linestyle='--', linewidth=0.8)

    plot_spatial_lines(
        axs[1, 1],
        spatial,
        'all_rois',
        'spatial_tau_rel_0_1',
        r'spatial $\tau$ / 0-1 s',
        )
    axs[1, 1].axhline(1, color='0.3', linestyle='--', linewidth=0.8)

    fig.tight_layout(w_pad=0.9, h_pad=1.0)
    for ext in ['.png', '.pdf']:
        out = save_stem / f'dopamine_norepinephrine_sensor_normalised_spatial_retention{ext}'
        fig.savefig(out, dpi=300, bbox_inches='tight')
        saved.append((ext[1:], out))
    plt.close(fig)

    print_statistics_section()
    for key in [
            'wholefield_peak_C_over_EC50',
            'roi_mean_C_over_EC50',
            'roi_tau_sensor_ratio',
            'wholefield_half_sensor_ratio',
            'late_early_ratio',
            ]:
        for row in summarise_by_transmitter(session, key):
            print(
                f'{row["sensor"]} {key}: '
                f'n={row["n"]}, median={row["median"]:.4g}, '
                f'mean={row["mean"]:.4g}, sem={row["sem"]:.4g}'
                )

    for source in ['releasing_rois', 'all_rois']:
        curr = spatial[
            (spatial['spatial_source'] == source)
            & (spatial['window'] == '3-4 s')
            ]
        for row in summarise_by_transmitter(curr, 'spatial_tau_rel_0_1'):
            print(
                f'{row["sensor"]} {source} spatial retention 3-4/0-1 s: '
                f'n={row["n"]}, median={row["median"]:.4g}, '
                f'mean={row["mean"]:.4g}'
                )
    print_files_saved(saved)

if __name__ == '__main__':
    main()
