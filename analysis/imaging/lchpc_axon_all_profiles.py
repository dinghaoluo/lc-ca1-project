# -*- coding: utf-8 -*-
'''
Created on Fri Mar 21 12:41:07 2025

organise properties of individual ROIs into a profile dataframe
    for axon-GCaMP LC-CA1 recordings

@author: Dinghao Luo
'''

#%% imports
from pathlib import Path
import sys

import numpy as np
import pandas as pd
from scipy.stats import sem
from tqdm import tqdm

repo_root = Path(__file__).resolve().parents[2]
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))
import peak_detection_functions as pdf

from common_functions import (
    get_GPU_availability,
    mpl_formatting,
    smooth_convolve,
)
from console_formatting import print_files_saved, print_session
import project_paths as pp
mpl_formatting()


#%% paths and parameters
LC_axon_stem  = pp.LCHPC_AXON_STEM
all_sess_stem = LC_axon_stem / 'all_sessions'
plot_payload_path = LC_axon_stem / 'LCHPC_axon_GCaMP_all_profiles_plot_payload.npy'


#%% dataframe initialisation
fname = LC_axon_stem / 'LCHPC_axon_GCaMP_all_profiles.pkl'
df = pd.DataFrame({
    'recname': [],  # Axxx-202xxxxx-0x
    'roi_type': [],  # str, 'primary' or 'constituent'
    'constituents': [],  # list, filled only if roi_type=='primary'
    'coord': [],  # list of tuples (x, y)
    'size': [],  # pixel count of the ROI
    'run_onset_peak': [],
    'run_onset_peak_dFF': [],  # float
    'run_onset_peak_ch2': [],  # as control
    'mean_profile': [],  # mean fluorescence profile (dFF)
    'sem_profile': [],
    'mean_profile_ch2': [],
    'sem_profile_ch2': [],
    'var': [],  # trialwise variability
    })


#%% load paths to recordings
import rec_list
paths = rec_list.pathLCHPCGCaMP


#%% parameters
SMOOTHED = 1
SIGMA = 3  # frames

SAMP_FREQ = 30  # in Hz
BEF = 3
MAX_TIME = 10  # collect (for each trial) a maximum of 10 s of profile
MAX_SAMPLES = SAMP_FREQ * MAX_TIME


#%% GPU acceleration
_, GPU_AVAILABLE, _ = get_GPU_availability()
plot_payloads = {'peak_detection': {}}


#%% main
for path in paths:
    recname = Path(path).name
    print_session(recname)

    # load data
    proc_path = all_sess_stem / recname
    proc_data_path = proc_path / 'processed_data'
    input_paths = {
        'F_dFF': proc_path / f'{recname}_all_run.npy',
        'F_dFF_const': proc_data_path / 'RO_aligned_const_dict.npy',
        'F2_dFF': proc_path / f'{recname}_all_run_ch2.npy',
        'F2_dFF_const': proc_data_path / 'RO_aligned_const_ch2_dict.npy',
        'valid_ROIs': proc_data_path / 'valid_ROIs_dict.npy',
        'valid_coords': proc_data_path / 'valid_ROIs_coord_dict.npy',
        'const_coords': proc_data_path / 'constituent_ROIs_coord_dict.npy',
        }
    F_dFF = (
        np.load(input_paths['F_dFF'], allow_pickle=True).item()
        | np.load(input_paths['F_dFF_const'], allow_pickle=True).item()
        )
    F2_dFF = (
        np.load(input_paths['F2_dFF'], allow_pickle=True).item()
        | np.load(input_paths['F2_dFF_const'], allow_pickle=True).item()
        )

    valid_ROIs_dict = np.load(input_paths['valid_ROIs'], allow_pickle=True).item()
    valid_ROIs_coord_dict = np.load(input_paths['valid_coords'], allow_pickle=True).item()
    constituent_ROIs_coord_dict = np.load(input_paths['const_coords'], allow_pickle=True).item()

    # get list of clunames
    primary_rois = set(
        [int(name.split(' ')[1]) for name in [*valid_ROIs_dict]]
        )
    constituent_rois = set(
        [roi
         for roi_list in valid_ROIs_dict.values()
         for roi in roi_list]
        )
    all_rois = primary_rois | constituent_rois

    for roi in tqdm(all_rois,
                    desc='collecting ROI profiles'):
        roiname = f'ROI {roi}'
        # print(roiname)
        dFF = F_dFF[roiname]
        dFF2 = F2_dFF[roiname]

        if SMOOTHED:
            dFF = smooth_convolve(dFF, sigma=SIGMA, axis=1)
            dFF2 = smooth_convolve(dFF2, sigma=SIGMA, axis=1)

        # get mean and sem spiking profiles
        # returned in lists for parqueting
        profile_dFF = np.array([trace[:MAX_SAMPLES] for trace in dFF])
        profile_dFF2 = np.array([trace[:MAX_SAMPLES] for trace in dFF2])
        mean_profile = list(np.mean(profile_dFF, axis=0))
        sem_profile = list(sem(profile_dFF, axis=0))
        mean_profile_ch2 = list(np.mean(profile_dFF2, axis=0))
        sem_profile_ch2 = list(sem(profile_dFF2, axis=0))

        # identity ('primary' or 'constituent')
        identity = 'primary' if roi in primary_rois else 'constituent'

        # constituent list
        if identity=='primary':
            constituents = valid_ROIs_dict[roiname]
        else:
            constituents = None

        # coord
        roi_idx = int(roiname.split(' ')[-1])
        coord = (
            valid_ROIs_coord_dict[roiname]
            if roi_idx in primary_rois
            else constituent_ROIs_coord_dict[roiname]
        )

        # run-onset peak detection
        peak, mean_prof, shuf_prof = pdf.peak_detection(
            dFF,
            around=4,  # check baseline on a higher threshold
            peak_width=2,  # 2 seconds to include some of the slightly offset peaks
            min_peak=.5,
            samp_freq=SAMP_FREQ,
            centre_bin=SAMP_FREQ*BEF,
            bootstrap=500,
            no_boundary=True,
            GPU_AVAILABLE=GPU_AVAILABLE,
            VERBOSE=False
            )
        if identity=='primary':
            plot_payloads['peak_detection'][f'{recname} {roiname}'] = {
                'recname': recname,
                'roiname': roiname,
                'mean_prof': mean_prof,
                'shuf_prof': shuf_prof,
                'peak': peak,
                'samp_freq': SAMP_FREQ
                }

        peak_ch2, _, _ = pdf.peak_detection(
            dFF2,
            around=4,  # check baseline on a higher threshold
            peak_width=2,  # 2 seconds to include some of the slightly offset peaks
            samp_freq=SAMP_FREQ,
            centre_bin=SAMP_FREQ*BEF,
            bootstrap=500,
            no_boundary=True,
            GPU_AVAILABLE=GPU_AVAILABLE,
            VERBOSE=False
            )

        # trial-wise variability
        max_length = max(len(trace) for trace in dFF)
        traces = [np.asarray(trace[:max_length], dtype=float) for trace in dFF]
        num_trials = len(traces)
        corr_matrix = np.full((num_trials, num_trials), np.nan)
        for i in range(num_trials):
            for j in range(i + 1, num_trials):
                corr_matrix[i, j] = np.corrcoef(traces[i], traces[j])[0, 1]
        corr_values = corr_matrix[np.triu_indices(num_trials, k=1)]
        variability = 1 - np.nanmedian(corr_values)

        # full ROI name
        full_roiname = f'{recname} {roiname}'

        # put into dataframe
        df.loc[full_roiname] = np.array(
            [recname,
             identity,
             constituents,
             coord,
             len(coord[0]),  # size of ROI
             peak,
             mean_prof[peak] if peak is not None else np.nan,
             peak_ch2,
             mean_profile,
             sem_profile,
             mean_profile_ch2,
             sem_profile_ch2,
             variability],
            dtype='object'
            )

## save dataframe
df.to_pickle(fname)
np.save(plot_payload_path, plot_payloads)
print_files_saved([
    ('dataframe', fname),
    ('peak-detection plotting data', plot_payload_path),
])
