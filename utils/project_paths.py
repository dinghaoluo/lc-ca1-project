# -*- coding: utf-8 -*-
'''
Created on 29 June 2026

project path stems for the LC-CA1 repository

Raw recording roots can still point to lab storage, while generated data and
figures resolve inside this repository. Edit the raw-input roots below when
using downloaded data outside the MPFI network-drive layout.

@author: Dinghao Luo
'''

#%% imports
from pathlib import Path

#%% locate the repository root
_THIS_DIR = Path(__file__).resolve().parent
_REPO_ROOT = _THIS_DIR.parent


#%% in-repo generated-output roots
REPO_ROOT = _REPO_ROOT
DATA_ROOT = _REPO_ROOT / 'data'
FIGURES_ROOT = _REPO_ROOT / 'figures'
OUTPUT_ROOT = _REPO_ROOT / 'outputs'
LOGS_ROOT = OUTPUT_ROOT / 'logs'

# Archive name retained for scripts that use ARCHIVE_ROOT as an
# analysis-ready-data stem. It points in-repo for full pipeline runs.
ARCHIVE_ROOT = DATA_ROOT
DINGHAO_ROOT = DATA_ROOT

# Edit these raw-input roots if the downloaded data are stored elsewhere.
MICEEXP_ROOT = Path(r'Z:/Dinghao/MiceExp')
RAPHAEL_ROOT = Path(r'Z:/Raphael_tests/mice_expdata')
BEHAV_DATA_ANALYSIS_ROOT = Path(r'Z:/Dinghao/Behav/DataAnalysis')
PUPIL_TRACKING_INPUT_STEM = Path(r'Z:/Dinghao/code_dinghao/pupil_tracking')
TWO_PHOTON_ROOT = Path(r'Z:/Dinghao/2p_recording')


#%% generated-data stems
BEHAVIOUR_STEM = DATA_ROOT / 'behaviour'
BEHAVIOUR_EXPERIMENTS_STEM = BEHAVIOUR_STEM / 'all_experiments'
BEHAVIOUR_SESSION_PROFILES_STEM = BEHAVIOUR_STEM / 'session_profiles'
BEHAVIOUR_EXPERIMENT_FOLDER_NAMES = {
    'LC': 'lc',
    'LCterm': 'lc_term',
    'HPC': 'hpc',
    'HPCRaphi': 'hpc_raphi',
    'HPCLC': 'hpc_lc',
    'HPCLCterm': 'hpc_lc_term',
    'HPCGRABNE': 'hpc_grabne',
    'LCHPCGCaMP': 'lc_hpc_gcamp',
    'LCHPCGCaMPImmobile': 'lc_hpc_gcamp_immobile',
    'HPCdLightLCOpto': 'hpc_dlight_lc_opto',
    'HPCdLightLCOptoInh': 'hpc_dlight_lc_opto_inh',
    'HPCnLightLCOpto': 'hpc_nlight_lc_opto',
    }

LC_EPHYS_STEM = DATA_ROOT / 'lc_ephys'
HPC_EPHYS_STEM = DATA_ROOT / 'hpc_ephys'
HPC_ALL_STEM = HPC_EPHYS_STEM / 'all_cells'
LC_OPTO_EPHYS_STEM = (
    LC_EPHYS_STEM / 'archive' / 'lc_tagging_optogenetics_behaviour'
    )
PHARMACOLOGY_STEM = DATA_ROOT / 'pharmacology'
PHARMACOLOGY_FIGURES_STEM = FIGURES_ROOT / 'pharmacology'
LCHPC_AXON_ROOT = DATA_ROOT / 'lchpc_axon_gcamp'
LCHPC_AXON_STEM = LCHPC_AXON_ROOT / 'active'
HPC_DLIGHT_STEM = DATA_ROOT / 'hpc_dlight'
HPC_DLIGHT_LC_OPTO_STEM = HPC_DLIGHT_STEM / 'lc_opto'
HPC_NLIGHT_STEM = DATA_ROOT / 'hpc_nlight'
HPC_NLIGHT_LC_OPTO_STEM = HPC_NLIGHT_STEM / 'lc_opto'
HPC_GRABNE_STEM = DATA_ROOT / 'hpc_grabne'
HPC_GRABNE_LC_OPTO_STEM = HPC_GRABNE_STEM / 'lc_opto'
BEHAVIOUR_LC_OPTO_STEM = BEHAVIOUR_STEM / 'lc_opto'
GRABNE_STEM = HPC_GRABNE_STEM / 'tone'
LC_FIGURES_STEM = FIGURES_ROOT / 'lc_ephys' / 'behaviour_examples'
RUN_BOUTS_STEM = DATA_ROOT / 'run_bouts'
LC_ALL_STEM = LC_EPHYS_STEM / 'all_cells'
LC_ALL_TAGGED_STEM = LC_EPHYS_STEM / 'tagging'
LC_BY_SESS_STEM = LC_EPHYS_STEM / 'by_session'
LC_TAGGED_BY_SESS_STEM = LC_EPHYS_STEM / 'tagged_by_session'
LCHPC_AXON_IMMOBILE_STEM = LCHPC_AXON_ROOT / 'immobile'
PUPIL_TRACKING_OUTPUT_STEM = BEHAVIOUR_STEM / 'pupil_tracking'
MATLAB_SPIKE_PIPELINE_STEM = DATA_ROOT / 'matlab_spike_pipeline'
SUITE2P_REGISTRATION_STEM = DATA_ROOT / 'suite2p_registration'


#%% generated-figure stems
BEHAVIOUR_FIGURES_STEM = FIGURES_ROOT / 'behaviour'
BEHAVIOUR_EXPERIMENTS_FIGURES_STEM = BEHAVIOUR_FIGURES_STEM / 'all_experiments'
BEHAVIOUR_SESSION_PROFILES_FIGURES_STEM = BEHAVIOUR_FIGURES_STEM / 'session_profiles'

LC_EPHYS_FIGURES_STEM = FIGURES_ROOT / 'lc_ephys'
HPC_EPHYS_FIGURES_STEM = FIGURES_ROOT / 'hpc_ephys'
HPC_ALL_FIGURES_STEM = HPC_EPHYS_FIGURES_STEM / 'all_cells'
LC_OPTO_EPHYS_FIGURES_STEM = (
    LC_EPHYS_FIGURES_STEM / 'archive' / 'lc_tagging_optogenetics_behaviour'
    )
HPC_OPTO_EPHYS_FIGURES_STEM = (
    HPC_EPHYS_FIGURES_STEM
    / 'archive'
    / 'lc_terminal_optogenetics_behaviour'
    )
EPHYS_HPC_LC_OPTO_FIGURES_STEM = HPC_EPHYS_FIGURES_STEM / 'lc_stim'
LC_STIM_CA1_EFFECT_FIGURES_STEM = (
    EPHYS_HPC_LC_OPTO_FIGURES_STEM / 'ca1_effect'
    )
LCHPC_AXON_FIGURES_ROOT = FIGURES_ROOT / 'lchpc_axon_gcamp'
LCHPC_AXON_FIGURES_STEM = LCHPC_AXON_FIGURES_ROOT / 'active'
HPC_DLIGHT_FIGURES_STEM = FIGURES_ROOT / 'hpc_dlight'
HPC_DLIGHT_LC_OPTO_FIGURES_STEM = HPC_DLIGHT_FIGURES_STEM / 'lc_opto'
HPC_DLIGHT_LC_OPTO_ALL_SESSIONS_FIGURES_STEM = (
    HPC_DLIGHT_LC_OPTO_FIGURES_STEM / 'all_sessions'
    )
HPC_DLIGHT_LC_OPTO_RELEASE_SITE_HETEROGENEITY_FIGURES_STEM = (
    HPC_DLIGHT_LC_OPTO_FIGURES_STEM / 'release_site_heterogeneity'
    )
HPC_NLIGHT_FIGURES_STEM = FIGURES_ROOT / 'hpc_nlight'
HPC_NLIGHT_LC_OPTO_FIGURES_STEM = HPC_NLIGHT_FIGURES_STEM / 'lc_opto'
HPC_NLIGHT_LC_OPTO_ALL_SESSIONS_FIGURES_STEM = (
    HPC_NLIGHT_LC_OPTO_FIGURES_STEM / 'all_sessions'
    )
HPC_NLIGHT_LC_OPTO_RELEASE_SITE_HETEROGENEITY_FIGURES_STEM = (
    HPC_NLIGHT_LC_OPTO_FIGURES_STEM / 'release_site_heterogeneity'
    )
HPC_GRABNE_FIGURES_STEM = FIGURES_ROOT / 'hpc_grabne'
HPC_GRABNE_LC_OPTO_FIGURES_STEM = HPC_GRABNE_FIGURES_STEM / 'lc_opto'
BEHAVIOUR_LC_OPTO_FIGURES_STEM = BEHAVIOUR_FIGURES_STEM / 'lc_opto'
GRABNE_FIGURES_STEM = HPC_GRABNE_FIGURES_STEM / 'tone'
RUN_BOUTS_FIGURES_STEM = FIGURES_ROOT / 'run_bouts'
LC_ALL_FIGURES_STEM = LC_EPHYS_FIGURES_STEM / 'all_cells'
LC_ALL_TAGGED_FIGURES_STEM = LC_EPHYS_FIGURES_STEM / 'tagging'
LC_BY_SESS_FIGURES_STEM = LC_EPHYS_FIGURES_STEM / 'by_session'
LC_TAGGED_BY_SESS_FIGURES_STEM = LC_EPHYS_FIGURES_STEM / 'tagged_by_session'
LCHPC_AXON_IMMOBILE_FIGURES_STEM = LCHPC_AXON_FIGURES_ROOT / 'immobile'
PUPIL_TRACKING_FIGURES_STEM = BEHAVIOUR_FIGURES_STEM / 'pupil_tracking'
MOVIE_FIGURES_STEM = FIGURES_ROOT / 'movies'
IMAGING_REF_CHANNEL_MAPS_FIGURES_STEM = FIGURES_ROOT / 'imaging_ref_channel_maps'


#%% path lookup
def behaviour_experiment_stem(dataset_name):
    return BEHAVIOUR_EXPERIMENTS_STEM / BEHAVIOUR_EXPERIMENT_FOLDER_NAMES[dataset_name]

def resolve_suite2p_session_stem(raw_session_path):
    '''Return raw Suite2p output if present, else repo-local output.

    Existing raw-side registrations remain the primary source on this machine.
    The in-repo data/suite2p_registration mirror is only used when no raw-side
    Suite2p result is present.
    '''
    raw_session_path = Path(raw_session_path)
    raw_suite2p_stem = raw_session_path / 'suite2p'
    if (raw_suite2p_stem / 'plane0' / 'ops.npy').exists():
        return raw_suite2p_stem

    archive_raw_stem = raw_session_path / 'processed' / 'suite2p'
    if (archive_raw_stem / 'plane0' / 'ops.npy').exists():
        return archive_raw_stem

    fallback_stem = SUITE2P_REGISTRATION_STEM / raw_session_path.name / 'suite2p'
    if (fallback_stem / 'plane0' / 'ops.npy').exists():
        return fallback_stem

    return raw_suite2p_stem

def resolve_matlab_pipeline_file(raw_file, recname):
    '''Return raw MATLAB output if present, else repo-local mirrored output.

    The lab/raw recording folder remains the primary source on this machine.
    The in-repo data/matlab_spike_pipeline mirror is only used when the raw
    MATLAB product is missing.
    '''
    raw_file = Path(raw_file)
    if raw_file.exists():
        return raw_file

    fallback_file = MATLAB_SPIKE_PIPELINE_STEM / recname / raw_file.name
    if fallback_file.exists():
        return fallback_file

    return raw_file
