# -*- coding: utf-8 -*-
'''
Created on 20 May 2026

plot LC-HPC axon peak-detection figures from the saved profile data

@author: Dinghao Luo
'''

#%% imports
from pathlib import Path
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')

repo_root = Path(__file__).resolve().parents[2]
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))

import peak_detection_functions as pdf
import project_paths as pp


#%% paths
payload_path = (
    pp.LCHPC_AXON_STEM
    / 'LCHPC_axon_GCaMP_all_profiles_plot_payload.npy'
    )


#%% main
payload = np.load(payload_path, allow_pickle=True).item()

for peak_payload in payload['peak_detection'].values():
    recname = peak_payload['recname']
    roiname = peak_payload['roiname']
    peak = peak_payload['peak']
    output_stem = (
        pp.LCHPC_AXON_FIGURES_STEM
        / 'peak_detection'
        / f'{recname} {roiname} {peak}'
        )
    output_stem.parent.mkdir(parents=True, exist_ok=True)
    pdf.plot_peak_v_shuf(
        roiname,
        peak_payload['mean_prof'],
        peak_payload['shuf_prof'],
        peak,
        peak_width=2,
        savepath=output_stem,
        samp_freq=peak_payload['samp_freq']
        )
