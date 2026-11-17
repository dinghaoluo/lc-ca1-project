# -*- coding: utf-8 -*-
'''
Created on Wed Nov 11 15:20:24 2026

early and late trial-count sensitivity and single-trial amplitude variance

@author: Dinghao Luo
'''

#%% imports
import sys
from pathlib import Path

repo_root = Path(__file__).resolve().parents[3]
if str(repo_root) not in sys.path:
    sys.path.insert(0, str(repo_root))
if str(repo_root / 'utils') not in sys.path:
    sys.path.insert(0, str(repo_root / 'utils'))

import project_paths as pp
import rec_list

import numpy as np
import pandas as pd
from threadpoolctl import threadpool_limits

from analyse_acr_consistency_aligned import (load_recording, event_frames,
                                              recnames)
from grid_response_test import quantify_grid_response

repo_root   = Path(__file__).resolve().parents[3]
OUTPUT_ROOT = pp.DLIGHT_HOTSPOTS_OUTPUT_ROOT
N_SHUFFLE   = 1000
SEED        = 42


#%% variance across kept trials
def amplitude_variance():
    rows = []
    for recname in list(recnames):
        behaviour, profile, traces = load_recording(recname)
        frames = event_frames(behaviour)
        result = quantify_grid_response(traces, frames, n_shuffle=N_SHUFFLE,
                                        seed=SEED, return_single_trial=True)
        keep = result['keep']
        kept = keep.sum(axis=1)
        up = profile['Up'].to_numpy(dtype=bool)
        observed_variance = np.full(len(profile), np.nan)
        shuffle_variance = np.full(len(profile), np.nan)
        for acr in range(len(profile)):
            if kept[acr] < 2:
                continue
            observed = result['trial_amplitude'][acr, keep[acr]]
            shuffled = result['shuffle_single_trial_amplitude'][acr, keep[acr], :]
            observed_variance[acr] = np.var(observed, ddof=1)
            shuffle_variance[acr] = np.median(np.var(shuffled, axis=0, ddof=1))
        ratio = np.divide(observed_variance, shuffle_variance,
                          out=np.full(len(profile), np.nan),
                          where=np.isfinite(shuffle_variance) & (shuffle_variance > 0))
        for label, group in (('Up', up), ('other', ~up)):
            finite = group & np.isfinite(ratio)
            median_ratio = float(np.median(ratio[finite])) if finite.any() else np.nan
            rows.append(dict(rec_id=recname, animal=recname.split('-')[0], group=label,
                             n_trials=len(frames), n_acr=int(group.sum()), n_ratio=int(finite.sum()),
                             n_fewer_two=int(np.sum(group & (kept < 2))),
                             n_zero_observed=int(np.sum(group & (observed_variance == 0))),
                             n_invalid_shuffle=int(np.sum(group & (kept >= 2)
                                                          & (~np.isfinite(shuffle_variance) | (shuffle_variance <= 0)))),
                             median_ratio=median_ratio))
        first, second = rows[-2], rows[-1]
        print(f'{recname}: {len(frames)} valid trials; Up median variance ratio={first["median_ratio"]:.6g} ({first["n_ratio"]} ACRs), other={second["median_ratio"]:.6g} ({second["n_ratio"]} ACRs)', flush=True)
        del result
    return pd.DataFrame(rows)


#%% run
if __name__ == '__main__':
    with threadpool_limits(limits=1):
        amplitude_variance()
