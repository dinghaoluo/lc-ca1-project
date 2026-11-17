# -*- coding: utf-8 -*-
'''
Created on Tue Nov 17 09:54:17 2026

plot the median saved Up template for each Dbh dLight animal

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

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter1d
from threadpoolctl import threadpool_limits
from common_functions import mpl_formatting

from check_template_variants import (load_recording, event_frames,
                                      align_trials, keep_aligned_trials,
                                      normalised_profiles, TIME)

mpl_formatting()


#%% templates and figure
def main():
    profiles_by_animal = {}
    recording_counts = {}
    saved_counts = {}
    for recname in rec_list.recDLightDbh:
        animal = recname.split('-')[0]
        behaviour, profile, dff = load_recording(recname)
        run_onset_frames = event_frames(behaviour)
        aligned, _ = keep_aligned_trials(align_trials(dff, run_onset_frames))
        up = profile['Up'].to_numpy(dtype=bool)
        profiles = normalised_profiles(np.nanmean(aligned[up], axis=1))
        profiles_by_animal.setdefault(animal, []).append(profiles)
        recording_counts[animal] = recording_counts.get(animal, 0) + 1
        saved_counts[animal] = saved_counts.get(animal, 0) + int(up.sum())
        print(f'{recname}: {len(profiles)} saved Up profiles', flush=True)
    fig, axes = plt.subplots(2, 2, figsize=(8, 5), sharex=True, sharey=True,
                              layout='constrained')
    with (pp.DLIGHT_HOTSPOTS_OUTPUT_ROOT / 'PROGRESS.md').open('a', encoding='utf-8') as log:
        log.write('\n## 17 November 2026: animal templates\n\n'
                  'For each animal I used every recording and its saved Up ACRs, '
                  'averaged all valid kept trials per ACR, subtracted the frames 30:60 mean '
                  'and divided by the frames 60:180 L2 norm with normalised_profiles unchanged. '
                  'The unsmoothed curve is the median of all contributing ACR profiles; '
                  'the smoothed curve uses sigma=3, mode=nearest, without renormalisation.\n\n')
        for axis, animal in zip(axes.flat, sorted(profiles_by_animal)):
            profiles = np.concatenate(profiles_by_animal[animal], axis=0)
            median = np.median(profiles, axis=0)
            smoothed = gaussian_filter1d(median, sigma=3, mode='nearest')
            axis.plot(TIME, median, color='0.55', linewidth=.6,
                       label='Unsmoothed median')
            axis.plot(TIME, smoothed, color='#228b22', linewidth=1.5,
                       label=f'{animal}, smoothed')
            axis.set(xlim=(-2, 4), xlabel='Time from run onset (s)',
                      ylabel='Normalised dF/F')
            axis.spines[['top', 'right']].set_visible(False)
            axis.legend(frameon=False, fontsize=8)
            log.write(f'{animal}: {recording_counts[animal]} recordings, '
                      f'{len(profiles)} contributing ACRs / {saved_counts[animal]} saved Up ACRs.\n')
        output = pp.DLIGHT_HOTSPOTS_FIGURE_ROOT / 'per_animal_templates.png'
        fig.savefig(output, dpi=200, facecolor='white')
        plt.close(fig)
        log.write('\nFigure: figures/dlight_hotspots/per_animal_templates.png; '
                  'four panels with shared axes, axis titles and legends only.\n')
    print(f'Saved {output}', flush=True)


#%% run
if __name__ == '__main__':
    with threadpool_limits(limits=1):
        main()
