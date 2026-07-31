# -*- coding: utf-8 -*-
'''
Created on Sat Aug 10 12:24:46 2024

fibre masks, Suite2p processing and imaging ROI measurements

@author: Dinghao Luo
'''

#%% imports
import numpy as np
import matplotlib.pyplot as plt

from scipy.ndimage import gaussian_filter, uniform_filter, minimum_filter
from skimage.morphology import remove_small_objects, remove_small_holes
from skimage.filters import threshold_local

import project_paths as pp


#%% filter pixels for dLight expression level
def generate_adaptive_membrane_mask(
    mean_img,
    gaussian_sigma=1.5,
    peak_min_distance=5,
    adaptive_block_size=21,  # odd window for uneven fields
    valley_radius=10,  # local soma-valley scale
    uniformity_thresh=3.5,            # higher -> stricter "cell-like"
    z_tau=3.8,
    min_region_size=150,
    hole_size_threshold=150,
    visualize=1
):
    '''
    combine cell-like membrane regions with locally bright neuropil

    with visualisation enabled, also return the broad mask and diagnostic figure
    '''
    img = mean_img.astype(np.float32)
    img = (img - np.percentile(img, 2)) / (np.percentile(img, 98) - np.percentile(img, 2) + 1e-6)
    img = np.clip(img, 0, 1)

    img_s = gaussian_filter(img, sigma=gaussian_sigma)
    # start broad because this mask contains both cell-like and neuropil-like signal
    base_mask = img_s > threshold_local(img_s, block_size=adaptive_block_size, method='gaussian')

    # local contrast separates uniform neuropil from dark soma regions
    r = valley_radius
    mu = uniform_filter(img_s, size=2*r+1, mode='reflect')
    mu2 = uniform_filter(img_s**2, size=2*r+1, mode='reflect')
    sigma = np.sqrt(np.maximum(mu2 - mu**2, 0.0))

    # local minimum to capture dark somas ("valleys")
    loc_min = minimum_filter(img_s, size=2*r+1, mode='reflect')
    # valley-contrast score: big when bright rims surround a dark soma
    uniformity = mu / (sigma + 1e-6)
    neuropil_like = (uniformity > uniformity_thresh)
    cell_like = ~neuropil_like

    # keep the broad local-threshold route in cell-like neighbourhoods
    cell_regions_mask = base_mask & cell_like
    cell_membrane_mask = cell_regions_mask

    # uniform neuropil needs a separate intensity/local-std threshold
    z = img_s / (sigma + 1e-6)
    neuropil_mask = (z > z_tau) & neuropil_like

    # avoid joining nearby structures during mask cleanup
    neuropil_mask = remove_small_holes(neuropil_mask, area_threshold=hole_size_threshold)
    neuropil_mask = remove_small_objects(neuropil_mask, min_size=min_region_size)

    final_mask = (cell_membrane_mask | neuropil_mask)
    final_mask = remove_small_objects(final_mask, min_size=min_region_size)

    debug = dict(img_s=img_s, mu=mu, sigma=sigma, loc_min=loc_min,
                 uniformity=uniformity, neuropil_like=neuropil_like,
                 cell_like=cell_like, z=z,
                 cell_membrane_mask=cell_membrane_mask, neuropil_mask=neuropil_mask)

    if visualize:
        neuropil_mask = neuropil_like.astype(bool)
        neuropil_zmap = debug['z']
        neuropil_zmap[~neuropil_mask] = np.nan
        fig, ax = plt.subplots(2, 3, figsize=(12, 8))
        ax[0,0].imshow(img, cmap='gray'); ax[0,0].set_title('Original'); ax[0,0].axis('off')
        ax[0,1].imshow(uniformity, cmap='magma'); ax[0,1].set_title('uniformity'); ax[0,1].axis('off')
        ax[0,2].imshow(neuropil_like, cmap='coolwarm'); ax[0,2].set_title('Neuropil-like'); ax[0,2].axis('off')
        ax[1,0].imshow(cell_membrane_mask, cmap='gray'); ax[1,0].set_title('Cell membranes'); ax[1,0].axis('off')
        ax[1,1].imshow(neuropil_zmap, cmap='viridis', vmax=np.nanpercentile(neuropil_zmap, 99));ax[1,1].set_title('Neuropil_Z_map'); ax[1,1].axis('off')
        ax[1,2].imshow(img, cmap='gray'); ax[1,2].imshow(final_mask, alpha=0.35, cmap='Reds')
        ax[1,2].set_title('Final'); ax[1,2].axis('off')
        plt.tight_layout()
        # plt.show()
        return base_mask, final_mask, fig

    else:
        return final_mask


#%% suite2p registration and roi extraction
DEFAULT_SUITE2P_REGISTRATION_OPS_PATH = pp.TWO_PHOTON_ROOT / 'registration_parameters.npy'

def run_suite2p_registration(path, ops_path=DEFAULT_SUITE2P_REGISTRATION_OPS_PATH):
    '''
    run Wang-lab Suite2p registration on one recording directory.

    parameters:
    - path: recording directory containing the input tiffs. output is written
      to raw imaging storage only when that session has no raw-side registration;
      new outputs are also mirrored to data/suite2p_registration/<recording>/
    - ops_path: path to the default registration parameter file
    '''
    from contextlib import redirect_stdout
    from datetime import timedelta
    from pathlib import Path
    from time import time
    import shutil
    import suite2p
    import sys

    path = Path(path)
    ops = np.load(ops_path, allow_pickle=True).item()

    ops['input_format'] = 'tif'
    ops['sparse_mode'] = True
    ops['roidetect'] = False
    ops['reg_tif'] = True
    ops['reg_tif_chan2'] = True
    ops['nonrigid'] = False

    fallback_stem = pp.SUITE2P_REGISTRATION_STEM / path.name
    raw_registration_exists = (
        (path / 'suite2p' / 'plane0' / 'ops.npy').exists()
        or (path / 'processed' / 'suite2p' / 'plane0' / 'ops.npy').exists()
        )
    if path.exists() and not raw_registration_exists:
        write_targets = [path, fallback_stem]
    else:
        write_targets = [fallback_stem]
    outdir = write_targets[0]
    outdir.mkdir(parents=True, exist_ok=True)
    ops['save_path0'] = str(outdir)

    print('registration starts')
    t0 = time()

    suite2p_dir = outdir / 'suite2p'
    suite2p_dir.mkdir(parents=True, exist_ok=True)
    pathlog = suite2p_dir / 'run_suite2p.log'

    db = {
        'data_path': [str(path)],
        'save_path0': str(outdir)
        }

    class Tee:
        def __init__(self, *streams):
            self.streams = streams

        def write(self, text):
            for stream in self.streams:
                stream.write(text)
                stream.flush()

        def flush(self):
            for stream in self.streams:
                stream.flush()

    with open(pathlog, 'w') as file:
        with redirect_stdout(Tee(sys.stdout, file)):
            print(f'running suite2p v{suite2p.version} from Spyder')
            suite2p.run_s2p(ops=ops, db=db)

    print(f'registration complete ({str(timedelta(seconds=int(time()-t0)))})\n')

    for mirror_root in write_targets[1:]:
        mirror_suite2p = mirror_root / 'suite2p'
        mirror_root.mkdir(parents=True, exist_ok=True)
        if mirror_suite2p.exists():
            shutil.rmtree(mirror_suite2p)
        shutil.copytree(outdir / 'suite2p', mirror_suite2p)
        print(f'registration mirrored to {mirror_suite2p}\n')

def run_suite2p_roi_extraction(path):
    '''
    run Wang-lab Suite2p ROI extraction after registration for one recording directory.

    parameters:
    - path: recording directory whose registered output can be resolved by project_paths
    '''
    from contextlib import redirect_stdout
    from datetime import timedelta
    from pathlib import Path
    from time import time
    import suite2p

    path = Path(path)
    register_path = path / 'processed'
    output_path = register_path

    ops = np.load(register_path / 'suite2p' / 'plane0' / 'ops.npy', allow_pickle=True).item()
    ops['sparse_mode'] = True
    ops['anatomical_only'] = False
    ops['roidetect'] = True
    ops['spatial_scale'] = 2
    ops['denoise'] = True

    ops['circular_neuropil'] = True
    ops['inner_neuropil_radius'] = 5

    ops['max_iterations'] = 1
    ops['high_pass'] = 200

    ops['wang:bin_size'] = 1
    ops['wang:high_pass_overlapping'] = True
    ops['wang:rolling_width'] = 30
    ops['wang:rolling_bin'] = 'max'
    ops['wang:use_alt_norm'] = True
    ops['wang:downsample_scale'] = 1
    ops['wang:thresh_act_pix'] = 0.04
    ops['wang:thresh_peak_default'] = 0.03

    ops['wang:save_roi_iterations'] = True
    ops['save_path_new'] = str(output_path)

    ops['wang:save_path_sparsedetect'] = str(output_path)
    ops['wang:neuropil_lam'] = True
    ops['wang:movie_chunk'] = 10000
    ops['wang:norm_method'] = 'max'

    output_path.mkdir(parents=True, exist_ok=True)
    log_path = output_path / 'run_suite2p-wang-lab.log'
    print(f'running log: {log_path}')

    db = {
        'data_path': [str(path)],
        'save_path0': str(output_path),
    }

    print('roi extraction starts...')
    t0 = time()
    with open(log_path, 'w') as file:
        with redirect_stdout(file):
            print(f'running suite2p v{suite2p.version} from Spyder')
            suite2p.run_s2p(ops=ops, db=db)

    print(f'roi extraction complete ({str(timedelta(seconds=int(time()-t0)))})\n')


#%% ROI and neuropil measurements
def build_roi_mask(roi_dict, shape=(512, 512)):
    '''
    build a combined boolean mask from ROI pixel coordinates.
    '''
    mask = np.zeros(shape, dtype=bool)
    for roi in roi_dict.values():
        mask[roi['ypix'], roi['xpix']] = True
    return mask

def identify_releasing_rois(pixel_ri_stim, roi_dict, alpha=0.05, min_ri=0.1):
    '''
    identify ROIs whose median stimulus response is significantly above zero.
    '''
    from scipy.stats import wilcoxon

    releasing = {}
    for roi_id, roi in roi_dict.items():
        values = pixel_ri_stim[roi['ypix'], roi['xpix'], :]
        means = np.nanmean(values, axis=0)
        means = [mean for mean in means if np.isfinite(mean)]
        if len(means) > 2:
            _, p_value = wilcoxon(means, alternative='greater')
            if p_value < alpha and np.mean(means) > min_ri:
                releasing[roi_id] = roi
    return releasing
