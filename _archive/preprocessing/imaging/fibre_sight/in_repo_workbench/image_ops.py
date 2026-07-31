'''
Created on 5 April 2026
normalisation, crop sampling, and augmentation for channel-2 images

@author: Dinghao Luo
'''

#%% imports
import numpy as np


#%% normalisation
def robust_normalise(image, low=1, high=99.7):
    image = np.asarray(image, dtype=np.float32)
    finite = np.isfinite(image)

    if not np.any(finite):
        return np.zeros_like(image, dtype=np.float32)

    # Scale each reference on its own percentiles; acquisition brightness should not become a label.
    lo, hi = np.nanpercentile(image[finite], [low, high])
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        lo = np.nanmin(image[finite])
        hi = np.nanmax(image[finite])

    if hi <= lo:
        return np.zeros_like(image, dtype=np.float32)

    image = (image - lo) / (hi - lo)
    return np.clip(image, 0, 1).astype(np.float32)


#%% crops
def sample_crop_bounds(mask, patch_size, rng, foreground_fraction=0.75):
    height, width = mask.shape
    patch_size = int(patch_size)

    if patch_size >= height or patch_size >= width:
        return 0, 0, min(patch_size, height), min(patch_size, width)

    # keep some background-only patches instead of teaching the model that
    # every crop contains an axon
    use_foreground = rng.random() < foreground_fraction and np.any(mask)
    if use_foreground:
        ypix, xpix = np.where(mask)
        idx = rng.integers(0, len(ypix))
        centre_y = int(ypix[idx])
        centre_x = int(xpix[idx])
        top = centre_y - patch_size // 2
        left = centre_x - patch_size // 2
    else:
        top = int(rng.integers(0, height - patch_size + 1))
        left = int(rng.integers(0, width - patch_size + 1))

    top = int(np.clip(top, 0, height - patch_size))
    left = int(np.clip(left, 0, width - patch_size))
    return top, left, patch_size, patch_size


def crop_pair(image, mask, bounds):
    top, left, height, width = bounds
    image_crop = image[top:top + height, left:left + width]
    mask_crop = mask[top:top + height, left:left + width]
    return image_crop, mask_crop


#%% augmentation
def augment_pair(image, mask, rng, intensity_jitter=0.15, noise_sd=0.02):
    if rng.random() < 0.5:
        image = np.flip(image, axis=0)
        mask = np.flip(mask, axis=0)

    if rng.random() < 0.5:
        image = np.flip(image, axis=1)
        mask = np.flip(mask, axis=1)

    k = int(rng.integers(0, 4))
    if k:
        image = np.rot90(image, k)
        mask = np.rot90(mask, k)

    if intensity_jitter > 0:
        scale = 1 + rng.uniform(-intensity_jitter, intensity_jitter)
        offset = rng.uniform(-intensity_jitter, intensity_jitter) * 0.25
        image = image * scale + offset

    if noise_sd > 0:
        image = image + rng.normal(0, noise_sd, size=image.shape)

    image = np.clip(image, 0, 1).astype(np.float32)
    mask = mask.astype(np.float32)
    return image.copy(), mask.copy()
