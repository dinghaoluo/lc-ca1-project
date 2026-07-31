'''
Created on 8 April 2026
Modified on 23 July 2026 to repair the Dataset import and keep target construction here
sample channel-2 image patches and their curated ROI masks

@author: Dinghao Luo
'''

#%% imports
import numpy as np
from scipy import ndimage as ndi

from .image_ops import augment_pair, crop_pair, robust_normalise, sample_crop_bounds
from .manifest import read_manifest
from .roi_io import load_roi_dict, roi_dict_to_mask

try:
    import torch
    from torch.utils.data import Dataset
except ImportError:
    torch = None
    Dataset = object


#%% targets
def make_target(mask, mode='foreground', support_radius=3):
    mask = np.asarray(mask).astype(bool)

    if mode == 'foreground':
        return mask[None, ...].astype(np.float32)

    if mode == 'foreground_support':
        # The wider channel was my recall-heavy trial; saved ROIs still use
        # foreground channel 0.
        support = ndi.binary_dilation(mask, iterations=int(support_radius))
        target = np.stack([mask, support], axis=0)
        return target.astype(np.float32)

    raise ValueError(f'unknown target mode: {mode}')


#%% dataset
class AxonROIDataset(Dataset):
    def __init__(
            self,
            manifest_path,
            split='train',
            patch_size=256,
            patches_per_image=24,
            foreground_fraction=0.75,
            normalise_percentiles=(1, 99.7),
            augment=True,
            cache_images=False,
            target_mode='foreground',
            support_radius=3,
            seed=7,
            ):
        if torch is None:
            raise ImportError('please install pytorch to use AxonROIDataset')

        self.rows = read_manifest(manifest_path, included_only=True, split=split)
        if len(self.rows) == 0:
            raise ValueError(f'no rows found for split {split}')

        self.split = split
        self.patch_size = int(patch_size)
        self.patches_per_image = int(patches_per_image)
        self.foreground_fraction = float(foreground_fraction)
        self.normalise_percentiles = tuple(normalise_percentiles)
        self.augment = bool(augment)
        self.cache_images = bool(cache_images)
        self.target_mode = target_mode
        self.support_radius = int(support_radius)
        self.seed = int(seed)
        self.epoch = 0
        self.cache = {}

    def __len__(self):
        return len(self.rows) * self.patches_per_image

    def set_epoch(self, epoch):
        self.epoch = int(epoch)

    def __getitem__(self, idx):
        row = self.rows[idx % len(self.rows)]
        image, mask = self.load_pair(row)

        # the epoch changes the crops while the seed keeps a rerun reproducible
        rng = np.random.default_rng(self.seed + self.epoch * len(self) + idx)
        bounds = sample_crop_bounds(
            mask,
            self.patch_size,
            rng,
            foreground_fraction=self.foreground_fraction,
            )
        image, mask = crop_pair(image, mask, bounds)

        if self.augment:
            image, mask = augment_pair(image, mask, rng)

        target = make_target(
            mask,
            mode=self.target_mode,
            support_radius=self.support_radius,
            )

        image = torch.from_numpy(image[None, ...].astype(np.float32))
        mask = torch.from_numpy(target.astype(np.float32))

        return {
            'image': image,
            'mask': mask,
            'session': row['session'],
            'animal': row['animal'],
            'image_path': row['image_path'],
            }

    def load_pair(self, row):
        cache_key = row['session']
        if self.cache_images and cache_key in self.cache:
            return self.cache[cache_key]

        image = np.load(row['image_path'])
        roi_dict = load_roi_dict(row['roi_path'])
        mask, _ = roi_dict_to_mask(roi_dict, image.shape)

        image = robust_normalise(
            image,
            low=self.normalise_percentiles[0],
            high=self.normalise_percentiles[1],
            )
        mask = mask.astype(np.float32)

        if self.cache_images:
            self.cache[cache_key] = (image, mask)

        return image, mask


#%% direct loading
def load_full_pair(row, normalise_percentiles=(1, 99.7)):
    image = np.load(row['image_path'])
    roi_dict = load_roi_dict(row['roi_path'])
    mask, _ = roi_dict_to_mask(roi_dict, image.shape)

    image = robust_normalise(
        image,
        low=normalise_percentiles[0],
        high=normalise_percentiles[1],
        )
    return image.astype(np.float32), mask.astype(np.float32)
