'''
Created on 29 April 2026
Modified on 3 June 2026
model loading and prediction used by the workbench

@author: Dinghao Luo
'''

#%% imports
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from ._repo import get_repo_root
from .config import load_config, resolve_path
from .fibre_segger_params import recommend_fibre_segger_params
from .postprocess import probability_to_roi_dict
from .predict_rois import load_trained_model, predict_probability
from .roi_io import save_roi_dict


#%% data structures
@dataclass
class AxonROIPrediction:
    roi_dict: dict
    labelled: np.ndarray
    probability: np.ndarray
    fibre_segger_params: dict
    checkpoint: Path
    threshold: float
    min_size: int
    tta: bool


#%% registry
def load_model_registry(path=None):
    if path is None:
        path = Path(__file__).resolve().parent / 'configs' / 'model_registry.yaml'
    return load_config(path)


def get_model_entry(model_name='dlight_hpc_lc_opto', registry_path=None, required=True):
    registry = load_model_registry(registry_path)
    models = registry.get('models', {})
    if model_name not in models:
        if required:
            raise KeyError(f'model not found in registry: {model_name}')
        return {}
    return models[model_name] or {}


def get_default_checkpoint(model_name='dlight_hpc_lc_opto', registry_path=None):
    entry = get_model_entry(model_name, registry_path=registry_path)
    return resolve_path(entry['checkpoint'], get_repo_root())


#%% predictor
class AxonROIPredictor:
    def __init__(self, checkpoint_path=None, model_name='dlight_hpc_lc_opto',
                 registry_path=None, device='auto', threshold=None, min_size=None,
                 tta=None):
        self.model_name = model_name
        self.registry_path = registry_path
        self.entry = get_model_entry(
            model_name,
            registry_path=registry_path,
            required=checkpoint_path is None,
            )
        self.checkpoint_path = Path(checkpoint_path) if checkpoint_path else get_default_checkpoint(
            model_name,
            registry_path=registry_path,
            )
        self.device = self._get_device(device)
        self.model = None
        self.checkpoint = None
        self.threshold = threshold
        self.min_size = min_size
        self.max_size = self.entry.get('max_size', None)
        self.tta = tta

    def load(self):
        if not self.checkpoint_path.exists():
            raise FileNotFoundError(f'model checkpoint not found: {self.checkpoint_path}')
        self.model, self.checkpoint = load_trained_model(self.checkpoint_path, self.device)
        post_cfg = dict(self.checkpoint.get('postprocess_config', {}))
        # the registry holds the operating point I settled on after tuning;
        # older checkpoints keep the defaults used during training
        if self.threshold is None:
            self.threshold = float(self.entry.get('threshold', post_cfg.get('threshold', 0.5)))
        else:
            self.threshold = float(self.threshold)
        if self.min_size is None:
            self.min_size = int(self.entry.get('min_size', post_cfg.get('min_size', 30)))
        else:
            self.min_size = int(self.min_size)
        if self.max_size is None:
            self.max_size = post_cfg.get('max_size', None)
        if self.tta is None:
            self.tta = bool(self.entry.get('tta', False))
        else:
            self.tta = bool(self.tta)
        return self

    def predict_image(self, image):
        if self.model is None:
            self.load()

        data_cfg = self.checkpoint.get('data_config', {})
        probability = predict_probability(
            image,
            self.model,
            self.device,
            normalise_percentiles=data_cfg.get('normalise_percentiles', [1, 99.7]),
            tta=self.tta,
            )
        # keep model output in the xpix/ypix dictionary used by the imaging pipeline
        roi_dict, labelled = probability_to_roi_dict(
            probability,
            threshold=self.threshold,
            min_size=self.min_size,
            max_size=self.max_size,
            )
        params = recommend_fibre_segger_params(image, roi_dict=roi_dict)

        return AxonROIPrediction(
            roi_dict=roi_dict,
            labelled=labelled,
            probability=probability,
            fibre_segger_params=params,
            checkpoint=self.checkpoint_path,
            threshold=self.threshold,
            min_size=self.min_size,
            tta=self.tta,
            )

    def predict_file(self, image_path, out_path=None):
        image_path = Path(image_path)
        image = np.load(image_path)
        prediction = self.predict_image(image)

        if out_path is not None:
            save_roi_dict(prediction.roi_dict, out_path)

        return prediction

    @staticmethod
    def _get_device(device):
        import torch

        if device == 'auto':
            device = 'cuda' if torch.cuda.is_available() else 'cpu'
        return torch.device(device)


#%% convenience functions
def predict_rois_for_gui(image, checkpoint_path=None, model_name='dlight_hpc_lc_opto',
                         device='auto', threshold=None, min_size=None, tta=None):
    predictor = AxonROIPredictor(
        checkpoint_path=checkpoint_path,
        model_name=model_name,
        device=device,
        threshold=threshold,
        min_size=min_size,
        tta=tta,
        )
    return predictor.predict_image(image)


def predict_roi_file_for_gui(image_path, out_path=None, checkpoint_path=None,
                             model_name='dlight_hpc_lc_opto', device='auto',
                             threshold=None, min_size=None, tta=None):
    predictor = AxonROIPredictor(
        checkpoint_path=checkpoint_path,
        model_name=model_name,
        device=device,
        threshold=threshold,
        min_size=min_size,
        tta=tta,
        )
    return predictor.predict_file(image_path, out_path=out_path)
