'''
Created on 22 April 2026
Modified on 23 June 2026
run a trained model on channel-2 references and save ROI dictionaries

@author: Dinghao Luo
'''

#%% imports
from pathlib import Path
import argparse
import sys

import numpy as np

if __package__ in {None, ''}:
    repo_root = Path(__file__).resolve().parents[3]
    if str(repo_root) not in sys.path:
        sys.path.insert(0, str(repo_root))

from preprocessing.imaging.fibre_sight._repo import add_repo_paths
from preprocessing.imaging.fibre_sight.image_ops import robust_normalise
from preprocessing.imaging.fibre_sight.model import build_model
from preprocessing.imaging.fibre_sight.postprocess import probability_to_roi_dict
from preprocessing.imaging.fibre_sight.roi_io import save_roi_dict

add_repo_paths()


#%% cli
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--image', type=Path, required=True)
    parser.add_argument('--checkpoint', type=Path, required=True)
    parser.add_argument('--out', type=Path, default=None)
    parser.add_argument('--threshold', type=float, default=None)
    parser.add_argument('--min-size', type=int, default=None)
    parser.add_argument('--tta', action='store_true')
    parser.add_argument('--device', default='auto')
    return parser.parse_args()


#%% prediction
def load_trained_model(checkpoint_path, device):
    import torch

    checkpoint = torch.load(checkpoint_path, map_location=device)
    model = build_model(checkpoint.get('model_config', {}))
    model.load_state_dict(checkpoint['model_state'])
    model.to(device)
    model.eval()
    return model, checkpoint


def predict_probability(image, model, device, normalise_percentiles=(1, 99.7), tta=False):
    import torch

    image = robust_normalise(
        image,
        low=normalise_percentiles[0],
        high=normalise_percentiles[1],
        )
    tensor = torch.from_numpy(image[None, None, ...].astype(np.float32)).to(device)

    with torch.no_grad():
        probability = predict_tensor_probability(tensor, model, tta=tta)

    return probability[0, 0].detach().cpu().numpy()


def predict_tensor_probability(tensor, model, tta=False):
    import torch

    if not tta:
        return torch.sigmoid(model(tensor)[:, :1])

    # average the four flip views without interpolating the thin ROI shapes
    predictions = []
    for dims in [(), (2,), (3,), (2, 3)]:
        if dims:
            input_tensor = torch.flip(tensor, dims=dims)
        else:
            input_tensor = tensor

        probability = torch.sigmoid(model(input_tensor)[:, :1])
        if dims:
            probability = torch.flip(probability, dims=dims)
        predictions.append(probability)

    return torch.mean(torch.stack(predictions), dim=0)


def predict_roi_dict(image_path, checkpoint_path, out_path=None, threshold=None,
                     min_size=None, device='auto', tta=False):
    import torch

    if device == 'auto':
        device = 'cuda' if torch.cuda.is_available() else 'cpu'
    device = torch.device(device)

    model, checkpoint = load_trained_model(checkpoint_path, device)
    data_cfg = checkpoint.get('data_config', {})
    post_cfg = dict(checkpoint.get('postprocess_config', {}))

    if threshold is not None:
        post_cfg['threshold'] = threshold
    if min_size is not None:
        post_cfg['min_size'] = min_size

    image = np.load(image_path)
    probability = predict_probability(
        image,
        model,
        device,
        normalise_percentiles=data_cfg.get('normalise_percentiles', [1, 99.7]),
        tta=tta,
        )
    roi_dict, labelled = probability_to_roi_dict(
        probability,
        threshold=post_cfg.get('threshold', 0.5),
        min_size=post_cfg.get('min_size', 30),
        max_size=post_cfg.get('max_size', None),
        )

    if out_path is None:
        recname = Path(image_path).name.split('_ref_mat')[0]
        out_path = Path(image_path).parent / f'{recname}_ROI_dict.npy'

    save_roi_dict(roi_dict, out_path)
    return roi_dict, labelled, probability


def main():
    args = parse_args()
    roi_dict, _, _ = predict_roi_dict(
        args.image,
        args.checkpoint,
        out_path=args.out,
        threshold=args.threshold,
        min_size=args.min_size,
        device=args.device,
        tta=args.tta,
        )
    print(f'predicted ROIs: {len(roi_dict)}')


if __name__ == '__main__':
    main()
