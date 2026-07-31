'''
Created on 27 April 2026
Modified on 2 June 2026
Modified on 23 June 2026
compare held-out model predictions with curated ROIs

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

from preprocessing.imaging.fibre_sight._repo import add_repo_paths, default_figure_root
from preprocessing.imaging.fibre_sight.api import AxonROIPredictor
from preprocessing.imaging.fibre_sight.image_ops import robust_normalise
from preprocessing.imaging.fibre_sight.manifest import read_manifest
from preprocessing.imaging.fibre_sight.plot_training_data import choose_rows, make_overlay
from preprocessing.imaging.fibre_sight.roi_io import load_roi_dict, roi_dict_to_label

add_repo_paths()
from utils.console_formatting import print_files_saved


#%% cli
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--split', default='test')
    parser.add_argument('--model-name', default='dlight_hpc_lc_opto')
    parser.add_argument('--checkpoint', type=Path, default=None)
    parser.add_argument('--threshold', type=float, default=None)
    parser.add_argument('--min-size', type=int, default=None)
    parser.add_argument('--no-tta', action='store_true')
    parser.add_argument('--n', type=int, default=4)
    parser.add_argument('--seed', type=int, default=17)
    parser.add_argument(
        '--out',
        type=Path,
        default=default_figure_root() / 'diagnostics' / 'model_prediction_overlays.png',
        )
    return parser.parse_args()


#%% plotting
def make_error_overlay(image, pred_mask, target_mask, alpha=0.6):
    base = robust_normalise(image)
    out = np.dstack([base, base, base])

    # keep the error colours fixed across sessions so each panel reads the same way
    true_positive = pred_mask & target_mask
    false_positive = pred_mask & ~target_mask
    false_negative = ~pred_mask & target_mask

    out[true_positive] = blend(out[true_positive], np.array([0.0, 0.9, 0.2]), alpha)
    out[false_positive] = blend(out[false_positive], np.array([1.0, 0.1, 0.8]), alpha)
    out[false_negative] = blend(out[false_negative], np.array([0.0, 0.45, 1.0]), alpha)
    return out


def blend(values, colour, alpha):
    return (1 - alpha) * values + alpha * colour


def plot_rows(rows, predictor, out_path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    from common_functions import mpl_formatting
    mpl_formatting()

    fig, axes = plt.subplots(
        len(rows),
        3,
        figsize=(10.5, 3.8 * len(rows)),
        constrained_layout=True,
        )
    axes = np.atleast_2d(axes)

    for row_idx, row in enumerate(rows):
        image = np.load(row['image_path'])
        target_dict = load_roi_dict(row['roi_path'])
        target_labelled, _, _ = roi_dict_to_label(target_dict, image.shape)
        prediction = predictor.predict_image(image)

        pred_mask = prediction.labelled > 0
        target_mask = target_labelled > 0

        axes[row_idx, 0].imshow(make_overlay(image, target_labelled), interpolation='nearest')
        axes[row_idx, 1].imshow(make_overlay(image, prediction.labelled), interpolation='nearest')
        axes[row_idx, 2].imshow(make_error_overlay(image, pred_mask, target_mask), interpolation='nearest')

        axes[row_idx, 0].set_ylabel(row['session'], fontsize=8)
        axes[row_idx, 0].set_title(f'curated ({int(row["roi_count"])} ROIs)', fontsize=9)
        axes[row_idx, 1].set_title(f'predicted ({len(prediction.roi_dict)} ROIs)', fontsize=9)
        axes[row_idx, 2].set_title('green TP, magenta FP, blue FN', fontsize=9)

    for ax in axes.ravel():
        ax.set_xticks([])
        ax.set_yticks([])

    fig.suptitle('held-out axon ROI model diagnostics', fontsize=12)
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=220, bbox_inches='tight')
    plt.close(fig)


def main():
    args = parse_args()
    rows = read_manifest(args.manifest, included_only=True, split=args.split)
    rows = choose_rows(rows, args.n, seed=args.seed)

    predictor = AxonROIPredictor(
        checkpoint_path=args.checkpoint,
        model_name=args.model_name,
        threshold=args.threshold,
        min_size=args.min_size,
        tta=False if args.no_tta else None,
    )
    plot_rows(rows, predictor, args.out)
    print_files_saved([
        ('plot', args.out),
    ], gap=1)


if __name__ == '__main__':
    main()
