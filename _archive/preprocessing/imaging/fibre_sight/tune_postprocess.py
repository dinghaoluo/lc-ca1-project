'''
Created on 25 April 2026
Modified on 2 June 2026
Modified on 23 June 2026
Archived on 23 July 2026 after settling on threshold 0.25 and minimum size 45
reuse saved probability maps while tuning threshold and ROI size

@author: Dinghao Luo
'''

#%% imports
from pathlib import Path
import argparse
import csv
import sys

import numpy as np

if __package__ in {None, ''}:
    repo_root = Path(__file__).resolve().parents[4]
    if str(repo_root) not in sys.path:
        sys.path.insert(0, str(repo_root))

from preprocessing.imaging.fibre_sight._repo import add_repo_paths
from preprocessing.imaging.fibre_sight.manifest import read_manifest
from preprocessing.imaging.fibre_sight.postprocess import probability_to_labels
from preprocessing.imaging.fibre_sight.predict_rois import load_trained_model, predict_probability
from preprocessing.imaging.fibre_sight.roi_io import load_roi_dict, roi_dict_to_label

add_repo_paths()
from utils.console_formatting import print_files_saved


#%% scoring
# This copy kept the threshold sweep usable after it left the active workbench.
def mask_scores(pred_mask, target_mask, eps=1e-8):
    pred_mask = np.asarray(pred_mask).astype(bool)
    target_mask = np.asarray(target_mask).astype(bool)

    tp = np.sum(pred_mask & target_mask)
    fp = np.sum(pred_mask & ~target_mask)
    fn = np.sum(~pred_mask & target_mask)

    precision = tp / (tp + fp + eps)
    recall = tp / (tp + fn + eps)
    dice = 2 * tp / (2 * tp + fp + fn + eps)
    iou = tp / (tp + fp + fn + eps)
    f2 = 5 * tp / (5 * tp + 4 * fn + fp + eps)

    return {
        'dice': float(dice),
        'f2': float(f2),
        'iou': float(iou),
        'precision': float(precision),
        'recall': float(recall),
        'true_positive_pixels': int(tp),
        'false_positive_pixels': int(fp),
        'false_negative_pixels': int(fn),
        }


def component_counts(pred_labelled, target_labelled):
    pred_ids = set(np.unique(pred_labelled)) - {0}
    target_ids = set(np.unique(target_labelled)) - {0}

    return {
        'predicted_components': len(pred_ids),
        'target_components': len(target_ids),
        }


#%% cli
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--checkpoint', type=Path, required=True)
    parser.add_argument('--split', default='val')
    parser.add_argument('--thresholds', default='0.25,0.30,0.35,0.40,0.45,0.50,0.55,0.60,0.65,0.70,0.75')
    parser.add_argument('--min-sizes', default='15,30,60,100,150')
    parser.add_argument('--metric', default='dice')
    parser.add_argument('--tta', action='store_true')
    parser.add_argument('--device', default='auto')
    parser.add_argument('--out', type=Path, default=None)
    return parser.parse_args()


def parse_values(text, dtype=float):
    return [dtype(value.strip()) for value in text.split(',') if value.strip()]


#%% prediction cache
def prepare_items(rows, checkpoint_path, device='auto', tta=False):
    import torch

    if device == 'auto':
        device = 'cuda' if torch.cuda.is_available() else 'cpu'
    device = torch.device(device)

    model, checkpoint = load_trained_model(checkpoint_path, device)
    data_cfg = checkpoint.get('data_config', {})

    # inference is the slow part, so every threshold and size setting reuses
    # the same probability maps
    items = []
    for row in rows:
        image = np.load(row['image_path'])
        roi_dict = load_roi_dict(row['roi_path'])
        target_labelled, _, _ = roi_dict_to_label(roi_dict, image.shape)
        probability = predict_probability(
            image,
            model,
            device,
            normalise_percentiles=data_cfg.get('normalise_percentiles', [1, 99.7]),
            tta=tta,
            )
        items.append({
            'row': row,
            'target_labelled': target_labelled,
            'probability': probability,
            })

    return items


#%% tuning
def score_setting(items, threshold, min_size):
    rows = []

    for item in items:
        pred_labelled = probability_to_labels(
            item['probability'],
            threshold=threshold,
            min_size=min_size,
            )
        target_labelled = item['target_labelled']

        scores = mask_scores(pred_labelled > 0, target_labelled > 0)
        scores.update(component_counts(pred_labelled, target_labelled))
        rows.append(scores)

    return summarise_scores(rows, threshold, min_size)


def summarise_scores(rows, threshold, min_size):
    summary = {
        'threshold': threshold,
        'min_size': min_size,
        'n_sessions': len(rows),
        }

    for key in ['dice', 'f2', 'iou', 'precision', 'recall',
                'predicted_components', 'target_components']:
        values = np.array([row[key] for row in rows], dtype=float)
        summary[f'mean_{key}'] = float(np.mean(values))
        summary[f'std_{key}'] = float(np.std(values))

    return summary


def tune(items, thresholds, min_sizes, metric='dice'):
    results = []

    for threshold in thresholds:
        for min_size in min_sizes:
            results.append(score_setting(items, threshold, min_size))

    metric_key = f'mean_{metric}'
    best = max(results, key=lambda row: row[metric_key])
    return results, best


def write_results(rows, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    keys = list(rows[0].keys())

    with open(path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def main():
    args = parse_args()
    rows = read_manifest(args.manifest, included_only=True, split=args.split)
    items = prepare_items(rows, args.checkpoint, device=args.device, tta=args.tta)

    results, best = tune(
        items,
        parse_values(args.thresholds, float),
        parse_values(args.min_sizes, int),
        metric=args.metric,
        )

    out_path = args.out
    if out_path is None:
        out_path = Path(args.checkpoint).parent / f'tune_postprocess_{args.split}.csv'
    write_results(results, out_path)

    print(
        f'best {args.metric}: '
        f'threshold={best["threshold"]}, '
        f'min_size={best["min_size"]}, '
        f'mean_{args.metric}={best[f"mean_{args.metric}"]:.4f}'
        )
    print_files_saved([
        ('tuning', out_path),
    ])


if __name__ == '__main__':
    main()
