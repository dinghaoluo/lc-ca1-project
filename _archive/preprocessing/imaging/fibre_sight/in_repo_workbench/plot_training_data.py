'''
Created on 11 April 2026
Modified on 2 June 2026
Modified on 23 June 2026
inspect channel-2 training images with their curated labels

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
from preprocessing.imaging.fibre_sight.image_ops import robust_normalise
from preprocessing.imaging.fibre_sight.manifest import read_manifest
from preprocessing.imaging.fibre_sight.roi_io import load_roi_dict, roi_dict_to_label

add_repo_paths()
from utils.console_formatting import print_files_saved


#%% cli
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--split', default=None)
    parser.add_argument('--n', type=int, default=6)
    parser.add_argument('--seed', type=int, default=7)
    parser.add_argument(
        '--out',
        type=Path,
        default=default_figure_root() / 'diagnostics' / 'training_label_overlays.png',
        )
    return parser.parse_args()


#%% plotting
def make_overlay(image, labelled, alpha=0.55):
    base = robust_normalise(image)
    out = np.dstack([base, base, base])

    if labelled.max() == 0:
        return out

    # stable colours make repeat previews comparable while the sampled sessions change
    rng = np.random.default_rng(20260512)
    colours = rng.uniform(0.1, 1.0, size=(int(labelled.max()) + 1, 3))
    overlay = colours[labelled]
    mask = labelled > 0
    out[mask] = (1 - alpha) * out[mask] + alpha * overlay[mask]

    boundary = mask & (
        (np.roll(labelled, 1, axis=0) != labelled) |
        (np.roll(labelled, -1, axis=0) != labelled) |
        (np.roll(labelled, 1, axis=1) != labelled) |
        (np.roll(labelled, -1, axis=1) != labelled)
        )
    out[boundary] = [1, 0.95, 0.05]
    return out


def choose_rows(rows, n_rows, seed=7):
    if len(rows) <= n_rows:
        return rows

    rng = np.random.default_rng(seed)
    idx = rng.choice(len(rows), size=n_rows, replace=False)
    return [rows[int(i)] for i in idx]


def plot_rows(rows, out_path):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    from common_functions import mpl_formatting
    mpl_formatting()

    n_rows = len(rows)
    n_cols = min(3, n_rows)
    n_fig_rows = int(np.ceil(n_rows / n_cols))

    fig, axes = plt.subplots(
        n_fig_rows,
        n_cols,
        figsize=(4.2 * n_cols, 4.5 * n_fig_rows),
        constrained_layout=True,
        )
    axes = np.atleast_1d(axes).ravel()

    for ax, row in zip(axes, rows):
        image = np.load(row['image_path'])
        roi_dict = load_roi_dict(row['roi_path'])
        labelled, _, _ = roi_dict_to_label(roi_dict, image.shape)

        ax.imshow(make_overlay(image, labelled), interpolation='nearest')
        ax.set_title(
            f'{row["session"]}\n'
            f'{row["roi_count"]} ROIs, {100 * row["positive_fraction"]:.2f}% labelled',
            fontsize=9,
            )
        ax.set_axis_off()

    for ax in axes[len(rows):]:
        ax.set_axis_off()

    fig.suptitle('channel-2 references with curated axon ROIs', fontsize=12)
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=220, bbox_inches='tight')
    plt.close(fig)


def main():
    args = parse_args()
    rows = read_manifest(args.manifest, included_only=True, split=args.split)
    rows = choose_rows(rows, args.n, seed=args.seed)
    plot_rows(rows, args.out)
    print_files_saved([
        ('plot', args.out),
    ], gap=1)


if __name__ == '__main__':
    main()
