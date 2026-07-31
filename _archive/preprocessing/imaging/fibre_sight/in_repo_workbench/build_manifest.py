'''
Created on 7 April 2026
Modified on 2 June 2026
Modified on 23 June 2026
scan labelled sessions and write the training table

@author: Dinghao Luo
'''

#%% imports
from pathlib import Path
import argparse
import sys

if __package__ in {None, ''}:
    repo_root = Path(__file__).resolve().parents[3]
    if str(repo_root) not in sys.path:
        sys.path.insert(0, str(repo_root))

from preprocessing.imaging.fibre_sight._repo import default_output_root, default_source_root
from preprocessing.imaging.fibre_sight.manifest import (
    assign_session_splits,
    scan_source_root,
    write_manifest,
    )

from utils.console_formatting import print_files_saved


#%% cli
def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--source-root', type=Path, default=default_source_root())
    parser.add_argument(
        '--out',
        type=Path,
        default=default_output_root() / 'manifests' / 'dlight_hpc_lc_opto_manifest.csv',
        )
    parser.add_argument('--val-fraction', type=float, default=0.15)
    parser.add_argument('--test-fraction', type=float, default=0.15)
    parser.add_argument('--seed', type=int, default=7)
    parser.add_argument('--no-splits', action='store_true')
    return parser.parse_args()


def print_summary(rows):
    included = [row for row in rows if row['included']]
    animals = sorted({row['animal'] for row in included})

    print(f'total sessions: {len(rows)}')
    print(f'included sessions: {len(included)}')
    print(f'animals: {animals}')
    print(f'total ROIs: {sum(row["roi_count"] for row in included)}')

    split_counts = {}
    for row in included:
        split_counts[row['split']] = split_counts.get(row['split'], 0) + 1
    if split_counts:
        print(f'splits: {split_counts}')


def main():
    args = parse_args()
    rows = scan_source_root(args.source_root)

    if not args.no_splits:
        rows = assign_session_splits(
            rows,
            val_fraction=args.val_fraction,
            test_fraction=args.test_fraction,
            seed=args.seed,
            )

    write_manifest(rows, args.out)
    print_summary(rows)
    print_files_saved([
        ('manifest', args.out),
    ])


if __name__ == '__main__':
    main()
