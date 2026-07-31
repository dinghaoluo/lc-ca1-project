# -*- coding: utf-8 -*-
'''
Created on Tue Jun  2 2026

stage selection and subprocess reporting for the pipeline wrappers

@author: Dinghao Luo
'''


#%% imports
from pathlib import Path
import subprocess
import sys
import time


#%% stage selection
def select_stage_defs(stage_defs, args, parser):
    if args.branch == 'all':
        selected = list(stage_defs)
    else:
        selected = [
            stage
            for stage in stage_defs
            if stage['branch'] == args.branch
            ]

    selected_keys = [stage['key'] for stage in selected]
    if args.from_stage is not None:
        if args.from_stage not in selected_keys:
            parser.error(f'{args.from_stage} is not in branch {args.branch}')
        selected = selected[selected_keys.index(args.from_stage):]

    if args.only:
        selected_keys = [stage['key'] for stage in selected]
        missing = [
            key
            for key in args.only
            if key not in selected_keys
            ]
        if missing:
            missing_text = ', '.join(missing)
            parser.error(f'{missing_text} not selected by branch/from-stage')
        wanted = set(args.only)
        selected = [
            stage
            for stage in selected
            if stage['key'] in wanted
            ]

    return selected


#%% formatting
def elapsed_text(seconds):
    if seconds < 60:
        return f'{seconds:.1f}s'
    minutes, seconds = divmod(seconds, 60)
    if minutes < 60:
        return f'{int(minutes)}m {seconds:.0f}s'
    hours, minutes = divmod(minutes, 60)
    return f'{int(hours)}h {int(minutes)}m {seconds:.0f}s'

def display_command(script_path, args=None):
    if args is None:
        args = []
    return subprocess.list2cmdline(
        [str(part) for part in ['python', '-u', str(script_path), *args]]
        )

def python_command(repo_root, script_path, args=None):
    if args is None:
        args = []
    script_path = Path(script_path)
    return [
        sys.executable,
        '-u',
        str(repo_root / script_path),
        *[str(arg) for arg in args],
        ]


#%% running
def run_pipeline(stages, repo_root, label, verbose=False, dry_run=False):
    pipeline_start = time.perf_counter()

    if verbose or dry_run:
        suffix = ' dry run' if dry_run else ''
        stage_text = 'stage' if len(stages) == 1 else 'stages'
        print(f'{label}{suffix} ({len(stages)} {stage_text})')

    total = len(stages)
    for number, stage in enumerate(stages, start=1):
        stage_key = stage['key']
        stage_label = stage['label']
        command = stage['command']
        display = stage.get('display') or subprocess.list2cmdline(
            [str(part) for part in command]
            )

        if dry_run:
            print(f'[{number}/{total}] {stage_key}')
            print(f'  {display}')
            continue

        if verbose:
            print(f'\n[{number}/{total}] {stage_label}')
            print(f'  {display}', flush=True)
        else:
            done = round(28 * number / total)
            bar = '#' * done + '-' * (28 - done)
            print(f'[{bar}] {number}/{total} {stage_label}', flush=True)

        stage_start = time.perf_counter()
        try:
            subprocess.run(command, cwd=repo_root, check=True)
        except subprocess.CalledProcessError:
            elapsed = elapsed_text(time.perf_counter() - stage_start)
            print(f'[{number}/{total}] failed after {elapsed}')
            raise

        elapsed = elapsed_text(time.perf_counter() - stage_start)
        if verbose:
            print(f'[{number}/{total}] done ({elapsed})')
        else:
            done = round(28 * number / total)
            bar = '#' * done + '-' * (28 - done)
            print(
                f'[{bar}] {number}/{total} done {stage_label} ({elapsed})',
                flush=True
                )

    total_elapsed = elapsed_text(time.perf_counter() - pipeline_start)
    if not dry_run:
        print(f'{label} complete ({total_elapsed})')
