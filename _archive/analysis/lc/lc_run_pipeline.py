# -*- coding: utf-8 -*-
'''
Created on Thu Feb 27 11:53:55 2025
Originally named lc_run_all.py

run MATLAB conversion, LC extraction, and profile building in recording order

MATLAB spike preprocessing must already have run; the current entry points are
    under 'preprocessing/matlab_spike_pipeline/'

@author: Dinghao Luo
'''


#%% import scripts
import argparse
import io
import importlib
import sys
import time
from contextlib import redirect_stdout, redirect_stderr
from pathlib import Path

repo_root = Path(__file__).resolve().parents[2]

lc_preprocessing_root = repo_root / 'preprocessing' / 'ephys' / 'lc'
lc_analysis_root = repo_root / 'analysis' / 'lc'
lc_plotting_root = repo_root / 'plotting' / 'lc'
for search_root in (lc_preprocessing_root, lc_analysis_root, lc_plotting_root):
    if str(search_root) not in sys.path:
        sys.path.insert(0, str(search_root))


class StageLog:
    def __init__(self, stream, prefix='  '):
        self.stream = stream
        self.prefix = prefix
        self.at_line_start = True

    def write(self, text):
        if text == '':
            return
        text = text.replace('\r', '\n')
        for chunk in text.splitlines(keepends=True):
            if chunk.strip() == '':
                self.at_line_start = chunk.endswith('\n')
                continue
            if self.at_line_start and chunk.strip():
                self.stream.write(self.prefix)
            self.stream.write(chunk)
            self.at_line_start = chunk.endswith('\n')
        self.stream.flush()

    def flush(self):
        self.stream.flush()


class QuietLog:
    def __init__(self, max_chars=20000):
        self.text = io.StringIO()
        self.max_chars = max_chars

    def write(self, text):
        self.text.write(text)
        if self.text.tell() > self.max_chars:
            current_text = self.text.getvalue()[-self.max_chars:]
            self.text = io.StringIO()
            self.text.write(current_text)

    def flush(self):
        pass

    def getvalue(self):
        return self.text.getvalue()


STAGES = [
    {
        'label': 'Spike times and ISIs',
        'module': 'lc_all_spikes_isis',
        'description': 'Read LC spike times and save per-cell ISI summaries.',
        },
    {
        'label': 'Waveforms, ACGs, and tagged identities',
        'module': 'lc_all_waveforms_acgs',
        'description': 'Save per-session waveforms, ACGs, and tagged-cell identities.',
        },
    {
        'label': 'Run-aligned spike extraction',
        'module': 'lc_all_extract',
        'argv': [],
        'description': 'Build per-session spike rasters and smoothed spike trains.',
        },
    {
        'label': 'ACG UMAP and k-means',
        'module': 'lc_all_identity_umap',
        'description': (
            'Use tagged Dbh+ units as references to classify the ACG '
            'UMAP/k-means split as putative Dbh+ or Dbh-.'
            ),
        },
    {
        'label': 'Plot identity UMAP',
        'module': 'plot_lc_all_identity_umap',
        'description': 'Render the UMAP identity figure from saved analysis output.',
        },
    {
        'label': 'Cell profiles',
        'module': 'lc_all_profiles',
        'argv': [],
        'profile_workers': True,
        'verbose_arg': True,
        'description': 'Compile cell properties into LC_all_cell_profiles.pkl.',
        },
    {
        'label': 'Plot cell profiles',
        'module': 'plot_lc_all_profiles',
        'description': 'Render profile summary figures from the saved profile table.',
        },
    {
        'label': 'Done',
        'module': None,
        'description': 'LC pipeline complete.',
        },
    ]


def build_arg_parser():
    parser = argparse.ArgumentParser(
        description='run LC processing pipeline'
        )
    parser.add_argument(
        '-v',
        '--verbose',
        '-verbose',
        action='store_true',
        help='print output from each stage',
        )
    parser.add_argument(
        '--n-workers',
        '--workers',
        dest='n_workers',
        type=int,
        default=None,
        help='number of workers for stages that support it',
        )
    return parser


def elapsed_text(seconds):
    if seconds < 60:
        return f'{seconds:.1f}s'
    minutes, seconds = divmod(seconds, 60)
    if minutes < 60:
        return f'{int(minutes)}m {seconds:.0f}s'
    hours, minutes = divmod(minutes, 60)
    return f'{int(hours)}h {int(minutes)}m {seconds:.0f}s'


def progress_bar(number, total, width=28):
    done = round(width * number / total)
    return '#' * done + '-' * (width - done)


def progress_line(stage, number, total):
    bar = progress_bar(number, total)
    return f'[{bar}] {number}/{total} {stage["label"]}'


def get_stage_args(stage, args):
    if 'argv' not in stage:
        return stage.get('args', ())

    argv = list(stage['argv'])
    if stage.get('verbose_arg') and args.verbose:
        argv.append('--verbose')
    if stage.get('profile_workers') and args.n_workers is not None:
        argv.extend(['--n-workers', str(args.n_workers)])
    return (argv,)


def run_stage(stage, number, total, args):
    if stage['module'] is None:
        return

    if args.verbose:
        print(f'\n[{number}/{total}] {stage["label"]}')
        print(f'  {stage["description"]}', flush=True)
        out_log = StageLog(sys.stdout)
        err_log = StageLog(sys.stderr)
    else:
        print(progress_line(stage, number, total), flush=True)
        out_log = QuietLog()
        err_log = QuietLog()

    start = time.perf_counter()
    try:
        with redirect_stdout(out_log), redirect_stderr(err_log):
            module = importlib.import_module(stage['module'])
            module.main(*get_stage_args(stage, args))
    except Exception:
        elapsed = elapsed_text(time.perf_counter() - start)
        if not args.verbose:
            saved_output = out_log.getvalue() + err_log.getvalue()
            if saved_output.strip():
                print(saved_output.rstrip())
        print(f'[{number}/{total}] failed after {elapsed}')
        raise

    elapsed = elapsed_text(time.perf_counter() - start)
    if args.verbose:
        print(f'[{number}/{total}] done ({elapsed})')


#%% pipeline proper
def main(argv=None):
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    start = time.perf_counter()
    runnable_stages = [stage for stage in STAGES if stage['module'] is not None]

    if args.verbose:
        print(f'LC pipeline ({len(runnable_stages)} stages)')

    for number, stage in enumerate(runnable_stages, start=1):
        run_stage(stage, number, len(runnable_stages), args)

    total_elapsed = elapsed_text(time.perf_counter() - start)
    if args.verbose:
        print(f'\nLC pipeline complete ({total_elapsed})')
    else:
        print(f'LC pipeline complete ({total_elapsed})')


if __name__ == '__main__':
    main(sys.argv[1:])
