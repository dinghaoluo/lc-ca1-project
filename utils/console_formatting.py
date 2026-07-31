# -*- coding: utf-8 -*-
'''
Created on 5 June 2026

console messages used during long recording and pipeline runs

@author: Dinghao Luo
'''

def print_session(name):
    print(f'\n{name}')

def print_status(label, message=None):
    label = str(label).strip().lower()
    if message is None or message == '':
        print(label)
    else:
        print(f'{label}: {message}')

def print_statistics_section(gap=2):
    bar = '-' * 11
    prefix = '\n' * gap
    print(f'{prefix}{bar}\nSTATISTICS\n{bar}')

def print_files_saved(entries, gap=1):
    bar = '-' * 12
    prefix = '\n' * gap
    print(f'{prefix}{bar}\nFILES SAVED\n{bar}')
    for label, path in entries:
        print(f'{label}: {path}')

def print_binwise_header(label, group0, group1, tests):
    print(f'\nbinwise stats for {label}:')
    print(
        f'bin | {group0} mean +/- SEM | {group1} mean +/- SEM | '
        + ' | '.join(tests)
    )

def print_binwise_row(bin_label, mean0, sem0, mean1, sem1, pvals):
    print(
        f'{bin_label:>15s} | '
        f'{mean0:.3f} +/- {sem0:.3f} | '
        f'{mean1:.3f} +/- {sem1:.3f} | '
        + ' '.join(f'{pval:.2e}' for pval in pvals)
    )
