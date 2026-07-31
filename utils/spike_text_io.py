# -*- coding: utf-8 -*-
'''
Created on Tue Aug  2 18:51:23 2022

spike text-file readers

@author: Dinghao Luo
'''

import numpy as np

def param2array(filename):
    '''read a line-based spike parameter file into a 1d array'''
    with open(filename, 'r', encoding='utf-8-sig') as file:
        string = file.read()

    return np.asarray(string.split('\n'))

def get_clu(n, clu):
    '''return spike indices for one cluster id'''
    if isinstance(n, int):
        n = str(n)

    return np.array(np.where(clu == n))
