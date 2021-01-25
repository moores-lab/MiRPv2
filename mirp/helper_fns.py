#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
MiRP - a microtubule RELION-based pipeline for cryo-EM image processing. 
MiRP v2 is designed to function with RELION v3.1, and is not compatible with earlier versions of RELION.
"""

__author__ = 'Alexander D. Cook & Joseph Atherton'
__license__ = 'GPLv3'
__version__ = '2.0'


from itertools import groupby
from operator import itemgetter
from collections import OrderedDict
import ast
import numpy as np


def get_window(index, lw, hi, size):
    lwin, hwin = index - lw, index + hi
    if lwin < 0 :
        lwin = 0
    if hwin > size:
        hwin = size
    return lwin, hwin


def group_dict_of_list(dict, *labels):
    keys = dict.keys()
    data = zip(*dict.values())
    grouped = []
    labels = [index_from_odict(dict, lbl) for lbl in labels]
    for _, grp in groupby(data, itemgetter(*labels)):
        group = OrderedDict()
        for k, v in zip(keys, zip(*grp)):
            group[k] = list(v)
        grouped.append(group)
    return grouped


def sort_dict_of_list(dict, *keys):
    data = trnsp_dict_of_lst(dict)
    data = sorted(data, key=itemgetter(*keys))
    return trnsp_lst_of_dict(data)


def split_dict_of_list(dict, indices):
    split = []
    split_vals = zip(*[np.split(val, indices) for val in dict.values()])
    for frac in split_vals:
        new = OrderedDict()
        for key, val in zip(dict.keys(), frac):
            new[key] = val
        split.append(new)
    return split


def index_from_odict(dict, key):
    for i, k in enumerate(dict.keys()):
        if k == key:
            return i
    

def trnsp_lst_of_dict(lod):
    dol = OrderedDict()
    for k in lod[0].keys():
        dol[k] = []
    for odict in lod:
        for key in odict:
            dol[key].append(odict[key])
    return dol


def trnsp_dict_of_lst(dict_of_list):
    return [OrderedDict(zip(dict_of_list, col))
            for col in zip(*dict_of_list.values())]


def literal_eval(var):
    try: var = ast.literal_eval(var)
    except: var = str(var)
    return var


def readfile(filepath):
    with open(filepath, 'r') as f:
        for line in f:
            fields = line.split()
            if not fields:
                pass
            elif fields[0].startswith('#'):
                pass
            else:
                yield fields