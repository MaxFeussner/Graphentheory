#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 12 11:53:27 2021

@author: TMC
"""

import os
import pandas as pd
from asymmetree.datastructures import PhyloTree
import asymmetree.treeevolve as te
import asymmetree.hgt as hgt
from asymmetree.cograph import Cotree
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
import networkx as netx


if '__file__' in vars():
    wk_dir = Path(os.path.dirname(os.path.realpath('__file__')))
else:
    print('We are running the script interactively')
    wk_dir = Path('/Users/TMC/Documents/01_Studium/02_Master/3_Semester/01_Graphentheorie/01_Prkatikum/Script/analyzeTrees.py').parent

# create output directory for data
Path(wk_dir / '01_Data').mkdir(parents=True, exist_ok=True)


# %% load parameter csv
parameter_Df = pd.read_csv(Path(wk_dir / '01_Data') / '01_Simulation_Parameters.csv')


# %% 2.1 Dependece between size of gene tree and numer of species and genes

# %% Iterate over all data

edges_ldt = []
leafes_s = []
leafes_tgt = []
fitch_edges = []
fraction_of_xenologs = []

ind = 0
for item in enumerate(parameter_Df.ID):
    path_s = str(item[1]) + '_species_tree.pickle'
    path_tgt = str(item[1]) + '_gene_tree.pickle'
    s = PhyloTree.load(Path(wk_dir / '01_Data' / path_s))
    tgt = PhyloTree.load(Path(wk_dir / '01_Data' / path_tgt))
    ogt = te.observable_tree(tgt)
    ldt = hgt.ldt_graph(ogt, s)
    transfer_edges = hgt.rs_transfer_edges(ogt, s)
    fitch = hgt.undirected_fitch(ogt, transfer_edges)
    fitch_edges.append(fitch.number_of_edges())
    edges_ldt.append(len(ldt.edges()))
    leafes_s.append(s.number_of_species)
    leafes_tgt.append(len(tgt.color_sorted_leaves()))
    ind += 1
    if ind == 100:
        print('Test')
    print(ind)

a = np.array(edges_ldt)
b = np.array(fitch_edges)
parameter_Df['LDT_Edges'] = edges_ldt
parameter_Df['Fitch_Edges'] = fitch_edges
parameter_Df['Fraction_of_Xenologs'] = np.divide(a, b, out = np.zeros_like(a), where=b != 0)
parameter_Df['Number_of_Species'] = leafes_s
parameter_Df['Number_of_leaves_tgt'] = leafes_tgt
parameter_Df.to_csv(Path(wk_dir / 'Tree_data.csv', index=False))
# %% Funktion to search in a pandas dataframe for an entry
def which(self)->int:
    '''
    Funktion to search in a pandas dataframe for an entry in an specific Column.
    Usage: which(nameOfDataFrame.ColumnName == "searchPattern")
    Raises
    ------
    Exception
        DESCRIPTION.
    Returns
    -------
    int
        Row Index of search pattern.
    '''
    try:
        self = list(iter(self))
    except TypeError as e:
        raise Exception("""'which' method can only be applied to iterables.
        {}""".format(str(e)))
    indices = [i for i, x in enumerate(self) if bool(x) == True]
    return(indices)

which(parameter_Df.ID == 'P0_0100')
