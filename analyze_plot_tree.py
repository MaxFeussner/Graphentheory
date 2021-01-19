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
import graphFunctions as gf



if '__file__' in vars():
    wk_dir = Path(os.path.dirname(os.path.realpath('__file__')))
else:
    print('We are running the script interactively')
    wk_dir = Path('/Users/TMC/Documents/01_Studium/02_Master/3_Semester/01_Graphentheorie/01_Prkatikum/Cloned/Graphentheory/analyze_plot_tree.py').parent

# create output directory for data
Path(wk_dir / '01_Data').mkdir(parents=True, exist_ok=True)


# %% load parameter csv
parameter_Df = pd.read_csv(Path(wk_dir / '01_Data') / '01_Simulation_Parameters.csv')


# %% Function for 2.2




# %% 2.1 Dependece between size of gene tree and numer of species and genes
# %% Iterate over all data

edges_ldt = []
leafes_s = []
leafes_tgt = []
fitch_rs_edges = []
fitch_true_edges = []
fraction_of_xenologs = []
Edges_rs_true = []
Edges_cd_true = []

ind = 0
for item in enumerate(parameter_Df.ID):
    path_s = str(item[1]) + '_species_tree.pickle'
    path_tgt = str(item[1]) + '_gene_tree.pickle'
    s = PhyloTree.load(Path(wk_dir / '01_Data' / path_s))
    tgt = PhyloTree.load(Path(wk_dir / '01_Data' / path_tgt))
    ogt = te.observable_tree(tgt)
    ldt = hgt.ldt_graph(ogt, s)

    constructor = hgt.RsScenarioConstructor(ldt)
    if constructor.run():
        S_rs = constructor.S
        T_rs = constructor.T
    transfer_edges_rs = hgt.rs_transfer_edges(T_rs, S_rs)
    fitch_rs = hgt.undirected_fitch(T_rs, transfer_edges_rs)
    fitch_rs_edges.append(fitch_rs.number_of_edges())

    cotree = Cotree.cotree(ldt)
    cotree_compl = cotree.complement(inplace=False)
    cd_list = gf.cluster_deletion(cotree_compl)
    fitch_cd = gf.build_graph(cd_list)

    transfer_edges_true = hgt.true_transfer_edges(ogt)
    fitch_true = hgt.undirected_fitch(ogt, transfer_edges_true)
    fitch_true_edges.append(fitch_true.number_of_edges())

    set_rs = set(fitch_rs.edges())
    set_cd = set(fitch_cd.edges())
    set_true = set(fitch_true.edges())
    tuple_list1 = set_true - set_cd
    tuple_list2 = set_true - set_rs
    change_tupel1 = gf.change_tupel(tuple_list1)
    change_tupel2 = gf.change_tupel(tuple_list2)
    Edges_cd_true.append(len((set_cd - change_tupel1) - set_true))
    Edges_rs_true.append(len((set_rs - change_tupel2) - set_true))

    edges_ldt.append(len(ldt.edges()))
    leafes_s.append(s.number_of_species)
    leafes_tgt.append(len(tgt.color_sorted_leaves()))
    ind += 1

    print(ind)

a = np.array(edges_ldt, dtype = np.float64)
b = np.array(fitch_true_edges, dtype = np.float64)
parameter_Df['LDT_Edges'] = edges_ldt
parameter_Df['Fitch_true_Edges'] = fitch_true_edges
parameter_Df['Fitch_rs_Edges'] = fitch_rs_edges
parameter_Df['Fraction_of_Xenologs'] = np.divide(a, b, out = np.zeros_like(a), where=b != 0)
parameter_Df['Number_of_Species'] = leafes_s
parameter_Df['Number_of_leaves_tgt'] = leafes_tgt
parameter_Df['Edges_cd_true'] = Edges_cd_true
parameter_Df['Edges_rs_true'] = Edges_rs_true
parameter_Df.to_csv(Path(wk_dir / 'Tree_data.csv', index=False))



