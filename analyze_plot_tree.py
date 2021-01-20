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

# %% Iterate over all data
ind = 0
for index, item in enumerate(parameter_Df.ID):
    print('Working on Tree # ', ind)
    path_s = str(item) + '_species_tree.pickle'
    path_tgt = str(item) + '_gene_tree.pickle'
    # load data
    s = PhyloTree.load(Path(wk_dir / '01_Data' / path_s))
    tgt = PhyloTree.load(Path(wk_dir / '01_Data' / path_tgt))
    # create graphs
    ogt = te.observable_tree(tgt)
    ldt = hgt.ldt_graph(ogt, s)

    # calculate some interesting parameters and store them in the dataframe
    transfer_edges_true = hgt.true_transfer_edges(ogt)
    fitch_true = hgt.undirected_fitch(ogt, transfer_edges_true)

    parameter_Df.loc[index, ('LDT_Edges')] = len(ldt.edges())
    parameter_Df.loc[index, ('Fitch_true_Edges')] = fitch_true.number_of_edges()

    a = np.array(len(ldt.edges()), dtype = np.float64)
    b = np.array(fitch_true.number_of_edges(), dtype = np.float64)

    parameter_Df.loc[index, ('Fraction_of_Xenologs')] = np.divide(a, b, out = np.zeros_like(a), where=b != 0)
    parameter_Df.loc[index, ('Number_of_Species')] = s.number_of_species
    parameter_Df.loc[index, ('Number_of_leaves_tgt')] = len(tgt.color_sorted_leaves())

    # create triples
    triples_T = set(ogt.get_triples(id_only=True))
    triples_S = set(s.get_triples(id_only=True))
    triple_ldt = set(gf.get_ldt_triples(ldt))
    triple_ldt_color = set(gf.get_ldt_triple_color(ldt))

    parameter_Df.loc[index, ('T_ldt_false_positive')] = len(gf.false_positive(triples_T, triple_ldt))
    parameter_Df.loc[index, ('T_ldt_true_positive')] = len(gf.true_positive(triples_T, triple_ldt))
    parameter_Df.loc[index, ('T_ldt_false_negative')] = len(gf.false_negative(triples_T, triple_ldt))
    parameter_Df.loc[index, ('S_ldt_false_positive')] = len(gf.false_positive(triples_S, triple_ldt_color))
    parameter_Df.loc[index, ('S_ldt_true_positive')] = len(gf.true_positive(triples_S, triple_ldt_color))
    parameter_Df.loc[index, ('S_ldt_false_negative')] = len(gf.false_negative(triples_S, triple_ldt_color))

    # %% Create subgraphs
    for percs in [1, 0.8, 0.6, 0.4, 0.2]:
        print('Subgraph: ' + str(percs * 100))
        # Generate Subgraphs
        ldtSub, fitch_trueSub = gf.buildSubgraph(ldt, fitch_true, percs)

        constructor = hgt.RsScenarioConstructor(ldt)

        if constructor.run():
            S_rs = constructor.S
            T_rs = constructor.T

        transfer_edges_rs = hgt.rs_transfer_edges(T_rs, S_rs)
        fitch_rs = hgt.undirected_fitch(T_rs, transfer_edges_rs)

        cotree = Cotree.cotree(ldt)
        cotree_compl = cotree.complement(inplace=False)
        cd_list = gf.cluster_deletion(cotree_compl)
        fitch_cd = gf.build_graph(cd_list)

        set_rs = set(fitch_rs.edges())
        set_cd = set(fitch_cd.edges())
        set_true = set(fitch_true.edges())

        parameter_Df.loc[index, ('Edges_rs_false_positive_' + str(percs * 100))] = len(gf.false_positive(set_true, set_rs))
        parameter_Df.loc[index, ('Edges_rs_false_negative_' + str(percs * 100))] = len(gf.false_negative(set_true, set_rs))
        parameter_Df.loc[index, ('Edges_rs_true_positive_' + str(percs * 100))] = len(gf.true_positive(set_true, set_rs))
        parameter_Df.loc[index, ('Fitch_rs_Edges_' + str(percs * 100))] = fitch_rs.number_of_edges()

        parameter_Df.loc[index, ('Edges_cd_false_positive_' + str(percs * 100))] = len(gf.false_positive(set_true, set_cd))
        parameter_Df.loc[index, ('Edges_cd_false_negative_' + str(percs * 100))] = len(gf.false_negative(set_true, set_cd))
        parameter_Df.loc[index, ('Edges_cd_true_positive_' + str(percs * 100))] = len(gf.true_positive(set_true, set_cd))

    ind += 1
    print('Done. .  .')
print('Realy Done!')

# %% Save Parameter as csv
parameter_Df.to_csv(Path(wk_dir / 'Tree_data.csv', index=False))
print('I also saved the results. .  .')






