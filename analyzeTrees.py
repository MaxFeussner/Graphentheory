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

# %% load trees
#s = PhyloTree.load(Path(wk_dir / '01_Data' / 'P0_0000_species_tree.pickle'))
#tgt = PhyloTree.load(Path(wk_dir / '01_Data' / 'P0_0000_gene_tree.pickle'))

s = te.simulate_species_tree(50,
                             model='innovation',
                             non_binary_prob=0.0,
                             planted=True,
                             remove_extinct=False,
                             rescale_to_height=1.0,
                             )

# true gene tree (contains losses) of type ’PhyloTree’
tgt = te.simulate_dated_gene_tree(s,
                                  dupl_rate=0.0,
                                  loss_rate=0.3,
                                  hgt_rate=0.5,
                                  dupl_polytomy=0.0,
                                  prohibit_extinction='per_species',
                                  replace_prob=1.0
                                  )
ogt = te.observable_tree(tgt)

# LDT and Fitch Graph
ldt = hgt.ldt_graph(ogt, s)
transfer_edges = hgt.rs_transfer_edges(ogt, s)
transfer_edges2 = hgt.true_transfer_edges(ogt)
fitch_true = hgt.undirected_fitch(ogt, transfer_edges)

# Build Cotree
cotree = Cotree.cotree(ldt)
cotree_compl = cotree.complement(inplace=False)


def cluster_deletion(tree):
    root = tree.root
    cliques = get_cliques(root)
    return cliques


def get_cliques(node):
    if node.label == "leaf":
        #liste mit liste
        q = [[node.ID]]
        return q
    # falls t(u) = 1 : merge Qs, in denen die Kinder von u sind
    elif node.label == "series":
        q_neu = [[]]
        for child in node.children:
            for i, child_cliques in enumerate(get_cliques(child)):
                if len(q_neu) > i:
                    q_neu[i].extend(child_cliques)
                else:
                    q_neu.append(child_cliques)
        return q_neu
    # falls t(u) = 0 : append Qs aneinander (nicht mergen), in  denen Kinder von u sind
    elif node.label == "parallel":
        cliques = []
        for child in node.children:
            #extend
            cliques.extend(get_cliques(child))
            #der groesse nach absteigend sortieren
        cliques.sort(key=len, reverse = True)
        return cliques
    else:
        print("Error - node ")


test1 = cluster_deletion(cotree_compl)

print(test1)


def build_graph(list_of_list):
    netx_graph = netx.Graph()
    for i in range(0, len(list_of_list)-1):
        counter = 1
        while counter <= (len(list_of_list)-1-i):
            for n in range(0, len(list_of_list[i])):
                for x in range(0, len(list_of_list[i+counter])):
                    print(list_of_list[i][n], list_of_list[i+counter][x])
                    netx_graph.add_edge(list_of_list[i][n], list_of_list[i+counter][x])
            counter += 1
    return netx_graph


fitch_cd = build_graph(test1)
fitch_cd.add_edge('a', 'b')
netx.draw(fitch_cd)
plt.savefig("test_graph")
test2 = fitch_cd.edges()
test3 = fitch_true.edges()
set_test2 = set(test2)
set_test3 = set(test3)
tuple_list = set_test3 - set_test2

changed_tuple_list = set()
for x in tuple_list:
    changed_tuple_list.add(((x[1]), (x[0])))

result = (set_test2 - changed_tuple_list) - set_test3

print("frisch")


if 0 == 1:
    # %% Iterate over all data
    for item in enumerate(parameter_Df.ID):
        s = PhyloTree.load(Path(wk_dir / '01_Data' / str(item + '_species_tree.pickle')))
        tgt = PhyloTree.load(Path(wk_dir / '01_Data' / str(item + '_gene_tree.pickle' )))




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
