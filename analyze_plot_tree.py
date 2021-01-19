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


# %% Function for 2.2

def change_tupel(tuple_list):
    changed_tuple_list = set()
    for x in tuple_list:
        changed_tuple_list.add(((x[1]), (x[0])))
    return changed_tuple_list


def build_graph(list_of_list):
    netx_graph = netx.Graph()
    for i in range(0, len(list_of_list)-1):
        counter = 1
        while counter <= (len(list_of_list)-1-i):
            for n in range(0, len(list_of_list[i])):
                for x in range(0, len(list_of_list[i+counter])):
                    #print(list_of_list[i][n], list_of_list[i+counter][x])
                    netx_graph.add_edge(list_of_list[i][n], list_of_list[i+counter][x])
            counter += 1
    return netx_graph


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
    cd_list = cluster_deletion(cotree_compl)
    fitch_cd = build_graph(cd_list)

    transfer_edges_true = hgt.true_transfer_edges(ogt)
    fitch_true = hgt.undirected_fitch(ogt, transfer_edges_true)
    fitch_true_edges.append(fitch_true.number_of_edges())

    set_rs = set(fitch_rs.edges())
    set_cd = set(fitch_cd.edges())
    set_true = set(fitch_true.edges())
    tuple_list1 = set_true - set_cd
    tuple_list2 = set_true - set_rs
    change_tupel1 = change_tupel(tuple_list1)
    change_tupel2 = change_tupel(tuple_list2)
    Edges_cd_true.append(len((set_cd - change_tupel1) - set_true))
    Edges_rs_true.append(len((set_rs - change_tupel2) - set_true))

    edges_ldt.append(len(ldt.edges()))
    leafes_s.append(s.number_of_species)
    leafes_tgt.append(len(tgt.color_sorted_leaves()))
    ind += 1

    print(ind)

a = np.array(edges_ldt)
b = np.array(fitch_true_edges)
parameter_Df['LDT_Edges'] = edges_ldt
parameter_Df['Fitch_true_Edges'] = fitch_true_edges
parameter_Df['Fitch_rs_Edges'] = fitch_rs_edges
parameter_Df['Fraction_of_Xenologs'] = np.divide(a, b, out = np.zeros_like(a), where=b != 0)
parameter_Df['Number_of_Species'] = leafes_s
parameter_Df['Number_of_leaves_tgt'] = leafes_tgt
parameter_Df['Edges_cd_true'] = Edges_cd_true
parameter_Df['Edges_rs_true'] = Edges_rs_true
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
