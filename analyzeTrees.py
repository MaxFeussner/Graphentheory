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
s = PhyloTree.load(Path(wk_dir / '01_Data' / 'P0_0000_species_tree.pickle'))
tgt = PhyloTree.load(Path(wk_dir / '01_Data' / 'P0_0000_gene_tree.pickle'))
ogt = te.observable_tree(tgt)

# LDT and Fitch Graph
ldt = hgt.ldt_graph(ogt, s)
transfer_edges = hgt.rs_transfer_edges(ogt, s)
transfer_edges2 = hgt.true_transfer_edges(ogt)
fitch = hgt.undirected_fitch(ogt, transfer_edges)

# %% 2.1 Dependece between size of gene tree and numer of species and genes 

#create cotree

cotree = Cotree.cotree(ldt)
cotree_compl =  cotree.complement(inplace=False)



# iterate in postorder over the tree



cliques = []
for node in cotree_compl.postorder():
     if node.label == 'leaf':
         # erstelle liste
         print('leaf')
     elif node.label == 'parallel':
         print('parallel')
         # 
     elif node.label == 'series':
         print('series')









def create_clique(item):
    if item.label == 'leaf':
        #cliques += set(node)
        x = 1
    elif item.label == 'parallel':
        create_clique(item.children)
 
            
print('Child', node.children.postorder)
















# %% Iterate over all data
for item in enumerate(parameter_Df.ID):
    s = PhyloTree.load(Path(wk_dir / '01_Data' / str(item + '_species_tree.pickle')))
    tgt = PhyloTree.load(Path(wk_dir / '01_Data' / str(item + '_gene_tree.pickle' )))
                       

counter_total = 0
counter_non_loss = 0

for i in s.leaves():
    counter_total +=1
    if not i.is_loss():
        counter_non_loss +=1


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
