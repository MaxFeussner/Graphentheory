#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 15:05:49 2021

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
import itertools as it
#import graphFunctions.py



s = te.simulate_species_tree(5,
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


trippleList = []

for tripples in it.combinations(ldt.nodes, 3):
    edgeCount = 0
    tempTrippleList = []
    
    for nodes in it.combinations(tripples, 2):

        if ldt.has_edge(nodes[0], nodes[1]):
            
            tempTrippleList = [nodes[0], nodes[1]]
            tempTrippleList.sort
            tempTrippleList.extend(set(nodes) - set(tripples))
            edgeCount += 1
            
    trippleList.append(tempTrippleList)


print(trippleList)
            


            

# (10, 13, 11)
    
# (3,6,2)           
#
            

            
            
            
            
            
            
            