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
import asymmetree.tools.GraphTools as gt
import graphFunctions as gf


s = te.simulate_species_tree(10,
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
transfer_edges_true = hgt.true_transfer_edges(ogt)
fitch_true = hgt.undirected_fitch(ogt, transfer_edges_true)

cotree = Cotree.cotree(ldt)
cotree_compl = cotree.complement(inplace=False)
cd_list = gf.cluster_deletion(cotree_compl)
fitch_cd = gf.build_graph(cd_list)
triples_T = ogt.get_triples(id_only=True)


def get_ldt_tripples(ldt):
    trippleList = []
    for tripples in it.combinations(ldt.nodes, 3):
        edgeCount = 0
        tempTrippleList = []
        for nodes in it.combinations(tripples, 2):
            if ldt.has_edge(nodes[0], nodes[1]):
                tempTrippleList = [nodes[0], nodes[1]]
                tempTrippleList.sort()
                tempTrippleList.extend(list({tripples[0], tripples[1], tripples[2]} - {nodes[0], nodes[1]}))
                edgeCount += 1
        if edgeCount == 1:
            trippleList.append(tempTrippleList[0])
    return trippleList


def get_ldt_tripple_color(ldt):
    trippleList = []
    for tripples in it.combinations(ldt.nodes, 3):
        if ldt.has_edge(tripples[0], tripples[1]) and ldt.has_edge(tripples[2], tripples[1]) and not ldt.has_edge(tripples[0], tripples[2]):
            if ldt.nodes[tripples[0]]["color"] != ldt.nodes[tripples[2]]["color"]:
                tempTrippleList = [tripples[0], tripples[2]]
                tempTrippleList.sort()
                tempTrippleList.extend([tripples[1]])
                trippleList.append(tripples)
        elif ldt.has_edge(tripples[0], tripples[2]) and ldt.has_edge(tripples[1], tripples[2]) and not ldt.has_edge(tripples[0], tripples[1]):
            if ldt.nodes[tripples[0]]["color"] != ldt.nodes[tripples[1]]["color"]:
                tempTrippleList = [tripples[0], tripples[1]]
                tempTrippleList.sort()
                tempTrippleList.extend([tripples[2]])
                trippleList.append(tripples)
        elif ldt.has_edge(tripples[1], tripples[0]) and ldt.has_edge(tripples[2], tripples[0]) and not ldt.has_edge(tripples[1], tripples[2]):
            if ldt.nodes[tripples[1]]["color"] != ldt.nodes[tripples[2]]["color"]:
                tempTrippleList = [tripples[1], tripples[2]]
                tempTrippleList.sort()
                tempTrippleList.extend([tripples[0]])
                trippleList.append(tripples)
    return trippleList


def sort_tripple(tripple_list):
    changed_tripple_list = []
    for tripples in tripple_list:
        tempTrippleList = [tripples[0], tripples[1]]
        tempTrippleList.sort()
        tempTrippleList.extend([tripples[2]])
        changed_tripple_list.append((tempTrippleList[0], tempTrippleList[1], tempTrippleList[2]))
    return changed_tripple_list


ct = gt.contingency_table(fitch_true, fitch_cd)
print(ct['tp'], ct['fp'])
print(gt.graphs_equal(fitch_true, fitch_cd))

#test = get_ldt_tripples(ldt)
#print(triples_T)
#print(sort_tripple(triples_T))

            

# (10, 13, 11)
    
# (3,6,2)           
#
            

            
            
            
            
            
            
            