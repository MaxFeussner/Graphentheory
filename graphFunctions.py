#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 15:11:24 2021

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


def change_tupel(tuple_list: list)-> list:
    '''
    Change the order of tupels

    Parameters
    ----------
    tuple_list : list
        List of tupels.

    Returns
    -------
    list
        A list with changed tupels.

    '''
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


def get_ldt_triples(ldt):
    tripleList = []
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
            tripleList.append((tempTrippleList[0], tempTrippleList[1], tempTrippleList[2]))
    return tripleList


def get_ldt_triple_color(ldt):
    tripleList = []
    for tripples in it.combinations(ldt.nodes, 3):
        if ldt.has_edge(tripples[0], tripples[1]) and ldt.has_edge(tripples[2], tripples[1]) and not ldt.has_edge(tripples[0], tripples[2]):
            if ldt.nodes[tripples[0]]["color"] != ldt.nodes[tripples[2]]["color"]:
                tempTrippleList = [tripples[0], tripples[2]]
                tempTrippleList.sort()
                tempTrippleList.extend([tripples[1]])
                tripleList.append(tripples)
        elif ldt.has_edge(tripples[0], tripples[2]) and ldt.has_edge(tripples[1], tripples[2]) and not ldt.has_edge(tripples[0], tripples[1]):
            if ldt.nodes[tripples[0]]["color"] != ldt.nodes[tripples[1]]["color"]:
                tempTrippleList = [tripples[0], tripples[1]]
                tempTrippleList.sort()
                tempTrippleList.extend([tripples[2]])
                tripleList.append(tripples)
        elif ldt.has_edge(tripples[1], tripples[0]) and ldt.has_edge(tripples[2], tripples[0]) and not ldt.has_edge(tripples[1], tripples[2]):
            if ldt.nodes[tripples[1]]["color"] != ldt.nodes[tripples[2]]["color"]:
                tempTrippleList = [tripples[1], tripples[2]]
                tempTrippleList.sort()
                tempTrippleList.extend([tripples[0]])
                tripleList.append(tripples)
    return tripleList


def sort_triple(tripple_list):
    changed_tripple_list = []
    for tripples in tripple_list:
        tempTrippleList = [tripples[0], tripples[1]]
        tempTrippleList.sort()
        tempTrippleList.extend([tripples[2]])
        changed_tripple_list.append((tempTrippleList[0], tempTrippleList[1], tempTrippleList[2]))
    return changed_tripple_list

        
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
