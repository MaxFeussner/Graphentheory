#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 19 15:11:24 2021

@author: Tobias Meißner
"""

import numpy as np
import networkx as netx
import itertools as it



def change_tupel(tuple_list: list) -> list:
    '''
    Change the order of tupels.

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


def build_graph(list_of_list: list):
    '''
    Builds a NetworkX from a list of Nodes

    Parameters
    ----------
    list_of_list : List
        A List with lists of nodes.

    Returns
    -------
    netx_graph : NetworkX.Graph
        A complete NetworkX.Graph.

    '''
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
    '''
    Creates List of Tupels with tree ID's with cliques

    Parameters
    ----------
    node : phylotree.node
        A Phylotree node.

    Returns
    -------
    TYPE
        A list of Tupels with cliques.

    '''
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


def buildSubgraph(ldt, fitch,  prozent: float):
    '''
    Generates a subgraphs with a given percentage of kept nodes

    Parameters
    ----------
    ldt : TYPE
        Graph objekt.
    prozent : float
        Percent to be kept.

    Returns
    -------
    Graph
        A Graph with some removed nodes.

    '''
    
    randomList = np.random.choice(a = list(ldt.nodes()),
                                  size = int((len(list(ldt.nodes())) * prozent)))
    if len(randomList) == 0:
        node_list = list(ldt.nodes())
        randomList = node_list[0]

    ldtSub = ldt.subgraph(randomList)
    fitchSub = fitch.subgraph(randomList)
    
    return (ldtSub, fitchSub)


def get_ldt_triples(ldt):
    tripleList = []
    for triples in it.combinations(ldt.nodes, 3):
        edgeCount = 0
        tempTripleList = []
        for nodes in it.combinations(triples, 2):
            if ldt.has_edge(nodes[0], nodes[1]):
                tempTripleList = [nodes[0], nodes[1]]
                tempTripleList.sort()
                tempTripleList.extend(list({triples[0], triples[1], triples[2]} - {nodes[0], nodes[1]}))
                edgeCount += 1
        if edgeCount == 1:
            tripleList.append((tempTripleList[0], tempTripleList[1], tempTripleList[2]))
    return tripleList


def get_ldt_triple_color(ldt):
    tripleList = []
    for triples in it.combinations(ldt.nodes, 3):
        if ldt.has_edge(triples[0], triples[1]) and ldt.has_edge(triples[2], triples[1]) and not ldt.has_edge(triples[0], triples[2]):
            if ldt.nodes[triples[0]]["color"] != ldt.nodes[triples[2]]["color"]:
                tempTripleList = [triples[0], triples[2]]
                tempTripleList.sort()
                tempTripleList.extend([triples[1]])
                tripleList.append(triples)
        elif ldt.has_edge(triples[0], triples[2]) and ldt.has_edge(triples[1], triples[2]) and not ldt.has_edge(triples[0], triples[1]):
            if ldt.nodes[triples[0]]["color"] != ldt.nodes[triples[1]]["color"]:
                tempTripleList = [triples[0], triples[1]]
                tempTripleList.sort()
                tempTripleList.extend([triples[2]])
                tripleList.append(triples)
        elif ldt.has_edge(triples[1], triples[0]) and ldt.has_edge(triples[2], triples[0]) and not ldt.has_edge(triples[1], triples[2]):
            if ldt.nodes[triples[1]]["color"] != ldt.nodes[triples[2]]["color"]:
                tempTripleList = [triples[1], triples[2]]
                tempTripleList.sort()
                tempTripleList.extend([triples[0]])
                tripleList.append(triples)
    return tripleList


def sort_triple(tripple_list: list) -> list:
    '''
    Sorts the first two elements of a Tupel in ascending order

    Parameters
    ----------
    tripple_list : list
        List of Tupels.

    Returns
    -------
    list
        Sorted lost of tupels.

    '''
    changed_tripple_list = []
    for triples in tripple_list:
        tempTripleList = [triples[0], triples[1]]
        tempTripleList.sort()
        tempTripleList.extend([triples[2]])
        changed_tripple_list.append((tempTripleList[0], tempTripleList[1], tempTripleList[2]))
    return changed_tripple_list


def false_positive(set_true, set_cd):
    
    tuple_list1 = set_true - set_cd
    change_tupel1 = change_tupel(tuple_list1)
    result_cd_fp = (set_cd - change_tupel1) - set_true
    return result_cd_fp


def false_negative(set_true, set_cd):
    tuple_list3 = set_cd - set_true
    change_tupel3 = change_tupel(tuple_list3)
    result_cd_fn = (set_true - change_tupel3) - set_cd
    return result_cd_fn


def true_positive(set_true, set_cd):
    set_true_double = set_true | change_tupel(set_true)
    tuple_list5 = set_true_double & set_cd
    return tuple_list5


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