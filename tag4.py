# -*- coding: utf-8 -*-
"""
Created on Thu Jan 14 11:09:37 2021

@author: jessi
"""
#Cluster Deletion
from  asymmetree  import  PhyloTree
import  asymmetree.treeevolve  as te
import  asymmetree.hgt as hgt
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from asymmetree.cograph import Cotree

wrk_dir = "C:/Users/jessi/.spyder-py3/PraktikumGT/"
dfl = pd.read_csv(wrk_dir + "parameter.csv")
dfID = pd.DataFrame(dfl)
IDs = dfID[dfID.columns[1]]

perc = [1, 0.8, 0.6, 0.4, 0.2]
a=1
#S, TGT aus Simulation zuordnen; für jede Simulation
for i, row in dfID.iterrows():
    ID = row[1]
    S = PhyloTree.load(wrk_dir + "simulation/S_" + ID + ".pickle")
    TGT = PhyloTree.load(wrk_dir + "simulation/TGT_" + ID + ".pickle")

# observable  gene  tree
    OGT = te.observable_tree(TGT)

# LDT laden
    for p in perc:
        ID = str(row[1]) + "_" + str(int(100*p))
        ldt = pd.read_pickle(wrk_dir +"LDT/LDT_"  + ID + ".pickle")



#LDT Graph -> Cotree
        ldtcotree = Cotree.cotree(ldt)

#LDT-Cotree -> Complement
        compldtcotree = Cotree.complement(ldtcotree, inplace=False)


        cliques = []
        for u in compldtcotree.postorder():
            # print("Label:")
            # print(u.label)
        # falls Blatt: Menge i ist {Blatt}
            #if u.is_leaf():
            if u.label == "leaf":
                q = [u]
                cliques.append(q)
        # falls t(u) = 1 : merge Qs, in denen die Kinder von u sind
            elif u.label == "series":
            #elif compldtcograph.nodes[u]["label"] == "series" : #t(v) = 1
                q_neu = []
                cliques_old = cliques
                for child in u.children:
                    #print("Cliques:")
                    #print(cliques)
                    for q in cliques_old:                            
                        #print(q_neu)
                        if [child] in q:
                            q_neu = q_neu.insert(0, q)
                            print([q])
                            print(cliques_old)
                            print(cliques)
                            cliques = cliques.remove(q)
                            print("neu Cliques:")
                            print(cliques)
                if cliques is None:
                    cliques = []
                #sollte vllt insert sein, damit groesstes q an erster stelle
                cliques.append(q_neu)
                #print("Cliques:")
                #print(cliques)
        # falls t(u) = 0 : append Qs aneinander (nicht mergen), in  denen Kinder von u sind
            elif u.label == "parallel":
        #   elif compldtcograph.nodes[u]["label"] == "parallel": #t(v) = 0
                print("TODO")
            else:
                print("Error")
        break
            