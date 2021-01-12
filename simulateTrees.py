#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 11 12:44:09 2021

@author: Tobias Meisner

Praktikum: Graphtheorie

"""

import os
import pandas as pd
import asymmetree.treeevolve as te
import numpy as np
from pathlib import Path

# %% Housekeeping
# check where the python script is located to store data in same directory
# adapt the 'else' part if you work interactiveley, it should work from terminal
if '__file__' in vars():
    wk_dir = Path(os.path.dirname(os.path.realpath('__file__')))
else:
    print('We are running the script interactively')
    wk_dir = Path('/Users/TMC/Documents/01_Studium/02_Master/3_Semester/01_Graphentheorie/01_Prkatikum/Script/simulateTrees.py').parent

# create output directory for data
Path(wk_dir / '01_Data').mkdir(parents=True, exist_ok=True)


# %% Define Parameter for Simulation     
# first we crate an empty pandas DataFrame where we save the parameters for 
# every simulated tree

simulation_rate = 10

colNames = ['ID',                                                              # Unique ID :int 
            'num_of_leaves',                                                   # random int range(10,50)                           
            'model',                                                           # innovation / Bird-Death (BDP)
            'non_binary_prob',                                                 # random contraction of inner verteces range(0.0,0.5)
            'planted',
            'remove_extinct',
            'rescale_to_height',                                               # Timespan of the tree 
            'dupl_rate',                                                       # values from script
            'loss_rate',
            'hgt_rate',
            'prohibit_extinction',                                             # avoids extinction 'per_family'
            'replace_prob']                                                    # 1.0

rowIndex = range(0,simulation_rate)

parameter_Df = pd.DataFrame(index = rowIndex, columns=colNames)


# list of parameters: duplication, loss, HGT rate (given in that order)
list_of_parameters = [(0.25, 0.25, 0.25),
                      (0.5, 0.5, 0.5),
                      (0.5, 0.5, 1.0),
                      (0.5, 0.5, 1.5),
                      (1.0, 1.0, 0.5),
                      (1.0, 1.0, 1.0),
                      (1.5, 1.5, 1.5)]

# create a dataframe with N x parameter length and fill it with parameters
ind = 0                                                                        # counter for the row index
for parameter in enumerate(list_of_parameters):
    for iteration in range(simulation_rate):
        parameter_Df.loc[ind, 'ID'] = 'P' + str(parameter[0]) + '_' + str(f'{iteration:04d}')
        parameter_Df.loc[ind, 'num_of_leaves'] = np.random.randint(low=10,high=50,dtype=int)
        parameter_Df.loc[ind, 'model'] = 'innovation'
        parameter_Df.loc[ind, 'non_binary_prob'] = np.random.uniform(low=0.0,high=0.5) 
        parameter_Df.loc[ind, 'planted'] = True 
        parameter_Df.loc[ind, 'remove_extinct'] = False
        parameter_Df.loc[ind, 'rescale_to_height'] = 1.0
        parameter_Df.loc[ind, 'dupl_rate'] = parameter[1][0]
        parameter_Df.loc[ind, 'loss_rate'] = parameter[1][1]
        parameter_Df.loc[ind, 'hgt_rate'] = parameter[1][2]
        parameter_Df.loc[ind, 'prohibit_extinction'] = 'per_family'
        parameter_Df.loc[ind, 'replace_prob'] = 1.0        
        ind += 1 
  
# save the parameter to file
parameter_Df.to_csv(Path(wk_dir / '01_Data') / '01_Simulation_Parameters.csv', index=False)

# load parameter csv
parameter_Df = pd.read_csv(Path(wk_dir / '01_Data') / '01_Simulation_Parameters.csv')

# %% Funktion to search in a pandas dataframe for an entry
def which(self):
    try:
        self = list(iter(self))
    except TypeError as e:
        raise Exception("""'which' method can only be applied to iterables.
        {}""".format(str(e)))
    indices = [i for i, x in enumerate(self) if bool(x) == True]
    return(indices)

# %% Simulation
# Simulate a species tree of type 'PhyoTree'
ind = 0
# build loop

len(parameter_Df.index)

for ind in range(len(parameter_Df.index)-1):
    # species tree of type ’PhyloTree’
    s = te.simulate_species_tree(parameter_Df.loc[ind, 'num_of_leaves'], 
                                 model = parameter_Df.loc[ind, 'model'],
                                 non_binary_prob = parameter_Df.loc[ind, 'non_binary_prob'],
                                 planted = parameter_Df.loc[ind, 'planted'],
                                 remove_extinct = parameter_Df.loc[ind, 'remove_extinct'],
                                 rescale_to_height = parameter_Df.loc[ind, 'rescale_to_height']
                                 )
    
    # true gene tree (contains losses) of type ’PhyloTree’
    tgt = te.simulate_dated_gene_tree(s,
                                      dupl_rate = parameter_Df.loc[ind, 'dupl_rate'],
                                      loss_rate = parameter_Df.loc[ind, 'loss_rate'],
                                      hgt_rate = parameter_Df.loc[ind, 'hgt_rate'],
                                      dupl_polytomy = 0.0,
                                      prohibit_extinction= parameter_Df.loc[ind, 'prohibit_extinction'],
                                      replace_prob = parameter_Df.loc[ind, 'replace_prob']
                                      )
    
    
    # serialization
    s.serialize(wk_dir / '01_Data' / str(parameter_Df.loc[ind, 'ID'] + '_species_tree.pickle'))
    tgt.serialize(wk_dir / '01_Data' / str(parameter_Df.loc[ind, 'ID'] + '_gene_tree.pickle'))

ogt = te.observable_tree(tgt)
