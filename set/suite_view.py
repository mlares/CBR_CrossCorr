# load_ext autoreload
# autoreload 2

# experiment, parsing and math
import cmfg
from Parser import Parser
from sys import argv
import numpy as np
import pickle
import math as m
import pandas as pd

ind_min = 12
ind_max = 22

inifile = '../set/config_TRNT_' + str(ind_min).zfill(3) + '.ini'
config = Parser(inifile)

cnames = config.p._fields
dfa = pd.DataFrame(columns=cnames)

for i, ids in enumerate(range(ind_min, ind_max+1)):
    inifile = '../set/config_TRNT_' + str(ids).zfill(3) + '.ini'
    print(inifile)
    config = Parser(inifile)
    dfa.loc[i] = config.p

sacar = ['n_jobs', 'verbose', 'run_parallel', 'showp', 'overwrite',
         'dir_output', 'dir_plots']
df2 = df = dfa.drop(columns=sacar)
 
df2.to_excel('../set/settings_TRNT.xlsx')