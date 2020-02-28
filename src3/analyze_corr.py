
import PixelSky as pxs
import pandas as pd
import healpy as hp
import numpy as np
import random
from astropy import units as u
import time
import pickle

#_________________________________________
# Parse parameters from configuration file
import configparser
import sys
config = configparser.ConfigParser()
if len(sys.argv) > 1:    
    filename = sys.argv[1]
else:
    filename = '../set/config_small.ini'
config.read(filename)
 
nside = int(config['maps']['filedata_cmb_nside'])
Nran = int(config['run']['Nran'])
njobs = int(config['run']['n_jobs']) 
Nexperiments = int(config['run']['Nexperiments'])

nbins = int(config['run']['corr_n_bins']) 
start = float(config['run']['corr_start']) 
stop = float(config['run']['corr_stop']) 
 

corr = pxs.Correlation(nside=nside, nran=1000, skymask=skymask, njobs=njobs)
corr.set_breaks(start=start, stop=stop, num=nbins+1)


f_root = config['out']['output_dir']
f_exp =  config['out']['pickle_name_root']
f_num =  config['out']['pickle_name_exp']
f_idx =  config['out']['pickle_name_idx'] 
f_ext =  config['out']['pickle_name_ext'] 

fname = f_root + f_exp + f_num + f_idx + f_ext

with open( fname, "rb" ) as f:
    res = pickle.load(f)


k = len(res)
n = len(res[0][0])

s_ac = np.zeros(n)
s_tt = np.zeros(n)

for r in res:
    s_ac = s_ac + r[1]
    s_tt = s_tt + r[0]

corr = s_tt / s_ac


#import matplotlib.pyplot as plt
#plt.plot(s_tt)



