#!/usr/bin/python/soft/miniconda3/bin/python
# %load_ext autoreload
# %autoreload 2
# para que ande el run recargando el modulo desde ipython

import PixelSky as pxs
import pandas as pd
import healpy as hp
import numpy as np
import random
from astropy import units as u
import time

import os
os.system("taskset -p 0xff %d" % os.getpid())

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

#_________________________________________
# Read CMB temperature map and mask
nside = int(config['maps']['filedata_cmb_nside'])
skymap = pxs.SkyMap(nside)
skymask = pxs.SkyMap(nside)

filedata = config['maps']['datadir_cmb']+ config['maps']['filedata_cmb_mapa']
skymap.load(filedata, field=( int(config['maps']['filedata_field_mapa']) ))

filedata = config['maps']['datadir_cmb']+ config['maps']['filedata_cmb_mask']
skymask.load(filedata, field=( int(config['maps']['filedata_field_mask']) ))
skymask.data = skymask.data.astype(int)

#________________________________
# Correlation

corr = pxs.Correlation(nside=nside, nran=1000, skymask=skymask, njobs=njobs)
corr.set_breaks(start=start, stop=stop, num=nbins+1)

# serial run test
# res = corr.correlation(1, skymap, skymask, nside, Nran)

# multithreading run
result = corr.correlation_II(range(Nexperiments), skymap, skymask)

import pickle
 
f_root = config['out']['output_dir']
f_exp =  config['out']['pickle_name_root']
f_num =  config['out']['pickle_name_exp'] 
f_idx =  config['out']['pickle_name_idx'] 
f_ext =  config['out']['pickle_name_ext'] 

fname = f_root + f_exp + f_num + f_idx + f_ext
 
with open( fname, "wb" ) as f:
    pickle.dump(result, f)


## cf = np.sum(res, 0)/Nexperiments/Nran
