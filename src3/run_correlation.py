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
import sys
from configparser import ConfigParser

import os
os.system("taskset -p 0xff %d" % os.getpid())

#_________________________________________
# Parse parameters from configuration file

filename = pxs.check_file(sys.argv)
config = ConfigParser()
config.read(filename) 

                
#_________________________________________
# Read CMB temperature map and mask

print('Reading CMB map and mask...')

nside = int(config['maps']['filedata_cmb_nside'])
mapa = pxs.SkyMap(nside)
mask = pxs.SkyMap(nside)

filedata = config['maps']['datadir_cmb']+ config['maps']['filedata_cmb_mapa']
mapa.load(filedata, field=( int(config['maps']['filedata_field_mapa']) ))

filedata = config['maps']['datadir_cmb']+ config['maps']['filedata_cmb_mask']
mask.load(filedata, field=( int(config['maps']['filedata_field_mask']) ))

                                    

#________________________________
# Correlation

print('Configuring correlation parameters...')
nbins = int(config['run']['rp_n_bins']) 
start = float(config['run']['rp_start']) 
stop = float(config['run']['rp_stop']) 
njobs = int(config['run']['n_jobs']) 
Nran = int(config['run']['nran']) 


corr = pxs.Correlation(nside=nside, nran=1000, skymask=mask, njobs=njobs)
corr.set_breaks(start=start, stop=stop, num=nbins+1)

# serial run test
# res = corr.correlation(1, mapa, mask) #, nside, Nran)
 
# multithreading run

Nexperiments = 100
result = corr.correlation_II(range(Nexperiments), mapa, mask)
 

#________________________________
# Save results
import pickle
         
if config['out']['save_pickle']:
    filedata = config['out']['output_dir']+\
               config['out']['pickle_name_root']+\
               config['out']['pickle_name_exp']+\
               config['out']['pickle_name_idx']+'.p'
     
    with open( filedata, "wb" ) as f:
        pickle.dump( result, open( filedata, "wb" ) )               



