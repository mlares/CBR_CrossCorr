
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

#_________________________________________
# Parse parameters from configuration file
import configparser
config = configparser.ConfigParser()
config.read('../set/config_small.ini')
 
nside = int(config['maps']['filedata_cmb_nside'])
Nran = int(config['run']['Nran'])

#_________________________________________
# Read CMB temperature map and mask
nside = int(config['maps']['filedata_cmb_nside'])
skymap = pxs.SkyMap(nside)
skymask = pxs.SkyMap(nside)

filedata = config['maps']['datadir_cmb']+ config['maps']['filedata_cmb_mapa']
skymap.load(filedata, field=( int(config['maps']['filedata_field_mapa']) ))

filedata = config['maps']['datadir_cmb']+ config['maps']['filedata_cmb_mask']
skymask.load(filedata, field=( int(config['maps']['filedata_field_mask']) ))
 
#________________________________
# Correlation

corr = pxs.Correlation()
nbins = int(config['run']['corr_n_bins']) 
start = float(config['run']['corr_start']) 
stop = float(config['run']['corr_stop']) 
corr.set_breaks(unit=u.arcmin, start=start, stop=stop, num=nbins+1)

njobs = int(config['run']['n_jobs']) 

# serial run test
res = corr.correlation_stack(skymap, skymask, nside, Nran)

# multithreading run
res = corr.correlation_II(range(100), skymap, skymask, nside, Nran, njobs)



print(res)
