
# %load_ext autoreload
# %autoreload 2
# para que ande el run recargando el modulo desde ipython

# python run_profile.py ../set/config_small.ini

import PixelSky as pxs
import pandas as pd
import healpy as hp
import numpy as np
from astropy import units as u
import time
import sys
from configparser import ConfigParser

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
# Galaxy catalog

print('Reading list of centers...')
# read...
glx_catalog = config['cats']['datadir_glx']+config['cats']['filedata_glx']
glx = pd.read_csv(glx_catalog, delim_whitespace=True, header=9)
# catalog: http://tdc-www.harvard.edu/2mrs/2mrs_readme.html

# positions...
phi_healpix = glx['RAdeg']*np.pi/180.
theta_healpix = (90. - glx['DECdeg'])*np.pi/180.
glx['vec'] = hp.ang2vec(theta_healpix, phi_healpix).tolist()

# filter galaxy catalog...
# l = ['A','X','B']
# spiral = [any([s in x for s in l]) for x in glx['type']]
# edgeon = glx['b/a'] < 0.8
# subset = spiral & edgeon
# centers = np.array(list(glx.vec[subset]))

centers = np.array(list(glx.vec))

#________________________________
# Radial profile

print('Configuring radial profile parameters...')
rp = pxs.RadialProfile()
nbins = int(config['run']['rp_n_bins']) 
start = float(config['run']['rp_start']) 
stop = float(config['run']['rp_stop']) 
rp.set_breaks(unit=u.arcmin, start=start, stop=stop, num=nbins+1)

njobs = int(config['run']['n_jobs']) 
max_centers = int(config['cats']['max_centers']) 
ncent = min(max_centers, len(centers))

use_parallel = config['run']['use_parallel'].lower() in\
               ['true', '1', 't', 'y', 'yes', 'yeah', 'yup',\
                       'certainly', 'uh-huh', 's', 'si', 'sip']
if use_parallel:
    res = rp.radialprofile_II(centers[:max_centers], mapa, mask, njobs)
else:
    res = rp.radialprofile(centers[:max_centers], mapa, mask, njobs)

#"""
#TEST:
#    probar como escala y como trabajan las versiones _II y normal
#"""
#
#rp.signal = np.mean(res, 1)
#rp.sigma = np.std(res, 1)


##________________________________
## Save results
#import pickle
#         
#if config['out']['save_pickle']:
#    filedata = config['out']['output_dir']+\
#               config['out']['pickle_name_root']+\
#               config['out']['pickle_name_exp']+\
#               config['out']['pickle_name_idx']+'.p'
#     
#    pickle.dump( rp, open( filedata, "wb" ) )









 



#  #________________________________
#  # run control sample (CS)
#  #
#  #N = len(glx)
#  #x = np.random.random((N,3))
#  #
#  #v = []
#  #for i in x:
#  #    v.append(np.dot(i,i))
#  #
#  #y = []
#  #for i in range(N):
#  #    y.append(np.division(x[i], v[i]))
#  #
#  #glx['vec'] = np.array(y)
#  
#  #def sample_spherical(npoints, ndim=3):
#  #    vec = np.random.randn(ndim, npoints)
#  #    vec /= np.linalg.norm(vec, axis=0)
#  #    return vec
#  
#  #CS = pxs.RadialProfile()
#  #CS.set_breaks(unit=u.arcmin, start=0., stop=30., num=5)
#  #
#  #res = CS.radialprofile_II(centers, mapa, mask)
#  #
#  #CS.signal = np.mean(res, 1)
#  #CS.sigma = np.std(res, 1)
