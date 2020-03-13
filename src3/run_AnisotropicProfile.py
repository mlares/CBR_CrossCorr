
# %load_ext autoreload
# %autoreload 2
# para que ande el run recargando el modulo desde ipython

# python run_profile.py ../set/config_small.ini

import PixelSky as pxs

###############################
import imp
pxs = imp.reload(pxs)
###############################


import pandas as pd
import healpy as hp
import numpy as np
from astropy import units as u
import time
import sys
from configparser import ConfigParser

#_________________________________________
# Parse parameters from configuration file

#filename = pxs.check_file(sys.argv)
filename = pxs.check_file(['', '../set/config_ani.ini'])
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
#glx['vec'] = hp.ang2vec(theta_healpix, phi_healpix).tolist()

glx['pa'] = len(glx['RAdeg'])* [1.]  #### HAY QUE LEERLO

# filter galaxy catalog...
# l = ['A','X','B']
# spiral = [any([s in x for s in l]) for x in glx['type']]
# edgeon = glx['b/a'] < 0.8
# subset = spiral & edgeon
# centers = np.array(list(glx.vec[subset]))

centers = pd.DataFrame(list(zip(phi_healpix, theta_healpix, glx['pa'])),
          columns=['phi','theta', 'pa'] )


#________________________________
# Radial profile


#================
print('Configuring radial profile parameters...')
ap = pxs.AnisotropicProfile()

nbins = int(config['run']['rp_n_bins']) 
start = float(config['run']['rp_start']) 
stop = float(config['run']['rp_stop']) 

ap.set_breaks_radial(unit=u.arcmin, start=start, stop=stop, num=nbins+1)

nbins = int(config['run']['theta_n_bins']) 
start = float(config['run']['theta_start']) 
stop = float(config['run']['theta_stop']) 

ap.set_breaks_angular(unit=u.rad, start=start, stop=stop, num=nbins+1)
#================

#res = ap.anisotropic_profile(centers[0:1], mapa, mask)
dists, thetas, temps, bins2d = ap.anisotropic_profile(centers[0:1], mapa, mask)




#max_centers = 10
#res = ap.anisotropic_profile_II(centers[:max_centers], ap[:max_centers], mapa, mask, njobs)





# njobs = int(config['run']['n_jobs']) 
# max_centers = int(config['cats']['max_centers']) 
# ncent = min(max_centers, len(centers))
# 
# if use_parallel:
#     res = rp.radialprofile_II(centers[:max_centers], mapa, mask, njobs)
# else:
#     res = rp.radialprofile(centers[:max_centers], mapa, mask, njobs)
# 
# """
# TEST:
#     probar como escala y como trabajan las versiones _II y normal
# """
# 
# rp.signal = np.mean(res, 1)
# rp.sigma = np.std(res, 1)
# 
# 
# #________________________________
# # Save results
# import pickle
#          
# if config['out']['save_pickle']:
#     filedata = config['out']['output_dir']+\
#                config['out']['pickle_name_root']+\
#                config['out']['pickle_name_exp']+\
#                config['out']['pickle_name_idx']+'.p'
#      
#     pickle.dump( rp, open( filedata, "wb" ) )



