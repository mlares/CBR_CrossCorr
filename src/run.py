
# %load_ext autoreload
# %autoreload 2
# para que ande el run recargando el modulo desde ipython

import PixelSky as pxs
import pandas as pd
import healpy as hp
import numpy as np
from astropy import units as u
import time


# para parsear los parametros de un experimento se pueden usar
# diferentes herramientas:
# https://docs.python.org/3/library/configparser.html
# https://tomassetti.me/parsing-in-python/#tools

#_________________________________________
# Parse parameters from configuration file

import configparser

config = configparser.ConfigParser()
config.read('../set/config.ini')


#_________________________________________
# Read CMB temperature map and mask

nside = int(config['maps']['filedata_cmb_nside'])
mapa = pxs.SkyMap(nside)
mask = pxs.SkyMap(nside)

# filename = '../dat/lensmap512_10arcmin_y2.fits'
# mapa.load(filename, field=(0))
# filename = '../dat/lensmask512_10arcmin_y2.fits'
# mask.load(filename, field=(0))


filedata = config['maps']['datadir_cmb']+ config['maps']['filedata_cmb_mapa']
mapa.load(filedata, field=( int(config['maps']['filedata_field_mapa']) ))

filedata = config['maps']['datadir_cmb']+ config['maps']['filedata_cmb_mask']
mask.load(filedata, field=( int(config['maps']['filedata_field_mask']) ))


#________________________________
# Galaxy catalog

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

rp = pxs.RadialProfile()
rp.set_breaks(unit=u.arcmin, start=0., stop=100., num=5)
res = rp.radialprofile_II(centers, mapa, mask, 15)
rp.signal = np.mean(res, 1)
rp.sigma = np.std(res, 1)



#________________________________
# run control sample (CS)
#
#N = len(glx)
#x = np.random.random((N,3))
#
#v = []
#for i in x:
#    v.append(np.dot(i,i))
#
#y = []
#for i in range(N):
#    y.append(np.division(x[i], v[i]))
#
#glx['vec'] = np.array(y)

#def sample_spherical(npoints, ndim=3):
#    vec = np.random.randn(ndim, npoints)
#    vec /= np.linalg.norm(vec, axis=0)
#    return vec

#CS = pxs.RadialProfile()
#CS.set_breaks(unit=u.arcmin, start=0., stop=30., num=5)
#
#res = CS.radialprofile_II(centers, mapa, mask)
#
#CS.signal = np.mean(res, 1)
#CS.sigma = np.std(res, 1)
 

#________________________________
# Save results
import pickle
         
if config['out']['save_pickle']:
    filedata = config['out']['output_dir']+\
               config['out']['pickle_name_root']+\
               config['out']['pickle_name_exp']+\
               config['out']['pickle_name_idx']+'.p'
     
    pickle.dump( rp, open( filedata, "wb" ) )

