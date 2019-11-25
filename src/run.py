
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


#________________________________
# Read CMB temperature map and mask

nside = 512
mapa = pxs.SkyMap(nside)
mask = pxs.SkyMap(nside)

# small map
# filename = '../dat/lensmap512_10arcmin_y2.fits'
# mapa.load(filename, field=(0))
# filename = '../dat/lensmask512_10arcmin_y2.fits'
# mask.load(filename, field=(0))


# large map
filename = '../dat/COM_CMB_IQU-smica_2048_R3.00_full.fits'
mapa.load(filename, field=(0))
mask.load(filename, field=(3))


#________________________________
# Galaxy catalog

# read...
glx_catalog = '../dat/2mrs_1175_done.dat'
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

pickle.dump( rp, open( "rp_run__tst.p", "wb" ) )
