import PixelSky as pxs
import pandas as pd
import healpy as hp
import numpy as np
from astropy import units as u

import joblib # import Parallel, delayed

import importlib
pxs = importlib.reload(pxs)

# Read CMB temperature map
nside = 512
mapa = pxs.SkyMap(nside)
mask = pxs.SkyMap(nside)

filename = '../dat/lensmap512_10arcmin_y2.fits'
mapa.load(filename, field=(0))
filename = '../dat/lensmask512_10arcmin_y2.fits'
mask.load(filename, field=(0))

# Read galaxy catalog
glx_catalog = '../dat/2mrs_1175_done.dat'
glx = pd.read_csv(glx_catalog, delim_whitespace=True, header=9)
# catalog: http://tdc-www.harvard.edu/2mrs/2mrs_readme.html

phi_healpix = glx['RAdeg']*np.pi/180.
theta_healpix = (90. - glx['DECdeg'])*np.pi/180.
glx['vec'] = hp.ang2vec(theta_healpix, phi_healpix).tolist()

# filter galaxy catalog
l = ['A','X','B']
spiral = [any([s in x for s in l]) for x in glx['type']]
edgeon = glx['b/a'] < 0.8
subset = spiral & edgeon

rp = pxs.RadialProfile()
rp.set_breaks(unit=u.arcmin, start=0., stop=30., num=10)

#rp.radialprofile(mapa, mask, glx)

# rp.radialprofile_II(glx.vec, mapa, mask)

pos = np.array(list(glx.vec))
rp.run_radialprofile(pos, mapa, mask)



print(rp.signal)
print(rp.sigma)

print('OK')

"""  
def f1(a, b, c=0, *args, **kwargs):
    # Esto es la documentación de la función 1
    partes = ["Parámetros pasados a f1", 
        f"a={a} b={b} args={args} kwargs={kwargs}"]
    return " -- ".join(partes)

class clase():
    def f2(c, d, *args, **kwargs):
        # Esto es la documentación de la función 2
        partes = ["Parámetros pasados a f2", 
            f"c={c} d={d} args={args} kwargs={kwargs}"]
        print(" -- ".join(partes))

        return(f1(c, d, *args, **kwargs))

print(clase.f2(1, 2, 4, l=9))
print('OK')
print('--------------')


# def f(a, b, c=1):
#     print(a+b+c)
""" 

# with Parallel(n_jobs=2) as parallel:
# ...    accumulator = 0.
# ...    n_iter = 0
# ...    while accumulator < 1000:
# ...        results = parallel(delayed(sqrt)(accumulator + i ** 2)
# ...                           for i in range(5))
# ...        accumulator += sum(results)  # synchronization barrier
# ...        n_iter += 1
