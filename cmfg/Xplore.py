import cmfg
from Parser import Parser
from sys import argv
import numpy as np
from astropy import units as u
from matplotlib import pyplot as plt
import pandas as pd
from astropy.cosmology import FlatLambdaCDM
import math


if len(argv) > 1:
    config = Parser(argv[1])
else:
    config = Parser()

X = cmfg.profile2d(config)
X.load_centers()
X.load_tracers()

a_r = np.array(X.centers['r_ext'])
plt.hist(a_r)
plt.xlabel('galaxy size in arcsec')
plt.show()

a_rad = np.array(X.centers['glx_size_rad'])
plt.hist(a_rad)
plt.xlabel('galaxy size in radians')
plt.show()

a_arcmin = (a_r*u.rad).to(u.arcmin).value
plt.hist(a_arcmin)
plt.xlabel('galaxy size in arcmins')
plt.show()

a_kpc = np.array(X.centers['glx_size_kpc'])
plt.hist(a_kpc)
plt.xlabel('galaxy size in kpc')
plt.show()

# Pasar de tamaño angular a kpc

c_light = 300000
theta_arcsec = 1.5
redshift = 0.01
H = 68
theta_rad = theta_arcsec/60./60.*math.pi/180.
theta_rad2 = (theta_arcsec*u.arcsec).to(u.rad).value

# metodo 1: (aprox.) rad to kpc
r_Mpc = math.tan(theta_rad/2)*2*c_light*redshift/H
r_kpc = r_Mpc * 1000.
print(r_kpc)

# metodo 2: (aprox.) rad to kpc
cosmo = FlatLambdaCDM(H0=H, Om0=0.308)


d_A = cosmo.angular_diameter_distance(redshift)
r_kpc = (theta_rad * d_A).to(u.kpc, u.dimensionless_angles())
print(r_kpc)

# Propiedades de las galaxias en 0.001 < z < 0.015

glx_catalog = '../dat/2mrs_1175_VAC.dat'
glx = pd.read_csv(glx_catalog, delim_whitespace=True, header=9)

redshift = np.array(glx['v']/300000.)
c_light = 300000
theta_arcsec = np.array(glx['r_ext'])
theta_arcsec = 10**theta_arcsec

H = 68
theta_rad = theta_arcsec/60./60.*math.pi/180.
r_Mpc = np.tan(theta_rad/2)*2*c_light*redshift/H
r_kpc = r_Mpc * 1000.

fig = plt.figure(figsize=(10, 10))
ax = fig.subplots(2, 2)
ax[0, 0].plot(redshift, r_kpc, marker='o', markersize=1., linestyle='None',
              alpha=0.2)
ax[0, 0].set_xlabel('redshift')
ax[0, 0].set_ylabel('galaxy size, kpc')
ax[0, 0].set_xlim(0., 0.015)
ax[0, 0].set_ylim(0., 20)

ax[0, 1].hist(r_kpc, bins=50)
ax[0, 1].set_xlabel('galaxy radius, r [kpc]')
ax[0, 1].set_ylabel('dN/dr')

ax[1, 0].hist(redshift, bins=50)
ax[1, 0].set_xlabel('redshift, z')
ax[1, 0].set_ylabel('dN/dz')
ax[1, 0].set_xlim(0., 0.015)

ax[1, 1].hist(theta_arcsec, bins=50)
ax[1, 1].set_xlabel(r'galaxy angular [arcsec], $\theta$')
ax[1, 1].set_ylabel(r'dN/d$\theta$')

plt.tight_layout()
plt.show()
