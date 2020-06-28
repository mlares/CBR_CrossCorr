
# experiment, parsing and math
import cmfg
from Parser import Parser
from sys import argv
import numpy as np
import pickle
import math as m
from matplotlib import pyplot as plt

Mycmap = plt.get_cmap('PuOr')

def rebin(M, R, T, rstart=0, tstart=None, cyclic=False):
    """perform a rebinning of a matrix.

    parameters
    ----------
    M : ndarray
        the original matrix
    R : array or list
        the bin groupings in the second index
    T : array or list
        the bin groupings in the first index

    returns
    -------
    Ar : ndarray
        A rebinned matrix
    """
    sx, sy = M.shape
    Ar = M.copy()
    d = 0
    i = rstart
    if tstart is None:
        J = [0]*M.shape[1]
    else:
        J = tstart

    for r, t, j0 in zip(R, T, J):
        j = j0
        while j < sx:
            w = Ar[j:(j+t+d), i:(i+r+d)]
            Ar[j:(j+t+d), i:(i+r+d)] = np.sum(w)/(w.shape[0]*w.shape[1])
            j = j + t                          
        if cyclic:
            print(i)
            p = np.sum(Ar[:j0, i:(i+r+d)] + Ar[-j0:, i:(i+r+d)])/(2*j0*r+d)
            Ar[:j0, i:(i+r+d)] = p
            Ar[-j0:, i:(i+r+d)] = p
        i = i + r
    return Ar


# plots
from pylab import*
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib import pyplot as plt
from matplotlib import ticker

if len(argv) > 1:
    config = Parser(argv[1])
else:
    config = Parser('../set/config_SYSH_020.ini')

f_input = (f"{config.p.dir_output}{config.p.experiment_id}"
           f"/R_{config.p.experiment_id}.pk")
print(f_input)

with open(f_input, 'rb') as f:
    H, K = pickle.load(f)



Ht = np.zeros([config.p.r_n_bins, config.p.theta_n_bins])
for h in H:
    Ht = Ht + h
Kt = np.zeros([config.p.r_n_bins, config.p.theta_n_bins])
for k in K:
    Kt = Kt + k
# profile: las galaxias mas grandes pesan mas
profile = Ht / np.maximum(Kt, 1)

h = profile.transpose()
#h = H[0].transpose()

R = [8, 8, 8, 8, 16, 16]
T = [8, 8, 8, 8, 8, 8]
J = [4]*6

A = rebin(h, R, T, tstart=J, cyclic=True)


fig = plt.figure()
ax = fig.subplots(2,1)
ax[0].imshow(h, origin='lower')
ax[1].imshow(A, origin='lower')
fig.savefig('plot_001.png')
plt.close()





Nrad = config.p.r_n_bins
Nth = config.p.theta_n_bins
rad = np.linspace(config.p.r_start, config.p.r_stop, Nrad+1)
rad = rad.value
rads = np.array((rad[1:]+rad[:-1])/2)

azm = np.linspace(config.p.theta_start, config.p.theta_stop, Nth+1)
azs = (azm[1:]+azm[:-1])/2
azs = azs.value
azs = azs - (azs[1]-azs[0])/2
r, th = np.meshgrid(rad, azm.value)

fig = plt.figure()
ax1 = fig.add_subplot(221, projection='polar')
ax1.pcolormesh(th, r, h, cmap=Mycmap,
                 linestyle='None')

ax2 = fig.add_subplot(222, projection='polar')
ax2.pcolormesh(th, r, A, cmap=Mycmap,
                 linestyle='None')
fig.savefig('plot_002.png')

