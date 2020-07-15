# experiment, parsing and math
import cmfg
from Parser import Parser
from sys import argv
import numpy as np
import pickle
import math as m
from matplotlib import pyplot as plt

# plots
from pylab import*
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib import pyplot as plt
from matplotlib import ticker

if len(argv) > 1:
    config = Parser(argv[1])
else:
    config = Parser()

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

# promedio de los perfiles (todas las galaxias pesan igual)
prof_avg = profile.sum(axis=1)

# perfil promedio del stacking
prof_tot = Ht.sum(axis=1) / np.maximum(Kt.sum(axis=1), 1)

h = profile.transpose()
#h = H[0].transpose()

R = [8, 8, 8, 8, 16, 16]
T = [8, 8, 4, 2, 2, 1]
J = [0, 0, 0, 0, 0, 0]

A = rebin(h, R, T, tstart=J, cyclic=False)

Nrad = config.p.r_n_bins
Nth = config.p.theta_n_bins
rad = np.linspace(config.p.r_start, config.p.r_stop, Nrad+1)
rad = rad.value
rads = np.array((rad[1:]+rad[:-1])/2)
h = profile.transpose()

azm = np.linspace(config.p.theta_start, config.p.theta_stop, Nth+1)
azs = (azm[1:]+azm[:-1])/2
azs = azs.value
azs = azs - (azs[1]-azs[0])/2
r, th = np.meshgrid(rad, azm.value)

fig = plt.figure(figsize=(8,5))
ax1 = fig.add_subplot(221, projection='polar')
ax1.pcolormesh(th, r, h, cmap=Mycmap,
                 linestyle='None')

ax2 = fig.add_subplot(222, projection='polar')
ax2.pcolormesh(th, r, A, cmap=Mycmap,
                 linestyle='None')
fig.savefig('plot_002.png')
plt.close()

fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(projection='polar')
ax.pcolormesh(th, r, h, cmap=Mycmap, linestyle='None')
fig.savefig('plot_003.png')
plt.close()

fig = plt.figure(figsize=(6,6))
ax = fig.add_subplot(projection='polar')
ax.pcolormesh(th, r, A, cmap=Mycmap, linestyle='None')
fig.savefig('plot_004.png')
plt.close()