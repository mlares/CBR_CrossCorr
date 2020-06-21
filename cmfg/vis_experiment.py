# load_ext autoreload
# autoreload 2

import cmfg
from Parser import Parser
from sys import argv
import numpy as np
import pickle
from matplotlib import pyplot as plt

if len(argv) > 1:
    config = Parser(argv[1])
else:
    config = Parser()

#X = cmfg.Correlation(config)
#X.load_centers()
#X.load_tracers()
#X.select_subsample_centers()

fout = f"R_{config.p.experiment_id}.pk"
print(fout)

with open(fout, 'rb') as f:
    res = pickle.load(f)

H, K = res
rvals = np.linspace(config.p.r_start, config.p.r_stop, config.p.r_n_bins)


Ht = np.zeros([config.p.r_n_bins, config.p.theta_n_bins])
for h in H:
    Ht = Ht + h
Kt = np.zeros([config.p.r_n_bins, config.p.theta_n_bins])
for k in K:
    Kt = Kt + k
profile = Ht.sum(axis=1) / np.maximum(Kt.sum(axis=1), 1)


 
fig = plt.figure()
ax = fig.subplots(1, 1)

# perfiles de la muestra total por angulos:
R = Ht / np.maximum(Kt, 1)
Rt = R.transpose()
m = Rt.mean()
profile_microk = profile * 1.e6

ax.plot(rvals[[0,-1]], [0, 0], color='grey', linewidth=2)
#for r in Rt[[0, 1, 6, 7, 8, 9, 14, 15]]:
#for r in Rt[[2, 3, 4, 5, 10, 11, 12, 13]]:
for i, r in enumerate(Rt):
    rmicrok = r * 1.e6
    if i in [0, 1, 6, 7, 8, 9, 14, 15]:
        col = 'orange'
    if i in [2, 3, 4, 5, 10, 11, 12, 13]:
        col = 'cadetblue'
    ax.plot(rvals, rmicrok, color=col, alpha=0.7)
ax.plot(rvals, profile_microk, color='black', alpha=1, linewidth=3)

ax.grid()
ax.legend()
#ax.set_ylim(-3.e-5, 3.e-5)
fig.savefig('plot0.png')
plt.close('all')
 



fig = plt.figure()
ax = fig.subplots(2, 1)

# perfiles de la muestra total y de cada centro,
# integrando todos los angulos:
for h, k in zip(H, K):
    hi = h.sum(axis=1)
    ki = k.sum(axis=1)
    #if np.any(ki<1): continue
    p = hi / np.maximum(ki, 1)
    ax[0].plot(rvals, p, color='cadetblue', alpha=0.1)
ax[0].plot(rvals, profile, color='orange', alpha=1)


# perfiles de la muestra total por angulos:

R = Ht / np.maximum(Kt, 1)
Rt = R.transpose()

m = Rt.mean()

#ax[1].plot(rvals[[0,-1]], [m,m], color='grey', linewidth=2)
for r in Rt:
    ax[1].plot(rvals, r, alpha=0.7)
ax[1].plot(rvals, profile, color='orange', alpha=1)

ax[0].grid()
ax[1].grid()
#ax.legend()
fig.savefig('plot.png')
plt.close('all')

