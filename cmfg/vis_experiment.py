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
ax = fig.subplots(2, 1)

# perfiles de la muestra total y de cada centro,
# integrando todos los angulos:
for h, k in zip(H, K):
    hi = h.sum(axis=1)
    ki = k.sum(axis=1)
    if np.any(ki<1): continue
    p = hi / np.maximum(ki, 1)
    ax[0].plot(rvals, p, color='cadetblue', alpha=0.7)
ax[0].plot(rvals, profile, color='orange', alpha=1)


# perfiles de la muestra total por angulos:

R = Ht / np.maximum(Kt, 1)
Rt = R.transpose()

m = Rt.mean()

ax[1].plot(rvals[[0,-1]], [m,m], color='grey', linewidth=2)
for r in Rt:
    ax[1].plot(rvals, r, color='cadetblue', alpha=0.7)

ax[0].grid()
ax[1].grid()
#ax.legend()
fig.savefig('plot.png')
plt.close('all')

