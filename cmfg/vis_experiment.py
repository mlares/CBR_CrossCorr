# load_ext autoreload
# autoreload 2

# experiment, parsing and math
import cmfg
from Parser import Parser
from sys import argv
import numpy as np
import pickle
import math as m

def rebin(A, R, T):

    sx, sy = A.shape

    i = 0
    for r, t in zip(R, T):

        j = 0
        while j < sy:
            A[i:(i+r+1), j:(j+t+1)] = np.sum(A[i:(i+r+1), j:(j+t+1)])
            j = j + t

        i = i + r
        Ar = A
        return Ar

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
    config = Parser()

# ---------------------note----------------------------
# run_experiment.py must be run before this script
# -----------------------------------------------------

f_input = (f"{config.p.dir_output}{config.p.experiment_id}"
           f"/R_{config.p.experiment_id}.pk")
print(f_input)

with open(f_input, 'rb') as f:
    res = pickle.load(f)

H, K = res
rvals = np.linspace(config.p.r_start, config.p.r_stop, config.p.r_n_bins)

# AGREGAR LAS BANDAS DE ERROR !!!

# Profile: computations
# ----------------------------------------------------------------------
Ht = np.zeros([config.p.r_n_bins, config.p.theta_n_bins])
for h in H:
    Ht = Ht + h
Kt = np.zeros([config.p.r_n_bins, config.p.theta_n_bins])
for k in K:
    Kt = Kt + k

# profile: las galaxias mas grandes pesan mas
profile = Ht.sum(axis=1) / np.maximum(Kt.sum(axis=1), 1)
 
# profile: las galaxias mas grandes pesan mas

P = H / K
profile = P.sum(axis=1)


t = np.linspace(config.p.theta_start, config.p.theta_stop, 
                config.p.theta_n_bins + 1)
tmean = (t[1:] + t[:-1])/2
para = []
perp = []
for i, tt in enumerate(tmean):
    f = abs(m.cos(tt.value)) > m.cos(m.pi/4)
    if f:
        para.append(i)
    else:
        perp.append(i)

R = Ht / np.maximum(Kt, 1)
Rt = R.transpose()
profile_microk = profile * 1.e6
R_para = Rt[para]
R_perp = Rt[perp]         

# Plot: settings
# ----------------------------------------------------------------------
Mycmap = plt.get_cmap('PuOr')

def fmt(x, pos):
    a, b = '{:.0e}'.format(x).split('e')
    b = int(b)
    if a == 0:
        return r'0'
    else:
        return r'${}\times10^{{{}}}$'.format(a, b)


# PLOT: Profile
# ----------------------------------------------------------------------
fig = plt.figure()
ax = fig.subplots(1, 1)

ax.plot(rvals[[0,-1]], [0, 0], color='grey', linewidth=2)

for i, r in enumerate(R_para):
    rmicrok = r * 1.e6
    ax.plot(rvals, rmicrok, color='orange', linewidth=5, alpha=0.3)

for i, r in enumerate(R_perp):
    rmicrok = r * 1.e6
    ax.plot(rvals, rmicrok, color='cadetblue', linewidth=5, alpha=0.5)


r_para_micro = R_para.sum(axis=0) / len(para) * 1.e6
ax.plot(rvals, r_para_micro, color='red', label='parallel to glx disk')

r_perp_micro = R_perp.sum(axis=0) / len(perp) * 1.e6
ax.plot(rvals, r_perp_micro, color='blue', 
        label='perpendicular to glx disk')

ax.plot(rvals, profile_microk, color='black', alpha=1, linewidth=2)

ax.grid()
ax.legend()
ax.set_xlabel('angular distance to glx center / glx ang size')
ax.set_ylabel('<$\Delta$T> [$\mu$K]')

fout = (f"{config.filenames.dir_plots}{config.p.experiment_id}/"
        f"{config.p.experiment_id}_align_angular.png")
fig.savefig(fout)
plt.close('all')
 

# PLOT: each center
# ----------------------------------------------------------------------
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
mn = Rt.mean()

for r in Rt:
    ax[1].plot(rvals, r, alpha=0.7)
ax[1].plot(rvals, profile, color='orange', alpha=1)

ax[0].grid()
ax[1].grid()
fout = (f"{config.filenames.dir_plots}{config.p.experiment_id}/"
        f"{config.p.experiment_id}_align_centers.png")
fig.savefig(fout)
plt.close('all')

# PLOT: each polar augmented
# ----------------------------------------------------------------------
fig = plt.figure()
ax = Axes3D(fig)

Nrad = config.p.r_n_bins
Nth = config.p.theta_n_bins
rad = np.linspace(config.p.r_start, config.p.r_stop, Nrad+1)
rad = rad.value
rads = np.array((rad[1:]+rad[:-1])/2)
Nrep = 100
rads = np.repeat(rads, Nrep)

azm = np.linspace(config.p.theta_start, config.p.theta_stop, Nrep*(Nth)+1)
azs = (azm[1:]+azm[:-1])/2
azs = azs.value

r, th = np.meshgrid(rads, azs)

Rt_mK = Rt * 1.e6
A = np.repeat(Rt_mK, 100, axis=0)
z = np.repeat(A, 100, axis=1)
zmin = z.min()
zmax = z.max()
zext = max(abs(zmin), abs(zmax))

plt.subplot(projection="polar")
cscale = plt.pcolormesh(th, r, z, cmap=Mycmap, 
                   vmin=-zext, vmax=zext,
                   linestyle='None',
                   shading='gouraud')
fig.colorbar(cscale, ax=ax, shrink=0.8, aspect=50,
             label='<$\Delta$T> [$\mu$K]')

plt.plot(azs, r, color='none') #, ls='none') 
fout = (f"{config.filenames.dir_plots}{config.p.experiment_id}/"
        f"{config.p.experiment_id}_align_polar.png")
fig.savefig(fout)

# PLOT: each polar matrix
# ----------------------------------------------------------------------
fig = plt.figure()
ax1 = Axes3D(fig)

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

plt.subplot(projection="polar")
a=plt.pcolormesh(th, r, Rt, cmap=Mycmap,
                 linestyle='None')
plt.colorbar(a)
#plt.plot(azs, r, color='k', ls='none') 

fout = (f"{config.filenames.dir_plots}{config.p.experiment_id}/"
        f"{config.p.experiment_id}_align_polar_original.png")
fig.savefig(fout)


# PLOT: each polar box
# ----------------------------------------------------------------------
fig = plt.figure()
ax = fig.add_subplot()

Nrad = config.p.r_n_bins
Nth = config.p.theta_n_bins
rad = np.linspace(config.p.r_start, config.p.r_stop, Nrad+1)
rad = rad.value
rads = np.array((rad[1:]+rad[:-1])/2)

azm = np.linspace(config.p.theta_start, config.p.theta_stop, Nth+1)
azs = (azm[1:]+azm[:-1])/2
azs = azs.value

r, th = np.meshgrid(rads, azs)
 
Rt_mK = Rt*1.e6
zmin = Rt_mK.min()
zmax = Rt_mK.max()
zext = max(abs(zmin), abs(zmax))
                    
cscale = ax.imshow(Rt_mK, aspect='auto', origin='lower',
                   extent=(rad[0], rad[-1],
                           azm[0].value, azm[-1].value),
                   vmin=-zext, vmax=zext,
                   cmap=Mycmap)

ax.set_xlabel('angular distance to glx. center / glx. size')
ax.set_ylabel('angle wrt glx. disk semi-major axis [rad]')

fig.colorbar(cscale, ax=ax, shrink=0.8, aspect=50,
             label='<$\Delta$T> [$\mu$K]')
             #format=ticker.FuncFormatter(fmt))
plt.tight_layout()

fout = (f"{config.filenames.dir_plots}{config.p.experiment_id}/"
        f"{config.p.experiment_id}_align_polar_box.png")
fig.savefig(fout)




## ----------------------------------------------------------------------
#from pylab import*
#from mpl_toolkits.mplot3d import Axes3D
#from matplotlib import cm
#Nrad = config.p.r_n_bins
#Nth = config.p.theta_n_bins
#rad = np.linspace(config.p.r_start, config.p.r_stop, Nrad+1)
#rads = np.array((rad[1:].value+rad.value[:-1])/2)
#Nrep = 100
#rads = np.repeat(rads, Nrep)
#azm = np.linspace(config.p.theta_start, config.p.theta_stop, Nrep*(Nth)+1)
#azs = (azm[1:]+azm[:-1])/2
#azs = azs.value
#r, th = np.meshgrid(rads, azs)
#A = np.repeat(Rt, 100, axis=0)
#z = np.repeat(A, 100, axis=1)
#fig = plt.figure()
#ax1 = fig.add_subplot(221, projection='polar')
#a = ax1.pcolormesh(th, r, z, cmap=plt.get_cmap('Spectral'))
#plt.thetagrids([theta * 15 for theta in range(360//15)])
#plt.grid()
#ax2 = fig.add_subplot(222, projection='polar')
#a = ax2.pcolormesh(th, r, z, cmap=plt.get_cmap('Spectral'))
#plt.thetagrids([theta * 15 for theta in range(360//15)])
#plt.grid()
#fout = f"{config.p.experiment_id}_align_polar2.png"
#fig.savefig(fout)



"""
TO DO:

- barras de error
- sumar bines para los chicos (adaptative binning)
- smooth matrix
- lista de settings

"""
