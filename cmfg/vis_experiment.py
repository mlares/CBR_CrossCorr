# load_ext autoreload
# autoreload 2

# experiment, parsing and math
from Parser import Parser
from sys import argv
import numpy as np
import pickle
from Process import rebin2d, rebin1d, profiles, rt_axes


# plots
# from pylab import*
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
# from matplotlib import cm
# from matplotlib import ticker

Mycmap = plt.get_cmap('PuOr')


if len(argv) > 1:
    config = Parser(argv[1])
else:
    config = Parser()

# ---------------------note----------------------------
# run_experiment.py must be run before this script
# -----------------------------------------------------

# ====================== LOAD DATA ===============================

f_input = (f"{config.p.dir_output}{config.p.experiment_id}"
           f"/profile_{config.p.experiment_id}.pk")
print(f_input)

with open(f_input, 'rb') as f:
    H, K = pickle.load(f)

r_breaks, r_means, t_breaks, t_means = rt_axes(config)

res = profiles(H, K, config)
mean_dT_cells, prof_avg, prof_stack, prof_para, prof_perp = res

N_r = config.p.r_n_bins
N_t = config.p.theta_n_bins


# ====================== LOAD CONTROL SAMPLE ====================

control = []
for r in range(config.p.control_n_samples):

    f_input = (f"{config.filenames.dir_output}{config.p.experiment_id}"
               f"/control_{config.p.experiment_id}_{r}.pk")
    print(f_input)

    with open(f_input, 'rb') as f:
        H, K = pickle.load(f)

    r_breaks, r_means, t_breaks, t_means = rt_axes(config)

    res = profiles(H, K, config)
    mean_dT_cells, prof_avg, prof_stack, prof_para, prof_perp = res

    control.append(prof_avg)


# ===================== PLOTS: PROFILE ==========================

# ____________
# perfil total y muestras de control

fig = plt.figure(figsize=(8, 5))
ax = fig.add_subplot()

ax.plot(r_means, prof_avg, linewidth=2, color='red')
ax.plot(r_means, prof_stack, linewidth=2, color='red')
for r in range(config.p.control_n_samples):
    ax.plot(r_means, control[r], linewidth=2, color='grey')

ax.set_xlabel('radial distance [glx size]')
ax.set_ylabel(r'<$\Delta$T> [$\mu$K]')

fout = (f"{config.filenames.dir_plots}{config.p.experiment_id}/"
        f"{config.p.experiment_id}_p000.png")
fig.savefig(fout)
plt.close()

# ____________
# perfil total y por dos regiones: perp y paralell

fig = plt.figure(figsize=(8, 5))
ax = fig.add_subplot()
ax.plot(r_means, prof_avg, linewidth=2, color='grey')
ax.plot(r_means, prof_stack, linewidth=1, color='k')

ax.plot(r_means, prof_para, color='red')
ax.plot(r_means, prof_perp, color='blue')

ax.set_xlabel('radial distance [glx size]')
ax.set_ylabel(r'<$\Delta$T> [$\mu$K]')

fout = (f"{config.filenames.dir_plots}{config.p.experiment_id}/"
        f"{config.p.experiment_id}_p001.png")
fig.savefig(fout)
plt.close()
# ____________
# perfil total y por dos regiones: perp y paralell REBIN

rb = 3
r_rb, rb_avg = rebin1d(r_means, prof_avg, rb)
r_rb, rb_stack = rebin1d(r_means, prof_stack, rb)
r_rb, rb_para = rebin1d(r_means, prof_para, rb)
r_rb, rb_perp = rebin1d(r_means, prof_perp, rb)

fig = plt.figure(figsize=(8, 5))
ax = fig.add_subplot()

ax.plot(r_means, prof_avg, linewidth=2, color='grey')
ax.plot(r_means, prof_stack, linewidth=1, color='k',
        alpha=0.3,
        marker='o', markerfacecolor='white', markeredgecolor='k')
ax.plot(r_means, prof_para, color='red',
        alpha=0.3,
        marker='o', markerfacecolor='white', markeredgecolor='red')
ax.plot(r_means, prof_perp, color='blue',
        alpha=0.3,
        marker='o', markerfacecolor='white', markeredgecolor='blue')

ax.set_xlabel('radial distance [glx size]')
ax.set_ylabel(r'<$\Delta$T> [$\mu$K]')

ax.plot(r_rb, rb_avg, linewidth=2, color='grey')
ax.plot(r_rb, rb_para, color='red')
ax.plot(r_rb, rb_perp, color='blue')

fout = (f"{config.filenames.dir_plots}{config.p.experiment_id}/"
        f"{config.p.experiment_id}_p002.png")
fig.savefig(fout)
plt.close()

# ____________
# matriz REBIN


R = [4, 4, 4, 4, 4]
T = [8, 4, 4, 2, 1]
J = [0, 0, 0, 0, 0]


Aorig = mean_dT_cells.transpose()

Arebin = rebin2d(Aorig, R, T, tstart=J, cyclic=False)
exts = (config.p.r_start.value, config.p.r_stop.value,
        config.p.theta_start.value, config.p.theta_stop.value)

fig = plt.figure(figsize=(8, 5))
ax = fig.subplots(2, 1)
ax[0].imshow(Aorig, extent=exts, origin='lower')
ax[0].set_xlabel('radial distance [glx size]')
ax[0].set_ylabel('angular distance [rad]')

cscale = ax[1].imshow(Arebin, extent=exts, origin='lower')
ax[1].set_xlabel('radial distance [glx size]')
ax[1].set_ylabel('angular distance [rad]')

fig.colorbar(cscale, ax=ax, shrink=0.8, aspect=50,
             label=r'<$\Delta$T> [$\mu$K]')

fout = (f"{config.filenames.dir_plots}{config.p.experiment_id}/"
        f"{config.p.experiment_id}_p003.png")
fig.savefig(fout)
plt.close()

# ____________
# augmented and rebinned polar matrix

fig = plt.figure(figsize=(8, 5))
ax = Axes3D(fig)

Nrep = 30
rads = np.repeat(r_means, Nrep)

azm = np.linspace(config.p.theta_start, config.p.theta_stop, Nrep*(N_t)+1)
azs = (azm[1:]+azm[:-1])/2
azs = azs.value

r, th = np.meshgrid(rads, azs)

A = np.repeat(Arebin, Nrep, axis=0)
z = np.repeat(A, Nrep, axis=1)
zmin = z.min()
zmax = z.max()
zext = max(abs(zmin), abs(zmax))

plt.subplot(projection="polar")
cscale = plt.pcolormesh(th, r, z, cmap=Mycmap,
                        vmin=-zext, vmax=zext,
                        linestyle='None')
fig.colorbar(cscale, ax=ax, shrink=0.8, aspect=50,
             label=r'<$\Delta$T> [$\mu$K]')
#
fout = (f"{config.filenames.dir_plots}{config.p.experiment_id}/"
        f"{config.p.experiment_id}_p007.png")
fig.savefig(fout)


# ____________
# matriz polar

r, th = np.meshgrid(r_breaks, t_breaks)

fig = plt.figure(figsize=(8, 5))
ax = fig.add_subplot(projection='polar')
cscale = ax.pcolormesh(th, r, mean_dT_cells.transpose(), cmap=Mycmap,
                       linestyle='None')
fig.colorbar(cscale, ax=ax, shrink=0.8, aspect=50,
             label=r'<$\Delta$T> [$\mu$K]')
fout = (f"{config.filenames.dir_plots}{config.p.experiment_id}/"
        f"{config.p.experiment_id}_p004.png")
fig.savefig(fout)
plt.close()

# ____________
# augmented polar matrix

fig = plt.figure(figsize=(8, 5))
ax = Axes3D(fig)

Nrep = 100
rads = np.repeat(r_means, Nrep)

azm = np.linspace(config.p.theta_start, config.p.theta_stop, Nrep*(N_t)+1)
azs = (azm[1:]+azm[:-1])/2
azs = azs.value

r, th = np.meshgrid(rads, azs)

A = np.repeat(mean_dT_cells, 100, axis=0)
z = np.repeat(A, 100, axis=1)
z = z.transpose()
zmin = z.min()
zmax = z.max()
zext = max(abs(zmin), abs(zmax))

plt.subplot(projection="polar")
cscale = plt.pcolormesh(th, r, z, cmap=Mycmap,
                        vmin=-zext, vmax=zext,
                        linestyle='None',
                        shading='gouraud')
fig.colorbar(cscale, ax=ax, shrink=0.8, aspect=50,
             label=r'<$\Delta$T> [$\mu$K]')
fout = (f"{config.filenames.dir_plots}{config.p.experiment_id}/"
        f"{config.p.experiment_id}_p005.png")
fig.savefig(fout)


# ____________
# augmented polar matrix

fig = plt.figure(figsize=(8, 5))
ax = Axes3D(fig)

Nrep = 100
rads = np.repeat(r_means, Nrep)

azm = np.linspace(config.p.theta_start, config.p.theta_stop, Nrep*(N_t)+1)
azs = (azm[1:]+azm[:-1])/2
azs = azs.value

r, th = np.meshgrid(rads, azs)

A = np.repeat(mean_dT_cells, 100, axis=0)
z = np.repeat(A, 100, axis=1)
z = z.transpose()
zmin = z.min()
zmax = z.max()
zext = max(abs(zmin), abs(zmax))

plt.subplot(projection="polar")
cscale = plt.pcolormesh(th, r, z, cmap=Mycmap,
                        vmin=-zext, vmax=zext,
                        linestyle='None')
fig.colorbar(cscale, ax=ax, shrink=0.8, aspect=50,
             label=r'<$\Delta$T> [$\mu$K]')

fout = (f"{config.filenames.dir_plots}{config.p.experiment_id}/"
        f"{config.p.experiment_id}_p006.png")
fig.savefig(fout)
plt.close('all')

# # PLOT: Profile
# # ----------------------------------------------------------------------
# fig = plt.figure()
# ax = fig.subplots(1, 1)
# ax.plot(rvals[[0,-1]], [0, 0], color='grey', linewidth=2)
# for i, r in enumerate(R_para):
#     rmicrok = r * 1.e6
#     ax.plot(rvals, rmicrok, color='orange', linewidth=5, alpha=0.3)
# for i, r in enumerate(R_perp):
#     rmicrok = r * 1.e6
#     ax.plot(rvals, rmicrok, color='cadetblue', linewidth=5, alpha=0.5)
# r_para_micro = R_para.sum(axis=0) / len(para) * 1.e6
# ax.plot(rvals, r_para_micro, color='red', label='parallel to glx disk')
# r_perp_micro = R_perp.sum(axis=0) / len(perp) * 1.e6
# ax.plot(rvals, r_perp_micro, color='blue',
#         label='perpendicular to glx disk')
# ax.plot(rvals, profile_microk, color='black', alpha=1, linewidth=2)
# ax.grid()
# ax.legend()
# ax.set_xlabel('angular distance to glx center / glx ang size')
# ax.set_ylabel('<$\Delta$T> [$\mu$K]')
# fout = (f"{config.filenames.dir_plots}{config.p.experiment_id}/"
#         f"{config.p.experiment_id}_align_angular.png")
# fig.savefig(fout)
# plt.close('all')
# # PLOT: each center
# # ----------------------------------------------------------------------
# fig = plt.figure()
# ax = fig.subplots(2, 1)
# # perfiles de la muestra total y de cada centro,
# # integrando todos los angulos:
# for h, k in zip(H, K):
#     hi = h.sum(axis=1)
#     ki = k.sum(axis=1)
#     #if np.any(ki<1): continue
#     p = hi / np.maximum(ki, 1)
#     ax[0].plot(rvals, p, color='cadetblue', alpha=0.1)
# ax[0].plot(rvals, profile, color='orange', alpha=1)
# # perfiles de la muestra total por angulos:
# R = Ht / np.maximum(Kt, 1)
# Rt = R.transpose()
# mn = Rt.mean()
# for r in Rt:
#     ax[1].plot(rvals, r, alpha=0.7)
# ax[1].plot(rvals, profile, color='orange', alpha=1)
# ax[0].grid()
# ax[1].grid()
# fout = (f"{config.filenames.dir_plots}{config.p.experiment_id}/"
#         f"{config.p.experiment_id}_align_centers.png")
# fig.savefig(fout)
# plt.close('all')
# # PLOT: each polar augmented
# # ----------------------------------------------------------------------
# fig = plt.figure()
# ax = Axes3D(fig)
# Nrad = config.p.r_n_bins
# Nth = config.p.theta_n_bins
# rad = np.linspace(config.p.r_start, config.p.r_stop, Nrad+1)
# rad = rad.value
# rads = np.array((rad[1:]+rad[:-1])/2)
# Nrep = 100
# rads = np.repeat(rads, Nrep)
# azm = np.linspace(config.p.theta_start, config.p.theta_stop, Nrep*(Nth)+1)
# azs = (azm[1:]+azm[:-1])/2
# azs = azs.value
# r, th = np.meshgrid(rads, azs)
# Rt_mK = Rt * 1.e6
# A = np.repeat(Rt_mK, 100, axis=0)
# z = np.repeat(A, 100, axis=1)
# zmin = z.min()
# zmax = z.max()
# zext = max(abs(zmin), abs(zmax))
# plt.subplot(projection="polar")
# cscale = plt.pcolormesh(th, r, z, cmap=Mycmap,
#                    vmin=-zext, vmax=zext,
#                    linestyle='None',
#                    shading='gouraud')
# fig.colorbar(cscale, ax=ax, shrink=0.8, aspect=50,
#              label='<$\Delta$T> [$\mu$K]')
# plt.plot(azs, r, color='none') #, ls='none')
# fout = (f"{config.filenames.dir_plots}{config.p.experiment_id}/"
#         f"{config.p.experiment_id}_align_polar.png")
# fig.savefig(fout)
# # PLOT: each polar matrix
# # ----------------------------------------------------------------------
# fig = plt.figure()
# ax1 = Axes3D(fig)
# Nrad = config.p.r_n_bins
# Nth = config.p.theta_n_bins
# rad = np.linspace(config.p.r_start, config.p.r_stop, Nrad+1)
# rad = rad.value
# rads = np.array((rad[1:]+rad[:-1])/2)
# azm = np.linspace(config.p.theta_start, config.p.theta_stop, Nth+1)
# azs = (azm[1:]+azm[:-1])/2
# azs = azs.value
# azs = azs - (azs[1]-azs[0])/2
# r, th = np.meshgrid(rad, azm.value)
# plt.subplot(projection="polar")
# a=plt.pcolormesh(th, r, Rt, cmap=Mycmap,
#                  linestyle='None')
# plt.colorbar(a)
# #plt.plot(azs, r, color='k', ls='none')
#
# fout = (f"{config.filenames.dir_plots}{config.p.experiment_id}/"
#         f"{config.p.experiment_id}_align_polar_original.png")
# fig.savefig(fout)
#
# # PLOT: each polar box
# # ----------------------------------------------------------------------
# fig = plt.figure()
# ax = fig.add_subplot()
#
# Nrad = config.p.r_n_bins
# Nth = config.p.theta_n_bins
# rad = np.linspace(config.p.r_start, config.p.r_stop, Nrad+1)
# rad = rad.value
# rads = np.array((rad[1:]+rad[:-1])/2)
# azm = np.linspace(config.p.theta_start, config.p.theta_stop, Nth+1)
# azs = (azm[1:]+azm[:-1])/2
# azs = azs.value
# r, th = np.meshgrid(rads, azs)
# Rt_mK = Rt*1.e6
# zmin = Rt_mK.min()
# zmax = Rt_mK.max()
# zext = max(abs(zmin), abs(zmax))
# cscale = ax.imshow(Rt_mK, aspect='auto', origin='lower',
#                    extent=(rad[0], rad[-1],
#                            azm[0].value, azm[-1].value),
#                    vmin=-zext, vmax=zext,
#                    cmap=Mycmap)
#
# ax.set_xlabel('angular distance to glx. center / glx. size')
# ax.set_ylabel('angle wrt glx. disk semi-major axis [rad]')
# fig.colorbar(cscale, ax=ax, shrink=0.8, aspect=50,
#              label='<$\Delta$T> [$\mu$K]')
#              #format=ticker.FuncFormatter(fmt))
# plt.tight_layout()
# fout = (f"{config.filenames.dir_plots}{config.p.experiment_id}/"
#         f"{config.p.experiment_id}_align_polar_box.png")
# fig.savefig(fout)
# ## ----------------------------------------------------------------------
# #from pylab import*
# #from mpl_toolkits.mplot3d import Axes3D
# #from matplotlib import cm
# #Nrad = config.p.r_n_bins
# #Nth = config.p.theta_n_bins
# #rad = np.linspace(config.p.r_start, config.p.r_stop, Nrad+1)
# #rads = np.array((rad[1:].value+rad.value[:-1])/2)
# #Nrep = 100
# #rads = np.repeat(rads, Nrep)
# #azm = np.linspace(config.p.theta_start, config.p.theta_stop, Nrep*(Nth)+1)
# #azs = (azm[1:]+azm[:-1])/2
# #azs = azs.value
# #r, th = np.meshgrid(rads, azs)
# #A = np.repeat(Rt, 100, axis=0)
# #z = np.repeat(A, 100, axis=1)
# #fig = plt.figure()
# #ax1 = fig.add_subplot(221, projection='polar')
# #a = ax1.pcolormesh(th, r, z, cmap=plt.get_cmap('Spectral'))
# #plt.thetagrids([theta * 15 for theta in range(360//15)])
# #plt.grid()
# #ax2 = fig.add_subplot(222, projection='polar')
# #a = ax2.pcolormesh(th, r, z, cmap=plt.get_cmap('Spectral'))
# #plt.thetagrids([theta * 15 for theta in range(360//15)])
# #plt.grid()
# #fout = f"{config.p.experiment_id}_align_polar2.png"
# #fig.savefig(fout)
#
# """
# TO DO:
#
# - barras de error
# - sumar bines para los chicos (adaptative binning)  DONE
# - smooth matrix
# - lista de settings
#
# """
