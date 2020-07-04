import numpy as np
import math as m

def rt_axes(config):
     
    N_r = config.p.r_n_bins
    r_breaks = np.linspace(config.p.r_start, config.p.r_stop, N_r+1)
    r_breaks = r_breaks.value
    r_means = np.array((r_breaks[1:]+r_breaks[:-1])/2)

    N_t = config.p.theta_n_bins
    t_breaks = np.linspace(config.p.theta_start, config.p.theta_stop, N_t+1)
    t_breaks = t_breaks.value
    t_means = (t_breaks[1:]+t_breaks[:-1])/2

    return r_breaks, r_means, t_breaks, t_means
                                                         


def profiles(H, K, config):

    r_breaks, r_means, t_breaks, t_means = rt_axes(config)
                                                       
    N_r = config.p.r_n_bins
    N_t = config.p.theta_n_bins

    Ht = np.zeros([N_r, N_t])
    for h in H:
        Ht = Ht + h * 1.e6

    Kt = np.zeros([N_r, N_t])
    for k in K:
        Kt = Kt + k

    # promedio de los perfiles (todas las galaxias pesan igual)
    mean_dT_cells = Ht / np.maximum(Kt, 1)
    prof_avg = mean_dT_cells.mean(axis=1)
    assert(len(prof_avg)==N_r)

    # perfil promedio del stacking
    prof_stack = Ht.sum(axis=1) / np.maximum(Kt.sum(axis=1), 1)
    assert(len(prof_stack)==N_r)

    para = []
    perp = []
    for i, tt in enumerate(t_means):
        f = abs(m.cos(tt)) > m.cos(m.pi/4)
        if f:
            para.append(i)
        else:
            perp.append(i)

    mean_dT_cells_trans = mean_dT_cells.transpose()
    prof_para = mean_dT_cells_trans[para].mean(axis=0)
    prof_perp = mean_dT_cells_trans[perp].mean(axis=0)         
    prof_all = mean_dT_cells_trans.mean(axis=0)
     
    assert(len(prof_para)==N_r)
    assert(len(prof_perp)==N_r)
    assert(len(prof_all)==N_r)

    return mean_dT_cells, prof_avg, prof_stack, prof_para, prof_perp


# ===================== PLOT UTILS ==========================

def rebin2d(M, R, T, rstart=0, tstart=None, cyclic=False):
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


def rebin1d(Mx, My, R, rstart=0, tstart=None, cyclic=False):
    """perform a rebinning of a list.

    parameters
    ----------
    M : list or ndarray
        the original array
    R : int
        the bin width

    returns
    -------
    Ar : ndarray
        A rebinned array
    """
    k = len(Mx) // R
    xrb = []
    yrb = []
    for i in range(k):
        ii = R*i
        xrb.append(Mx[ii:(ii+R)].mean())
        yrb.append(My[ii:(ii+R)].mean())
    return(xrb, yrb)

def fmt(x, pos):
    """Set format for numbers in scale
    """
    a, b = '{:.0e}'.format(x).split('e')
    b = int(b)
    if a == 0:
        return r'0'
    else:
        return r'${}\times10^{{{}}}$'.format(a, b)
 
