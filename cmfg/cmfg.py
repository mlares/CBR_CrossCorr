# FCMB

from PixelSky import SkyMap
# from SkyUnits import SkyUnits
import pandas as pd
import healpy as hp
import numpy as np
from astropy import units as u
# import time
from scipy.spatial.transform import Rotation as R
from math import atan2, pi
from tqdm import tqdm


def unwrap_run(correlation, c, **kwarg):
    """Wrap the serial function for parallel run.

    This function just call the serialized version, but allows to run
    it concurrently.
    """
    return correlation.run(c, **kwarg)


class Correlation:
    '''
    class Correlation
    methods for computing angular correlations in the CMB.

    Methods
    -------
    load_centers:
    load_tracers:
    corr_II:
    run: computes the radial profile.
    '''

    def __init__(self, config):
        self.config = config
        self.centers = None
        self.map = None
        self.mask = None

    def load_centers(self):
        """load centers.
        """
        conf = self.config.filenames
        # read Galaxy catalog
        glx_catalog = conf.datadir_glx + conf.filedata_glx
        glx = pd.read_csv(glx_catalog, delim_whitespace=True, header=9)

        phi_healpix = glx['RAdeg']*np.pi/180.
        theta_healpix = (90. - glx['DECdeg'])*np.pi/180.
        glx['phi'] = phi_healpix
        glx['theta'] = theta_healpix
        glx['vec'] = hp.ang2vec(theta_healpix, phi_healpix).tolist()
        glx['pa'] = np.zeros(glx.shape[0])

        self.centers = glx

    def load_tracers(self):
        """load tracers.
        """
        conf = self.config.filenames
        print(conf.datadir_cmb + conf.filedata_cmb_mapa)
        nside = int(self.config['cmb']['filedata_cmb_nside'])
        mapa = SkyMap(nside)
        mask = SkyMap(nside)

        filedata = conf.datadir_cmb + conf.filedata_cmb_mapa
        column = int(self.config['cmb']['filedata_field_mapa'])
        mapa.load(filedata, field=(column))
        filedata = conf.datadir_cmb + conf.filedata_cmb_mask
        column = int(self.config['cmb']['filedata_field_mask'])
        mask.load(filedata, field=(column))

        # averiguar si esto duplica la memoria
        self.map = mapa
        self.mask = mask

    def corr_II(self):
        """Run an experiment, parallel version.

        Paralelization is made on the basis of centers

        Parameters
        ----------
        params: the parameters
        njobs: number of jobs
        """
        from joblib import Parallel, delayed

        njobs = self.conf.p.njobs
        Pll = Parallel(n_jobs=njobs, verbose=5, prefer="processes")

        params = self.params.values.tolist()
        ids = np.array(range(len(params))) + 1

        interactive = self.conf.p.interactive
        ntr = [interactive]*len(params)

        z = zip([self]*len(params), params, ids, ntr)

        d_experiment = delayed(unwrap_run)

        results = Pll(d_experiment(i) for i in z)

        # ojo a esto hay que procesarlo para que sea una correlacion
        return results

    def run(self, parallel=False, njobs=1):
        """Compute correlations in CMB data.

        When centers are fixed, it returns the stacked radial profile.


        Parameters
        ----------
        centers : list or array
            List of centers
        tracers : list or array
            List of tracers

        Returns
        -------
        H : array like
            Cross product of temperatures
        K : array like
            Counts of pairs contributing to each bin
        profile : array like
            Array containing the mean temperature per radial bin.
            Angular bins are also returned if required from
            configuration.  The scale of the radial coordinate depends
            on the configuration. All configuration parameters are
            stored in self.config
        """
        if self.config.p.verbose:
            print('starting computations...')
        skymap = self.map
        centers = self.centers

        Nb_r = self.config.p.r_n_bins
        Nb_a = self.config.p.theta_n_bins

        # breaks for angular distance, in radians
        rmin = self.config.p.r_start
        rmax = self.config.p.r_stop
        rmin = rmin.to(u.rad).value
        rmax = rmax.to(u.rad).value
        radii = np.linspace(rmin, rmax, Nb_r+1)

        # breaks for angle, in radians
        amin = self.config.p.theta_start
        amax = self.config.p.theta_stop
        amin = amin.to(u.rad).value
        amax = amax.to(u.rad).value
        angles = np.linspace(amin, amax, Nb_a+1)

        # initialize 2d histogram
        bins2d = [radii, angles]
        Ht = np.zeros([Nb_r, Nb_a])
        Kt = np.zeros([Nb_r, Nb_a])

        if self.config.p.showp:
            bf1 = "{desc}: {percentage:.4f}% | "
            bf2 = "{n_fmt}/{total_fmt} ({elapsed}/{remaining})"
            bf = ''.join([bf1, bf2])
            total = centers.shape[0]
            iterator = tqdm(centers.itertuples(), total=total, bar_format=bf)
        else:
            iterator = centers.itertuples()

        for center in iterator:

            # compute rotation matrix
            phi = float(center.phi)
            theta = float(center.theta)
            pa = float(center.pa)

            vector = hp.ang2vec(center.theta, center.phi)
            rotate_pa = R.from_euler('zyz', [-phi, -theta, pa])

            # querydisc
            listpixs = hp.query_disc(skymap.nside,
                                     vector,
                                     rmax,
                                     inclusive=True,
                                     fact=4,
                                     nest=False)

            dists = []
            thetas = []
            temps = skymap.data[listpixs]

            for ipix in listpixs:

                v = hp.pix2vec(skymap.nside, ipix)
                w = rotate_pa.apply(v)

                # each center is in position [0,0,1] wrt the new system
                dist = hp.rotator.angdist(w, [0, 0, 1])
                """el angulo wrt el plano de la galaxia (dado por PA) sale
                de la proyeccion en los nuevos ejes X (o sea w[0])
                e Y (o sea w[1])"""
                theta = atan2(w[1], w[0])
                if theta < 0:
                    theta = theta + 2*pi

                dists.append(dist[0])
                thetas.append(theta)

            H = np.histogram2d(dists, thetas, bins=bins2d,
                               weights=temps, density=False)
            H = H[0]
            K = np.histogram2d(dists, thetas, bins=bins2d,
                               density=False)
            K = K[0]
            print(H)

        Ht = Ht + H
        Kt = Kt + K

        return([Ht, Kt])

    def testrot(self):
        centers = self.centers
        skymap = self.map
        if self.config.p.showp:
            bf1 = "{desc}: {percentage:.4f}% | "
            bf2 = "{n_fmt}/{total_fmt} ({elapsed}/{remaining})"
            bf = ''.join([bf1, bf2])
            total = centers.shape[0]
            iterator = tqdm(centers.itertuples(), total=total, bar_format=bf)
        else:
            iterator = centers.itertuples()

        for center in iterator:

            # compute rotation matrix
            phi = float(center.phi)
            theta = float(center.theta)
            pa = float(center.pa)

            vector = hp.ang2vec(center.theta, center.phi)
            rotate_pa = R.from_euler('zyz', [-phi, -theta, pa])
            listpixs = hp.query_disc(
                        skymap.nside,
                        vector,
                        0.15,
                        inclusive=True,
                        fact=4,
                        nest=False)

            for ipix in listpixs:
                v = hp.pix2vec(skymap.nside, ipix)
                w = rotate_pa.apply(v)
                print(v, w)
                # w must be ~ [0,0,1]


def test_rotation(N):
    from math import pi, acos
    from random import random
    alphas = [random()*360 for _ in range(N)]
    deltas = [90 - acos(random()*2 - 1)*180/pi for _ in range(N)]
    for alpha, delta in zip(alphas, deltas):
        colatitude = 90-delta
        theta_healpix = (colatitude)*pi/180
        phi_healpix = alpha*pi/180
        v = hp.ang2vec(theta_healpix, phi_healpix)
        rotate_pa = R.from_euler('zyz', [-alpha, -colatitude, 0], degrees=True)
        w = rotate_pa.apply(v)
        print(f"alpha={alpha}, delta={delta}\n v={v}\nw={w}\n")
