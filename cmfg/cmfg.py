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


def unwrap_run(arg, **kwarg):
    """Wrap the serial function for parallel run.

    This function just call the serialized version, but allows to run
    it concurrently.
    """
    return Correlation.run_batch(*arg, **kwarg)


class Correlation:
    '''
    Correlation (class): compute angular correlations in CMB maps.

    Methods
    -------
    load_centers : load centers from file using config.
    load_tracers : load tracers from file using config.
    initialize_counters : initialize counters
    run : computes the radial profile for a sample.
    run_single : computes the radial profile for a single center. 
    run_batch : serial computation of the radial profile for a sample. 
    run_batch_II : computes in parallel the radial profile for a sample.
    '''

    def __init__(self, config):
        self.config = config
        self.centers = None
        self.map = None
        self.mask = None

    def load_centers(self):
        """load centers from galaxy catalogue.
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

        self.centers = glx

    def load_tracers(self):
        """load tracers from Healpix map of the CMBR.
        """
        conf = self.config.filenames
        if self.config.p.verbose:
            print(conf.datadir_cmb + conf.filedata_cmb_mapa)

        nside = int(self.config['cmb']['filedata_cmb_nside'])
        mapa = SkyMap(nside)
        mask = SkyMap(nside)

        filedata = conf.datadir_cmb + conf.filedata_cmb_mapa
        column = int(self.config['cmb']['filedata_field_mapa'])
        mapa.load(filedata, field=(column), verbose=self.config.p.verbose)
        filedata = conf.datadir_cmb + conf.filedata_cmb_mask
        column = int(self.config['cmb']['filedata_field_mask'])
        mask.load(filedata, field=(column), verbose=self.config.p.verbose)
        # averiguar si esto duplica la memoria
        self.map = mapa
        self.mask = mask

        npixs = hp.nside2npix(nside)
        self.masked_indices = [i for i in range(npixs) if self.mask.data[i]]

    def initialize_counters(self):
        """Initialize counters for temperature map aroung centers.

        This is kept separate for organization of the code.
        """
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

        return bins2d, rmax, Ht, Kt

    def select_subsample_centers(self):
        """Make a selection of centers.

        The selection is made by filtering on galxy type and distance
        to the galaxy (redshift). This makes use of the configuration
        paramenters from the .ini file.
        """
        zmin = self.config.p.redshift_min
        zmax = self.config.p.redshift_max

        Sa_lbl = ['1','2']
        Sb_lbl = ['3','4']
        Sc_lbl = ['5','6','7','8']
        Sd_lbl = ['7','8']
        Sno_lbl = ['10', '11', '12', '15', '16', '19', '20']

        Stypes = []
        gtypes = self.config.p.galaxy_types
        for s in gtypes:
            if s=='a':
                for k in Sa_lbl:
                    Stypes.append(k)
            if s=='b':
                for k in Sb_lbl:
                    Stypes.append(k)
            if s=='c':
                for k in Sc_lbl:
                    Stypes.append(k)
            if s=='d':
                for k in Sd_lbl:
                    Stypes.append(k)

        filt = []
        for t in self.centers['type']:
            f = t[0] in Stypes and not (t[:2] in Sno_lbl)
            filt.append(f)

        self.centers = self.centers[filt]

        return None

    def run_single(self, center):
        """Compute the temperature in CMB data around a center.

        Parameters
        ----------
        center : list or array
            List of centers

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
        skymap = self.map
        bins2d, rmax, Ht, Kt = self.initialize_counters()

        # compute rotation matrix
        phi = float(center[1].phi)
        theta = float(center[1].theta)
        pa = float(center[1].pa)

        vector = hp.ang2vec(center[1].theta, center[1].phi)
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
        K = np.histogram2d(dists, thetas, bins=bins2d,
                           density=False)

        return H[0], K[0]

    def run(self, parallel=None, njobs=1):
        """Compute (stacked) temperature map in CMB data around centers.

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
        if isinstance(parallel, bool):
            run_parallel = parallel
        else:
            run_parallel = self.config.p.run_parallel

        if self.config.p.verbose:
            print('starting computations...')

        centers = self.centers
        bins2d, rmax, Ht, Kt = self.initialize_counters()

        if run_parallel:
            Ht, Kt = self.run_batch_II()
        else:
            centers_ids = range(len(centers))
            Ht, Kt = self.run_batch(centers, centers_ids)

        R = Ht / np.maximum(Kt, 1)

        return R

    def run_batch(self, centers, index):
        """Compute (stacked) temperature map in CMB data around centers.

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

        Notes
        -----
        This function admits a pandas dataframe or a row, with the
        type as traversed with df.iterrows()
        When called from unwrap_run, it receives a "row"
        """
        bins2d, rmax, Ht, Kt = self.initialize_counters()

        if isinstance(centers, tuple):
            # a single center
            Ht, Kt = self.run_single(centers)
        else:
            # a dataframe
            if self.config.p.showp:
                bf1 = "{desc}: {percentage:.4f}% | "
                bf2 = "{n_fmt}/{total_fmt} ({elapsed}/{remaining})"
                bf = ''.join([bf1, bf2])
                total = centers.shape[0]
                iterator = tqdm(centers.iterrows(),
                                total=total, bar_format=bf)
            else:
                total = centers.shape[0]
                iterator = centers.iterrows()

            for center in iterator:
                H, K = self.run_single(center)
                Ht = Ht + H
                Kt = Kt + K

        return Ht, Kt

    def run_batch_II(self):
        """Compute (stacked, parallel) temperature map around centers.

        Paralelization is made on the basis of centers

        Parameters
        ----------
        njobs: number of jobs
        """
        from joblib import Parallel, delayed

        njobs = self.config.p.n_jobs

        if self.config.p.verbose:
            vlevel = 5
        else:
            vlevel = 0
        Pll = Parallel(n_jobs=njobs, verbose=vlevel, prefer="processes")
        #centers = self.centers.values.tolist()
        centers = self.centers
        Ncenters = centers.shape[0]
        ids = np.array(range(Ncenters)) + 1

        cntrs = []
        for c in centers.iterrows():
            cntrs.append(c)

        z = zip([self]*Ncenters, cntrs, ids)

        d_experiment = delayed(unwrap_run)

        results = Pll(d_experiment(i) for i in z)

        Ht, Kt = np.array(results).sum(axis=0)
        return Ht, Kt

    def testrot(self):
        """Test rotation.

        Loads data and extract a small disc centered on each centers,
        in order to erify if the angles between selected pixels and
        their respective centers are small
        """
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
                # w must be ~ [0,0,1]
                print(v, w)


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
        return None
