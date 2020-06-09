# FCMB

from PixelSky import SkyMap
# from SkyUnits import SkyUnits
import pandas as pd
import healpy as hp
import numpy as np
import astropy.units as u
# import time
from scipy.spatial.transform import Rotation as R
from math import atan2


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
        self.tracer = None

    def load_centers(self):
        """load centers.
        """
        conf = self.config.filenames
        # read Galaxy catalog
        glx_catalog = conf.datadir_glx + conf.filedata_glx
        glx = pd.read_csv(glx_catalog, delim_whitespace=True, header=9)

        phi_healpix = glx['RAdeg']*np.pi/180.
        theta_healpix = (90. - glx['DECdeg'])*np.pi/180.
        glx['vec'] = hp.ang2vec(theta_healpix, phi_healpix).tolist()
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

    def run(self, centers, tracers,
            parallel=False, njobs=1):
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

        skymap = self.tracers

        for center in centers:
            # calcular la matriz de rotacion
            phi = float(center.phi)
            theta = float(center.theta)
            pa = float(center.pa)

            r = R.from_euler('zxz', [phi, theta, pa], degrees=True)
            vector = hp.ang2vec(center.theta, center.phi)[0]

            radius = max(self.breaks_rad).to(u.rad).value

            # querydisc
            listpixs = hp.query_disc(
                    skymap.nside,
                    vector,
                    radius,
                    inclusive=True,
                    fact=4,
                    nest=False)

            dists = []
            thetas = []

            temps = skymap.data[listpixs]

            for ipix in listpixs:

                v = hp.pix2vec(skymap.nside, ipix)
                w = r.apply(v)

                # la galaxia en el nuevo sistema de coordenadas esta en [0,0,1]
                dist = hp.rotator.angdist(w, [0, 0, 1])
                """el angulo wrt el plano de la galaxia (dado por PA) sale
                de la proyeccion en los nuevos ejes X (o sea w[0])
                e Y (o sea w[1])"""
                theta = atan2(w[1], w[0])

                dists.append(dist[0])
                thetas.append(theta)

            dists = dists * u.rad
            thetas = thetas * u.rad

            rr = self.breaks_rad.to(u.rad)
            # aa = self.breaks_ang.to(u.rad)
            bins2d = [rr, self.breaks_ang]

        H = np.histogram2d(dists, thetas, bins=bins2d,
                           weights=temps, density=False)[0]
        K = np.histogram2d(dists, thetas, bins=bins2d, density=False)[0]

        # return([dists, thetas, temps, bins2d])
        return([H, K])
