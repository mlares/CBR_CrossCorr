from Parser import Parser
from PixelSky import SkyMap, PixelTools
from SkyUnits import SkyUnits
import pandas as pd
import healpy as hp
import numpy as np
import astropy.units as u
import time
from scipy.spatial.transform import Rotation as R
from math import atan2


def unwrap_run(correlation, c, **kwarg):
    """Wrap the serial function for parallel run.

    This function just call the serialized version, but allows to run
    it concurrently.
    """
    return correlation.run(c, **kwarg)


class Correlation:
    #{{{
    '''
    class Correlation
    methods for computing angular correlations in the CMB
    methods:
        set_breaks: select bin scheme for the profile
        radialprofile: computes the radial profilee
    ''' 

    # def __init__(self, breaks=[0], skymask=None, nside=256, nran=0, njobs=1):
    def __init__(self, config):
        #{{{
        import numpy as np
        import healpy as hp

        #self.breaks = breaks
        #self.N = len(breaks)-1
        #self.nside = nside
        #npixs = hp.nside2npix(self.nside)
        #self.masked_indices = [i for i in range(npixs) if skymask.data[i]]
        #self.N_ma_IDs = len(self.masked_indices)
        #self.nran = nran
        #self.njobs = njobs
        self.config = config
        #}}}
 

    def load_centers(self):

        conf = self.config.filenames
        # read Galaxy catalog
        glx_catalog = conf.datadir_glx + conf.filedata_glx
        glx = pd.read_csv(glx_catalog, delim_whitespace=True, header=9)

        phi_healpix = glx['RAdeg']*np.pi/180.
        theta_healpix = (90. - glx['DECdeg'])*np.pi/180.
        glx['vec'] = hp.ang2vec(theta_healpix, phi_healpix).tolist()

        self.centers = 0


    def load_tracers(self):

        conf = self.config.filenames
        print(conf.datadir_cmb + conf.filedata_cmb_mapa)
        nside = int(self.config['cmb']['filedata_cmb_nside'])
        mapa = SkyMap(nside)
        mask = SkyMap(nside)

        filedata = conf.datadir_cmb + conf.filedata_cmb_mapa 
        mapa.load(filedata, field=( int(self.config['cmb']['filedata_field_mapa']) ))
        filedata = conf.datadir_cmb + conf.filedata_cmb_mask
        mask.load(filedata, field=( int(self.config['cmb']['filedata_field_mask']) ))

        # averiguar si esto duplica la memoria
        self.map = mapa
        self.mask = mask


#########  LA IDEA ES HACER UNA FUNCION CORR QUE REEMPLACE A LAS 3
#########  FUNCIONES: RADIALPROFILE, ANISOTROPICPROFILE Y CORRELATION

    def corr_II():
        """Run an experiment, parallel version.

        Paralelization is made on the basis of centers

        Parameters
        ----------
        params: the parameters
        njobs: number of jobs
        """
        from joblib import Parallel, delayed

        Pll = Parallel(n_jobs=njobs, verbose=5, prefer="processes")


        params = self.params.values.tolist()
        ids = np.array(range(len(params))) + 1
        ntr = [interactive]*len(params)

        z = zip([self]*len(params), params, ids, ntr)

        d_experiment = delayed(unwrap_run)

        results = Pll(d_experiment(i) for i in z)                 

        # ojo a esto hay que procesarlo para que sea una correlacion
        return results

    def run(self, centers, tracers, 
             parallel=False, njobs=1):

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
            dist = hp.rotator.angdist(w, [0,0,1])
            # el angulo wrt el plano de la galaxia (dado por PA) sale
            # de la proyeccion en los nuevos ejes X (o sea w[0]) e Y (o sea w[1])
            theta = atan2(w[1], w[0])

            dists.append(dist[0])
            thetas.append(theta)

        dists = dists * u.rad
        thetas = thetas * u.rad

        rr = self.breaks_rad.to(u.rad)
        aa = self.breaks_ang.to(u.rad)
        bins2d = [rr, self.breaks_ang]

        H = np.histogram2d(dists, thetas, bins=bins2d, weights=temps, density=False)[0]
        K = np.histogram2d(dists, thetas, bins=bins2d, density=False)[0]
 
        return([dists, thetas, temps, bins2d])
        #return([H, K])


    def run_initial_(self, centers, tracers, 
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
        import numpy as np
        import healpy as hp
        import astropy.units as u
        import time
        from scipy.spatial.transform import Rotation as R
        from math import atan2

        #if autocorr: (montecarlo)
        #    listpixs_C = 
        #    listpixs_T = 
        #    pix2vec
        #    pix2vec
        #    vec_C =
        #    vec_T = 

        #if profile: (querydisc)
        #    listpixs_C = 
        #    listpixs_T = 
        #    pix2vec
        #    pix2vec
        #    vec_C =
        #    vec_T = 
        #    #rotation

        ## pasar pixs a vecs

        #for i in vec_C:
        #    for i in vec_T:
        #        dist =
        #        theta =
        #        dists.append(dist)
        #        thetas.append(theta)
        #               
        #dists = dists * u.rad
        #thetas = thetas * u.rad

        #rr = self.breaks_rad.to(u.rad)
        #aa = self.breaks_ang.to(u.rad)
        #bins2d = [rr, self.breaks_ang]

        #H = np.histogram2d(dists, thetas, bins=bins2d,
        #                   weights=temps, density=False)[0]
        #K = np.histogram2d(dists, thetas, bins=bins2d,
        #                   density=False)[0]
                       


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
            dist = hp.rotator.angdist(w, [0,0,1])
            # el angulo wrt el plano de la galaxia (dado por PA) sale
            # de la proyeccion en los nuevos ejes X (o sea w[0]) e Y (o sea w[1])
            theta = atan2(w[1], w[0])

            dists.append(dist[0])
            thetas.append(theta)

        dists = dists * u.rad
        thetas = thetas * u.rad

        rr = self.breaks_rad.to(u.rad)
        aa = self.breaks_ang.to(u.rad)
        bins2d = [rr, self.breaks_ang]

        H = np.histogram2d(dists, thetas, bins=bins2d, weights=temps, density=False)[0]
        K = np.histogram2d(dists, thetas, bins=bins2d, density=False)[0]
 
        return([dists, thetas, temps, bins2d])
        #return([H, K])




 


