from fcmb import Parser
from fcmb import PixelSky
from fcmb import SkyUnits


from joblib import Parallel, delayed

def unwrap_profile_self(arg, **kwarg):
    return RadialProfile.radialprofile(*arg, **kwarg)

def unwrap_anisotropicprofile_self(arg, **kwarg):
    return AnisotropicProfile.anisotropic_profile(c, **kwarg)

def unwrap_correlation_self(correlation, c, **kwarg):
    return correlation.correlation(c, **kwarg)


class RadialProfile:
    # {{{
    '''
    class RadialProfile
    methods for computing angular correlations in the CMB
    methods:
        set_breaks: select bin scheme for the profile
        radialprofile: computes the radial profile
        radialprofile_II: computes the radial profile in parallel
    ''' 

    def __init__(self, breaks=[0], Nran=0):
        #{{{
        """init(self, breaks, Nran) : sets the partition (binning
        scheme) for the computing of the radial profile.

        Tasks:
        1. this works as a combination of np.linspace and units.

        Args:
            unit:
            selected unit for the distance to the center

        Raises:
            errors?

        Returns:
        """               
        import numpy as np

        self.breaks = breaks
        N = len(breaks)-1
        self.N = N
        self.max_centers = 0
        self.signal = np.zeros(N)
        self.sigma = np.zeros(N)
        self.controlsample_mean = np.zeros(N)
        self.controlsample_sigma = np.zeros(N)
        #}}} 
 
    def set_breaks(self, unit, *args, **kwargs):
        #{{{
        """set_breaks(self, unit) : sets the breaks for the binned 
        profile.

        Tasks:
        1. this works as a combination of np.linspace and units.

        Args:
            unit:
            selected unit for the distance to the center

        Raises:
            errors?

        Returns:
        """                           

        import numpy as np
        self.breaks = np.linspace(*args, **kwargs)
        self.breaks = self.breaks * unit
        self.N = len(self.breaks)-1
        self.signal = np.zeros(self.N)
        self.sigma = np.zeros(self.N)
        #}}} 

    def radialprofile(self, center, skymap, skymask):
        #{{{
        """radialprofile(self, skymap) : computes the radial profile of
        CMB pixels around a selected center

        Tasks:
        1. traverse all centers (paralalize here)
        2. traverse all radial bins
        3. traverse all pixels in the ring
        4. compute the mean
        5. store the mean values for all the rings

        Args:
            skymap (class SkyMap):
            Map of the cosmic background, including scalar and mask

            centers_catalog (class Centers):
            Catalog of the centers, including (x, y, z) position
            in Healpix convention and position angle of the galaxy
            disk.

        Raises:
            errors?

        Returns:
            profdata:
            proferror:
            uncertaintydata:
            uncertaintyerror:
        """ 

        # en la version paralela hace un solo centro cada vez
        # que estra a esta funcion
        import numpy as np
        import healpy as hp
        import astropy.units as u
        import time

        radiifloat = self.breaks.to(u.rad)
        listpixs_internal = []
        listpixs_mask = []
        profile = []
        first = True

        for radiusfloat in radiifloat:
            listpixs_external = hp.query_disc(
                skymap.nside,
                center,
                radiusfloat.value,
                inclusive=True,
                fact=4,
                nest=False)

            if(not first):
                listpixs_ring = list(set(listpixs_external) -
                                    set(listpixs_internal))
                listpixs_mask = skymask.data[listpixs_ring]
                mean_ring = np.nanmean(skymap.data[listpixs_ring])
                profile.append(mean_ring)
            first = False
            listpixs_internal = listpixs_external.copy()

        return(profile)
        #}}}
     
    def radialprofile_II(self, centers, skymap, skymask, njobs):
        #{{{
        """radialprofile_II(self, skymap) : computes the radial profile of
        CMB pixels around selected centers in parallel. Uses a wrapper
        and the joblib library.

        Tasks:
        1. traverse all centers (paralalize here)
        2. traverse all radial bins
        3. traverse all pixels in the ring
        4. compute the mean
        5. store the mean values for all the rings

        Args:
            skymap (class SkyMap):
            Map of the cosmic background, including scalar and mask

            centers_catalog (class Centers):
            Catalog of the centers, including (x, y, z) position
            in Healpix convention and position angle of the galaxy
            disk.

        Raises:
            errors?

        Returns:
            profdata:
            proferror:
            uncertaintydata:
            uncertaintyerror:
        """                            
        results = []

        # threading? multiprocessing?
        results = Parallel(n_jobs=njobs, verbose=5, backend="multiprocessing")\
            (delayed(unwrap_profile_self)(i, skymap=skymap, skymask=skymask) 
                    for i in zip([self]*len(centers), centers))

        return(results)
        #}}} 
    #}}}

class AnisotropicProfile:
    # {{{
    '''
    class AnisotropicProfile
    methods for computing angular correlations in the CMB, as a
    function of the angle wrt the position angle of the galaxy and the
    radial distance
    methods:
        set_breaks: select bin scheme for the profile
    ''' 

    def __init__(self, breaks=[0], Nran=0):
        #{{{
        """init(self, breaks, Nran) : sets the partition (binning
        scheme) for the computing of the radial profile.

        Tasks:
        1. this works as a combination of np.linspace and units.

        Args:
            unit:
            selected unit for the distance to the center

        Raises:
            errors?

        Returns:
        """               
        import numpy as np

        self.breaks_rad = breaks
        self.breaks_ang = breaks
        self.Nrad = len(self.breaks_rad) - 1
        self.Nang = len(self.breaks_ang) - 1
        self.max_centers = 0
        #}}} 
 
    def set_breaks_radial(self, unit, *args, **kwargs):
        #{{{
        """set_breaks(self, unit) : sets the breaks for the binned 
        profile.

        Tasks:
        1. this works as a combination of np.linspace and units.

        Args:
            unit:
            selected unit for the distance to the center

        Raises:
            errors?

        Returns:
        """                           

        import numpy as np
        self.breaks_rad = np.linspace(*args, **kwargs)
        self.breaks_rad = self.breaks_rad * unit
        self.Nrad = len(self.breaks_rad)-1
        #}}} 

    def set_breaks_angular(self, unit, *args, **kwargs):
        #{{{
        """set_breaks(self, unit) : sets the breaks for the binned 
        profile.

        Tasks:
        1. this works as a combination of np.linspace and units.

        Args:
            unit:
            selected unit for the distance to the center

        Raises:
            errors?

        Returns:
        """                           

        import numpy as np
        self.breaks_ang = np.linspace(*args, **kwargs)
        self.breaks_ang = self.breaks_ang * unit
        self.Nang = len(self.breaks_ang)-1
        #}}} 

    def anisotropic_profile(self, center, skymap, skymask):
        #{{{
        """radialprofile(self, skymap) : computes the radial profile of
        CMB pixels around a selected center
 
       armar la matriz del perfil (bineado en r y o)

       recorrer la lista de centros

          rotar: armar la matriz de rotacion con angulos de Euler

          recorrer lista de pixels

              calcular el angulo entre la direccion del disco y la direccion al pixel
              calcular la distancia angular entre el centro de la glx y el pixel
              binear esas dos angulos y guardar                           


        """ 

        # en la version paralela hace un solo centro cada vez
        # que estra a esta funcion
        import numpy as np
        import healpy as hp
        import astropy.units as u
        import time
        from scipy.spatial.transform import Rotation as R
        from math import atan2

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

        #}}}
     
    def anisotropic_profile_II(self, centers, skymap, skymask, njobs):
        #{{{
        """radialprofile_II(self, skymap) : computes the radial profile of
        CMB pixels around selected centers in parallel. Uses a wrapper
        and the joblib library.

        Tasks:
        1. traverse all centers (paralalize here)
        2. traverse all radial bins
        3. traverse all pixels in the ring
        4. compute the mean
        5. store the mean values for all the rings

        Args:
            skymap (class SkyMap):
            Map of the cosmic background, including scalar and mask

            centers_catalog (class Centers):
            Catalog of the centers, including (x, y, z) position
            in Healpix convention and position angle of the galaxy
            disk.

        Raises:
            errors?

        Returns:
            profdata:
            proferror:
            uncertaintydata:
            uncertaintyerror:
        """                            
        results = []

        # threading? multiprocessing?
        results = Parallel(n_jobs=njobs, verbose=5, backend="multiprocessing")\
            (delayed(unwrap_anisotropicprofile_self)(i, skymap=skymap, skymask=skymask) 
                    for i in zip([self]*len(centers), centers))

        return(results)
        #}}} 
    #}}}

class Correlation:
    #{{{
    '''
    class Correlation
    methods for computing angular correlations in the CMB
    methods:
        set_breaks: select bin scheme for the profile
        radialprofile: computes the radial profilee
    ''' 

    def __init__(self, breaks=[0], skymask=None, nside=256, nran=0, njobs=1):
        #{{{
        import numpy as np
        import healpy as hp

        self.breaks = breaks
        self.N = len(breaks)-1
        self.nside = nside
        npixs = hp.nside2npix(self.nside)
        self.masked_indices = [i for i in range(npixs) if skymask.data[i]]
        self.N_ma_IDs = len(self.masked_indices)
        self.nran = nran
        self.njobs = njobs

        #}}}
 
    def set_breaks(self, *args, **kwargs):
        #{{{
        import numpy as np
        self.breaks = np.linspace(*args, **kwargs)
        self.N = len(self.breaks)-1
        self.signal = np.zeros(self.N)
        self.sigma = np.zeros(self.N)
        #}}}

    def correlation(self, k, skymap, skymask):
        #{{{
        # en la version paralela hace un solo centro cada vez
        # que estra a esta funcion
        import numpy as np
        import healpy as hp
        import astropy.units as u
        import time
        import random

        random.seed((k+42)*3)
        # draw a set of nran pixel pairs
        x = np.array([random.random() for _ in range(self.nran*2)])
        x = x*self.N_ma_IDs
        x = x.astype(int)
        x = x.reshape((-1, 2))

        # compute frequency
        coso = np.zeros(self.nran)
        tt = np.zeros(self.nran)
        for j, idxs in enumerate(x):
            v1 = hp.pix2vec(self.nside, self.masked_indices[idxs[0]])
            v2 = hp.pix2vec(self.nside, self.masked_indices[idxs[1]])
            coso[j] = np.dot(v1,v2)
            tt[j] = skymap.data[idxs[0]]*skymap.data[idxs[1]]

        H = np.histogram(coso, bins=self.breaks, weights=tt, density=False)[0]
        K = np.histogram(coso, bins=self.breaks, density=False)[0]

        return(H, K)
        #}}}
 
    def correlation_II(self, centers, skymap, skymask):
        #{{{
        results = []

        with Parallel(n_jobs=self.njobs, verbose=5, backend="threading") as P:
            results = P(
                delayed(unwrap_correlation_self)(self, c, skymap=skymap,skymask=skymask)
                for c in centers 
             )

        return(results)
        #}}}
    #}}}

















    def corr(self, centers, tracers, 
             parallel=False, njobs=1):
        """Compute correlations in CMB data.

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
        """
        #{{{
        import numpy as np
        import healpy as hp
        import astropy.units as u
        import time
        from scipy.spatial.transform import Rotation as R
        from math import atan2

        if autocorr: (montecarlo)
            listpixs_C = 
            listpixs_T = 
            pix2vec
            pix2vec
            vec_C =
            vec_T = 

        if profile: (querydisc)
            listpixs_C = 
            listpixs_T = 
            pix2vec
            pix2vec
            vec_C =
            vec_T = 
            #rotation

        # pasar pixs a vecs

        for i in vec_C:
            for i in vec_T:
                dist =
                theta =
                dists.append(dist)
                thetas.append(theta)
                       
        dists = dists * u.rad
        thetas = thetas * u.rad

        rr = self.breaks_rad.to(u.rad)
        aa = self.breaks_ang.to(u.rad)
        bins2d = [rr, self.breaks_ang]

        H = np.histogram2d(dists, thetas, bins=bins2d,
                           weights=temps, density=False)[0]
        K = np.histogram2d(dists, thetas, bins=bins2d,
                           density=False)[0]
                       












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

        #}}}
 
