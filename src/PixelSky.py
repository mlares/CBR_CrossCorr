class SkyUnits():
    #{{{
    '''
    Transform units in the sky for a given cosmology
    '''

    from astropy.cosmology import Planck15 as cosmodel
    from astropy import units as u

    def __init__(self, theta=0, z=0., distance=0., glxsize=0.):
        self.theta = theta
        self.z = z
        self.distance = distance
        self.glxsize = glxsize

    def ang2proj(self, theta, z):

        from astropy.cosmology import Planck15 as cosmodel
        from astropy import units as u
        d_A = cosmodel.angular_diameter_distance(z)
        distance_Mpc = (theta * d_A).to(u.kpc, u.dimensionless_angles())
        return(distance_Mpc)

    def ang2glxsize(self, theta, z, glxsize):
        # in development (always compute??)
        rglxsize=1.
        return(rglxsize)

    def glxsize2ang(self, glxsize):
        # in development (always compute??)
        theta=1.
        return(theta)

    def glxsize2proj(self, size):
        # in development (always compute??)
        Rp = 1.
        return(Rp)

    def proj2ang(self, proj):
        # in development (always compute??)
        theta=1.
        return(theta)

    def proj2glxsize(self, proj):
        # in development (always compute??)
        rglxsize=1.
        return(rglxsize)
    #}}}


class SkyMap():
    # {{{
    '''
    class SkyMap: methods for computing angular correlations in the CMB
    methods:
        load: loads a CMB map
    '''

    import healpy as hp

    def __init__(self, nside=256, ordering='ring', frame='equatorial'):

        import healpy as hp

        self.nside = nside
        self.ordering = 'ring'
        self.frame = 'equatorial'
        self.npixs = hp.nside2npix(self.nside)

        # fac =
        # Dsel: sample of galaxies
        # self, skymap, nside,fac,rprof,Dsel,vec,hp_data_sel,hp_mask):

    def __len__(self):
        return self.npixs

    def __repr__(self):
        return 'Sky Map with {!s} pixels in {!s} order'.format(\
            self.npixs, self.ordering)

    def __str__(self):
        return 'Sky Map with {!s} pixels in {!s} order'.format(\
            self.npixs, self.ordering)

    def load(self, filename, *args, **kwargs):
        '''
        Reads the CMB map

        Args:
            filename (str): the file name of the map to be read

        Raises:

        Returns:
            readmap: a healpix map, class ?

        '''

        import healpy as hp

        d = hp.read_map(filename, h=True, **kwargs)
        self.data = d[0]
        self.header = d[1]

        return(True)
    
    def apply_mask(self, mask):

        import healpy as hp

        m = self[0].copy()
        k = mask[0].copy()
        m[k<0.5] = hp.UNSEEN
        masked_map = hp.ma(m)
        return(masked_map)

    # }}}
         
from joblib import Parallel, delayed

def unwrap_profile_self(arg, **kwarg):
    return RadialProfile.radialprofile(*arg, **kwarg)

def unwrap_correlation_self(arg, **kwarg):
    return Correlation.correlation(*arg, **kwarg)


class RadialProfile:
    # {{{
    '''
    class RadialProfile
    methods for computing angular correlations in the CMB
    methods:
        set_breaks: select bin scheme for the profile
        radialprofile: computes the radial profilee
    ''' 

    def __init__(self, breaks=[0], Nran=0):
        #{{{
        import numpy as np

        self.breaks = breaks
        N = len(breaks)-1
        self.N = N
        self.signal = np.zeros(N)
        self.sigma = np.zeros(N)
        self.controlsample_mean = np.zeros(N)
        self.controlsample_sigma = np.zeros(N)
        #}}} 
 
    def set_breaks(self, unit, *args, **kwargs):
        #{{{

        import numpy as np
        self.breaks = np.linspace(*args, **kwargs)
        self.breaks = self.breaks * unit
        self.N = len(self.breaks)-1
        self.signal = np.zeros(self.N)
        self.sigma = np.zeros(self.N)
        #}}} 

    def radialprofile_stack(self, center, skymap, skymask):
        #{{{
        """radialprofile_stack(self, skymap) : computes the stacked radial profile of
        CMB pixels around selected centers

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
        results = []

        results = Parallel(n_jobs=njobs, verbose=5, backend="threading")\
            (delayed(unwrap_profile_self)(i, skymap=skymap, skymask=skymask) 
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
 
    def set_breaks(self, unit, *args, **kwargs):
        #{{{
        import numpy as np
        self.breaks = np.linspace(*args, **kwargs)
        self.breaks = self.breaks * unit
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

        H = np.histogram(coso, bins=self.breaks, weights=tt, density=True)[0]
        return(H)
        #}}}
 
    def correlation_II(self, centers, skymap, skymask):
        #{{{
        results = []

        results = Parallel(n_jobs=self.njobs, verbose=5, backend="loky")\
            (delayed(unwrap_correlation_self)(i, skymap=skymap,skymask=skymask)
                    for i in zip([self]*len(centers), centers) )

        return(results)
        #}}}

class PixelTools:
    #{{{
    
    def spread_pixels(self, Nside_low, Nside_high, ID):
        #{{{
        """
        returns a list of pixel IDs in the Nside_high resolution
        from a pixel ID in the Nside_low resolution.
        """

        from math import log
        Llow = int(log(Nside_low, 2))
        Lhigh = int(log(Nside_high, 2))

        print(Llow, Lhigh)
        b = bin(ID)
        DN = Lhigh-Llow
        a = [bin(i)[2:].zfill(2**DN) for i in range(4**DN)]
        pix_IDs = []
        for i in a:
            x = (b[2:].zfill(Llow) + i)
            pix_IDs.append(int(x, 2))
        return(pix_IDs)
        #}}}

    def downgrade_pixels(Nside_low, Nside_high, ID):
        #{{{
        """
        Function to be used for testing purposes
        It accomplish the same task than spread_pixels, 
        but with a less efficient, brute force method.
        """
        from math import log
        import healpy as hp
      
        Npix_low = hp.nside2npix(Nside_low)
        Npix_high = hp.nside2npix(Nside_high)
        l = []
        a = []
        for i in range(Npix_high):
            v = hp.pix2vec(Nside_high, i, nest=True)
            k = hp.vec2pix(Nside_low, v[0], v[1], v[2], nest=True)
            if (k==ID):
                l.append(i)
                a.append(bin(i))
        return([l, a])
        #}}}

    #}}}

