class SkyUnits():
    '''
    Transform units in the sky for a given cosmology
    '''

    from astropy.cosmology import Planck15 as cosmodel
    from astropy import units as u

    def ang2proj(theta, z):

        from astropy.cosmology import Planck15 as cosmodel
        from astropy import units as u
        d_A = cosmodel.angular_diameter_distance(z)
        distance_Mpc = (theta * d_A).to(u.kpc, u.dimensionless_angles())
        return(distance_Mpc)

    def ang2glxsize(theta, z, glxsize):
        return(rglxsize)

    def glxsize2ang():
        return(theta)

    def glxsize2proj():
        return(Rp)

    def proj2ang():
        return(theta)

    def proj2glxsize():
        return(rglxsize)


class Centers():
    '''
    class SkyMap: methods for computing angular correlations in the CMB
    methods:
        load: loads a CMB map
    '''

    import healpy as hp

    def __init__(self, pos):
        self.pos = pos


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
        return 'Sky Map with {!s} pixels in {!s} order'.format(self.npixs, self.ordering)

    def __str__(self):
        return 'Sky Map with {!s} pixels in {!s} order'.format(self.npixs, self.ordering)

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
        return(d)
    
    def apply_mask(self, mask):

        import healpy as hp

        m = self[0].copy()
        k = mask[0].copy()
        m[k<0.5] = hp.UNSEEN
        masked_map = hp.ma(m)
        return(masked_map)


    # }}}


class RadialProfile():
    # {{{
    '''
    class RadialProfile
    methods for computing angular correlations in the CMB
    methods:
        set_breaks: select bin scheme for the profile
        radialprofile: computes the radial profilee
    '''
    def __init__(self, breaks=[0], Nran=0):

        import numpy as np

        self.breaks = breaks
        N = len(breaks)-1
        self.signal = np.zeros(N)
        self.sigma = np.zeros(N)

        self.controlsample_mean = np.zeros(N)
        self.controlsample_sigma = np.zeros(N)

    def print(self):
        print(self.breaks)

    def set_breaks(self, unit, *args, **kwargs):

        import numpy as np
        self.breaks = np.linspace(*args, **kwargs)
        self.breaks = self.breaks * unit

    def radialprofile(self, skymap, skymask, centers_catalog):
        """radialprofile(self, skymap) : computes the stacked radial profile of
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
        import numpy as np
        import healpy as hp
        from scipy.stats import sem, nanmean

        # convert radii to radians
        radiifloat = self.breaks.to(u.rad)

        # initialize
        profdata = []

        # despues tratamos de paralelizar esto
        for center in range(glxcat.vec):

            listpixs_internal = []
            listpixs_mask = []

            for radiusfloat in radiifloat:

                listpixs_external = hp.query_disc(
                    skymap.nside,
                    center,
                    radiusfloat,
                    inclusive=True,
                    fact=4,
                    nest=False)

                listpixs_ring = list(set(listpixs_external) -
                                     set(listpixs_internal))
                listpixs_mask = skymask[0][listpixs_ring]
                mean_ring = np.nanmean(skymap[listpix_ring])
                profdata.append(mean_ring)
                listpixs_internal = listpixs_external.copy()

    # }}}
