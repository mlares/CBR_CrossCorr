class SkyUnits():
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
        rglxsize = 1.
        return(rglxsize)

    def glxsize2ang(self, glxsize):
        # in development (always compute??)
        theta = 1.
        return(theta)

    def glxsize2proj(self, size):
        # in development (always compute??)
        Rp = 1.
        return(Rp)

    def proj2ang(self, proj):
        # in development (always compute??)
        theta = 1.
        return(theta)

    def proj2glxsize(self, proj):
        # in development (always compute??)
        rglxsize = 1.
        return(rglxsize)

    def check_file(sys_args):
        from os.path import isfile

        if len(sys_args) == 2:
            filename = sys_args[1]
            if isfile(filename):
                msg = "Loading configuration parameters from {}"
                print(msg.format(filename))
            else:
                print("Input argument is not a valid file")
                raise SystemExit(1)

        else:
            print('Configuration file expected (just 1 argument)')
            print('example:  python run_experiment.py ../set/config.ini')
            raise SystemExit(1)

        return filename

# https://stackoverflow.com/questions/3609852/which-is-the-best-way-to
#  -allow-configuration-options-be-overridden-at-the-comman
# https://docs.python.org/3/library/configparser.html
# https://tomassetti.me/parsing-in-python/#tools


class SkyMap():
    '''
    class SkyMap: utilities to work with pixelized maps. 
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
        string = 'Sky Map with {!s} pixels in {!s} order'
        return string.format(self.npixs, self.ordering)

    def __str__(self):
        string = 'Sky Map with {!s} pixels in {!s} order'
        return string.format(
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
        m[k < 0.5] = hp.UNSEEN
        masked_map = hp.ma(m)
        return(masked_map)
