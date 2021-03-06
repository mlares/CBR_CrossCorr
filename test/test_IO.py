
# procedure for testing:
# go to root directory and run $ pip install -e .
# go to tst/ and run pytest


import pytest
import PixelSky as pxs
import numpy as np

def test_radial_profile_N():
    #{{{
    """
    test_radial_profile(self):
        test the construction of the partition for the radial
        profile

    Tasks:

    Args:

    Raises:
        errors?

    Returns:
    """
    import astropy.units as u

    rp = pxs.RadialProfile()
    N = 10
    rp.set_breaks(unit=u.arcmin, start=0., stop=30., num=N)
    assert rp.N == (N-1)
    #}}}

def test_radial_profile_linspace():
    #{{{
    """
    test_radial_profile_linspace(self):
        test the construction of the partition for the radial
        profile using the linspace function from the numpy
        package.

    Tasks:

    Args:

    Raises:
        errors?

    Returns:
    """
    import astropy.units as u

    rp = pxs.RadialProfile()
    N = 10
    rp.set_breaks(unit=u.arcmin, start=0., stop=30., num=N, endpoint=True)
    assert rp.N == (N-1)
    #}}}

def test_load_config_1():
    #{{{
    sys_args = ['python','../set/config_small.ini']
    filename = pxs.check_file(sys_args)
    assert isinstance(filename, str)

def test_load_config_2():
    #{{{
    sys_args = ['../set/config_small.ini']
    with pytest.raises(SystemExit):
        pxs.check_file(sys_args)
 
def test_load_config_3():
    #{{{
    sys_args = ['../set/nonexistentfile.ini']
    with pytest.raises(SystemExit):
        pxs.check_file(sys_args)



    #config = pxs.check_file(['python','../set/config_small.ini']) 
    #with pytest.raises(SystemExit):
    #    config.load_config()
    #}}}

#def test_load_config_2():
#    #{{{
#    config = pxs.check_file(['python','../set/non_existent_file.ini']) 
#    with pytest.raises(SystemExit):
#        config.load_config()
#    #}}}
# 
#def test_load_config_3():
#    #{{{
#    config = pxs.check_file([]) 
#    with pytest.raises(SystemExit):
#        config.load_config()
#    #}}}
# 
#def test_load_config_4():
#    #{{{
#    config = pxs.check_file(["hola", "file", "another"]) 
#    with pytest.raises(SystemExit):
#        config.load_config()
#    #}}}
# 
#def test_load_config_pars():
#    #{{{
#    config = pxs.check_file(['','../set/config_small.ini']) 
#    config.load_config()

    #}}}
 


#
#def test_radial_profile_average(self):
#    #{{{
#    import pandas as pd
#    import healpy as hp
#    import numpy as np
#    from astropy import units as u
#
#    import time
#
#    # Read CMB temperature map
#    nside = 512
#    mapa = pxs.SkyMap(nside)
#    mask = pxs.SkyMap(nside)
#
#    filename = '../dat/lensmap512_10arcmin_y2.fits'
#    mapa.load(filename, field=(0))
#    filename = '../dat/lensmask512_10arcmin_y2.fits'
#    mask.load(filename, field=(0))
#
#    # set the map at a fixed value of 1.
#    value = 1.
#    mapa.data = mapa.data*0. + value
#
#    # Read galaxy catalog
#    glx_catalog = '../dat/2mrs_1175_done.dat'
#    glx = pd.read_csv(glx_catalog, delim_whitespace=True, header=9)
#    # catalog: http://tdc-www.harvard.edu/2mrs/2mrs_readme.html
#
#    phi_healpix = glx['RAdeg']*np.pi/180.
#    theta_healpix = (90. - glx['DECdeg'])*np.pi/180.
#    glx['vec'] = hp.ang2vec(theta_healpix, phi_healpix).tolist()
#
#    # filter galaxy catalog
#    l = ['A','X','B']
#    spiral = [any([s in x for s in l]) for x in glx['type']]
#    edgeon = glx['b/a'] < 0.8
#    subset = spiral & edgeon
#    centers = np.array(list(glx.vec[subset]))
#
#    rp = pxs.RadialProfile()
#    rp.set_breaks(unit=u.arcmin, start=0., stop=30., num=5)
#
#    res = rp.radialprofile_II(centers, mapa, mask) 
#
#    rp.signal = np.mean(res, 1)
#
#    self.assertEqual(rp.signal, value)
#    #}}}
#
#
