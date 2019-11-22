import PixelSky as pxs
import numpy as np

def test_radial_profile_N(self):
    rp = pxs.RadialProfile()
    N = 10
    rp.set_breaks(unit=u.arcmin, start=0., stop=30., num=N)
    self.assertEqual(len(rp.N), N)

   
"""
TO DO:

- test passing different arguments for linspace
- test passing different units
- test assignment of RadialProfile attributes
- ...

"""
