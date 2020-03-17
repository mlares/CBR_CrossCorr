import numpy as np
from astropy import coordinates as coo
from astropy.modeling import rotations as Rot
from astropy import units as u
from scipy.spatial.transform import Rotation as R


# coo.spherical_to_cartesian(1., delta, alpha)
#Note that the input angles should be in latitude/longitude or
#elevation/azimuthal form.  I.e., the origin is along the equator
#rather than at the north pole.
 
# Rot.EulerAngleRotation(lon, lat, phi, 'zyz')
#Rotates one coordinate system into another (fixed) coordinate system.
#All coordinate systems are right-handed. The sign of the angles is
#determined by the right-hand rule.

# Rot.spherical2cartesian(lon, lat) ?
 
# alpha in [0, 360] deg
alpha = np.random.rand()*360.*u.deg

# delta in [-90, 90] deg
delta = math.acos(random()*2.-1.)*u.rad
delta = 90.*u.deg - delta.to(u.deg) 
phi = 90.*u.deg

v = coo.spherical_to_cartesian(1., delta, alpha)
r1 = R.from_euler('z', -alpha.value, degrees=True)
r2 = R.from_euler('y', -90.+delta.value, degrees=True)
r3 = R.from_euler('z', -phi.value, degrees=True)


v1 = r1.apply(v)
v2 = r2.apply(v1)
v_new = r3.apply(v2)


 


#r = R.from_euler('zxz', [lon, lat, phi], degrees=True)

r1 = R.from_euler('z', lon, degrees=True)  
r2 = R.from_euler('x', lat, degrees=True)  
r = r2 * r1

v_new = r.apply(v)
print(v_new)



#from math import atan2
#
## calcular la matriz de rotacion    
#phi = float(center.phi)
#theta = float(center.theta)
#pa = float(center.pa)
#
#r = R.from_euler('zxz', [phi, theta, pa], degrees=True)
 
 







