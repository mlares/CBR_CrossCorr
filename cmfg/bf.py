import cmfg
from Parser import Parser
import numpy as np
import healpy as hp
from scipy.spatial.transform import Rotation as R

config = Parser()

X = cmfg.Correlation(config)
X.load_centers() 

center = [0]
center.append(X.centers.iloc[12])

phi = float(center[1].phi) 
theta = float(center[1].theta) 
pa = float(center[1].pa)

rmax = 0.12
vector = hp.ang2vec(center[1].theta, center[1].phi)
rotate_pa = R.from_euler('zyz', [-phi, -theta, pa])

listpixs = hp.query_disc(512, vector, rmax, inclusive=False, 
                         fact=4, nest=False)

for ipix in listpixs:
    v = hp.pix2vec(512, ipix)
    w = rotate_pa.apply(v)
    dist = hp.rotator.angdist(w, [0, 0, 1])
    #print(math.acos(np.dot(w, [0,0,1]))/rmax, dist[0]/rmax)
    print(dist[0]/rmax)

