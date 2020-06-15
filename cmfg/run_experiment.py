# load_ext autoreload
# autoreload 2

import cmfg
from Parser import Parser
from sys import argv
import numpy as np
import pickle

if len(argv) > 1:
    config = Parser(argv[1])
else:
    config = Parser()

X = cmfg.Correlation(config)
X.load_centers()
X.load_tracers()
X.select_subsample_centers()
R = X.run()

fout = f"R_{config.p.experiment_id}.pk"
pickle.dump(R, open(fout, 'wb'))


#from matplotlib import pyplot as plt
#fig = plt.figure()
#ax = fig.add_subplot()
#
#Rt = R.transpose()
#for r in Rt:
#    ax.plot(r)
#plt.show()


   
#   # filter galaxy catalog...
#   # l = ['A','X','B']
#   # spiral = [any([s in x for s in l]) for x in glx['type']]
#   # edgeon = glx['b/a'] < 0.8
#   # subset = spiral & edgeon
#   # centers = np.array(list(glx.vec[subset]))
#   
#   centers = pd.DataFrame(list(zip(phi_healpix, theta_healpix, glx['pa'])),
#             columns=['phi','theta', 'pa'] )
