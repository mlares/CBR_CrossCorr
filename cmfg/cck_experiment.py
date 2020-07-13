# load_ext autoreload
# autoreload 2

from Parser import Parser
import pickle
from matplotlib import pyplot as plt

config = Parser('../set/Ts05.ini')
f_input = (f"{config.p.dir_output}{config.p.experiment_id}"
           f"/R_{config.p.experiment_id}.pk")
with open(f_input, 'rb') as f:
    H1, K1 = pickle.load(f)

config = Parser('../set/Ts07.ini')
f_input = (f"{config.p.dir_output}{config.p.experiment_id}"
           f"/R_{config.p.experiment_id}.pk")
with open(f_input, 'rb') as f:
    H2, K2 = pickle.load(f)

plt.plot(H1[0], color='firebrick', linewidth=2, label='H1')
plt.plot(H2[0], color='cadetblue', linewidth=2, label='H2')
plt.plot(H1[0]-H2[0], color='slategrey', linewidth=3, label='H1-H2')

plt.legend()
plt.show()


# from PixelSky import PixelTools
# import numpy as np
# import healpy as hp
# px = PixelTools()
# # esto funciona:
# l = hp.query_disc(16, [1, 0, 0], 0.07, nest=True)
# p16 = px.spread_pixels(16, 32, l[0], order='nest')
# p32 =hp.query_disc(32, [1, 0, 0], 0.07, nest=True)
# # ahora en RING:
# l = hp.query_disc(16, [1, 0, 0], 0.07, nest=False)
# p16 = px.spread_pixels(16, 32, l[0], order='ring')
# p32 =hp.query_disc(32, [1, 0, 0], 0.07, nest=False)
