# para generar un archivo de una catálogo random (de "galaxias")

import numpy as np
from math import acos
from random import random
import pandas as pd

N = 100000

x = [random()*360. for _ in range(N)]

y = [acos(random()*2.-1.) for _ in range(N)]

df = pd.DataFrame(zip(x,y), columns=['RAdeg','DECdeg'] )

df.to_csv('../dat/random_catalog.dat', index=False)


