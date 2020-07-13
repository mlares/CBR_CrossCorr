# load_ext autoreload
# autoreload 2

import cmfg
from Parser import Parser
from sys import argv
import pickle

if len(argv) > 1:
    config = Parser(argv[1])
else:
    config = Parser()

X = cmfg.profile2d(config)
X.load_centers()
X.select_subsample_centers()

print(f"Number of centers in subsample: {X.centers.shape[0]}")

X.load_tracers()
res = X.run()
H, K = res

fout = (f"{config.p.dir_output}{config.p.experiment_id}"
        f"/R_{config.p.experiment_id}.pk")
pickle.dump(res, open(fout, 'wb'))
