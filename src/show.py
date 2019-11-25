import pickle
import configparser

config = configparser.ConfigParser()
config.read('../set/config.ini')

filedata = config['out']['pickle_name_root']+
           config['out']['pickle_name_exp']+
           config['out']['pickle_name_idx']+'.p'


pickle.dump( rp, open( filedata, "r" ) )


