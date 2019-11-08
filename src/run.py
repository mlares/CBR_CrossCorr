import PixelSky as pxs

import pandas as pd

# Read CMB temperature map
mapa = pxs.SkyMap(nside)
mask = pxs.SkyMap(nside)

filename = '../dat/lensmap512_10arcmin_y2.fits'
mapa = mp.load(filename, field=(0))
filename = '../dat/lensmask512_10arcmin_y2.fits'
mask = mp.load(filename, field=(0))

# mapa.apply_mask(mask) # no funciona??



# Read galaxy catalog
glx_catalog = '../dat/2mrs_1175_done.dat'
glx = pd.read_csv(glx_catalog, delim_whitespace=True, header=9)
# catalog: http://tdc-www.harvard.edu/2mrs/2mrs_readme.html

# filter galaxy catalog
l = ['A','X','B']
spiral = [any([s in x for s in l]) for x in glx['type']]
edgeon = glx['b/a'] < 0.8
subset = spiral & edgeon

rp = pxs.RadialProfile()
rp.set_breaks(unit=1., start=0., stop=30., num=10)









print('OK')



#  
# def f1(a, b, c=0, *args, **kwargs):
#     """Esto es la documentación de la función 1"""
#     partes = ["Parámetros pasados a f1", 
#         f"a={a} b={b} args={args} kwargs={kwargs}"]
#     return " -- ".join(partes)
# 
# class clase():
#     def f2(c, d, *args, **kwargs):
#         """Esto es la documentación de la función 2"""
#         partes = ["Parámetros pasados a f2", 
#             f"c={c} d={d} args={args} kwargs={kwargs}"]
#         print(" -- ".join(partes))
# 
#         return(f1(c, d, *args, **kwargs))
# 
# print(clase.f2(1, 2, 4, l=9))
# print('OK')
# print('--------------')
# 
# 
# # def f(a, b, c=1):
# #     print(a+b+c)
# 
#  
