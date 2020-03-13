


def AsymmetricProfile(centros, angulos):
   """
   armar la matriz del perfil (bineado en r y o)

   recorrer la lista de centros

      rotar: armar la matriz de rotacion con angulos de Euler

      recorrer lista de pixels

          calcular el angulo entre la direccion del disco y la direccion al pixel
          calcular la distancia angular entre el centro de la glx y el pixel
          binear esas dos angulos y guardar
   """
   from scipy.spatial.transform import Rotation as R

   # recorrer la lista de centros
   Ht = np.zeros(nbins_x, nbins_y)
   Kt = np.zeros(nbins_x, nbins_y)

   for i in centros:

       # calcular la matriz de rotacion

       
        alpha = 45.
        delta = 0.
        q=-45.
        r = R.from_euler('zxz', [alpha, delta, q], degrees=True)


       # querydisc
       listpixs = hpx.query_dics(...)

       dists = []
       thetas = []
       temps = []

       for ipix in listpixs:

           v = hpx.pix2vec(ipix)
           w = r.apply(v)

           dist = hpx.angdist(w, [0,0,1])
           theta = np.atan2(w[1], w[0])

           dists.append(dist)
           thetas.append(theta)
           temps.append(T[ipix])

       H = np.histogram2D(dists, thetas, bins=bins, weights=temps, density=False)[0]
       K = np.histogram2D(dists, thetas, bins=bins, density=False)[0]

       # stacking
       Ht = Ht + H.H
       Kt = Kt + K.H


    # mapa estackeado:
    mapa = Ht / Kt

    # para graficar:
    H.xedges, H.yedges, mapa





