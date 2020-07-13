***********************************
source codes by ML (former version)
***********************************

These code files are stored in *src3* directory

La documentación se genera sola a partir de los comentarios, y se
puede acceder mediante el link al modulo :class:`PixelSky`


La idea es que PixelSky.py define una serie de clases para trabajar
con orientación a objetos.

Por ejemplo, está la clase :class:`PixelSky.RadialProfile` y la clase :class:`PixelSky.Correlation`.

Una vez que se define una instancia de un objeto de una determinada
clase, se pueden aplicar los métodos de esa clase.  Por ejemplo, un
objeto de tipo :class:`PixelSky.RadialProfile` puede ejecutar un
método que fija la partición (i.e., inicio y fin, y cantidad de bines
para el perfil), :meth:`PixelSky.RadialProfile.set_breaks`, y puede
ejecutar la cuenta con :meth:`PixelSky.RadialProfile.radialprofile`.

Para trabajar con el perfil radial, por ejemplo, se puede usar esta
estructura básica (basado en run_profile.py):


.. code-block::

   import PixelSky as pxs
   rp = pxs.RadialProfile()
   rp.set_breaks(unit=u.arcmin, start=start, stop=stop, num=nbins+1)
   res = rp.radialprofile_II(centers, mapa, mask, njobs)


Además de esto, hay otras tareas, a saber:

* parsing del archivo de configuración
* lectura de datos
* escritura de los resultados


Hay dos programas, run_correlation.py y run_profile.py que corren la
correlación y el perfil, respectivamente.

El análisis de los resultados se hace con analyze_corr.py.  Los demás
archivos son de desarrollo (y por lo tanto, en realidad no deberían
estar en control de versión, pero bueh...).


Ahora veamos cada uno en detalle, por ejemplo para el perfil:


Parsing del archivo de configuración
====================================

Todos los parámetros a los que se les asigna un valor en el archivo de
configuración se deben leer usando el módulo `configparser <https://docs.python.org/3/library/configparser.html>`_.


.. code-block::

   import configparser
   config = configparser.ConfigParser()
   config.read('../set/config_small.ini')    

Las variables se acceden a través de las secciones y los nombres
asignados en el archivo de configuración.  Por ejemplo, para acceder
a la variable 'datadir_cmb' de la sección 'maps', *config['maps']['datadir_cmb']*


Lectura de datos
====================================

datos del CMB:

.. code-block::

   filedata = config['maps']['datadir_cmb']+ config['maps']['filedata_cmb_mapa']
   mapa.load(filedata, field=( int(config['maps']['filedata_field_mapa']) ))

   filedata = config['maps']['datadir_cmb']+ config['maps']['filedata_cmb_mask']
   mask.load(filedata, field=( int(config['maps']['filedata_field_mask']) ))

datos de los catálogos de galaxias:

.. code-block::

   glx_catalog = config['cats']['datadir_glx']+config['cats']['filedata_glx']
   glx = pd.read_csv(glx_catalog, delim_whitespace=True, header=9)         


Curación de datos
====================================

.. code-block::

   # positions...
   phi_healpix = glx['RAdeg']*np.pi/180.
   theta_healpix = (90. - glx['DECdeg'])*np.pi/180.
   glx['vec'] = hp.ang2vec(theta_healpix, phi_healpix).tolist()

   filter galaxy catalog...
   l = ['A','X','B']
   spiral = [any([s in x for s in l]) for x in glx['type']]
   edgeon = glx['b/a'] < 0.8
   subset = spiral & edgeon
   centers = np.array(list(glx.vec[subset]))

   centers = np.array(list(glx.vec))

Cómputo del perfil
====================================

Para el cómputo del perfil se asignan los valores de los parámetros
usando el archivo de configuración.

.. code-block::

   # crear el objeto tipo "perfil radial"
   rp = pxs.RadialProfile()

   # configurar el bineado usando los parámetos 
   # del archivo de configuración:
   nbins = int(config['run']['rp_n_bins']) 
   start = float(config['run']['rp_start']) 
   stop = float(config['run']['rp_stop']) 
   rp.set_breaks(unit=u.arcmin, start=start, stop=stop, num=nbins+1)

   # fijar los parámetros de paralelismo
   njobs = int(config['run']['n_jobs']) 

   # hacer el cómputo (función rp.radialprofile_II)
   res = rp.radialprofile_II(centers, mapa, mask, njobs)

   # los resultados están en el objeto "rp"
   rp.signal = np.mean(res, 1)
   rp.sigma = np.std(res, 1)


Escritura de los resultados
====================================

Los resultados se escriben si config['out']['save_pickle'] es True.
El nombre del archivo de salida se construye a partir de los valores
guardados en en archivo de configuración.

.. code-block::

   import pickle
            
   if config['out']['save_pickle']:
       filedata = config['out']['output_dir']+\
                  config['out']['pickle_name_root']+\
                  config['out']['pickle_name_exp']+\
                  config['out']['pickle_name_idx']+'.p'
        
       pickle.dump( rp, open( filedata, "wb" ) )



Paralelismo
====================================

El paralelismo está implementado en el método
:meth:`PixelSky.RadialProfile.radialprofile_II`, mediante un wrapper
de la función serial :meth:`PixelSky.RadialProfile.radialprofile`.
El wrapper es el método :meth:`PixelSky.RadialProfile.unwrap_profile_self`, que usa el paquete `joblib <https://joblib.readthedocs.io/en/latest/>`_.

.. code-block::

   class RadialProfile:
      ...
      from joblib import Parallel, delayed

      def unwrap_profile_self(arg, **kwarg):
          return RadialProfile.radialprofile(*arg, **kwarg)
                                                            

      def radialprofile_II(self, centers, skymap, skymask, njobs):
        results = []
        results = Parallel(n_jobs=njobs, verbose=5, backend="threading")\
            (delayed(unwrap_profile_self)(i, skymap=skymap, skymask=skymask) 
                    for i in zip([self]*len(centers), centers))

        return(results)




