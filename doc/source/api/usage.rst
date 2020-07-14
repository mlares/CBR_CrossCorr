*********
API Usage
*********

- instalation through pypy not yet implemented
- make setup.py installer
- from a python script, call import PixelSky


This project is organized as an API to be used from a python prompt.

Steps:

- Complete the configuration of the experiment
- All the settings of the experimets are parsed from the configuration
  files using configparser.


Prerequisites
=============

* Put data files on the ``dat`` directory.
* Complete the names of the data files in the configuration file

Data
====

Data is stored in the *dat* directory.


=========================================  =================================================
 filename                                   contents
=========================================  =================================================
2mrs_1175_done.dat                          CMB temperature map
COM_CMB_IQU-smica_2048_R3.00_full.fits      CMB temperature map
lensmap512_10arcmin_y2.fits*                CMB lensing map
lensmask512_10arcmin_y2.fits*               CMB lensing mask map
=========================================  =================================================





Configuration files
===================


.. code-block::
   [DEFAULT]

   [maps]

   # CMB MAPS
   datadir_cmb = ../dat/
   filedata_cmb_nside = 512

   filedata_cmb_mapa = lensmap512_10arcmin_y2.fits
   filedata_field_mapa = 0

   # masks
   filedata_cmb_mask = lensmask512_10arcmin_y2.fits
   filedata_field_mask = 0


   [cats]
   # CATALOGUES
   datadir_glx = ../dat/
   filedata_glx = 2mrs_1175_done.dat


   [run]
   # CONFIGURATIONS FOR EXPERIMENT AND COMPUTATIONS

   # to be passed to joblib 
   n_jobs = 2
   # breaks for radial profile
   rp_n_bins = 10
   rp_start = 0.
   rp_stop = 100.
   # breaks for correlation
   corr_n_bins = 77
   corr_start = -1.
   corr_stop = 1.
   # MonteCarlo for correlation
   Nran = 1000
   Nexperiments = 10

   [out]
   # OUTPUT SETTINGS

   save_pickle = True
   output_dir = ../out/
   pickle_name_root = rp_run_
   pickle_name_exp = nobjs15_
   pickle_name_idx = 01




Interactive usage
=================

For a simple test, go to cmfg and run:

.. code-block::

   $ python run_profile.py ../set/config_small.ini


Run experiments at IATE
=======================

In order to use the `HPC services at IATE <https://wiki.oac.uncor.edu/doku.php>`_ the following steps shoul be followed:


1. log in into a cluster (e.g., ``ssh clemente``)
2. git clone or pull the `CBR_correlation <https://github.com/mlares/CBR_CrossCorr>`_ project.
3. prepare a SLURM script (src/submit_python_jobs.sh)
4. launch the script: ``sbatch submit_python_jobs.sh``


SLURM script example for *clemente* running python in parallel:

.. code-block::
   #!/bin/bash

   # SLURM script for: CLEMENTE
    
   ## Las líneas #SBATCH configuran los recursos de la tarea
   ## (aunque parezcan estar comentadas)

   # More info:
   # http://homeowmorphism.com/articles/17/Python-Slurm-Cluster-Five-Minutes


   ## Nombre de la tarea
   #SBATCH --job-name=CMB_corr

   ## Cola de trabajos a la cual enviar.
   #SBATCH --partition=small

   ## tasks requested
   #SBATCH --ntasks=1
   #SBATCH --cpus-per-task=20

   ## STDOUT
   #SBATCH -o submit_python_jobs.out

   ## STDOUT
   #SBATCH -e submit_python_jobs.err

   ## Tiempo de ejecucion. Formato dias-horas:minutos.
   #SBATCH --time 0-1:00

   ## Script que se ejecuta al arrancar el trabajo

   ## Cargar el entorno del usuario incluyendo la funcionalidad de modules
   ## No tocar
   . /etc/profile

   # conda init bash
   # source /home/${USER}/.bashrc

   module load gcc/8.2.0
   conda activate
   # por las dudas activar conda antes de correr el sbatch

   ## Launch program

   srun python /home/mlares/CBR_CrossCorr/src/run_correlation.py ../set/config_big.ini

   ## launch script
   ## $>sbatch submit_python_jobs.sh







