.. _parsing:

***********************************
Parsing
***********************************
 
All parameters in an experiment have a value assigned through a
configuration or settings file, which is read using the module
:mod:`Parse`, which in turn is based on `configparser <https://docs.python.org/3/library/configparser.html>`_.

There are two main strategies to load the configuration file:

1. From a command line, ``python run_experiment.py config.ini``
2. From the python interpreter

Both strategies can be used with the same code, as follows:


.. code:: python

   from Parser import Parser
   from sys import argv

   if len(argv) > 1:
       config = Parser(argv[1])
   else:
       config = Parser()  

Which loads the default configuration file, set in the variable
``DEFAULT_INI`` in the :mod:``Parser`` module.

Also, a configuration file can be loaded from the python interpreter
directly

.. code:: python

    from Parser import Parser
    config = Parser('custom_file.ini')

Finally, once the configuration object has been instatiated woth the
default variable values, they can be changed with the
:meth:``Parser.load_config`` method, for example:

.. code:: python

    from Parser import Parser
    config = Parser('custom_file.ini')
    n_jobs = 4
    run_parallel = 'yes'
    config.load_config(keys=['n_jobs', 'run_parallel'], 
                       values=[n_jobs, run_parallel])

Variables are accessed by sections.  For example, in order to 
access the variable 'datadir_cmb' in the section 'maps':

.. code:: python

    print(config['maps']['datadir_cmb'])
 

The sections in the configuration file are:

* **experiment**: unique ID for the experiment
* **cmb**: data for the CMB maps
* **glx**: data for the galaxy catalogue
* **run**: computing configuration, options and optimizations
* **out**: output data
* **UX**: user experience


Configuration of required fields
=====================================

The following are the required fields for a basic experiment:

.. code-block::

   [experiment]

   experiment_ID = TST

   [cmb]

   datadir_cmb = ../dat/
   filedata_cmb_nside = 512
   filedata_cmb_mapa = lensmap512_10arcmin_y2.fits
   filedata_field_mapa = 0
   filedata_cmb_mask = lensmask512_10arcmin_y2.fits
   filedata_field_mask = 0

   [glx]

   datadir_glx = ../dat/
   filedata_glx = 2mrs_1175_VAC.dat

   [run]

   theta_start = 0.
   theta_stop = 2*pi
   theta_n_bins = 8
   theta_units = rad      

   r_start = 0.
   r_stop = 40
   r_n_bins = 20
   r_units = angular

   [out]

   save_pickle = False
   dir_output = ../out/
   pickle_name_root = run_
   pickle_name_ext = .pk

   dir_plots = ../plt/
   plot_name_root = corr_
   plot_format = pdf
   clobber = y
   plot_fname = plot
   plot_ftype = PNG

Here we present the detailed description of these fields. 
For all variables which admit Y/N options, the following values are
accepted: 

   * YES (y, yes, s, si, true). Case insensitive.
   * NO (n, no, false). Case insensitive.
    

**experiment_ID**
   An identifier for each experiment.  When running a new experiment,
   directories will be created to store the output results and plots.
   Examples: EXP_001, 01, TEST, etc. (without spaces).
**datadir_cmb**
   Directory where the data files with the CMB temperature maps and
   masks are located.  It accepts absolute paths or relative paths to
   the directory where the ``Parse`` object is executed.
**filedata_cmb_nside**
   Healpix ``nside`` corresponding to the map. If a wrong value is
   set, an error is raised. If a wrong value is set, an error is
   raised.
**filedata_cmb_mapa**
   File with the Healpix map with the temperatures.
**filedata_field_mapa**
   Field in the FITS file with the Healpix map with the temperatures.
**filedata_cmb_mask**
   File with the Healpix mask with the temperatures.
**filedata_field_mask**
   Field in the FITS file with the Healpix mask with the temperatures.
**datadir_glx**
   Directory where the data files with the galaxy catalogues are located.
   It accepts absolute paths or relative paths to
   the directory where the ``Parse`` object is executed.
**filedata_glx**
   File name with the galaxy catalogue.
**theta_start**
   Starting value of the angle with respect to the galaxy disk
**theta_stop**
   Ending value of the angle with respect to the galaxy disk
**theta_n_bins**
   Number of bins in the range [theta_start, theta_stop]
**theta_units**
   Units for the values of theta_start and theta_stop.
   
   Options:

   * rad,
   * arcmin,
   * arcsec.

**r_start**
   Starting value of the angular distance to the center
**r_stop**
   Ending value of the angular distance to the center
**r_n_bins**
   Number of bins in the range [r_start, r_stop]
**r_units**
   Units for the values of r_start and r_stop.
   
   Options:

   * rad: radians
   * arcmin: arc minutes
   * arcsec: arc seconds
   * angular: distance is normalized to the angular size of each
     galaxy
   * physical: distance is normalized to the physical size of each
     galaxy    
 
**save_pickle**
   Wether to save the results in pickle files.
   Options: Y/N

**dir_output**
   Directory of output data files.

**pickle_name_root**
   Root of pickle filename
**pickle_name_ext**
   Extension of pickle filename (e.g., 'pk')
**dir_plots**
   Directory of output data files.

**plot_name_root**
   Root of plot filename
**plot_format**
   Format of plot filename
**clobber**
   Wether to overwrite the output files when repeating experiments.
   Options: Y/N.
**plot_fname**
   root name for the plot
**plot_ftype**
   filetype for the plot


Configuration of optional fields
=====================================

The following are optional fields for a given experiment:

.. code-block::

   [glx]

   max_centers = 150
   control_sample = no
   control_n_samples = 10
   control_ranmap = no
   control_angles = no

   [run]

   n_jobs = 4
   run_parallel = n
    
   r_avg_cuts = 99
   r_avg_fact = 1.

   adaptative_resolution = yes
   adaptative_res_nside = 128
   adaptative_res_dilut = 0 8 15

   disk_align = yes

   galaxy_types = Sb Sc
   redshift_min = 0.001
   redshift_max = 0.020

   ellipt_min = 0.
   ellipt_max = 0.2

   glx_angsize_min = no
   glx_angsize_max = no
   glx_angsize_unit = no

   glx_physize_min = 8.
   glx_physize_max = 99.
   glx_physize_unit = kpc

   [UX]

   show_progress = y
   verbose = y
   interactive = n


Detailed description of the optional fields:

max_centers
   description
control_sample
   description
control_n_samples
   description
control_ranmap
   description
control_angles
   description
n_jobs
   description
run_parallel
   description
r_avg_cuts
   description
r_avg_fact
   description
adaptative_resolution
   description
adaptative_res_nside
   description
adaptative_res_dilut
   description
disk_align
   description
galaxy_types
   description
redshift_min
   description
redshift_max
   description
ellipt_min
   description
ellipt_max
   description
glx_angsize_min
   description
glx_angsize_max
   description
glx_angsize_unit
   description
glx_physize_min
   description
glx_physize_max
   description
glx_physize_unit
   description
show_progress
   description
verbose
   description
interactive
   description

