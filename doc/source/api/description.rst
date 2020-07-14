***********************************
Description
***********************************

This software offers tools to compute the stacked temperature profile
around a set of selected centers.

A complete experiment can be run with the following interactive code:
 
.. code:: python
                
   >>> import cmfg
   >>> from Parser import Parser
   >>> from sys import argv
   >>> import pickle
   >>> config = Parser('config.ini')
   >>> X = cmfg.profile2d(config)
   >>> X.load_centers()
   >>> X.select_subsample_centers()
   >>> X.load_tracers()
   >>> res = X.run()
                            
or, alternatively, from the command line with:
 
.. code-block::
 
   $ python run_experiment.py settings.ini

where the file ``settings.ini`` contains the code indicated before,
and the ``settings.ini`` file contains the settings for the
experiment.

In what follows, we explain the details of both the code and the
settings.

Minimal example code
=====================

Utilities to compute the profile are contained in the following
modules:

- :mod:`cmfg`
- :mod:`PixelSky`
- :mod:`Parser`
- :mod:`Process`

We also provide example codes:

- for running an experiment: `run_experiment.py`
- for visualizing the results of an experiment: `vis_experiment.py`


Step by step description
========================

First, we need to load the settings from a configuration file:

.. code:: python
                
   >>> from Parser import Parser
   >>> config = Parser('config.ini')

Parameters can also be changed using the function
:meth:`Parse.load_config`.  For more details on loading or changing
the parameters, see :ref:`parsing`.

We first need to instantiate an object of the class
:class:`cmfg.profile2d`.  The initialization requires a configuration
object as argument.

.. code:: python
                
   >>> import cmfg
   >>> X = cmfg.profile2d(config)

An object like this inherits the methods of the class, which include:

- :meth:`cmfg.profile2d.load_centers`
- :meth:`cmfg.profile2d.select_subsample_centers`
- :meth:`cmfg.profile2d.load_tracers`
- :meth:`cmfg.profile2d.run`

First, we need to load the centers and select a subsample, if needed:

.. code:: python
                
   >>> X.load_centers()
   >>> X.select_subsample_centers()

We also need to load the temperature map.  In this case the pixels are
the "tracers" since the profile is equivalent to a cross correlation.

.. code:: python

   >>> X.load_tracers()

Finally, the experiment can be run with the :meth:`cmfg.profile2d.run` method:

.. code:: python

   >>> res = X.run()

In this example, ``res`` is a list that contains:

1. A list of arrays containing the sum of temperatures, each array correspond to one center
2. A list of arrays containing the total number of pixels contributing to the sum of temperatures, each array correspond to one center



..
   NOTES
   -----

   : How to make links in sphinx
   https://sublime-and-sphinx-guide.readthedocs.io/en/latest/references.html

