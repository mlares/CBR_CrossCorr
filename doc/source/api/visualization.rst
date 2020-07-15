***********************************
Visualization
***********************************


Post processing and visualizing
===============================

.. code:: python

   from Parser import Parser
   import pickle
   from Process import rebin2d, rebin1d, profiles, rt_axes

   config = Parser('settings.ini')
   f_input = (f"{config.p.dir_output}{config.p.experiment_id}"
              f"/profile_{config.p.experiment_id}.pk")

   with open(f_input, 'rb') as f:
       H, K = pickle.load(f)

   r_breaks, r_means, t_breaks, t_means = rt_axes(config)
   res = profiles(H, K, config)
   mean_dT_cells, prof_avg, prof_stack, prof_para, prof_perp = res
     


..
   NOTES
   -----

   : How to make links in sphinx
   https://sublime-and-sphinx-guide.readthedocs.io/en/latest/references.html

