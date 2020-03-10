***********
Testing
***********

Make several tests to reach a good code coverage and verify if results are as expected.

Proposed tests to develop
=========================

  * check parsing
  * read a synthetic map and verify a flat profile.
  * Check if CMB maps are read correctly
  * Check paths and files
  * read a sample of random centers and compute an averaged flat profile
  * compute the profile of just one center (fixed?)
  * checl size of CMB map and Nside
  * compare results of serial and parallel versions.

  * test passing different arguments for linspace
  * test passing different units
  * test assignment of RadialProfile attributes


Testing tools and procedures
============================

In order to make testing, we should use any of the following tools:

* `pytest <https://docs.pytest.org/en/latest/>`_
* hypothesis


**pytest examples**

"pytest will run all files of the form test_\*.py or \*_test.py in the current directory and its subdirectories."
So, simply go to the tst directory, and run pytest.

In the environment:

.. code-block::

   pip install pytest


In the code (example from pytest documentation):

.. code-block::

   def inc(x):
       return x + 1


   def test_answer():
       assert inc(3) == 5

How to run test:

From the CLI, write:

.. code-block::

   pytest

** coverage **

Desde el entorno, instalar los paquetes coverage y pytest-cov:

.. code-block::

   pip install coverage pytest-cov

Para calcular la cobertura de código, correr:


.. code-block::

   coverage run tst/test_IO.py
   coverage report


Se puede integrar el pytest con el codecov:

.. code-block::

   pytest tst/test_IO.py --cov test_cov --cov-fail-under 90







