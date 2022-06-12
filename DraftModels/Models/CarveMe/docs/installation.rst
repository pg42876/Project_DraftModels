.. highlight:: shell

============
Installation
============

**CarveMe** currently supports Python 3.6 and 3.7.
We recommend using the `Anaconda python
distribution <https://www.continuum.io/downloads>`_.

CarveMe can be easily installed using the **pip** package manager:

.. code-block:: console

    $ pip install carveme

Additionally, you must manually install two external dependencies:

- diamond_ (conda install -c bioconda diamond)
- IBM CPLEX_ Optimizer

.. _diamond: https://github.com/bbuchfink/diamond
.. _CPLEX: https://www.ibm.com/analytics/cplex-optimizer

Note that you will need to register with IBM to obtain an `academic license for CPLEX <https://www.ibm.com/academic/home>`_.

**IMPORTANT:** After installing CPLEX, do not forget to install the CPLEX python API (see the CPLEX documentation_ for details).

.. _documentation: https://www.ibm.com/support/knowledgecenter/SSSA5P_12.7.1/ilog.odms.cplex.help/CPLEX/GettingStarted/topics/set_up/Python_setup.html

Everything should be ready now! See the next section for instructions on how to start *carving*.
