.. image:: VASCA_icon.png
  :width: 1024
  :alt: VASCA icon

Variable Source Cluster Analysis
================================

The VAriable Source Cluster Analysis (VASCA) is an astronomy pipeline
for time-variable sources. Its main purpose is to create a source catalog
based on GALEX data. A pipeline description can be found in the 
`publication <https://arxiv.org/abs/2405.14269>`__.


Installation
------------

Typically, you want to create a new python 3.10 environment, e.g. with conda or mamba:

.. code:: bash

   mamba create -n vasca python=3.10
   mamba activate vasca
   
The installation steps are:

1. Clone this repository:

.. code:: bash

   git clone git@github.com:rbuehler/vasca.git
 
2. Install the ``VASCA`` package:

.. code:: bash

  cd ./vasca
  pip install -e .

3. Setup the resource manager by including your cloud path location in the the ``.env_template`` file and rename it to ``.env``.


Running the pipeline and post processing
----------------------------------------


We use `Jupyter lab <https://github.com/jupyterlab/jupyterlab>`__ for post processing, functional examples given in ``vasca/examples``.

Coding guidelines
-----------------

We use the `PEP 8 <https://realpython.com/python-pep8/>`__ coding conventions.
Before contributing please consider the use of automatic code formatting
tools like `isort <https://github.com/pycqa/isort>`__,
`flake8 <https://github.com/PyCQA/flake8>`__ and
`black <https://black.readthedocs.io/en/stable/#>`__. We set `88 characters <https://black.readthedocs.io/en/stable/the_black_code_style/current_style.html?highlight=88%20#line-length>`__ as the default line width. The recommended Python
version to use is 3.10.x . For docstrings we use the
`numpy <https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_numpy.html>`__ 
format.

Build the documentation
-----------------------

For documentation we use `SPHINX <https://www.sphinx-doc.org/en/master/>`__.
To build it run the following:

.. code:: bash

	sphinx-apidoc -f -o docs vasca
	cd docs/
	make html

To create Unified Modeling Language diagramms install `pyreverse <https://pylint.pycqa.org/en/latest/pyreverse.html>`__ and `graphviz <https://graphviz.org/>`__ and run:

.. code:: bash

	pyreverse vasca -o png -d ./docs/
