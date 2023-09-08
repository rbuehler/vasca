Variable Source Cluster Analysis
================================

The VAriable Source Cluster Analysis (VASCA) is an astronomy pipeline
for time-variable sources. Its main purpose is to create a source catalog
based on GALEX data.

.. image:: https://gitlab.desy.de/ultrasat-camera/vasca/badges/main/pipeline.svg
    :target: https://gitlab.desy.de/ultrasat-camera/vasca/-/commits/main
    :alt: pipeline status
    
.. image:: https://gitlab.desy.de/ultrasat-camera/vasca/badges/main/coverage.svg
    :target: https://gitlab.desy.de/ultrasat-camera/vasca/-/commits/main
    :alt: coverage report

Installation
------------

Typically, you want to create a new python 3.10 environment, e.g. with conda or mamba:

.. code:: bash

   mamba create -n vasca python=3.10
   mamba activate vasca
   
The installation steps are:

1. Clone this repository:

.. code:: bash

   git clone git@gitlab.desy.de:ultrasat-camera/vasca.git
 
2. Install the ``VASCA`` package:

.. code:: bash

  cd ./vasca
  pip install -e .

3. Setup the resource manager by including your cloud path location in the the ``.env_template`` file and rename it to ``.env``.


**Prototyping and functional examples**

We use `Jupyter lab <https://github.com/jupyterlab/jupyterlab>`__ for prototyping and for functional examples given in ``vasca/examples``.
Follow these `instructions <https://albertauyeung.github.io/2020/08/17/pyenv-jupyter.html/>`__ to add  a pyenv-generated virtual environment as a Jupyter kernel. Jupyter lab extensions can be used to enable interactive Matoplotlib figures: First install the `ipywidgets <https://github.com/jupyter-widgets/ipywidgets>`__ extension and then the `ipympl <https://github.com/matplotlib/ipympl>`__ extension. For plotting of sky maps in notebooks we often use `Imviz <https://jdaviz.readthedocs.io/en/latest/imviz/index.html>`__.

Coding guidelines
-----------------

We use the `PEP 8 <https://realpython.com/python-pep8/>`__ coding conventions.
Before contributing please consider the use of automatic code formatting
tools like `isort <https://github.com/pycqa/isort>`__,
`flake8 <https://github.com/PyCQA/flake8>`__ and
`black <https://black.readthedocs.io/en/stable/#>`__. We set `88 characters <https://black.readthedocs.io/en/stable/the_black_code_style/current_style.html?highlight=88%20#line-length>`__ as the default line width. The recommended Python
version to use is 3.9.x . For docstrings we use the
`ņumpy <https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_numpy.html>`__ 
format.

**Documentation**

For documentation we use `SPHINX <https://www.sphinx-doc.org/en/master/>`__. To make them yourself be 
sure to have the ``sphinx-rtd-theme``, ``sphinx.ext.autodoc``
and ``sphinx.ext.napoleon``  installed (see 
`here <https://betterprogramming.pub/auto-documenting-a-python-project-using-sphinx-8878f9ddc6e9>`__ 
for a quick guide).
Here the steps to follow:

.. code:: bash

	sphinx-apidoc -f -o docs vasca
	cd docs/
	make html

To create Unified Modeling Language diagramms install `pyreverse <https://pylint.pycqa.org/en/latest/pyreverse.html>`__ and `graphviz <https://graphviz.org/>`__ and run:

.. code:: bash

	pyreverse vasca -o png -d ./docs/
