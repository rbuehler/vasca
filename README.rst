Ultraviolet Variability Analysis
================================

The Ultraviolet Variability Analysis (UVVA) is an astronomy pipeline
for time-variable sources. Its main purpose is to create a source catalog
based on GALEX data.

.. image:: https://gitlab.desy.de/ultrasat-camera/uc_uvva/badges/main/pipeline.svg
    :target: https://gitlab.desy.de/ultrasat-camera/uc_uvva/-/commits/main
    :alt: pipeline status
    
.. image:: https://gitlab.desy.de/ultrasat-camera/uc_uvva/badges/main/coverage.svg
    :target: https://gitlab.desy.de/ultrasat-camera/uc_uvva/-/commits/main
    :alt: coverage report



Installation
------------

1. Clone this repository:

.. code:: bash

   git clone https://gitlab.desy.de/ultrasat-camera/uc_uvva.git
   
2. Setup your environment, either using ``pyenv`` or ``Anaconda`` (see below).

3. Install the ``UVVA`` package:

.. code:: bash

  cd ./uc_uvva
  pip install -e .

4. Setup the resource manager by including your cloud path location in the the ``.env_template`` file and rename it to ``.env``.

**Environment setup with pyenv**

Setup a virtual environment
(`pyenv <https://github.com/pyenv/pyenv>`__ and
`pyenv-virtualenv <https://github.com/pyenv/pyenv-virtualenv>`__, install
the requirements and general purpose modules using ``pip``.

Install Python version 3.9.9 and setup a virtual environment (e.g. named
``uc_uvva_venv39``) after following the pyenv  installation described
`here <https://github.com/pyenv/pyenv#installation>`__ and
`here <https://github.com/pyenv/pyenv-virtualenv#installation>`__:

.. code:: bash

   pyenv virtualenv 3.9.9 uc_uvva_venv39 

Make sure that the virtual environment is automatically
`enabled <https://github.com/pyenv/pyenv/blob/master/COMMANDS.md#pyenv-local>`__ 
for commands executed from the repository root directory:

.. code:: bash

   cd <your_install_folder>/uv_uvva   
   pyenv local uc_uvva_venv39 

Install the necessary Python modules:

.. code:: bash

   pip install -r requirements.txt   
   pip install -e .

**Environment setup with Anaconda**

You can setup the environment with
`Anaconda <https://www.anaconda.com/products/individual>`__ with the
``conda_env.yml`` file (see
`here <https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#create-env-from-file>`__):

.. code:: bash

   conda env create -f conda_env.yml
   conda activate uc_uvva

**Prototyping and functional examples**

We use `Jupyter lab <https://github.com/jupyterlab/jupyterlab>`__ for prototyping and for functional examples given in ``uvva/examples``.
Follow these `instructions <https://albertauyeung.github.io/2020/08/17/pyenv-jupyter.html/>`__ to add  a pyenv-generated virtual environment as a Jupyter kernel. Jupyter lab extensions can be used to enable interactive Matoplotlib figures: First install the `ipywidgets <https://github.com/jupyter-widgets/ipywidgets>`__ extension and then the `ipympl <https://github.com/matplotlib/ipympl>`__ extension.

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

For documentation we use `SPHINX <https://www.sphinx-doc.org/en/master/>`__. To make them yourself be 
sure to have the ``sphinx-rtd-theme``, ``sphinx.ext.autodoc``
and ``sphinx.ext.napoleon``  installed (see 
`here <https://betterprogramming.pub/auto-documenting-a-python-project-using-sphinx-8878f9ddc6e9>`__ 
for a quick guide).
