Ultraviolet Variability Analysis (UVVA)
========================================

UVVA is an astronomy analysis pipeline of ultraviolet sources. Its main purpose is to create catalogs of variable sources. Primarily, this will be done based on GALEX data. 

Installation
------------

1. Clone this repository:

.. code:: bash

   git clone https://gitlab.desy.de/ultrasat-camera/uc_uvva.git
   
2. Setup your environment, either using ``pyenv`` or ``Anaconda`` (see below).

3. Install the ``UVVO`` package:

.. code:: bash

  cd ./uc_uvva
  pip install -e .

4. Setup the resource manager by including your cloud path location in the the ``.env_template`` file and rename it to ``.env``.

Environment setup with pyenv
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Setup a virtual environment
(`pyenv <https://github.com/pyenv/pyenv>`__ and
`pyenv-virtualenv <https://github.com/pyenv/pyenv-virtualenv>`__, install
the requirements and general purpose modules using ``pip``.

Install Python version 3.9.1 and setup a virtual environment (e.g. named
``uc_science_venv39``) after following the pyenv  installation described
`here <https://github.com/pyenv/pyenv#installation>`__ and
`here <https://github.com/pyenv/pyenv-virtualenv#installation>`__:

.. code:: bash

   pyenv virtualenv 3.9.1 uc_science_venv39 

Make sure that the virtual environment is automatically
`enabled <https://github.com/pyenv/pyenv/blob/master/COMMANDS.md#pyenv-local>`__ 
for commands executed from the repository root directory:

.. code:: bash

   cd <your_install_folder>/uc_science   
   pyenv local uc_science_venv39 

Install the necessary Python modules:

.. code:: bash

   pip install -r requirements.txt   
   pip install -e .

Environment setup with Anaconda
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can setup the environment with
`Anaconda <https://www.anaconda.com/products/individual>`__ with the
``conda_env.yml`` file (see
`here <https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#create-env-from-file>`__):

.. code:: bash

   conda env create -f conda_env.yml
   conda activate uc_uvva

Coding guidelines
-----------------

We use the `PEP 8 <https://realpython.com/python-pep8/>`__ coding conventions.
Before contributing please consider the use of automatic code formatting
tools like `isort <https://github.com/pycqa/isort>`__,
`flake8 <https://github.com/PyCQA/flake8>`__ and
`black <https://black.readthedocs.io/en/stable/#>`__. The recommended Python
version to use is 3.9.x . For docstrings we use the
`ņumpy <https://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_numpy.html>`__ 
format.

For documentation we use `SPHINX <https://www.sphinx-doc.org/en/master/>`__. To make them yourself be 
sure to have the ``sphinx-rtd-theme``, ``sphinx.ext.autodoc``
and ``sphinx.ext.napoleon``  installed (see 
`here <https://betterprogramming.pub/auto-documenting-a-python-project-using-sphinx-8878f9ddc6e9>`__ 
for a quick guide).
