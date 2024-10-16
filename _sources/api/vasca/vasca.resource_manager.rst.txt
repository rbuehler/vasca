:py:mod:`vasca.resource_manager`
================================

.. py:module:: vasca.resource_manager

.. autodoc2-docstring:: vasca.resource_manager
   :allowtitles:

Module Contents
---------------

Classes
~~~~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`ResourceManager <vasca.resource_manager.ResourceManager>`
     - .. autodoc2-docstring:: vasca.resource_manager.ResourceManager
          :summary:

Data
~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`CLASS_DIR <vasca.resource_manager.CLASS_DIR>`
     - .. autodoc2-docstring:: vasca.resource_manager.CLASS_DIR
          :summary:
   * - :py:obj:`PACKAGE_DIR <vasca.resource_manager.PACKAGE_DIR>`
     - .. autodoc2-docstring:: vasca.resource_manager.PACKAGE_DIR
          :summary:

API
~~~

.. py:data:: CLASS_DIR
   :canonical: vasca.resource_manager.CLASS_DIR
   :value: 'dirname(...)'

   .. autodoc2-docstring:: vasca.resource_manager.CLASS_DIR

.. py:data:: PACKAGE_DIR
   :canonical: vasca.resource_manager.PACKAGE_DIR
   :value: None

   .. autodoc2-docstring:: vasca.resource_manager.PACKAGE_DIR

.. py:class:: ResourceManager(verbose=False)
   :canonical: vasca.resource_manager.ResourceManager

   .. autodoc2-docstring:: vasca.resource_manager.ResourceManager

   .. rubric:: Initialization

   .. autodoc2-docstring:: vasca.resource_manager.ResourceManager.__init__

   .. py:method:: __enter__()
      :canonical: vasca.resource_manager.ResourceManager.__enter__

      .. autodoc2-docstring:: vasca.resource_manager.ResourceManager.__enter__

   .. py:method:: __exit__(exc_type, exc_value, traceback)
      :canonical: vasca.resource_manager.ResourceManager.__exit__

      .. autodoc2-docstring:: vasca.resource_manager.ResourceManager.__exit__

   .. py:method:: _load_metadata()
      :canonical: vasca.resource_manager.ResourceManager._load_metadata

      .. autodoc2-docstring:: vasca.resource_manager.ResourceManager._load_metadata

   .. py:method:: _load_env(config_path=PACKAGE_DIR + '/.env', overwrite=False)
      :canonical: vasca.resource_manager.ResourceManager._load_env

      .. autodoc2-docstring:: vasca.resource_manager.ResourceManager._load_env

   .. py:method:: _log_env_status()
      :canonical: vasca.resource_manager.ResourceManager._log_env_status

      .. autodoc2-docstring:: vasca.resource_manager.ResourceManager._log_env_status

   .. py:method:: _check_resource_catalog()
      :canonical: vasca.resource_manager.ResourceManager._check_resource_catalog

      .. autodoc2-docstring:: vasca.resource_manager.ResourceManager._check_resource_catalog

   .. py:method:: _check_resoruce_env_info()
      :canonical: vasca.resource_manager.ResourceManager._check_resoruce_env_info

      .. autodoc2-docstring:: vasca.resource_manager.ResourceManager._check_resoruce_env_info

   .. py:method:: get_path(resource, storage)
      :canonical: vasca.resource_manager.ResourceManager.get_path

      .. autodoc2-docstring:: vasca.resource_manager.ResourceManager.get_path
