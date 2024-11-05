:py:mod:`vasca.vasca_pipe`
==========================

.. py:module:: vasca.vasca_pipe

.. autodoc2-docstring:: vasca.vasca_pipe
   :allowtitles:

Module Contents
---------------

Functions
~~~~~~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`set_pipe_dir <vasca.vasca_pipe.set_pipe_dir>`
     - .. autodoc2-docstring:: vasca.vasca_pipe.set_pipe_dir
          :summary:
   * - :py:obj:`set_logger <vasca.vasca_pipe.set_logger>`
     - .. autodoc2-docstring:: vasca.vasca_pipe.set_logger
          :summary:
   * - :py:obj:`keep_base_field <vasca.vasca_pipe.keep_base_field>`
     - .. autodoc2-docstring:: vasca.vasca_pipe.keep_base_field
          :summary:
   * - :py:obj:`run_field <vasca.vasca_pipe.run_field>`
     - .. autodoc2-docstring:: vasca.vasca_pipe.run_field
          :summary:
   * - :py:obj:`run_field_docs <vasca.vasca_pipe.run_field_docs>`
     - .. autodoc2-docstring:: vasca.vasca_pipe.run_field_docs
          :summary:
   * - :py:obj:`run_cluster_fields <vasca.vasca_pipe.run_cluster_fields>`
     - .. autodoc2-docstring:: vasca.vasca_pipe.run_cluster_fields
          :summary:
   * - :py:obj:`run <vasca.vasca_pipe.run>`
     - .. autodoc2-docstring:: vasca.vasca_pipe.run
          :summary:
   * - :py:obj:`run_from_file <vasca.vasca_pipe.run_from_file>`
     - .. autodoc2-docstring:: vasca.vasca_pipe.run_from_file
          :summary:

API
~~~

.. py:function:: set_pipe_dir(vasca_cfg: dict) -> pathlib.Path
   :canonical: vasca.vasca_pipe.set_pipe_dir

   .. autodoc2-docstring:: vasca.vasca_pipe.set_pipe_dir

.. py:function:: set_logger(vasca_cfg)
   :canonical: vasca.vasca_pipe.set_logger

   .. autodoc2-docstring:: vasca.vasca_pipe.set_logger

.. py:function:: keep_base_field(field)
   :canonical: vasca.vasca_pipe.keep_base_field

   .. autodoc2-docstring:: vasca.vasca_pipe.keep_base_field

.. py:function:: run_field(obs_nr, field_id, rg, vasca_cfg)
   :canonical: vasca.vasca_pipe.run_field

   .. autodoc2-docstring:: vasca.vasca_pipe.run_field

.. py:function:: run_field_docs(obs_nr: int, field_id: str, rg: vasca.region.Region, vasca_cfg: dict) -> vasca.field.BaseField | None
   :canonical: vasca.vasca_pipe.run_field_docs

   .. autodoc2-docstring:: vasca.vasca_pipe.run_field_docs

.. py:function:: run_cluster_fields(meanshift_cfg, tt_fd_src, tt_fd_det=None, cluster_coadd=False)
   :canonical: vasca.vasca_pipe.run_cluster_fields

   .. autodoc2-docstring:: vasca.vasca_pipe.run_cluster_fields

.. py:function:: run(vasca_cfg)
   :canonical: vasca.vasca_pipe.run

   .. autodoc2-docstring:: vasca.vasca_pipe.run

.. py:function:: run_from_file()
   :canonical: vasca.vasca_pipe.run_from_file

   .. autodoc2-docstring:: vasca.vasca_pipe.run_from_file
