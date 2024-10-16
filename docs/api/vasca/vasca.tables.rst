:py:mod:`vasca.tables`
======================

.. py:module:: vasca.tables

.. autodoc2-docstring:: vasca.tables
   :allowtitles:

Module Contents
---------------

Classes
~~~~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`TableCollection <vasca.tables.TableCollection>`
     - .. autodoc2-docstring:: vasca.tables.TableCollection
          :summary:

Data
~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`dimless <vasca.tables.dimless>`
     - .. autodoc2-docstring:: vasca.tables.dimless
          :summary:
   * - :py:obj:`FILE_DIR <vasca.tables.FILE_DIR>`
     - .. autodoc2-docstring:: vasca.tables.FILE_DIR
          :summary:
   * - :py:obj:`ROOT_DIR <vasca.tables.ROOT_DIR>`
     - .. autodoc2-docstring:: vasca.tables.ROOT_DIR
          :summary:

API
~~~

.. py:data:: dimless
   :canonical: vasca.tables.dimless
   :value: None

   .. autodoc2-docstring:: vasca.tables.dimless

.. py:data:: FILE_DIR
   :canonical: vasca.tables.FILE_DIR
   :value: 'dirname(...)'

   .. autodoc2-docstring:: vasca.tables.FILE_DIR

.. py:data:: ROOT_DIR
   :canonical: vasca.tables.ROOT_DIR
   :value: None

   .. autodoc2-docstring:: vasca.tables.ROOT_DIR

.. py:class:: TableCollection()
   :canonical: vasca.tables.TableCollection

   Bases: :py:obj:`object`

   .. autodoc2-docstring:: vasca.tables.TableCollection

   .. rubric:: Initialization

   .. autodoc2-docstring:: vasca.tables.TableCollection.__init__

   .. py:method:: table_from_template(dd_data, template_name)
      :canonical: vasca.tables.TableCollection.table_from_template
      :staticmethod:

      .. autodoc2-docstring:: vasca.tables.TableCollection.table_from_template

   .. py:method:: remove_tables(ll_table_names)
      :canonical: vasca.tables.TableCollection.remove_tables

      .. autodoc2-docstring:: vasca.tables.TableCollection.remove_tables

   .. py:method:: add_table(data, template_name)
      :canonical: vasca.tables.TableCollection.add_table

      .. autodoc2-docstring:: vasca.tables.TableCollection.add_table

   .. py:method:: remove_unselected(table_name)
      :canonical: vasca.tables.TableCollection.remove_unselected

      .. autodoc2-docstring:: vasca.tables.TableCollection.remove_unselected

   .. py:method:: write_to_fits(file_name='tables.fits', overwrite=True, fits_verify='fix')
      :canonical: vasca.tables.TableCollection.write_to_fits

      .. autodoc2-docstring:: vasca.tables.TableCollection.write_to_fits

   .. py:method:: load_from_fits(file_name)
      :canonical: vasca.tables.TableCollection.load_from_fits

      .. autodoc2-docstring:: vasca.tables.TableCollection.load_from_fits

   .. py:method:: info()
      :canonical: vasca.tables.TableCollection.info

      .. autodoc2-docstring:: vasca.tables.TableCollection.info

   .. py:method:: __str__()
      :canonical: vasca.tables.TableCollection.__str__

      .. autodoc2-docstring:: vasca.tables.TableCollection.__str__

   .. py:method:: select_from_config(dd_selections)
      :canonical: vasca.tables.TableCollection.select_from_config

      .. autodoc2-docstring:: vasca.tables.TableCollection.select_from_config

   .. py:method:: select_rows(selections, remove_unselected=False)
      :canonical: vasca.tables.TableCollection.select_rows

      .. autodoc2-docstring:: vasca.tables.TableCollection.select_rows

   .. py:method:: get_light_curve(fd_src_ids=None, rg_src_ids=None, flux_var='flux')
      :canonical: vasca.tables.TableCollection.get_light_curve

      .. autodoc2-docstring:: vasca.tables.TableCollection.get_light_curve

   .. py:method:: cluster_meanshift(**ms_kw)
      :canonical: vasca.tables.TableCollection.cluster_meanshift

      .. autodoc2-docstring:: vasca.tables.TableCollection.cluster_meanshift

   .. py:method:: set_src_stats(src_id_name='fd_src_id')
      :canonical: vasca.tables.TableCollection.set_src_stats

      .. autodoc2-docstring:: vasca.tables.TableCollection.set_src_stats

   .. py:method:: set_hardness_ratio(obs_filter_id1=1, obs_filter_id2=2)
      :canonical: vasca.tables.TableCollection.set_hardness_ratio

      .. autodoc2-docstring:: vasca.tables.TableCollection.set_hardness_ratio

   .. py:method:: add_column(table_name, col_name, col_data=None)
      :canonical: vasca.tables.TableCollection.add_column

      .. autodoc2-docstring:: vasca.tables.TableCollection.add_column

   .. py:method:: copy_table_columns(tab_name_to, tab_name_from, copy_vars, match_var='rg_src_id', select_matched=False)
      :canonical: vasca.tables.TableCollection.copy_table_columns

      .. autodoc2-docstring:: vasca.tables.TableCollection.copy_table_columns

   .. py:method:: cross_match(dist_max=1.5 * uu.arcsec, dist_s2n_max=3, cat_table_name='tt_coadd_sources', cat_id_name='coadd_src_id', cat_name='coadd', src_table_name='tt_sources')
      :canonical: vasca.tables.TableCollection.cross_match

      .. autodoc2-docstring:: vasca.tables.TableCollection.cross_match
