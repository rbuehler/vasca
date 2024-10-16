:py:mod:`vasca.field`
=====================

.. py:module:: vasca.field

.. autodoc2-docstring:: vasca.field
   :allowtitles:

Module Contents
---------------

Classes
~~~~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`BaseField <vasca.field.BaseField>`
     - .. autodoc2-docstring:: vasca.field.BaseField
          :summary:
   * - :py:obj:`GALEXField <vasca.field.GALEXField>`
     - .. autodoc2-docstring:: vasca.field.GALEXField
          :summary:
   * - :py:obj:`GALEXDSField <vasca.field.GALEXDSField>`
     - .. autodoc2-docstring:: vasca.field.GALEXDSField
          :summary:

Data
~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`FILE_DIR <vasca.field.FILE_DIR>`
     - .. autodoc2-docstring:: vasca.field.FILE_DIR
          :summary:
   * - :py:obj:`ROOT_DIR <vasca.field.ROOT_DIR>`
     - .. autodoc2-docstring:: vasca.field.ROOT_DIR
          :summary:
   * - :py:obj:`dimless <vasca.field.dimless>`
     - .. autodoc2-docstring:: vasca.field.dimless
          :summary:

API
~~~

.. py:data:: FILE_DIR
   :canonical: vasca.field.FILE_DIR
   :value: 'dirname(...)'

   .. autodoc2-docstring:: vasca.field.FILE_DIR

.. py:data:: ROOT_DIR
   :canonical: vasca.field.ROOT_DIR
   :value: None

   .. autodoc2-docstring:: vasca.field.ROOT_DIR

.. py:data:: dimless
   :canonical: vasca.field.dimless
   :value: None

   .. autodoc2-docstring:: vasca.field.dimless

.. py:class:: BaseField()
   :canonical: vasca.field.BaseField

   Bases: :py:obj:`vasca.tables.TableCollection`

   .. autodoc2-docstring:: vasca.field.BaseField

   .. rubric:: Initialization

   .. autodoc2-docstring:: vasca.field.BaseField.__init__

   .. py:method:: load_sky_map(file_name, img_attr='ref_img')
      :canonical: vasca.field.BaseField.load_sky_map

      .. autodoc2-docstring:: vasca.field.BaseField.load_sky_map

   .. py:method:: get_upper_limits()
      :canonical: vasca.field.BaseField.get_upper_limits

      .. autodoc2-docstring:: vasca.field.BaseField.get_upper_limits

   .. py:method:: set_light_curve(add_upper_limits=True)
      :canonical: vasca.field.BaseField.set_light_curve

      .. autodoc2-docstring:: vasca.field.BaseField.set_light_curve

   .. py:method:: remove_double_visit_detections()
      :canonical: vasca.field.BaseField.remove_double_visit_detections

      .. autodoc2-docstring:: vasca.field.BaseField.remove_double_visit_detections

   .. py:method:: set_field_attr(dd_names=None)
      :canonical: vasca.field.BaseField.set_field_attr

      .. autodoc2-docstring:: vasca.field.BaseField.set_field_attr

   .. py:method:: get_field_par(par_name, table_name)
      :canonical: vasca.field.BaseField.get_field_par

      .. autodoc2-docstring:: vasca.field.BaseField.get_field_par

   .. py:method:: get_sky_region()
      :canonical: vasca.field.BaseField.get_sky_region

      .. autodoc2-docstring:: vasca.field.BaseField.get_sky_region

   .. py:method:: load_from_fits(file_name='tables.fits')
      :canonical: vasca.field.BaseField.load_from_fits

      .. autodoc2-docstring:: vasca.field.BaseField.load_from_fits

.. py:class:: GALEXField(obs_id, obs_filter=None, data_path=None, visits_data_path=None)
   :canonical: vasca.field.GALEXField

   Bases: :py:obj:`vasca.field.BaseField`

   .. autodoc2-docstring:: vasca.field.GALEXField

   .. rubric:: Initialization

   .. autodoc2-docstring:: vasca.field.GALEXField.__init__

   .. py:method:: from_VASCA(obs_id, obs_filter='NUV', fits_path=None, **kwargs)
      :canonical: vasca.field.GALEXField.from_VASCA
      :classmethod:

      .. autodoc2-docstring:: vasca.field.GALEXField.from_VASCA

   .. py:method:: from_MAST(obs_id, obs_filter='NUV', refresh=False, load_products='TABLES', write=True, **kwargs)
      :canonical: vasca.field.GALEXField.from_MAST
      :classmethod:

      .. autodoc2-docstring:: vasca.field.GALEXField.from_MAST

   .. py:method:: load(gfield_id, obs_filter='NUV', method='MAST_LOCAL', load_products='TABLES', **field_kwargs)
      :canonical: vasca.field.GALEXField.load
      :staticmethod:

      .. autodoc2-docstring:: vasca.field.GALEXField.load

   .. py:method:: get_visit_upper_limits(tt_visits)
      :canonical: vasca.field.GALEXField.get_visit_upper_limits
      :staticmethod:

      .. autodoc2-docstring:: vasca.field.GALEXField.get_visit_upper_limits

   .. py:method:: _load_galex_field_info(obs_id, obs_filter, col_names=None, refresh=False)
      :canonical: vasca.field.GALEXField._load_galex_field_info

      .. autodoc2-docstring:: vasca.field.GALEXField._load_galex_field_info

   .. py:method:: _load_galex_visits_info(obs_id, obs_filter, col_names=None)
      :canonical: vasca.field.GALEXField._load_galex_visits_info

      .. autodoc2-docstring:: vasca.field.GALEXField._load_galex_visits_info

   .. py:method:: _load_galex_archive_products(obs_id, obs_filter, col_names=None, dd_products=None, ref_maps_only=False, refresh=False)
      :canonical: vasca.field.GALEXField._load_galex_archive_products

      .. autodoc2-docstring:: vasca.field.GALEXField._load_galex_archive_products

.. py:class:: GALEXDSField(field_name, data_path=None, visits_data_path=None, **kwargs)
   :canonical: vasca.field.GALEXDSField

   Bases: :py:obj:`vasca.field.BaseField`

   .. autodoc2-docstring:: vasca.field.GALEXDSField

   .. rubric:: Initialization

   .. autodoc2-docstring:: vasca.field.GALEXDSField.__init__

   .. py:method:: name_id_map(a)
      :canonical: vasca.field.GALEXDSField.name_id_map

      .. autodoc2-docstring:: vasca.field.GALEXDSField.name_id_map

   .. py:method:: _load_info(field_name)
      :canonical: vasca.field.GALEXDSField._load_info

      .. autodoc2-docstring:: vasca.field.GALEXDSField._load_info

   .. py:method:: _load_data_products(field_name, ref_maps_only=False)
      :canonical: vasca.field.GALEXDSField._load_data_products

      .. autodoc2-docstring:: vasca.field.GALEXDSField._load_data_products

   .. py:method:: from_VASCA(field_name, fits_path=None, **kwargs)
      :canonical: vasca.field.GALEXDSField.from_VASCA
      :classmethod:

      .. autodoc2-docstring:: vasca.field.GALEXDSField.from_VASCA

   .. py:method:: from_MAST(field_name, load_products='TABLES', write=True, **kwargs)
      :canonical: vasca.field.GALEXDSField.from_MAST
      :classmethod:

      .. autodoc2-docstring:: vasca.field.GALEXDSField.from_MAST

   .. py:method:: load(field_name, method='MAST_LOCAL', load_products='TABLES', **field_kwargs)
      :canonical: vasca.field.GALEXDSField.load
      :staticmethod:

      .. autodoc2-docstring:: vasca.field.GALEXDSField.load
