:py:mod:`vasca.region`
======================

.. py:module:: vasca.region

.. autodoc2-docstring:: vasca.region
   :allowtitles:

Module Contents
---------------

Classes
~~~~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`Region <vasca.region.Region>`
     - .. autodoc2-docstring:: vasca.region.Region
          :summary:

API
~~~

.. py:class:: Region()
   :canonical: vasca.region.Region

   Bases: :py:obj:`vasca.tables.TableCollection`

   .. autodoc2-docstring:: vasca.region.Region

   .. rubric:: Initialization

   .. autodoc2-docstring:: vasca.region.Region.__init__

   .. py:method:: load_from_config(vasca_cfg)
      :canonical: vasca.region.Region.load_from_config
      :classmethod:

      .. autodoc2-docstring:: vasca.region.Region.load_from_config

   .. py:method:: add_table_from_fields(table_name, only_selected=False)
      :canonical: vasca.region.Region.add_table_from_fields

      .. autodoc2-docstring:: vasca.region.Region.add_table_from_fields

   .. py:method:: add_coverage_hp(nside=4096, coord_sys='galactic')
      :canonical: vasca.region.Region.add_coverage_hp

      .. autodoc2-docstring:: vasca.region.Region.add_coverage_hp

   .. py:method:: load_from_fits(file_name, load_fields=False)
      :canonical: vasca.region.Region.load_from_fits

      .. autodoc2-docstring:: vasca.region.Region.load_from_fits

   .. py:method:: get_src_from_id(rg_src_id, load_from_file=True, write_to_file=True, add_sed=True, add_gphoton=True, add_spectrum=True)
      :canonical: vasca.region.Region.get_src_from_id

      .. autodoc2-docstring:: vasca.region.Region.get_src_from_id

   .. py:method:: set_src_id_info()
      :canonical: vasca.region.Region.set_src_id_info

      .. autodoc2-docstring:: vasca.region.Region.set_src_id_info

   .. py:method:: get_src_from_sky_pos(coordx, coordy, frame='icrs')
      :canonical: vasca.region.Region.get_src_from_sky_pos

      .. autodoc2-docstring:: vasca.region.Region.get_src_from_sky_pos

   .. py:method:: get_field(field_id=None, rg_fd_id=None, load_method='FITS', add_field=False, mast_products='TABLES', field_kwargs=dict())
      :canonical: vasca.region.Region.get_field

      .. autodoc2-docstring:: vasca.region.Region.get_field

   .. py:method:: cross_match_cds(query_radius=1.5 * uu.arcsec, query_table='I/355/gaiadr3', vizier_columns=['*', 'PQSO', 'PGal', 'PSS', 'RPlx', 'VarFlag', 'o_Gmag', 'RFRP', 'RFBP', 'AG', 'E(BP-RP)'], overwrite=False)
      :canonical: vasca.region.Region.cross_match_cds

      .. autodoc2-docstring:: vasca.region.Region.cross_match_cds

   .. py:method:: add_simbad_otype_info()
      :canonical: vasca.region.Region.add_simbad_otype_info

      .. autodoc2-docstring:: vasca.region.Region.add_simbad_otype_info

   .. py:method:: synch_src_sel(remove_unselected=False)
      :canonical: vasca.region.Region.synch_src_sel

      .. autodoc2-docstring:: vasca.region.Region.synch_src_sel

   .. py:method:: get_region_catalog()
      :canonical: vasca.region.Region.get_region_catalog

      .. autodoc2-docstring:: vasca.region.Region.get_region_catalog

   .. py:method:: set_LombScargle(obs_filters=['NUV', 'FUV'], nbins_min=20)
      :canonical: vasca.region.Region.set_LombScargle

      .. autodoc2-docstring:: vasca.region.Region.set_LombScargle

   .. py:method:: redo_src_selection(cfg_file_name='./vasca_cfg.yaml')
      :canonical: vasca.region.Region.redo_src_selection

      .. autodoc2-docstring:: vasca.region.Region.redo_src_selection
