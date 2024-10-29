:py:mod:`vasca.source`
======================

.. py:module:: vasca.source

.. autodoc2-docstring:: vasca.source
   :allowtitles:

Module Contents
---------------

Classes
~~~~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`Source <vasca.source.Source>`
     - .. autodoc2-docstring:: vasca.source.Source
          :summary:

Data
~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`rm <vasca.source.rm>`
     - .. autodoc2-docstring:: vasca.source.rm
          :summary:

API
~~~

.. py:data:: rm
   :canonical: vasca.source.rm
   :value: 'ResourceManager(...)'

   .. autodoc2-docstring:: vasca.source.rm

.. py:class:: Source()
   :canonical: vasca.source.Source

   Bases: :py:obj:`vasca.tables.TableCollection`

   .. autodoc2-docstring:: vasca.source.Source

   .. rubric:: Initialization

   .. autodoc2-docstring:: vasca.source.Source.__init__

   .. py:method:: add_vizier_SED(vizier_radius: astropy.units.Quantity = 1 * uu.arcsec) -> None
      :canonical: vasca.source.Source.add_vizier_SED

      .. autodoc2-docstring:: vasca.source.Source.add_vizier_SED

   .. py:method:: add_gphoton_lc(s2n_min: float = 3.0, tbin: int = -1) -> None
      :canonical: vasca.source.Source.add_gphoton_lc

      .. autodoc2-docstring:: vasca.source.Source.add_gphoton_lc

   .. py:method:: add_spectrum(search_radius: astropy.units.Quantity = 2 * uu.arcsec) -> None
      :canonical: vasca.source.Source.add_spectrum

      .. autodoc2-docstring:: vasca.source.Source.add_spectrum
