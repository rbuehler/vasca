:py:mod:`vasca.visualization`
=============================

.. py:module:: vasca.visualization

.. autodoc2-docstring:: vasca.visualization
   :allowtitles:

Module Contents
---------------

Functions
~~~~~~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`plot_sky_sources <vasca.visualization.plot_sky_sources>`
     - .. autodoc2-docstring:: vasca.visualization.plot_sky_sources
          :summary:
   * - :py:obj:`plot_field_sky_map <vasca.visualization.plot_field_sky_map>`
     - .. autodoc2-docstring:: vasca.visualization.plot_field_sky_map
          :summary:
   * - :py:obj:`plot_region_sky_gnomeview <vasca.visualization.plot_region_sky_gnomeview>`
     - .. autodoc2-docstring:: vasca.visualization.plot_region_sky_gnomeview
          :summary:
   * - :py:obj:`plot_region_sky_mollview <vasca.visualization.plot_region_sky_mollview>`
     - .. autodoc2-docstring:: vasca.visualization.plot_region_sky_mollview
          :summary:
   * - :py:obj:`plot_table_hist <vasca.visualization.plot_table_hist>`
     - .. autodoc2-docstring:: vasca.visualization.plot_table_hist
          :summary:
   * - :py:obj:`scatter_hist <vasca.visualization.scatter_hist>`
     - .. autodoc2-docstring:: vasca.visualization.scatter_hist
          :summary:
   * - :py:obj:`plot_table_scatter <vasca.visualization.plot_table_scatter>`
     - .. autodoc2-docstring:: vasca.visualization.plot_table_scatter
          :summary:
   * - :py:obj:`plot_pipe_diagnostic <vasca.visualization.plot_pipe_diagnostic>`
     - .. autodoc2-docstring:: vasca.visualization.plot_pipe_diagnostic
          :summary:
   * - :py:obj:`plot_light_curves <vasca.visualization.plot_light_curves>`
     - .. autodoc2-docstring:: vasca.visualization.plot_light_curves
          :summary:
   * - :py:obj:`plot_light_curve <vasca.visualization.plot_light_curve>`
     - .. autodoc2-docstring:: vasca.visualization.plot_light_curve
          :summary:
   * - :py:obj:`plot_lombscargle <vasca.visualization.plot_lombscargle>`
     - .. autodoc2-docstring:: vasca.visualization.plot_lombscargle
          :summary:
   * - :py:obj:`plot_sed <vasca.visualization.plot_sed>`
     - .. autodoc2-docstring:: vasca.visualization.plot_sed
          :summary:

API
~~~

.. py:function:: plot_sky_sources(tt_src, tt_det=None, only_selected=True, ax=None, src_id='rg_src_id', sky_region_wcs=None, draw_labels=True, src_kwargs=None, det_kwargs=None)
   :canonical: vasca.visualization.plot_sky_sources

   .. autodoc2-docstring:: vasca.visualization.plot_sky_sources

.. py:function:: plot_field_sky_map(field, fig=None, ax=None, img_idx=-1, sky_region=None, **img_kwargs)
   :canonical: vasca.visualization.plot_field_sky_map

   .. autodoc2-docstring:: vasca.visualization.plot_field_sky_map

.. py:function:: plot_region_sky_gnomeview(region, ra, dec, sel_srcs=True, gw_kwargs=None, ps_kwargs=None)
   :canonical: vasca.visualization.plot_region_sky_gnomeview

   .. autodoc2-docstring:: vasca.visualization.plot_region_sky_gnomeview

.. py:function:: plot_region_sky_mollview(region, var='nr_vis', mw_kwargs=None)
   :canonical: vasca.visualization.plot_region_sky_mollview

   .. autodoc2-docstring:: vasca.visualization.plot_region_sky_mollview

.. py:function:: plot_table_hist(tt, var, ax=None, logx=False, obs_filter_id=None, **hist_kwargs)
   :canonical: vasca.visualization.plot_table_hist

   .. autodoc2-docstring:: vasca.visualization.plot_table_hist

.. py:function:: scatter_hist(x, y, ax, ax_histx, ax_histy)
   :canonical: vasca.visualization.scatter_hist

   .. autodoc2-docstring:: vasca.visualization.scatter_hist

.. py:function:: plot_table_scatter(tt, varx, vary, ax=None, xlim=None, ylim=None, invert_xaxis=None, invert_yaxis=None, xscale='linear', yscale='linear', obs_filter_id=None, grp_var='sel', grp_vals=None, add_projection=False, **scatter_kwargs)
   :canonical: vasca.visualization.plot_table_scatter

   .. autodoc2-docstring:: vasca.visualization.plot_table_scatter

.. py:function:: plot_pipe_diagnostic(tc, table_name, plot_type, fig_size=(12, 8), obs_filter_id=None)
   :canonical: vasca.visualization.plot_pipe_diagnostic

   .. autodoc2-docstring:: vasca.visualization.plot_pipe_diagnostic

.. py:function:: plot_light_curves(tc, fd_src_ids=None, rg_src_ids=None, fig=None, ax=None, ylim=None, plot_upper_limits=True, flux_var='flux', **errorbar_kwargs)
   :canonical: vasca.visualization.plot_light_curves

   .. autodoc2-docstring:: vasca.visualization.plot_light_curves

.. py:function:: plot_light_curve(tc_src, fig=None, ax=None, show_gphoton=True, add_axes=True, **errorbar_kwargs)
   :canonical: vasca.visualization.plot_light_curve

   .. autodoc2-docstring:: vasca.visualization.plot_light_curve

.. py:function:: plot_lombscargle(tt_lc, fig=None, ax=None, ax_phase=None, ax_lc=None, obs_filter='NUV', nbins_min=10, logy=False, freq_range=[0.03, 2] / uu.d, plot_dtbins=True)
   :canonical: vasca.visualization.plot_lombscargle

   .. autodoc2-docstring:: vasca.visualization.plot_lombscargle

.. py:function:: plot_sed(tc_src, fig=None, ax=None, plot_spec_lines=False, plot_spec=False, **errorbar_kwargs)
   :canonical: vasca.visualization.plot_sed

   .. autodoc2-docstring:: vasca.visualization.plot_sed
