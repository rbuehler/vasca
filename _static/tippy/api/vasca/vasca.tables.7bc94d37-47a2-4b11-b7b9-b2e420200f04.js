selector_to_html = {"a[href=\"#vasca.tables.TableCollection.table_from_template\"]": "<dt class=\"sig sig-object py\" id=\"vasca.tables.TableCollection.table_from_template\">\n<em class=\"property\"><span class=\"pre\">static</span><span class=\"w\"> </span></em><span class=\"sig-name descname\"><span class=\"pre\">table_from_template</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">dd_data</span></span><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"n\"><a class=\"reference external\" href=\"https://docs.python.org/3/library/stdtypes.html#dict\" title=\"(in Python v3.13)\"><span class=\"pre\">dict</span></a></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">template_name</span></span><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"n\"><a class=\"reference external\" href=\"https://docs.python.org/3/library/stdtypes.html#str\" title=\"(in Python v3.13)\"><span class=\"pre\">str</span></a></span></em><span class=\"sig-paren\">)</span> <span class=\"sig-return\"><span class=\"sig-return-icon\">\u2192</span> <span class=\"sig-return-typehint\"><a class=\"reference external\" href=\"https://docs.astropy.org/en/stable/api/astropy.table.Table.html#astropy.table.Table\" title=\"(in Astropy v6.1)\"><span class=\"pre\">astropy.table.Table</span></a></span></span><a class=\"reference internal\" href=\"../../_modules/vasca/tables.html#TableCollection.table_from_template\"><span class=\"viewcode-link\"><span class=\"pre\">[source]</span></span></a></dt><dd><p>Creates a new astropy table.</p></dd>", "a[href=\"#vasca.tables.TableCollection.set_hardness_ratio\"]": "<dt class=\"sig sig-object py\" id=\"vasca.tables.TableCollection.set_hardness_ratio\">\n<span class=\"sig-name descname\"><span class=\"pre\">set_hardness_ratio</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">obs_filter_id1</span></span><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"n\"><a class=\"reference external\" href=\"https://docs.python.org/3/library/functions.html#int\" title=\"(in Python v3.13)\"><span class=\"pre\">int</span></a></span><span class=\"w\"> </span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"w\"> </span><span class=\"default_value\"><span class=\"pre\">1</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">obs_filter_id2</span></span><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"n\"><a class=\"reference external\" href=\"https://docs.python.org/3/library/functions.html#int\" title=\"(in Python v3.13)\"><span class=\"pre\">int</span></a></span><span class=\"w\"> </span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"w\"> </span><span class=\"default_value\"><span class=\"pre\">2</span></span></em><span class=\"sig-paren\">)</span> <span class=\"sig-return\"><span class=\"sig-return-icon\">\u2192</span> <span class=\"sig-return-typehint\"><a class=\"reference external\" href=\"https://docs.python.org/3/library/constants.html#None\" title=\"(in Python v3.13)\"><span class=\"pre\">None</span></a></span></span><a class=\"reference internal\" href=\"../../_modules/vasca/tables.html#TableCollection.set_hardness_ratio\"><span class=\"viewcode-link\"><span class=\"pre\">[source]</span></span></a></dt><dd><p>Calculated hardness ratio from detections flux(filter_2)/ flux(filter_1).\nOnly simultaneous detections are considered</p></dd>", "a[href=\"#module-vasca.tables\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\"><a class=\"reference internal\" href=\"#module-vasca.tables\" title=\"vasca.tables\"><code class=\"xref py py-mod docutils literal notranslate\"><span class=\"pre\">vasca.tables</span></code></a><a class=\"headerlink\" href=\"#module-vasca.tables\" title=\"Link to this heading\">\u00b6</a></h1><h2>Module Contents<a class=\"headerlink\" href=\"#module-contents\" title=\"Link to this heading\">\u00b6</a></h2><h3>Classes<a class=\"headerlink\" href=\"#classes\" title=\"Link to this heading\">\u00b6</a></h3>", "a[href=\"#vasca.tables.TableCollection.cross_match\"]": "<dt class=\"sig sig-object py\" id=\"vasca.tables.TableCollection.cross_match\">\n<span class=\"sig-name descname\"><span class=\"pre\">cross_match</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">dist_max</span></span><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"n\"><a class=\"reference external\" href=\"https://docs.astropy.org/en/stable/api/astropy.units.Quantity.html#astropy.units.Quantity\" title=\"(in Astropy v6.1)\"><span class=\"pre\">astropy.units.Quantity</span></a></span><span class=\"w\"> </span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"w\"> </span><span class=\"default_value\"><span class=\"pre\">1.5</span> <span class=\"pre\">*</span> <span class=\"pre\">uu.arcsec</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">dist_s2n_max</span></span><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"n\"><a class=\"reference external\" href=\"https://docs.python.org/3/library/functions.html#float\" title=\"(in Python v3.13)\"><span class=\"pre\">float</span></a></span><span class=\"w\"> </span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"w\"> </span><span class=\"default_value\"><span class=\"pre\">3</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">cat_table_name</span></span><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"n\"><a class=\"reference external\" href=\"https://docs.python.org/3/library/stdtypes.html#str\" title=\"(in Python v3.13)\"><span class=\"pre\">str</span></a></span><span class=\"w\"> </span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"w\"> </span><span class=\"default_value\"><span class=\"pre\">'tt_coadd_sources'</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">cat_id_name</span></span><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"n\"><a class=\"reference external\" href=\"https://docs.python.org/3/library/stdtypes.html#str\" title=\"(in Python v3.13)\"><span class=\"pre\">str</span></a></span><span class=\"w\"> </span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"w\"> </span><span class=\"default_value\"><span class=\"pre\">'coadd_src_id'</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">cat_name</span></span><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"n\"><a class=\"reference external\" href=\"https://docs.python.org/3/library/stdtypes.html#str\" title=\"(in Python v3.13)\"><span class=\"pre\">str</span></a></span><span class=\"w\"> </span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"w\"> </span><span class=\"default_value\"><span class=\"pre\">'coadd'</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">src_table_name</span></span><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"n\"><a class=\"reference external\" href=\"https://docs.python.org/3/library/stdtypes.html#str\" title=\"(in Python v3.13)\"><span class=\"pre\">str</span></a></span><span class=\"w\"> </span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"w\"> </span><span class=\"default_value\"><span class=\"pre\">'tt_sources'</span></span></em><span class=\"sig-paren\">)</span> <span class=\"sig-return\"><span class=\"sig-return-icon\">\u2192</span> <span class=\"sig-return-typehint\"><a class=\"reference external\" href=\"https://docs.python.org/3/library/constants.html#None\" title=\"(in Python v3.13)\"><span class=\"pre\">None</span></a></span></span><a class=\"reference internal\" href=\"../../_modules/vasca/tables.html#TableCollection.cross_match\"><span class=\"viewcode-link\"><span class=\"pre\">[source]</span></span></a></dt><dd><p>Cross match sources to a catalog by position. Typically this is the coadd\ncatalog.</p></dd>", "a[href=\"#api\"]": "<h3 class=\"tippy-header\" style=\"margin-top: 0;\">API<a class=\"headerlink\" href=\"#api\" title=\"Link to this heading\">\u00b6</a></h3>", "a[href=\"#vasca.tables.TableCollection.remove_tables\"]": "<dt class=\"sig sig-object py\" id=\"vasca.tables.TableCollection.remove_tables\">\n<span class=\"sig-name descname\"><span class=\"pre\">remove_tables</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">ll_table_names</span></span><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"n\"><a class=\"reference external\" href=\"https://docs.python.org/3/library/stdtypes.html#list\" title=\"(in Python v3.13)\"><span class=\"pre\">list</span></a><span class=\"p\"><span class=\"pre\">[</span></span><a class=\"reference external\" href=\"https://docs.python.org/3/library/stdtypes.html#str\" title=\"(in Python v3.13)\"><span class=\"pre\">str</span></a><span class=\"p\"><span class=\"pre\">]</span></span></span></em><span class=\"sig-paren\">)</span> <span class=\"sig-return\"><span class=\"sig-return-icon\">\u2192</span> <span class=\"sig-return-typehint\"><a class=\"reference external\" href=\"https://docs.python.org/3/library/constants.html#None\" title=\"(in Python v3.13)\"><span class=\"pre\">None</span></a></span></span><a class=\"reference internal\" href=\"../../_modules/vasca/tables.html#TableCollection.remove_tables\"><span class=\"viewcode-link\"><span class=\"pre\">[source]</span></span></a></dt><dd><p>Removes table from collection</p></dd>", "a[href=\"#vasca.tables.TableCollection.select_from_config\"]": "<dt class=\"sig sig-object py\" id=\"vasca.tables.TableCollection.select_from_config\">\n<span class=\"sig-name descname\"><span class=\"pre\">select_from_config</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">dd_selections</span></span><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"n\"><a class=\"reference external\" href=\"https://docs.python.org/3/library/stdtypes.html#dict\" title=\"(in Python v3.13)\"><span class=\"pre\">dict</span></a></span></em><span class=\"sig-paren\">)</span> <span class=\"sig-return\"><span class=\"sig-return-icon\">\u2192</span> <span class=\"sig-return-typehint\"><a class=\"reference external\" href=\"https://docs.python.org/3/library/constants.html#None\" title=\"(in Python v3.13)\"><span class=\"pre\">None</span></a></span></span><a class=\"reference internal\" href=\"../../_modules/vasca/tables.html#TableCollection.select_from_config\"><span class=\"viewcode-link\"><span class=\"pre\">[source]</span></span></a></dt><dd><p>Apply multiple selections at once.</p></dd>", "a[href=\"#vasca.tables.TableCollection.add_table\"]": "<dt class=\"sig sig-object py\" id=\"vasca.tables.TableCollection.add_table\">\n<span class=\"sig-name descname\"><span class=\"pre\">add_table</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">data</span></span><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"n\"><a class=\"reference external\" href=\"https://docs.python.org/3/library/stdtypes.html#list\" title=\"(in Python v3.13)\"><span class=\"pre\">list</span></a><span class=\"w\"> </span><span class=\"p\"><span class=\"pre\">|</span></span><span class=\"w\"> </span><a class=\"reference external\" href=\"https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray\" title=\"(in NumPy v2.1)\"><span class=\"pre\">numpy.ndarray</span></a></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">template_name</span></span><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"n\"><a class=\"reference external\" href=\"https://docs.python.org/3/library/stdtypes.html#str\" title=\"(in Python v3.13)\"><span class=\"pre\">str</span></a></span></em><span class=\"sig-paren\">)</span> <span class=\"sig-return\"><span class=\"sig-return-icon\">\u2192</span> <span class=\"sig-return-typehint\"><a class=\"reference external\" href=\"https://docs.python.org/3/library/constants.html#None\" title=\"(in Python v3.13)\"><span class=\"pre\">None</span></a></span></span><a class=\"reference internal\" href=\"../../_modules/vasca/tables.html#TableCollection.add_table\"><span class=\"viewcode-link\"><span class=\"pre\">[source]</span></span></a></dt><dd><p>Add a VASCA table to the colection</p></dd>", "a[href=\"#vasca.tables.dimless\"]": "<dt class=\"sig sig-object py\" id=\"vasca.tables.dimless\">\n<span class=\"sig-prename descclassname\"><span class=\"pre\">vasca.tables.</span></span><span class=\"sig-name descname\"><span class=\"pre\">dimless</span></span><em class=\"property\"><span class=\"w\"> </span><span class=\"p\"><span class=\"pre\">=</span></span><span class=\"w\"> </span><span class=\"pre\">None</span></em></dt><dd></dd>", "a[href=\"vasca.tables_dict.html#module-vasca.tables_dict\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\"><a class=\"reference internal\" href=\"#module-vasca.tables_dict\" title=\"vasca.tables_dict\"><code class=\"xref py py-mod docutils literal notranslate\"><span class=\"pre\">vasca.tables_dict</span></code></a><a class=\"headerlink\" href=\"#module-vasca.tables_dict\" title=\"Link to this heading\">\u00b6</a></h1><p>Defines dictionary for the tables used by vasca.tables.TableCollection</p>", "a[href=\"#classes\"]": "<h3 class=\"tippy-header\" style=\"margin-top: 0;\">Classes<a class=\"headerlink\" href=\"#classes\" title=\"Link to this heading\">\u00b6</a></h3>", "a[href=\"#vasca.tables.TableCollection.set_src_stats\"]": "<dt class=\"sig sig-object py\" id=\"vasca.tables.TableCollection.set_src_stats\">\n<span class=\"sig-name descname\"><span class=\"pre\">set_src_stats</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">src_id_name</span></span><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"n\"><a class=\"reference external\" href=\"https://docs.python.org/3/library/stdtypes.html#str\" title=\"(in Python v3.13)\"><span class=\"pre\">str</span></a></span><span class=\"w\"> </span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"w\"> </span><span class=\"default_value\"><span class=\"pre\">'fd_src_id'</span></span></em><span class=\"sig-paren\">)</span> <span class=\"sig-return\"><span class=\"sig-return-icon\">\u2192</span> <span class=\"sig-return-typehint\"><a class=\"reference external\" href=\"https://docs.python.org/3/library/constants.html#None\" title=\"(in Python v3.13)\"><span class=\"pre\">None</span></a></span></span><a class=\"reference internal\" href=\"../../_modules/vasca/tables.html#TableCollection.set_src_stats\"><span class=\"viewcode-link\"><span class=\"pre\">[source]</span></span></a></dt><dd><p>Calculates source parameters from detections and stores them\nin the source table (tt_source).</p></dd>", "a[href=\"#vasca.tables.TableCollection.load_from_fits\"]": "<dt class=\"sig sig-object py\" id=\"vasca.tables.TableCollection.load_from_fits\">\n<span class=\"sig-name descname\"><span class=\"pre\">load_from_fits</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">file_name</span></span><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"n\"><a class=\"reference external\" href=\"https://docs.python.org/3/library/stdtypes.html#str\" title=\"(in Python v3.13)\"><span class=\"pre\">str</span></a></span></em><span class=\"sig-paren\">)</span> <span class=\"sig-return\"><span class=\"sig-return-icon\">\u2192</span> <span class=\"sig-return-typehint\"><a class=\"reference external\" href=\"https://docs.python.org/3/library/constants.html#None\" title=\"(in Python v3.13)\"><span class=\"pre\">None</span></a></span></span><a class=\"reference internal\" href=\"../../_modules/vasca/tables.html#TableCollection.load_from_fits\"><span class=\"viewcode-link\"><span class=\"pre\">[source]</span></span></a></dt><dd><p>Loads collection from a fits file</p></dd>", "a[href=\"#vasca.tables.TableCollection.write_to_fits\"]": "<dt class=\"sig sig-object py\" id=\"vasca.tables.TableCollection.write_to_fits\">\n<span class=\"sig-name descname\"><span class=\"pre\">write_to_fits</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">file_name</span></span><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"n\"><a class=\"reference external\" href=\"https://docs.python.org/3/library/stdtypes.html#str\" title=\"(in Python v3.13)\"><span class=\"pre\">str</span></a></span><span class=\"w\"> </span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"w\"> </span><span class=\"default_value\"><span class=\"pre\">'tables.fits'</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">overwrite</span></span><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"n\"><a class=\"reference external\" href=\"https://docs.python.org/3/library/functions.html#bool\" title=\"(in Python v3.13)\"><span class=\"pre\">bool</span></a></span><span class=\"w\"> </span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"w\"> </span><span class=\"default_value\"><span class=\"pre\">True</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">fits_verify</span></span><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"n\"><a class=\"reference external\" href=\"https://docs.python.org/3/library/stdtypes.html#str\" title=\"(in Python v3.13)\"><span class=\"pre\">str</span></a></span><span class=\"w\"> </span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"w\"> </span><span class=\"default_value\"><span class=\"pre\">'fix'</span></span></em><span class=\"sig-paren\">)</span> <span class=\"sig-return\"><span class=\"sig-return-icon\">\u2192</span> <span class=\"sig-return-typehint\"><a class=\"reference external\" href=\"https://docs.python.org/3/library/constants.html#None\" title=\"(in Python v3.13)\"><span class=\"pre\">None</span></a></span></span><a class=\"reference internal\" href=\"../../_modules/vasca/tables.html#TableCollection.write_to_fits\"><span class=\"viewcode-link\"><span class=\"pre\">[source]</span></span></a></dt><dd><p>Write tables and image of a collection to a fits file.</p></dd>", "a[href=\"#vasca.tables.ROOT_DIR\"]": "<dt class=\"sig sig-object py\" id=\"vasca.tables.ROOT_DIR\">\n<span class=\"sig-prename descclassname\"><span class=\"pre\">vasca.tables.</span></span><span class=\"sig-name descname\"><span class=\"pre\">ROOT_DIR</span></span><em class=\"property\"><span class=\"w\"> </span><span class=\"p\"><span class=\"pre\">=</span></span><span class=\"w\"> </span><span class=\"pre\">None</span></em></dt><dd></dd>", "a[href=\"#vasca.tables.TableCollection.cluster_meanshift\"]": "<dt class=\"sig sig-object py\" id=\"vasca.tables.TableCollection.cluster_meanshift\">\n<span class=\"sig-name descname\"><span class=\"pre\">cluster_meanshift</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"o\"><span class=\"pre\">**</span></span><span class=\"n\"><span class=\"pre\">ms_kw</span></span></em><span class=\"sig-paren\">)</span> <span class=\"sig-return\"><span class=\"sig-return-icon\">\u2192</span> <span class=\"sig-return-typehint\"><a class=\"reference external\" href=\"https://docs.python.org/3/library/functions.html#int\" title=\"(in Python v3.13)\"><span class=\"pre\">int</span></a></span></span><a class=\"reference internal\" href=\"../../_modules/vasca/tables.html#TableCollection.cluster_meanshift\"><span class=\"viewcode-link\"><span class=\"pre\">[source]</span></span></a></dt><dd><p>Apply <a class=\"reference external\" href=\"https://scikit-learn.org/stable/modules/generated/sklearn.cluster.MeanShift.html\">MeanShift</a> clustering algorithm using to derive sources. Runs only on\nselected detections or sources.</p></dd>", "a[href=\"#vasca.tables.TableCollection.get_light_curve\"]": "<dt class=\"sig sig-object py\" id=\"vasca.tables.TableCollection.get_light_curve\">\n<span class=\"sig-name descname\"><span class=\"pre\">get_light_curve</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">fd_src_ids</span></span><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"n\"><a class=\"reference external\" href=\"https://docs.python.org/3/library/functions.html#int\" title=\"(in Python v3.13)\"><span class=\"pre\">int</span></a><span class=\"w\"> </span><span class=\"p\"><span class=\"pre\">|</span></span><span class=\"w\"> </span><a class=\"reference external\" href=\"https://docs.python.org/3/library/stdtypes.html#list\" title=\"(in Python v3.13)\"><span class=\"pre\">list</span></a><span class=\"p\"><span class=\"pre\">[</span></span><a class=\"reference external\" href=\"https://docs.python.org/3/library/functions.html#int\" title=\"(in Python v3.13)\"><span class=\"pre\">int</span></a><span class=\"p\"><span class=\"pre\">]</span></span><span class=\"w\"> </span><span class=\"p\"><span class=\"pre\">|</span></span><span class=\"w\"> </span><a class=\"reference external\" href=\"https://docs.python.org/3/library/constants.html#None\" title=\"(in Python v3.13)\"><span class=\"pre\">None</span></a></span><span class=\"w\"> </span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"w\"> </span><span class=\"default_value\"><span class=\"pre\">None</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">rg_src_ids</span></span><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"n\"><a class=\"reference external\" href=\"https://docs.python.org/3/library/functions.html#int\" title=\"(in Python v3.13)\"><span class=\"pre\">int</span></a><span class=\"w\"> </span><span class=\"p\"><span class=\"pre\">|</span></span><span class=\"w\"> </span><a class=\"reference external\" href=\"https://docs.python.org/3/library/stdtypes.html#list\" title=\"(in Python v3.13)\"><span class=\"pre\">list</span></a><span class=\"p\"><span class=\"pre\">[</span></span><a class=\"reference external\" href=\"https://docs.python.org/3/library/functions.html#int\" title=\"(in Python v3.13)\"><span class=\"pre\">int</span></a><span class=\"p\"><span class=\"pre\">]</span></span><span class=\"w\"> </span><span class=\"p\"><span class=\"pre\">|</span></span><span class=\"w\"> </span><a class=\"reference external\" href=\"https://docs.python.org/3/library/constants.html#None\" title=\"(in Python v3.13)\"><span class=\"pre\">None</span></a></span><span class=\"w\"> </span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"w\"> </span><span class=\"default_value\"><span class=\"pre\">None</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">flux_var</span></span><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"n\"><a class=\"reference external\" href=\"https://docs.python.org/3/library/stdtypes.html#str\" title=\"(in Python v3.13)\"><span class=\"pre\">str</span></a></span><span class=\"w\"> </span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"w\"> </span><span class=\"default_value\"><span class=\"pre\">'flux'</span></span></em><span class=\"sig-paren\">)</span> <span class=\"sig-return\"><span class=\"sig-return-icon\">\u2192</span> <span class=\"sig-return-typehint\"><a class=\"reference external\" href=\"https://docs.python.org/3/library/stdtypes.html#dict\" title=\"(in Python v3.13)\"><span class=\"pre\">dict</span></a></span></span><a class=\"reference internal\" href=\"../../_modules/vasca/tables.html#TableCollection.get_light_curve\"><span class=\"viewcode-link\"><span class=\"pre\">[source]</span></span></a></dt><dd><p>Get a light curves for one or list of sources, for regions or fields.</p></dd>", "a[href=\"#module-contents\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Module Contents<a class=\"headerlink\" href=\"#module-contents\" title=\"Link to this heading\">\u00b6</a></h2><h3>Classes<a class=\"headerlink\" href=\"#classes\" title=\"Link to this heading\">\u00b6</a></h3>", "a[href=\"#vasca.tables.TableCollection.remove_unselected\"]": "<dt class=\"sig sig-object py\" id=\"vasca.tables.TableCollection.remove_unselected\">\n<span class=\"sig-name descname\"><span class=\"pre\">remove_unselected</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">table_name</span></span><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"n\"><a class=\"reference external\" href=\"https://docs.python.org/3/library/stdtypes.html#str\" title=\"(in Python v3.13)\"><span class=\"pre\">str</span></a></span></em><span class=\"sig-paren\">)</span> <span class=\"sig-return\"><span class=\"sig-return-icon\">\u2192</span> <span class=\"sig-return-typehint\"><a class=\"reference external\" href=\"https://docs.python.org/3/library/constants.html#None\" title=\"(in Python v3.13)\"><span class=\"pre\">None</span></a></span></span><a class=\"reference internal\" href=\"../../_modules/vasca/tables.html#TableCollection.remove_unselected\"><span class=\"viewcode-link\"><span class=\"pre\">[source]</span></span></a></dt><dd><p>Remove unselected rows from given table</p></dd>", "a[href=\"#vasca.tables.TableCollection.select_rows\"]": "<dt class=\"sig sig-object py\" id=\"vasca.tables.TableCollection.select_rows\">\n<span class=\"sig-name descname\"><span class=\"pre\">select_rows</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">selections</span></span><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"n\"><a class=\"reference external\" href=\"https://docs.python.org/3/library/stdtypes.html#dict\" title=\"(in Python v3.13)\"><span class=\"pre\">dict</span></a></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">remove_unselected</span></span><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"n\"><a class=\"reference external\" href=\"https://docs.python.org/3/library/functions.html#bool\" title=\"(in Python v3.13)\"><span class=\"pre\">bool</span></a></span><span class=\"w\"> </span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"w\"> </span><span class=\"default_value\"><span class=\"pre\">False</span></span></em><span class=\"sig-paren\">)</span> <span class=\"sig-return\"><span class=\"sig-return-icon\">\u2192</span> <span class=\"sig-return-typehint\"><a class=\"reference external\" href=\"https://docs.python.org/3/library/constants.html#None\" title=\"(in Python v3.13)\"><span class=\"pre\">None</span></a></span></span><a class=\"reference internal\" href=\"../../_modules/vasca/tables.html#TableCollection.select_rows\"><span class=\"viewcode-link\"><span class=\"pre\">[source]</span></span></a></dt><dd><p>Apply selection to a passed table.</p></dd>", "a[href=\"#vasca.tables.FILE_DIR\"]": "<dt class=\"sig sig-object py\" id=\"vasca.tables.FILE_DIR\">\n<span class=\"sig-prename descclassname\"><span class=\"pre\">vasca.tables.</span></span><span class=\"sig-name descname\"><span class=\"pre\">FILE_DIR</span></span><em class=\"property\"><span class=\"w\"> </span><span class=\"p\"><span class=\"pre\">=</span></span><span class=\"w\"> </span><span class=\"pre\">'dirname(...)'</span></em></dt><dd></dd>", "a[href=\"#vasca.tables.TableCollection.__str__\"]": "<dt class=\"sig sig-object py\" id=\"vasca.tables.TableCollection.__str__\">\n<span class=\"sig-name descname\"><span class=\"pre\">__str__</span></span><span class=\"sig-paren\">(</span><span class=\"sig-paren\">)</span> <span class=\"sig-return\"><span class=\"sig-return-icon\">\u2192</span> <span class=\"sig-return-typehint\"><a class=\"reference external\" href=\"https://docs.python.org/3/library/stdtypes.html#str\" title=\"(in Python v3.13)\"><span class=\"pre\">str</span></a></span></span><a class=\"reference internal\" href=\"../../_modules/vasca/tables.html#TableCollection.__str__\"><span class=\"viewcode-link\"><span class=\"pre\">[source]</span></span></a></dt><dd><p>Return string with information about the collection.</p></dd>", "a[href=\"#data\"]": "<h3 class=\"tippy-header\" style=\"margin-top: 0;\">Data<a class=\"headerlink\" href=\"#data\" title=\"Link to this heading\">\u00b6</a></h3>", "a[href=\"#vasca.tables.TableCollection.copy_table_columns\"]": "<dt class=\"sig sig-object py\" id=\"vasca.tables.TableCollection.copy_table_columns\">\n<span class=\"sig-name descname\"><span class=\"pre\">copy_table_columns</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">tab_name_to</span></span><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"n\"><a class=\"reference external\" href=\"https://docs.python.org/3/library/stdtypes.html#str\" title=\"(in Python v3.13)\"><span class=\"pre\">str</span></a></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">tab_name_from</span></span><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"n\"><a class=\"reference external\" href=\"https://docs.python.org/3/library/stdtypes.html#str\" title=\"(in Python v3.13)\"><span class=\"pre\">str</span></a></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">copy_vars</span></span><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"n\"><a class=\"reference external\" href=\"https://docs.python.org/3/library/stdtypes.html#list\" title=\"(in Python v3.13)\"><span class=\"pre\">list</span></a><span class=\"p\"><span class=\"pre\">[</span></span><a class=\"reference external\" href=\"https://docs.python.org/3/library/stdtypes.html#str\" title=\"(in Python v3.13)\"><span class=\"pre\">str</span></a><span class=\"p\"><span class=\"pre\">]</span></span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">match_var</span></span><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"n\"><a class=\"reference external\" href=\"https://docs.python.org/3/library/stdtypes.html#str\" title=\"(in Python v3.13)\"><span class=\"pre\">str</span></a></span><span class=\"w\"> </span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"w\"> </span><span class=\"default_value\"><span class=\"pre\">'rg_src_id'</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">select_matched</span></span><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"n\"><a class=\"reference external\" href=\"https://docs.python.org/3/library/functions.html#bool\" title=\"(in Python v3.13)\"><span class=\"pre\">bool</span></a></span><span class=\"w\"> </span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"w\"> </span><span class=\"default_value\"><span class=\"pre\">False</span></span></em><span class=\"sig-paren\">)</span> <span class=\"sig-return\"><span class=\"sig-return-icon\">\u2192</span> <span class=\"sig-return-typehint\"><a class=\"reference external\" href=\"https://docs.python.org/3/library/constants.html#None\" title=\"(in Python v3.13)\"><span class=\"pre\">None</span></a></span></span><a class=\"reference internal\" href=\"../../_modules/vasca/tables.html#TableCollection.copy_table_columns\"><span class=\"viewcode-link\"><span class=\"pre\">[source]</span></span></a></dt><dd><p>Copy a column from one table to the other, for those columns that have a\nmatching variable value \u2018match_variable\u2019</p></dd>", "a[href=\"#vasca.tables.TableCollection\"]": "<dt class=\"sig sig-object py\" id=\"vasca.tables.TableCollection\">\n<em class=\"property\"><span class=\"pre\">class</span><span class=\"w\"> </span></em><span class=\"sig-prename descclassname\"><span class=\"pre\">vasca.tables.</span></span><span class=\"sig-name descname\"><span class=\"pre\">TableCollection</span></span><a class=\"reference internal\" href=\"../../_modules/vasca/tables.html#TableCollection\"><span class=\"viewcode-link\"><span class=\"pre\">[source]</span></span></a></dt><dd><p>Collection of <a class=\"reference external\" href=\"https://docs.astropy.org/en/stable/api/astropy.table.Table.html#astropy.table.Table\" title=\"(in Astropy v6.1)\"><code class=\"xref py py-class docutils literal notranslate\"><span class=\"pre\">Table</span></code></a> objects. Base calss for data storage\nclasses in VASCA.</p><p class=\"rubric\">Initialization</p></dd>", "a[href=\"#vasca.tables.TableCollection.info\"]": "<dt class=\"sig sig-object py\" id=\"vasca.tables.TableCollection.info\">\n<span class=\"sig-name descname\"><span class=\"pre\">info</span></span><span class=\"sig-paren\">(</span><span class=\"sig-paren\">)</span> <span class=\"sig-return\"><span class=\"sig-return-icon\">\u2192</span> <span class=\"sig-return-typehint\"><a class=\"reference external\" href=\"https://docs.python.org/3/library/constants.html#None\" title=\"(in Python v3.13)\"><span class=\"pre\">None</span></a></span></span><a class=\"reference internal\" href=\"../../_modules/vasca/tables.html#TableCollection.info\"><span class=\"viewcode-link\"><span class=\"pre\">[source]</span></span></a></dt><dd><p>Print out information on tables in the collection.</p></dd>", "a[href=\"#vasca.tables.TableCollection.add_column\"]": "<dt class=\"sig sig-object py\" id=\"vasca.tables.TableCollection.add_column\">\n<span class=\"sig-name descname\"><span class=\"pre\">add_column</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">table_name</span></span><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"n\"><a class=\"reference external\" href=\"https://docs.python.org/3/library/stdtypes.html#str\" title=\"(in Python v3.13)\"><span class=\"pre\">str</span></a></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">col_name</span></span><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"n\"><a class=\"reference external\" href=\"https://docs.python.org/3/library/stdtypes.html#str\" title=\"(in Python v3.13)\"><span class=\"pre\">str</span></a></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">col_data</span></span><span class=\"p\"><span class=\"pre\">:</span></span><span class=\"w\"> </span><span class=\"n\"><a class=\"reference external\" href=\"https://docs.python.org/3/library/stdtypes.html#dict\" title=\"(in Python v3.13)\"><span class=\"pre\">dict</span></a><span class=\"w\"> </span><span class=\"p\"><span class=\"pre\">|</span></span><span class=\"w\"> </span><a class=\"reference external\" href=\"https://numpy.org/doc/stable/reference/generated/numpy.ndarray.html#numpy.ndarray\" title=\"(in NumPy v2.1)\"><span class=\"pre\">numpy.ndarray</span></a><span class=\"w\"> </span><span class=\"p\"><span class=\"pre\">|</span></span><span class=\"w\"> </span><a class=\"reference external\" href=\"https://docs.python.org/3/library/constants.html#None\" title=\"(in Python v3.13)\"><span class=\"pre\">None</span></a></span><span class=\"w\"> </span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"w\"> </span><span class=\"default_value\"><span class=\"pre\">None</span></span></em><span class=\"sig-paren\">)</span> <span class=\"sig-return\"><span class=\"sig-return-icon\">\u2192</span> <span class=\"sig-return-typehint\"><a class=\"reference external\" href=\"https://docs.python.org/3/library/constants.html#None\" title=\"(in Python v3.13)\"><span class=\"pre\">None</span></a></span></span><a class=\"reference internal\" href=\"../../_modules/vasca/tables.html#TableCollection.add_column\"><span class=\"viewcode-link\"><span class=\"pre\">[source]</span></span></a></dt><dd><p>Adds column in a table, using the predefined VASCA columns.\nIf column exists already replace it.</p></dd>"}
skip_classes = ["headerlink", "sd-stretched-link"]

window.onload = function () {
    for (const [select, tip_html] of Object.entries(selector_to_html)) {
        const links = document.querySelectorAll(`div.content ${select}`);
        for (const link of links) {
            if (skip_classes.some(c => link.classList.contains(c))) {
                continue;
            }

            tippy(link, {
                content: tip_html,
                allowHTML: true,
                arrow: true,
                placement: 'auto-start', maxWidth: 500, interactive: false,
                onShow(instance) {MathJax.typesetPromise([instance.popper]).then(() => {});},
            });
        };
    };
    console.log("tippy tips loaded!");
};
