selector_to_html = {"a[href=\"#classes\"]": "<h3 class=\"tippy-header\" style=\"margin-top: 0;\">Classes<a class=\"headerlink\" href=\"#classes\" title=\"Link to this heading\">\u00b6</a></h3>", "a[href=\"#vasca.region.Region.add_simbad_otype_info\"]": "<dt class=\"sig sig-object py\" id=\"vasca.region.Region.add_simbad_otype_info\">\n<span class=\"sig-name descname\"><span class=\"pre\">add_simbad_otype_info</span></span><span class=\"sig-paren\">(</span><span class=\"sig-paren\">)</span></dt><dd><p>Add table explaing SIMBAD object groups</p></dd>", "a[href=\"#vasca.region.Region.set_src_id_info\"]": "<dt class=\"sig sig-object py\" id=\"vasca.region.Region.set_src_id_info\">\n<span class=\"sig-name descname\"><span class=\"pre\">set_src_id_info</span></span><span class=\"sig-paren\">(</span><span class=\"sig-paren\">)</span></dt><dd><p>Stores the mapping of rg_src_id to rg_fd_id and fd_src_id\ninto tt_src_id_map table.</p><p>None</p></dd>", "a[href=\"#api\"]": "<h3 class=\"tippy-header\" style=\"margin-top: 0;\">API<a class=\"headerlink\" href=\"#api\" title=\"Link to this heading\">\u00b6</a></h3>", "a[href=\"#vasca.region.Region.redo_src_selection\"]": "<dt class=\"sig sig-object py\" id=\"vasca.region.Region.redo_src_selection\">\n<span class=\"sig-name descname\"><span class=\"pre\">redo_src_selection</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">cfg_file_name</span></span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"default_value\"><span class=\"pre\">'./vasca_cfg.yaml'</span></span></em><span class=\"sig-paren\">)</span></dt><dd><p>Redo source selection, set in the tt_sources[\u201csel\u201d] column\nof the region based on the passed configuration file.</p><p>None</p></dd>", "a[href=\"#vasca.region.Region.add_table_from_fields\"]": "<dt class=\"sig sig-object py\" id=\"vasca.region.Region.add_table_from_fields\">\n<span class=\"sig-name descname\"><span class=\"pre\">add_table_from_fields</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">table_name</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">only_selected</span></span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"default_value\"><span class=\"pre\">False</span></span></em><span class=\"sig-paren\">)</span></dt><dd><p>Add tables from the fields to the region by stacking them,\nadding the rg_fd_id column.</p><p>None</p></dd>", "a[href=\"#vasca.region.Region\"]": "<dt class=\"sig sig-object py\" id=\"vasca.region.Region\">\n<em class=\"property\"><span class=\"pre\">class</span><span class=\"w\"> </span></em><span class=\"sig-prename descclassname\"><span class=\"pre\">vasca.region.</span></span><span class=\"sig-name descname\"><span class=\"pre\">Region</span></span></dt><dd><p>Bases: <a class=\"reference internal\" href=\"vasca.tables.html#vasca.tables.TableCollection\" title=\"vasca.tables.TableCollection\"><code class=\"xref py py-obj docutils literal notranslate\"><span class=\"pre\">vasca.tables.TableCollection</span></code></a></p><p>Defines a region in the sky as a\nlist of vasca.field objects. It provides functionality to\nloop over fields to derive source lists, etc.</p><p class=\"rubric\">Initialization</p><p>Many class attributes are stored in <a class=\"reference external\" href=\"https://docs.astropy.org/en/stable/api/astropy.table.Table.html\">astropy.table.Table</a>. To see a\ndescription of each of their columns run :meth: <cite>~vasca.Regions.info</cite>.</p><p>None</p></dd>", "a[href=\"#vasca.region.Region.cross_match_cds\"]": "<dt class=\"sig sig-object py\" id=\"vasca.region.Region.cross_match_cds\">\n<span class=\"sig-name descname\"><span class=\"pre\">cross_match_cds</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">query_radius</span></span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"default_value\"><span class=\"pre\">1.5</span> <span class=\"pre\">*</span> <span class=\"pre\">uu.arcsec</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">query_table</span></span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"default_value\"><span class=\"pre\">'I/355/gaiadr3'</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">vizier_columns</span></span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"default_value\"><span class=\"pre\">['*',</span> <span class=\"pre\">'PQSO',</span> <span class=\"pre\">'PGal',</span> <span class=\"pre\">'PSS',</span> <span class=\"pre\">'RPlx',</span> <span class=\"pre\">'VarFlag',</span> <span class=\"pre\">'o_Gmag',</span> <span class=\"pre\">'RFRP',</span> <span class=\"pre\">'RFBP',</span> <span class=\"pre\">'AG',</span> <span class=\"pre\">'E(BP-RP)']</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">overwrite</span></span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"default_value\"><span class=\"pre\">False</span></span></em><span class=\"sig-paren\">)</span></dt><dd><p>Match sources in region with SIMBAD-catalogs or Vizier database catalog. Runs\nonly over selected sources.</p><p>None</p></dd>", "a[href=\"#vasca.region.Region.load_from_config\"]": "<dt class=\"sig sig-object py\" id=\"vasca.region.Region.load_from_config\">\n<em class=\"property\"><span class=\"pre\">classmethod</span><span class=\"w\"> </span></em><span class=\"sig-name descname\"><span class=\"pre\">load_from_config</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">vasca_cfg</span></span></em><span class=\"sig-paren\">)</span></dt><dd><p>Loads region from configuration dictionary.</p><p>None</p></dd>", "a[href=\"vasca.tables.html#vasca.tables.TableCollection\"]": "<dt class=\"sig sig-object py\" id=\"vasca.tables.TableCollection\">\n<em class=\"property\"><span class=\"pre\">class</span><span class=\"w\"> </span></em><span class=\"sig-prename descclassname\"><span class=\"pre\">vasca.tables.</span></span><span class=\"sig-name descname\"><span class=\"pre\">TableCollection</span></span><a class=\"reference internal\" href=\"../../_modules/vasca/tables.html#TableCollection\"><span class=\"viewcode-link\"><span class=\"pre\">[source]</span></span></a></dt><dd><p>Bases: <a class=\"reference external\" href=\"https://docs.python.org/3/library/functions.html#object\" title=\"(in Python v3.13)\"><code class=\"xref py py-obj docutils literal notranslate\"><span class=\"pre\">object</span></code></a></p><p>Collection of <a href=\"#id3\"><span class=\"problematic\" id=\"id4\">astropy.table.Table_</span></a> objects. Base calss for data\nstorage classes in VASCA.</p><p class=\"rubric\">Initialization</p></dd>", "a[href=\"#module-vasca.region\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\"><a class=\"reference internal\" href=\"#module-vasca.region\" title=\"vasca.region\"><code class=\"xref py py-mod docutils literal notranslate\"><span class=\"pre\">vasca.region</span></code></a><a class=\"headerlink\" href=\"#module-vasca.region\" title=\"Link to this heading\">\u00b6</a></h1><h2>Module Contents<a class=\"headerlink\" href=\"#module-contents\" title=\"Link to this heading\">\u00b6</a></h2><h3>Classes<a class=\"headerlink\" href=\"#classes\" title=\"Link to this heading\">\u00b6</a></h3>", "a[href=\"#vasca.region.Region.load_from_fits\"]": "<dt class=\"sig sig-object py\" id=\"vasca.region.Region.load_from_fits\">\n<span class=\"sig-name descname\"><span class=\"pre\">load_from_fits</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">file_name</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">load_fields</span></span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"default_value\"><span class=\"pre\">False</span></span></em><span class=\"sig-paren\">)</span></dt><dd><p>Loads field from a fits file.</p><p>None</p></dd>", "a[href=\"#vasca.region.Region.add_coverage_hp\"]": "<dt class=\"sig sig-object py\" id=\"vasca.region.Region.add_coverage_hp\">\n<span class=\"sig-name descname\"><span class=\"pre\">add_coverage_hp</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">nside</span></span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"default_value\"><span class=\"pre\">4096</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">coord_sys</span></span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"default_value\"><span class=\"pre\">'galactic'</span></span></em><span class=\"sig-paren\">)</span></dt><dd><p>Creates healpix arrays of Nr visits, fields and total exposure.</p></dd>", "a[href=\"#vasca.region.Region.set_LombScargle\"]": "<dt class=\"sig sig-object py\" id=\"vasca.region.Region.set_LombScargle\">\n<span class=\"sig-name descname\"><span class=\"pre\">set_LombScargle</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">obs_filters</span></span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"default_value\"><span class=\"pre\">['NUV',</span> <span class=\"pre\">'FUV']</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">nbins_min</span></span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"default_value\"><span class=\"pre\">20</span></span></em><span class=\"sig-paren\">)</span></dt><dd><p>Apply LombScargle analysis to selected sources. Results\nare stored into the tt_lombscargle table of the region.</p><p>None</p></dd>", "a[href=\"#vasca.region.Region.get_src_from_id\"]": "<dt class=\"sig sig-object py\" id=\"vasca.region.Region.get_src_from_id\">\n<span class=\"sig-name descname\"><span class=\"pre\">get_src_from_id</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">rg_src_id</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">load_from_file</span></span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"default_value\"><span class=\"pre\">True</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">write_to_file</span></span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"default_value\"><span class=\"pre\">True</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">add_sed</span></span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"default_value\"><span class=\"pre\">True</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">add_gphoton</span></span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"default_value\"><span class=\"pre\">True</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">add_spectrum</span></span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"default_value\"><span class=\"pre\">True</span></span></em><span class=\"sig-paren\">)</span></dt><dd><p>Get Source object containing all region table entries\nrelevant for the passed rg_src_id.</p></dd>", "a[href=\"#vasca.region.Region.synch_src_sel\"]": "<dt class=\"sig sig-object py\" id=\"vasca.region.Region.synch_src_sel\">\n<span class=\"sig-name descname\"><span class=\"pre\">synch_src_sel</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">remove_unselected</span></span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"default_value\"><span class=\"pre\">False</span></span></em><span class=\"sig-paren\">)</span></dt><dd><p>Synchronize selections among tables. Select rows for tables containing\n\u201crg_src_id\u201d for sources selected in tt_sources.</p><p>None</p></dd>", "a[href=\"#vasca.region.Region.get_region_catalog\"]": "<dt class=\"sig sig-object py\" id=\"vasca.region.Region.get_region_catalog\">\n<span class=\"sig-name descname\"><span class=\"pre\">get_region_catalog</span></span><span class=\"sig-paren\">(</span><span class=\"sig-paren\">)</span></dt><dd><p>Create a reduced region, which only contains info on selected sources</p></dd>", "a[href=\"#vasca.region.Region.get_field\"]": "<dt class=\"sig sig-object py\" id=\"vasca.region.Region.get_field\">\n<span class=\"sig-name descname\"><span class=\"pre\">get_field</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">field_id</span></span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"default_value\"><span class=\"pre\">None</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">rg_fd_id</span></span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"default_value\"><span class=\"pre\">None</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">load_method</span></span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"default_value\"><span class=\"pre\">'FITS'</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">add_field</span></span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"default_value\"><span class=\"pre\">False</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">mast_products</span></span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"default_value\"><span class=\"pre\">'TABLES'</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">field_kwargs</span></span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"default_value\"><span class=\"pre\">dict()</span></span></em><span class=\"sig-paren\">)</span></dt><dd><p>Load a field from a region, tt_fields table needs to include this field.</p></dd>", "a[href=\"#vasca.region.Region.get_src_from_sky_pos\"]": "<dt class=\"sig sig-object py\" id=\"vasca.region.Region.get_src_from_sky_pos\">\n<span class=\"sig-name descname\"><span class=\"pre\">get_src_from_sky_pos</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">coordx</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">coordy</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">frame</span></span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"default_value\"><span class=\"pre\">'icrs'</span></span></em><span class=\"sig-paren\">)</span></dt><dd><p>Get Source object containing all region table entries\nrelevant for the passed source position. The nearest source is matched.</p></dd>", "a[href=\"#module-contents\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Module Contents<a class=\"headerlink\" href=\"#module-contents\" title=\"Link to this heading\">\u00b6</a></h2><h3>Classes<a class=\"headerlink\" href=\"#classes\" title=\"Link to this heading\">\u00b6</a></h3>"}
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
