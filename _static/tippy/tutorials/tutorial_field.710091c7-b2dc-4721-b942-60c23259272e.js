selector_to_html = {"a[href=\"../api/vasca/vasca.tables.html#vasca.tables.TableCollection.info\"]": "<dt class=\"sig sig-object py\" id=\"vasca.tables.TableCollection.info\">\n<span class=\"sig-name descname\"><span class=\"pre\">info</span></span><span class=\"sig-paren\">(</span><span class=\"sig-paren\">)</span> <span class=\"sig-return\"><span class=\"sig-return-icon\">\u2192</span> <span class=\"sig-return-typehint\"><a class=\"reference external\" href=\"https://docs.python.org/3/library/constants.html#None\" title=\"(in Python v3.13)\"><span class=\"pre\">None</span></a></span></span><a class=\"reference internal\" href=\"../_modules/vasca/tables.html#TableCollection.info\"><span class=\"viewcode-link\"><span class=\"pre\">[source]</span></span></a></dt><dd><p>Print out information on tables in the collection.</p></dd>", "a[href=\"../api/vasca/vasca.field.html#vasca.field.BaseField\"]": "<dt class=\"sig sig-object py\" id=\"vasca.field.BaseField\">\n<em class=\"property\"><span class=\"pre\">class</span><span class=\"w\"> </span></em><span class=\"sig-prename descclassname\"><span class=\"pre\">vasca.field.</span></span><span class=\"sig-name descname\"><span class=\"pre\">BaseField</span></span><a class=\"reference internal\" href=\"../_modules/vasca/field.html#BaseField\"><span class=\"viewcode-link\"><span class=\"pre\">[source]</span></span></a></dt><dd><p>Bases: <a class=\"reference internal\" href=\"../api/vasca/vasca.tables.html#vasca.tables.TableCollection\" title=\"vasca.tables.TableCollection\"><code class=\"xref py py-obj docutils literal notranslate\"><span class=\"pre\">vasca.tables.TableCollection</span></code></a></p><p>Class that defines the basic\ndata structure for field-based analysis. One <em>field</em> is generally\nthe area in the sky covered by a telescope in one observation.\nA field is generally composed of several <em>visits</em> of the telescope\nat different times.</p><p>This class contains the main functionality for source\ndetection and drawing. To be inherited by field analysis classes,\nwhich can then be tailored to the needs of the observatories supported\nby the VASCA pipeline.</p><p class=\"rubric\">Initialization</p><p>Many class attributes are stored in <a class=\"reference external\" href=\"https://docs.astropy.org/en/stable/api/astropy.table.Table.html\">astropy.table.Table</a>. To see a\ndescription of each of their columns run :meth: <cite>~vasca.field.BaseField.info</cite>.</p><p>None</p></dd>", "a[href=\"../api/vasca/vasca.field.html#vasca.field.GALEXField\"]": "<dt class=\"sig sig-object py\" id=\"vasca.field.GALEXField\">\n<em class=\"property\"><span class=\"pre\">class</span><span class=\"w\"> </span></em><span class=\"sig-prename descclassname\"><span class=\"pre\">vasca.field.</span></span><span class=\"sig-name descname\"><span class=\"pre\">GALEXField</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">obs_id</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">obs_filter</span></span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"default_value\"><span class=\"pre\">None</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">data_path</span></span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"default_value\"><span class=\"pre\">None</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">visits_data_path</span></span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"default_value\"><span class=\"pre\">None</span></span></em><span class=\"sig-paren\">)</span><a class=\"reference internal\" href=\"../_modules/vasca/field.html#GALEXField\"><span class=\"viewcode-link\"><span class=\"pre\">[source]</span></span></a></dt><dd><p>Bases: <a class=\"reference internal\" href=\"#vasca.field.BaseField\" title=\"vasca.field.BaseField\"><code class=\"xref py py-obj docutils literal notranslate\"><span class=\"pre\">vasca.field.BaseField</span></code></a></p><p>Instance of one GALEX field</p><p class=\"rubric\">Initialization</p><p>Initializes a new GALEXField instance with\nskeleton VASCA data structure.</p></dd>", "a[href=\"#galexfield\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">GALEXField<a class=\"headerlink\" href=\"#galexfield\" title=\"Link to this heading\">\u00b6</a></h1><p>This is a short tutorial on the <a class=\"reference internal\" href=\"../api/vasca/vasca.field.html#vasca.field.GALEXField\" title=\"vasca.field.GALEXField\"><code class=\"xref myst py py-class docutils literal notranslate\"><span class=\"pre\">GALEXField</span></code></a> class. This is an implementation of\nVASCA\u2019s <a class=\"reference internal\" href=\"../api/vasca/vasca.field.html#vasca.field.BaseField\" title=\"vasca.field.BaseField\"><code class=\"xref myst py py-class docutils literal notranslate\"><span class=\"pre\">BaseField</span></code></a> tailored to the specific needs of GALEX archival data. It\nserves first and foremost as the interface between the instrument specifics, like raw\ndata handling and nomenclature, and the rest of VASCA\u2019s functionality. Important are\nthe loading functions that download raw data from <a class=\"reference external\" href=\"https://astroquery.readthedocs.io/en/latest/mast/mast.html\">MAST</a>\nand load already processed field data into memory.</p>", "a[href=\"#example\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Example<a class=\"headerlink\" href=\"#example\" title=\"Link to this heading\">\u00b6</a></h2><p>Let\u2019s load GALEX field data that is already downloaded (used for unit testing)</p>"}
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
