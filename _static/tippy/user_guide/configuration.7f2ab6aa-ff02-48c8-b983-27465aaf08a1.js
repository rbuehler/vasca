selector_to_html = {"a[href=\"#configuration\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">Configuration<a class=\"headerlink\" href=\"#configuration\" title=\"Link to this heading\">\u00b6</a></h1><p>VASCA\u2019s pipeline can be configured using yaml files. See <a class=\"reference external\" href=\"https://github.com/rbuehler/vasca/blob/main/vasca/vasca_cfg.yaml\"><code class=\"docutils literal notranslate\"><span class=\"pre\">vasca/vasca_cfg.yml</span></code></a>\nas an example below. For a more concrete example see the pipeline tutorial <a class=\"reference internal\" href=\"../tutorials/tutorial_pipe.html\"><span class=\"std std-doc\">here</span></a>.</p>", "a[href=\"../tutorials/tutorial_pipe.html\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">Pipeline<a class=\"headerlink\" href=\"#pipeline\" title=\"Link to this heading\">\u00b6</a></h1><p>This is a tutorial showcasing VASCA\u2019s pipeline flow on a simple example. We will go\nthrough all the steps equivalent to what is done in <a class=\"reference internal\" href=\"../api/vasca/vasca.vasca_pipe.html#vasca.vasca_pipe.run_from_file\" title=\"vasca.vasca_pipe.run_from_file\"><code class=\"xref myst py py-func docutils literal notranslate\"><span class=\"pre\">vasca_pipe.run_from_file</span></code></a>.\nThis is the same function that is called when starting the pipeline from the CLI using <code class=\"docutils literal notranslate\"><span class=\"pre\">vasca-pipe</span></code>.</p><p>The goal is to create a VASCA <a class=\"reference internal\" href=\"../api/vasca/vasca.region.html#vasca.region.Region\" title=\"vasca.region.Region\"><code class=\"xref myst py py-class docutils literal notranslate\"><span class=\"pre\">Region</span></code></a> from multiple <a class=\"reference internal\" href=\"../api/vasca/vasca.field.html#vasca.field.GALEXField\" title=\"vasca.field.GALEXField\"><code class=\"xref myst py py-class docutils literal notranslate\"><span class=\"pre\">GALEXField</span></code></a> for which we\ndownload the raw data online from <a class=\"reference external\" href=\"https://astroquery.readthedocs.io/en/latest/mast/mast.html\">MAST</a>.\nWe apply quality cuts and do source clustering followed by variability analysis.</p>"}
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
