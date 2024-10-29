selector_to_html = {"a[href=\"#intro\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Intro<a class=\"headerlink\" href=\"#intro\" title=\"Link to this heading\">\u00b6</a></h2><p>All tutorials are jupyter-based. This documentation uses <a class=\"reference external\" href=\"https://myst-nb.readthedocs.io/en/latest/index.html\"><code class=\"docutils literal notranslate\"><span class=\"pre\">myst-nb</span></code></a>\nand <a class=\"reference external\" href=\"https://jupytext.readthedocs.io/en/latest/index.html\"><code class=\"docutils literal notranslate\"><span class=\"pre\">jupytext</span></code></a> to execute and\nrender the content.</p><p>It is quite remarkable since I can use inline code cells as well:</p>", "a[href=\"simple_example.html\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">Test Tutorial<a class=\"headerlink\" href=\"#test-tutorial\" title=\"Link to this heading\">\u00b6</a></h1><p>The contents of this page are edited in a python file which is converted to a markdown\nfile prior to the sphinx build and then executed during build time. See how long it\ntook to run this notebook <a class=\"reference internal\" href=\"#execution-statistics\">below</a>.</p>", "a[href=\"#tutorials\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">Tutorials<a class=\"headerlink\" href=\"#tutorials\" title=\"Link to this heading\">\u00b6</a></h1><p>In this section various tutorials are provided. This is a markdown file and this section\nis still under development.</p>", "a[href=\"test_table.html\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">Table Display<a class=\"headerlink\" href=\"#table-display\" title=\"Link to this heading\">\u00b6</a></h1><p>This test tutorial shows a Pandas DataFrame as an interactive table.</p>", "a[href=\"#list-of-tutorials\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">List of Tutorials<a class=\"headerlink\" href=\"#list-of-tutorials\" title=\"Link to this heading\">\u00b6</a></h2>"}
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
