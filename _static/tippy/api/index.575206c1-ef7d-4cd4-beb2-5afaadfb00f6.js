selector_to_html = {"a[href=\"vasca/vasca.region.html\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\"><a class=\"reference internal\" href=\"#module-vasca.region\" title=\"vasca.region\"><code class=\"xref py py-mod docutils literal notranslate\"><span class=\"pre\">vasca.region</span></code></a><a class=\"headerlink\" href=\"#module-vasca.region\" title=\"Link to this heading\">\u00b6</a></h1><h2>Module Contents<a class=\"headerlink\" href=\"#module-contents\" title=\"Link to this heading\">\u00b6</a></h2><h3>Classes<a class=\"headerlink\" href=\"#classes\" title=\"Link to this heading\">\u00b6</a></h3>", "a[href=\"#api-reference\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">API Reference<a class=\"headerlink\" href=\"#api-reference\" title=\"Link to this heading\">\u00b6</a></h1><p>This page contains auto-generated API reference documentation <a class=\"footnote-reference brackets\" href=\"#f1\" id=\"id1\" role=\"doc-noteref\"><span class=\"fn-bracket\">[</span>1<span class=\"fn-bracket\">]</span></a>.</p>", "a[href=\"vasca/vasca.tables_dict.html\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\"><a class=\"reference internal\" href=\"#module-vasca.tables_dict\" title=\"vasca.tables_dict\"><code class=\"xref py py-mod docutils literal notranslate\"><span class=\"pre\">vasca.tables_dict</span></code></a><a class=\"headerlink\" href=\"#module-vasca.tables_dict\" title=\"Link to this heading\">\u00b6</a></h1><p>Defines dictionary for the tables used by vasca.tables.TableCollection</p>", "a[href=\"vasca/vasca.html\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\"><a class=\"reference internal\" href=\"#module-vasca\" title=\"vasca\"><code class=\"xref py py-mod docutils literal notranslate\"><span class=\"pre\">vasca</span></code></a><a class=\"headerlink\" href=\"#module-vasca\" title=\"Link to this heading\">\u00b6</a></h1><h2>Submodules<a class=\"headerlink\" href=\"#submodules\" title=\"Link to this heading\">\u00b6</a></h2>", "a[href=\"vasca/vasca.tables.html\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\"><a class=\"reference internal\" href=\"#module-vasca.tables\" title=\"vasca.tables\"><code class=\"xref py py-mod docutils literal notranslate\"><span class=\"pre\">vasca.tables</span></code></a><a class=\"headerlink\" href=\"#module-vasca.tables\" title=\"Link to this heading\">\u00b6</a></h1><h2>Module Contents<a class=\"headerlink\" href=\"#module-contents\" title=\"Link to this heading\">\u00b6</a></h2><h3>Classes<a class=\"headerlink\" href=\"#classes\" title=\"Link to this heading\">\u00b6</a></h3>", "a[href=\"vasca/vasca.vasca_pipe.html\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\"><a class=\"reference internal\" href=\"#module-vasca.vasca_pipe\" title=\"vasca.vasca_pipe\"><code class=\"xref py py-mod docutils literal notranslate\"><span class=\"pre\">vasca.vasca_pipe</span></code></a><a class=\"headerlink\" href=\"#module-vasca.vasca_pipe\" title=\"Link to this heading\">\u00b6</a></h1><p>Script that runs the VASCA pipeline.</p>", "a[href=\"#f1\"]": "<aside class=\"footnote brackets\" id=\"f1\" role=\"doc-footnote\">\n<span class=\"label\"><span class=\"fn-bracket\">[</span><a href=\"#id1\" role=\"doc-backlink\">1</a><span class=\"fn-bracket\">]</span></span>\n<p>Created with <a class=\"reference external\" href=\"https://github.com/chrisjsewell/sphinx-autodoc2\">sphinx-autodoc2</a></p>\n</aside>", "a[href=\"vasca/vasca.utils.html\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\"><a class=\"reference internal\" href=\"#module-vasca.utils\" title=\"vasca.utils\"><code class=\"xref py py-mod docutils literal notranslate\"><span class=\"pre\">vasca.utils</span></code></a><a class=\"headerlink\" href=\"#module-vasca.utils\" title=\"Link to this heading\">\u00b6</a></h1><p>Utilities for VASCA</p>", "a[href=\"vasca/vasca.field.html\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\"><a class=\"reference internal\" href=\"#module-vasca.field\" title=\"vasca.field\"><code class=\"xref py py-mod docutils literal notranslate\"><span class=\"pre\">vasca.field</span></code></a><a class=\"headerlink\" href=\"#module-vasca.field\" title=\"Link to this heading\">\u00b6</a></h1><p>Field classes for VASCA</p>", "a[href=\"vasca/vasca.visualization.html\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\"><a class=\"reference internal\" href=\"#module-vasca.visualization\" title=\"vasca.visualization\"><code class=\"xref py py-mod docutils literal notranslate\"><span class=\"pre\">vasca.visualization</span></code></a><a class=\"headerlink\" href=\"#module-vasca.visualization\" title=\"Link to this heading\">\u00b6</a></h1><p>Visualization related methods for VASCA</p>", "a[href=\"vasca/vasca.resource_manager.html\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\"><a class=\"reference internal\" href=\"#module-vasca.resource_manager\" title=\"vasca.resource_manager\"><code class=\"xref py py-mod docutils literal notranslate\"><span class=\"pre\">vasca.resource_manager</span></code></a><a class=\"headerlink\" href=\"#module-vasca.resource_manager\" title=\"Link to this heading\">\u00b6</a></h1><p>Resource manager for VASCA</p>", "a[href=\"vasca/vasca.source.html\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\"><a class=\"reference internal\" href=\"#module-vasca.source\" title=\"vasca.source\"><code class=\"xref py py-mod docutils literal notranslate\"><span class=\"pre\">vasca.source</span></code></a><a class=\"headerlink\" href=\"#module-vasca.source\" title=\"Link to this heading\">\u00b6</a></h1><p>Created on Thu Jan 19 10:17:54 2023</p><p>@author: buehler</p>"}
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
