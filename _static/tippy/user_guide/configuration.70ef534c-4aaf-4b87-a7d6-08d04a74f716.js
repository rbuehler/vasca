selector_to_html = {"a[href=\"#configuration\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">Configuration<a class=\"headerlink\" href=\"#configuration\" title=\"Link to this heading\">\u00b6</a></h1><p>VASCA\u2019s pipeline can be configured using yaml files. See <a class=\"reference external\" href=\"https://github.com/rbuehler/vasca/blob/main/vasca/vasca_cfg.yaml\"><code class=\"docutils literal notranslate\"><span class=\"pre\">vasca/vasca_cfg.yml</span></code></a>\nas an example:</p>"}
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
