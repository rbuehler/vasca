selector_to_html = {"a[href=\"#processing-flow\"]": "<figure class=\"align-default\" id=\"processing-flow\">\n<a class=\"bg-primary mb-1 reference internal image-reference\" href=\"../_images/VASCA_processing_flow_v2.jpg\"><img alt=\"data_model\" class=\"bg-primary mb-1\" src=\"../_images/VASCA_processing_flow_v2.jpg\" style=\"width: 600px;\"/></a>\n<figcaption>\n<p><span class=\"caption-text\">The pipeline processing flow in VASCA</span><a class=\"headerlink\" href=\"#processing-flow\" title=\"Link to this image\">\u00b6</a></p>\n</figcaption>\n</figure>", "a[href=\"#pipeline\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">Pipeline<a class=\"headerlink\" href=\"#pipeline\" title=\"Link to this heading\">\u00b6</a></h1><p>VASCA\u2019s pipeline flow is composed of individual components for which\nparallel processing is enabled.</p>"}
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
