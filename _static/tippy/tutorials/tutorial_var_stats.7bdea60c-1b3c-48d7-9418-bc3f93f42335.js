selector_to_html = {"a[href=\"#run-example\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Run Example<a class=\"headerlink\" href=\"#run-example\" title=\"Link to this heading\">\u00b6</a></h2>", "a[href=\"#variability-statistics\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">Variability Statistics<a class=\"headerlink\" href=\"#variability-statistics\" title=\"Link to this heading\">\u00b6</a></h1><p>This tutorial showcases VASCA\u2019s main statistics computation (<a class=\"reference internal\" href=\"../api/vasca/vasca.utils.html#vasca.utils.get_var_stat\" title=\"vasca.utils.get_var_stat\"><code class=\"xref myst py py-func docutils literal notranslate\"><span class=\"pre\">get_var_stat</span></code></a>) based\non synthetic data. The simple simulation framework allows for the variation of data\npoints, flux uncertainties, number of flare events as well as intrinsic variability\nparameters like amplitude and period.</p>", "a[href=\"#synthetic-light-curves\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Synthetic Light Curves<a class=\"headerlink\" href=\"#synthetic-light-curves\" title=\"Link to this heading\">\u00b6</a></h2>", "a[href=\"#custom-light-curves\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Custom Light Curves<a class=\"headerlink\" href=\"#custom-light-curves\" title=\"Link to this heading\">\u00b6</a></h2>", "a[href=\"#visualization\"]": "<h2 class=\"tippy-header\" style=\"margin-top: 0;\">Visualization<a class=\"headerlink\" href=\"#visualization\" title=\"Link to this heading\">\u00b6</a></h2>", "a[href=\"../api/vasca/vasca.utils.html#vasca.utils.get_var_stat\"]": "<dt class=\"sig sig-object py\" id=\"vasca.utils.get_var_stat\">\n<span class=\"sig-prename descclassname\"><span class=\"pre\">vasca.utils.</span></span><span class=\"sig-name descname\"><span class=\"pre\">get_var_stat</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">vals</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">vals_err</span></span></em><span class=\"sig-paren\">)</span><a class=\"reference internal\" href=\"../_modules/vasca/utils.html#get_var_stat\"><span class=\"viewcode-link\"><span class=\"pre\">[source]</span></span></a></dt><dd><p>Calculate variability parameters</p></dd>"}
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
