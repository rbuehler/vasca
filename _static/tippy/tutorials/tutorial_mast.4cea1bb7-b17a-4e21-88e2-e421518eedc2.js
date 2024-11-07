selector_to_html = {"a[href=\"#mast-download\"]": "<h1 class=\"tippy-header\" style=\"margin-top: 0;\">MAST Download<a class=\"headerlink\" href=\"#mast-download\" title=\"Link to this heading\">\u00b6</a></h1><p>In this short tutorial the MAST query functions of <a class=\"reference internal\" href=\"../api/vasca/vasca.field.html#vasca.field.GALEXField\" title=\"vasca.field.GALEXField\"><code class=\"xref myst py py-class docutils literal notranslate\"><span class=\"pre\">GALEXField</span></code></a> are used to\ndownload fresh data.</p>", "a[href=\"../api/vasca/vasca.field.html#vasca.field.GALEXField\"]": "<dt class=\"sig sig-object py\" id=\"vasca.field.GALEXField\">\n<em class=\"property\"><span class=\"pre\">class</span><span class=\"w\"> </span></em><span class=\"sig-prename descclassname\"><span class=\"pre\">vasca.field.</span></span><span class=\"sig-name descname\"><span class=\"pre\">GALEXField</span></span><span class=\"sig-paren\">(</span><em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">obs_id</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">obs_filter</span></span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"default_value\"><span class=\"pre\">None</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">data_path</span></span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"default_value\"><span class=\"pre\">None</span></span></em>, <em class=\"sig-param\"><span class=\"n\"><span class=\"pre\">visits_data_path</span></span><span class=\"o\"><span class=\"pre\">=</span></span><span class=\"default_value\"><span class=\"pre\">None</span></span></em><span class=\"sig-paren\">)</span><a class=\"reference internal\" href=\"../_modules/vasca/field.html#GALEXField\"><span class=\"viewcode-link\"><span class=\"pre\">[source]</span></span></a></dt><dd><p>Bases: <a class=\"reference internal\" href=\"#vasca.field.BaseField\" title=\"vasca.field.BaseField\"><code class=\"xref py py-obj docutils literal notranslate\"><span class=\"pre\">vasca.field.BaseField</span></code></a></p><p>Instance of one GALEX field</p><p class=\"rubric\">Initialization</p><p>Initializes a new GALEXField instance with\nskeleton VASCA data structure.</p></dd>"}
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
