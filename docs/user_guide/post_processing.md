# Post-Processing
The value provided by VASCA not only comes from its data model and the source clustering
pipeline, but also from the wide-ranging suite of post-processing functionality. The most
important ones are highlighted below:

Variability Detection
: Statistical variability is measured by testing the flux variation as a function of time
against the hypothesis of constant flux. See the [{meth}`get_var_stat`](#get_var_stat)
method for implementation details.

Catalog Cross-Matching
: VASCA implements cross-matching searches for SIMBAD and Gaia catalogs out-of-the box.
These two catalogs allow for limited source classification of VASCA sources. Users may
even search in local catalog files for which an easy integration exists. A chance-
coincidence analysis can be performed to determine the probability of a matched source
to be associated to a catalog object by pure chance. See [{meth}`cross_match_cds`](#cross_match_cds)
for implementation details.

Lomb-Scargele Variability
: Using the [Lomb-Scargele](https://docs.astropy.org/en/stable/timeseries/lombscargle.html)
algorithm users can asses a light curve for periodic variability. This may help to
classify sources where no result is found in any of the cross-matching catalog searches.
See [{meth}`run_LombScargle`](#run_LombScargle) for implementation details.

Spectral Energy Distribution
: Using Vizier, the SED is queried for associated sources. A black body fit gives further
insight into correctness of the classification the potential physical processes driving
the variability. See [{meth}`query_vizier_sed`](#query_vizier_sed) for implementation
details.

Publication-Ready Source Catalog
: The final step in the post-processing chain is the creation of the final source
catalog. Publishing a requires to follow certain protocols and naming conventions of the
various tables, their columns and not unimportantly the source IDs. By design, this is
all covered in VASCA automatically. For more info see this Jupyter [example](https://github.com/rbuehler/vasca/blob/main/vasca/examples/vasca_pipe_post_process.ipynb). 

Data Visualization
: Most of the above will require some portion of manual work on the pipeline results.
This is streamlined in VASCA by a sequential set of Jupyter notebooks that guide users
through the post-processing and support by providing useful data visualization functions. 