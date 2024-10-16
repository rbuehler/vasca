---
myst:
  html_meta:
    "description lang=en": |
      Top-level documentation for VASCA package, with links to the rest of the site..
---
<!-- start Intro -->
```{include} ../README.md
:relative-docs: docs/
:relative-images:
```
<!-- end Intro -->

<!-- start Docs Overview -->
# Documentation!
```{toctree}
:maxdepth: 2
:caption: Contents:
   
getting_started
user_guide/index
```
<!-- end Docs Overview -->

# Modules
<!-- ```{eval-rst}
.. autosummary::

  vasca.resource_manager.ResourceManager
  vasca.field.GALEXField
``` -->
```{autodoc2-summary}
:renderer: myst

~vasca.resource_manager.ResourceManager
~vasca.field.GALEXField
```
