---
file_format: mystnb
kernelspec:
    name: vasca-github
---
# Tutorials
In this section various tutorials are provided. This is a markdown file and this section
is still under development.

## Intro
All tutorials are jupyter-based. This documentation uses [`myst-nb`](https://myst-nb.readthedocs.io/en/latest/index.html)
and [`jupytext`](https://jupytext.readthedocs.io/en/latest/index.html) to execute and
render the content.

It is quite remarkable since I can use inline code cells as well:

```{code-cell}
# Set up a random generator with fixed seed
import numpy as np

seed = 123
rng = np.random.default_rng(seed)
```

Now let's see what we can do with it:
```{code-cell}
# Print a random number. This should be the same number everytime this cell is evaluated.
rng.random()
```

## List of Tutorials
```{toctree}
:maxdepth: 2
:titlesonly:

simple_example
```