---
jupytext:
  text_represenation:
    extension: .py
    format_name: percent
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.16.4
kernelspec:
  display_name: vasca-github
  language: python
  name: vasca-github
---

```{code-cell}
:tags: [remove-cell]

# ruff: noqa: T201
```

# Test Tutorial

The contents of this page are edited in a python file which is converted to a markdown
file prior to the sphinx build and then executed during build time. See how long it
took to run this notebook [below](#execution-statistics).

+++

## Simple Test
This is a simple test function
Let's try what happens to myst-style sphinx admonitions:
:::{hint}
Pairing `.py` files with `.ipynb` files can make version control of notebooks better.
The python file is used for version control the notebook file is excluded from it, but
can be used for developing the notebook contents in a proper Jupyter environment. This
allows to seamlessly use all the tools and extensions in my code editor the help me
develop Python code.
:::

```{code-cell}
# Just an inline code comment explaining that the function below is very simple.
def f(x: float) -> float:
    return 3 * x + 1
```

```{code-cell}
:tags: [raises-exception]

assert f(4) == 12
print(f(4))
```

```{code-cell}
:tags: [hide-output]

# The output of this cell will be collapsed by default
# using the `hide-output` cell metadata tag
print(1 == 1)
for i in range(20):
    print(f"Number: {i}")
```

```{code-cell}
%%time
import numpy as np

print(f"{np.pi:1.5f}")
```

```{code-cell}
# Simple plot example
import matplotlib.pyplot as plt

Fs = 8000
f = 5
sample = 8000
x = np.arange(sample)
y = np.sin(2 * np.pi * f * x / Fs)
plt.plot(x, y)
plt.xlabel("samples")
plt.ylabel("amplitude")
```

## Test Markdown Features

Lets put a relative link to one of the sections in the user guide in this documentation
[here](../user_guide/index.md#using-vasca)

This is a simple and elegant definitions list:

**First Item**
: The first item in this list is bold

[Second Item](https://en.wikipedia.org/wiki/Tidal_disruption_event)
: The second item links to TDEs


## Execution Statistics
```{nb-exec-table}
```
