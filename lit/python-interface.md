---
title: Python interface
author: Johan Hidding
---

We extended the original code with functions to read a `input.cfg` so that we can run multiple instances from Python.

``` {.python file=marlpde/__init__.py}

```

``` {.python file=marlpde/marlpde.py}
import numpy

@dataclass
class Scenario:
    pass

@dataclass
class Solver:
    pass


```
