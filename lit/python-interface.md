---
title: Python interface
author: Johan Hidding
---

We extended the original code with functions to read a `input.cfg` so that we can run multiple instances from Python.

``` {.python file=marlpde/__init__.py}

```

``` {.python file=marlpde/__main__.py}
import argh
from pathlib import Path

from .marlpde import (write_input_cfg, Solver, Scenario)

def main(path: Path = Path(".")):
    write_input_cfg(path, Solver(), Scenario())

if __name__ == "__main__":
    argh.dispatch_command(main)
```

``` {.python file=marlpde/marlpde.py}
# import numpy
import configparser
from pathlib import Path
from dataclasses import (dataclass, asdict)
from subprocess import (run)
import h5py as h5

@dataclass
class Scenario:
    mua: float    = 100.09
    rhoa: float   = 2.95
    rhoc: float   = 2.71
    rhot: float   = 2.8
    rhow: float   = 1.023
    D0ca: float   = 131.9
    D0co3: float  = 272.6
    Ka: float     = 6.4565e-7
    Kc: float     = 4.2658e-7
    beta: float   = 0.01
    k1: float     = 0.0
    k2: float     = 0.01
    k3: float     = 0.001
    nn: float     = 2.8
    m: float      = 2.48
    S: float      = 0.01
    cAthy: float  = 0.1
    phiinf: float = 0.01
    phi0: float   = 0.7
    ca0: float    = 6.5313e-4
    co30: float   = 6.5313e-4
    ccal0: float  = 0.3
    cara0: float  = 0.6

@dataclass
class Solver:
    dt: float     = 5.e-6
    xdis: float   = 50.0
    xcem: float   = -100.0
    xcemf: float  = 1000.0
    length: float = 500.0
    eps: float    = 1.e-2
    Th: float     = 100.0
    tmax: int     = 20_000
    outt: int     =  1_000
    outx: int     = 50_000
    N: int        = 200

def write_input_cfg(path: Path, solver: Solver, scenario: Scenario):
    cfg = configparser.ConfigParser()
    cfg.optionxform = str
    cfg["Solver"] = asdict(solver)
    cfg["Scenario"] = asdict(scenario)
    path.mkdir(parents=True, exist_ok=True)
    with open(path / "input.cfg", "w") as f_cfg:
        cfg.write(f_cfg)

def run_marl_pde(path: Path, exe_dir: Path = Path(".")):
    run(exe_dir / "marl-pde", cwd=path, check=True)

def output_data(path: Path):
    return h5.File(path / "output.h5", mode="r")
```
