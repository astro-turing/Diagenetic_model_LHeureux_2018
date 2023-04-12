# ~\~ language=Python filename=marlpde/marlpde.py
# ~\~ begin <<lit/python-interface.md|marlpde/marlpde.py>>[init]
# import numpy
import configparser
from pathlib import Path
from dataclasses import (dataclass, asdict)
from subprocess import (run)
import h5py as h5
from pint import UnitRegistry


u = UnitRegistry()
quantity = u.Quantity

@dataclass
class Scenario:
    mua: quantity    = 100.09 * u.g/u.mol
    rhoa: quantity   = 2.95 * u.g/u.cm**3
    rhoc: quantity   = 2.71 * u.g/u.cm**3
    rhot: quantity   = 2.8 * u.g/u.cm**3
    rhow: quantity   = 1.023 * u.g/u.cm**3
    D0ca: quantity   = 131.9 * u.cm**2/u.a
    D0co3: quantity  = 272.6 * u.cm**2/u.a
    Ka: quantity     = 10**(-6.19) * u.M**2
    Kc: quantity     = 10**(-6.37) * u.M**2
    beta: quantity   = 0.1 * u.cm / u.a
    b: quantity      = 10.0 / u.kPa
    k1: quantity     = 1.0 / u.a
    k2: quantity     = 1.0 / u.a
    k3: quantity     = 0.1 / u.a
    k4: quantity     = 0.1 / u.a
    nn: quantity     = 2.8 * u.dimensionless
    m: quantity      = 2.48 * u.dimensionless
    S: quantity      = 0.1 * u.cm / u.a
    # cAthy: quantity  = 0.1 * u.dimensionless
    phiinf: quantity = 0.01 * u.dimensionless
    phi0: quantity   = 0.8 * u.dimensionless
    phiin: quantity  = 0.8 * u.dimensionless
    ca0: quantity    = 0.326e-3 * u.M
    co30: quantity   = 0.326e-3 * u.M
    ccal0: quantity  = 0.3 * u.dimensionless
    cara0: quantity  = 0.6 * u.dimensionless
    xdis: quantity   = 50.0 * u.cm       # x_d   (start of dissolution zone)
    xcem: quantity   = -100.0 * u.cm
    xcemf: quantity  = 1000.0 * u.cm
    length: quantity = 500.0 * u.cm
    Th: quantity     = 100.0 * u.cm      # h_d   (height of dissolution zone)
    phi00: quantity  = 0.8 * u.dimensionless
    ca00: quantity   = 0.326e-3 * u.M    # sqrt(Kc) / 2
    co300: quantity  = 0.326e-3 * u.M    # sqrt(Kc) / 2
    ccal00: quantity = 0.3 * u.dimensionless
    cara00: quantity = 0.6 * u.dimensionless

@dataclass
class Solver:
    dt: float     = 1.e-5
    eps: float    = 1.e-2
    tmax: int     = 100_000
    outt: int     =   1_000      # timesteps inbetween writing
    outx: int     =  25_000
    N: int        = 200

def cAthy(s: Scenario):
    return ((1 - s.phi0) * s.b * 9.81 * u['m/s²'] * s.rhow).to('cm⁻¹')

def write_input_cfg(path: Path, solver: Solver, scenario: Scenario):
    cfg = configparser.ConfigParser()
    cfg.optionxform = str
    cfg["Solver"] = asdict(solver)

    units = { k: v.units for (k, v) in asdict(Scenario()).items() }
    magnitudes = { k: v.to(units[k]).magnitude
                   for (k, v) in asdict(scenario).items() }
    magnitudes["cAthy"] = cAthy(scenario).magnitude
    cfg["Scenario"] = magnitudes
    path.mkdir(parents=True, exist_ok=True)
    with open(path / "input.cfg", "w") as f_cfg:
        cfg.write(f_cfg)

def run_marl_pde(path: Path, exe_dir: Path = Path(".")):
    run(exe_dir / "marl-pde", cwd=path, check=True)

def output_data(path: Path):
    return h5.File(path / "output.h5", mode="r")
# ~\~ end
