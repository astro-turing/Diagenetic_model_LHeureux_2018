# ~\~ language=Python filename=marlpde/__main__.py
# ~\~ begin <<lit/python-interface.md|marlpde/__main__.py>>[init]
import argh
from pathlib import Path

from .marlpde import (write_input_cfg, Solver, Scenario)

def main(path: Path = Path(".")):
    write_input_cfg(path, Solver(), Scenario())

if __name__ == "__main__":
    argh.dispatch_command(main)
# ~\~ end
