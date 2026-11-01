"""
GUI front-end for the c back-end.
"""

import json
import math
import sys
import time
from contextlib import contextmanager

import numpy as np
import matplotlib.ticker

from . import bridge
from . import geez
from . import paths
from . import interface

__all__ = ["run"]


def now_this_is_brother():
    intf = interface.Interface()
    intf.load_json(paths.ROOT / "systems/tenkay.json")

    @contextmanager
    def time_me(label):
        start = time.perf_counter()
        try:
            yield
        finally:
            end = time.perf_counter()
            print(f"{label} took {end - start:.4g} s.")
    print("Simulation starting...")
    with time_me("Simulation"):
        ret = intf.sim()

    print(intf.state)
    if ret is not None:
        print("FAILED:", ret)
        return 1

    intf.view(what="important")
    return 0



def run():
    with geez.instance():
        return now_this_is_brother()

if __name__ == "__main__":
    if len(sys.argv) > 1:
        raise RuntimeError(f"{paths.pyname()} does not have command line args")
    sys.exit(run())
