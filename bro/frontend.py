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

__all__ = ["run"]


def get_interpretation():
    interp = bridge.Interpretation()
    IN = interp.INPUT
    OUT = interp.OUTPUT

    interp.append("x", interp.F64, IN)

    interp.finalise()
    return interp

def get_state(interp):
    state = bridge.State(interp)

    state["x"] = 0.8

    return state



def now_this_is_brother():
    interp = get_interpretation()
    state = get_state(interp)
    print("python")
    print(state)


    @contextmanager
    def time_me(label):
        start = time.perf_counter()
        try:
            yield
        finally:
            end = time.perf_counter()
            print(f"{label} took {end - start:.4g} s.")
    with time_me("Simulation"):
        print("executing")
        ret = state.execute()
    if ret is not None:
        print("FAILED:", ret)
        return 1

    print("end")
    print(state)
    return 0


def run():
    with geez.instance():
        return now_this_is_brother()


if __name__ == "__main__":
    if len(sys.argv) > 1:
        raise RuntimeError("bruv.frontend does not have command line args")
    sys.exit(run())
