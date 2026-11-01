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

    print("python")
    print(intf.state)


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
        ret = intf.state.execute()
    print()
    print(intf.state)
    if ret is not None:
        print("FAILED:", ret)
        return 1




    import matplotlib.pyplot as plt

    get_out = lambda s: intf.state[f"out_{s}"].view(intf.state["out_count"])
    n = int(intf.state["out_upto"])
    get = lambda k: get_out(k)[:n]
    t = get("t")

    names = ["T_o", "m_o", "rr_f", "alt_r", "vel_r"]
    fig, axs = plt.subplots(2, 3, figsize=(14, 12), sharex=True)
    axs = axs.ravel()
    for ax, k in zip(axs, names):
        ax.plot(t, get(k))
        ax.set_title(k)
        ax.grid(alpha=0.3)

    for ax in axs[-3:]:
        ax.set_xlabel("t [s]")
    fig.tight_layout()
    plt.show()


    print("end")
    return 0




def run():
    with geez.instance():
        return now_this_is_brother()

if __name__ == "__main__":
    if len(sys.argv) > 1:
        raise RuntimeError(f"{paths.pyname()} does not have command line args")
    sys.exit(run())
