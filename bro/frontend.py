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
    INDATA = interp.INPUT_DATA
    OUTDATA = interp.OUTPUT_DATA

    interp.append("static_fire", interp.I64, IN)

    interp.append("alt0_r", interp.F64, IN)
    interp.append("D_r", interp.F64, IN)
    interp.append("A_r", interp.F64, OUT)
    interp.append("mdry_r", interp.F64, IN)
    interp.append("mwet0_r", interp.F64, OUT)

    interp.append("V_t", interp.F64, IN)
    interp.append("C_tw", interp.F64, IN)

    interp.append("T0_o", interp.F64, IN)
    interp.append("m0_o", interp.F64, IN | OUT)
    interp.append("vff0_l", interp.F64, IN | OUT)

    interp.append("model_inj", interp.I64, IN)
    interp.append("Cd_inj", interp.F64, IN)
    interp.append("no_injo", interp.F64, IN | OUT)
    interp.append("D_injo", interp.F64, IN | OUT)
    interp.append("Deq_injo", interp.F64, IN | OUT)
    interp.append("A_injo", interp.F64, OUT)
    interp.append("CdA_inj", interp.F64, OUT)

    interp.append("etaCstar", interp.F64, IN)
    interp.append("etaCf", interp.F64, IN)
    interp.append("Vempty_c", interp.F64, IN)

    interp.append("L0_f", interp.F64, IN)
    interp.append("th0_f", interp.F64, IN)
    interp.append("OD0_f", interp.F64, IN)
    interp.append("rho_f", interp.F64, IN | OUT)
    interp.append("m0_f", interp.F64, IN | OUT)

    interp.append("AEAT", interp.F64, IN)
    interp.append("NLF", interp.F64, IN)
    interp.append("D_tht", interp.F64, IN)
    interp.append("A_tht", interp.F64, OUT)
    interp.append("D_exit", interp.F64, OUT)
    interp.append("A_exit", interp.F64, OUT)
    interp.append("theta_div", interp.F64, OUT)
    interp.append("theta_exit", interp.F64, OUT)

    interp.append("out_upto", interp.I64, IN | OUT)
    interp.append("out_count", interp.I64, IN)
    interp.append("out_t", interp.PTR_F64, IN | OUTDATA)
    interp.append("out_T_o", interp.PTR_F64, IN | OUTDATA)
    interp.append("out_m_o", interp.PTR_F64, IN | OUTDATA)
    interp.append("out_rr_f", interp.PTR_F64, IN | OUTDATA)
    interp.append("out_alt_r", interp.PTR_F64, IN | OUTDATA)
    interp.append("out_vel_r", interp.PTR_F64, IN | OUTDATA)

    interp.finalise()
    return interp

def get_state(interp):
    nan = float("nan")
    state = bridge.State(interp)

    state["static_fire"] = 0

    state["alt0_r"] = 0.0
    state["D_r"] = 154e-3
    state["mdry_r"] = 42.5

    state["V_t"] = 11e-3
    state["C_tw"] = 9.059*896/8

    state["T0_o"] = 273.15 + 14.0
    state["m0_o"] = nan
    state["vff0_l"] = 1.0

    state["model_inj"] = 3
    state["Cd_inj"] = 0.65
    state["no_injo"] = 85.0
    state["D_injo"] = nan
    state["Deq_injo"] = 14.0e-3

    state["etaCstar"] = 0.75
    state["etaCf"] = 0.97
    state["Vempty_c"] = 3.61e-3

    state["L0_f"] = 385e-3
    state["th0_f"] = 18e-3
    state["OD0_f"] = 95e-3
    state["rho_f"] = nan
    state["m0_f"] = 1.47035

    state["AEAT"] = 3.0
    state["NLF"] = 1.0
    state["D_tht"] = 40e-3

    state["out_upto"] = 0
    state["out_count"] = 1_000_000
    new_out = lambda: np.empty(shape=(state["out_count"],), dtype=np.float64)
    state["out_t"] = new_out()
    state["out_T_o"] = new_out()
    state["out_m_o"] = new_out()
    state["out_rr_f"] = new_out()
    state["out_alt_r"] = new_out()
    state["out_vel_r"] = new_out()

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
    print("Simulation starting...")
    with time_me("Simulation"):
        ret = state.execute()
    print()
    print(state)
    if ret is not None:
        print("FAILED:", ret)
        return 1




    import matplotlib.pyplot as plt

    get_out = lambda s: state[f"out_{s}"].view(state["out_count"])
    n = int(state["out_upto"])
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
