# da rocket.

import numpy as np

from .. import interface
from .. import geez
from .. import paths

from . import search
from . import units as u
from .system import *



class Rocket(System):
    # --- inputs ---
    V_t: u.ratio
    T0_o: u.degC
    Cd_inj: u.ratio
    etaCstar: u.ratio
    Deq_injo: u.mm
    D_tht: u.mm

    # --- god ---
    @computed(None)
    def deity(self):
        intf = interface.Interface()
        intf.load_json(paths.ROOT / "systems/tenkay.json")
        intf.input_direct("V_t", self.V_t)
        intf.input_direct("vff0_l", float("nan"))
        intf.input_direct("m0_o", (self.V_t/11e-3)
                                  * 9.100114194284641) # 11L 14dC satliq fill.
        intf.input_direct("T0_o", self.T0_o)
        intf.input_direct("Cd_inj", self.Cd_inj)
        intf.input_direct("etaCstar", self.etaCstar)
        intf.input_direct("Deq_injo", self.Deq_injo)
        intf.input_direct("D_tht", self.D_tht)
        intf.sim()
        return intf
    @computed(None)
    def god(self):
        intf = self.deity
        array = intf.snapshots()

        t = array["t"]
        idxof = lambda mask: next(iter(np.flatnonzero(mask)), mask.size)
        Nliq = idxof(array["m_l"] == 0)
        Nburn = idxof(array["Thrust"] != array["Thrust"])

        def otr_velocity():
            diff = array["alt_r"] - intf["alt0_r"] - intf.ROD_HEIGHT
            sign_changes = np.where(np.diff(np.signbit(diff)))[0]
            if sign_changes.size == 0:
                return nan
            idx = sign_changes[0]

            t0, t1 = t[idx], t[idx+1]
            alt0, alt1 = diff[idx], diff[idx+1]
            t_cross = t0 - alt0 * (t1 - t0) / (alt1 - alt0)

            return np.interp(t_cross, t, array["vel_r"])

        # security thru obscurity.
        class God:
            __init__ = lambda s, x: setattr(s, "x", x)
            __getitem__ = lambda s, x: s.x[x]
            __repr__ = lambda s: "??"
        return God(dict(
            apogee=array["alt_r"][-1],
            peak_thrust=array["Thrust"][0],
            impulse=np.trapezoid(array["Thrust"][:Nburn], x=t[:Nburn]),
            min_stiffness=array["Pr_inj"][Nliq - 1] - 1,
            otr_velocity=otr_velocity(),
        ))

    # --- outputs (sim results) ---
    apogee = calc(u.kft, lambda s: s.god["apogee"])
    peak_thrust = calc(u.kN, lambda s: s.god["peak_thrust"])
    impulse = calc(u.kNs, lambda s: s.god["impulse"])
    min_stiffness = calc(u.percent, lambda s: s.god["min_stiffness"])
    otr_velocity = calc(u.m_s, lambda s: s.god["otr_velocity"])



    # --- tolerances ---
    def tolerances(s):
        return dict(
            V_t=0.08,
            T0_o=(s.T0_o, s.T0_o + 2.5),
            Cd_inj=0.18,
            etaCstar=0.12,
        )
    worst_apogee = calc(u.kft,
        lambda s: search.corners(s, "apogee", s.tolerances()).value)
    worst_peak_thrust = calc(u.kN,
        lambda s: search.corners(s, "peak_thrust", s.tolerances()).value)
    worst_impulse = calc(u.kNs,
        lambda s: search.corners(s, "impulse", s.tolerances()).value)
    worst_min_stiffness = calc(u.percent,
        lambda s: search.corners(s, "min_stiffness", s.tolerances()).value)
    worst_otr_velocity = calc(u.m_s,
        lambda s: search.corners(s, "otr_velocity", s.tolerances()).value)


    def plot(self):
        return self.deity.view()




def check_candidate():
    rocket = Rocket(
        V_t=11.0e-3,
        T0_o=273.15 + 14,
        Cd_inj=0.65,
        etaCstar=0.75,
        Deq_injo=10e-3,
        D_tht=36e-3,
    )
    print("===== BASELINE =====")
    rocket.report()
    rocket.plot()

    print("===== SENSITIVITY =====")
    sens = search.sensitivity(rocket, [
        "apogee",
        "peak_thrust",
        "impulse",
        "min_stiffness",
        "otr_velocity",
    ], ["V_t", "T0_o", "Cd_inj", "etaCstar", "Deq_injo", "D_tht"])
    sens.report()
    print()

    print("===== TOLERANCES =====")
    tol = rocket.tolerances()
    for metric in ["apogee", "peak_thrust", "impulse", "min_stiffness",
            "otr_velocity"]:
        print()
        search.corners(rocket, metric, tol, worst="min").report()


    """>> cheeky output dump:
===== BASELINE =====
Rocket
 (- dependent)
    V_t               0.011
    T0_o              14 °C
    Cd_inj            0.65
    etaCstar          0.75
    Deq_injo          10 mm
    D_tht             36 mm
  - deity             <bro.interface.Interface object at 0x0000024255F77A10>
  - god               ??
  - apogee            13.23 kft
  - peak_thrust       4.005 kN
  - impulse           16.25 kNs
  - min_stiffness     33.19 %
  - otr_velocity      25.59 m/s
  - worst_apogee      8.57 kft
  - worst_peak_thrust  3.163 kN
  - worst_impulse     12.77 kNs
  - worst_min_stiffness  17.87 %
  - worst_otr_velocity  22.15 m/s

===== SENSITIVITY =====
Sensitivity (% change in metric per % change in param):
  param     apogee  peak_thrust  impulse  min_stiffness  otr_velocity
  -------------------------------------------------------------------
  V_t       +1.560       +0.000   +1.012         -0.055        -0.092
  T0_o      +0.325      +22.420   +0.017         +4.643        +1.713
  Cd_inj    +0.155       +0.646   +0.050         -1.799        +0.367
  etaCstar  +1.967       +0.679   +1.119         -1.949        +0.386
  Deq_injo  +0.312       +1.291   +0.110         -3.622        +0.733
  D_tht     -0.314       +0.643   -0.224         +3.893        +0.379

===== TOLERANCES =====

Worst-case (min) apogee:
  nominal  13.23 kft
  worst    8.57 kft
  at:
    V_t               0.01012
    T0_o              14 °C
    Cd_inj            0.533
    etaCstar          0.66

Worst-case (min) peak_thrust:
  nominal  4.005 kN
  worst    3.163 kN
  at:
    V_t               0.01012
    T0_o              14 °C
    Cd_inj            0.533
    etaCstar          0.66

Worst-case (min) impulse:
  nominal  16.25 kNs
  worst    12.77 kNs
  at:
    V_t               0.01012
    T0_o              14 °C
    Cd_inj            0.533
    etaCstar          0.66

Worst-case (min) min_stiffness:
  nominal  33.19 %
  worst    17.87 %
  at:
    V_t               0.01188
    T0_o              14 °C
    Cd_inj            0.767
    etaCstar          0.84

Worst-case (min) otr_velocity:
  nominal  25.59 m/s
  worst    22.15 m/s
  at:
    V_t               0.01188
    T0_o              14 °C
    Cd_inj            0.533
    etaCstar          0.66

    """



def search_space():
    import matplotlib.pyplot as plt

    rocket = Rocket(
        T0_o=273.15 + 14,
        Cd_inj=0.65,
        Deq_injo=14.0e-3,
        etaCstar=0.75,
        D_tht=36e-3,
    )
    all_metrics = [
        "apogee",
        "worst_apogee",
        "peak_thrust",
        "worst_peak_thrust",
        "impulse",
        "worst_impulse",
        "min_stiffness",
        "worst_min_stiffness",
        "otr_velocity",
        "worst_otr_velocity",
    ]
    metrics = [
        "min_stiffness",
        "apogee",
        "otr_velocity",
    ]
    metrics_plusworst = [
        "min_stiffness",
        "worst_min_stiffness",
        "apogee",
        "worst_apogee",
        "otr_velocity",
        "worst_otr_velocity",
    ]


    # --- 0. nominal baseline ---
    rocket.report()

    # --- 1. sensitivity: which inputs matter, and for what? ---
    inputs = ["V_t", "T0_o", "Cd_inj", "etaCstar", "Deq_injo", "D_tht"]
    sens = search.sensitivity(rocket, all_metrics, inputs)
    sens.report()
    for m in metrics:
        print(f"\nmost influential on {m}:")
        for p, e in sens.ranked(m):
            print(f"  {p:<10}{e:+.3f}")

    # --- 2. 1D sweeps: throat size and injector area ---
    sw_tht = search.sweep(rocket, "D_tht", np.linspace(28e-3, 44e-3, 17),
            *metrics)
    sw_tht.report()
    sw_tht.plot()

    sw_T = search.sweep(rocket, "Deq_injo", np.linspace(9e-3, 16e-3, 15),
            *metrics)
    sw_T.report()
    sw_T.plot()

    # --- 3. 2D grid: injector vs throat, where does stiffness run out? ---
    g = search.grid(rocket,
            "D_tht", np.linspace(28e-3, 44e-3, 9),
            "Deq_injo", np.linspace(9e-3, 16e-3, 9),
            "min_stiffness")
    g.plot(contours=(0.25,))

    # --- 5. tolerance stack-up: worst case of every output ---
    # (2^4 = 16 sims per metric.) `worst` is the bad direction for each.
    tol = rocket.tolerances()
    for metric in metrics:
        print()
        search.corners(rocket, metric, tol, worst="min").report()

    # --- 6. optimise: highest apogee that stays safe and stable ---
    opt = search.optimise(rocket, "apogee",
            params={"D_tht": (28e-3, 44e-3),
                    "Deq_injo": (9e-3, 16e-3)},
            constraints={"min_stiffness": (0.3, None),
                         "otr_velocity": (28.0, None)},
            maximise=True, seed=0, maxiter=20, popsize=8)
    opt.report()

    # does the optimum survive tolerances?
    print()
    tol = opt.system.tolerances()
    search.corners(opt.system, "min_stiffness", tol, worst="min").report()
    search.corners(opt.system, "otr_velocity", tol, worst="min").report()

    plt.show()


def run():
    with geez.instance():
        # search_space()
        check_candidate()
