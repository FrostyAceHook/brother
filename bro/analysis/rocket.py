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
                                  * 8.799545319460602) # 11L 18dC satliq fill.
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
    impulse = calc(u.kNs, lambda s: s.god["impulse"])
    peak_thrust = calc(u.kN, lambda s: s.god["peak_thrust"])
    otr_velocity = calc(u.m_s, lambda s: s.god["otr_velocity"])
    min_stiffness = calc(u.percent, lambda s: s.god["min_stiffness"])

    # --- tolerances ---
    def tolerances(s):
        return dict(
            V_t=0.1,
            T0_o=(s.T0_o - 2.5, s.T0_o + 2.5),
            Cd_inj=0.18,
            etaCstar=0.12,
        )
    worst_apogee = calc(u.kft,
        lambda s: search.corners(s, "apogee", s.tolerances()).value)
    worst_impulse = calc(u.kNs,
        lambda s: search.corners(s, "impulse", s.tolerances()).value)
    worst_peak_thrust = calc(u.kN,
        lambda s: search.corners(s, "peak_thrust", s.tolerances(), worst="max")
                        .value)
    worst_otr_velocity = calc(u.m_s,
        lambda s: search.corners(s, "otr_velocity", s.tolerances()).value)
    worst_min_stiffness = calc(u.percent,
        lambda s: search.corners(s, "min_stiffness", s.tolerances()).value)


    def plot(self):
        return self.deity.view()




def check_candidate():
    rocket = Rocket(
        V_t=11.0e-3,
        T0_o=273.15 + 18,
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
        "impulse",
        "peak_thrust",
        "otr_velocity",
        "min_stiffness",
    ], ["V_t", "T0_o", "Cd_inj", "etaCstar", "Deq_injo", "D_tht"])
    sens.report()
    print()

    print("===== TOLERANCES =====")
    tol = rocket.tolerances()
    for metric in ["apogee", "impulse", "peak_thrust", "otr_velocity",
            "min_stiffness"]:
        print()
        search.corners(rocket, metric, tol,
                worst="max" if metric == "peak_thrust" else "min").report()


    """>> cheeky output dump:
===== BASELINE =====
Rocket
 (- dependent)
    V_t               0.011
    T0_o              18 °C
    Cd_inj            0.65
    etaCstar          0.75
    Deq_injo          10 mm
    D_tht             36 mm
  - deity             <bro.interface.Interface object at 0x00000209F34378C0>
  - god               ??
  - apogee            12.2 kft
  - impulse           15.67 kNs
  - peak_thrust       4.211 kN
  - otr_velocity      26.43 m/s
  - min_stiffness     36.76 %
  - worst_apogee      7.657 kft
  - worst_impulse     12.04 kNs
  - worst_peak_thrust  6.523 kN
  - worst_otr_velocity  22.4 m/s
  - worst_min_stiffness  18.81 %

===== SENSITIVITY =====
Sensitivity (% change in metric per % change in param):
  param     apogee  impulse  peak_thrust  otr_velocity  min_stiffness
  -------------------------------------------------------------------
  V_t       +1.530   +1.005       -0.000        -0.089         -0.037
  T0_o      +0.213   -0.029      +20.423        +1.489         +4.621
  Cd_inj    +0.149   +0.050       +0.674        +0.379         -1.761
  etaCstar  +1.931   +1.121       +0.718        +0.404         -1.910
  Deq_injo  +0.296   +0.109       +1.349        +0.759         -3.546
  D_tht     -0.336   -0.232       +0.565        +0.334         +3.818

===== TOLERANCES =====

Worst-case (min) apogee:
  nominal  12.2 kft
  worst    7.657 kft
  at:
    V_t               0.0099
    T0_o              15.5 °C
    Cd_inj            0.533
    etaCstar          0.66

Worst-case (min) impulse:
  nominal  15.67 kNs
  worst    12.04 kNs
  at:
    V_t               0.0099
    T0_o              20.5 °C
    Cd_inj            0.533
    etaCstar          0.66

Worst-case (max) peak_thrust:
  nominal  4.211 kN
  worst    6.523 kN
  at:
    V_t               0.0121
    T0_o              20.5 °C
    Cd_inj            0.767
    etaCstar          0.84

Worst-case (min) otr_velocity:
  nominal  26.43 m/s
  worst    22.4 m/s
  at:
    V_t               0.0121
    T0_o              15.5 °C
    Cd_inj            0.533
    etaCstar          0.66

Worst-case (min) min_stiffness:
  nominal  36.76 %
  worst    18.81 %
  at:
    V_t               0.0121
    T0_o              15.5 °C
    Cd_inj            0.767
    etaCstar          0.84
    """



def search_space():
    import matplotlib.pyplot as plt

    rocket = Rocket(
        V_t=11.0e-3,
        T0_o=273.15 + 18,
        Cd_inj=0.65,
        etaCstar=0.75,
        Deq_injo=14e-3,
        D_tht=36e-3,
    )
    all_metrics = [
        "apogee",
        "worst_apogee",
        "impulse",
        "worst_impulse",
        "peak_thrust",
        "worst_peak_thrust",
        "otr_velocity",
        "worst_otr_velocity",
        "min_stiffness",
        "worst_min_stiffness",
    ]
    metrics = [
        "apogee",
        "otr_velocity",
        "min_stiffness",
    ]
    metrics_plusworst = [
        "apogee",
        "worst_apogee",
        "otr_velocity",
        "worst_otr_velocity",
        "min_stiffness",
        "worst_min_stiffness",
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
            constraints={"otr_velocity": (25.0, None),
                         "min_stiffness": (0.3, None)},
            maximise=True, seed=0, maxiter=20, popsize=8)
    opt.report()

    # does the optimum survive tolerances?
    print()
    tol = opt.system.tolerances()
    search.corners(opt.system, "otr_velocity", tol, worst="min").report()
    search.corners(opt.system, "min_stiffness", tol, worst="min").report()

    plt.show()


def run():
    with geez.instance():
        # search_space()
        check_candidate()
