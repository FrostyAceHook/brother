# da rocket.

"""
Push-fit pneumatically-actuator quick-disconnect.
"""

import numpy as np

from .. import interface
from .. import geez
from .. import paths

from . import search
from . import units as u
from .system import *



class Rocket(System):
    # --- inputs ---
    T0_o: u.K
    vff0_l: u.ratio
    Cd_inj: u.ratio
    Deq_injo: u.mm
    etaCstar: u.ratio
    D_tht: u.mm

    # --- god ---
    @computed(None)
    def deity(self):
        intf = interface.Interface()
        intf.load_json(paths.ROOT / "systems/tenkay.json")
        intf.input_direct("T0_o", self.T0_o)
        intf.input_direct("vff0_l", self.vff0_l)
        intf.input_direct("Cd_inj", self.Cd_inj)
        intf.input_direct("Deq_injo", self.Deq_injo)
        intf.input_direct("etaCstar", self.etaCstar)
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

        return dict(
            apogee=array["alt_r"][-1],
            peak_thrust=array["Thrust"][0],
            impulse=np.trapezoid(array["Thrust"][:Nburn], x=t[:Nburn]),
            max_stiffness=array["Pr_inj"][0] - 1,
            min_stiffness=array["Pr_inj"][Nliq - 1] - 1,
            otr_velocity=otr_velocity(),
        )

    # --- outputs (sim results) ---
    apogee = calc(u.m, lambda s: s.god["apogee"])
    peak_thrust = calc(u.N, lambda s: s.god["peak_thrust"])
    impulse = calc(u.Ns, lambda s: s.god["impulse"])
    max_stiffness = calc(u.percent, lambda s: s.god["max_stiffness"])
    min_stiffness = calc(u.percent, lambda s: s.god["min_stiffness"])
    otr_velocity = calc(u.m_s, lambda s: s.god["otr_velocity"])



    # --- tolerances ---
    def tolerances(s):
        return dict(
            T0_o=(s.T0_o - 5, s.T0_o + 5),
            vff0_l=(s.vff0_l - 0.2, s.vff0_l + 0.0),
            Cd_inj=0.1,
            etaCstar=0.1,
        )
    worst_apogee = calc(u.m,
        lambda s: search.corners(s, "apogee", s.tolerances()).value)
    worst_peak_thrust = calc(u.N,
        lambda s: search.corners(s, "peak_thrust", s.tolerances()).value)
    worst_impulse = calc(u.Ns,
        lambda s: search.corners(s, "impulse", s.tolerances()).value)
    worst_max_stiffness = calc(u.percent,
        lambda s: search.corners(s, "max_stiffness", s.tolerances()).value)
    worst_min_stiffness = calc(u.percent,
        lambda s: search.corners(s, "min_stiffness", s.tolerances()).value)
    worst_otr_velocity = calc(u.m_s,
        lambda s: search.corners(s, "otr_velocity", s.tolerances()).value)


    def plot(self):
        return self.deity.view()



def main():
    rocket = Rocket(
        T0_o=273.15 + 14,
        vff0_l=1,
        Cd_inj=0.65,
        Deq_injo=11.0e-3,
        etaCstar=0.75,
        D_tht=36e-3,
    )
    rocket.report()
    rocket.plot()
