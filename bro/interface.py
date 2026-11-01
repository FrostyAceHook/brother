"""
Interface object which understands the state struct and interacts directly with
the c back-end.
"""

import ast
import json
import functools
import math
from dataclasses import dataclass
from math import pi, nan, inf

import numpy as np
import matplotlib

from . import bridge
from . import geez



### DESCRIPTORS #################################################################

class Caster:
    def __init__(self, display=None):
        self.display = display
    def to_value(self, string):
        raise NotImplementedError()
    def to_string(self, value):
        if value is None:
            return ""
        return str(value)

    _ALLOWED_NAMES = {
        # Constants.
        "pi": math.pi,
        "e": math.e,
        # Functions.
        "pow": pow,
        "sqrt": math.sqrt,
        "cbrt": math.cbrt,
        # Trig.
        "torad": math.radians,
        "todeg": math.degrees,
        "sin": math.sin,
        "cos": math.cos,
        "tan": math.tan,
        "asin": math.asin,
        "acos": math.acos,
        "atan": math.atan,
        "sind": lambda x: math.sin(math.radians(x)),
        "cosd": lambda x: math.cos(math.radians(x)),
        "tand": lambda x: math.tan(math.radians(x)),
        "asind": lambda x: math.degrees(math.asin(x)),
        "acosd": lambda x: math.degrees(math.acos(x)),
        "atand": lambda x: math.degrees(math.atan(x)),
        # Lists.
        "len": len,
        "min": min,
        "max": max,
        "sum": sum,
        "prod": math.prod,
    }
    def atof(self, string):
        if not string:
            return nan
        # Use ast parsing to allow general math expressions.
        node = ast.parse(string, mode="eval")
        # Validate each node.
        for subnode in ast.walk(node):
            # Type must be legal.
            if isinstance(subnode, (ast.Assign, ast.AnnAssign, ast.AugAssign,
                    ast.NamedExpr)):
                raise ValueError("illegal assigning node "
                                f"{repr(type(subnode).__name__)}")
            # Name must be legal.
            if isinstance(subnode, ast.Name):
                if subnode.id not in self._ALLOWED_NAMES:
                    raise ValueError(f"illegal name {repr(subnode.id)}")
        # Worlds safest eval (trust).
        compiled = compile(node, "<string>", "eval")
        v = eval(compiled, {"__builtins__": {}}, self._ALLOWED_NAMES)
        if not isinstance(v, (bool, int, float)):
            raise ValueError(f"must be float, got {repr(type(v).__name__)}")
        v = float(v)
        if v != v:
            raise ValueError("must be nonnan")
        if abs(v) == inf:
            raise ValueError("must be noninf")
        return v
    def ftoa(self, value):
        if value is None:
            return ""
        string = f"{value:.8g}"
        if string.endswith(".0"):
            string = string[:-2]
        return string

class Dimensioned(Caster):
    def __init__(self, display, scale=1.0, offset=0.0):
        super().__init__(display)
        self.scale = scale
        self.offset = offset
    def to_value(self, string):
        return self.atof(string) * self.scale + self.offset
    def to_string(self, value):
        return self.ftoa(None if value is None else \
                         (value - self.offset) / self.scale)
Dimensioned.none = Dimensioned("-")
Dimensioned.count = Dimensioned("no.")
Dimensioned.percent = Dimensioned("%", scale=1e-2)
Dimensioned.deg = Dimensioned("°", scale=pi/180.0)
Dimensioned.m = Dimensioned("m")
Dimensioned.mm = Dimensioned("mm", scale=1e-3)
Dimensioned.ft = Dimensioned("ft", scale=1.0 / 3.281)
Dimensioned.kft = Dimensioned("kft", scale=1e3 / 3.281)
Dimensioned.mm2 = Dimensioned("mm²", scale=1e-6)
Dimensioned.L = Dimensioned("L", scale=1e-3)
Dimensioned.s = Dimensioned("s")
Dimensioned.m_s = Dimensioned("m/s")
Dimensioned.m_s2 = Dimensioned("m/s²")
Dimensioned.N = Dimensioned("N")
Dimensioned.Pa = Dimensioned("Pa")
Dimensioned.MPa = Dimensioned("MPa", scale=1e6)
Dimensioned.bar = Dimensioned("bar", scale=1e5)
Dimensioned.K = Dimensioned("K")
Dimensioned.dC = Dimensioned("°C", offset=273.15)
Dimensioned.kg = Dimensioned("kg")
Dimensioned.kg_m3 = Dimensioned("kg/m³")
Dimensioned.kg_s = Dimensioned("kg/s")
Dimensioned.J_kg_K = Dimensioned("J/kg/K")
Dimensioned.J_K = Dimensioned("J/K")
Dimensioned.J_kg = Dimensioned("J/kg")
Dimensioned.g_mol = Dimensioned("g/mol", scale=1e3)


class Bool(Caster):
    def __init__(self):
        super().__init__("(y/n)")
    def to_value(self, string):
        if string.strip().casefold() in {"y", "yes", "yeah", "yeagh", "true"}:
            return True
        if string.strip().casefold() in {"n", "no", "nah", "nogh", "false"}:
            return False
        raise ValueError(f"excepted 'true' or 'false', got {repr(string)}")
    def to_string(self, value):
        if value is None:
            return ""
        return "ny"[not not value]



@dataclass(frozen=True, kw_only=True)
class Parameter:
    c_name: str
    c_type: int
    c_flags: int
    display_name: str = None
    caster: Caster = None
    in_snapshot: bool = False

    @classmethod
    def input(cls, c_name, display_name, caster, **kw):
        opt = dict(display_name=display_name, caster=caster,
                   c_name=c_name, c_type=bridge.Interpretation.F64,
                   c_flags=bridge.Interpretation.INPUT)
        opt.update(kw)
        return cls(**opt)
    @classmethod
    def output(cls, c_name, display_name, caster, **kw):
        opt = dict(display_name=display_name, caster=caster,
                   c_name=c_name, c_type=bridge.Interpretation.F64,
                   c_flags=bridge.Interpretation.OUTPUT)
        opt.update(kw)
        return cls(**opt)
    @classmethod
    def silent(cls, c_name, **kw):
        opt = dict(c_name=c_name, c_type=bridge.Interpretation.F64,
                   c_flags=bridge.Interpretation.OUTPUT)
        opt.update(kw)
        return cls(**opt)
    @classmethod
    def either(cls, c_name, display_name, caster, **kw):
        opt = dict(display_name=display_name, caster=caster,
                   c_name=c_name, c_type=bridge.Interpretation.F64,
                   c_flags=bridge.Interpretation.INPUT
                         | bridge.Interpretation.OUTPUT)
        opt.update(kw)
        return cls(**opt)

    @classmethod
    def snapshotted(cls, c_name, display_name, dimensioned, **kw):
        opt = dict(display_name=display_name, caster=dimensioned,
                   in_snapshot=True,
                   c_name=c_name, c_type=bridge.Interpretation.F64,
                   c_flags=bridge.Interpretation.OUTPUT)
        assert isinstance(dimensioned, Dimensioned)
        opt.update(kw)
        return cls(**opt)

    @property
    def is_displayed(self):
        return self.display_name is not None
    @property
    def is_input(self):
        return self.is_displayed and \
               bool(self.c_flags & bridge.Interpretation.INPUT)
    @property
    def is_output(self):
        return self.is_displayed and \
               bool(self.c_flags & bridge.Interpretation.OUTPUT)



### STATE INTERPRETATION ########################################################

@functools.cache
def params_grouped():
    inp = Parameter.input
    out = Parameter.output
    sil = Parameter.silent
    eit = Parameter.either
    snap = Parameter.snapshotted

    F64 = bridge.Interpretation.F64
    I64 = bridge.Interpretation.I64
    PTR_F64 = bridge.Interpretation.PTR_F64

    IN  = bridge.Interpretation.INPUT
    OUT = bridge.Interpretation.OUTPUT
    INDATA  = bridge.Interpretation.INPUT_DATA
    OUTDATA = bridge.Interpretation.OUTPUT_DATA

    return {
        "Setup": [
            inp("static_fire", "Static fire?", Bool()),
            inp("alt0_r", "Initial altitude", Dimensioned.m),
        ],
        "Rocket": [
            inp("D_r", "Outer diameter", Dimensioned.mm),
            sil("A_r"),
            inp("mdry_r", "Dry mass", Dimensioned.kg),
            out("mwet0_r", "Wet mass", Dimensioned.kg),
        ],
        "Oxidiser Tank": [
            inp("V_t", "Volume", Dimensioned.L),
            inp("C_tw", "Heat capacity", Dimensioned.J_K),
        ],
        "Oxidiser Pool": [
            inp("T0_o", "Initial temperature", Dimensioned.dC),
            eit("m0_o", "Initial mass", Dimensioned.kg),
            eit("vff0_l", "Initial liq. vol. fraction", Dimensioned.none),
        ],
        "Injector": [
            inp("Cd_inj", "Discharge coefficient", Dimensioned.none),
            eit("no_injo", "Orifice count", Dimensioned.count),
            eit("D_injo", "Orifice diameter", Dimensioned.mm),
            eit("Deq_injo", "Orifice eq. single diameter", Dimensioned.mm),
            out("A_injo", "Orifice net area", Dimensioned.mm2),
            sil("CdA_inj"),
        ],
        "Combustion Chamber": [
            inp("etaCstar", "C* efficiency", Dimensioned.percent),
            inp("etaCf", "Cf efficiency", Dimensioned.percent),
            inp("Vempty_c", "Empty volume", Dimensioned.L),
        ],
        "Fuel Grain": [
            inp("L0_f", "Length", Dimensioned.mm),
            inp("OD0_f", "Outer diameter", Dimensioned.mm),
            inp("th0_f", "Thickness", Dimensioned.mm),
            eit("rho_f", "Density", Dimensioned.kg_m3),
            eit("m0_f", "Mass", Dimensioned.kg),
        ],
        "Nozzle": [
            eit("D_tht", "Throat diameter", Dimensioned.mm),
            eit("A_tht", "Throat area", Dimensioned.mm2),
            eit("D_exit", "Exit diameter", Dimensioned.mm),
            eit("A_exit", "Exit area", Dimensioned.mm2),
            eit("AEAT", "Exit-throat area ratio", Dimensioned.none),
            inp("NLF", "Length fraction", Dimensioned.percent),
            out("theta_div", "Diverging half-angle", Dimensioned.deg),
            out("theta_exit", "Exit half-angle", Dimensioned.deg),
        ],
        "__snapshot__": [
            sil("snap_upto", c_type=I64, c_flags=IN | OUT),
            sil("snap_count", c_type=I64, c_flags=IN),
            sil("snap_hash", c_type=I64, c_flags=IN),
            sil("snap_array", c_type=PTR_F64, c_flags=IN | OUTDATA),
        ],
        "__snapshot_interp__": [
            snap("t", None,
                Dimensioned.s),
            snap("alt_r", "Altitude",
                Dimensioned.kft),
            snap("vel_r", "Velocity",
                Dimensioned.m_s),
            snap("acc_r", "Acceleration",
                Dimensioned.m_s2),
            snap("mwet_r", "Rocket wet mass",
                Dimensioned.kg),
            snap("P_t", "Tank pressure",
                Dimensioned.bar),
            snap("P_c", "Chamber pressure",
                Dimensioned.bar),
            snap("T_o", "Oxidiser temperature",
                Dimensioned.dC),
            snap("m_o", "Oxidiser mass",
                Dimensioned.kg),
            snap("m_l", "Liquid oxidiser mass",
                Dimensioned.kg),
            snap("m_v", "Vapour oxidiser mass",
                Dimensioned.kg),
            snap("rho_l", "Liquid oxidiser density",
                Dimensioned.kg_m3),
            snap("rho_v", "Vapour oxidiser density",
                Dimensioned.kg_m3),
            snap("vff_l", "Tank liquid volume fraction",
                Dimensioned.none),
            snap("X_o", "Oxidiser vapour quality",
                Dimensioned.none),
            snap("s_o", "Oxidiser specific entropy",
                Dimensioned.J_kg_K),
            snap("h_o", "Oxidiser specific enthalpy",
                Dimensioned.J_kg),
            snap("u_o", "Oxidiser specific internal energy",
                Dimensioned.J_kg),
            snap("cp_o", "Oxidiser isobaric specific heat capacity",
                Dimensioned.J_kg_K),
            snap("cv_o", "Oxidiser isochoric specific heat capacity",
                Dimensioned.J_kg_K),
            snap("Y_o", "Oxidiser specific heat ratio",
                Dimensioned.none),
            snap("m_f", "Fuel mass",
                Dimensioned.kg),
            snap("rr_f", "Fuel regression depth",
                Dimensioned.mm),
            snap("Pr_inj", "Injector pressure ratio",
                Dimensioned.none),
            snap("dm_inj", "Oxidiser mass flow",
                Dimensioned.kg_s),
            snap("dm_reg", "Fuel mass flow",
                Dimensioned.kg_s),
            snap("dm_out", "Exhaust mass flow",
                Dimensioned.kg_s),
            snap("T_cc", "Chamber gas temperature",
                Dimensioned.K),
            snap("rho_cc", "Chamber gas density",
                Dimensioned.kg_m3),
            snap("a_cc", "Chamber gas speed of sound",
                Dimensioned.m_s),
            snap("Y_cc", "Chamber gas specific heat ratio",
                Dimensioned.none),
            snap("cp_cc", "Chamber gas isobaric specific heat capacity",
                Dimensioned.J_kg_K),
            snap("cv_cc", "Chamber gas isochoric specific heat capacity",
                Dimensioned.J_kg_K),
            snap("Mw_cc", "Chamber gas molar mass",
                Dimensioned.g_mol),
            snap("R_cc", "Chamber gas specific gas constant",
                Dimensioned.J_kg_K),
            snap("T_tht", "Throat gas temperature",
                Dimensioned.K),
            snap("P_tht", "Throat pressure",
                Dimensioned.bar),
            snap("rho_tht", "Throat gas density",
                Dimensioned.kg_m3),
            snap("a_tht", "Throat gas speed of sound",
                Dimensioned.m_s),
            snap("T_exit", "Exit gas temperature",
                Dimensioned.K),
            snap("P_exit", "Exit pressure",
                Dimensioned.bar),
            snap("rho_exit", "Exit gas density",
                Dimensioned.kg_m3),
            snap("a_exit", "Exit gas speed of sound",
                Dimensioned.m_s),
            snap("M_exit", "Exit gas mach number",
                Dimensioned.none),
            snap("ofr", "Oxidiser-fuel ratio",
                Dimensioned.none),
            snap("Isp", "Specific impulse",
                Dimensioned.s),
            snap("Thrust", "Thrust",
                Dimensioned.N),
            snap("Drag", "Drag force",
                Dimensioned.N),
            snap("Gravity", "Gravity force",
                Dimensioned.N),
        ],
    }

@functools.cache
def params():
    return {p.c_name: p for group_name, group in params_grouped().items()
                        for p in group
                        if not group_name.endswith("interp__")}

@functools.cache
def params_snap_interp():
    return {p.c_name: p for p in params_grouped()["__snapshot_interp__"]}

@functools.cache
def snapshot_dtype():
    snap_interp = params_grouped()["__snapshot_interp__"]
    dtype = np.dtype([(p.c_name, np.float64) for p in snap_interp])
    assert dtype.itemsize == 8*len(snap_interp)
    return dtype



### INTERFACE ###################################################################

class Interface:
    SNAPSHOT_COUNT = 100_000

    def __init__(self):
        # Create the bridging object.
        self.interp = bridge.Interpretation()
        for p in params().values():
            self.interp.append(p.c_name, p.c_type, p.c_flags)
        self.interp.finalise()
        self.state = bridge.State(self.interp)

        # Setup snapshotting, first by hashing our interpretation.
        snapshot_interp = bridge.Interpretation()
        for p in params_grouped()["__snapshot_interp__"]:
            snapshot_interp.append(p.c_name, p.c_type, p.c_flags)
        snapshot_interp.finalise()
        # `snap_upto` set before each run.
        self.state["snap_count"] = self.SNAPSHOT_COUNT
        self.state["snap_hash"] = snapshot_interp.get_hash()
        self.state["snap_array"] = np.empty(
            shape=(self.SNAPSHOT_COUNT * snapshot_dtype().itemsize,),
            dtype=np.float64,
        )

        # Setup inputs.
        self.boxes = {} # c_name -> user-input-string.
        for name, p in params().items():
            if not p.is_input:
                continue
            self.input_parsed(name, None)

    def input_direct(self, name, value):
        string = str(value)
        self.state[name] = value
        self.boxes[name] = string
    def input_parsed(self, name, string):
        if string is None:
            self.state[name] = nan
        else:
            self.state[name] = params()[name].caster.to_value(string)
        self.boxes[name] = string

    def save_json(self, path):
        with open(path, "w", encoding="utf-8") as file:
            json.dump(self.boxes, file)

    def load_json(self, path):
        with open(path, "r", encoding="utf-8") as file:
            boxes = json.load(file)
        for name, string in boxes.items():
            self.input_parsed(name, string)

    def sim(self):
        self.state["snap_upto"] = 0
        ret = self.state.execute()
        if ret is not None:
            raise RuntimeError(ret)

    def view(self, what="important"):

        # Extract from bridge as a flat numpy float array, then re-view numpy
        # array in correct dtype.
        snap_interp = params_grouped()["__snapshot_interp__"]
        dtype = snapshot_dtype()
        array = self.state["snap_array"] \
                .sized(self.state["snap_upto"] * len(snap_interp)) \
                .view(dtype)
        def number(name):
            return self.state[name]
        def trace(name):
            caster = params_snap_interp()[name].caster
            return (array[name] - caster.offset) / caster.scale

        t = trace("t")

        if what == "all":
            # Plot each thing against time.
            window = geez.new_window()
            for i in range(1, len(snap_interp)):
                j = ((i - 1) % 6)
                if j == 0:
                    _, axes = window.new_plots(rows=2, cols=3)
                param = snap_interp[i]
                data = trace(param.c_name)
                label = f"{param.display_name} [{param.caster.display}]"
                axes[j // 3, j % 3].plot(t, data)
                axes[j // 3, j % 3].set_title(label)
                if (j // 3) == 1:
                    axes[j // 3, j % 3].set_xlabel("Time [s]")

        elif what != "important":
            assert False, f"unknown 'what': {what}"


        # headline numbers:
        #   apogee
        #   impulse
        #   burn time (liq/total)
        #   prop masses (ox/fuel)
        #   otr velocity
        # headline traces:
        #   thrust
        #   isp
        #   P_t & P_c
        #   Pr_inj
        #   mfr (ox/fuel)
        #   ofr
        ROD_HEIGHT = 5 # [m]

        PALETTE = {
            "red":    "C3",
            "blue":   "C0",
            "amber":  "C1",
            "teal":   "C2",
            "purple": "C4",
            "pink":   "C6",
            "maroon": "#8a0000",
        }

        Nliq = np.argmax(array["m_l"] == 0)
        Nburn = np.argmax(array["Thrust"] != array["Thrust"])
        assert Nliq > 0
        assert Nburn > 0
        Tburn = t[Nburn]

        def otr_velocity():
            diff = array["alt_r"] - number("alt0_r") - ROD_HEIGHT
            sign_changes = np.where(np.diff(np.signbit(diff)))[0]
            if sign_changes.size == 0:
                return nan
            idx = sign_changes[0]

            t0, t1 = t[idx], t[idx+1]
            alt0, alt1 = diff[idx], diff[idx+1]
            t_cross = t0 - alt0 * (t1 - t0) / (alt1 - alt0)

            return np.interp(t_cross, t, array["vel_r"])

        numbers = {
            "apogee": trace("alt_r")[-1],
            "impulse": (np.trapezoid(array["Thrust"][:Nburn], x=t[:Nburn]),
                        np.trapezoid(array["Thrust"][:Nliq], x=t[:Nliq])),
            "burn_time": (t[Nburn], t[Nliq]),
            "prop_masses": (number("m0_o"), number("m0_f")),
            "otr_velocity": otr_velocity(),
        }
        traces = {
            "thrust": 1e-3*trace("Thrust"),
            "isp": trace("Isp"),
            "pressures": (trace("P_t"), trace("P_c")),
            "stiffness": 100*(trace("Pr_inj") - 1),
            "mfr": (trace("dm_inj"), trace("dm_reg")),
            "ofr": trace("ofr"),
        }


        def tile_numbers(ax):
            from matplotlib.patches import FancyBboxPatch
            from matplotlib.colors import to_rgba
            ax.axis("off")
            stats = [  # (label, colour, value)
                ("apogee", PALETTE["amber"],
                    f'{numbers["apogee"]:.2f} kft'
                ),
                ("impulse (liq.)", PALETTE["pink"],
                    f'{numbers["impulse"][0]*1e-3:.1f} '
                    f'({numbers["impulse"][1]*1e-3:.1f}) kNs'
                ),
                ("burn time (liq.)", PALETTE["teal"],
                    f'{numbers["burn_time"][0]:.2f} '
                    f'({numbers["burn_time"][1]:.2f}) s'
                ),
                ("prop. masses (ox/fuel)", PALETTE["purple"],
                    f'{numbers["prop_masses"][0]:.2f} / '
                    f'{numbers["prop_masses"][1]:.2f} kg'
                ),
                ("otr. velocity", PALETTE["blue"],
                    f'{numbers["otr_velocity"]:.1f} m/s'
                ),
            ]
            # patch coords are axes fractions, and this axes is ~12x wider than tall,
            # so isotropic (round-in-pixels) corners need mutation_aspect = w/h
            pos = ax.get_position()
            figw, figh = ax.figure.get_size_inches()
            aspect = (pos.width * figw) / (pos.height * figh)
            n = len(stats)
            w = 0.85 / n
            for i, (label, col, value) in enumerate(stats):
                x = (i + 0.5) / n
                ax.add_patch(FancyBboxPatch(
                    (x - w / 2, 0.06), w, 0.88, transform=ax.transAxes,
                    boxstyle="round,pad=0.006,rounding_size=0.018",
                    mutation_aspect=aspect,
                    facecolor=to_rgba(col, 0.08), edgecolor=to_rgba(col, 0.55),
                    lw=1.3, clip_on=False))
                ax.text(x, 0.72, label, ha="center", va="center", fontsize=11,
                        color="0.45", fontweight="bold", transform=ax.transAxes)
                ax.text(x, 0.28, value, ha="center", va="center", fontsize=15,
                        color=col, fontweight="bold", transform=ax.transAxes)


        def clamp_lower(ax, min_val, axis="y"):
            get, set_ = (ax.get_ylim, ax.set_ylim) \
                        if axis == "y" else \
                        (ax.get_xlim, ax.set_xlim)
            lo, hi = get()
            set_(max(lo, min_val), hi)

        def clamp_run(ax):
            ax.set_xlim(-0.05*t[-1], 1.05*t[-1])
        def clamp_burn(ax):
            ax.set_xlim(-0.05*Tburn, 1.05*Tburn)

        def plot_callout(ax, y, col, pos=1.0, fmt="{:.2f}"):
            ax.axhline(y, color=col, ls="--", lw=1)
            ax.annotate(fmt.format(y), xy=(pos, y),
                        xycoords=ax.get_yaxis_transform(), # x: axes frac, y: data
                        ha="left" if pos >= 1.0 else "center", va="center",
                        color=col, fontsize=10, fontweight="bold")
        def plot_signal(ax, t, y, ycallout, col, pos=1.0, fmt="{:.2f}",
                **plot_kw):
            ax.plot(t, y, color=col, **plot_kw)
            plot_callout(ax, ycallout, col, pos, fmt)

        def tile_thrust(ax):
            plot_signal(ax, t, traces["thrust"], np.nanmax(traces["thrust"]),
                    PALETTE["amber"], pos=1.01, fmt="{:.2f} kN")
            ax.set_title("thrust [kN]", fontweight="bold")
            clamp_burn(ax)
            clamp_lower(ax, 0.0)

        def tile_isp(ax):
            plot_signal(ax, t, traces["isp"], np.nanmax(traces["isp"]),
                    PALETTE["pink"], pos=1.01, fmt="{:.1f} s")
            ax.set_title("isp [s]", fontweight="bold")
            clamp_burn(ax)
            clamp_lower(ax, 0.0)

        def tile_pressures(ax):
            plot_signal(ax, t, traces["pressures"][0],
                    np.nanmax(traces["pressures"][0]),
                    PALETTE["blue"], pos=1.01, fmt="{:.1f} bar",
                    label="tank")
            plot_signal(ax, t, traces["pressures"][1],
                    np.nanmax(traces["pressures"][1]),
                    PALETTE["red"], pos=1.01, fmt="{:.1f} bar",
                    label="cc")
            ax.legend()
            ax.set_title("pressures [bar]", fontweight="bold")
            clamp_burn(ax)
            clamp_lower(ax, 0.0)

        def tile_stiffness(ax):
            ax.plot(t, traces["stiffness"], color=PALETTE["teal"])
            plot_callout(ax, traces["stiffness"][0],
                    PALETTE["teal"], pos=1.01, fmt="{:.1f}%")
            plot_callout(ax, np.nanmin(traces["stiffness"][:Nburn]),
                    PALETTE["teal"], pos=1.01, fmt="{:.1f}%")
            ax.set_title("stiffness [%]", fontweight="bold")
            clamp_burn(ax)
            ax.set_ylim(0.0, 50.0)

        def tile_mfr(ax):
            plot_signal(ax, t, traces["mfr"][0],
                    np.nanmax(traces["mfr"][0]),
                    PALETTE["blue"], pos=1.01, fmt="{:.2f} kg/s",
                    label="ox")
            plot_signal(ax, t, traces["mfr"][1],
                    np.nanmax(traces["mfr"][1]),
                    PALETTE["red"], pos=1.01, fmt="{:.2f} kg/s",
                    label="fuel")
            ax.legend()
            ax.set_title("mfr [kg/s]", fontweight="bold")
            clamp_burn(ax)
            clamp_lower(ax, 0.0)

        def tile_ofr(ax):
            plot_signal(ax, t, traces["ofr"], traces["ofr"][:Nliq].sum() / Nliq,
                    PALETTE["purple"], pos=1.01, fmt="{:.1f} s")
            ax.set_title("ofr [-]", fontweight="bold")
            clamp_burn(ax)
            clamp_lower(ax, 0.0)

        def clean_axes(ax):
            ax.title.set_color("0.15")
            ax.xaxis.label.set_color("0.15")
            ax.yaxis.label.set_color("0.15")
            ax.tick_params(colors="0.15") # tick marks + tick labels
            for spine in ax.spines.values():
                spine.set_color("0.15")

            ax.grid(False, axis="both", which="both")
            autminloc = matplotlib.ticker.AutoMinorLocator
            ax.xaxis.set_minor_locator(autminloc())
            ax.yaxis.set_minor_locator(autminloc())
            ax.grid(True, axis="both", which="major", lw=1.0, ls="-", alpha=0.6)
            ax.grid(True, axis="both", which="minor", lw=0.6, ls="--",
                    alpha=0.35)

        TILES = [
            [tile_numbers],
            [tile_thrust, tile_pressures, tile_mfr],
            [tile_isp,    tile_stiffness, tile_ofr],
        ]
        HEIGHT_RATIOS = [0.22, 1, 0.7]

        nrows, ncols = len(TILES), len(TILES[0])
        fig = geez.new_figure(figsize=(12.0, 7.0), constrained_layout={
            "w_pad": 0.1,
            "h_pad": 0.1,
            "wspace": 0.0,
            "hspace": 0.1,
        })
        gs = fig.add_gridspec(len(TILES), 1, height_ratios=HEIGHT_RATIOS)
        for r, row in enumerate(TILES):
            if r == 0:
                # Row 0 (tile_numbers) spans normally
                sub = gs[r].subgridspec(1, len(row))
                for c, tile in enumerate(row):
                    ax = fig.add_subplot(sub[c])
                    tile(ax)
                    clean_axes(ax)
            else:
                # Rows 1 & 2: Add 3 active columns + 1 tiny empty padding column on the right
                sub = gs[r].subgridspec(1, 4, width_ratios=[1, 1, 1, 0.08])
                for c, tile in enumerate(row):
                    ax = fig.add_subplot(sub[c])  # Only plot in columns 0, 1, and 2
                    tile(ax)
                    clean_axes(ax)
                # Column 3 (sub[3]) is deliberately lef
