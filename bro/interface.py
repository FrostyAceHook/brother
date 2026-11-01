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

        window = geez.new_window()

        # Extract from bridge as a flat numpy float array, then re-view numpy
        # array in correct dtype.
        snap_interp = params_grouped()["__snapshot_interp__"]
        dtype = snapshot_dtype()
        array = self.state["snap_array"] \
                .sized(self.state["snap_upto"] * len(snap_interp)) \
                .view(dtype)
        def get(name):
            caster = params_snap_interp()[name].caster
            return (array[name] - caster.offset) / caster.scale

        # Plot each thing against time.
        t = get("t")
        for i in range(1, len(snap_interp)):
            j = ((i - 1) % 6)
            if j == 0:
                _, axes = window.new_plots(rows=2, cols=3)
            param = snap_interp[i]
            data = get(param.c_name)
            label = f"{param.display_name} [{param.caster.display}]"
            axes[j // 3, j % 3].plot(t, data)
            axes[j // 3, j % 3].set_title(label)
            if (j // 3) == 1:
                axes[j // 3, j % 3].set_xlabel("Time [s]")
