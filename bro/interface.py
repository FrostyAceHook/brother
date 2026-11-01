"""
Interface object which understands the state struct and interacts directly with
the c back-end.
"""

import ast
import json
import math
import sys
import time
from dataclasses import dataclass
from math import pi, nan, inf

import numpy as np

from . import bridge
from . import geez
from . import paths

INTERP_F64 = bridge.Interpretation.F64
INTERP_I64 = bridge.Interpretation.I64
INTERP_PTR_F32 = bridge.Interpretation.PTR_F32
INTERP_PTR_F64 = bridge.Interpretation.PTR_F64
INTERP_PTR_I8  = bridge.Interpretation.PTR_I8
INTERP_PTR_I16 = bridge.Interpretation.PTR_I16
INTERP_PTR_I32 = bridge.Interpretation.PTR_I32
INTERP_PTR_I64 = bridge.Interpretation.PTR_I64
INTERP_PTR_U8  = bridge.Interpretation.PTR_U8
INTERP_PTR_U16 = bridge.Interpretation.PTR_U16
INTERP_PTR_U32 = bridge.Interpretation.PTR_U32
INTERP_PTR_U64 = bridge.Interpretation.PTR_U64

INTERP_IN  = bridge.Interpretation.INPUT
INTERP_OUT = bridge.Interpretation.OUTPUT
INTERP_INDATA  = bridge.Interpretation.INPUT_DATA
INTERP_OUTDATA = bridge.Interpretation.OUTPUT_DATA



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
            if isinstance(subnode, (ast.Assign, ast.AnnAssign, ast.AugAssign, ast.NamedExpr)):
                raise ValueError(f"illegal assigning node {repr(type(subnode).__name__)}")
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
Dimensioned.m = Dimensioned("m")
Dimensioned.mm = Dimensioned("mm", scale=1e-3)
Dimensioned.ft = Dimensioned("ft", scale=1.0 / 3.281)
Dimensioned.kft = Dimensioned("kft", scale=1e3 / 3.281)
Dimensioned.mm2 = Dimensioned("mm^2", scale=1e-6)
Dimensioned.L = Dimensioned("L", scale=1e-3)
Dimensioned.s = Dimensioned("s")
Dimensioned.N = Dimensioned("N")
Dimensioned.Pa = Dimensioned("Pa")
Dimensioned.MPa = Dimensioned("MPa", scale=1e6)
Dimensioned.bar = Dimensioned("bar", scale=1e5)
Dimensioned.K = Dimensioned("K")
Dimensioned.dC = Dimensioned("°C", offset=273.15)
Dimensioned.kg = Dimensioned("kg")
Dimensioned.kg_m3 = Dimensioned("kg/m^3")
Dimensioned.kg_s = Dimensioned("kg/s")
Dimensioned.J_kg_K = Dimensioned("J/kg/K")
Dimensioned.J_K = Dimensioned("J/K")
Dimensioned.deg = Dimensioned("°", scale=pi/180.0)
Dimensioned.percent = Dimensioned("%", scale=1e-2)


class Bool(Caster):
    def __init__(self):
        super().__init__("(y/n)")
    def to_value(self, string):
        if string.strip().casefold() in {"y", "yes", "yeah", "yeagh", "t", "true"}:
            return True
        if string.strip().casefold() in {"n", "no", "nah", "nogh", "f", "false"}:
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

    @classmethod
    def input(cls, c_name, display_name, caster, **kw):
        opt = dict(c_name=c_name, display_name=display_name, caster=caster,
                   c_type=INTERP_F64, c_flags=INTERP_IN)
        opt.update(kw)
        return cls(**opt)
    @classmethod
    def output(cls, c_name, display_name, caster, **kw):
        opt = dict(c_name=c_name, display_name=display_name, caster=caster,
                   c_type=INTERP_F64, c_flags=INTERP_OUT)
        opt.update(kw)
        return cls(**opt)
    @classmethod
    def silent(cls, c_name, **kw):
        opt = dict(c_name=c_name, c_type=INTERP_F64, c_flags=INTERP_OUT)
        opt.update(kw)
        return cls(**opt)
    @classmethod
    def either(cls, c_name, display_name, caster, **kw):
        opt = dict(c_name=c_name, display_name=display_name, caster=caster,
                   c_type=INTERP_F64, c_flags=INTERP_IN | INTERP_OUT)
        opt.update(kw)
        return cls(**opt)

    @classmethod
    def output_array(cls, c_name, **kw):
        opt = dict(c_name=c_name, c_type=INTERP_PTR_F64,
                   c_flags=INTERP_IN | INTERP_OUTDATA)
        opt.update(kw)
        return cls(**opt)

    @property
    def is_displayed(self):
        return self.display_name is not None
    @property
    def is_input(self):
        return self.is_displayed and bool(self.c_flags & INTERP_IN)
    @property
    def is_output(self):
        return self.is_displayed and bool(self.c_flags & INTERP_OUT)



### STATE INTERPRETATION ########################################################

PARAMS_GROUPED = {
    "Setup": [
        Parameter.input("static_fire", "Static fire?", Bool()),
        Parameter.input("alt0_r", "Initial altitude", Dimensioned.m),
    ],
    "Rocket": [
        Parameter.input("D_r", "Outer diameter", Dimensioned.mm),
        Parameter.silent("A_r"),
        Parameter.input("mdry_r", "Dry mass", Dimensioned.kg),
        Parameter.output("mwet0_r", "Wet mass", Dimensioned.kg),
    ],
    "Oxidiser Tank": [
        Parameter.input("V_t", "Volume", Dimensioned.L),
        Parameter.input("C_tw", "Heat capacity", Dimensioned.J_K),
    ],
    "Oxidiser Pool": [
        Parameter.input("T0_o", "Initial temperature", Dimensioned.dC),
        Parameter.either("m0_o", "Initial mass", Dimensioned.kg),
        Parameter.either("vff0_l", "Initial liq. vol. fraction", Dimensioned.none),
    ],
    "Injector": [
        Parameter.input("Cd_inj", "Discharge coefficient", Dimensioned.none),
        Parameter.either("no_injo", "Orifice count", Dimensioned.count),
        Parameter.either("D_injo", "Orifice diameter", Dimensioned.mm),
        Parameter.either("Deq_injo", "Orifice eq. single diameter", Dimensioned.mm),
        Parameter.output("A_injo", "Orifice net area", Dimensioned.mm2),
        Parameter.silent("CdA_inj"),
    ],
    "Combustion Chamber": [
        Parameter.input("etaCstar", "C* efficiency", Dimensioned.percent),
        Parameter.input("etaCf", "Cf efficiency", Dimensioned.percent),
        Parameter.input("Vempty_c", "Empty volume", Dimensioned.L),
    ],
    "Fuel Grain": [
        Parameter.input("L0_f", "Length", Dimensioned.mm),
        Parameter.input("OD0_f", "Outer diameter", Dimensioned.mm),
        Parameter.input("th0_f", "Thickness", Dimensioned.mm),
        Parameter.either("rho_f", "Density", Dimensioned.kg_m3),
        Parameter.either("m0_f", "Mass", Dimensioned.kg),
    ],
    "Nozzle": [
        Parameter.either("D_tht", "Throat diameter", Dimensioned.mm),
        Parameter.either("A_tht", "Throat area", Dimensioned.mm2),
        Parameter.either("D_exit", "Exit diameter", Dimensioned.mm),
        Parameter.either("A_exit", "Exit area", Dimensioned.mm2),
        Parameter.either("AEAT", "Exit-throat area ratio", Dimensioned.none),
        Parameter.input("NLF", "Length fraction", Dimensioned.percent),
        Parameter.output("theta_div", "Diverging half-angle", Dimensioned.deg),
        Parameter.output("theta_exit", "Exit half-angle", Dimensioned.deg),
    ],
    "__outputs__": [
        Parameter.silent("out_upto", c_type=INTERP_I64,
                c_flags=INTERP_IN | INTERP_OUT),
        Parameter.silent("out_count", c_type=INTERP_I64, c_flags=INTERP_IN),
        Parameter.output_array("out_t"),
        Parameter.output_array("out_T_o"),
        Parameter.output_array("out_m_o"),
        Parameter.output_array("out_rr_f"),
        Parameter.output_array("out_alt_r"),
        Parameter.output_array("out_vel_r"),
    ]
}
PARAMS = {p.c_name: p for group in PARAMS_GROUPED.values() for p in group}



### INTERFACE ###################################################################

class Interface:
    OUTPUT_ARRAY_COUNT = 1_000_000

    def __init__(self):
        # Create the bridging object.
        self.interp = bridge.Interpretation()
        for param in PARAMS.values():
            self.interp.append(param.c_name, param.c_type, param.c_flags)
        self.interp.finalise()
        self.state = bridge.State(self.interp)

        # Initialise output arrays.
        self.state["out_upto"] = 0
        self.state["out_count"] = self.OUTPUT_ARRAY_COUNT
        for p in PARAMS_GROUPED["__outputs__"][2:]:
            assert p.c_type == INTERP_PTR_F64
            assert p.c_flags == (INTERP_IN | INTERP_OUTDATA)
            self.state[p.c_name] = np.empty(
                shape=(self.OUTPUT_ARRAY_COUNT,),
                dtype=np.float64,
            )

        # Setup inputs.
        self.boxes = {} # c_name -> user-input-string.
        for name, p in PARAMS.items():
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
            self.state[name] = PARAMS[name].caster.to_value(string)
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
        ret = self.state.execute()
        if ret is not None:
            raise RuntimeError(ret)
