"""
Runs the static approximator, outputting all found approximations.
"""

import argparse
from math import prod

import numpy as np
import matplotlib

from .. import paths
from .. import geez
from ..geez import new_figure, new_plots, new_window, no_window

from . import n2o
from . import cea

__all__ = ["run", "main"]



# ============================================================================
# Shared LUT axis (mirrors lutAxis in lut.h).


class LutAxis:
    """Table axis evenly spaced in a transformed coordinate g(x):
        "lin"   g(x) = x
        "log"   g(x) = log2(x)
        "sqrt"  g(x) = sign(x - c) * |x - c|^(1/2)  (nodes dense near c)
        "cbrt"  g(x) = cbrt(x - c)                    (denser still near c)
    Interpolation is linear in g, exactly as lut.c does it."""

    _C_SCALE = {"lin": "LUT_LIN", "log": "LUT_LOG", "sqrt": "LUT_SQRT",
                "cbrt": "LUT_CBRT"}

    def __init__(self, lo, hi, n, scale="lin", c=0.0):
        assert scale in self._C_SCALE, scale
        self.lo = float(lo)
        self.hi = float(hi)
        self.n = int(n)
        self.scale = scale
        self.c = float(c)
        self.glo = float(self.g(self.lo))
        self.ghi = float(self.g(self.hi))

    def g(self, x):
        x = np.asarray(x, dtype=np.float64)
        if self.scale == "log":
            return np.log2(x)
        if self.scale == "sqrt":
            d = x - self.c
            return np.sign(d) * np.sqrt(np.abs(d))
        if self.scale == "cbrt":
            return np.cbrt(x - self.c)
        return x

    def ginv(self, y):
        y = np.asarray(y, dtype=np.float64)
        if self.scale == "log":
            return np.exp2(y)
        if self.scale == "sqrt":
            return self.c + np.sign(y) * y*y
        if self.scale == "cbrt":
            return self.c + y*y*y
        return y

    def nodes(self):
        x = self.ginv(np.linspace(self.glo, self.ghi, self.n))
        x[0] = self.lo # exact endpoints (no round-trip error).
        x[-1] = self.hi
        return x

    def frac(self, x):
        return (self.g(x) - self.glo) / (self.ghi - self.glo) * (self.n - 1)

    def from_frac(self, f):
        t = np.asarray(f, dtype=np.float64) / (self.n - 1)
        return self.ginv(self.glo + t*(self.ghi - self.glo))

    def mids(self):
        """Cell centres in g (the worst case for linear interpolation)."""
        return self.from_frac(np.arange(self.n - 1) + 0.5)

    def c_init(self):
        return (f"{{ "
                f".lo = {self.lo!r}, "
                f".hi = {self.hi!r}, "
                f".glo = {self.glo!r}, "
                f".ghi = {self.ghi!r}, "
                f".n = {self.n}, "
                f".scale = {self._C_SCALE[self.scale]}, "
                f".c = {self.c!r}, "
            f"}}")

    def describe(self):
        about = f" about {self.c:.6g}" if self.scale in ("sqrt", "cbrt") else ""
        return (f"{self.n} nodes over [{self.lo:.6g}, {self.hi:.6g}], "
                f"{self.scale}{about}")


def _interp_axes(values, axes, coords):
    """Multilinear interpolation matching lut.c. `values` has one axis per
    entry of `axes`; `coords` is a list of physical-coordinate arrays."""
    fr = [ax.frac(c) for ax, c in zip(axes, coords)]
    idx = [np.clip(np.floor(f).astype(int), 0, ax.n - 2)
           for f, ax in zip(fr, axes)]
    t = [f - i for f, i in zip(fr, idx)]
    v = np.asarray(values).astype(np.float32).astype(np.float64)
    out = 0.0
    dims = len(axes)
    for corner in range(1 << dims):
        w = 1.0
        sel = []
        for d in range(dims):
            bit = (corner >> d) & 1
            w = w * (t[d] if bit else (1.0 - t[d]))
            sel.append(idx[d] + bit)
        out = out + w * v[tuple(sel)]
    return out


# ============================================================================
# N2O: saturation-aligned lookup tables.
#
# The (T, rho) plane is split along the saturation dome (see n2o.py) so no
# table cell straddles the dome boundary:
#   - 1D saturated tables (liquid and vapour) over T in [Ttrip, Tsat_hi].
#   - 1D dense-table ceiling rho_dense_max(T) over T in [Ttrip, Tmax].
#   - 2D dense table (liquid / dense supercritical) over (T, psi).
#   - 2D light table (vapour / light supercritical) over (T, psi).
# Inside the dome nothing is tabulated; the generated C builds those states
# from the saturated tables with the lever rule.

N2O_AXES = dict(
    sat=LutAxis(n2o.Ttrip, n2o.Tsat_hi, 2000, scale="cbrt", c=n2o.Tcrit),
    bound=LutAxis(n2o.Ttrip, n2o.Tmax, 400),
    T2d=LutAxis(n2o.Ttrip, n2o.Tmax, 220, scale="sqrt", c=n2o.Tcrit),
    psi=LutAxis(0.0, 1.0, 100),
)

# Fields tabulated in the 2D single-phase tables. T and rho are inputs, X is
# computed analytically, and the saturation slopes only exist on the
# saturation line.
N2O_FIELDS_2D = [f for f in n2o.FIELDS
                 if f not in ("T", "rho", "X", "Y", "IPLE")
                 and f not in n2o.FIELDS_SAT_ONLY]
# Fields the C dispatcher derives from tabulated ones rather than tabulating:
#   Y    = cp/cv     (keeps it exactly consistent with cp and cv)
#   IPLE = rho*a^2/P (IPLE ~ 1/P varies hyperbolically along psi at low T,
#                     which bilinear interpolation handles badly; a and P are
#                     both smooth)
N2O_FIELDS_DERIVED = ["Y", "IPLE"]
# Fields tabulated in the 1D saturated tables.
N2O_FIELDS_SAT = [f for f in n2o.FIELDS if f not in ("T", "X")]

# Errors are also reported excluding this band around Tcrit, where saturated
# properties are singular and no table does well.
N2O_CRIT_BAND = 3.0 # [K]


def _fmt_f32_table(tbl):
    # x.01234567ennnf,
    # =17 diggies.
    entrylen = 17
    display_cols = (80 - 4) // entrylen
    lines = []
    for start in range(0, len(tbl), display_cols):
        entries = []
        for val in tbl[start:start + display_cols]:
            s = f"{val:+.8e}f," if val == val else "fNAN,"
            entries.append(" "*(entrylen - len(s)) + s)
        lines.append(" " * 4 + "".join(entries))
    return lines

def _err_summary(truth, approx, keep=None):
    """Max |approx/truth - 1| and max |approx - truth| / range, ignoring NaN."""
    truth = np.asarray(truth, dtype=np.float64)
    approx = np.asarray(approx, dtype=np.float64)
    if keep is not None:
        truth = truth[keep]
        approx = approx[keep]
    ok = np.isfinite(truth) & np.isfinite(approx)
    if not ok.any():
        return float("nan"), float("nan")
    t = truth[ok]
    a = approx[ok]
    with np.errstate(divide="ignore", invalid="ignore"):
        rel_to_value = np.nanmax(np.abs(a / t - 1.0))
    span = np.nanmax(t) - np.nanmin(t)
    rel_to_range = np.nanmax(np.abs(a - t)) / span if span > 0 else 0.0
    return rel_to_value, rel_to_range


class N2OTable1D:
    """Several 1D LUTs (one per field) sharing one T axis."""

    def __init__(self, prefix, func, fields, axis):
        self.prefix = prefix
        self.func = func          # T -> array ordered as `fields`
        self.fields = fields
        self.axis = axis
        self.T = axis.nodes()
        self.values = np.array([func(T) for T in self.T]) # (n, nfields)

    def __call__(self, k, T):
        return _interp_axes(self.values[:, k], [self.axis], [T])

    def c_source(self, group, axis=True):
        """C text: the axis (`<GROUP>_AXIS`) and one named f32 array per
        field."""
        G = group.upper()
        lines = [f"/* {self.prefix}: 1D LUTs in T, {self.axis.describe()}. */"]
        if axis:
            lines += [f"static const lutAxis {G}_AXIS = {self.axis.c_init()};",
                      ""]
        for k, field in enumerate(self.fields):
            lines.append(f"static const f32 {self.prefix}_{field}_tbl"
                         f"[{self.axis.n}] = {{")
            lines.extend(_fmt_f32_table(self.values[:, k].astype(np.float32)))
            lines.append("};")
        lines.append("")
        return lines

    def run(self):
        Tmid = self.axis.mids()
        truth = np.array([self.func(T) for T in Tmid])
        keep = np.abs(Tmid - n2o.Tcrit) > N2O_CRIT_BAND
        print(f"    /* {self.prefix}: {self.axis.describe()} */")
        for k, field in enumerate(self.fields):
            approx = self(k, Tmid)
            e_all = _err_summary(truth[:, k], approx)
            e_op = _err_summary(truth[:, k], approx, keep)
            print(f"    /*   {field:>11}: max err {e_all[0]*100:7.3f}% "
                  f"(of range {e_all[1]*100:7.3f}%), excl. Tcrit+-"
                  f"{N2O_CRIT_BAND:g}K {e_op[0]*100:8.5f}% */")


class N2OTable2D:
    """Several 2D LUTs (one per field) sharing (T, psi) axes."""

    def __init__(self, prefix, func, fields, field_idx, T_axis, psi_axis):
        self.prefix = prefix
        self.func = func          # (T, psi) -> array ordered as n2o.FIELDS
        self.fields = fields
        self.field_idx = field_idx
        self.T_axis = T_axis
        self.psi_axis = psi_axis
        self.T = T_axis.nodes()
        self.psi = psi_axis.nodes()
        vals = np.empty((T_axis.n, psi_axis.n, len(fields)))
        for i, T in enumerate(self.T):
            for j, psi in enumerate(self.psi):
                vals[i, j] = func(T, psi)[field_idx]
        self.values = vals

    def __call__(self, k, T, psi):
        return _interp_axes(self.values[:, :, k], [self.T_axis, self.psi_axis],
                [T, psi])

    def c_source(self, group):
        """C text: the axes (`<GROUP>_T_AXIS`, `<GROUP>_PSI_AXIS`) and one
        named f32 array per field, C-ordered (row = T, col = psi)."""
        G = group.upper()
        lines = [
            f"/* {self.prefix}: 2D LUTs in (T, psi). */",
            f"/*   T:   {self.T_axis.describe()} */",
            f"/*   psi: {self.psi_axis.describe()} */",
            f"static const lutAxis {G}_T_AXIS = {self.T_axis.c_init()};",
            f"static const lutAxis {G}_PSI_AXIS = {self.psi_axis.c_init()};",
            "",
        ]
        n = self.T_axis.n * self.psi_axis.n
        for k, field in enumerate(self.fields):
            lines.append(f"static const f32 {self.prefix}_{field}_tbl[{n}] = {{")
            lines.extend(_fmt_f32_table(
                self.values[:, :, k].astype(np.float32).flatten("C")))
            lines.append("};")
        lines.append("")
        return lines

    def run(self):
        # Cell centres are the worst case for bilinear interpolation.
        TT, PP = np.meshgrid(self.T_axis.mids(), self.psi_axis.mids(),
                             indexing="ij")
        truth = np.empty(TT.shape + (len(self.fields),))
        for i in range(TT.shape[0]):
            for j in range(TT.shape[1]):
                truth[i, j] = self.func(TT[i, j], PP[i, j])[self.field_idx]
        keep = np.abs(TT - n2o.Tcrit) > N2O_CRIT_BAND
        print(f"    /* {self.prefix}: T {self.T_axis.describe()}; "
              f"psi {self.psi_axis.describe()} */")
        for k, field in enumerate(self.fields):
            name = f"{self.prefix}_{field}"
            approx = self(k, TT, PP)
            e_all = _err_summary(truth[..., k], approx)
            e_op = _err_summary(truth[..., k], approx, keep)
            print(f"    /*   {field:>11}: max err {e_all[0]*100:7.3f}% "
                  f"(of range {e_all[1]*100:7.3f}%), excl. Tcrit+-"
                  f"{N2O_CRIT_BAND:g}K {e_op[0]*100:8.5f}% */")
            if summary_window is not None:
                fig, axes = summary_window.new_plots(rows=1, cols=2,
                        title=name)
                axes = np.asarray(axes).ravel()
                cont = axes[0].contourf(TT, PP, truth[..., k], levels=300,
                        cmap="viridis")
                fig.colorbar(cont, ax=axes[0])
                axes[0].set_title("actual function (T, psi)")
                with np.errstate(divide="ignore", invalid="ignore"):
                    err = 100*np.abs(approx / truth[..., k] - 1)
                cont = axes[1].contourf(TT, PP, err, levels=300,
                        cmap="viridis")
                fig.colorbar(cont, ax=axes[1])
                fig.suptitle(name)
                axes[1].set_title("abs error %")


def _write_tbl(name, lines):
    paths.APPROXIMATOR_TBLS.mkdir(parents=True, exist_ok=True)
    with open(paths.APPROXIMATOR_TBLS / f"{name}.i", "w") as f:
        f.write("\n".join([f"/* Generated by approximator/run.py -- do not "
                           f"edit. */"] + lines) + "\n")

def _n2o_consts():
    def xmacro(name, fields):
        return [f"#define {name}(X) \\"] \
             + [f"    X({f}) \\" for f in fields] + [""]
    lines = [
        "/* N2O constants (CoolProp) and table domain. */",
        f"#define N2O_Ttrip   ({n2o.Ttrip!r})",
        f"#define N2O_Ptrip   ({n2o.Ptrip!r})",
        f"#define N2O_Tcrit   ({n2o.Tcrit!r})",
        f"#define N2O_Pcrit   ({n2o.Pcrit!r})",
        f"#define N2O_rhocrit ({n2o.rhocrit!r})",
        f"#define N2O_R       ({n2o.R!r})",
        f"#define N2O_Tmin    ({n2o.Tmin!r})",
        f"#define N2O_Tmax    ({n2o.Tmax!r})",
        f"#define N2O_Pmax    ({n2o.Pmax!r})",
        f"#define N2O_rhomin  ({n2o.rhomin!r})",
        f"#define N2O_Tsat_hi ({n2o.Tsat_hi!r})",
        "",
        "// Fields present in the saturated (1D) tables.",
    ] + xmacro("N2O_TBL_SAT_FIELDS", N2O_FIELDS_SAT) + [
        "// Fields present in the dense/light (2D) tables.",
    ] + xmacro("N2O_TBL_2D_FIELDS", N2O_FIELDS_2D)
    _write_tbl("N2O_consts", lines)


def gen_n2o_tables(axes=None):
    axes = N2O_AXES | (axes or {})

    idx_sat = [n2o.FIELDS.index(f) for f in N2O_FIELDS_SAT]
    idx_2d = [n2o.FIELDS.index(f) for f in N2O_FIELDS_2D]

    satliq = N2OTable1D("N2O_satliq",
            lambda T: n2o.query_sat(T, 0)[idx_sat], N2O_FIELDS_SAT, axes["sat"])
    satvap = N2OTable1D("N2O_satvap",
            lambda T: n2o.query_sat(T, 1)[idx_sat], N2O_FIELDS_SAT, axes["sat"])
    bound = N2OTable1D("N2O_bound", lambda T: np.array([n2o.rho_dense_max(T)]),
            ["rho_dense_max"], axes["bound"])
    dense = N2OTable2D("N2O_dense", n2o.query_dense, N2O_FIELDS_2D, idx_2d,
            axes["T2d"], axes["psi"])
    light = N2OTable2D("N2O_light", n2o.query_light, N2O_FIELDS_2D, idx_2d,
            axes["T2d"], axes["psi"])

    for tbl in (satliq, satvap, bound, dense, light):
        tbl.run()

    _n2o_consts()
    # Both saturated sides share one grid, so they share one file/macro set.
    _write_tbl("N2O_sat", satliq.c_source("N2O_SAT")
                        + satvap.c_source("N2O_SAT", axis=False))
    _write_tbl("N2O_bound", bound.c_source("N2O_BOUND"))
    _write_tbl("N2O_dense", dense.c_source("N2O_DENSE"))
    _write_tbl("N2O_light", light.c_source("N2O_LIGHT"))

    print("    /* N2O tables written: tbl/N2O_{consts,sat,bound,dense,light}.i */")
    print("")






# ============================================================================
# CEA: per-propellant combustion lookup tables.
#
# Every propellant combination gets its own tables over its own domain, written
# to tbl/cea_<key>.i, plus one generated index (tbl/cea_index.i) that maps
# (propellant, station, field) to a table so cea.c can dispatch on a
# `ceaProp` argument.
#
# Axes are evenly spaced in either x or log(x). Chamber pressure is log-spaced
# for the hybrids (it spans ~2 decades there, and equilibrium properties vary
# roughly with log P), as is the supersonic area ratio (Ivac ~ log AR).
#
# Stations, as selected by the C-side `AR` argument:
#   AR == 0      chamber (stagnation)       2D in (P, ofr)
#   AR == 1      throat                     2D in (P, ofr)
#   AR >  1      supersonic, A/At = AR      3D in (P, ofr, AR)
#   0 < AR < 1   subsonic,   A/At = 1/AR    3D in (P, ofr, AR)
#
# Values come from cea.CEA (configured per propellant), whose cacher grid is far
# finer than these tables, so its interpolation is negligible here.


# All (station, field) pairs the C API can ask for. Pairs a propellant doesn't
# tabulate get NULL entries in the index (and assert if requested).
CEA_STATIONS = ["cc", "tht", "sup", "sub"]
CEA_FIELDS = [
    # name,   CEA station name, units (as returned by cea.CEA)
    ("T",     "t",     "temperature [K]"),
    ("P",     "p",     "pressure"),
    ("rho",   "rho",   "density [kg/m^3]"),
    ("M",     "mach",  "Mach number [-]"),
    ("a",     "son",   "speed of sound [m/s]"),
    ("Y",     "gamma", "ratio of specific heats [-]"),
    ("cp",    "cp",    "isobaric heat capacity"),
    ("mu",    "visc",  "dynamic viscosity"),
    ("Pr",    "pran",  "Prandtl number [-]"),
    ("Mw",    "mw",    "molar mass"),
    ("Ivac",  "ivac",  "vacuum specific impulse [s]"),
]
CEA_FIELD_NAMES = [f[0] for f in CEA_FIELDS]
_CEA_NAME = {f[0]: f[1] for f in CEA_FIELDS}

def _cea_source_name(station, field):
    # Chamber values come from the "c_*" outputs; every other station is the
    # exit-plane output of a problem with the matching area ratio.
    return ("c_" if station == "cc" else "") + _CEA_NAME[field]

# Fields per station. Mach is exact at the chamber (0) and throat (1), and Ivac
# has no meaning at the chamber.
_ALL_STN = [f for f in CEA_FIELD_NAMES]
_CC = [f for f in CEA_FIELD_NAMES if f not in ("M", "Ivac")]
_THT = [f for f in CEA_FIELD_NAMES if f != "M"]

# Propellants. `enum` must match `ceaProp` in cea.h.
#   P    : chamber pressure axis [Pa]
#   ofr  : oxidiser/fuel ratio axis
#   sup/sub : area-ratio axes (None to skip that station)
#   fields  : which fields to tabulate at each station
CEA_PROPELLANTS = [
    dict(
        enum="CEA_N2O_PARAFFIN", key="n2o_paraffin", ox="n2o", fuel="paraffin",
        P=LutAxis(0.5e5, 90e5, 30, scale="log"), ofr=LutAxis(1.0, 15.0, 351),
        P3=LutAxis(0.5e5, 90e5, 12, scale="log"), ofr3=LutAxis(1.0, 15.0, 141),
        sup=LutAxis(1.2, 12.0, 24, scale="log"), sub=None,
        fields=dict(cc=_CC, tht=_THT, sup=["Ivac"], sub=[]),
    ),
]

CEA_N_CHECK = 400 # random points per table group for the error check.


def _interp_nd(values, axes, coords):
    return _interp_axes(values, axes, coords)


class CeaTableGroup:
    """All fields of one (propellant, station), sharing one grid."""

    def __init__(self, prop, station, axes, fields, evaluate):
        self.prop = prop
        self.station = station
        self.axes = axes          # list of LutAxis (2 or 3)
        self.fields = fields      # list of field names
        self.evaluate = evaluate  # (field, *coords) -> array
        self.name = f"cea_{prop['key']}_{station}"
        grids = np.meshgrid(*[ax.nodes() for ax in axes], indexing="ij")
        self.values = {f: np.asarray(evaluate(f, *grids), dtype=np.float64)
                       for f in fields}

    def check(self, rng):
        # Random points, uniform in table-index space (so log axes are sampled
        # evenly in log).
        coords = [ax.from_frac(rng.uniform(0, ax.n - 1, CEA_N_CHECK))
                  for ax in self.axes]
        print(f"    /* {self.name}: "
              + " x ".join(f"{ax.n}({ax.scale})"
                           for ax in self.axes) + " */")
        for f in self.fields:
            truth = np.asarray(self.evaluate(f, *coords), dtype=np.float64)
            approx = _interp_nd(self.values[f], self.axes, coords)
            e = _err_summary(truth, approx)
            print(f"    /*   {f:>6}: max err {e[0]*100:.3g}% "
                  f"(of range {e[1]*100:.3g}%) */")
            if summary_window is not None and len(self.axes) == 2:
                X, Y = np.meshgrid(*[ax.nodes() for ax in self.axes],
                                   indexing="ij")
                fig, axs = summary_window.new_plots(rows=1, cols=1,
                        title=f"LUT {self.name}_{f}")
                axs = np.asarray(axs).ravel()
                cont = axs[0].contourf(X, Y, self.values[f], levels=300,
                        cmap="viridis")
                fig.colorbar(cont, ax=axs[0])
                axs[0].set_title(f"{self.name}_{f} (table nodes)")

    def c_source(self):
        lines = [
            f"static const ceaGrid_ {self.name}_grid_ = {{ {len(self.axes)}, {{",
        ] + [f"    {ax.c_init()}," for ax in self.axes] + ["} };"]
        for f in self.fields:
            tbl = self.values[f].astype(np.float32).flatten("C")
            lines.append(f"static const f32 {self.name}_{f}_tbl[{len(tbl)}] = {{")
            lines.extend(_fmt_f32_table(tbl))
            lines.append("};")
        lines.append("")
        return lines


def gen_cea_tables(propellants=None, seed=0):
    """Builds, checks and writes the CEA tables for each propellant, then the
    dispatch index. Writes tbl/cea_<key>.i and tbl/cea_index.i."""
    propellants = propellants or CEA_PROPELLANTS
    rng = np.random.default_rng(seed)
    present = {} # (station, field) -> {enum: group name}

    for prop in propellants:
        print(f"    /* CEA {prop['ox']}/{prop['fuel']} ({prop['enum']}) */")
        with cea.CEA.configure(prop["ox"], prop["fuel"]):

            def station_eval(station):
                def ev(field, P, ofr, AR=None):
                    if station == "cc":
                        AR = 1.0 # chamber outputs don't depend on AR.
                    elif station == "tht":
                        AR = 1.0
                    getter = cea.CEA[_cea_source_name(station, field)]
                    return getter(P, ofr, AR)
                return ev

            groups = []
            for station in ("cc", "tht"):
                fields = prop["fields"][station]
                if fields:
                    groups.append(CeaTableGroup(prop, station,
                            [prop["P"], prop["ofr"]], fields,
                            station_eval(station)))
            for station in ("sup", "sub"):
                fields = prop["fields"][station]
                if prop[station] is not None and fields:
                    groups.append(CeaTableGroup(prop, station,
                            [prop["P3"], prop["ofr3"], prop[station]], fields,
                            station_eval(station)))

            for g in groups:
                g.check(rng)

        lines = [f"/* CEA {prop['ox']}/{prop['fuel']} ({prop['enum']}). */", ""]
        for g in groups:
            lines += g.c_source()
            for f in g.fields:
                present.setdefault((g.station, f), {})[prop["enum"]] = g.name
        _write_tbl(f"cea_{prop['key']}", lines)

    # Index: one array per (station, field) the C API knows about, indexed by
    # ceaProp. Missing entries are zero (NULL), which cea.c asserts on.
    lines = ["/* (station, field) -> per-propellant table. */"]
    lines += [f'#include "cea_{p["key"]}.i"' for p in propellants]
    lines.append("")
    keys = [(s, f) for s in CEA_STATIONS for f in CEA_FIELD_NAMES] \
         + [("perfexp", "AEAT")]
    for station, f in keys:
        lines.append(f"static const ceaTable_ cea_{station}_{f}_tbls_"
                     f"[CEA_PROP_COUNT] = {{")
        # Every propellant gets an explicit entry (NULL if not tabulated), so
        # no array is ever empty-initialised (ISO C forbids `= {}` pre-C23).
        have = present.get((station, f), {})
        for p in propellants:
            gname = have.get(p["enum"])
            if gname is None:
                lines.append(f"    [{p['enum']}] = {{ NULL, NULL }},")
            else:
                lines.extend([f"    [{p['enum']}] = {{",
                              f"        {gname}_{f}_tbl,",
                              f"        &{gname}_grid_",
                              f"    }},"])
        lines.append("};")
    lines.append("")
    _write_tbl("cea_index", lines)
    print("    /* CEA tables written: tbl/cea_index.i, "
          + ", ".join(f"tbl/cea_{p['key']}.i" for p in propellants) + " */")
    print("")







def _run():

    from contextlib import contextmanager
    import time
    @contextmanager
    def timer(title, count=None):
        start = time.perf_counter()
        try:
            yield
        finally:
            dt = time.perf_counter() - start
            if count is None:
                print(f"{title}: {dt:.6f} s")
            elif dt > 0:
                per = dt / count
                n_3h = int(count * (3 * 3600) / dt)
                print(f"{title}: {dt:.6f} s "
                      f"({count} iters; {per*1e3:.1f} ms/iter; "
                      f"{n_3h:,} iters in 3 hr)")
            else:
                print(f"{title}: {dt:.6f} s "
                      f"({count} iters; too fast to extrapolate)")

    # with timer("N2O"):
    #     gen_n2o_tables()

    with timer("CEA"):
        gen_cea_tables()




def run(view=True, save=False):
    global summary_window
    # N2O tables are built directly from CoolProp at the table nodes (no
    # cacher), so there's no double interpolation and nothing to configure.
    with geez.instance(view, save):
        # stoud redirect nested to not capture meta output.
        with paths.splice_stdout(paths.APPROXIMATOR_OUTPUT, view, save):
            summary_window = geez.new_window("summary", emptyok=True)
            _run()


def main():
    parser = argparse.ArgumentParser(
        description="Static function approximator, generates c source code "
                    "which can be used as an approximation of some mathematical"
                    "/data-derived function."
    )

    group = parser.add_mutually_exclusive_group()
    group.add_argument("-v", "--view", action="store_true",
            help="view output (default, unless saving)")
    group.add_argument("-n", "--no-view", action="store_true",
            help="do not view output")

    parser.add_argument("-s", "--save", action="store_true",
            help="save output to file, inc. figures")

    args = parser.parse_args()

    save = args.save
    view = (args.view | (not args.save)) & (not args.no_view)
    run(view, save)

if __name__ == "__main__":
    main()
