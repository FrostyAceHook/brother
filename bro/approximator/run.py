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

__all__ = ["run", "main"]


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

N2O_SIZE_SAT = 2000          # 1D saturated tables (cheap, and steep near Tcrit)
N2O_SIZE_DENSE_MAX = 400     # 1D dense-table ceiling
N2O_SIZE_DENSE = (150, 100)  # (T, psi)
N2O_SIZE_LIGHT = (150, 100)  # (T, psi)

# Fields tabulated in the 2D single-phase tables. T and rho are inputs, X is
# computed analytically in C, and the saturation slopes only exist on the
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
    """Several 1D LUTs (one per field) sharing a grid in T."""

    def __init__(self, prefix, func, fields, Tlo, Thi, size):
        self.prefix = prefix
        self.func = func          # T -> array ordered as `fields`
        self.fields = fields
        self.Tlo = Tlo
        self.Thi = Thi
        self.size = size
        self.T = np.linspace(Tlo, Thi, size)
        self.values = np.array([func(T) for T in self.T]) # (size, nfields)

    def __call__(self, k, T):
        x = (np.asarray(T, dtype=np.float64) - self.Tlo) \
          / (self.Thi - self.Tlo) * (self.size - 1)
        i = np.clip(np.floor(x).astype(int), 0, self.size - 2)
        t = x - i
        v = self.values[:, k].astype(np.float32).astype(np.float64)
        return v[i] + t*(v[i + 1] - v[i])

    def c_source(self, group):
        """C text: grid macros (prefixed by `group`) and one named f32 array
        per field."""
        G = group.upper()
        lines = [
            f"/* {self.prefix}: {self.size}-point 1D LUTs in T. */",
            f"#define {G}_TLO ({self.Tlo!r})",
            f"#define {G}_THI ({self.Thi!r})",
            f"#define {G}_LEN ({self.size})",
            "",
        ]
        for k, field in enumerate(self.fields):
            lines.append(f"static const f32 {self.prefix}_{field}_tbl"
                         f"[{G}_LEN] = {{")
            lines.extend(_fmt_f32_table(self.values[:, k].astype(np.float32)))
            lines.append("};")
        lines.append("")
        return lines

    def run(self):
        Tmid = 0.5*(self.T[:-1] + self.T[1:])
        truth = np.array([self.func(T) for T in Tmid])
        keep = np.abs(Tmid - n2o.Tcrit) > N2O_CRIT_BAND
        print(f"    /* {self.prefix}: {self.size}-point 1D LUTs over "
              f"T in [{self.Tlo:.6g}, {self.Thi:.6g}] K */")
        for k, field in enumerate(self.fields):
            name = f"{self.prefix}_{field}"
            approx = self(k, Tmid)
            e_all = _err_summary(truth[:, k], approx)
            e_op = _err_summary(truth[:, k], approx, keep)
            print(f"    /*   {field:>11}: max err {e_all[0]*100:7.3f}% "
                  f"(of range {e_all[1]*100:7.3f}%), excl. Tcrit+-"
                  f"{N2O_CRIT_BAND:g}K {e_op[0]*100:8.5f}% */")


class N2OTable2D:
    """Several 2D LUTs (one per field) sharing a grid in (T, psi)."""

    def __init__(self, prefix, func, fields, field_idx, Tlo, Thi, size):
        self.prefix = prefix
        self.func = func          # (T, psi) -> array ordered as n2o.FIELDS
        self.fields = fields
        self.field_idx = field_idx
        self.Tlo = Tlo
        self.Thi = Thi
        self.shape = size
        self.T = np.linspace(Tlo, Thi, size[0])
        self.psi = np.linspace(0.0, 1.0, size[1])
        vals = np.empty((size[0], size[1], len(fields)))
        for i, T in enumerate(self.T):
            for j, psi in enumerate(self.psi):
                vals[i, j] = func(T, psi)[field_idx]
        self.values = vals

    def __call__(self, k, T, psi):
        nx, ny = self.shape
        x = (np.asarray(T, dtype=np.float64) - self.Tlo) \
          / (self.Thi - self.Tlo) * (nx - 1)
        y = np.asarray(psi, dtype=np.float64) * (ny - 1)
        i = np.clip(np.floor(x).astype(int), 0, nx - 2)
        j = np.clip(np.floor(y).astype(int), 0, ny - 2)
        t = x - i
        s = y - j
        v = self.values[:, :, k].astype(np.float32).astype(np.float64)
        c0 = v[i, j] + s*(v[i, j + 1] - v[i, j])
        c1 = v[i + 1, j] + s*(v[i + 1, j + 1] - v[i + 1, j])
        return c0 + t*(c1 - c0)

    def c_source(self, group):
        """C text: grid macros (prefixed by `group`) and one named f32 array
        per field, C-ordered (row = T index, col = psi index)."""
        G = group.upper()
        nx, ny = self.shape
        lines = [
            f"/* {self.prefix}: {nx}x{ny} 2D LUTs in (T, psi), psi in [0, 1]. */",
            f"#define {G}_TLO ({self.Tlo!r})",
            f"#define {G}_THI ({self.Thi!r})",
            f"#define {G}_NT ({nx})",
            f"#define {G}_NPSI ({ny})",
            "",
        ]
        for k, field in enumerate(self.fields):
            lines.append(f"static const f32 {self.prefix}_{field}_tbl"
                         f"[{G}_NT*{G}_NPSI] = {{")
            lines.extend(_fmt_f32_table(
                self.values[:, :, k].astype(np.float32).flatten("C")))
            lines.append("};")
        lines.append("")
        return lines

    def run(self):
        # Cell centres are the worst case for bilinear interpolation.
        Tm = 0.5*(self.T[:-1] + self.T[1:])
        Pm = 0.5*(self.psi[:-1] + self.psi[1:])
        TT, PP = np.meshgrid(Tm, Pm, indexing="ij")
        truth = np.empty(TT.shape + (len(self.fields),))
        for i in range(TT.shape[0]):
            for j in range(TT.shape[1]):
                truth[i, j] = self.func(TT[i, j], PP[i, j])[self.field_idx]
        keep = np.abs(TT - n2o.Tcrit) > N2O_CRIT_BAND
        nx, ny = self.shape
        print(f"    /* {self.prefix}: {nx}x{ny} 2D LUTs over T in "
              f"[{self.Tlo:.6g}, {self.Thi:.6g}] K, psi in [0, 1] */")
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
    """Constants and field lists shared by n2o.h/n2o.c. Keeps the C side in
    lock-step with the coordinate mapping used to build the tables."""
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


def find_n2o_tables(sizes=None):
    """Builds all N2O tables, reports their accuracy and writes:
        tbl/N2O_consts.i  constants + X-macro field lists (included by n2o.h)
        tbl/N2O_sat.i     1D saturated liquid/vapour tables
        tbl/N2O_bound.i   1D dense-table ceiling rho_dense_max(T)
        tbl/N2O_dense.i   2D dense (liquid/supercritical) tables in (T, psi)
        tbl/N2O_light.i   2D light (vapour/supercritical) tables in (T, psi)
    """
    sizes = sizes or {}
    n_sat = sizes.get("sat", N2O_SIZE_SAT)
    n_max = sizes.get("dense_max", N2O_SIZE_DENSE_MAX)
    s_dense = sizes.get("dense", N2O_SIZE_DENSE)
    s_light = sizes.get("light", N2O_SIZE_LIGHT)

    idx_sat = [n2o.FIELDS.index(f) for f in N2O_FIELDS_SAT]
    idx_2d = [n2o.FIELDS.index(f) for f in N2O_FIELDS_2D]

    satliq = N2OTable1D("N2O_satliq",
            lambda T: n2o.query_sat(T, 0)[idx_sat],
            N2O_FIELDS_SAT, n2o.Ttrip, n2o.Tsat_hi, n_sat)
    satvap = N2OTable1D("N2O_satvap",
            lambda T: n2o.query_sat(T, 1)[idx_sat],
            N2O_FIELDS_SAT, n2o.Ttrip, n2o.Tsat_hi, n_sat)
    bound = N2OTable1D("N2O_bound", lambda T: np.array([n2o.rho_dense_max(T)]),
            ["rho_dense_max"], n2o.Ttrip, n2o.Tmax, n_max)
    dense = N2OTable2D("N2O_dense", n2o.query_dense, N2O_FIELDS_2D, idx_2d,
            n2o.Ttrip, n2o.Tmax, s_dense)
    light = N2OTable2D("N2O_light", n2o.query_light, N2O_FIELDS_2D, idx_2d,
            n2o.Ttrip, n2o.Tmax, s_light)

    for tbl in (satliq, satvap, bound, dense, light):
        tbl.run()

    _n2o_consts()
    # Both saturated sides share one grid, so they share one file/macro set.
    sat_lines = satliq.c_source("N2O_SAT")
    sat_lines += satvap.c_source("N2O_SAT")[5:] # skip duplicate grid macros
    _write_tbl("N2O_sat", sat_lines)
    _write_tbl("N2O_bound", bound.c_source("N2O_BOUND"))
    _write_tbl("N2O_dense", dense.c_source("N2O_DENSE"))
    _write_tbl("N2O_light", light.c_source("N2O_LIGHT"))

    print("    /* N2O tables written: tbl/N2O_{consts,sat,bound,dense,light}.i */")
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

    with timer("N2O"):
        print("N2O_Ttrip:   ", n2o.Ttrip)
        print("N2O_Ptrip:   ", n2o.Ptrip)
        print("N2O_Tcrit:   ", n2o.Tcrit)
        print("N2O_Pcrit:   ", n2o.Pcrit)
        print("N2O_rhocrit: ", n2o.rhocrit)
        print("N2O_Tmax:    ", n2o.Tmax)
        print("N2O_Pmax:    ", n2o.Pmax)
        print("N2O_rhomin:  ", n2o.rhomin)
        print()
        find_n2o_tables()




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
