"""
CoolProp nitrous oxide lookups.

All properties are mass-specific SI. Every query returns the same `FIELDS`
struct (unused fields are NaN).

Query kinds:
  - N2O()       : cached general state at (T, rho). Inside the dome, cv and
                  the partial derivatives are the equilibrium-mixture values,
                  while cp, Y, a and IPLE (undefined for a mixture) are NaN.
                  Saturation-line slopes are always NaN.
  - N2Osatliq() : cached saturated liquid at T, evaluated as the single-phase
  - N2Osatvap() : cached saturated vapour at T, (one-sided limit), so every
                  property is defined.
  - query_dense(T, psi) / query_light(T, psi) : uncached single-phase states in
                  saturation-aligned coordinates (see below). These are what
                  the lookup tables are built from.

Saturation-aligned coordinates
------------------------------
The (T, rho) plane is split along the saturation dome so that no lookup-table
cell ever straddles the dome boundary (where P has a kink and single-phase
properties become undefined):

  dense (liquid / dense supercritical):
      rho = rho_split_l(T) + psi * (rho_dense_max(T) - rho_split_l(T))
  light (vapour / light supercritical):
      rho = rhomin * (rho_split_v(T) / rhomin) ** psi

with psi in [0, 1]. rho_split_l/v(T) are the saturated densities below Tcrit,
and both equal rhocrit above it (continuous at the critical point).
rho_dense_max(T) is the density at Pmax. Inside the dome nothing is tabulated:
states there are built from the saturated queries via the lever rule.
"""

import math

import CoolProp
import CoolProp.CoolProp as CP
import numpy as np

from . import cacher

__all__ = [
    "FIELDS", "FIELDS_THERMO", "FIELDS_SAT_ONLY", "FIELDS_SINGLE_PHASE_ONLY", "N2O_Properties",
    "N2O", "N2Osatliq", "N2Osatvap",
    "Ttrip", "Ptrip", "Tcrit", "Pcrit", "rhocrit", "Mw", "R",
    "Tmin", "Tmax", "Pmax", "rhomin", "Tsat_hi", "Psat_hi", "Them_lo",
    "rho_split_l", "rho_split_v", "rho_dense_max",
    "rho_from_dense_psi", "rho_from_light_psi",
    "psi_dense", "psi_light",
    "query", "query_sat", "query_dense", "query_light",
    "Tsat", "hem_choke_T_sat",
]


# Tags cache files so stale fields aren't used. Bump on any field change.
_VERSION = "v3"


FIELDS = [
    "T",    # Temperature [K]
    "P",    # Pressure [Pa]
    "rho",  # Density [kg/m^3]
    "s",    # Entropy [J/kg/K]
    "h",    # Enthalpy [J/kg]
    "u",    # Internal energy [J/kg]
    "X",    # Vapour quality [-] (<0 compressed liquid, >1 superheated vapour,
            #   NaN supercritical; exactly 0/1 for saturated queries)
    "Z",    # Compressibility factor [-]
]
FIELDS_THERMO = [
    # Single-phase values outside the dome. Inside it, cp, Y, a and IPLE are
    # NaN (undefined for an equilibrium mixture: cp is infinite, and the
    # mixture sound speed is a different quantity), while cv and the partials
    # are the equilibrium-mixture values (see `_mixture_partials`).
    "cp",   # Isobaric heat capacity [J/kg/K]
    "cv",   # Isochoric heat capacity [J/kg/K]
    "Y",    # Heat capacity ratio cp/cv [-]
    "a",    # Speed of sound [m/s]
    "IPLE", # Isentropic power law exponent rho*a^2/P [-]
    "dP_dT_rho",  # (dP/dT) at const rho [Pa/K]
    "dP_drho_T",  # (dP/drho) at const T [Pa*m^3/kg]
    "du_drho_T",  # (du/drho) at const T [J*m^3/kg^2]
    "dh_dT_rho",  # (dh/dT) at const rho [J/kg/K]
    "dh_drho_T",  # (dh/drho) at const T [J*m^3/kg^2]
    "ds_drho_T",  # (ds/drho) at const T [J*m^3/kg^2/K]
]
# Subset of FIELDS_THERMO that is NaN inside the dome.
FIELDS_SINGLE_PHASE_ONLY = ["cp", "Y", "a", "IPLE"]
FIELDS_SAT_ONLY = [
    # Saturated queries only (NaN for all others). Total derivatives along the
    # saturation curve, not partials:
    "dPsat_dT",   # Saturation pressure slope [Pa/K]
    "drhosat_dT", # Density slope [kg/m^3/K]
    "dssat_dT",   # Entropy slope [J/kg/K^2]
    "dhsat_dT",   # Enthalpy slope [J/kg/K]
    "dusat_dT",   # Internal energy slope [J/kg/K]
]
FIELDS = FIELDS + FIELDS_THERMO + FIELDS_SAT_ONLY
N2O_Properties = cacher.cls_Result(FIELDS)


_FLUID = "NitrousOxide"

Ttrip = CP.PropsSI("Ttriple", _FLUID)       # 182.33 K
Ptrip = CP.PropsSI("ptriple", _FLUID)       # 87.8 kPa
Tcrit = CP.PropsSI("Tcrit", _FLUID)         # 309.52 K
Pcrit = CP.PropsSI("pcrit", _FLUID)         # 7.245 MPa
rhocrit = CP.PropsSI("rhocrit", _FLUID)     # 452.9 kg/m^3
Mw = CP.PropsSI("molar_mass", _FLUID)       # 0.0440128 kg/mol
R = CP.PropsSI("gas_constant", _FLUID) / Mw # specific gas constant [J/kg/K]

_T_EOS_MAX = CP.PropsSI("Tmax", _FLUID)     # 525 K

# Domain of the lookup tables.
Tmin = Ttrip
Tmax = 333.15      # [K]
Pmax = 15e6        # [Pa], keeps rho_dense_max > rhocrit up to Tmax.
rhomin = 1.0       # [kg/m^3]
# Saturated tables stop just short of the critical point, where saturated
# properties are singular.
Tsat_hi = Tcrit - 0.01
Psat_hi = CP.PropsSI("P", "T", Tsat_hi, "Q", 0, _FLUID) # Psat(Tsat_hi) [Pa]
# Lower end of the tabulated HEM choke: the choke sits ~1-20 K below the inlet
# temperature, and the root needs some room above the triple point to bracket.
Them_lo = Ttrip + 2.0


def _new_state():
    return CoolProp.AbstractState("HEOS", _FLUID)

# Module-level scratch states for the uncached helpers.
_st = _new_state()
_st_sat = _new_state()


# -------------------------------------------------------------------------------
# Coordinate mappings.

def _sat_rho(T, Q):
    _st_sat.unspecify_phase()
    _st_sat.update(CoolProp.QT_INPUTS, Q, T)
    return _st_sat.rhomass()

def rho_split_l(T):
    """Dense/dome boundary density: saturated liquid below Tcrit, rhocrit above."""
    return _sat_rho(T, 0.0) if T < Tcrit else rhocrit

def rho_split_v(T):
    """Light/dome boundary density: saturated vapour below Tcrit, rhocrit above."""
    return _sat_rho(T, 1.0) if T < Tcrit else rhocrit

def rho_dense_max(T):
    """Upper density bound of the dense table: density at Pmax."""
    _st_sat.unspecify_phase()
    _st_sat.update(CoolProp.PT_INPUTS, Pmax, T)
    return _st_sat.rhomass()

def rho_from_dense_psi(T, psi):
    rl = rho_split_l(T)
    return rl + psi * (rho_dense_max(T) - rl)

def rho_from_light_psi(T, psi):
    return rhomin * (rho_split_v(T) / rhomin) ** psi

def psi_dense(T, rho):
    rl = rho_split_l(T)
    return (rho - rl) / (rho_dense_max(T) - rl)

def psi_light(T, rho):
    return math.log(rho / rhomin) / math.log(rho_split_v(T) / rhomin)


# -------------------------------------------------------------------------------
# Property evaluation.

def _single_phase(st, out):
    cp = st.cpmass()
    cv = st.cvmass()
    a = st.speed_sound()
    out["cp"] = cp
    out["cv"] = cv
    out["Y"] = cp / cv
    out["a"] = a
    out["IPLE"] = st.rhomass() * a * a / st.p()
    d = st.first_partial_deriv
    out["dP_dT_rho"] = d(CoolProp.iP, CoolProp.iT, CoolProp.iDmass)
    out["dP_drho_T"] = d(CoolProp.iP, CoolProp.iDmass, CoolProp.iT)
    out["du_drho_T"] = d(CoolProp.iUmass, CoolProp.iDmass, CoolProp.iT)
    out["dh_dT_rho"] = d(CoolProp.iHmass, CoolProp.iT, CoolProp.iDmass)
    out["dh_drho_T"] = d(CoolProp.iHmass, CoolProp.iDmass, CoolProp.iT)
    out["ds_drho_T"] = d(CoolProp.iSmass, CoolProp.iDmass, CoolProp.iT)

def _common(st, out):
    out["T"] = st.T()
    out["P"] = st.p()
    out["rho"] = st.rhomass()
    out["s"] = st.smass()
    out["h"] = st.hmass()
    out["u"] = st.umass()
    out["Z"] = st.compressibility_factor()

def _mixture_partials(sat, T, rho, out):
    """Equilibrium-mixture cv and partial derivatives inside the dome, from the
    saturated states and their slopes along the saturation curve. With
    v = 1/rho, X = (v - v_l)/v_fg and phi = phi_l + X*phi_fg (phi in u, h, s):
        (dphi/drho)_T = -phi_fg / (rho^2 * v_fg)
        (dphi/dT)_rho = phi_l' + X*phi_fg' + phi_fg * dX/dT,
            dX/dT = -(v_l' + X*v_fg') / v_fg
        (dP/dT)_rho = dPsat/dT,  (dP/drho)_T = 0,  cv = (du/dT)_rho
    where ' is the derivative along the saturation curve."""
    fsd = sat.first_saturation_deriv
    sides = []
    for Q in (0.0, 1.0):
        sat.unspecify_phase()
        sat.update(CoolProp.QT_INPUTS, Q, T)
        sides.append(dict(
            v=1.0/sat.rhomass(),
            dv=-fsd(CoolProp.iDmass, CoolProp.iT) / sat.rhomass()**2,
            u=sat.umass(), du=fsd(CoolProp.iUmass, CoolProp.iT),
            h=sat.hmass(), dh=fsd(CoolProp.iHmass, CoolProp.iT),
            s=sat.smass(), ds=fsd(CoolProp.iSmass, CoolProp.iT),
            dP=fsd(CoolProp.iP, CoolProp.iT),
        ))
    l, v = sides
    vfg = v["v"] - l["v"]
    X = (1.0/rho - l["v"]) / vfg
    dXdT = -(l["dv"] + X*(v["dv"] - l["dv"])) / vfg
    def dT_rho(k):
        return l["d" + k] + X*(v["d" + k] - l["d" + k]) + (v[k] - l[k])*dXdT
    def drho_T(k):
        return -(v[k] - l[k]) / (rho*rho*vfg)
    out["cv"] = dT_rho("u")
    out["dP_dT_rho"] = l["dP"]
    out["dP_drho_T"] = 0.0
    out["du_drho_T"] = drho_T("u")
    out["dh_dT_rho"] = dT_rho("h")
    out["dh_drho_T"] = drho_T("h")
    out["ds_drho_T"] = drho_T("s")

def _quality(sat, T, rho):
    if T >= Tcrit:
        return math.nan
    sat.unspecify_phase()
    sat.update(CoolProp.QT_INPUTS, 0.0, T)
    rho_l = sat.rhomass()
    sat.update(CoolProp.QT_INPUTS, 1.0, T)
    rho_v = sat.rhomass()
    return (1.0/rho - 1.0/rho_l) / (1.0/rho_v - 1.0/rho_l)

def _query(st, sat, T, rho):
    out = {}
    if not (Ttrip <= T <= _T_EOS_MAX) or not (rho > 0.0):
        return out
    X = _quality(sat, T, rho)
    out["X"] = X
    st.unspecify_phase()
    st.update(CoolProp.DmassT_INPUTS, rho, T)
    _common(st, out)
    inside_dome = (T < Tcrit) and (0.0 < X < 1.0)
    if inside_dome:
        _mixture_partials(sat, T, rho, out)
    else:
        _single_phase(st, out)
    return out

def _query_sat(st, T, Q):
    out = {}
    if not (Ttrip <= T < Tcrit):
        return out
    # Saturation-line slopes need the saturated state itself.
    st.unspecify_phase()
    st.update(CoolProp.QT_INPUTS, Q, T)
    rho_sat = st.rhomass()
    fsd = st.first_saturation_deriv
    out["dPsat_dT"] = fsd(CoolProp.iP, CoolProp.iT)
    out["drhosat_dT"] = fsd(CoolProp.iDmass, CoolProp.iT)
    out["dssat_dT"] = fsd(CoolProp.iSmass, CoolProp.iT)
    out["dhsat_dT"] = fsd(CoolProp.iHmass, CoolProp.iT)
    out["dusat_dT"] = fsd(CoolProp.iUmass, CoolProp.iT)
    # Everything else from the single-phase EOS evaluated exactly at the
    # saturated density, with the phase forced so CoolProp doesn't classify
    # the boundary point as two-phase.
    try:
        st.specify_phase(CoolProp.iphase_liquid if Q == 0
                         else CoolProp.iphase_gas)
        st.update(CoolProp.DmassT_INPUTS, rho_sat, T)
        _common(st, out)
        _single_phase(st, out)
    finally:
        st.unspecify_phase()
    out["X"] = float(Q)
    return out

def _query_side(st, sat, T, rho, liquid_side):
    """Single-phase state on a known side of the dome (phase forced below
    Tcrit so points exactly on the boundary stay single-phase)."""
    out = {}
    if not (Ttrip <= T <= _T_EOS_MAX) or not (rho > 0.0):
        return out
    out["X"] = _quality(sat, T, rho)
    try:
        if T < Tcrit:
            st.specify_phase(CoolProp.iphase_liquid if liquid_side
                             else CoolProp.iphase_gas)
        st.update(CoolProp.DmassT_INPUTS, rho, T)
        _common(st, out)
        _single_phase(st, out)
    finally:
        st.unspecify_phase()
    return out


def _pack_arr(out):
    return np.array([out.get(name, math.nan) for name in FIELDS],
                    dtype=np.float64)

def _pack(out):
    return N2O_Properties.from_arr(_pack_arr(out))

def _safe(fn, *args):
    try:
        return _pack(fn(*args))
    except ValueError: # CoolProp raises ValueError on failure.
        return N2O_Properties.NAN

def _safe_arr(fn, *args):
    try:
        return _pack_arr(fn(*args))
    except ValueError:
        return np.full(len(FIELDS), math.nan)


# -------------------------------------------------------------------------------
# Uncached public queries (return arrays ordered as `FIELDS`).

def query(T, rho):
    return _safe_arr(_query, _st, _st_sat, float(T), float(rho))

def query_sat(T, Q):
    return _safe_arr(_query_sat, _st, float(T), int(Q))

def query_dense(T, psi):
    T = float(T)
    try:
        rho = rho_from_dense_psi(T, float(psi))
    except ValueError:
        return np.full(len(FIELDS), math.nan)
    return _safe_arr(_query_side, _st, _st_sat, T, rho, True)

def query_light(T, psi):
    T = float(T)
    try:
        rho = rho_from_light_psi(T, float(psi))
    except ValueError:
        return np.full(len(FIELDS), math.nan)
    return _safe_arr(_query_side, _st, _st_sat, T, rho, False)


# -------------------------------------------------------------------------------
# Saturation temperature and the saturated-inlet HEM choke.

def Tsat(P):
    """Saturation temperature [K] at pressure P [Pa], Ptrip <= P < Pcrit."""
    _st_sat.unspecify_phase()
    _st_sat.update(CoolProp.PQ_INPUTS, float(P), 0.0)
    return _st_sat.T()

def _sonic_resid(s_u, h_u, T):
    """V^2 - a_eq^2 for an isentropic HEM expansion from (s_u, h_u) to
    saturation temperature T, with
        V^2 = 2*(h_u - h),  a_eq^2 = -v^2 * Psat'(T) / (dv/dT)_s,
        v = v_l + X*v_fg,  X = (s_u - s_l)/s_fg,  X' = -(s_l' + X*s_fg')/s_fg,
        (dv/dT)_s = v_l' + X*v_fg' + v_fg*X'
    (' = d/dT along saturation). Its root is the HEM choke."""
    sides = []
    for Q in (0.0, 1.0):
        _st_sat.unspecify_phase()
        _st_sat.update(CoolProp.QT_INPUTS, Q, T)
        fsd = _st_sat.first_saturation_deriv
        rho = _st_sat.rhomass()
        sides.append(dict(
            v=1.0/rho, dv=-fsd(CoolProp.iDmass, CoolProp.iT)/rho**2,
            s=_st_sat.smass(), ds=fsd(CoolProp.iSmass, CoolProp.iT),
            h=_st_sat.hmass(), dP=fsd(CoolProp.iP, CoolProp.iT),
        ))
    l, v = sides
    sfg = v["s"] - l["s"]
    X = (s_u - l["s"]) / sfg
    dX = -(l["ds"] + X*(v["ds"] - l["ds"])) / sfg
    vm = l["v"] + X*(v["v"] - l["v"])
    dvm = l["dv"] + X*(v["dv"] - l["dv"]) + (v["v"] - l["v"])*dX
    a2 = -vm*vm * l["dP"] / dvm
    h = l["h"] + X*(v["h"] - l["h"])
    return 2.0*(h_u - h) - a2

def hem_choke_T_sat(T_u):
    """HEM choke exit (saturation) temperature [K] for a saturated-liquid inlet
    at T_u: where the isentropic flow speed reaches the equilibrium sound
    speed. The same point as the maximum mass flux, but found as a root (well-
    conditioned) rather than as a maximum (flat, so a search there wanders over
    ~1 K and drags the exit quality with it). NaN if no root is bracketed."""
    T_u = float(T_u)
    _st_sat.unspecify_phase()
    _st_sat.update(CoolProp.QT_INPUTS, 0.0, T_u)
    s_u = _st_sat.smass()
    h_u = _st_sat.hmass()
    # Illinois false position; f > 0 far below T_u (V large, a_eq small),
    # f < 0 just below T_u (V -> 0).
    lo, hi = Ttrip + 0.01, T_u - 1e-4
    flo, fhi = _sonic_resid(s_u, h_u, lo), _sonic_resid(s_u, h_u, hi)
    if not (flo > 0.0 > fhi):
        return math.nan
    side = 0
    for _ in range(200):
        T = (lo*fhi - hi*flo) / (fhi - flo)
        f = _sonic_resid(s_u, h_u, T)
        if f > 0.0:
            lo, flo = T, f
            if side == +1:
                fhi *= 0.5
            side = +1
        else:
            hi, fhi = T, f
            if side == -1:
                flo *= 0.5
            side = -1
        if hi - lo < 1e-9:
            break
    return 0.5*(lo + hi)


# -------------------------------------------------------------------------------
# Cached queries.

@cacher.cls_Fetcher
def N2O():
    st = _new_state()
    sat = _new_state()
    spacing = np.array([0.1, 0.25]) # T [K], rho [kg/m^3]
    def f(T, rho):
        return _safe(_query, st, sat, float(T), float(rho))
    return cacher.Cacher(f"N2O-{_VERSION}", spacing, f,
            N2O_Properties.from_arr, N2O_Properties.to_arr, N2O_Properties.NAN)


def _sat_fetcher(name, Q):
    st = _new_state()
    spacing = np.array([0.1]) # T [K]
    def f(T):
        # Clamp so grid nodes just outside the valid range return the boundary
        # value rather than NaN (which would poison interpolation).
        T = min(max(float(T), Ttrip), Tsat_hi)
        return _safe(_query_sat, st, T, Q)
    return cacher.Cacher(f"{name}-{_VERSION}", spacing, f,
            N2O_Properties.from_arr, N2O_Properties.to_arr, N2O_Properties.NAN)

@cacher.cls_Fetcher
def N2Osatliq():
    return _sat_fetcher("N2Osatliq", 0)

@cacher.cls_Fetcher
def N2Osatvap():
    return _sat_fetcher("N2Osatvap", 1)
