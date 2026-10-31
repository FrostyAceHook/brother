#include "n2o.h"

#include "assertion.h"
#include "lut.h"
#include "maths.h"

// Generated tables. Each defines its axes (lutAxis) and one named f32 array per
// field (N2O_<group>_<field>_tbl):
//   N2O_sat.i    N2O_SAT_AXIS (T, cbrt about Tcrit)
//   N2O_bound.i  N2O_BOUND_AXIS (T, linear)
//   N2O_dense.i  N2O_DENSE_T_AXIS (T, sqrt about Tcrit), N2O_DENSE_PSI_AXIS
//   N2O_light.i  N2O_LIGHT_T_AXIS (T, sqrt about Tcrit), N2O_LIGHT_PSI_AXIS
//   N2O_tsat.i   N2O_TSAT_AXIS (P, log)
//   N2O_hem.i    N2O_HEM_AXIS (T, cbrt about Tcrit)
#include "tbl/N2O_sat.i"
#include "tbl/N2O_bound.i"
#include "tbl/N2O_dense.i"
#include "tbl/N2O_light.i"
#include "tbl/N2O_tsat.i"
#include "tbl/N2O_hem.i"



// =========================================================================== //
// = SATURATED =============================================================== //
// =========================================================================== //

#define N2O_DEFINE_SAT_(side, name)                                         \
    f64 N2O_##side##_##name(f64 T) {                                        \
        assert(N2O_Ttrip <= T && T < N2O_Tcrit,                             \
                "saturated N2O input oob: T=%g", T);                        \
        /* Untabulated sliver [Tsat_hi, Tcrit) takes the last value. */     \
        T = min(T, N2O_SAT_AXIS.hi);                                        \
        f64 X = lut1(N2O_##side##_##name##_tbl, &N2O_SAT_AXIS, T);          \
        assert(notnan(X), "saturated N2O lookup nan: T=%g", T);             \
        return X;                                                           \
    }
#define N2O_DEFINE_SATLIQ_(name) N2O_DEFINE_SAT_(satliq, name)
#define N2O_DEFINE_SATVAP_(name) N2O_DEFINE_SAT_(satvap, name)
N2O_TBL_SAT_FIELDS(N2O_DEFINE_SATLIQ_)
N2O_TBL_SAT_FIELDS(N2O_DEFINE_SATVAP_)
#undef N2O_DEFINE_SATLIQ_
#undef N2O_DEFINE_SATVAP_
#undef N2O_DEFINE_SAT_


f64 N2O_Psat(f64 T) {
    return N2O_satliq_P(T);
}

f64 N2O_Tsat(f64 P) {
    assert(N2O_Ptrip <= P && P < N2O_Pcrit, "N2O Tsat input oob: P=%g", P);
    if (P <= N2O_TSAT_AXIS.hi)
        return lut1(N2O_tsat_T_tbl, &N2O_TSAT_AXIS, P);
    // Untabulated sliver [Psat_hi, Pcrit): blend to the critical point.
    f64 T_hi = (f64)N2O_tsat_T_tbl[N2O_TSAT_AXIS.n - 1];
    f64 x = (P - N2O_TSAT_AXIS.hi) / (N2O_Pcrit - N2O_TSAT_AXIS.hi);
    return lerp(T_hi, N2O_Tcrit, x);
}


f64 N2O_satliq_T_chokeHEM(f64 T_u) {
    if (!lut_inside(&N2O_HEM_AXIS, T_u, 0.0))
        return NAN;
    return lut1(N2O_hem_Tchoke_tbl, &N2O_HEM_AXIS, T_u);
}


// =========================================================================== //
// = REGION LOCATION ========================================================= //
// =========================================================================== //

// Dome boundary densities: saturated below Tcrit, rhocrit above, linearly
// blended over the untabulated sliver [Tsat_hi, Tcrit] (continuous at Tcrit).
static f64 N2O_rho_split_(const f32* tbl, f64 T) {
    if (T >= N2O_Tcrit)
        return N2O_rhocrit;
    if (T <= N2O_SAT_AXIS.hi)
        return lut1(tbl, &N2O_SAT_AXIS, T);
    f64 rho_hi = (f64)tbl[N2O_SAT_AXIS.n - 1];
    f64 x = (T - N2O_SAT_AXIS.hi) / (N2O_Tcrit - N2O_SAT_AXIS.hi);
    return lerp(rho_hi, N2O_rhocrit, x);
}
static f64 N2O_rho_split_l_(f64 T) {
    return N2O_rho_split_(N2O_satliq_rho_tbl, T);
}
static f64 N2O_rho_split_v_(f64 T) {
    return N2O_rho_split_(N2O_satvap_rho_tbl, T);
}
static f64 N2O_rho_dense_max_(f64 T) {
    return lut1(N2O_bound_rho_dense_max_tbl, &N2O_BOUND_AXIS, T);
}

enum {
    N2O_REGION_DENSE_,
    N2O_REGION_LIGHT_,
    N2O_REGION_DOME_,
};
typedef struct N2O_Loc_ {
    i32 region;
    f64 psi; // table coordinate (dense/light only)
    f64 rl;  // dense/dome boundary density
    f64 rv;  // light/dome boundary density
} N2O_Loc_;

static N2O_Loc_ N2O_locate_(f64 T, f64 rho) {
    assert(N2O_Tmin <= T && T <= N2O_Tmax, "N2O input oob: T=%g", T);
    N2O_Loc_ l;
    l.rl = N2O_rho_split_l_(T);
    l.rv = N2O_rho_split_v_(T);
    if (rho >= l.rl) {
        f64 rho_max = N2O_rho_dense_max_(T);
        assert(rho <= rho_max*(1.0 + 1e-6),
                "N2O input oob (above Pmax): T=%g, rho=%g", T, rho);
        l.region = N2O_REGION_DENSE_;
        l.psi = (rho - l.rl) / (rho_max - l.rl);
    } else if (rho <= l.rv) {
        assert(rho >= N2O_rhomin*(1.0 - 1e-6),
                "N2O input oob (below rhomin): T=%g, rho=%g", T, rho);
        l.region = N2O_REGION_LIGHT_;
        l.psi = log2(rho / N2O_rhomin) / log2(l.rv / N2O_rhomin);
    } else {
        l.region = N2O_REGION_DOME_;
        l.psi = NAN;
    }
    return l;
}

static f64 N2O_quality_(f64 rho, N2O_Loc_ l) {
    return (1.0/rho - 1.0/l.rl) / (1.0/l.rv - 1.0/l.rl);
}



// =========================================================================== //
// = GENERAL ================================================================= //
// =========================================================================== //

// Inside-the-dome values (lever rule on the saturated tables). Every field in
// N2O_TBL_2D_FIELDS needs one; a new tabulated field without one is a compile
// error rather than a silently wrong value.
static f64 N2O_dome_P_(f64 T, f64 rho, N2O_Loc_ l) {
    (void)rho; (void)l;
    return N2O_satliq_P(T);
}
#define N2O_DOME_LEVER_(name)                                               \
    static f64 N2O_dome_##name##_(f64 T, f64 rho, N2O_Loc_ l) {             \
        f64 X = N2O_quality_(rho, l);                                       \
        return lerp(N2O_satliq_##name(T), N2O_satvap_##name(T), X);         \
    }
N2O_DOME_LEVER_(s)
N2O_DOME_LEVER_(h)
N2O_DOME_LEVER_(u)
#undef N2O_DOME_LEVER_
static f64 N2O_dome_Z_(f64 T, f64 rho, N2O_Loc_ l) {
    (void)l;
    return N2O_satliq_P(T) / (rho * N2O_R * T);
}
// Undefined for an equilibrium mixture (cp is infinite; the mixture sound speed
// is a different quantity).
#define N2O_DOME_NAN_(name)                                                 \
    static f64 N2O_dome_##name##_(f64 T, f64 rho, N2O_Loc_ l) {             \
        (void)T; (void)rho; (void)l;                                        \
        return NAN;                                                         \
    }
N2O_DOME_NAN_(cp)
N2O_DOME_NAN_(a)
#undef N2O_DOME_NAN_

// Equilibrium-mixture partials, from the saturated tables and their slopes
// along the saturation curve ('). With v = 1/rho, X = (v - v_l)/v_fg and
// phi = phi_l + X*phi_fg for phi in {u, h, s}:
//   (dphi/drho)_T = -phi_fg / (rho^2 * v_fg)
//   (dphi/dT)_rho = phi_l' + X*phi_fg' + phi_fg * dX/dT,
//       dX/dT = -(v_l' + X*v_fg') / v_fg
//   (dP/dT)_rho = dPsat/dT,  (dP/drho)_T = 0,  cv = (du/dT)_rho
typedef struct N2O_Mix_ {
    f64 X;    // quality
    f64 vfg;  // v_v - v_l
    f64 dXdT; // (dX/dT)_rho
} N2O_Mix_;
static N2O_Mix_ N2O_mix_(f64 T, f64 rho, N2O_Loc_ l) {
    f64 vl = 1.0/l.rl;
    f64 vv = 1.0/l.rv;
    f64 vfg = vv - vl;
    f64 X = (1.0/rho - vl) / vfg;
    f64 dvl = -N2O_satliq_drhosat_dT(T) * vl*vl;
    f64 dvv = -N2O_satvap_drhosat_dT(T) * vv*vv;
    return (N2O_Mix_){ X, vfg, -(dvl + X*(dvv - dvl)) / vfg };
}
static f64 N2O_mix_dT_rho_(N2O_Mix_ m, f64 fl, f64 fv, f64 dfl, f64 dfv) {
    return dfl + m.X*(dfv - dfl) + (fv - fl)*m.dXdT;
}
static f64 N2O_mix_drho_T_(N2O_Mix_ m, f64 rho, f64 fl, f64 fv) {
    return -(fv - fl) / (rho*rho*m.vfg);
}

static f64 N2O_dome_cv_(f64 T, f64 rho, N2O_Loc_ l) {
    return N2O_mix_dT_rho_(N2O_mix_(T, rho, l),
            N2O_satliq_u(T), N2O_satvap_u(T),
            N2O_satliq_dusat_dT(T), N2O_satvap_dusat_dT(T));
}
static f64 N2O_dome_dh_dT_rho_(f64 T, f64 rho, N2O_Loc_ l) {
    return N2O_mix_dT_rho_(N2O_mix_(T, rho, l),
            N2O_satliq_h(T), N2O_satvap_h(T),
            N2O_satliq_dhsat_dT(T), N2O_satvap_dhsat_dT(T));
}
static f64 N2O_dome_dP_dT_rho_(f64 T, f64 rho, N2O_Loc_ l) {
    (void)rho; (void)l;
    return N2O_satliq_dPsat_dT(T);
}
static f64 N2O_dome_dP_drho_T_(f64 T, f64 rho, N2O_Loc_ l) {
    (void)T; (void)rho; (void)l;
    return 0.0;
}
#define N2O_DOME_DRHO_(fld)                                                 \
    static f64 N2O_dome_d##fld##_drho_T_(f64 T, f64 rho, N2O_Loc_ l) {      \
        return N2O_mix_drho_T_(N2O_mix_(T, rho, l), rho,                    \
                N2O_satliq_##fld(T), N2O_satvap_##fld(T));                  \
    }
N2O_DOME_DRHO_(u)
N2O_DOME_DRHO_(h)
N2O_DOME_DRHO_(s)
#undef N2O_DOME_DRHO_

#define N2O_DEFINE_GENERAL_(name)                                           \
    f64 N2O_##name(f64 T, f64 rho) {                                        \
        N2O_Loc_ l = N2O_locate_(T, rho);                                   \
        switch (l.region) {                                                 \
            case N2O_REGION_DENSE_:                                         \
                return lut2(N2O_dense_##name##_tbl,                         \
                        &N2O_DENSE_T_AXIS, &N2O_DENSE_PSI_AXIS, T, l.psi);  \
            case N2O_REGION_LIGHT_:                                         \
                return lut2(N2O_light_##name##_tbl,                         \
                        &N2O_LIGHT_T_AXIS, &N2O_LIGHT_PSI_AXIS, T, l.psi);  \
            default:                                                        \
                return N2O_dome_##name##_(T, rho, l);                       \
        }                                                                   \
    }
N2O_TBL_2D_FIELDS(N2O_DEFINE_GENERAL_)
#undef N2O_DEFINE_GENERAL_


// Derived rather than tabulated: Y stays exactly consistent with cp and cv, and
// IPLE (~1/P along psi at low T) interpolates badly while a and P don't.
f64 N2O_Y(f64 T, f64 rho) {
    return N2O_cp(T, rho) / N2O_cv(T, rho); // NaN inside the dome (cp is).
}
f64 N2O_IPLE(f64 T, f64 rho) {
    f64 a = N2O_a(T, rho); // NaN inside the dome.
    return rho * a * a / N2O_P(T, rho);
}

f64 N2O_X(f64 T, f64 rho) {
    if (T >= N2O_Tcrit)
        return NAN;
    N2O_Loc_ l = N2O_locate_(T, rho);
    return N2O_quality_(rho, l);
}