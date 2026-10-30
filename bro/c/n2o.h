#pragma once
#include "br.h"

// Nitrous oxide properties, approximated by lookup tables generated from
// CoolProp (see approximator/n2o.py and approximator/run.py).
//
// Tables are split along the saturation dome so no table cell straddles it:
//   - saturated liquid/vapour: 1D in T, over [N2O_Ttrip, N2O_Tsat_hi].
//   - dense (liquid, dense supercritical): 2D in (T, psi), from the saturated
//     liquid density (or rhocrit above Tcrit) up to the density at N2O_Pmax.
//   - light (vapour, light supercritical): 2D in (T, psi), from N2O_rhomin up
//     to the saturated vapour density (or rhocrit above Tcrit).
//   - inside the dome: nothing tabulated; built from the saturated tables with
//     the lever rule.
//
// General (T, rho) queries require T in [N2O_Tmin, N2O_Tmax] and rho in
// [N2O_rhomin, rho(T, N2O_Pmax)]. Saturated queries require T in
// [N2O_Ttrip, N2O_Tcrit). Out-of-domain inputs assert.

#include "tbl/N2O_consts.i"


// ========================= //
//          GENERAL          //
// ========================= //

// Always defined. Inside the dome, P is the saturation pressure, s/h/u/Z are
// lever-rule mixture values, and cv and the partial derivatives are the
// equilibrium-mixture values ((dP/drho)_T = 0, (dP/dT)_rho = dPsat/dT).
f64 N2O_P(f64 T, f64 rho);
f64 N2O_s(f64 T, f64 rho);
f64 N2O_h(f64 T, f64 rho);
f64 N2O_u(f64 T, f64 rho);
f64 N2O_Z(f64 T, f64 rho);
f64 N2O_X(f64 T, f64 rho); // Vapour quality [-] (<0 compressed liquid,
                           //   >1 superheated vapour, NaN supercritical)
f64 N2O_cv(f64 T, f64 rho);
f64 N2O_dP_dT_rho(f64 T, f64 rho);
f64 N2O_dP_drho_T(f64 T, f64 rho);
f64 N2O_du_drho_T(f64 T, f64 rho);
f64 N2O_dh_dT_rho(f64 T, f64 rho);
f64 N2O_dh_drho_T(f64 T, f64 rho);
f64 N2O_ds_drho_T(f64 T, f64 rho);

// Single-phase only (NaN inside the dome, where they're undefined for an
// equilibrium mixture: cp is infinite, and the mixture sound speed is a
// different quantity):
f64 N2O_cp(f64 T, f64 rho);
f64 N2O_Y(f64 T, f64 rho);
f64 N2O_a(f64 T, f64 rho);
f64 N2O_IPLE(f64 T, f64 rho);


// ========================= //
//         SATURATED         //
// ========================= //

// Saturation pressure/temperature:
f64 N2O_Psat(f64 T);
f64 N2O_Tsat(f64 P);

// Saturated liquid/vapour at T, as the one-sided single-phase limit (so every
// property is defined), requiring T in [N2O_Ttrip, N2O_Tcrit). The `*sat_dT`
// fields are total derivatives along the saturation curve, not partials.

// Saturated liquid:
f64 N2O_satliq_P(f64 T);
f64 N2O_satliq_rho(f64 T);
f64 N2O_satliq_s(f64 T);
f64 N2O_satliq_h(f64 T);
f64 N2O_satliq_u(f64 T);
f64 N2O_satliq_Z(f64 T);
f64 N2O_satliq_cp(f64 T);
f64 N2O_satliq_cv(f64 T);
f64 N2O_satliq_Y(f64 T);
f64 N2O_satliq_a(f64 T);
f64 N2O_satliq_IPLE(f64 T);
f64 N2O_satliq_dP_dT_rho(f64 T);
f64 N2O_satliq_dP_drho_T(f64 T);
f64 N2O_satliq_du_drho_T(f64 T);
f64 N2O_satliq_dh_dT_rho(f64 T);
f64 N2O_satliq_dh_drho_T(f64 T);
f64 N2O_satliq_ds_drho_T(f64 T);
f64 N2O_satliq_dPsat_dT(f64 T);
f64 N2O_satliq_drhosat_dT(f64 T);
f64 N2O_satliq_dssat_dT(f64 T);
f64 N2O_satliq_dhsat_dT(f64 T);
f64 N2O_satliq_dusat_dT(f64 T);

// Saturated vapour:
f64 N2O_satvap_P(f64 T);
f64 N2O_satvap_rho(f64 T);
f64 N2O_satvap_s(f64 T);
f64 N2O_satvap_h(f64 T);
f64 N2O_satvap_u(f64 T);
f64 N2O_satvap_Z(f64 T);
f64 N2O_satvap_cp(f64 T);
f64 N2O_satvap_cv(f64 T);
f64 N2O_satvap_Y(f64 T);
f64 N2O_satvap_a(f64 T);
f64 N2O_satvap_IPLE(f64 T);
f64 N2O_satvap_dP_dT_rho(f64 T);
f64 N2O_satvap_dP_drho_T(f64 T);
f64 N2O_satvap_du_drho_T(f64 T);
f64 N2O_satvap_dh_dT_rho(f64 T);
f64 N2O_satvap_dh_drho_T(f64 T);
f64 N2O_satvap_ds_drho_T(f64 T);
f64 N2O_satvap_dPsat_dT(f64 T);
f64 N2O_satvap_drhosat_dT(f64 T);
f64 N2O_satvap_dssat_dT(f64 T);
f64 N2O_satvap_dhsat_dT(f64 T);
f64 N2O_satvap_dusat_dT(f64 T);
