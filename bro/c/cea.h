#pragma once
#include "br.h"

// NASA-CEA outputs, approximated by lookup tables statically generated (see
// approximator/cea.py and approximator/run.py).
//
// Each propellant combination has its own tables over its own domain. Out-of-
// domain inputs, and fields not tabulated for a propellant, assert.
//
// Stations are selected by the area ratio `AR`:
//   AR == CEA_AR_stag (0)  chamber (stagnation).
//   AR == CEA_AR_tht  (1)  throat.
//   AR >  1                supersonic, A/A_throat = AR.
//   0 < AR < 1             subsonic,   A/A_throat = 1/AR.

enum {
    CEA_N2O_PARAFFIN,
    CEA_PROP_COUNT,
};

#define CEA_AR_stag (0.0)
#define CEA_AR_tht (1.0)

f64 cea_T(i32 prop, f64 P0_cc, f64 ofr, f64 AR);
f64 cea_P(i32 prop, f64 P0_cc, f64 ofr, f64 AR);
f64 cea_rho(i32 prop, f64 P0_cc, f64 ofr, f64 AR);
f64 cea_M(i32 prop, f64 P0_cc, f64 ofr, f64 AR);
f64 cea_a(i32 prop, f64 P0_cc, f64 ofr, f64 AR);
f64 cea_gamma(i32 prop, f64 P0_cc, f64 ofr, f64 AR);
f64 cea_cp(i32 prop, f64 P0_cc, f64 ofr, f64 AR);
f64 cea_mu(i32 prop, f64 P0_cc, f64 ofr, f64 AR);
f64 cea_Pr(i32 prop, f64 P0_cc, f64 ofr, f64 AR);
f64 cea_Mw(i32 prop, f64 P0_cc, f64 ofr, f64 AR);
f64 cea_Ivac(i32 prop, f64 P0_cc, f64 ofr, f64 AR); // AR>=1.
