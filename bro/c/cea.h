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

// TODO: gen these from approximator. currently hardcoded.
#define CEA_min_P (0.5e5)
#define CEA_max_P (90e5)
#define CEA_min_ofr (1.0)
#define CEA_max_ofr (15.0)
#define CEA_min_AR (1.2)
#define CEA_max_AR (12.0)

#define CEA_AR_stag (0.0)
#define CEA_AR_tht (1.0)

f64 cea_T(f64 P0_cc, f64 ofr, f64 AR);
f64 cea_P(f64 P0_cc, f64 ofr, f64 AR);
f64 cea_rho(f64 P0_cc, f64 ofr, f64 AR);
f64 cea_M(f64 P0_cc, f64 ofr, f64 AR);
f64 cea_a(f64 P0_cc, f64 ofr, f64 AR);
f64 cea_Y(f64 P0_cc, f64 ofr, f64 AR);
f64 cea_cp(f64 P0_cc, f64 ofr, f64 AR);
f64 cea_mu(f64 P0_cc, f64 ofr, f64 AR);
f64 cea_Pr(f64 P0_cc, f64 ofr, f64 AR);
f64 cea_Mw(f64 P0_cc, f64 ofr, f64 AR);
f64 cea_Ivac(f64 P0_cc, f64 ofr, f64 AR); // AR>=1.
