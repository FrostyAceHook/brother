#pragma once
#include "br.h"

#include "relations.h"

// Ambient properties from the International Standard Atmosphere, approximated
// via rational polynomials by 'bro/approximator.py'.

#define STANDARD_GRAVITY (9.80665) // straight off the dome.

// These properties are assumed, as "standard conditions".
#define AMB_STANCOND_Mw (28.9647e-3)
#define AMB_STANCOND_R (GAS_CONSTANT / AMB_STANCOND_Mw)
#define AMB_STANCOND_T (293.15)
#define AMB_STANCOND_cp (1004.7)
#define AMB_STANCOND_cv (AMB_STANCOND_cp - AMB_STANCOND_R)

f64 amb_T(f64 alt);
f64 amb_P(f64 alt);
f64 amb_rho(f64 alt);
f64 amb_a(f64 alt);
f64 amb_g(f64 alt);
