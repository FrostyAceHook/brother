#pragma once
#include "br.h"

#include "relations.h"

// Ambient properties from the International Standard Atmosphere, approximated
// via rational polynomials by 'bro/approximator.py'.
//
// These functions accept any altitude, but default to assuming no atmosphere
// above 60km and a return as-if at 0km when given below 0km. Note that all
// inputs and outputs are base si.

#define STANDARD_GRAVITY (9.80665) // straight off the dome.

#define AMB_softmax_alt (60e3)
#define AMB_softmin_alt (0.0)

// These properties are assumed, as "standard conditions".
#define AMB_STANCOND_Mw (28.9647e-3)
#define AMB_STANCOND_R (GAS_CONSTANT / AMB_STANCOND_Mw)
#define AMB_STANCOND_T (293.15)
#define AMB_STANCOND_cp (1005.0)
#define AMB_STANCOND_cv (AMB_STANCOND_cp - AMB_STANCOND_R)

f64 amb_P(f64 alt);
f64 amb_rho(f64 alt);
f64 amb_a(f64 alt);
f64 amb_g(f64 alt);
