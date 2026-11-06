#pragma once
#include "br.h"

// Grain geometry + regression functions.

// (Empirical) regression rate.
f64 grain_drr(f64 G_o);

typedef struct Grain {
    f64 L;
    f64 ID;
    f64 OD;
    f64 rr;
} Grain;

f64 grain_max_rr(const Grain* grain);
f64 grain_least_th(const Grain* grain);

f64 grain_Aneg(const Grain* grain);
f64 grain_dAneg(const Grain* grain, f64 Aneg, f64 drr);
f64 grain_V(const Grain* grain, f64 Aneg);
f64 grain_dV(const Grain* grain, f64 Aneg, f64 dAneg);
