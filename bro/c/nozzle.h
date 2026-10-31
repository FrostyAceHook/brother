#pragma once
#include "br.h"

// Rao nozzle.

// Nozzle diverging and exit angles extrapolated from data in nasa paper:
// https://ntrs.nasa.gov/citations/19770009165

#define NZL_min_AEAT (2.0)
#define NZL_max_AEAT (12.0)
#define NZL_min_NLF (0.6)
#define NZL_max_NLF (1.0)

f64 nzl_theta_div(f64 AEAT, f64 NLF);
f64 nzl_theta_exit(f64 AEAT, f64 NLF);