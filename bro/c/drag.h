#pragma once
#include "br.h"

// Rocket drag.

f64 rocket_CD(f64 alt);
f64 rocket_Drag(f64 alt, f64 vel, f64 A);
