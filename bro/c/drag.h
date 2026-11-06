#pragma once
#include "br.h"

// Rocket drag.

f64 drag_CD(f64 alt);
f64 drag_force(f64 alt, f64 vel, f64 A);
