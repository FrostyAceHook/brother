#pragma once
#include "br.h"


// Geomtric formulas.

f64 circ_diam(f64 A);
f64 circ_area(f64 D);
f64 circ_perim(f64 D);
f64 cyl_vol(f64 L, f64 D);
f64 cyl_L(f64 V, f64 D);
f64 tube_vol(f64 L, f64 ID, f64 OD);


// Structual pipe/tube.

typedef struct Pipe {
    f64 L;  // cannot be nan.
    f64 ID; // may be nan.
    f64 OD; // may be nan.
} Pipe;
f64 pipe_th(const Pipe* p);
f64 pipe_V(const Pipe* p);
void pipe_set_th(Pipe* p, f64 th);
void pipe_set_th_for_P(Pipe* p, f64 P, f64 Ys, f64 sf);
f64 pipe_sf_for_P(const Pipe* p, f64 P, f64 Ys);
