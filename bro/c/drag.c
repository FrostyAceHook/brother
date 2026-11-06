#include "drag.h"

#include "assertion.h"
#include "maths.h"

#include "ambient.h"


#define DRAG_TBL_Mmin (0.05)
#define DRAG_TBL_Mmax (0.9 + 0.05*5)
static const f64 DRAG_TBL_M[] = {
    // Drag curve from structures for Ody V2 2026/11/06, clean model.
    0.391603058,
    0.357034263,
    0.340667007,
    0.333157043,
    0.325230972,
    0.302864458,
    0.325845682,
    0.321057597,
    0.319039084,
    0.316793308,
    0.314393559,
    0.313004018,
    0.309672971,
    0.308727786,
    0.307576269,
    0.306617945,
    0.307402859,
    0.308465884,

    /* MADE UP TO NOT FAULT SIM AT EXTREMES: */
    0.33,
    0.36,
    0.40,
    0.42,
    0.42, // flat forever.
};


f64 drag_CD(f64 M) {
    // Use lookup table.
    assert(M >= 0.0, "M=%f", M);
    f64 x = invlerp(DRAG_TBL_Mmin, DRAG_TBL_Mmax, M);
    x *= numel(DRAG_TBL_M) - 1;
    i32 i = min(max((i32)x, 0), numel(DRAG_TBL_M) - 2);
    f64 t = x - i;
    return lerp(DRAG_TBL_M[i], DRAG_TBL_M[i + 1], t);
}

f64 drag_force(f64 alt, f64 vel, f64 A) {
    // Mach-dependant drag.
    f64 M = abs(vel) / amb_a(alt);
    f64 CD = drag_CD(M);
    f64 Drag = CD * A * 0.5 * amb_rho(alt) * vel*vel;
    // "did you remember to make the drag negative?"
    if (vel > 0.0)
        Drag = -Drag;
    return Drag;
}
