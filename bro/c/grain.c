#include "grain.h"

#include "geometry.h"
#include "maths.h"


f64 grain_drr(f64 G_o) {
    // Paraffin+N2O regression rate constants from Hybrid Rocket Propulsion
    // Handbook, Karp & Jens.
    return 1.55e-4 * sqrt(G_o);
}


f64 grain_max_rr(const Grain* grain) {
    return 0.5*(grain->OD - grain->ID);
}
f64 grain_least_th(const Grain* grain) {
    f64 th = 0.5*(grain->OD - grain->ID) - grain->rr;
    return max(th, 0);
}


f64 grain_Aneg(const Grain* grain) {
    f64 ID = grain->ID + 2*grain->rr;
    ID = min(ID, grain->OD);
    return circ_area(ID);
}

f64 grain_dAneg(const Grain* grain, f64 Aneg, f64 drr) {
    (void)Aneg;
    // A = pi/4 * ID^2
    // d/dt (A) = d/dt (pi/4 * ID^2)
    // dA = pi/4 * d/dt (ID^2)
    // dA = pi/4 * 2 * ID * dID
    // dA = pi/4 * 2 * ID * (2 * drr)
    // dA = pi * ID * drr
    f64 ID = grain->ID + 2*grain->rr;
    if (ID >= grain->OD)
        return 0.0;
    return PI * ID * drr;
}

f64 grain_V(const Grain* grain, f64 Aneg) {
    return grain->L * (circ_area(grain->OD) - Aneg);
}

f64 grain_dV(const Grain* grain, f64 Aneg, f64 dAneg) {
    (void)Aneg;
    return grain->L * (-dAneg);
}


void grain_radial_contour(const Grain* grain, i32 count, f64* rstr xs,
        f64* rstr ys) {
    assert(count >= 50, "count=%d", count);
    assert(xs, "null pointer");
    assert(ys, "null pointer");

    f64 Ir = 0.5*grain->ID + grain->rr;
    f64 Or = 0.5*grain->OD;

    Ir = min(Ir, Or);

    i32 i = 0;

    // Half the points for the interior, half for the exterior.
    i32 count_i = count / 2;
    i32 count_nan = 1;
    i32 count_o = count - i - count_i - count_nan;

    i32 start;

    // Inner circle me.
    start = i;
    while (i < start + count_i) {
        f64 a = TWOPI * (i - start) / (f64)(count_i - 1);
        f64 x = Ir * cos(a);
        f64 y = Ir * sin(a);
        xs[i] = x;
        ys[i] = y;
        ++i;
    }

    // Nan me (to disconnect the two thangs).
    xs[i] = NAN;
    ys[i] = NAN;
    ++i;

    // Outer circle me.
    start = i;
    while (i < start + count_o) {
        f64 a = TWOPI * (i - start) / (f64)(count_o - 1);
        f64 x = Or * cos(a);
        f64 y = Or * sin(a);
        xs[i] = x;
        ys[i] = y;
        ++i;
    }

    assert(i == count, "i=%d, count=%d", i, count);
    return;
}


void grain_axial_contour(const Grain* grain, f64 x0, i32 count, f64* rstr xs,
        f64* rstr ys) {
    assert(xs, "null pointer");
    assert(ys, "null pointer");

    assert(count >= 5, "count=%d", count);
    while (count > 5) {
        *xs++ = NAN;
        *ys++ = NAN;
        --count;
    }
    f64 Ir = 0.5*grain->ID + grain->rr;
    f64 Or = 0.5*grain->OD;
    Ir = min(Ir, Or);

    xs[0] = x0;
    ys[0] = Or;
    xs[1] = x0 + grain->L;
    ys[1] = Or;
    xs[2] = x0 + grain->L;
    ys[2] = Ir;
    xs[3] = x0;
    ys[3] = Ir;
    xs[4] = x0;
    ys[4] = Or;
    return;
}
