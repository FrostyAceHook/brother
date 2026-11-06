#include "grain.h"

#include "assertion.h"
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
