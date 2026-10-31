#include "drag.h"

#include "assertion.h"
#include "maths.h"

#include "ambient.h"


// Rocket properties from structures for current cycle.
// https://www.desmos.com/calculator/00d7pcn1ec
#define ROCKET_CD_max_M (1.25)

f64 rocket_CD(f64 M) {
    assert(M >= 0.0, "M=%f", M);
    assert(M <= ROCKET_CD_max_M, "M=%f", M);
  #if 1
    f64 n4 = -2.68015197546;
    f64 n3 = +4.32516313836;
    f64 n2 = -4.05555843745;
    f64 n1 = +1.31009964679;
    f64 n0 = +0.130126875872;
    f64 d3 = +6.45904203456;
    f64 d2 = -11.540104789;
    f64 d1 = +4.86614101695;
    f64 d0 = +0.289084735643;
    return (n0 + M*(n1 + M*(n2 + M*(n3 + M*(n4 + /* 1.0* */M)))))
         / (d0 + M*(d1 + M*(d2 + M*d3)));
  #else
    // https://www.desmos.com/calculator/f2mvcvsrd3
    f64 n0 = -16.2988133908;
    f64 n1 = -48.1659386797;
    f64 n2 = +143.674493616;
    f64 n3 = -80.1863491945;
    f64 d0 = -35.0797586734;
    f64 d1 = -193.441052967;
    f64 d2 = +540.844019766;
    f64 d3 = -370.66858712;
    f64 d4 = +56.1759212259;
    return (n0 + M*(n1 + M*(n2 + M*n3)))
         / (d0 + M*(d1 + M*(d2 + M*(d3 + M*d4))));
  #endif
}

f64 rocket_Drag(f64 alt, f64 vel, f64 A) {
    if (alt > AMB_softmax_alt)
        return 0.0;
    // Mach-dependant drag.
    f64 M = abs(vel) / amb_a(alt);
    f64 CD = rocket_CD(M);
    f64 Drag = CD * A * 0.5 * amb_rho(alt) * vel*vel;
    // "did you remember to make the drag negative?"
    if (vel > 0.0)
        Drag = -Drag;
    return Drag;
}
