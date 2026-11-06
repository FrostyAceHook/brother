#include "geometry.h"

#include "assertion.h"
#include "maths.h"


f64 circ_diam(f64 A) {
    assert(A >= 0.0, "A=%f", A);
    return sqrt(A / PI_4);
}
f64 circ_area(f64 D) {
    assert(D >= 0.0, "D=%f", D);
    return PI_4 * D*D;
}
f64 circ_perim(f64 D) {
    assert(D >= 0.0, "D=%f", D);
    return PI * D;
}
f64 cyl_vol(f64 L, f64 D) {
    assert(L >= 0.0, "L=%f", L);
    assert(D >= 0.0, "D=%f", D);
    return L * PI_4 * D*D;
}
f64 cyl_L(f64 V, f64 D) {
    assert(V >= 0.0, "V=%f", V);
    assert(D >= 0.0, "D=%f", D);
    return V / PI_4 / D / D;
}
f64 tube_vol(f64 L, f64 ID, f64 OD) {
    assert(L >= 0.0, "L=%f", L);
    assert(OD >= ID, "ID=%f, OD=%f", ID, OD);
    assert(ID >= 0.0, "ID=%f", ID);
    return L * PI_4 * (OD*OD - ID*ID);
}
