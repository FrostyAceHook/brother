#include "geometry.h"

#include "assertion.h"
#include "maths.h"


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


f64 pipe_th(const Pipe* p) {
    assert(notnan(p->ID), "must give");
    assert(notnan(p->OD), "must give");
    assert(p->OD >= p->ID, "ID=%f, OD=%f", p->ID, p->OD);
    return 0.5 * (p->OD - p->ID);
}
f64 pipe_V(const Pipe* p) {
    assert(notnan(p->ID), "must give");
    assert(notnan(p->OD), "must give");
    assert(p->L >= 0.0, "L=%f", p->L);
    assert(p->OD >= p->ID, "ID=%f, OD=%f", p->ID, p->OD);
    return tube_vol(p->L, p->ID, p->OD);
}
void pipe_set_th(Pipe* p, f64 th) {
    assert(notnan(p->ID) + notnan(p->OD) != 2, "must give exactly one");
    assert(notnan(p->ID) + notnan(p->OD) != 0, "must give exactly one");
    if (notnan(p->ID)) {
        assert(th >= 0.0, "ID=%f, th=%f", p->ID, th);
        p->OD = p->ID + 2.0*th;
    } else {
        assert(2.0*th <= p->OD, "OD=%f, th=%f", p->OD, th);
        p->ID = p->OD - 2.0*th;
    }
}
void pipe_set_th_for_P(Pipe* p, f64 P, f64 Ys, f64 sf) {
    // Finding min thickness for hoop stress, using thin-walled approximation:
    //  th = sf * P * ((ID + OD) / 2) / (2 * Ys)
    // Consider:
    // OD = ID + 2*th  [constrained inner diameter]
    //  th = sf * P * (2*ID + 2*th) / 4 / Ys
    //  th = sf * P * (ID + th) / 2 / Ys
    //  th = sf * P/2/Ys * ID + sf * P/2/Ys * th
    //  th - sf * P/2/Ys * th = sf * P/2/Ys * ID
    //  th * (1 - sf * P/2/Ys) = sf * P/2/Ys * ID
    //  th * (2 * Ys - sf * P) = sf * P * ID
    // => th = sf * P * ID / (2 * Ys - sf * P)
    // Consider:
    // ID = OD - 2*th  [constrained outer diameter]
    //  th = sf * P * (2*OD - 2*th) / 4 / Ys
    //  th = sf * P * (OD - th) / 2 / Ys
    //  th = sf * P/2/Ys * OD - sf * P/2/Ys * th
    //  th + sf * P/2/Ys * th = sf * P/2/Ys * OD
    //  th * (1 + sf * P/2/Ys) = sf * P/2/Ys * OD
    //  th * (2 * Ys + sf * P) = sf * P * OD
    // => th = sf * P * OD / (2 * Ys + sf * P)
    // Note that these relations mean it is not always possible to construct a
    // thick enough tank (in either sitation).
    assert(notnan(p->ID) + notnan(p->OD) != 2, "must give exactly one");
    assert(notnan(p->ID) + notnan(p->OD) != 0, "must give exactly one");
    if (notnan(p->ID)) {
        f64 th = sf * P * p->ID / (2.0 * Ys - sf * P);
        assert(th > 0.0, "ID=%f, th=%f, P=%f, Ys=%f, sf=%f", p->ID, th, P, Ys,
                sf);
        p->OD = p->ID + 2.0*th;
    } else {
        f64 th = sf * P * p->OD / (2.0 * Ys + sf * P);
        assert(2.0*th < p->OD, "OD=%f, th=%f, P=%f, Ys=%f, sf=%f", p->OD, th, P,
                Ys, sf);
        p->ID = p->OD - 2.0*th;
    }
}
f64 pipe_sf_for_P(const Pipe* p, f64 P, f64 Ys) {
    // Finding sf from hoop stress, using thin-walled approximation:
    //  th = sf * P * ((ID + OD) / 2) / (2 * Ys)
    // => sf = 4 * th * Ys / (P * (ID + OD))
    return 4.0 * pipe_th(p) * Ys / P / (p->ID + p->OD);
}
