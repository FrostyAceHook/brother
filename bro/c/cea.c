#include "cea.h"

#include "assertion.h"
#include "maths.h"


// Table axis: evenly spaced in x, or in log(x) if `log`.
typedef struct ceaAxis_ {
    f64 lo;
    f64 hi;
    i32 n;
    i32 log;
} ceaAxis_;
typedef struct ceaGrid_ {
    i32 dims; // 2: (P, ofr), 3: (P, ofr, AR)
    ceaAxis_ ax[3];
} ceaGrid_;
typedef struct ceaTable_ {
    const f32* tbl;       // C-ordered, NULL if not tabulated.
    const ceaGrid_* grid;
} ceaTable_;

// Generated: per-propellant tables, and cea_<station>_<field>_tbls_[] arrays
// indexed by ceaProp.
#include "tbl/cea_index.i"


static const char* cea_prop_name_(i32 prop) {
    switch (prop) {
        case CEA_N2O_PARAFFIN: return "N2O/paraffin";
        default: return "?";
    }
}

// Fractional index of `x` along `ax`, asserting it's within the axis.
static f64 cea_frac_(const ceaAxis_* ax, f64 x, const char* what) {
    f64 tol = 1e-9 * max(abs(ax->lo), abs(ax->hi));
    assert(ax->lo - tol <= x && x <= ax->hi + tol,
            "CEA input oob: %s=%g not in [%g, %g]", what, x, ax->lo, ax->hi);
    f64 t = (ax->log) ? log2(x / ax->lo) / log2(ax->hi / ax->lo)
                      : (x - ax->lo) / (ax->hi - ax->lo);
    return t * (ax->n - 1);
}

static f64 cea_lookup_(const ceaTable_* tables, i32 prop, const char* what,
        f64 P0_cc, f64 ofr, f64 AR) {
    assert(0 <= prop && prop < CEA_PROP_COUNT, "invalid CEA prop: %d", prop);
    const ceaTable_* t = &tables[prop];
    assert(t->tbl != NULL, "CEA %s not tabulated for %s (AR=%g)", what,
            cea_prop_name_(prop), AR);
    const ceaGrid_* g = t->grid;

    f64 x[3] = { P0_cc, ofr, AR };
    const char* names[3] = { "P0_cc", "ofr", "AR" };
    i32 idx[3];
    f64 w[3];
    for (i32 d=0; d<g->dims; ++d) {
        f64 f = cea_frac_(&g->ax[d], x[d], names[d]);
        idx[d] = min(max((i32)f, 0), g->ax[d].n - 2);
        w[d] = f - idx[d];
    }

    f64 X;
    if (g->dims == 2) {
        i32 ny = g->ax[1].n;
        const f32* r0 = t->tbl + idx[0]*ny + idx[1];
        const f32* r1 = r0 + ny;
        f64 X0 = lerp((f64)r0[0], (f64)r0[1], w[1]);
        f64 X1 = lerp((f64)r1[0], (f64)r1[1], w[1]);
        X = lerp(X0, X1, w[0]);
    } else {
        i32 ny = g->ax[1].n;
        i32 nz = g->ax[2].n;
        const f32* c00 = t->tbl + (idx[0]*ny + idx[1])*nz + idx[2];
        const f32* c01 = c00 + nz;
        const f32* c10 = c00 + ny*nz;
        const f32* c11 = c10 + nz;
        f64 X00 = lerp((f64)c00[0], (f64)c00[1], w[2]);
        f64 X01 = lerp((f64)c01[0], (f64)c01[1], w[2]);
        f64 X10 = lerp((f64)c10[0], (f64)c10[1], w[2]);
        f64 X11 = lerp((f64)c11[0], (f64)c11[1], w[2]);
        f64 X0 = lerp(X00, X01, w[1]);
        f64 X1 = lerp(X10, X11, w[1]);
        X = lerp(X0, X1, w[0]);
    }
    assert(notnan(X), "CEA %s nan for %s: P0_cc=%g, ofr=%g, AR=%g", what,
            cea_prop_name_(prop), P0_cc, ofr, AR);
    return X;
}



// Station dispatch on AR.
#define CEA_DEFINE_STATION_(name)                                           \
    f64 cea_##name(i32 prop, f64 P0_cc, f64 ofr, f64 AR) {                  \
        const ceaTable_* tables;                                            \
        if (AR == CEA_AR_stag)                                              \
            tables = cea_cc_##name##_tbls_;                                 \
        else if (AR == CEA_AR_tht)                                          \
            tables = cea_tht_##name##_tbls_;                                \
        else if (AR > 1.0)                                                  \
            tables = cea_sup_##name##_tbls_;                                \
        else {                                                              \
            assert(AR > 0.0, "CEA input oob: AR=%g", AR);                   \
            tables = cea_sub_##name##_tbls_;                                \
        }                                                                   \
        return cea_lookup_(tables, prop, #name, P0_cc, ofr, AR);            \
    }
CEA_DEFINE_STATION_(T)
CEA_DEFINE_STATION_(P)
CEA_DEFINE_STATION_(rho)
CEA_DEFINE_STATION_(a)
CEA_DEFINE_STATION_(gamma)
CEA_DEFINE_STATION_(cp)
CEA_DEFINE_STATION_(mu)
CEA_DEFINE_STATION_(Pr)
CEA_DEFINE_STATION_(Mw)
CEA_DEFINE_STATION_(Ivac)
#undef CEA_DEFINE_STATION_

f64 cea_M(i32 prop, f64 P0_cc, f64 ofr, f64 AR) {
    assert(0 <= prop && prop < CEA_PROP_COUNT, "invalid CEA prop: %d", prop);
    if (AR == CEA_AR_stag)
        return 0.0;
    if (AR == CEA_AR_tht)
        return 1.0;
    const ceaTable_* tables = (AR > 1.0) ? cea_sup_M_tbls_ : cea_sub_M_tbls_;
    assert(AR > 0.0, "CEA input oob: AR=%g", AR);
    return cea_lookup_(tables, prop, "M", P0_cc, ofr, AR);
}
