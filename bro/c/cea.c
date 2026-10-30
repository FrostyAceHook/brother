#include "cea.h"

#include "assertion.h"
#include "lut.h"
#include "maths.h"


typedef struct ceaGrid_ {
    i32 dims; // 2: (P, ofr), 3: (P, ofr, AR)
    lutAxis ax[3];
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

static f64 cea_lookup_(const ceaTable_* tables, i32 prop, const char* what,
        f64 P0_cc, f64 ofr, f64 AR) {
    assert(0 <= prop && prop < CEA_PROP_COUNT, "invalid CEA prop: %d", prop);
    const ceaTable_* t = &tables[prop];
    assert(t->tbl != NULL, "CEA %s not tabulated for %s (AR=%g)", what,
            cea_prop_name_(prop), AR);
    const ceaGrid_* g = t->grid;

    f64 x[3] = { P0_cc, ofr, AR };
    const char* names[3] = { "P0_cc", "ofr", "AR" };
    for (i32 d=0; d<g->dims; ++d) {
        assert(lut_inside(&g->ax[d], x[d], 1e-9),
                "CEA input oob for %s: %s=%g not in [%g, %g]",
                cea_prop_name_(prop), names[d], x[d], g->ax[d].lo,
                g->ax[d].hi);
    }

    f64 X = (g->dims == 2)
          ? lut2(t->tbl, &g->ax[0], &g->ax[1], x[0], x[1])
          : lut3(t->tbl, &g->ax[0], &g->ax[1], &g->ax[2], x[0], x[1], x[2]);
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
