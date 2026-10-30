#pragma once
#include "br.h"

// Evenly-spaced lookup tables over a transformed coordinate.
//
// Each table axis is evenly spaced in g(x) rather than in x, so nodes can be
// concentrated where a function varies fastest:
//   LUT_LIN   g(x) = x
//   LUT_LOG   g(x) = log2(x)                    (x > 0; spans decades evenly)
//   LUT_SQRT  g(x) = sign(x - c) * |x - c|^(1/2) (nodes dense near x = c)
//   LUT_CBRT  g(x) = cbrt(x - c)                 (nodes denser still near c)
//
// Tables are flattened, C-ordered f32 arrays. Lookups interpolate (multi-)
// linearly in g; inputs outside an axis extrapolate from the end cell,
// so callers check bounds (lut_inside) where that matters.

enum {
    LUT_LIN,
    LUT_LOG,
    LUT_SQRT,
    LUT_CBRT,
};
typedef i32 lutScale;

typedef struct lutAxis {
    f64 lo;    // physical lower bound (first node)
    f64 hi;    // physical upper bound (last node)
    f64 glo;   // g(lo), precomputed by the table generator
    f64 ghi;   // g(hi), precomputed by the table generator
    i32 n;     // number of nodes, >=2
    lutScale scale; // how the function g(x) is shaped.
    f64 c;     // centre for LUT_SQRT/LUT_CBRT (unused otherwise)
} lutAxis;

// Transformed coordinate g(x) of an axis.
f64 lut_g(const lutAxis* ax, f64 x);

// Fractional node index of x along an axis (unclamped).
f64 lut_frac(const lutAxis* ax, f64 x);

// Whether x lies within [lo, hi], with a relative tolerance on the span.
i32 lut_inside(const lutAxis* ax, f64 x, f64 rtol);

// Interpolated lookups. `tbl` holds prod(n) values, last axis fastest.
f64 lut1(const f32* tbl, const lutAxis* ax, f64 x);
f64 lut2(const f32* tbl, const lutAxis* ax0, const lutAxis* ax1,
        f64 x0, f64 x1);
f64 lut3(const f32* tbl, const lutAxis* ax0, const lutAxis* ax1,
        const lutAxis* ax2, f64 x0, f64 x1, f64 x2);
