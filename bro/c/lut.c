#include "lut.h"

#include "maths.h"


f64 lut_g(const lutAxis* ax, f64 x) {
    switch (ax->scale) {
      case LUT_LOG:
        return log2(x);
      case LUT_SQRT: {
        f64 d = x - ax->c;
        return (d < 0.0) ? -sqrt(-d) : sqrt(d);
      }
      case LUT_CBRT:
        return cbrt(x - ax->c);
      default:
        return x;
    }
}

f64 lut_frac(const lutAxis* ax, f64 x) {
    return (lut_g(ax, x) - ax->glo) / (ax->ghi - ax->glo) * (ax->n - 1);
}

i32 lut_inside(const lutAxis* ax, f64 x, f64 rtol) {
    f64 tol = rtol * (ax->hi - ax->lo);
    return (ax->lo - tol <= x) && (x <= ax->hi + tol);
}

// Cell index and in-cell weight for x (end cells extrapolate).
static i32 lut_cell_(const lutAxis* ax, f64 x, f64* w) {
    f64 f = lut_frac(ax, x);
    i32 i = (f <= 0.0) ? 0 : (i32)f;
    i = min(i, ax->n - 2);
    *w = f - i;
    return i;
}

f64 lut1(const f32* tbl, const lutAxis* ax, f64 x) {
    f64 t;
    i32 i = lut_cell_(ax, x, &t);
    return lerp((f64)tbl[i], (f64)tbl[i + 1], t);
}

f64 lut2(const f32* tbl, const lutAxis* ax0, const lutAxis* ax1,
        f64 x0, f64 x1) {
    f64 t, u;
    i32 i = lut_cell_(ax0, x0, &t);
    i32 j = lut_cell_(ax1, x1, &u);
    i32 n1 = ax1->n;
    const f32* r0 = tbl + i*n1 + j;
    const f32* r1 = r0 + n1;
    f64 X0 = lerp((f64)r0[0], (f64)r0[1], u);
    f64 X1 = lerp((f64)r1[0], (f64)r1[1], u);
    return lerp(X0, X1, t);
}

f64 lut3(const f32* tbl, const lutAxis* ax0, const lutAxis* ax1,
        const lutAxis* ax2, f64 x0, f64 x1, f64 x2) {
    f64 t, u, v;
    i32 i = lut_cell_(ax0, x0, &t);
    i32 j = lut_cell_(ax1, x1, &u);
    i32 k = lut_cell_(ax2, x2, &v);
    i32 n1 = ax1->n;
    i32 n2 = ax2->n;
    const f32* c00 = tbl + (i*n1 + j)*n2 + k;
    const f32* c01 = c00 + n2;
    const f32* c10 = c00 + n1*n2;
    const f32* c11 = c10 + n2;
    f64 X00 = lerp((f64)c00[0], (f64)c00[1], v);
    f64 X01 = lerp((f64)c01[0], (f64)c01[1], v);
    f64 X10 = lerp((f64)c10[0], (f64)c10[1], v);
    f64 X11 = lerp((f64)c11[0], (f64)c11[1], v);
    f64 X0 = lerp(X00, X01, u);
    f64 X1 = lerp(X10, X11, u);
    return lerp(X0, X1, t);
}
