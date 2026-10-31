#include "ambient.h"

#include "maths.h"


f64 amb_P(f64 alt) {
    // Just say properties below sea level are constant.
    alt = max(alt, 0.0);
    // Above 60k, its about 0 pressure (21.5 Pa at 60k, 0.021% of sea level).
    if (alt > 60e3)
        return 0.0;
    // Different approxs for different regions.
    if (alt > 20.063e3) {
        f64 x1 = alt;
        f64 x2 = x1*x1;
        f64 n0 = +2839481483.3687034;
        f64 n1 = -103143.5518901285;
        f64 d0 = +239165.4892657625;
        f64 d1 = -24.24291511989479;
        f64 d2 = +0.00114675186987331;
        f64 Num = n0 + n1*x1 + x2;
        f64 Den = d0 + d1*x1 + d2*x2;
        return Num / Den;
    }
    f64 x1 = alt + 10e3; // need some offset to avoid /0.
    f64 x2 = x1*x1;
    f64 n0 = +1326889135.6347075;
    f64 n1 = -70527.4205624988;
    f64 d0 = +5522.836249663364;
    f64 d2 = +1.6121842727811546e-05;
    f64 Num = n0 + n1*x1 + x2;
    f64 Den = d0 + d2*x2;
    return Num / Den;
}

f64 amb_rho(f64 alt) {
    // Just say properties below sea level are constant.
    alt = max(alt, 0.0);
    // Above 60k, its about 0 density (0.303 g/m^3 at 60k, 0.025% of sea level).
    if (alt > 60e3)
        return 0.0;
    // Different approxs for different regions.
    if (alt > 20.063e3) {
        f64 x1 = alt;
        f64 x2 = x1*x1;
        f64 n0 = +2239506917.973228;
        f64 n1 = -92441.9508707563;
        f64 d0 = -65681468.84943648;
        f64 d2 = +22.324464921748664;
        f64 Num = n0 + n1*x1 + x2;
        f64 Den = d0 + d2*x2;
        return Num / Den;
    }
    f64 x1 = alt + 10e3; // offset to avoid /0.
    f64 x2 = x1*x1;
    f64 n0 = +1080314808.6263719;
    f64 n1 = -64847.90278878072;
    f64 d0 = +345190107.69729966;
    f64 d1 = +12653.90131058989;
    f64 d2 = -0.37409509100171673;
    f64 Num = n0 + n1*x1 + x2;
    f64 Den = d0 + d1*x1 + d2*x2;
    return Num / Den;
}

f64 amb_a(f64 alt) {
    // Just say properties below sea level are constant.
    alt = max(alt, 0.0);
    // Above 60k, the speed of sound keeps changing but we dont keep calcing it,
    // since we'd only ever use it together with density and that goes to zero so
    // it becomes irrelevant.
    if (alt > 60e3)
        alt = 60e3;
    f64 x1 = alt + 10e3; // offset to avoid /0.
    f64 x2 = x1*x1;
    f64 x3 = x2*x1;
    f64 x4 = x3*x1;
    f64 n1 = +2453911045.412046;
    f64 d0 = +27961612000.005264;
    f64 d1 = -46762.970478368705;
    f64 d2 = +565.8014318466305;
    f64 d3 = -0.00988835027584721;
    f64 d4 = +9.180518407728973e-08;
    f64 Num = n1*x1 + x3;
    f64 Den = d0 + d1*x1 + d2*x2 + d3*x3 + d4*x4;
    return Num / Den;
}

// oh yeah we modelling gravity dropoff in this bitch.
f64 amb_g(f64 alt) {
    // Newtons law of universal gravitation:
    //  F = G*M*m / r^2
    //  g = G*M / r^2
    // let:
    // R = 6.3568e6  [radius of earth]
    // g0 = G*M / R^2  [standard gravity]
    //  g = G*M / r^2 * R^2 / R^2
    //  g = R^2 / r^2 * G*M / R^2
    //  g = R^2 / r^2 * g0
    //  g = R^2 / (R + alt)^2 * g0
    f64 R = 6.3568e6;
    f64 factor = R / (R + alt);
    return STANDARD_GRAVITY * factor*factor;
}
