#include "ambient.h"

#include "maths.h"


// International Standard Atmosphere layers. Lapse rate is positive for
// temperature falling with altitude.
typedef struct ambLayer {
    f64 alt; // base altitude [m]
    f64 T;   // base temperature [K]
    f64 P;   // base pressure [Pa]
    f64 L;   // lapse rate [K/m]
} ambLayer;
static const ambLayer amb_layers[] = {
    {     0.0, 288.15,  101325.0, +0.0065 },
    { 11019.0, 216.65,  22632.1,   0.0    },
    { 20063.0, 216.65,  5474.89,  -0.001  },
    { 32162.0, 228.65,  868.02,   -0.0028 },
    { 47350.0, 270.65,  110.91,    0.0    },
    { 51412.0, 270.65,  66.94,    +0.0028 },
    { 71802.0, 214.65,  3.96,     +0.002  },
    { 86000.0, 186.946, 0.3734,    0.0    },
};


// Temperature and pressure at `alt`, within whichever layer contains it (the
// last layer extends upwards forever).
static void amb_T_P(f64 alt, f64* rstr T, f64* rstr P) {
    i32 i = 0;
    while (i+1 < numel(amb_layers) && alt >= amb_layers[i+1].alt)
        ++i;
    const ambLayer* l = amb_layers + i;
    f64 Dalt = alt - l->alt;
    if (l->L == 0.0) {
        // Isothermal layer.
        *T = l->T;
        *P = l->P * exp2(LOG2E * -STANDARD_GRAVITY*Dalt/AMB_STANCOND_R/l->T);
    } else {
        // Gradient layer.
        *T = l->T - l->L * Dalt;
        *P = l->P * pow(*T / l->T, STANDARD_GRAVITY/AMB_STANCOND_R/l->L);
    }
}



f64 amb_T(f64 alt) {
    f64 T, P;
    amb_T_P(alt, &T, &P);
    return T;
}

f64 amb_P(f64 alt) {
    f64 T, P;
    amb_T_P(alt, &T, &P);
    return P;
}

f64 amb_rho(f64 alt) {
    f64 T, P;
    amb_T_P(alt, &T, &P);
    return P / AMB_STANCOND_R / T; // ideal gas.
}

f64 amb_a(f64 alt) {
    f64 T, P;
    amb_T_P(alt, &T, &P);
    f64 Y = AMB_STANCOND_cp / AMB_STANCOND_cv;
    return sqrt(Y * AMB_STANCOND_R * T);
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
