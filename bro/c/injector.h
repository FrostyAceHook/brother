#pragma once
#include "br.h"



// ========================= //
//        GAS ORIFICE        //
// ========================= //

typedef struct Hole {
    f64 CdA;
    f64 P_u; // saturated.
    f64 P_d;
    f64 Y;
    f64 R;
    f64 T; // saturated.
    f64 Z;
} Hole;
// Isentropic compressible flow through an orifice, choked or unchoked [kg/s].
// No flow for P_d >= P_u.
f64 hole_dm(const Hole* hole);


// Same as above, but split into an upstream init (everything independent of
// P_d) and a cheap flow query.
typedef struct HoleU {
    f64 P_u;
    f64 Y;
    f64 Prinv_choke;
    f64 dm_choke;
    f64 factor0;
} HoleU;
HoleU holeU(const Hole* hole); // `hole->P_d` is ignored.
f64 holeU_dm(const HoleU* up, f64 P_d);


// ========================== //
//     LIQUID N2O ORIFICE     //
// ========================== //

enum {
    HOLE_N2O_SPC,  // single-phase compressible liquid
    HOLE_N2O_HEM,  // homogeneous equilibrium (two-phase), choke-limited
    HOLE_N2O_DYER, // Dyer/NHNE blend of SPC and HEM, weighted by kappa
    HOLE_N2O_FML,  // FML blend of SPC and HEM, weighted by exit void fraction
};
typedef struct HoleN2O {
    i32 model; // =`HOLE_N2O_*`
    f64 CdA;
    f64 P_u; // stagnation.
    f64 P_d;
    f64 T_u; // stagnation.
    f64 rho_u; // should be saturated or compressed liquid, altho clamped.
} HoleN2O;
// Liquid nitrous orifice mass flow rate. No flow for P_d >= P_u.
f64 hole_N2O_dm(const HoleN2O* hole);


// Same as above, but split into an upstream init (everything independent of
// P_d) and a cheap flow query.
typedef struct HoleN2OU {
    i32 model;          // =`HOLE_N2O_*`
    f64 CdA;            // Cd*A [m^2]
    f64 P_u;            // Upstream pressure [Pa]
    f64 Psat_u;         // Saturation pressure at T_u [Pa] (Dyer)
    f64 s_u;            // Upstream entropy [J/kg/K]
    f64 h_u;            // Upstream enthalpy [J/kg]
    f64 IPLE_u;         // Upstream isentropic power law exponent [-]
    f64 Prinv_chokeSPC; // SPC critical pressure ratio [-]
    f64 factor0_SPC;    // CdA*sqrt(2*rho_u*P_u*n/(n - 1)) [kg/s]
    f64 T_d_chokeHEM;   // HEM choke exit temperature [K] (NaN for SPC)
    f64 P_d_chokeHEM;   // HEM choke exit pressure [Pa] (NaN for SPC)
    f64 X_chokeHEM;     // HEM exit quality when choked [-] (NaN for SPC)
    f64 dm_chokeHEM;    // HEM mass flow rate when choked [kg/s] (NaN for SPC)
    f64 k_chokeFML;     // FML blend weight when choked [-] (NaN unless FML)
} HoleN2OU;
HoleN2OU hole_N2OU(const HoleN2O* hole); // `hole->P_d` is ignored.
f64 hole_N2OU_dm(const HoleN2OU* up, f64 P_d);
