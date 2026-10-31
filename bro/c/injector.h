#pragma once
#include "br.h"



// ========================= //
//        GAS ORIFICE        //
// ========================= //

typedef struct Hole {
    f64 Cd;  // Discharge coefficient [-]
    f64 A;   // Flow area [m^2]
    f64 P_u; // Upstream (stagnation) pressure [Pa]
    f64 P_d; // Downstream pressure [Pa]
    f64 Y;   // Ratio of specific heats [-]
    f64 R;   // Specific gas constant [J/kg/K]
    f64 T;   // Upstream (stagnation) temperature [K]
    f64 Z;   // Upstream compressibility factor [-]
} Hole;
// Isentropic compressible flow through an orifice, choked or unchoked [kg/s].
// No flow for P_d >= P_u.
f64 hole_dm(const Hole* hole);


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
    f64 Cd;    // Discharge coefficient [-]
    f64 A;     // Flow area [m^2]
    f64 P_u;   // Upstream pressure [Pa]
    f64 P_d;   // Downstream pressure [Pa]
    f64 T_u;   // Upstream temperature [K]
    f64 rho_u; // Upstream density [kg/m^3]; saturated or compressed liquid
               //   (values slightly below saturated are treated as saturated)
} HoleN2O;
// Liquid nitrous orifice mass flow rate. No flow for P_d >= P_u.
f64 hole_N2O_dm(const HoleN2O* hole);
