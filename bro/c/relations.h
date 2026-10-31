#pragma once
#include "br.h"


// Straight off the dome.
#define GAS_CONSTANT (8.31446261815324)
#define STEFAN_BOLTZMAN_CONSTANT (5.670374419184429e-8)


typedef struct SpecificHeatRatio {
    f64 Y;
    f64 n; /* (Y + 1)/(Y - 1)/2 */
    f64 sup_M_seed_n; /* (1 + Y - Y^0.2)/3 */
    f64 sup_M_seed_m; /* Y^3 / 64 */
    f64 sub_M_seed_n; /* 0.55 + 0.7/(Y + 2) */
    f64 sub_M_seed_m; /* 0.8*Y + 2 */
} SpecificHeatRatio;
#define get_shr(args...) ( init_shr(&(SpecificHeatRatio){0}, args) )
SpecificHeatRatio* init_shr(SpecificHeatRatio* shr, f64 gamma);


f64 isentropic_sup_M(f64 A_on_Astar, const SpecificHeatRatio* shr);
f64 isentropic_sub_M(f64 A_on_Astar, const SpecificHeatRatio* shr);
f64 isentropic_M(i32 subsonic, f64 A_on_Astar, const SpecificHeatRatio* shr);

void isentropic_shr_M(SpecificHeatRatio* shr, f64* rstr M, i32 subsonic,
        f64 A_on_Astar, f64 P0_cc, f64 ofr, f64 gamma_cc, f64 gamma_tht);

f64 isentropic_A_on_Astar(f64 M, const SpecificHeatRatio* shr);

f64 isentropic_T_on_T0(f64 M, const SpecificHeatRatio* shr);
f64 isentropic_P_on_P0(f64 M, const SpecificHeatRatio* shr);
f64 isentropic_rho_on_rho0(f64 M, const SpecificHeatRatio* shr);

f64 get_Y1M22(f64 M, const SpecificHeatRatio* shr);
f64 isentropicx_T_on_T0(f64 Y1M22, const SpecificHeatRatio* shr);
f64 isentropicx_P_on_P0(f64 Y1M22, const SpecificHeatRatio* shr);
f64 isentropicx_rho_on_rho0(f64 Y1M22, const SpecificHeatRatio* shr);

f64 isentropic_M_from_P_on_P0(f64 P_on_P0, const SpecificHeatRatio* shr);
