#include "injector.h"

#include "assertion.h"
#include "maths.h"
#include "n2o.h"



// =========================================================================== //
// = GAS ORIFICE ============================================================= //
// =========================================================================== //

f64 hole_dm(const Hole* hole) {
    // Let:
    //  Prinv = P_d / P_u
    //  beta = (2 / (Y + 1)) ** (1 / (Y - 1))
    //  Prinv_choke = (2 / (Y + 1)) * beta
    //             [= (2/(Y+1))**(Y/(Y-1))]
    // Choked (Prinv <= Prinv_choke):
    //  alpha = Prinv_choke * beta
    //       [= (2/(Y+1))**((Y+1)/(Y-1))]
    // Unchoked, with g = Prinv**(1/Y):
    //  alpha = 2/(Y - 1) * g * (g - Prinv)
    //       [= 2/(Y-1)*(Pr**(2/Y) - Pr**((Y+1)/Y))]
    // Mass flow rate:
    //  dm = P_u * Cd * A * sqrt(alpha * Y / (Z * R * T))
    f64 P_u = hole->P_u;
    f64 Y = hole->Y;
    if (P_u <= 0.0)
        return 0.0;
    f64 Prinv = hole->P_d / P_u;
    if (Prinv >= 1.0)
        return 0.0;
    f64 beta = pow(2.0 / (Y + 1.0), 1.0 / (Y - 1.0));
    f64 Prinv_choke = 2.0 / (Y + 1.0) * beta;
    f64 alpha;
    if (Prinv <= Prinv_choke) {
        alpha = Prinv_choke * beta;
    } else {
        f64 g = pow(Prinv, 1.0 / Y);
        alpha = 2.0 / (Y - 1.0) * g * (g - Prinv);
    }
    return P_u * hole->Cd * hole->A
         * sqrt(alpha * Y / (hole->Z * hole->R * hole->T));
}



// =========================================================================== //
// = LIQUID NITROUS ORIFICE ================================================== //
// =========================================================================== //

// Homogeneous-equilibrium exit state for an isentropic expansion from (s_u, h_u)
// to saturation temperature T_d (pressure Psat(T_d)). Parametrising by T_d
// rather than P_d avoids inverting Psat on every evaluation.
typedef struct HemExit_ {
    f64 X;   // exit quality [-]
    f64 G;   // mass flux [kg/m^2/s] (0 if not in the dome or no enthalpy drop)
    f64 rho_l, rho_v; // exit saturated densities [kg/m^3]
} HemExit_;
static HemExit_ hem_exit_(f64 s_u, f64 h_u, f64 T_d) {
    HemExit_ e;
    f64 s_l = N2O_satliq_s(T_d);
    f64 s_v = N2O_satvap_s(T_d);
    e.X = (s_u - s_l) / (s_v - s_l);
    e.rho_l = N2O_satliq_rho(T_d);
    e.rho_v = N2O_satvap_rho(T_d);
    if (e.X <= 0.0) {
        e.G = 0.0;
        return e;
    }
    f64 h_d = lerp(N2O_satliq_h(T_d), N2O_satvap_h(T_d), e.X);
    f64 v_d = lerp(1.0/e.rho_l, 1.0/e.rho_v, e.X); // mix in specific volume.
    f64 dh = h_u - h_d;
    // A meaningfully negative drop means the property tables are inconsistent
    // (cf. the old compressed-liquid/saturated mismatch), not physics.
    assert(dh > -100.0, "HEM enthalpy rise: dh=%g, T_d=%g", dh, T_d);
    e.G = sqrt(2.0*max(dh, 0.0)) / v_d;
    return e;
}

// HEM choke: the exit saturation temperature maximising the mass flux (where the
// flow speed equals the equilibrium sound speed). G is unimodal in T_d, so
// golden-section search; for compressed inlets G is zero above the saturation
// crossing and the search converges onto the crossing, which is HEM's choke
// there.
static f64 hem_choke_T_(f64 s_u, f64 h_u, f64 T_u) {
    f64 a = N2O_Ttrip + 0.01;
    f64 b = T_u - 1e-6;
    const f64 g = 0.6180339887498949;
    f64 c = b - g*(b - a);
    f64 d = a + g*(b - a);
    f64 Gc = hem_exit_(s_u, h_u, c).G;
    f64 Gd = hem_exit_(s_u, h_u, d).G;
    while (b - a > 1e-5) { // [K]; ~30 iterations from a ~125 K bracket.
        if (Gc > Gd) {
            b = d; d = c; Gd = Gc;
            c = b - g*(b - a);
            Gc = hem_exit_(s_u, h_u, c).G;
        } else {
            a = c; c = d; Gc = Gd;
            d = a + g*(b - a);
            Gd = hem_exit_(s_u, h_u, d).G;
        }
    }
    return 0.5*(a + b);
}

f64 hole_N2O_dm(const HoleN2O* hole) {
    // FML nitrous oxide injector mass flow rate model.
    // https://doi.org/10.3390/aerospace9120828
    // Describes SPI/SPC/HEM/Dyer well also.
    // https://doi.org/10.13009/EUCASS2023-728
    // follow-up goes into further/different detail also.
    i32 model = hole->model;
    f64 Cd = hole->Cd;
    f64 A = hole->A;
    f64 P_u = hole->P_u;
    f64 P_d = hole->P_d;
    f64 T_u = hole->T_u;

    assert(model == HOLE_N2O_SPC
        || model == HOLE_N2O_HEM
        || model == HOLE_N2O_DYER
        || model == HOLE_N2O_FML, "invalid injector model: %d", model);

    // No flow without a forward pressure difference.
    if (P_u <= 0.0)
        return 0.0;
    f64 Prinv = P_d / P_u;
    if (Prinv >= 0.99)
        return 0.0;

    // Upstream state: saturated liquid, or compressed liquid at (T_u, rho_u).
    // Densities a hair below saturated (round-off from the caller) count as
    // saturated; the general lookup would otherwise put them in the dome.
    //
    // s_u and h_u are anchored to the saturated tables, plus the compressed-
    // liquid offset from the dense table. HEM's enthalpy drop (h_u - h_d, with
    // h_d from the saturated tables) is a small difference of large numbers
    // (~250 J/kg at P_d/P_u = 0.95), so h_u and h_d must come from the same,
    // most accurate tables: the dense table's absolute h error (~10 J/kg) would
    // otherwise be a few % of the drop near P_u. The offset form cancels the
    // dense table's error at psi = 0 exactly.
    f64 rho_sat = N2O_satliq_rho(T_u);
    f64 rho_u = max(hole->rho_u, rho_sat);
    f64 s_u = N2O_satliq_s(T_u);
    f64 h_u = N2O_satliq_h(T_u);
    f64 n_u = N2O_satliq_IPLE(T_u);
    if (rho_u > rho_sat) {
        s_u += N2O_s(T_u, rho_u) - N2O_s(T_u, rho_sat);
        h_u += N2O_h(T_u, rho_u) - N2O_h(T_u, rho_sat);
        n_u = N2O_IPLE(T_u, rho_u);
    }
    assert(n_u > 1.0, "IPLE_u=%g", n_u);


    // Single-phase compressible (SPC), choked at the isentropic power-law
    // critical ratio (in practice never reached by liquid: n ~ 30).
    //  x**(2/n) - x**((n + 1)/n)  =  x**(1/n) * (x**(1/n) - x)
    f64 Prinv_choke_SPC = pow(2.0/(n_u + 1.0), n_u/(n_u - 1.0));
    f64 Prinv_SPC = max(Prinv, Prinv_choke_SPC);
    f64 x_1_n = pow(Prinv_SPC, 1.0/n_u);
    f64 inner = x_1_n * (x_1_n - Prinv_SPC);
    f64 dm_SPC = Cd*A*sqrt(2.0*rho_u*P_u*n_u/(n_u - 1.0)*inner);

    if (model == HOLE_N2O_SPC)
        return dm_SPC;


    // Homogeneous equilibrium (HEM), evaluated at the exit-plane state: the
    // downstream pressure, or the choke if that's higher.
    f64 T_choke = hem_choke_T_(s_u, h_u, T_u);
    f64 P_choke = N2O_Psat(T_choke);
    f64 T_e = (max(P_d, N2O_Ptrip) <= P_choke)
            ? T_choke
            : N2O_Tsat(max(P_d, N2O_Ptrip));
    HemExit_ e = hem_exit_(s_u, h_u, T_e);

    // Still single-phase (liquid) at the exit: SPC is the right answer.
    if (e.X <= 0.0)
        return dm_SPC;
    assert(e.X <= 1.0, "HEM exit quality X_e=%g", e.X);
    f64 dm_HEM = Cd*A*e.G;

    if (model == HOLE_N2O_HEM)
        return dm_HEM;


    // Blend SPC and HEM: dm = (1 - k)*dm_SPC + k*dm_HEM.
    f64 k;
    if (model == HOLE_N2O_DYER) {
        // Dyer/Solomon: kappa = sqrt((P_u - P_d) / (Psat(T_u) - P_d)),
        // k = 1/(1 + kappa). Uses the actual downstream pressure.
        f64 kappa = sqrt((P_u - P_d) / (N2O_Psat(T_u) - P_d));
        k = 1.0 / (1.0 + kappa);
    } else /* HOLE_N2O_FML */ {
        // FML: k = exit void fraction with Zivi slip S = (rho_l/rho_v)^(1/3),
        //  alpha = 1 / (1 + (1 - X)/X * S * rho_v/rho_l)
        //        = 1 / (1 + (1 - X)/X * (rho_v/rho_l)^(2/3)),
        // at the same (choke-limited) exit state as the HEM term.
        f64 r = pow(e.rho_v / e.rho_l, 2.0/3.0);
        k = 1.0 / (1.0 + (1.0 - e.X)/e.X * r);
    }
    assert(within(k, 0.0, 1.0), "injector blend weight k=%g", k);

    return dm_SPC*(1.0 - k) + dm_HEM*k;
}
