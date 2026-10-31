#include "injector.h"

#include "assertion.h"
#include "maths.h"
#include "n2o.h"



// =========================================================================== //
// = GAS ORIFICE ============================================================= //
// =========================================================================== //

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
//  dm = P_u * CdA * sqrt(alpha * Y / (Z * R * T))

HoleU holeU(const Hole* hole) {
    HoleU up;
    f64 Y = hole->Y;
    up.P_u = hole->P_u;
    up.Y = Y;
    up.factor0 = hole->P_u * hole->CdA
                * sqrt(Y / hole->Z / hole->R / hole->T);
    f64 beta = pow(2.0 / (Y + 1.0), 1.0 / (Y - 1.0));
    up.Prinv_choke = 2.0 / (Y + 1.0) * beta;
    up.dm_choke = up.factor0 * sqrt(up.Prinv_choke * beta);
    return up;
}

f64 holeU_dm(const HoleU* up, f64 P_d) {
    if (up->P_u <= 0.0)
        return 0.0;
    f64 Prinv = P_d / up->P_u;
    if (Prinv >= 1.0)
        return 0.0;
    if (Prinv <= up->Prinv_choke)
        return up->dm_choke;
    f64 Y = up->Y;
    f64 g = pow(Prinv, 1.0 / Y);
    return up->factor0 * sqrt(2.0 / (Y - 1.0) * g * (g - Prinv));
}

f64 hole_dm(const Hole* hole) {
    HoleU up = holeU(hole);
    return holeU_dm(&up, hole->P_d);
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

// FML blend weight: exit void fraction with Zivi slip S = (rho_l/rho_v)^(1/3),
//  alpha = 1 / (1 + (1 - X)/X * S * rho_v/rho_l)
//        = 1 / (1 + (1 - X)/X * (rho_v/rho_l)^(2/3)).
static f64 get_k_FML(HemExit_ e) {
    f64 r = pow(e.rho_v / e.rho_l, 2.0/3.0);
    return 1.0 / (1.0 + (1.0 - e.X)/e.X * r);
}

// HEM choke: the exit saturation temperature maximising the mass flux (where the
// flow speed equals the equilibrium sound speed). G is unimodal in T_d, so
// golden-section search; for compressed inlets G is zero above the saturation
// crossing and the search converges onto the crossing, which is HEM's choke
// there.
static f64 get_T_chokeHEM(f64 s_u, f64 h_u, f64 T_u) {
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

HoleN2OU hole_N2OU(const HoleN2O* hole) {
    HoleN2OU up;

    // FML nitrous oxide injector mass flow rate model.
    // https://doi.org/10.3390/aerospace9120828
    // Describes SPI/SPC/HEM/Dyer well also.
    // https://doi.org/10.13009/EUCASS2023-728
    // follow-up goes into further/different detail also.

    f64 T_u = hole->T_u;
    up.model = hole->model;
    up.CdA = hole->CdA;
    up.P_u = hole->P_u;
    up.Psat_u = N2O_Psat(T_u);

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
    f64 rho_satliq = N2O_satliq_rho(T_u);
    f64 rho_u = max(hole->rho_u, rho_satliq);
    up.s_u = N2O_satliq_s(T_u);
    up.h_u = N2O_satliq_h(T_u);
    up.IPLE_u = N2O_satliq_IPLE(T_u);
    if (rho_u > rho_satliq) {
        up.s_u += N2O_s(T_u, rho_u) - N2O_s(T_u, rho_satliq);
        up.h_u += N2O_h(T_u, rho_u) - N2O_h(T_u, rho_satliq);
        up.IPLE_u = N2O_IPLE(T_u, rho_u);
    }
    f64 n = up.IPLE_u;
    assert(n > 1.0, "IPLE_u=%g", n);

    // SPC: critical ratio and the P_d-independent prefactor.
    up.Prinv_chokeSPC = pow(2.0/(n + 1.0), n/(n - 1.0));
    up.factor0_SPC = up.CdA * sqrt(2.0*rho_u*up.P_u*n/(n - 1.0));

    // HEM: the choke, and the whole HEM term (and FML weight) at the choke,
    // which are what every choked query returns.
    if (up.model == HOLE_N2O_SPC) {
        up.T_d_chokeHEM = NAN;
        up.P_d_chokeHEM = NAN;
        up.X_chokeHEM = NAN;
        up.dm_chokeHEM = NAN;
        up.k_chokeFML = NAN;
    } else {
        // Use the lut for saturated inlets, otherwise search.
        f64 T_d_sat_chokeHEM = (rho_u == rho_satliq)
                             ? N2O_satliq_T_chokeHEM(T_u)
                             : NAN;
        up.T_d_chokeHEM = (notnan(T_d_sat_chokeHEM))
                        ? T_d_sat_chokeHEM
                        : get_T_chokeHEM(up.s_u, up.h_u, T_u);
        up.P_d_chokeHEM = N2O_Psat(up.T_d_chokeHEM);
        HemExit_ e = hem_exit_(up.s_u, up.h_u, up.T_d_chokeHEM);
        up.X_chokeHEM = e.X;
        up.dm_chokeHEM = up.CdA * e.G;
        up.k_chokeFML = (up.model == HOLE_N2O_FML && e.X > 0.0)
                       ? get_k_FML(e)
                       : NAN;
    }
    return up;
}

f64 hole_N2OU_dm(const HoleN2OU* up, f64 P_d) {
    // No flow without a forward pressure difference.
    if (up->P_u <= 0.0)
        return 0.0;
    f64 Prinv = P_d / up->P_u;
    if (Prinv >= 0.99)
        return 0.0;


    // Single-phase compressible (SPC), choked at the isentropic power-law
    // critical ratio.
    //  x**(2/IPLE) - x**((IPLE + 1)/IPLE)  =  x**(1/IPLE) * (x**(1/IPLE) - x)
    f64 Prinv_SPC = max(Prinv, up->Prinv_chokeSPC);
    f64 x_1_n = pow(Prinv_SPC, 1.0/up->IPLE_u);
    f64 dm_SPC = up->factor0_SPC * sqrt(x_1_n * (x_1_n - Prinv_SPC));

    if (up->model == HOLE_N2O_SPC)
        return dm_SPC;


    // Homogeneous equilibrium (HEM), evaluated at the exit-plane state: the
    // downstream pressure, or the choke if that's higher. Choked queries use
    // the cached choke values.
    f64 P_e = max(P_d, N2O_Ptrip);
    i32 choked = (P_e <= up->P_d_chokeHEM);
    f64 dm_HEM;
    f64 k_FML = NAN;
    if (choked) {
        // Still single-phase (liquid) at the exit: SPC is the right answer.
        // (Possible only for compressed inlets, choking at the saturation
        // crossing.)
        if (up->X_chokeHEM <= 0.0)
            return dm_SPC;
        dm_HEM = up->dm_chokeHEM;
        k_FML = up->k_chokeFML;
    } else {
        HemExit_ e = hem_exit_(up->s_u, up->h_u, N2O_Tsat(P_e));
        // Still single-phase (liquid) at the exit: SPC is the right answer.
        if (e.X <= 0.0)
            return dm_SPC;
        assert(e.X <= 1.0, "HEM exit quality X_e=%g", e.X);
        dm_HEM = up->CdA * e.G;
        if (up->model == HOLE_N2O_FML)
            k_FML = get_k_FML(e);
    }

    if (up->model == HOLE_N2O_HEM)
        return dm_HEM;


    // Blend SPC and HEM: dm = (1 - k)*dm_SPC + k*dm_HEM.
    f64 k;
    if (up->model == HOLE_N2O_DYER) {
        // Dyer/Solomon: kappa = sqrt((P_u - P_d) / (Psat(T_u) - P_d)),
        // k = 1/(1 + kappa). Uses the actual downstream pressure.
        f64 kappa = sqrt((up->P_u - P_d) / (up->Psat_u - P_d));
        k = 1.0 / (1.0 + kappa);
    } else /* HOLE_N2O_FML */ {
        // At the same (choke-limited) exit state as the HEM term.
        k = k_FML;
    }
    assert(within(k, 0.0, 1.0), "k=%g", k);

    return dm_SPC*(1.0 - k) + dm_HEM*k;
}

f64 hole_N2O_dm(const HoleN2O* hole) {
    HoleN2OU up = hole_N2OU(hole);
    return hole_N2OU_dm(&up, hole->P_d);
}
