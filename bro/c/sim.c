#include "sim.h"

#include "assertion.h"
#include "maths.h"
#include "optim.h"

#include "ambient.h"
#include "cea.h"
#include "drag.h"
#include "geometry.h"
#include "grain.h"
#include "injector.h"
#include "n2o.h"
#include "nozzle.h"
#include "relations.h"



// =========================================================================== //
// = CONSTANTS =============================================================== //
// =========================================================================== //

// Some (tweakable) system constants.

static const f64 comb_Dt = 1.0e-3;
static const f64 gotime_Dt = 50e-3;
static const f64 comb_max_t = 50.0; // failsafe.
static const f64 max_t = 1000.0; // failsafe.
static const f64 negligible_P = 1.2e5; // assume no tank pressure if P <= this.
static const f64 negligible_m = 1e-3; // assume nothing if m <= this.



// =========================================================================== //
// = STATE HELPERS =========================================================== //
// =========================================================================== //

typedef struct simRunning { // time-dependant constantly-updating parameters.
    i32 onfire;
    f64 t;
    f64 T_o;
    f64 m_o;
    f64 rr_f;
    f64 alt_r;
    f64 vel_r;

    f64 P_c_hint; // cc pressure solved each step, use prev as seed.
} simRunning;


#define get_grain(args...) ( init_grain(&(Grain){0}, args) )
static Grain* init_grain(Grain* g, simState* s, f64 rr) {
    *g = (Grain){
        .L=s->L0_f,
        .ID=s->OD0_f - 2*s->th0_f,
        .OD=s->OD0_f,
        .rr=rr,
    };
    return g;
}




// =========================================================================== //
// = INJECTOR FLOW =========================================================== //
// =========================================================================== //

// Cached injector flow state (jus like me rn).
typedef struct InjFlow {
    i32 liquid_draining;
    union {
        HoleN2OU liq;
        HoleU gas;
    };
} InjFlow;

static f64 inj_flow(const InjFlow* inj, f64 P_d) {
    return (inj->liquid_draining)
         ? hole_N2OU_dm(&inj->liq, P_d)
         : holeU_dm(&inj->gas, P_d);
}




// =========================================================================== //
// = CC PRESSURE ============================================================= //
// =========================================================================== //

// Quasi-steady chamber: its filling time (~2-3 ms) is ~1000x below the burn
// time, so the chamber is always in mass balance. P_c is the root of
//  resid(P_c) = dm_out(P_c) - dm_inj(P_c) - dm_reg(dm_inj(P_c)),
// which increases with P_c, is <= 0 at P_a (no nozzle flow) and > 0 at P_t (no
// injector flow), so exactly one root lies in [P_a, P_t].

// Everything the balance needs that's fixed for one state.
typedef struct CCctx {
    const simState* s;
    const InjFlow* inj; // null if no oxidiser left.
    f64 P_a;
    i32 regressing;
    f64 Aneg_f;
    f64 dm_per_drr;
} CCctx;

// Chamber state at a candidate pressure.
typedef struct CCstate {
    f64 P_c;
    f64 dm_inj;
    f64 dm_reg;
    f64 drr_f;
    f64 ofr;
    f64 P_cea;   // P_c clamped into the CEA domain
    f64 ofr_cea; // ofr clamped into the CEA domain
    f64 T_c;
    f64 R_c;
    f64 Y_c;
    f64 dm_out;
    f64 resid; // dm_out - dm_inj - dm_reg
} CCstate;

static CCstate cc_eval(const CCctx* c, f64 P_c) {
    const simState* s = c->s;
    CCstate e;
    e.P_c = P_c;
    e.dm_inj = (c->inj) ? inj_flow(c->inj, P_c) : 0.0;

    // Regression from the oxidiser mass flux through the port. Fuel mass flow
    // is linear in the regression rate at fixed geometry, so only the
    // regression law itself is evaluated per candidate.
    if (c->regressing) {
        e.drr_f = grain_drr(e.dm_inj / c->Aneg_f);
        e.dm_reg = c->dm_per_drr * e.drr_f;
    } else {
        e.drr_f = 0.0;
        e.dm_reg = 0.0;
    }
    e.ofr = (e.dm_reg > 0.0) ? (e.dm_inj / e.dm_reg) : 0.0;

    // Products: adiabatic flame state at (P_c, ofr), clamped into the table
    // domain, with C* efficiency as a temperature deficit (C* ~ sqrt(R*T)).
    e.P_cea = min(max(P_c, CEA_min_P), CEA_max_P);
    e.ofr_cea = min(max(e.ofr, CEA_min_ofr), CEA_max_ofr);
    e.T_c = cea_T(e.P_cea, e.ofr_cea, CEA_AR_stag)
          * s->etaCstar * s->etaCstar;
    e.R_c = GAS_CONSTANT / cea_Mw(e.P_cea, e.ofr_cea, CEA_AR_stag);
    e.Y_c = cea_Y(e.P_cea, e.ofr_cea, CEA_AR_stag);

    e.dm_out = hole_dm(&(Hole){
        .CdA = s->A_tht,
        .P_u = P_c,
        .P_d = c->P_a,
        .Y = e.Y_c,
        .R = e.R_c,
        .T = e.T_c,
        .Z = 1.0,
    });
    e.resid = e.dm_out - e.dm_inj - e.dm_reg;
    return e;
}

// Solves for P_c in [P_a, P_t]. `P_hint` (previous solution, NaN for none) only
// affects the iteration count.
static CCstate cc_solve(const CCctx* c, f64 P_t, f64 P_hint) {
    // Bracket over [P_a, P_t].
    CCstate lo = cc_eval(c, c->P_a);
    if (lo.resid >= 0.0 || P_t <= c->P_a)
        return lo; // No injector mass flow, leave cc at P_a.
    CCstate hi = cc_eval(c, P_t);
    assert(hi.resid > 0.0, "cc bracket: resid(P_t=%g)=%g", P_t, hi.resid);

    // Tighten the bracket with the hint.
    if (notnan(P_hint) && c->P_a < P_hint && P_hint < P_t) {
        CCstate m = cc_eval(c, P_hint);
        if (m.resid == 0.0)
            return m;
        if (m.resid < 0.0)
            lo = m;
        else
            hi = m;
    }

    // Illinois false position solving.
    f64 flo = lo.resid;
    f64 fhi = hi.resid;
    enum { SIDE_NONE, SIDE_LO, SIDE_HI };
    i32 side = SIDE_NONE; // no end kept yet, so no halving on the first step.
    for (i32 iter=0; ; ++iter) {
        enum { MAX_ITERS = 60 };

        f64 P = (lo.P_c*fhi - hi.P_c*flo) / (fhi - flo);
        CCstate m = cc_eval(c, P);
        f64 inflow = m.dm_inj + m.dm_reg;

        // Converged when accurate mass balance or tiny bracket.
        if (abs(m.resid) <= 1e-5*inflow || hi.P_c - lo.P_c < 1.0)
            return m;

        if (m.resid > 0.0) {
            hi = m;
            fhi = m.resid;
            if (side == SIDE_HI)
                flo *= 0.5;
            side = SIDE_HI;
        } else {
            lo = m;
            flo = m.resid;
            if (side == SIDE_LO)
                fhi *= 0.5;
            side = SIDE_LO;
        }
        assert(iter < MAX_ITERS, "cc solve failed to converge: P_c in [%g, %g]",
                lo.P_c, hi.P_c);
    }
    __builtin_unreachable();
}




// =========================================================================== //
// = INITIALISE ============================================================== //
// =========================================================================== //

static void sim_initial(simState* s, simRunning* rn) {

    assert(oneof(s->static_fire, 0, 1), "static_fire=%g", s->static_fire);


    assert(notnan(s->alt0_r), "alt0_r=%g", s->alt0_r);
    assert(s->D_r > 0, "D_r=%g", s->D_r);
    s->A_r = circ_area(s->D_r);
    assert(s->mdry_r > 0, "mdry_r=%g", s->mdry_r);
    /* mwet0_r calced later. */


    assert(s->V_t > 0, "V_t=%g", s->V_t);
    assert(s->C_tw >= 0, "C_tw=%g", s->C_tw); // 0 = adiabatic, no wall.


    assert(N2O_Ttrip < s->T0_o && s->T0_o < N2O_Tcrit, "T0_o=%g", s->T0_o);
    assert(notnan(s->m0_o) + notnan(s->vff0_l) == 1,
            "need exactly one of m0_o=%g, vff0_l=%g", s->m0_o, s->vff0_l);
    f64 rho_l = N2O_satliq_rho(s->T0_o);
    f64 rho_v = N2O_satvap_rho(s->T0_o);
    // Note we allow vff outside 0--1.
    if (isnan(s->m0_o)) {
        s->m0_o = s->V_t * lerp(rho_v, rho_l, s->vff0_l);
    } else {
        s->vff0_l = invlerp(rho_v, rho_l, s->m0_o/s->V_t);
    }
    assert(s->m0_o >= 0, "m0_o=%g", s->m0_o);


    assert(0 < s->Cd_inj && s->Cd_inj <= 1, "Cd_inj=%g", s->Cd_inj);
    assert(notnan(s->no_injo) + notnan(s->D_injo) + notnan(s->Deq_injo) == 2,
            "need exactly two of no_injo=%g, D_injo=%g, Deq_injo=%g", s->no_injo,
            s->D_injo, s->Deq_injo);
    if (notnan(s->no_injo))
        assert(s->no_injo >= 1 && s->no_injo == floor(s->no_injo),
                "no_injo=%g", s->no_injo);
    if (notnan(s->D_injo))
        assert(s->D_injo > 0, "D_injo=%g", s->D_injo);
    if (notnan(s->Deq_injo))
        assert(s->Deq_injo > 0, "Deq_injo=%g", s->Deq_injo);
    if (isnan(s->no_injo))
        s->no_injo = sqed(s->Deq_injo/s->D_injo);
    else if (isnan(s->D_injo))
        s->D_injo = s->Deq_injo / sqrt(s->no_injo);
    else
        s->Deq_injo = s->D_injo * sqrt(s->no_injo);
    s->A_injo = circ_area(s->Deq_injo);
    s->CdA_inj = s->Cd_inj*s->A_injo;


    assert(0 < s->etaCstar && s->etaCstar <= 1, "etaCstar=%g", s->etaCstar);
    assert(0 < s->etaCf && s->etaCf <= 1, "etaCf=%g", s->etaCf);


    assert(s->L0_f > 0, "L0_f=%g", s->L0_f);
    assert(s->th0_f > 0 && 2*s->th0_f < s->OD0_f, "th0_f=%g, OD0_f=%g", s->th0_f,
            s->OD0_f);
    assert(notnan(s->rho_f) + notnan(s->m0_f) == 1,
            "need exactly one of rho_f=%g, m0_f=%g", s->rho_f, s->m0_f);
    Grain* grain = get_grain(s, 0);
    f64 V0_f = grain_V(grain, grain_Aneg(grain));
    if (isnan(s->rho_f))
        s->rho_f = s->m0_f / V0_f;
    else
        s->m0_f = s->rho_f * V0_f;
    assert(s->rho_f > 0, "rho_f=%g", s->rho_f);
    assert(s->m0_f > 0, "m0_f=%g", s->m0_f);


    assert(notnan(s->D_tht) + notnan(s->A_tht) == 1,
            "need exactly one of D_tht=%g, A_tht=%g", s->D_tht, s->A_tht);
    if (isnan(s->A_tht))
        s->A_tht = circ_area(s->D_tht);
    else
        s->D_tht = circ_diam(s->A_tht);
    assert(s->D_tht > 0, "D_tht=%g", s->D_tht);
    assert(s->A_tht > 0, "A_tht=%g", s->A_tht);
    assert(notnan(s->D_exit) + notnan(s->A_exit) + notnan(s->AEAT) == 1,
            "need exactly one of D_exit=%g, A_exit=%g, AEAT=%g", s->D_exit,
            s->A_exit, s->AEAT);
    if (notnan(s->D_exit)) {
        s->A_exit = circ_area(s->D_exit);
        s->AEAT = s->A_exit / s->A_tht;
    } else if (notnan(s->A_exit)) {
        s->D_exit = circ_diam(s->A_exit);
        s->AEAT = s->A_exit / s->A_tht;
    } else {
        s->A_exit = s->AEAT * s->A_tht;
        s->D_exit = circ_diam(s->A_exit);
    }
    assert(s->D_exit > 0, "D_exit=%g", s->D_exit);
    assert(s->A_exit > 0, "A_exit=%g", s->A_exit);
    assert(s->AEAT > 1, "AEAT=%g", s->AEAT);
    assert(0 < s->NLF && s->NLF <= 1, "NLF=%g", s->NLF);
    assert(within(s->AEAT, NZL_min_AEAT, NZL_max_AEAT), "AEAT=%g", s->AEAT);
    assert(within(s->AEAT, CEA_min_AR, CEA_max_AR), "AEAT=%g", s->AEAT);
    assert(within(s->NLF, NZL_min_NLF, NZL_max_NLF), "NLF=%g", s->NLF);
    s->theta_div = nzl_theta_div(s->AEAT, s->NLF);
    s->theta_exit = nzl_theta_exit(s->AEAT, s->NLF);
    assert(0 < s->theta_div && s->theta_div < PI_2,
            "theta_div=%g", s->theta_div);
    assert(0 <= s->theta_exit && s->theta_exit < PI_2,
            "theta_exit=%g", s->theta_exit);


    assert(s->snap_count >= 0, "snap_count=%lld", s->snap_count);
    if (s->snap_count) {
        assert(0 <= s->snap_upto, "snap_upto=%lld", s->snap_upto);
        assert(s->snap_array, "null pointer");
    }


    // Waiting checks:

    s->mwet0_r = s->mdry_r + s->m0_o + s->m0_f;



    // Initialise the running state:
    rn->onfire = 1; // start ignited.
    rn->t = 0; // start at the start.
    rn->T_o = s->T0_o; // initial temp given.
    rn->m_o = s->m0_o;
    rn->rr_f = 0;
    rn->alt_r = s->alt0_r; // launch altitude given.
    rn->vel_r = 0; // start stationary (dont model the onboard we sacrificed
                   // giving it a push).
    rn->P_c_hint = NAN; // no hint initially.
}




// =========================================================================== //
// = DERIVATIVES ============================================================= //
// =========================================================================== //

static void sim_diff(simState* s, simRunning* rstr rn, simRunning* rstr d_dt,
        simSnapshot* ss /* may be null */) {
    // Unpack current state variables.
    i32 onfire = rn->onfire;
    f64 t     = rn->t; (void)t; // not used lmao.
    f64 T_o   = rn->T_o;
    f64 m_o   = rn->m_o;
    f64 rr_f  = rn->rr_f;
    f64 alt_r = rn->alt_r;
    f64 vel_r = rn->vel_r;

    // If sim is unstable, generally manifests as cooked tank properties.
    assert(T_o > 0.0, "T_o=%f", T_o);


    // Tank state, and the cached injector flow state for it.
    f64 P_t = 0.0;
    f64 rho_o = NAN;
    f64 X_o = NAN;
    f64 h_out = NAN;
    InjFlow inj;
    // Still oxidiser left (note it asymptotes to 0kg so there should always be
    // oxidiser left but anyway).
    i32 have_ox = (m_o > 0.0);
    if (have_ox) {
        rho_o = m_o / s->V_t;
        X_o = N2O_X(T_o, rho_o);
        P_t = N2O_P(T_o, rho_o);

        // What leaves the tank: liquid while any remains, else vapour.
        inj.liquid_draining = (T_o < N2O_Tcrit && X_o < 1.0);
        if (inj.liquid_draining) {
            i32 compressed = (X_o < 0.0);
            f64 rho_u = (compressed) ? rho_o : N2O_satliq_rho(T_o);
            h_out = (compressed) ? N2O_h(T_o, rho_o) : N2O_satliq_h(T_o);
            inj.liq = hole_N2OU(&(HoleN2O){
                .model = HOLE_N2O_FML, // TODO: dont hardcode.
                .CdA = s->CdA_inj,
                .P_u = P_t,
                .T_u = T_o,
                .rho_u = rho_u,
            });
        } else {
            h_out = N2O_h(T_o, rho_o);
            inj.gas = holeU(&(Hole){
                .CdA = s->CdA_inj,
                .P_u = P_t,
                .Y = N2O_Y(T_o, rho_o),
                .R = N2O_R,
                .T = T_o,
                .Z = N2O_Z(T_o, rho_o),
            });
        }
    }


    // Reconstruct some fuel state.
    Grain* grain = get_grain(s, rr_f);
    f64 Aneg_f = grain_Aneg(grain);
    f64 V_f = grain_V(grain, Aneg_f);
    f64 m_f = s->rho_f * V_f;

    // Fuel mass flow per unit regression rate [kg/m] = rho_f * burning area.
    // The grain functions are linear in the regression rate at fixed geometry,
    // so evaluate the mass flow [kg/s] at a reference rate [m/s] and divide it
    // back out: kg/s / (m/s) = kg/m. (drr_ref's value is arbitrary.)
    f64 drr_f_ref = 1.0; // [m/s]
    f64 dm_per_drr = s->rho_f
                   * -grain_dV(grain, Aneg_f,
                               grain_dAneg(grain, Aneg_f, drr_f_ref))
                   / drr_f_ref;


    // Ambient pressure changes with altitude.
    f64 P_a = amb_P(alt_r);


    // Chamber: quasi-steady mass balance while burning, else ambient.
    CCctx cc_ctx = {
        .s = s,
        .inj = (have_ox) ? &inj : NULL,
        .P_a = P_a,
        .regressing = onfire,
        .Aneg_f = Aneg_f,
        .dm_per_drr = dm_per_drr,
    };
    CCstate cc;
    if (onfire) {
        cc = cc_solve(&cc_ctx, P_t, rn->P_c_hint);
        rn->P_c_hint = cc.P_c; // hint future steps.
    } else {
        // Assume chamber at ambient.
        cc = cc_eval(&cc_ctx, P_a);
    }
    f64 dm_inj = cc.dm_inj;



    // Tank first law (closed form, valid in every phase).
    f64 dT_o = 0.0;
    if (have_ox) {
        f64 Cv = m_o * N2O_cv(T_o, rho_o) + s->C_tw;
        dT_o = dm_inj * (N2O_u(T_o, rho_o) - h_out
                       + rho_o * N2O_du_drho_T(T_o, rho_o)) / Cv;
    }
    f64 dm_o = -dm_inj;



    // Calc me thrust.
    f64 Thrust;
    if (onfire) {
        // Thrust from CEA, with back-pressure loss term and C* efficiency.
        f64 Ivac = cea_Ivac(cc.P_cea, cc.ofr_cea, s->AEAT);
        Thrust = cc.dm_out * Ivac * STANDARD_GRAVITY;
        Thrust *= s->etaCstar * s->etaCf; // efficiency penalty.
        Thrust -= P_a*s->A_exit; // back-pressure force balance.
        // Pretend it doesnt suck the rocket downwards. This happens at extremely
        // poor combustion conditions, when the exhaust is very over-expanded and
        // the higher atmospheric pressure dominates.
        Thrust = max(Thrust, 0.0);
    } else {
        // Assume no thrust without combustion (i.e. ignore any residual thrust
        // from still-pressurised tank/cc).
        Thrust = 0.0;
    }



    // Calculate drag force.
    f64 Drag = drag_force(alt_r, vel_r, s->A_r);

    // Calculate weight force.
    f64 mwet_r = s->mdry_r + m_o + m_f;
    f64 Gravity = -mwet_r * amb_g(alt_r);

    // Calculate acceleration.
    //  F = m*a  [a classic, newton i believe]
    // => a = F/m
    f64 dvel_r = (Thrust + Drag + Gravity) / mwet_r;



    // Got all state time-dependant derivatives, lets pack it into the return
    // struct.
    d_dt->t = 1; // extremely complex derivation, left as an exercise for the
                 // reader.
    d_dt->T_o   = dT_o;
    d_dt->m_o   = dm_o;
    d_dt->rr_f  = cc.drr_f;
    d_dt->alt_r = vel_r;
    d_dt->vel_r = dvel_r;


    // Fill out full snapshot if requested.
    if (!ss)
        return;

    // Tank phase split. Liquid and vapour don't exist supercritically (X is
    // NaN there), so those fields are NaN.
    f64 rho_l = NAN;
    f64 rho_v = NAN;
    f64 m_l = NAN;
    f64 vff_l = NAN;
    if (!have_ox) {
        m_l = 0.0;
        vff_l = 0.0;
    } else if (notnan(X_o)) {
        if (X_o <= 0.0) { // liquid.
            rho_l = rho_o;
            m_l = m_o;
            vff_l = 1.0;
        } else if (X_o >= 1.0) { // vapour.
            rho_v = rho_o;
            m_l = 0.0;
            vff_l = 0.0;
        } else { // saturated mixture.
            rho_l = N2O_satliq_rho(T_o);
            rho_v = N2O_satvap_rho(T_o);
            m_l = (1.0 - X_o) * m_o;
            vff_l = m_l / rho_l / s->V_t;
        }
    }

    // Chamber gas: ideal gas (frozen Y, R, with C* efficiency in T), so the
    // throat follows from isentropic relations at choke.
    f64 Y = cc.Y_c;
    f64 R = cc.R_c;
    f64 T_cc = (onfire) ? cc.T_c : NAN;
    f64 P_cc = (onfire) ? cc.P_c : NAN;
    f64 rho_cc = P_cc/R/T_cc; // TODO: add to relations
    f64 a_cc = sqrt(Y*R*T_cc); // TODO: ^
    f64 cv_g = R / (Y - 1.0);
    SpecificHeatRatio* shr = get_shr(Y);

    // Throat plane.
    f64 T_tht = NAN;
    f64 P_tht = NAN;
    f64 rho_tht = NAN;
    f64 a_tht = NAN;
    if (onfire) {
        T_tht   = T_cc   * isentropic_T_on_T0(1.0, shr);
        P_tht   = P_cc   * isentropic_P_on_P0(1.0, shr);
        rho_tht = rho_cc * isentropic_rho_on_rho0(1.0, shr);
        a_tht   = sqrt(Y*R*T_tht);
    }

    // Exit plane.
    f64 T_exit = NAN;
    f64 P_exit = NAN;
    f64 rho_exit = NAN;
    f64 a_exit = NAN;
    f64 M_exit = NAN;
    if (onfire) {
        // Use isentropic relations smile.
        f64 M = isentropic_sup_M(s->AEAT, shr);
        T_exit   = T_cc   * isentropic_T_on_T0(M, shr);
        P_exit   = P_cc   * isentropic_P_on_P0(M, shr);
        rho_exit = rho_cc * isentropic_rho_on_rho0(M, shr);
        a_exit   = sqrt(Y*R*T_exit);
        M_exit   = M;
    }


    ss->t = rn->t;
    ss->alt_r = rn->alt_r;
    ss->vel_r = rn->vel_r;
    ss->acc_r = (s->static_fire) ? 0 : dvel_r;
    ss->mwet_r = mwet_r;
    ss->P_t = P_t;
    ss->P_c = cc.P_c;
    ss->T_o = T_o;
    ss->m_o = m_o;
    ss->m_l = m_l;
    ss->m_v = m_o - m_l;
    ss->rho_l = rho_l;
    ss->rho_v = rho_v;
    ss->vff_l = vff_l;
    ss->X_o = X_o;
    ss->s_o = (have_ox) ? N2O_s(T_o, rho_o) : NAN;
    ss->h_o = (have_ox) ? N2O_h(T_o, rho_o) : NAN;
    ss->u_o = (have_ox) ? N2O_u(T_o, rho_o) : NAN;
    ss->cp_o = (have_ox) ? N2O_cp(T_o, rho_o) : NAN; // NaN in the dome.
    ss->cv_o = (have_ox) ? N2O_cv(T_o, rho_o) : NAN; // mixture cv in the dome.
    ss->Y_o = (have_ox) ? N2O_Y(T_o, rho_o) : NAN; // NaN in the dome.
    ss->m_f = m_f;
    ss->rr_f = rr_f;
    ss->Pr_inj = (cc.P_c > 0.0) ? P_t / cc.P_c : NAN; // tank/chamber.
    ss->dm_inj = cc.dm_inj;
    ss->dm_reg = cc.dm_reg;
    ss->dm_out = cc.dm_out;
    ss->T_cc = T_cc;
    ss->rho_cc = rho_cc;
    ss->a_cc = a_cc;
    ss->Y_cc = (onfire) ? Y : NAN;
    ss->cp_cc = (onfire) ? cv_g + R : NAN;
    ss->cv_cc = (onfire) ? cv_g : NAN;
    ss->Mw_cc = (onfire) ? GAS_CONSTANT / R : NAN;
    ss->R_cc = (onfire) ? R : NAN;
    ss->T_tht = T_tht;
    ss->P_tht = P_tht;
    ss->rho_tht = rho_tht;
    ss->a_tht = a_tht;
    ss->T_exit = T_exit;
    ss->P_exit = P_exit;
    ss->rho_exit = rho_exit;
    ss->a_exit = a_exit;
    ss->M_exit = M_exit;
    ss->ofr = (onfire) ? cc.ofr : NAN;
    ss->Isp = (cc.dm_out > 0.0) ? Thrust/cc.dm_out/STANDARD_GRAVITY : NAN;
    ss->Thrust = (onfire) ? Thrust : NAN;
    ss->Drag = Drag;
    ss->Gravity = Gravity;
    ss->Cstar = (onfire) ? cc.P_c * s->A_tht / cc.dm_out : NAN;
    ss->Cf = (onfire) ? Thrust / cc.P_c / s->A_tht : NAN;
}



// =========================================================================== //
// = SIMULATION ============================================================== //
// =========================================================================== //

enum {
    SIM_KEEP_GOING = 0,
    SIM_STOP = 1,
};
static __attribute((__hot__)) i32 sim_ulate(simState* s, simRunning* rn) {
    // Different time step sizes for different periods.
    f64 Dt = (rn->onfire) ? comb_Dt : gotime_Dt; // gotime.
    f64 this_max_t = (rn->onfire) ? comb_max_t : max_t;
    assert(rn->t < this_max_t, "t=%g, this_max_t=%g", rn->t, this_max_t);


    // Stop sim when combustion stops.
    i32 stop = (!rn->onfire);
    // Don't stop if non-empty tank.
    if (rn->m_o > 0) {
        f64 rho_o = rn->m_o / s->V_t;
        f64 P_t = N2O_P(rn->T_o, rho_o);
        stop &= (P_t < negligible_P);
    }
    // Unless not static fire, then also wait for travelling downwards.
    if (!s->static_fire)
        stop &= (rn->vel_r <= 0.0);
    if (stop)
        return SIM_STOP;


    // Fill out snapshot if needed.
    simSnapshot snap;
    i32 recording = (s->snap_count > 0);

    // Do the four state derivatives that are needed for runge-kutta 4
    // integration.
    simRunning derivatives[4]; // all members are derivatives.
    simRunning prev = *rn;
    for (i32 i=0; i<4; ++i) {
        if (i > 0) {
            // Setup new input, based on rk4 formula.
            f64 rk4_Dt = (i == 3) ? Dt : 0.5 * Dt;
            simRunning* d = derivatives + i - 1;
            *rn = (simRunning){
                .onfire = prev.onfire,
                .t     = prev.t     + rk4_Dt,
                .T_o   = prev.T_o   + rk4_Dt * d->T_o,
                .m_o   = prev.m_o   + rk4_Dt * d->m_o,
                .rr_f  = prev.rr_f  + rk4_Dt * d->rr_f,
                .alt_r = prev.alt_r + rk4_Dt * d->alt_r,
                .vel_r = prev.vel_r + rk4_Dt * d->vel_r,
                .P_c_hint = rn->P_c_hint, // keep latest solve as hint.
            };
        }
        sim_diff(s, rn, derivatives + i, (i == 0 && recording) ? &snap : NULL);

        // Never moves under static fire.
        if (s->static_fire)
            derivatives[i].vel_r = 0.0;
    }

    // Apply rk4 formula to find a weighted average of the four derivatives we
    // just calculated and use that for this entire step.
    simRunning* k1 = derivatives + 0;
    simRunning* k2 = derivatives + 1;
    simRunning* k3 = derivatives + 2;
    simRunning* k4 = derivatives + 3;
    f64 dT_o   = (k1->T_o   + 2*k2->T_o   + 2*k3->T_o   + k4->T_o  ) / 6;
    f64 dm_o   = (k1->m_o   + 2*k2->m_o   + 2*k3->m_o   + k4->m_o  ) / 6;
    f64 drr_f  = (k1->rr_f  + 2*k2->rr_f  + 2*k3->rr_f  + k4->rr_f ) / 6;
    f64 dalt_r = (k1->alt_r + 2*k2->alt_r + 2*k3->alt_r + k4->alt_r) / 6;
    f64 dvel_r = (k1->vel_r + 2*k2->vel_r + 2*k3->vel_r + k4->vel_r) / 6;

    // Step me.
    rn->t     = prev.t     + Dt;
    rn->T_o   = prev.T_o   + Dt * dT_o;
    rn->m_o   = prev.m_o   + Dt * dm_o;
    rn->rr_f  = prev.rr_f  + Dt * drr_f;
    rn->alt_r = prev.alt_r + Dt * dalt_r;
    rn->vel_r = prev.vel_r + Dt * dvel_r;


    // Ensure we don't enter the ground.
    if (rn->alt_r <= s->alt0_r) {
        rn->alt_r = s->alt0_r;
        rn->vel_r = max(rn->vel_r, 0.0); // pushes against ground.
    }


    // Collapse negligible masses to 0.
    if (rn->m_o <= negligible_m) {
        rn->m_o = 0;
        rn->onfire = 0;
    }
    Grain* grain = get_grain(s, rn->rr_f); // re-calc.
    if (s->rho_f * grain_V(grain, grain_Aneg(grain)) <= negligible_m) {
        rn->rr_f = grain_max_rr(grain);
        rn->onfire = 0;
    }


    // Extinguish once tank pressure barely above ambient.
    if (rn->m_o > 0 && rn->onfire) {
        f64 rho_o = rn->m_o / s->V_t;
        f64 P_t = N2O_P(rn->T_o, rho_o);
        if (P_t < negligible_P)
            rn->onfire = 0;
    }


    // Output snapshot if requested.
    if (s->snap_count > 0) {
        assert(s->snap_upto < s->snap_count, "snap_upto=%lld, snap_count=%lld",
                s->snap_upto, s->snap_count);
        memcpy(s->snap_array + s->snap_upto*SIM_SNAPSHOT_N, &snap, sizeof(snap));
        ++s->snap_upto;
    }


    // Not finished.
    return SIM_KEEP_GOING;
}


void sim_execute(simState* rstr s) {
    assert(s->snap_hash == sim_snapshot_hash(),
            "snapshot hash does not match, proposal is dismissed");

    // Integrate the system differential, tracking a running state.
    simRunning* rn = &(simRunning){0};
    sim_initial(s, rn);
    while (sim_ulate(s, rn) != SIM_STOP);
}




// =========================================================================== //
// = HASHES ================================================================== //
// =========================================================================== //

c_IH sim_snapshot_hash(void) {
    c_IH h = c_ih_initial();
    #define X(name) \
        h = c_ih_append(h, #name, C_F64 | C_OUTPUT);
    SIM_SNAPSHOT
    #undef X
    return h;
}

c_IH sim_interpretation_hash(void) {
    c_IH h = c_ih_initial();
    #define X(name, type, flags)                                \
        h = c_ih_append(h, #name, flags | generic(objof(type)   \
            , f64:  C_F64                                       \
            , i64:  C_I64                                       \
            , f32*: C_PTR_F32                                   \
            , f64*: C_PTR_F64                                   \
            , i8*:  C_PTR_I8                                    \
            , i16*: C_PTR_I16                                   \
            , i32*: C_PTR_I32                                   \
            , i64*: C_PTR_I64                                   \
            , u8*:  C_PTR_U8                                    \
            , u16*: C_PTR_U16                                   \
            , u32*: C_PTR_U32                                   \
            , u64*: C_PTR_U64                                   \
        ));
    SIM_INTERPRETATION
    #undef X
    return h;
}
