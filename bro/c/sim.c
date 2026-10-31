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



// Some (tweakable) system constants.

static const f64 startup_Dt = 0.1e-3;
static const f64 comb_Dt = 0.4e-3;
static const f64 gotime_Dt = 50e-3;
static const f64 startup_max_t = 0.2; // startup length of all systems.
static const f64 comb_max_t = 50.0; // failsafe.
static const f64 max_t = 1000.0; // failsafe.
static const f64 negligible_m = 2e-3; // assume nothing if mass <= this.
static const f64 comb_min_ofr = 2.0; // combustion extinguishes below.
static const f64 comb_max_ofr = 15.0; // combustion extinguishes above.
static const f64 min_peak_P_c = 5e5; // ballpark.



typedef struct simRunning { // time-dependant constantly-updating parameters.
    i64 onfire;
    f64 t;
    f64 T_o;
    f64 m_o;
    f64 rr_f;
    f64 m_g;
    f64 N_g;
    f64 Cv_g;
    f64 U_g;
    f64 alt_r;
    f64 vel_r;
} simRunning;

typedef struct simPrev { // previous-tick values, for integration.
    f64 Thrust;
    f64 ofr;
} simPrev;


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




static void sim_initial(simState* s, simRunning* rn) {

    assert(oneof(s->static_fire, 0, 1), "static_fire=%lld", s->static_fire);


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


    assert(s->Ot_mov >= 0, "Ot_mov=%g", s->Ot_mov);


    assert(oneof(s->model_inj, HOLE_N2O_SPC, HOLE_N2O_HEM, HOLE_N2O_DYER,
            HOLE_N2O_FML), "model_inj=%lld", s->model_inj);
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


    assert(0 < s->etaCstar && s->etaCstar <= 1, "etaCstar=%g", s->etaCstar);
    assert(0 < s->etaCf && s->etaCf <= 1, "etaCf=%g", s->etaCf);
    /* Vempty_c checked later */


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


    assert(within(s->AEAT, NZL_min_AEAT, NZL_max_AEAT), "AEAT=%g", s->AEAT);
    assert(within(s->NLF, NZL_min_NLF, NZL_max_NLF), "NLF=%g", s->NLF);
    assert(within(s->AEAT, CEA_min_AR, CEA_max_AR), "AEAT=%g", s->AEAT);
    assert(s->D_tht > 0, "D_tht=%g", s->D_tht);
    s->A_tht = circ_area(s->D_tht);
    s->A_exit = s->A_tht * s->AEAT;
    s->D_exit = s->D_tht * sqrt(s->AEAT);
    s->theta_div = nzl_theta_div(s->AEAT, s->NLF);
    s->theta_exit = nzl_theta_exit(s->AEAT, s->NLF);
    assert(0 < s->theta_div && s->theta_div < PI_2,
            "theta_div=%g", s->theta_div);
    assert(0 <= s->theta_exit && s->theta_exit < PI_2,
            "theta_exit=%g", s->theta_exit);


    assert(s->out_count >= 0, "out_count=%lld", s->out_count);
    if (s->out_count == 0) {
        // No storage: all buffers are treated as null.
        s->out_upto = intmin(i64);
        s->out_t = NULL;
        s->out_T_o = NULL;
        s->out_m_o = NULL;
        s->out_rr_f = NULL;
        s->out_m_g = NULL;
        s->out_N_g = NULL;
        s->out_Cv_g = NULL;
        s->out_U_g = NULL;
        s->out_alt_r = NULL;
        s->out_vel_r = NULL;
    } else {
        assert(0 <= s->out_upto, "out_upto=%lld", s->out_upto);
    }


    // Waiting checks:

    // Empty chamber volume must hold the grain with gas space to spare.
    assert(s->Vempty_c > V0_f, "Vempty_c=%g, V0_f=%g", s->Vempty_c, V0_f);

    s->mwet0_r = s->mdry_r + s->m0_o + s->m0_f;



    // Initialise the running state:
    rn->onfire = 1; // start ignited.
    rn->t = 0; // start at the start.
    rn->T_o = s->T0_o; // initial temp given.
    rn->m_o = s->m0_o;
    rn->rr_f = 0;

    // Chamber initially full of ambient air.
    f64 P_a = amb_P(s->alt0_r);
    f64 V0_g = s->Vempty_c - V0_f;
    rn->m_g = P_a * V0_g / AMB_STANCOND_R / AMB_STANCOND_T;
    rn->N_g = rn->m_g / AMB_STANCOND_Mw;
    rn->Cv_g = rn->m_g * AMB_STANCOND_cv;
    rn->U_g = rn->Cv_g * AMB_STANCOND_T;

    // Starts stationary where specified.
    rn->alt_r = s->alt0_r;
    rn->vel_r = 0;
}




// Sets the members of `derivatives` to the current time-derivative of each
// variable, `onfire` is set to its new state (since it is discrete and does not
// have a derivative).
static void sim_diff(simState* s, simRunning* rstr rn,
        simRunning* rstr derivatives) {
    // Unpack current state variables.
    i64 onfire = rn->onfire;
    f64 t     = rn->t;
    f64 T_o   = rn->T_o;
    f64 m_o   = rn->m_o;
    f64 rr_f  = rn->rr_f;
    f64 m_g   = rn->m_g;
    f64 N_g   = rn->N_g;
    f64 Cv_g  = rn->Cv_g;
    f64 U_g   = rn->U_g;
    f64 alt_r = rn->alt_r;
    f64 vel_r = rn->vel_r;

    // Ambient pressure changes with altitude.
    f64 P_a = amb_P(alt_r);

    // Reconstruct some cc/fuel state.
    Grain* grain = get_grain(s, rr_f);
    f64 Aneg_f = grain_Aneg(grain);
    f64 V_f = grain_V(grain, Aneg_f);
    f64 V_c = s->Vempty_c - V_f;
    f64 m_f = s->rho_f * V_f;

    // Reconstruct some cc gas state.
    f64 cv_g = Cv_g / m_g;
    f64 R_g  = GAS_CONSTANT * N_g / m_g;
    f64 cp_g = cv_g + R_g;
    f64 Y_g  = cp_g / cv_g;
    f64 T_g  = U_g / Cv_g;
    f64 P_c  = (onfire) ? (N_g * GAS_CONSTANT * T_g / V_c) : P_a;



    // If its still simming, somethings gone wrong.
    assert(t < ((onfire) ? comb_max_t : max_t), "onfire=%f", onfire);

    // If sim is unstable, generally manifests as cooked cc gas properties.
    assert(T_o > 0.0, "T_o=%f", T_o);
    assert(m_g >= 0.0, "m_g=%f", m_g);
    assert(T_g > 0.0, "T_g=%f", T_g);
    assert(Y_g > 0.0, "Y_g=%f", Y_g);
    assert(R_g > 0.0, "R_g=%f", R_g);


    // To model mov opening, we apply some proportion value to the injector area,
    // where it starts at 0 and ends at 1 at the opening time.
    f64 A_injo_frac = 1.0;
    if (t < s->Ot_mov) {
        // https://www.desmos.com/calculator/3oetvtaoxb
        f64 x = t / s->Ot_mov;
        A_injo_frac = 1.0 - (x - 1.0)*(x - 1.0);
    }
    f64 A_injo = s->A_injo * A_injo_frac;



    // Properties determined around injector:
    f64 P_t, dm_inj, dT_o;

    // Still oxidiser left (note it asymptotes to 0kg so there should always be
    // oxidiser left but anyway).
    if (m_o > 0.0) {
        f64 rho_o = m_o / s->V_t;
        f64 X_o = N2O_X(T_o, rho_o);
        P_t = N2O_P(T_o, rho_o);

        f64 P_d = max(P_c, N2O_Ptrip);

        // What leaves the tank: liquid while any remains, else vapour.
        f64 h_out;
        if (T_o < N2O_Tcrit && X_o < 1.0) {
            f32 compressed = (X_o < 0.0);
            f64 rho_u = (compressed) ? rho_o : N2O_satliq_rho(T_o);
            h_out = (compressed) ? N2O_h(T_o, rho_o) : N2O_satliq_h(T_o);
            dm_inj = hole_N2O_dm(&(HoleN2O){
                .model=s->model_inj,
                .Cd=s->Cd_inj,
                .A=A_injo,
                .P_u=P_t,
                .P_d=P_d,
                .T_u=T_o,
                .rho_u=rho_u,
            });
        } else {
            h_out = N2O_h(T_o, rho_o);
            dm_inj = hole_dm(&(Hole){
                .Cd=s->Cd_inj,
                .A=A_injo,
                .P_u=P_t,
                .P_d=P_d,
                .Y=N2O_Y(T_o, rho_o),
                .R=N2O_R,
                .T=T_o,
                .Z=N2O_Z(T_o, rho_o),
            });
        }

        // First law.
        f64 Cv = m_o * N2O_cv(T_o, rho_o) + s->C_tw;
        dT_o = dm_inj * (N2O_u(T_o, rho_o) - h_out
                       + rho_o * N2O_du_drho_T(T_o, rho_o)) / Cv;
    } else {
        P_t = 0.0;
        dm_inj = 0.0;
        dT_o = 0.0;
    }

    f64 dm_o = -dm_inj;



    // Do fuel regression.
    f64 drr_f;
    f64 dm_reg;
    // Gotta be fire and fuel left to regress.
    if (onfire && (m_f > 0.0)) {

        // Get oxidiser mass flux through the fuel grain.
        f64 G_o = dm_inj / Aneg_f;
        // Use empircally determined regression rate.
        drr_f = grain_drr(G_o);

        // Determine the amount of fuel regressing.
        f64 dAneg_f = grain_dAneg(grain, Aneg_f, drr_f);
        f64 dV_f = grain_dV(grain, Aneg_f, dAneg_f);
        dm_reg = s->rho_f * (-dV_f);

    } else {
        // No fuel left.
        drr_f = 0.0;
        dm_reg = 0.0;
    }



    // Do nozzle flow.
    f64 dm_out = hole_dm(&(Hole){
        .Cd = 1.0,
        .A = s->A_tht,
        .P_u = P_c,
        .P_d = P_a,
        .Y = Y_g,
        .R = R_g,
        .T = T_g,
        .Z = 1.0, // assume ideal gas.
    });



    // Do combustion things.

    f64 Thrust;
    f64 dm_g;
    f64 dN_g;
    f64 dCv_g;
    f64 dU_g;

    // Instantaneous oxidiser-fuel ratio.
    f64 ofr = (dm_reg > 0.0) ? (dm_inj / dm_reg) : 0.0;

    if (onfire) {
        // Tables assert out of domain; clamp to the propellant's domain here.
        f64 P_cea = min(max(P_c, CEA_min_P), CEA_max_P);
        f64 ofr_cea = min(max(ofr, CEA_min_ofr), CEA_max_ofr);

        // Fresh products: adiabatic flame state, with C* efficiency applied as
        // a temperature deficit (C* ~ sqrt(R*T) at fixed Y, Mw).
        f64 dm_n = dm_inj + dm_reg;
        f64 T_n  = cea_T(P_cea, ofr_cea, CEA_AR_stag);
        T_n *= s->etaCstar*s->etaCstar; // apply C* efficiency.
        f64 Mw_n = cea_Mw(P_cea, ofr_cea, CEA_AR_stag);
        f64 Y_n  = cea_Y(P_cea, ofr_cea, CEA_AR_stag);
        f64 R_n  = GAS_CONSTANT / Mw_n;
        f64 cv_n = R_n / (Y_n - 1.0);
        f64 h_n  = (cv_n + R_n) * T_n;

        // Open fixed-volume reservoir; outflow leaves with the mixture's own
        // intensive properties.
        dm_g  = dm_n - dm_out;
        dN_g  = dm_n / Mw_n - dm_out * N_g / m_g;
        dCv_g = dm_n * cv_n - dm_out * cv_g;
        dU_g  = dm_n * h_n - dm_out * cp_g * T_g;

        // Thrust from CEA, with back-pressure loss term and C* efficiency.
        f64 Ivac = cea_Ivac(P_cea, ofr_cea, s->AEAT);
        Thrust = dm_out * Ivac * STANDARD_GRAVITY;
        Thrust *= s->etaCstar * s->etaCf; // efficiency penalty.
        Thrust -= P_a*s->A_exit; // back-pressure force balance.
        // Pretend it doesnt suck the rocket downwards. This happens at extremely
        // poor combustion conditions, when the exhaust is very over-expanded and
        // the higher atmospheric pressure dominates.
        Thrust = max(Thrust, 0.0);

        // Assume that combustion stops when ofr becomes too extreme (but only
        // after startup).
        if (t > startup_max_t)
            onfire = within(ofr, comb_min_ofr, comb_max_ofr);

    } else {
        // Assume no thrust without combustion (i.e. ignore any residual thrust
        // from still-pressurised tank/cc).
        Thrust = 0.0;

        // Ignore cc gas properties if not combusting.
        dm_g = 0.0;
        dN_g = 0.0;
        dCv_g = 0.0;
        dU_g = 0.0;
    }



    // Calculate drag force.
    f64 Drag = rocket_Drag(alt_r, vel_r, s->A_r);

    // Calculate weight force.
    f64 mwet_r = s->mdry_r + m_o + m_f;
    f64 Gravity = -mwet_r * amb_g(alt_r);

    // Calculate acceleration.
    //  F = m*a  [a classic, newton i believe]
    // => a = F/m
    f64 dvel_r = (Thrust + Drag + Gravity) / mwet_r;



    // Got all state time-dependant derivatives, lets pack it into the return
    // struct.
    derivatives->onfire = onfire;
    derivatives->t = 1; // extremely complex derivation, left as an exercise for
                        // the reader.
    derivatives->T_o   = dT_o;
    derivatives->m_o   = dm_o;
    derivatives->rr_f  = drr_f;
    derivatives->m_g   = dm_g;
    derivatives->N_g   = dN_g;
    derivatives->Cv_g  = dCv_g;
    derivatives->U_g   = dU_g;
    derivatives->alt_r = vel_r;
    derivatives->vel_r = dvel_r;
}



enum {
    SIM_KEEP_GOING = 0,
    SIM_STOP = 1,
};
static __attribute((__hot__)) i32 sim_ulate(simState* s, simRunning* rn) {
    // Different time step sizes for different periods.
    f64 Dt = (rn->t <= startup_max_t) ? startup_Dt
           : (rn->onfire) ? comb_Dt
           : gotime_Dt; // gotime.

    // Stop sim when combustion stops.
    i32 stop = (!rn->onfire);
    // Unless not static fire, then also wait for travelling downwards.
    if (!s->static_fire)
        stop &= (rn->vel_r <= 0.0);
    if (stop)
        return SIM_STOP;

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
                .m_g   = prev.m_g   + rk4_Dt * d->m_g,
                .N_g   = prev.N_g   + rk4_Dt * d->N_g,
                .Cv_g  = prev.Cv_g  + rk4_Dt * d->Cv_g,
                .U_g   = prev.U_g   + rk4_Dt * d->U_g,
                .alt_r = prev.alt_r + rk4_Dt * d->alt_r,
                .vel_r = prev.vel_r + rk4_Dt * d->vel_r,
            };
        }
        sim_diff(s, rn, derivatives + i);

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
    f64 Donfire = (k1->onfire || k2->onfire || k3->onfire || k4->onfire);
    f64 dT_o   = (k1->T_o   + 2*k2->T_o   + 2*k3->T_o   + k4->T_o  ) / 6;
    f64 dm_o   = (k1->m_o   + 2*k2->m_o   + 2*k3->m_o   + k4->m_o  ) / 6;
    f64 drr_f  = (k1->rr_f  + 2*k2->rr_f  + 2*k3->rr_f  + k4->rr_f ) / 6;
    f64 dm_g   = (k1->m_g   + 2*k2->m_g   + 2*k3->m_g   + k4->m_g  ) / 6;
    f64 dN_g   = (k1->N_g   + 2*k2->N_g   + 2*k3->N_g   + k4->N_g  ) / 6;
    f64 dCv_g  = (k1->Cv_g  + 2*k2->Cv_g  + 2*k3->Cv_g  + k4->Cv_g ) / 6;
    f64 dU_g   = (k1->U_g   + 2*k2->U_g   + 2*k3->U_g   + k4->U_g  ) / 6;
    f64 dalt_r = (k1->alt_r + 2*k2->alt_r + 2*k3->alt_r + k4->alt_r) / 6;
    f64 dvel_r = (k1->vel_r + 2*k2->vel_r + 2*k3->vel_r + k4->vel_r) / 6;

    // Step me.
    rn->onfire = Donfire;
    rn->t     = prev.t     + Dt;
    rn->T_o   = prev.T_o   + Dt * dT_o;
    rn->m_o   = prev.m_o   + Dt * dm_o;
    rn->rr_f  = prev.rr_f  + Dt * drr_f;
    rn->m_g   = prev.m_g   + Dt * dm_g;
    rn->N_g   = prev.N_g   + Dt * dN_g;
    rn->Cv_g  = prev.Cv_g  + Dt * dCv_g;
    rn->U_g   = prev.U_g   + Dt * dU_g ;
    rn->alt_r = prev.alt_r + Dt * dalt_r;
    rn->vel_r = prev.vel_r + Dt * dvel_r;


    // Ensure we don't enter the ground.
    if (rn->alt_r <= s->alt0_r) {
        rn->alt_r = s->alt0_r;
        rn->vel_r = max(rn->vel_r, 0.0); // pushes against ground.
    }


    // Collapse negligible masses to 0.
    if (rn->m_o <= negligible_m)
        rn->m_o = 0.0;
    Grain* grain = get_grain(s, rn->rr_f);
    if (s->rho_f * grain_V(grain, grain_Aneg(grain)) <= negligible_m) {
        rn->rr_f = grain_max_rr(grain);
        rn->onfire = 0;
        // Note that this onfire culling outside sim_diff is kiiinda dodgy, but
        // it kept having a trailing zero in the ofr without it and like she'll
        // be right.
    }
    // Don't let gas masses slip below 0.
    if (rn->m_g < 0.0)
        rn->m_g = 0.0;


    // TODO: tmp outputs
    assert(s->out_upto < s->out_count, "out_upto=%d, out_count=%d",
            s->out_upto, s->out_count);
    if (s->out_t)      s->out_t[s->out_upto]      = rn->t;
    if (s->out_T_o)    s->out_T_o[s->out_upto]    = rn->T_o;
    if (s->out_m_o)    s->out_m_o[s->out_upto]    = rn->m_o;
    if (s->out_rr_f)   s->out_rr_f[s->out_upto]   = rn->rr_f;
    if (s->out_m_g)    s->out_m_g[s->out_upto]    = rn->m_g;
    if (s->out_N_g)    s->out_N_g[s->out_upto]    = rn->N_g;
    if (s->out_Cv_g)   s->out_Cv_g[s->out_upto]   = rn->Cv_g;
    if (s->out_U_g)    s->out_U_g[s->out_upto]    = rn->U_g;
    if (s->out_alt_r)  s->out_alt_r[s->out_upto]  = rn->alt_r;
    if (s->out_vel_r)  s->out_vel_r[s->out_upto]  = rn->vel_r;
    ++s->out_upto;


    // Not finished.
    return SIM_KEEP_GOING;
}


void sim_execute(simState* rstr s) {
    // Integrate the system differential.
    simRunning* rn = &(simRunning){0};
    sim_initial(s, rn);
    while (sim_ulate(s, rn) != SIM_STOP);
}











// =========================================================================== //
// = INTERPRETATION HASH ===================================================== //
// =========================================================================== //

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
