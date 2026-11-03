#pragma once
#include "br.h"

#include "../bridge/bridge.h"



/*
Cheeky rocket ascii:

   .     .
   :     :
   |     | ----- locked (by rest of team)
   | ___ |
   ||   ||
   ||   || ----- tank (contains oxidiser)
   ||___||
   |  X  | ----- mov (main oxidiser valve)
   | _#_ | ----- injector
   ||, ,||
   ||| ||| ----- cc (combustion chamber, contains fuel grain)
  /||' '||\
 / | \ / | \
/  |_/ \_|  \
|_/   |   \_| -- fins
      |
      '--------- nozzle

Only the independant parameters of the system are taken an inputs, while all
outputs are dependant. Note that this choice of "basis" variables (the
independant parameters) is arbitrary, but some make significantly more sense than
others (i.e. the relationship is difficult to inverse).

Optional output buffers (`snap_*`): must be supplied by the user, and all must
contain at-least `snap_count` elements. If `snap_count` is <=0, all buffers are
assumed to be null and will not be filled. `snap_upto` contains the index of the
next write into the buffers, aka the number of elements already written provided
it was set to 0 before simming.
*/


#define SIM_SNAPSHOT    \
    X(t)                \
    X(alt_r)            \
    X(vel_r)            \
    X(acc_r)            \
    X(mwet_r)           \
    X(P_t)              \
    X(P_c)              \
    X(T_o)              \
    X(m_o)              \
    X(m_l)              \
    X(m_v)              \
    X(rho_l)            \
    X(rho_v)            \
    X(vff_l)            \
    X(X_o)              \
    X(s_o)              \
    X(h_o)              \
    X(u_o)              \
    X(cp_o)             \
    X(cv_o)             \
    X(Y_o)              \
    X(m_f)              \
    X(rr_f)             \
    X(Pr_inj)           \
    X(dm_inj)           \
    X(dm_reg)           \
    X(dm_out)           \
    X(T_cc)             \
    X(rho_cc)           \
    X(a_cc)             \
    X(Y_cc)             \
    X(cp_cc)            \
    X(cv_cc)            \
    X(Mw_cc)            \
    X(R_cc)             \
    X(T_tht)            \
    X(P_tht)            \
    X(rho_tht)          \
    X(a_tht)            \
    X(T_exit)           \
    X(P_exit)           \
    X(rho_exit)         \
    X(a_exit)           \
    X(M_exit)           \
    X(ofr)              \
    X(Isp)              \
    X(Thrust)           \
    X(Drag)             \
    X(Gravity)          \
    X(Cstar)            \
    X(Cf)               \


typedef struct simSnapshot {
    #define X(name) f64 name;
    SIM_SNAPSHOT
    #undef X
} simSnapshot;

enum { SIM_SNAPSHOT_N = sizeof(simSnapshot) / sizeof(f64) };
static_assert(sizeof(simSnapshot) == SIM_SNAPSHOT_N*sizeof(f64));



// Canonical interpretation of the state array.
#define SIM_INTERPRETATION                                                      \
    X(static_fire, f64, C_INPUT)                                                \
                                                                                \
    X(alt0_r, f64, C_INPUT)                                                     \
    X(D_r, f64, C_INPUT)                                                        \
    X(A_r, f64, C_OUTPUT)                                                       \
    X(mdry_r, f64, C_INPUT)                                                     \
    X(mwet0_r, f64, C_OUTPUT)                                                   \
                                                                                \
    X(V_t, f64, C_INPUT)                                                        \
    X(C_tw, f64, C_INPUT)                                                       \
                                                                                \
    X(T0_o, f64, C_INPUT)                                                       \
    X(m0_o, f64, C_INPUT | C_OUTPUT)                                            \
    X(vff0_l, f64, C_INPUT | C_OUTPUT)                                          \
                                                                                \
    X(Cd_inj, f64, C_INPUT)                                                     \
    X(no_injo, f64, C_INPUT | C_OUTPUT)                                         \
    X(D_injo, f64, C_INPUT | C_OUTPUT)                                          \
    X(Deq_injo, f64, C_INPUT | C_OUTPUT)                                        \
    X(A_injo, f64, C_OUTPUT)                                                    \
    X(CdA_inj, f64, C_OUTPUT)                                                   \
                                                                                \
    X(etaCstar, f64, C_INPUT)                                                   \
    X(etaCf, f64, C_INPUT)                                                      \
                                                                                \
    X(L0_f, f64, C_INPUT)                                                       \
    X(OD0_f, f64, C_INPUT)                                                      \
    X(th0_f, f64, C_INPUT)                                                      \
    X(rho_f, f64, C_INPUT | C_OUTPUT)                                           \
    X(m0_f, f64, C_INPUT | C_OUTPUT)                                            \
                                                                                \
    X(D_tht, f64, C_INPUT | C_OUTPUT)                                           \
    X(A_tht, f64, C_INPUT | C_OUTPUT)                                           \
    X(D_exit, f64, C_INPUT | C_OUTPUT)                                          \
    X(A_exit, f64, C_INPUT | C_OUTPUT)                                          \
    X(AEAT, f64, C_INPUT | C_OUTPUT)                                            \
    X(NLF, f64, C_INPUT)                                                        \
    X(theta_div, f64, C_OUTPUT)                                                 \
    X(theta_exit, f64, C_OUTPUT)                                                \
                                                                                \
    X(snap_upto, i64, C_INPUT | C_OUTPUT)                                       \
    X(snap_count, i64, C_INPUT)                                                 \
    X(snap_hash, i64, C_INPUT)                                                  \
    X(snap_array, f64 MAYALIAS*, C_INPUT | C_OUTPUT_DATA)                       \



// Da state array (RIP bro sis).
typedef struct simState {
    #define X(name, type, flags) type name;
    SIM_INTERPRETATION
    #undef X
} simState;

// Returns the hash of the canonical interpretation of the state array.
c_IH sim_interpretation_hash(void);
// Returns the hash of the canonical interpretation of the snapshot.
c_IH sim_snapshot_hash(void);

// Simulation entrypoint. Errors are handled via asserts, caller is required to
// setup assertion failed handling.
void sim_execute(simState* rstr s);
