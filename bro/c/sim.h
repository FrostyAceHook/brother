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

Optional output buffers (`out_*`): must be supplied by the user, and each must
be either null or contain at-least `out_count` elements. If `out_count` is <=0,
all buffers are assumed to be null. `out_upto` contains the index of the next
write into the buffers, aka the number of elements already written provided it
was set to 0 before simming.

Contour/curve buffers (`*_contour*`, `drag_curve_*`): always exactly filled if
their count is >0.
*/

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
    X(Vempty_c, f64, C_INPUT)                                                   \
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
    X(out_upto, i64, C_INPUT | C_OUTPUT)                                        \
    X(out_count, i64, C_INPUT)                                                  \
    X(out_t, f64*, C_INPUT | C_OUTPUT_DATA)                                     \
    X(out_T_o, f64*, C_INPUT | C_OUTPUT_DATA)                                   \
    X(out_m_o, f64*, C_INPUT | C_OUTPUT_DATA)                                   \
    X(out_rr_f, f64*, C_INPUT | C_OUTPUT_DATA)                                  \
    X(out_alt_r, f64*, C_INPUT | C_OUTPUT_DATA)                                 \
    X(out_vel_r, f64*, C_INPUT | C_OUTPUT_DATA)                                 \


// Da state array (RIP bro sis).
typedef struct simState {
    #define X(name, type, flags) type name;
    SIM_INTERPRETATION
    #undef X
} simState;

// Returns the hash of the canonical interpretation of the state array.
c_IH sim_interpretation_hash(void);

// Simulation entrypoint. Errors are handled via asserts, caller is required to
// setup assertion failed handling.
void sim_execute(simState* rstr s);
