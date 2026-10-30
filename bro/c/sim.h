#pragma once
#include "br.h"

#include "../bridge/bridge.h"


// Canonical interpretation of the state array.
#define SIM_INTERPRETATION                                      \
    X(x, f64, C_INPUT)


// Da state array.
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
