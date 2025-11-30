/**
 * @file rk4_solver.c
 * @brief C-RIVE RK4 Adaptive Solver Implementation
 * 
 * 4th-order Runge-Kutta solver for CGEM biogeochemistry.
 * Ported from C-RIVE manage_simulation.c, calc_new_C_comp(), calc_final_conc()
 * 
 * CITATION:
 * Wang, S., Flipo, N., Romary, T., 2018. Time-dependent global sensitivity analysis
 * of the C-RIVE biogeochemical model. Water Research 144, 341-355.
 * 
 * COPYRIGHT: Original C-RIVE code (c) 2023 Contributors to the librive library.
 * Eclipse Public License v2.0
 */

#include "rk4_solver.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ===========================================================================
 * Configuration and Allocation
 * ===========================================================================*/

void rk4_init_config(RK4Config *config) {
    if (!config) return;
    
    strcpy(config->name, "Runge-Kutta 4");
    config->kmax = RK4_STAGES;
    config->max_div_dt = RK4_MAX_DIV_DT;
    
    /* Copy Butcher tableau coefficients */
    for (int i = 0; i < RK4_STAGES; i++) {
        config->dt_RK[i] = RK4_DT_FRAC[i];
        config->coef_RK[i] = RK4_COEF[i];
    }
}

RK4State *rk4_alloc_state(int num_species) {
    if (num_species <= 0) return NULL;
    
    RK4State *state = (RK4State *)calloc(1, sizeof(RK4State));
    if (!state) return NULL;
    
    state->num_species = num_species;
    
    state->C = (double *)calloc((size_t)num_species, sizeof(double));
    state->newC = (double *)calloc((size_t)num_species, sizeof(double));
    state->C_old = (double *)calloc((size_t)num_species, sizeof(double));
    state->mb = (RK4MassBalance *)calloc((size_t)num_species, sizeof(RK4MassBalance));
    
    if (!state->C || !state->newC || !state->C_old || !state->mb) {
        rk4_free_state(state);
        return NULL;
    }
    
    return state;
}

void rk4_free_state(RK4State *state) {
    if (!state) return;
    
    if (state->C) free(state->C);
    if (state->newC) free(state->newC);
    if (state->C_old) free(state->C_old);
    if (state->mb) free(state->mb);
    
    free(state);
}

void rk4_reset_mb(RK4State *state) {
    if (!state || !state->mb) return;
    
    for (int s = 0; s < state->num_species; s++) {
        for (int k = 0; k < RK4_STAGES; k++) {
            state->mb[s].dC[k] = 0.0;
            for (int r = 0; r < CGEM_NUM_REACTIONS; r++) {
                state->mb[s].dmb[k][r] = 0.0;
            }
        }
        state->mb[s].deltaC = 0.0;
        for (int r = 0; r < CGEM_NUM_REACTIONS; r++) {
            state->mb[s].deltamb[r] = 0.0;
        }
    }
}

/* ===========================================================================
 * RK4 Stage Updates - Direct port from C-RIVE calc_new_C_comp()
 * ===========================================================================*/

void rk4_update_newC(RK4State *state, const RK4Config *config, double dt, int k) {
    if (!state || !config || k < 0 || k >= RK4_STAGES) return;
    
    if (k == 0) {
        /* First stage: newC = C_old (no modification) */
        for (int s = 0; s < state->num_species; s++) {
            state->newC[s] = state->C_old[s];
        }
    } else {
        /* Subsequent stages: newC = C_old + dt * dt_RK[k] * dC[k-1] */
        double dt_frac = dt * config->dt_RK[k];
        for (int s = 0; s < state->num_species; s++) {
            state->newC[s] = state->C_old[s] + dt_frac * state->mb[s].dC[k - 1];
        }
    }
}

/* ===========================================================================
 * RK4 Finalization - Direct port from C-RIVE calc_final_conc()
 * ===========================================================================*/

void rk4_finalize(RK4State *state, const RK4Config *config, double dt) {
    if (!state || !config) return;
    
    for (int s = 0; s < state->num_species; s++) {
        /* Accumulate weighted stage derivatives */
        double deltaC = 0.0;
        for (int k = 0; k < config->kmax; k++) {
            deltaC += config->coef_RK[k] * state->mb[s].dC[k];
        }
        state->mb[s].deltaC = dt * deltaC;
        
        /* Update final concentration */
        state->C[s] = state->C_old[s] + state->mb[s].deltaC;
        
        /* Accumulate mass balance terms */
        for (int r = 0; r < CGEM_NUM_REACTIONS; r++) {
            double deltamb = 0.0;
            for (int k = 0; k < config->kmax; k++) {
                deltamb += config->coef_RK[k] * state->mb[s].dmb[k][r];
            }
            state->mb[s].deltamb[r] = dt * deltamb;
        }
    }
}

/* ===========================================================================
 * Adaptive Time Step - From C-RIVE calc_new_dt()
 * ===========================================================================*/

double rk4_calc_adaptive_dt(const RK4State *state, double dt, const RK4Config *config) {
    if (!state || !config) return dt;
    
    double dt_new = dt;
    
    /* Check each species for potential instability */
    for (int s = 0; s < state->num_species; s++) {
        double C = state->newC[s];
        double dC = state->mb[s].dC[0];  /* Use first stage derivative */
        
        if (fabs(dC) > RK4_EPS && C > RK4_EPS) {
            /* CFL-like condition: dt should not deplete concentration */
            double dt_species = fabs(C / dC) * 0.5;  /* Safety factor of 0.5 */
            
            if (dt_species < dt_new) {
                dt_new = dt_species;
            }
        }
    }
    
    /* Apply minimum subdivision limit */
    if (dt_new < dt / config->max_div_dt) {
        dt_new = dt / config->max_div_dt;
    }
    
    return dt_new;
}

/* ===========================================================================
 * Negative Concentration Correction - From C-RIVE procedure_exceptions.c
 * ===========================================================================*/

int rk4_check_negative(RK4State *state, int k) {
    if (!state) return 0;
    
    int corrections = 0;
    
    for (int s = 0; s < state->num_species; s++) {
        if (state->newC[s] < 0.0) {
            /* Simple correction: set to small positive value */
            state->newC[s] = RK4_EPS;
            corrections++;
            
            /* Adjust derivative to prevent further negativity */
            if (k >= 0 && k < RK4_STAGES) {
                state->mb[s].dC[k] = 0.0;
            }
        }
    }
    
    return corrections;
}

/* ===========================================================================
 * Stoichiometric Constraints
 * ===========================================================================*/

void rk4_apply_stoichiometry(RK4State *state, int k) {
    if (!state || k < 0 || k >= RK4_STAGES) return;
    
    /* Example: Ensure O2 consumption doesn't exceed available O2
     * This would be species-specific and depends on model setup
     * 
     * For now, just ensure non-negative values
     */
    for (int s = 0; s < state->num_species; s++) {
        if (state->newC[s] < 0.0) {
            state->newC[s] = 0.0;
        }
    }
}

/* ===========================================================================
 * Complete RK4 Integration
 * ===========================================================================*/

int rk4_integrate(RK4State *state, const RK4Config *config, double dt,
                  RK4RateFunc rate_func, void *user_data) {
    if (!state || !config || !rate_func) return -1;
    
    /* 1. Store initial concentrations */
    for (int s = 0; s < state->num_species; s++) {
        state->C_old[s] = state->C[s];
    }
    
    /* 2. Reset mass balance accumulators */
    rk4_reset_mb(state);
    
    /* 3. Perform RK4 stages */
    for (int k = 0; k < config->kmax; k++) {
        /* a. Update newC from previous stage */
        rk4_update_newC(state, config, dt, k);
        
        /* b. Calculate rates at this stage */
        int ret = rate_func(state, k, user_data);
        if (ret != 0) {
            fprintf(stderr, "RK4: Rate function failed at stage %d\n", k);
            return -1;
        }
        
        /* c. Check for negative concentrations */
        rk4_check_negative(state, k);
    }
    
    /* 4. Accumulate final concentrations */
    rk4_finalize(state, config, dt);
    
    /* 5. Final non-negativity check */
    for (int s = 0; s < state->num_species; s++) {
        if (state->C[s] < 0.0) {
            state->C[s] = 0.0;
        }
    }
    
    return 0;
}

/* ===========================================================================
 * Branch Data Transfer
 * ===========================================================================*/

void rk4_load_from_branch(RK4State *state, double **conc, int idx, int num_species) {
    if (!state || !conc || idx < 0) return;
    
    for (int s = 0; s < num_species && s < state->num_species; s++) {
        if (conc[s]) {
            state->C[s] = conc[s][idx];
            state->newC[s] = conc[s][idx];
            state->C_old[s] = conc[s][idx];
        }
    }
}

void rk4_save_to_branch(const RK4State *state, double **conc, int idx, int num_species) {
    if (!state || !conc || idx < 0) return;
    
    for (int s = 0; s < num_species && s < state->num_species; s++) {
        if (conc[s]) {
            conc[s][idx] = state->C[s];
        }
    }
}
