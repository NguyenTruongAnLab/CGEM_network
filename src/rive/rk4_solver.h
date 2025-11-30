/**
 * @file rk4_solver.h
 * @brief C-RIVE RK4 Adaptive Solver for CGEM Biogeochemistry
 * 
 * Implements the 4th-order Runge-Kutta solver from C-RIVE for biogeochemical
 * reactions. Includes adaptive time-stepping for numerical stability.
 * 
 * Ported from: C-RIVE struct_RIVE.h (s_numerical), manage_simulation.c
 * 
 * The RK4 scheme:
 *   k1 = f(t, y)
 *   k2 = f(t + dt/2, y + dt*k1/2)
 *   k3 = f(t + dt/2, y + dt*k2/2)
 *   k4 = f(t + dt, y + dt*k3)
 *   y_new = y + dt*(k1 + 2*k2 + 2*k3 + k4)/6
 * 
 * References:
 * - Wang et al. (2018) Water Research 144, 341-355
 * - C-RIVE struct_RIVE.h, calc_new_C_comp()
 * 
 * CITATION:
 * Wang, S., Flipo, N., Romary, T., 2018. Time-dependent global sensitivity analysis
 * of the C-RIVE biogeochemical model. Water Research 144, 341-355.
 */

#ifndef CGEM_RK4_SOLVER_H
#define CGEM_RK4_SOLVER_H

#include "../define.h"

/* ===========================================================================
 * RK4 Solver Constants (from C-RIVE struct_RIVE.h)
 * ===========================================================================*/

#define RK4_STAGES 4           /* Number of RK4 stages */
#define RK4_MAX_DIV_DT 100.0   /* Maximum time step subdivision factor */
#define RK4_EPS 1e-10          /* Small number for numerical stability */

/* RK4 Butcher tableau coefficients
 * From C-RIVE: dt_RK[4] and coef_RK[4]
 */
static const double RK4_DT_FRAC[RK4_STAGES] = {0.0, 0.5, 0.5, 1.0};
static const double RK4_COEF[RK4_STAGES] = {1.0/6.0, 2.0/6.0, 2.0/6.0, 1.0/6.0};

/* ===========================================================================
 * RK4 Solver Configuration
 * ===========================================================================*/

/**
 * RK4 solver configuration (from C-RIVE s_numerical)
 */
typedef struct {
    char name[32];              /* Solver name ("Runge-Kutta 4") */
    int kmax;                   /* Max stages (4 for RK4, 1 for Euler) */
    double max_div_dt;          /* Max time step subdivision */
    double dt_RK[RK4_STAGES];   /* Fractional time steps */
    double coef_RK[RK4_STAGES]; /* Stage coefficients */
} RK4Config;

/**
 * Mass balance structure for a single species (from C-RIVE s_mb_species)
 * Tracks changes from each biogeochemical process
 */
typedef struct {
    /* RK4 stage derivatives [4 stages][NVAR_MB processes] */
    double dmb[RK4_STAGES][CGEM_NUM_REACTIONS];
    
    /* RK4 stage concentration changes [4 stages] */
    double dC[RK4_STAGES];
    
    /* Total concentration change over timestep */
    double deltaC;
    
    /* Final mass balance terms */
    double deltamb[CGEM_NUM_REACTIONS];
    
    /* Specific process rates (for diagnostics) */
    double phot;        /* Photosynthesis */
    double resp;        /* Respiration */
    double mort;        /* Mortality */
    double growth;      /* Growth */
    double yield;       /* Growth yield */
} RK4MassBalance;

/**
 * RK4 solver state for a single grid cell
 */
typedef struct {
    int num_species;
    double *C;              /* Current concentrations [num_species] */
    double *newC;           /* Intermediate concentrations [num_species] */
    double *C_old;          /* Concentrations at start of step [num_species] */
    RK4MassBalance *mb;     /* Mass balance per species [num_species] */
} RK4State;

/* ===========================================================================
 * Function Prototypes
 * ===========================================================================*/

/**
 * Initialize RK4 configuration with default values
 * @param config Pointer to configuration structure
 */
void rk4_init_config(RK4Config *config);

/**
 * Allocate RK4 state for a given number of species
 * @param num_species Number of species to track
 * @return Pointer to allocated state, or NULL on failure
 */
RK4State *rk4_alloc_state(int num_species);

/**
 * Free RK4 state memory
 * @param state Pointer to state to free
 */
void rk4_free_state(RK4State *state);

/**
 * Reset RK4 mass balance accumulators to zero
 * @param state Pointer to solver state
 */
void rk4_reset_mb(RK4State *state);

/**
 * Perform one RK4 stage calculation
 * Updates newC based on current stage derivatives
 * 
 * From C-RIVE calc_new_C_comp():
 *   newC = C_old + dt * dt_RK[k] * dC[k-1]   (for k > 0)
 *   newC = C_old                              (for k = 0)
 * 
 * @param state Solver state
 * @param config Solver configuration
 * @param dt Full time step [s]
 * @param k RK4 stage (0-3)
 */
void rk4_update_newC(RK4State *state, const RK4Config *config, double dt, int k);

/**
 * Accumulate final concentrations after all RK4 stages
 * 
 * From C-RIVE calc_final_conc():
 *   deltaC = dt * sum_k(coef_RK[k] * dC[k])
 *   C_new = C_old + deltaC
 * 
 * @param state Solver state
 * @param config Solver configuration
 * @param dt Time step [s]
 */
void rk4_finalize(RK4State *state, const RK4Config *config, double dt);

/**
 * Calculate adaptive time step based on species rates
 * Prevents instabilities from fast reactions
 * 
 * From C-RIVE calc_new_dt():
 *   dt_new = min(dt, C / |dC| * factor)
 * 
 * @param state Solver state
 * @param dt Requested time step [s]
 * @param config Solver configuration
 * @return Adjusted time step [s]
 */
double rk4_calc_adaptive_dt(const RK4State *state, double dt, const RK4Config *config);

/**
 * Check and correct negative concentrations
 * From C-RIVE procedure_exceptions.c
 * 
 * @param state Solver state
 * @param k Current RK4 stage
 * @return Number of corrections made
 */
int rk4_check_negative(RK4State *state, int k);

/**
 * Apply stoichiometric constraints between coupled species
 * (e.g., O2 consumption balanced with organic matter oxidation)
 * 
 * @param state Solver state
 * @param k Current RK4 stage
 */
void rk4_apply_stoichiometry(RK4State *state, int k);

/* ===========================================================================
 * RK4 Integration Helper
 * ===========================================================================*/

/**
 * Callback type for rate calculation function
 * 
 * @param state Current solver state (newC contains intermediate values)
 * @param k Current RK4 stage (0-3)
 * @param user_data User-provided context (e.g., Branch pointer)
 * @return 0 on success, -1 on error
 */
typedef int (*RK4RateFunc)(RK4State *state, int k, void *user_data);

/**
 * Perform complete RK4 integration over one timestep
 * 
 * Algorithm:
 *   1. Store initial concentrations in C_old
 *   2. For k = 0 to 3:
 *      a. Update newC from previous stage (if k > 0)
 *      b. Call rate_func to calculate dC[k]
 *      c. Check for negative concentrations
 *   3. Accumulate final concentrations
 * 
 * @param state Solver state
 * @param config Solver configuration
 * @param dt Time step [s]
 * @param rate_func Function to calculate reaction rates
 * @param user_data User context passed to rate_func
 * @return 0 on success, -1 on error
 */
int rk4_integrate(RK4State *state, const RK4Config *config, double dt,
                  RK4RateFunc rate_func, void *user_data);

/* ===========================================================================
 * Utility Functions
 * ===========================================================================*/

/**
 * Copy concentrations from branch arrays to RK4 state
 * @param state RK4 state to populate
 * @param conc 2D concentration array [species][grid]
 * @param idx Grid index
 * @param num_species Number of species
 */
void rk4_load_from_branch(RK4State *state, double **conc, int idx, int num_species);

/**
 * Copy concentrations from RK4 state back to branch arrays
 * @param state RK4 state to read from
 * @param conc 2D concentration array [species][grid]
 * @param idx Grid index
 * @param num_species Number of species
 */
void rk4_save_to_branch(const RK4State *state, double **conc, int idx, int num_species);

#endif /* CGEM_RK4_SOLVER_H */
