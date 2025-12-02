/**
 * @file calibration.h
 * @brief Generic Data-Driven Calibration Module for C-GEM Network
 * 
 * This module provides a flexible, config-file-based calibration framework
 * that allows parameter tuning without recompiling. It supports:
 * 
 * 1. MULTI-STAGE CALIBRATION:
 *    - Stage 1: Hydrodynamics (tidal range, salinity intrusion)
 *    - Stage 2: Sediment transport (SPM, ETM location)
 *    - Stage 3: Biogeochemistry (O2, nutrients, carbon)
 * 
 * 2. PARAMETER TARGETING:
 *    - GLOBAL: Apply to entire network
 *    - GROUP:  Apply to branch group (e.g., all estuaries)
 *    - BRANCH: Apply to specific branch
 * 
 * 3. OBJECTIVE FUNCTIONS:
 *    - Salinity intrusion length (4 PSU isohaline)
 *    - Mean/min/max of state variables
 *    - Tidal range at gauging stations
 *    - RMSE against observation time series
 * 
 * Configuration files (CSV format):
 *   - calibration_params.csv: Parameters to tune
 *   - calibration_targets.csv: Observations to match
 * 
 * Optimization: Uses NLopt library (optional) or simple grid search
 * 
 * Reference: 
 * - Savenije (2005) - Estuarine salt intrusion
 * - Wang et al. (2018) - C-RIVE sensitivity analysis
 * 
 * @author Nguyen Truong An
 * @date December 2024
 */

#ifndef CGEM_CALIBRATION_H
#define CGEM_CALIBRATION_H

#include "../network.h"
#include "../define.h"

/* NLopt library for optimization */
#include <nlopt.h>

/* ============================================================================
 * CONFIGURATION CONSTANTS
 * ============================================================================*/

#define CALIB_MAX_PARAMS        64      /* Maximum tunable parameters */
#define CALIB_MAX_OBJECTIVES    32      /* Maximum calibration targets */
#define CALIB_MAX_NAME_LEN      64      /* Parameter/branch name length */
#define CALIB_MAX_OBS_POINTS    1000    /* Max observation time points */
#define CALIB_SALINITY_THRESH   4.0     /* Default intrusion threshold [PSU] */

/* ============================================================================
 * ENUMERATIONS
 * ============================================================================*/

/**
 * @brief Target scope for parameter application
 */
typedef enum {
    TARGET_GLOBAL = 0,      /* Apply to entire network (BiogeoParams) */
    TARGET_GROUP = 1,       /* Apply to all branches in a group */
    TARGET_BRANCH = 2       /* Apply to specific branch */
} CalibTargetType;

/**
 * @brief Variable types that can be calibrated
 * 
 * Organized by calibration stage:
 * - HYDRO_*: Stage 1 (Hydrodynamics)
 * - SED_*:   Stage 2 (Sediment)
 * - BIO_*:   Stage 3 (Biogeochemistry)
 */
typedef enum {
    /* Stage 1: Hydrodynamics */
    VAR_CHEZY = 0,          /* Chezy friction coefficient [m^0.5/s] */
    VAR_LC_CONV,            /* Convergence length [m] */
    VAR_VDB_COEF,           /* Van den Burgh dispersion coefficient [-] */
    VAR_D0,                 /* Dispersion at mouth [m²/s] */
    VAR_MIXING_ALPHA,       /* Fischer mixing efficiency α for D0 calculation [-] */
    VAR_STORAGE_RATIO,      /* Storage width ratio RS [-] */
    VAR_RS_CHANNEL,         /* In-channel storage ratio */
    VAR_RS_FLOODPLAIN,      /* Floodplain storage ratio */
    VAR_MANNING_N,          /* Manning friction coefficient [s/m^1/3] */
    
    /* Stage 2: Sediment Transport */
    VAR_WS,                 /* Settling velocity [m/s] */
    VAR_FLOC_SAL_SCALE,     /* Flocculation salinity scale [PSU] */
    VAR_FLOC_FACTOR_MAX,    /* Maximum flocculation factor [-] */
    VAR_TAU_ERO,            /* Erosion threshold stress [Pa] */
    VAR_TAU_DEP,            /* Deposition threshold stress [Pa] */
    VAR_MERO,               /* Erosion coefficient [kg/m²/s] */
    VAR_POC_SPM_RATIO,      /* POC/SPM coupling ratio [-] */
    
    /* Stage 3: Biogeochemistry */
    VAR_KOX,                /* TOC oxidation rate [1/day] */
    VAR_KNIT,               /* Nitrification rate [1/day] */
    VAR_KDENIT,             /* Denitrification rate [1/day] */
    VAR_PBMAX1,             /* Phy1 max photosynthesis [1/day] */
    VAR_PBMAX2,             /* Phy2 max photosynthesis [1/day] */
    VAR_KMORT1,             /* Phy1 mortality rate [1/day] */
    VAR_KMORT2,             /* Phy2 mortality rate [1/day] */
    VAR_WIND_COEFF,         /* Wind-driven gas exchange coefficient */
    VAR_CURRENT_K,          /* Current reaeration factor */
    VAR_BENTHIC_RESP,       /* Benthic respiration rate [µmol/m²/day] */
    VAR_KPADS,              /* P-adsorption constant [L/µmol] */
    VAR_PAC,                /* P-adsorption capacity [µmol/mg] */
    
    /* Tropical-specific */
    VAR_SAL_STRESS_THRESH,  /* Salinity stress threshold [PSU] */
    VAR_SAL_STRESS_COEF,    /* Salinity stress coefficient */
    
    VAR_COUNT               /* Total number of variable types */
} CalibVarType;

/**
 * @brief Objective function types (what we compare to reality)
 */
typedef enum {
    OBJ_INTRUSION_LEN = 0,  /* Salt intrusion length (km to threshold) */
    OBJ_TIDAL_RANGE,        /* Tidal range at location [m] */
    OBJ_MEAN,               /* Mean value over domain/time */
    OBJ_MIN,                /* Minimum value */
    OBJ_MAX,                /* Maximum value */
    OBJ_TIMESERIES_RMSE,    /* RMSE against observation time series */
    OBJ_SPATIAL_RMSE,       /* RMSE against spatial profile */
    OBJ_ETM_LOCATION,       /* Location of Estuarine Turbidity Maximum [km] */
    OBJ_GRADIENT,           /* Gradient (e.g., dS/dx at location) */
    OBJ_COUNT
} CalibObjType;

/**
 * @brief Calibration stage
 */
typedef enum {
    STAGE_HYDRO = 1,        /* Stage 1: Hydrodynamics */
    STAGE_SEDIMENT = 2,     /* Stage 2: Sediment transport */
    STAGE_BIOGEOCHEM = 3    /* Stage 3: Biogeochemistry */
} CalibStage;

/* ============================================================================
 * DATA STRUCTURES
 * ============================================================================*/

/**
 * @brief Single calibration parameter definition
 * 
 * Loaded from calibration_params.csv:
 * Name, TargetType, TargetID, VarType, Min, Max, Initial, Stage
 */
typedef struct {
    char name[CALIB_MAX_NAME_LEN];  /* User-defined parameter name */
    CalibTargetType target_type;    /* GLOBAL, GROUP, or BRANCH */
    int target_id;                  /* Group ID or Branch ID (0 for global) */
    char target_name[CALIB_MAX_NAME_LEN]; /* Branch/group name (alternative to ID) */
    CalibVarType var_type;          /* Which variable this controls */
    
    double min_val;                 /* Lower bound for optimization */
    double max_val;                 /* Upper bound for optimization */
    double initial_val;             /* Starting value */
    double current_val;             /* Current value during optimization */
    double best_val;                /* Best value found so far */
    
    CalibStage stage;               /* Which calibration stage */
    int enabled;                    /* 1 = include in optimization */
} CalibParam;

/**
 * @brief Single calibration objective (target to match)
 * 
 * Loaded from calibration_targets.csv:
 * BranchName, Variable, Statistic, TargetValue, Weight, Tolerance
 */
typedef struct {
    char branch_name[CALIB_MAX_NAME_LEN]; /* Branch to evaluate */
    int branch_idx;                 /* Resolved branch index (-1 = all) */
    
    int species_idx;                /* Species index (-1 for hydro variables) */
    char species_name[32];          /* Species name (e.g., "SALINITY", "O2") */
    
    CalibObjType obj_type;          /* Type of objective function */
    double target_value;            /* Reality we want to match */
    double weight;                  /* Relative importance (higher = more important) */
    double tolerance;               /* Acceptable deviation (optional) */
    
    /* For intrusion length: salinity threshold */
    double threshold;               /* e.g., 4 PSU for salt intrusion */
    
    /* For spatial location (km from mouth) */
    double location_km;             /* Location for point observation [km] */
    
    /* For time series comparison (SEASONAL CALIBRATION) */
    double *obs_time;               /* Observation times [days from start] */
    double *obs_value;              /* Observation values */
    int num_obs;                    /* Number of observations */
    int use_seasonal;               /* 1 = use time-series RMSE, 0 = use mean */
    
    CalibStage stage;               /* Which calibration stage */
    int enabled;                    /* 1 = include in objective */
    
    /* Computed during calibration */
    double model_value;             /* Last computed model value */
    double residual;                /* (model - target) */
    double contribution;            /* Weighted contribution to total RMSE */
} CalibObjective;

/**
 * @brief Main calibration engine state
 */
typedef struct {
    /* Configuration */
    CalibParam *params;             /* Array of parameters to tune */
    int num_params;                 /* Number of parameters */
    
    CalibObjective *objectives;     /* Array of objectives */
    int num_objectives;             /* Number of objectives */
    
    /* Current calibration stage */
    CalibStage current_stage;       /* 1, 2, or 3 */
    
    /* Optimization settings */
    int max_iterations;             /* Max optimizer iterations */
    double tolerance;               /* Convergence tolerance */
    int algorithm;                  /* Optimization algorithm:
                                       0 = NLOPT_LN_NELDERMEAD (derivative-free simplex)
                                       1 = NLOPT_LN_BOBYQA (recommended, derivative-free)
                                       2 = NLOPT_GN_DIRECT (global optimization)
                                       3 = NLOPT_LN_COBYLA (constrained optimization) */
    double rel_tol;                 /* Relative tolerance for NLopt */
    double abs_tol;                 /* Absolute tolerance for NLopt */
    
    /* Results */
    double best_rmse;               /* Best RMSE achieved */
    int num_evaluations;            /* Total model evaluations */
    int converged;                  /* 1 if optimization converged */
    
    /* Network reference (set before optimization) */
    Network *network;               /* Pointer to network being calibrated */
    CaseConfig *config;             /* Case configuration */
    
    /* Output */
    char output_dir[CGEM_MAX_PATH]; /* Directory for calibration output */
    FILE *log_file;                 /* Log file handle */
    int verbose;                    /* Verbosity level (0-3) */
} CalibrationEngine;

/* ============================================================================
 * FUNCTION DECLARATIONS
 * ============================================================================*/

/* --- Initialization and Configuration --- */

/**
 * @brief Create and initialize a calibration engine
 * @return Pointer to new engine, or NULL on failure
 */
CalibrationEngine* calib_create_engine(void);

/**
 * @brief Free calibration engine and all resources
 * @param engine Engine to free
 */
void calib_free_engine(CalibrationEngine *engine);

/**
 * @brief Load calibration parameters from CSV file
 * @param engine Calibration engine
 * @param filepath Path to calibration_params.csv
 * @return 0 on success, -1 on failure
 */
int calib_load_params(CalibrationEngine *engine, const char *filepath);

/**
 * @brief Load calibration objectives from CSV file
 * @param engine Calibration engine
 * @param filepath Path to calibration_targets.csv
 * @return 0 on success, -1 on failure
 */
int calib_load_objectives(CalibrationEngine *engine, const char *filepath);

/**
 * @brief Load observation time series for an objective
 * @param obj Objective to update
 * @param filepath Path to observation CSV (time, value columns)
 * @return 0 on success, -1 on failure
 */
int calib_load_observations(CalibObjective *obj, const char *filepath);

/**
 * @brief Load seasonal targets from CSV file
 * 
 * Seasonal format: Branch, Variable, Location_km, Time_Day, Value
 * This allows calibration against both dry season and wet season observations.
 * 
 * @param engine Calibration engine
 * @param filepath Path to seasonal_targets.csv
 * @return 0 on success, -1 on failure
 */
int calib_load_seasonal_targets(CalibrationEngine *engine, const char *filepath);

/**
 * @brief Calculate seasonal RMSE for time-series objectives
 * 
 * Interpolates model output at observation times and computes RMSE.
 * 
 * @param engine Calibration engine
 * @param obj Objective with seasonal data
 * @param network Network with simulation results
 * @return Seasonal RMSE for this objective
 */
double calib_calc_seasonal_rmse(CalibrationEngine *engine, CalibObjective *obj, 
                                 Network *network);

/* --- Parameter Application --- */

/**
 * @brief Apply current parameter values to network
 * 
 * Iterates through enabled parameters and updates corresponding
 * Network/Branch/BiogeoParams values based on target type.
 * 
 * @param engine Calibration engine with parameters
 * @param network Network to update
 * @return 0 on success, -1 on failure
 */
int calib_apply_params_to_network(CalibrationEngine *engine, Network *network);

/**
 * @brief Apply parameters for a specific stage only
 * @param engine Calibration engine
 * @param network Network to update
 * @param stage Stage (1=hydro, 2=sed, 3=biogeo)
 * @return 0 on success
 */
int calib_apply_stage_params(CalibrationEngine *engine, Network *network, CalibStage stage);

/* --- Objective Calculation --- */

/**
 * @brief Calculate total objective function (RMSE)
 * 
 * Runs the model and computes weighted RMSE against all enabled objectives.
 * 
 * @param engine Calibration engine
 * @param network Network with current parameters
 * @param config Case configuration
 * @return Weighted RMSE (lower is better)
 */
double calib_calculate_objective(CalibrationEngine *engine, Network *network, 
                                  CaseConfig *config);

/**
 * @brief Calculate objective for current stage only
 * @param engine Calibration engine
 * @param network Network
 * @param config Configuration
 * @param stage Stage to evaluate
 * @return Stage RMSE
 */
double calib_calculate_stage_objective(CalibrationEngine *engine, Network *network,
                                       CaseConfig *config, CalibStage stage);

/**
 * @brief Calculate salt intrusion length for a branch
 * 
 * Finds the distance from ocean where salinity drops below threshold.
 * Uses linear interpolation between grid points.
 * 
 * @param branch Branch to analyze
 * @param threshold Salinity threshold [PSU] (default 4.0)
 * @return Intrusion length [km], or -1 if threshold not crossed
 */
double calib_calc_intrusion_length(Branch *branch, double threshold);

/**
 * @brief Calculate tidal range at a location
 * @param branch Branch
 * @param km_location Distance from mouth [km]
 * @return Tidal range [m] (max - min water level)
 */
double calib_calc_tidal_range(Branch *branch, double km_location);

/**
 * @brief Calculate ETM (Estuarine Turbidity Maximum) location
 * @param branch Branch with SPM data
 * @return Distance from mouth [km] where SPM is maximum
 */
double calib_calc_etm_location(Branch *branch);

/* --- Optimization --- */

/**
 * @brief Run optimization for current stage
 * 
 * Uses NLopt (if available) or simple grid search to find optimal parameters.
 * 
 * @param engine Calibration engine
 * @param network Network to calibrate
 * @param config Case configuration
 * @return Final RMSE, or -1 on failure
 */
double calib_run_optimization(CalibrationEngine *engine, Network *network,
                              CaseConfig *config);

/**
 * @brief Run multi-stage calibration workflow
 * 
 * Executes the three-stage calibration:
 * 1. Hydrodynamics (chezy, LC, RS, vdb) -> match tidal range, salinity intrusion
 * 2. Sediment (ws, floc, tau) -> match SPM, ETM location
 * 3. Biogeochemistry (kox, knit, etc.) -> match O2, nutrients
 * 
 * @param engine Calibration engine
 * @param network Network to calibrate
 * @param config Case configuration
 * @return 0 on success, -1 on failure
 */
int calib_run_multistage(CalibrationEngine *engine, Network *network,
                         CaseConfig *config);

/* --- Utility Functions --- */

/**
 * @brief Parse variable type from string
 * @param str String like "CHEZY", "KOX", etc.
 * @return CalibVarType enum value, or -1 if unknown
 */
CalibVarType calib_parse_var_type(const char *str);

/**
 * @brief Parse objective type from string
 * @param str String like "INTRUSION_LEN", "MEAN", etc.
 * @return CalibObjType enum value, or -1 if unknown
 */
CalibObjType calib_parse_obj_type(const char *str);

/**
 * @brief Parse target type from string
 * @param str String like "GLOBAL", "GROUP", "BRANCH"
 * @return CalibTargetType enum value
 */
CalibTargetType calib_parse_target_type(const char *str);

/**
 * @brief Get variable name string
 * @param var Variable type
 * @return String representation
 */
const char* calib_var_type_name(CalibVarType var);

/**
 * @brief Find branch index by name
 * @param network Network
 * @param name Branch name
 * @return Branch index, or -1 if not found
 */
int calib_find_branch(Network *network, const char *name);

/**
 * @brief Find species index by name
 * @param name Species name (e.g., "SALINITY", "O2")
 * @return Species index, or -1 for hydro variables
 */
int calib_find_species(const char *name);

/* --- Output and Reporting --- */

/**
 * @brief Write calibration results to file
 * @param engine Calibration engine
 * @param filepath Output file path
 * @return 0 on success
 */
int calib_write_results(CalibrationEngine *engine, const char *filepath);

/**
 * @brief Print calibration summary to console
 * @param engine Calibration engine
 */
void calib_print_summary(CalibrationEngine *engine);

/**
 * @brief Log a calibration iteration
 * @param engine Engine with log file
 * @param iteration Iteration number
 * @param rmse Current RMSE
 */
void calib_log_iteration(CalibrationEngine *engine, int iteration, double rmse);

#endif /* CGEM_CALIBRATION_H */
