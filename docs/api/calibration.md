# Calibration API

## Data Structures

### CalibrationEngine

```c
typedef struct {
    CalibParam *params;
    int num_params;
    int active_params;          // Parameters for current stage
    
    CalibObjective *objectives;
    int num_objectives;
    
    // Seasonal objectives
    CalibObjective *seasonal_objs;
    int num_seasonal_objs;
    
    // State
    int current_stage;
    int max_iterations;
    double tolerance;
    
    // Results
    double best_rmse;
    double *best_params;
    int num_evaluations;
} CalibrationEngine;
```

### CalibParam

```c
typedef struct {
    char name[64];
    
    // Target
    TargetType target_type;     // GLOBAL or BRANCH
    char target_id[64];         // Branch name or "-"
    
    // Variable
    VarType var_type;           // CHEZY, VDB_COEF, WS, etc.
    
    // Bounds
    double min_val;
    double max_val;
    double initial;
    double current;
    double best;
    
    // Stage
    int stage;
} CalibParam;
```

### CalibObjective

```c
typedef struct {
    char branch_name[64];
    int branch_idx;
    
    // Variable
    char variable[32];
    int species_idx;
    
    // Type
    StatisticType stat_type;    // MEAN, MIN, MAX, INTRUSION_LEN, ETM_LOC
    double target_value;
    double weight;
    double threshold;           // For intrusion length
    
    // Seasonal (NEW)
    double location_km;         // Location for point observation
    double *obs_time;           // Observation times [days]
    double *obs_value;          // Observation values
    int num_obs;                // Number of observations
    int use_seasonal;           // 1 = use time-series RMSE
} CalibObjective;
```

## Main Functions

### calib_init

Initialize calibration engine.

```c
CalibrationEngine *calib_init(const char *case_dir);
```

### calib_load_params

Load parameters from CSV.

```c
int calib_load_params(CalibrationEngine *engine, const char *filepath);
```

### calib_load_objectives

Load objectives from CSV.

```c
int calib_load_objectives(CalibrationEngine *engine, const char *filepath);
```

### calib_load_seasonal_targets

Load seasonal observations.

```c
int calib_load_seasonal_targets(CalibrationEngine *engine, const char *filepath);
```

### calib_run_optimization

Run NLopt optimization.

```c
double calib_run_optimization(
    CalibrationEngine *engine,
    Network *network,
    CaseConfig *config
);
```

**Returns:** Best RMSE achieved

### calib_calculate_objective

Calculate total objective function.

```c
double calib_calculate_objective(
    CalibrationEngine *engine,
    Network *network
);
```

### calib_calc_seasonal_rmse

Calculate RMSE for seasonal objective.

```c
double calib_calc_seasonal_rmse(
    CalibrationEngine *engine,
    CalibObjective *obj,
    Network *network
);
```

### calib_apply_params_to_network

Apply current parameters to network.

```c
void calib_apply_params_to_network(
    CalibrationEngine *engine,
    Network *network
);
```

### calib_save_results

Save calibration results.

```c
int calib_save_results(
    CalibrationEngine *engine,
    const char *output_dir
);
```

### calib_free

Free calibration engine.

```c
void calib_free(CalibrationEngine *engine);
```

## NLopt Integration

### nlopt_objective_wrapper

Wrapper for NLopt objective function.

```c
double nlopt_objective_wrapper(
    unsigned n,
    const double *x,
    double *grad,
    void *data
);
```

### Supported Algorithms

| ID | NLopt Algorithm | Description |
|----|-----------------|-------------|
| 0 | `NLOPT_LN_NELDERMEAD` | Nelder-Mead simplex |
| 1 | `NLOPT_LN_BOBYQA` | BOBYQA (default) |
| 2 | `NLOPT_GN_DIRECT` | DIRECT global |
| 3 | `NLOPT_LN_COBYLA` | COBYLA constrained |
| 4 | `NLOPT_GN_CRS2_LM` | CRS2 random search |
| 5 | `NLOPT_LN_SBPLX` | Subplex |

## Example Usage

```c
// Initialize
CalibrationEngine *engine = calib_init(case_dir);

// Load configuration
calib_load_params(engine, "calibration_params.csv");
calib_load_objectives(engine, "calibration_targets.csv");
calib_load_seasonal_targets(engine, "seasonal_targets.csv");

// Set stage
engine->current_stage = 1;

// Run optimization
double best_rmse = calib_run_optimization(engine, network, config);

// Save results
calib_save_results(engine, output_dir);

// Cleanup
calib_free(engine);
```
