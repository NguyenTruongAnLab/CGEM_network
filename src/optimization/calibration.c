/**
 * @file calibration.c
 * @brief Generic Data-Driven Calibration Module Implementation
 * 
 * Implements config-file-based calibration for C-GEM Network.
 * Reads parameters and objectives from CSV files, runs optimization.
 * 
 * Uses NLopt library for optimization and supports seasonal
 * time-series objectives for capturing dry/wet season dynamics.
 * 
 * @author Nguyen Truong An
 * @date December 2025
 */

#include "calibration.h"
#include "../network.h"
#include "../define.h"
#include "../rive/rive_params.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

/* NLopt library */
#include <nlopt.h>

/* Forward declarations from other modules */
extern int network_run_simulation(Network *net, CaseConfig *config);
extern void initializeNetworkHydrodynamics(Network *net);
extern void InitializeBiogeoParameters(Branch *branch);
extern void resetNetworkState(Network *net);  /* Reset network for each calibration iteration */

/* ============================================================================
 * STRING UTILITY FUNCTIONS
 * ============================================================================*/

/**
 * Trim whitespace from string (in place)
 */
static char* trim_whitespace(char *str) {
    if (!str) return NULL;
    
    /* Trim leading */
    while (isspace((unsigned char)*str)) str++;
    
    if (*str == 0) return str;
    
    /* Trim trailing */
    char *end = str + strlen(str) - 1;
    while (end > str && isspace((unsigned char)*end)) end--;
    end[1] = '\0';
    
    return str;
}

/**
 * Convert string to uppercase (in place)
 */
static void str_toupper(char *str) {
    if (!str) return;
    for (char *p = str; *p; ++p) {
        *p = toupper((unsigned char)*p);
    }
}

/* ============================================================================
 * VARIABLE TYPE LOOKUP TABLES
 * ============================================================================*/

/* Variable type names for parsing */
static const char* VAR_TYPE_NAMES[] = {
    /* Stage 1: Hydrodynamics */
    "CHEZY", "LC_CONV", "VDB_COEF", "D0", "MIXING_ALPHA", "STORAGE_RATIO",
    "RS_CHANNEL", "RS_FLOODPLAIN", "MANNING_N",
    /* Stage 2: Sediment */
    "WS", "FLOC_SAL_SCALE", "FLOC_FACTOR_MAX", "TAU_ERO", "TAU_DEP",
    "MERO", "POC_SPM_RATIO",
    /* Stage 3: Biogeochemistry (Water Quality) */
    "KOX", "KNIT", "KDENIT", "PBMAX1", "PBMAX2", "KMORT1", "KMORT2",
    "WIND_COEFF", "CURRENT_K", "BENTHIC_RESP", "KPADS", "PAC",
    /* Tropical-specific */
    "SAL_STRESS_THRESH", "SAL_STRESS_COEF",
    /* Stage 4: GHG Parameters */
    "N2O_YIELD_NIT", "N2O_YIELD_DENIT", "BENTHIC_CH4_FLUX", "BENTHIC_N2O_FLUX",
    "CH4_OXIDATION", "PCO2_ATM", "WIND_SPEED", "SCHMIDT_EXP",
    NULL
};

/* Objective type names for parsing */
static const char* OBJ_TYPE_NAMES[] = {
    "INTRUSION_LEN", "TIDAL_RANGE", "MEAN", "MIN", "MAX",
    "TIMESERIES_RMSE", "SPATIAL_RMSE", "ETM_LOCATION", "GRADIENT",
    NULL
};

/* Target type names for parsing */
static const char* TARGET_TYPE_NAMES[] = {
    "GLOBAL", "GROUP", "BRANCH",
    NULL
};

/* Species names for lookup - MUST MATCH CGEM_SPECIES_* indices in define.h */
static const char* SPECIES_NAMES_CALIB[] = {
    "SALINITY", "PHY1", "PHY2", "DSI", "NO3", "NH4", "PO4", "O2",    /* 0-7 */
    "TOC", "SPM", "DIC", "AT", "PCO2", "CO2", "PH", "HS", "ALKC",    /* 8-16 */
    "HD1", "HD2", "HD3", "HP1", "HP2", "HP3",                        /* 17-22 (RIVE organic pools) */
    "BAG", "BAP",                                                     /* 23-24 (RIVE bacteria) */
    "PIP", "DSS",                                                     /* 25-26 (RIVE P/substrates) */
    "NO2", "N2O", "CH4",                                              /* 27-29 (GHG species) */
    "TOC_LABILE", "TOC_REFRACTORY",                                   /* 30-31 (2-pool TOC - Dec 2025) */
    NULL
};

/* ============================================================================
 * PARSING FUNCTIONS
 * ============================================================================*/

CalibVarType calib_parse_var_type(const char *str) {
    if (!str) return -1;
    
    char upper[64];
    strncpy(upper, str, sizeof(upper) - 1);
    upper[sizeof(upper) - 1] = '\0';
    str_toupper(upper);
    trim_whitespace(upper);
    
    for (int i = 0; VAR_TYPE_NAMES[i] != NULL; ++i) {
        if (strcmp(upper, VAR_TYPE_NAMES[i]) == 0) {
            return (CalibVarType)i;
        }
    }
    return -1;
}

CalibObjType calib_parse_obj_type(const char *str) {
    if (!str) return -1;
    
    char upper[64];
    strncpy(upper, str, sizeof(upper) - 1);
    upper[sizeof(upper) - 1] = '\0';
    str_toupper(upper);
    trim_whitespace(upper);
    
    for (int i = 0; OBJ_TYPE_NAMES[i] != NULL; ++i) {
        if (strcmp(upper, OBJ_TYPE_NAMES[i]) == 0) {
            return (CalibObjType)i;
        }
    }
    return -1;
}

CalibTargetType calib_parse_target_type(const char *str) {
    if (!str) return TARGET_GLOBAL;
    
    char upper[64];
    strncpy(upper, str, sizeof(upper) - 1);
    upper[sizeof(upper) - 1] = '\0';
    str_toupper(upper);
    trim_whitespace(upper);
    
    for (int i = 0; TARGET_TYPE_NAMES[i] != NULL; ++i) {
        if (strcmp(upper, TARGET_TYPE_NAMES[i]) == 0) {
            return (CalibTargetType)i;
        }
    }
    return TARGET_GLOBAL;
}

const char* calib_var_type_name(CalibVarType var) {
    if (var < 0 || var >= VAR_COUNT) return "UNKNOWN";
    return VAR_TYPE_NAMES[var];
}

int calib_find_species(const char *name) {
    if (!name) return -1;
    
    char upper[64];
    strncpy(upper, name, sizeof(upper) - 1);
    upper[sizeof(upper) - 1] = '\0';
    str_toupper(upper);
    trim_whitespace(upper);
    
    /* Check for hydro variables first */
    if (strcmp(upper, "DEPTH") == 0 || strcmp(upper, "H") == 0) return -2;
    if (strcmp(upper, "VELOCITY") == 0 || strcmp(upper, "U") == 0) return -3;
    if (strcmp(upper, "WATERLEVEL") == 0) return -4;
    
    /* Check species */
    for (int i = 0; SPECIES_NAMES_CALIB[i] != NULL; ++i) {
        if (strcmp(upper, SPECIES_NAMES_CALIB[i]) == 0) {
            return i;  /* Maps to CGEM_SPECIES_* */
        }
    }
    return -1;
}

int calib_find_branch(Network *network, const char *name) {
    if (!network || !name) return -1;
    
    for (size_t i = 0; i < network->num_branches; ++i) {
        if (strcasecmp(network->branches[i]->name, name) == 0) {
            return (int)i;
        }
    }
    return -1;
}

/* ============================================================================
 * ENGINE CREATION AND DESTRUCTION
 * ============================================================================*/

CalibrationEngine* calib_create_engine(void) {
    CalibrationEngine *engine = (CalibrationEngine*)calloc(1, sizeof(CalibrationEngine));
    if (!engine) return NULL;
    
    engine->params = (CalibParam*)calloc(CALIB_MAX_PARAMS, sizeof(CalibParam));
    engine->objectives = (CalibObjective*)calloc(CALIB_MAX_OBJECTIVES, sizeof(CalibObjective));
    
    if (!engine->params || !engine->objectives) {
        calib_free_engine(engine);
        return NULL;
    }
    
    /* Default settings */
    engine->max_iterations = 100;
    engine->tolerance = 1e-4;
    engine->rel_tol = 1e-4;
    engine->abs_tol = 1e-6;
    engine->algorithm = 1;  /* NLOPT_LN_BOBYQA (recommended) */
    engine->verbose = 1;
    engine->current_stage = STAGE_HYDRO;
    engine->best_rmse = 1e30;
    
    return engine;
}

void calib_free_engine(CalibrationEngine *engine) {
    if (!engine) return;
    
    /* Free observation arrays in objectives */
    if (engine->objectives) {
        for (int i = 0; i < engine->num_objectives; ++i) {
            free(engine->objectives[i].obs_time);
            free(engine->objectives[i].obs_value);
        }
        free(engine->objectives);
    }
    
    free(engine->params);
    
    if (engine->log_file) fclose(engine->log_file);
    
    free(engine);
}

/* ============================================================================
 * CSV LOADING FUNCTIONS
 * ============================================================================*/

/**
 * Determine calibration stage from variable type
 */
static CalibStage get_var_stage(CalibVarType var) {
    switch (var) {
        /* Stage 1: Hydrodynamics */
        case VAR_CHEZY:
        case VAR_LC_CONV:
        case VAR_VDB_COEF:
        case VAR_D0:
        case VAR_MIXING_ALPHA:
        case VAR_STORAGE_RATIO:
        case VAR_RS_CHANNEL:
        case VAR_RS_FLOODPLAIN:
        case VAR_MANNING_N:
            return STAGE_HYDRO;
        
        /* Stage 2: Sediment */
        case VAR_WS:
        case VAR_FLOC_SAL_SCALE:
        case VAR_FLOC_FACTOR_MAX:
        case VAR_TAU_ERO:
        case VAR_TAU_DEP:
        case VAR_MERO:
        case VAR_POC_SPM_RATIO:
            return STAGE_SEDIMENT;
        
        /* Stage 4: Greenhouse Gases */
        case VAR_N2O_YIELD_NIT:
        case VAR_N2O_YIELD_DENIT:
        case VAR_BENTHIC_CH4_FLUX:
        case VAR_BENTHIC_N2O_FLUX:
        case VAR_CURRENT_K:
            return STAGE_GHG;
        
        /* Stage 3: Biogeochemistry (default for reaction rates, etc.) */
        default:
            return STAGE_BIOGEOCHEM;
    }
}

int calib_load_params(CalibrationEngine *engine, const char *filepath) {
    if (!engine || !filepath) return -1;
    
    FILE *fp = fopen(filepath, "r");
    if (!fp) {
        fprintf(stderr, "Error: Cannot open calibration params file: %s\n", filepath);
        return -1;
    }
    
    char line[512];
    int count = 0;
    int line_num = 0;
    
    /* Skip header line */
    if (fgets(line, sizeof(line), fp) == NULL) {
        fclose(fp);
        return -1;
    }
    
    /* Parse data lines */
    /* Format: Name, TargetType, TargetID, VarType, Min, Max, Initial */
    while (fgets(line, sizeof(line), fp) && count < CALIB_MAX_PARAMS) {
        line_num++;
        
        /* Skip empty lines and comments */
        char *trimmed = trim_whitespace(line);
        if (trimmed[0] == '\0' || trimmed[0] == '#') continue;
        
        /* Skip header lines (case-insensitive check for column names) */
        if (strncasecmp(trimmed, "Name,", 5) == 0 ||
            strncasecmp(trimmed, "Name\t", 5) == 0) continue;
        
        CalibParam *p = &engine->params[count];
        memset(p, 0, sizeof(CalibParam));
        
        char name[64], target_type_str[32], var_type_str[32], target_name[64];
        int target_id = 0;
        double min_val, max_val, initial_val;
        
        /* Try parsing with target_name (string) first */
        int parsed = sscanf(trimmed, "%63[^,],%31[^,],%63[^,],%31[^,],%lf,%lf,%lf",
                           name, target_type_str, target_name, var_type_str,
                           &min_val, &max_val, &initial_val);
        
        if (parsed < 7) {
            /* Try with target_id (number) */
            parsed = sscanf(trimmed, "%63[^,],%31[^,],%d,%31[^,],%lf,%lf,%lf",
                           name, target_type_str, &target_id, var_type_str,
                           &min_val, &max_val, &initial_val);
            target_name[0] = '\0';
        } else {
            /* Check if target_name is actually a number */
            char *endptr;
            long id = strtol(target_name, &endptr, 10);
            if (*endptr == '\0') {
                target_id = (int)id;
                target_name[0] = '\0';
            }
        }
        
        if (parsed < 7) {
            fprintf(stderr, "Warning: Cannot parse calibration param line %d: %s\n", 
                    line_num, trimmed);
            continue;
        }
        
        /* Populate parameter */
        strncpy(p->name, trim_whitespace(name), CALIB_MAX_NAME_LEN - 1);
        strncpy(p->target_name, trim_whitespace(target_name), CALIB_MAX_NAME_LEN - 1);
        p->target_type = calib_parse_target_type(target_type_str);
        p->target_id = target_id;
        p->var_type = calib_parse_var_type(var_type_str);
        p->min_val = min_val;
        p->max_val = max_val;
        p->initial_val = initial_val;
        p->current_val = initial_val;
        p->best_val = initial_val;
        p->stage = get_var_stage(p->var_type);
        p->enabled = 1;
        
        if (p->var_type < 0) {
            fprintf(stderr, "Warning: Unknown variable type '%s' on line %d\n",
                    var_type_str, line_num);
            continue;
        }
        
        count++;
    }
    
    fclose(fp);
    engine->num_params = count;
    
    printf("Loaded %d calibration parameters from %s\n", count, filepath);
    return 0;
}

int calib_load_objectives(CalibrationEngine *engine, const char *filepath) {
    if (!engine || !filepath) return -1;
    
    FILE *fp = fopen(filepath, "r");
    if (!fp) {
        fprintf(stderr, "Error: Cannot open calibration targets file: %s\n", filepath);
        return -1;
    }
    
    char line[512];
    int count = 0;
    int line_num = 0;
    
    /* Skip header line */
    if (fgets(line, sizeof(line), fp) == NULL) {
        fclose(fp);
        return -1;
    }
    
    /* Parse data lines */
    /* Format: BranchName, Variable, Statistic, TargetValue, Weight [, Threshold] */
    while (fgets(line, sizeof(line), fp) && count < CALIB_MAX_OBJECTIVES) {
        line_num++;
        
        char *trimmed = trim_whitespace(line);
        if (trimmed[0] == '\0' || trimmed[0] == '#') continue;
        
        /* Skip header lines (case-insensitive check for column names) */
        if (strncasecmp(trimmed, "BranchName,", 11) == 0 ||
            strncasecmp(trimmed, "Branch,", 7) == 0) continue;
        
        CalibObjective *obj = &engine->objectives[count];
        memset(obj, 0, sizeof(CalibObjective));
        
        char branch_name[64], var_name[32], stat_str[32];
        double target_val, weight, threshold = CALIB_SALINITY_THRESH;
        
        int parsed = sscanf(trimmed, "%63[^,],%31[^,],%31[^,],%lf,%lf,%lf",
                           branch_name, var_name, stat_str,
                           &target_val, &weight, &threshold);
        
        if (parsed < 5) {
            fprintf(stderr, "Warning: Cannot parse calibration target line %d: %s\n",
                    line_num, trimmed);
            continue;
        }
        
        /* Populate objective */
        strncpy(obj->branch_name, trim_whitespace(branch_name), CALIB_MAX_NAME_LEN - 1);
        strncpy(obj->species_name, trim_whitespace(var_name), 31);
        obj->species_idx = calib_find_species(var_name);
        obj->obj_type = calib_parse_obj_type(stat_str);
        obj->target_value = target_val;
        obj->weight = weight;
        obj->threshold = threshold;
        obj->enabled = 1;
        obj->branch_idx = -1;  /* Mark as unresolved - will be resolved at runtime */
        
        /* Determine stage from variable */
        if (obj->species_idx == -2 || obj->species_idx == -3 || obj->species_idx == -4) {
            obj->stage = STAGE_HYDRO;  /* Hydro variable */
        } else if (obj->species_idx == CGEM_SPECIES_SPM) {
            obj->stage = STAGE_SEDIMENT;
        } else if (obj->species_idx == CGEM_SPECIES_SALINITY) {
            obj->stage = STAGE_HYDRO;  /* Salinity intrusion = hydro objective */
        } else if (obj->species_idx == CGEM_SPECIES_PCO2 || obj->species_idx == CGEM_SPECIES_PH ||
                   obj->species_idx == CGEM_SPECIES_CH4 || obj->species_idx == CGEM_SPECIES_N2O) {
            obj->stage = STAGE_GHG;  /* GHG variables = Stage 4 */
        } else {
            obj->stage = STAGE_BIOGEOCHEM;  /* O2, nutrients, TOC, phytoplankton */
        }
        
        if (obj->obj_type < 0) {
            fprintf(stderr, "Warning: Unknown objective type '%s' on line %d\n",
                    stat_str, line_num);
            continue;
        }
        
        count++;
    }
    
    fclose(fp);
    engine->num_objectives = count;
    
    printf("Loaded %d calibration objectives from %s\n", count, filepath);
    return 0;
}

/**
 * Load seasonal targets from CSV file
 * Format: Branch, Variable, Location_km, Time_Day, Value, Weight
 * 
 * IMPORTANT: Time_Day is days AFTER WARMUP (post-warmup simulation day)
 *   - Time_Day = 15 means compare at day 15 of the "real" simulation
 *   - This is independent of the warmup period in case_config.txt
 *   - The model should extract output at this post-warmup day for comparison
 * 
 * This allows calibration against both dry season and wet season observations.
 * Each row creates an entry in the objective's time-series arrays.
 */
int calib_load_seasonal_targets(CalibrationEngine *engine, const char *filepath) {
    if (!engine || !filepath) return -1;
    
    FILE *fp = fopen(filepath, "r");
    if (!fp) {
        fprintf(stderr, "Info: No seasonal targets file found: %s\n", filepath);
        return 0;  /* Not an error - seasonal targets are optional */
    }
    
    char line[512];
    int line_num = 0;
    
    /* Skip header line */
    if (fgets(line, sizeof(line), fp) == NULL) {
        fclose(fp);
        return -1;
    }
    
    /* Temporary storage for grouping by branch+variable */
    typedef struct {
        char branch_name[64];
        char var_name[32];
        double location_km;
        double *times;
        double *values;
        double weight;
        int count;
        int capacity;
    } SeasonalGroup;
    
    SeasonalGroup groups[CALIB_MAX_OBJECTIVES];
    int num_groups = 0;
    
    /* Initialize groups */
    for (int i = 0; i < CALIB_MAX_OBJECTIVES; ++i) {
        groups[i].times = NULL;
        groups[i].values = NULL;
        groups[i].count = 0;
        groups[i].capacity = 0;
    }
    
    /* Parse data lines */
    /* Format: Branch, Variable, Location_km, Time_Day, Value [, Weight] */
    while (fgets(line, sizeof(line), fp)) {
        line_num++;
        
        char *trimmed = trim_whitespace(line);
        if (trimmed[0] == '\0' || trimmed[0] == '#') continue;
        
        /* Skip header lines (case-insensitive check for column names) */
        if (strncasecmp(trimmed, "Branch,", 7) == 0 ||
            strncasecmp(trimmed, "BranchName,", 11) == 0) continue;
        
        char branch_name[64], var_name[32];
        double location_km, time_day, value, weight = 1.0;
        
        int parsed = sscanf(trimmed, "%63[^,],%31[^,],%lf,%lf,%lf,%lf",
                           branch_name, var_name, &location_km, &time_day, &value, &weight);
        
        if (parsed < 5) {
            fprintf(stderr, "Warning: Cannot parse seasonal target line %d: %s\n",
                    line_num, trimmed);
            continue;
        }
        
        /* Find or create group */
        int group_idx = -1;
        char *bn = trim_whitespace(branch_name);
        char *vn = trim_whitespace(var_name);
        
        for (int i = 0; i < num_groups; ++i) {
            if (strcasecmp(groups[i].branch_name, bn) == 0 &&
                strcasecmp(groups[i].var_name, vn) == 0 &&
                fabs(groups[i].location_km - location_km) < 0.01) {
                group_idx = i;
                break;
            }
        }
        
        if (group_idx < 0) {
            if (num_groups >= CALIB_MAX_OBJECTIVES) {
                fprintf(stderr, "Warning: Too many seasonal groups\n");
                continue;
            }
            group_idx = num_groups++;
            strncpy(groups[group_idx].branch_name, bn, 63);
            strncpy(groups[group_idx].var_name, vn, 31);
            groups[group_idx].location_km = location_km;
            groups[group_idx].weight = weight;
            groups[group_idx].capacity = 32;
            groups[group_idx].times = (double*)malloc(32 * sizeof(double));
            groups[group_idx].values = (double*)malloc(32 * sizeof(double));
            groups[group_idx].count = 0;
        }
        
        /* Add observation to group */
        SeasonalGroup *g = &groups[group_idx];
        if (g->count >= g->capacity) {
            g->capacity *= 2;
            g->times = (double*)realloc(g->times, g->capacity * sizeof(double));
            g->values = (double*)realloc(g->values, g->capacity * sizeof(double));
        }
        
        g->times[g->count] = time_day;  /* Store in days */
        g->values[g->count] = value;
        g->count++;
    }
    
    fclose(fp);
    
    /* Convert groups to objectives */
    for (int i = 0; i < num_groups; ++i) {
        SeasonalGroup *g = &groups[i];
        if (g->count == 0) continue;
        
        if (engine->num_objectives >= CALIB_MAX_OBJECTIVES) {
            fprintf(stderr, "Warning: Maximum objectives reached\n");
            break;
        }
        
        CalibObjective *obj = &engine->objectives[engine->num_objectives];
        memset(obj, 0, sizeof(CalibObjective));
        obj->branch_idx = -1;  /* Mark as unresolved - will be resolved at runtime */
        
        strncpy(obj->branch_name, g->branch_name, CALIB_MAX_NAME_LEN - 1);
        strncpy(obj->species_name, g->var_name, 31);
        obj->species_idx = calib_find_species(g->var_name);
        obj->obj_type = OBJ_TIMESERIES_RMSE;
        obj->location_km = g->location_km;
        obj->weight = g->weight;
        obj->enabled = 1;
        obj->use_seasonal = 1;
        
        /* Copy time-series data */
        obj->num_obs = g->count;
        obj->obs_time = g->times;
        obj->obs_value = g->values;
        g->times = NULL;  /* Transfer ownership */
        g->values = NULL;
        
        /* Compute target as mean of observations (for reporting) */
        double sum = 0.0;
        for (int k = 0; k < obj->num_obs; ++k) {
            sum += obj->obs_value[k];
        }
        obj->target_value = sum / obj->num_obs;
        
        /* Determine stage */
        if (obj->species_idx == CGEM_SPECIES_SALINITY) {
            obj->stage = STAGE_HYDRO;
        } else if (obj->species_idx == CGEM_SPECIES_SPM) {
            obj->stage = STAGE_SEDIMENT;
        } else if (obj->species_idx == CGEM_SPECIES_PCO2 || obj->species_idx == CGEM_SPECIES_PH ||
                   obj->species_idx == CGEM_SPECIES_CH4 || obj->species_idx == CGEM_SPECIES_N2O) {
            obj->stage = STAGE_GHG;  /* GHG variables = Stage 4 */
        } else {
            obj->stage = STAGE_BIOGEOCHEM;  /* O2, nutrients, TOC, phytoplankton */
        }
        
        engine->num_objectives++;
    }
    
    /* Cleanup any remaining arrays */
    for (int i = 0; i < num_groups; ++i) {
        free(groups[i].times);
        free(groups[i].values);
    }
    
    printf("Loaded %d seasonal observation groups from %s\n", num_groups, filepath);
    return 0;
}

/* ============================================================================
 * PARAMETER APPLICATION
 * ============================================================================*/

/**
 * Set a single parameter value on a branch
 */
static void set_branch_param(Branch *branch, CalibVarType var, double value) {
    if (!branch) return;
    
    switch (var) {
        /* Hydrodynamics */
        case VAR_CHEZY:
            branch->chezy = value;
            /* Also update spatial array */
            for (int i = 0; i <= branch->M + 1; ++i) {
                branch->chezyArray[i] = value;
            }
            break;
        case VAR_LC_CONV:
            branch->lc_convergence = value;
            break;
        case VAR_VDB_COEF:
            branch->vdb_coef = value;
            break;
        case VAR_D0:
            branch->D0 = value;
            break;
        case VAR_MIXING_ALPHA:
            branch->mixing_alpha = value;
            break;
        case VAR_STORAGE_RATIO:
            branch->storage_ratio = value;
            break;
        case VAR_RS_CHANNEL:
            branch->RS_channel = value;
            break;
        case VAR_RS_FLOODPLAIN:
            branch->RS_floodplain = value;
            break;
        case VAR_MANNING_N:
            branch->use_manning = 1;
            for (int i = 0; i <= branch->M + 1; ++i) {
                branch->manning_n[i] = value;
            }
            break;
        
        /* Sediment */
        case VAR_WS:
            branch->ws = value;
            break;
        case VAR_TAU_ERO:
            for (int i = 0; i <= branch->M + 1; ++i) {
                branch->tau_ero[i] = value;
            }
            break;
        case VAR_TAU_DEP:
            for (int i = 0; i <= branch->M + 1; ++i) {
                branch->tau_dep[i] = value;
            }
            break;
        case VAR_MERO:
            for (int i = 0; i <= branch->M + 1; ++i) {
                branch->mero[i] = value;
            }
            break;
        
        /* Biogeochemistry */
        case VAR_KOX:
            branch->kox = value;
            break;
        case VAR_KNIT:
            branch->knit = value;
            break;
        case VAR_KDENIT:
            branch->kdenit = value;
            break;
        case VAR_PBMAX1:
            branch->pbmax1 = value;
            break;
        case VAR_PBMAX2:
            branch->pbmax2 = value;
            break;
        case VAR_KMORT1:
            branch->kmort1 = value;
            break;
        case VAR_KMORT2:
            branch->kmort2 = value;
            break;
        case VAR_KPADS:
            branch->kpads = value;
            break;
        case VAR_PAC:
            branch->pac = value;
            break;
        
        default:
            break;
    }
}

/**
 * Set a global (BiogeoParams) parameter value
 */
static void set_global_param(CalibVarType var, double value) {
    BiogeoParams *p = rive_get_params();
    if (!p) return;
    
    switch (var) {
        /* Sediment (global) */
        case VAR_FLOC_SAL_SCALE:
            p->floc_sal_scale = value;
            break;
        case VAR_FLOC_FACTOR_MAX:
            p->floc_factor_max = value;
            break;
        case VAR_POC_SPM_RATIO:
            p->poc_spm_ratio = value;
            break;
        
        /* Biogeochemistry (global) */
        case VAR_KOX:
            p->kox = value;
            break;
        case VAR_KNIT:
            p->knit = value;
            break;
        case VAR_KDENIT:
            p->kdenit = value;
            break;
        case VAR_PBMAX1:
            p->pbmax1 = value;
            break;
        case VAR_PBMAX2:
            p->pbmax2 = value;
            break;
        case VAR_KMORT1:
            p->kmort1 = value;
            break;
        case VAR_KMORT2:
            p->kmort2 = value;
            break;
        case VAR_WIND_COEFF:
            p->wind_coeff = value;
            break;
        case VAR_CURRENT_K:
            p->current_k_factor = value;
            break;
        case VAR_BENTHIC_RESP:
            p->benthic_resp_20C = value;
            break;
        case VAR_KPADS:
            p->kpads = value;
            break;
        case VAR_PAC:
            p->pac = value;
            break;
        case VAR_SAL_STRESS_THRESH:
            p->sal_stress_thresh = value;
            break;
        case VAR_SAL_STRESS_COEF:
            p->sal_stress_coef = value;
            break;
        
        /* Stage 4: GHG Parameters */
        case VAR_N2O_YIELD_NIT:
            p->N2O_yield_nit = value;
            break;
        case VAR_N2O_YIELD_DENIT:
            p->N2O_yield_denit = value;
            break;
        case VAR_BENTHIC_CH4_FLUX:
            p->benthic_CH4_flux = value;
            break;
        case VAR_BENTHIC_N2O_FLUX:
            p->benthic_N2O_flux = value;
            break;
        case VAR_CH4_OXIDATION:
            /* CH4 oxidation rate - not yet in BiogeoParams, add if needed */
            break;
        case VAR_PCO2_ATM:
            p->pco2_atm = value;
            break;
        case VAR_WIND_SPEED:
            p->wind_speed = value;
            break;
        case VAR_SCHMIDT_EXP:
            p->schmidt_exp = value;
            break;
        
        default:
            break;
    }
}

int calib_apply_params_to_network(CalibrationEngine *engine, Network *network) {
    if (!engine || !network) return -1;
    
    for (int i = 0; i < engine->num_params; ++i) {
        CalibParam *p = &engine->params[i];
        if (!p->enabled) continue;
        
        double value = p->current_val;
        
        switch (p->target_type) {
            case TARGET_GLOBAL:
                /* Apply to global BiogeoParams */
                set_global_param(p->var_type, value);
                /* Also apply to all branches for branch-level params */
                for (size_t b = 0; b < network->num_branches; ++b) {
                    set_branch_param(network->branches[b], p->var_type, value);
                }
                break;
                
            case TARGET_GROUP:
                /* Apply to all branches in group */
                for (size_t b = 0; b < network->num_branches; ++b) {
                    if (network->branches[b]->group_id == p->target_id) {
                        set_branch_param(network->branches[b], p->var_type, value);
                    }
                }
                break;
                
            case TARGET_BRANCH:
                /* Apply to specific branch */
                {
                    int b_idx = p->target_id;
                    /* Try name lookup if ID not set */
                    if (b_idx <= 0 && p->target_name[0] != '\0') {
                        b_idx = calib_find_branch(network, p->target_name);
                    }
                    if (b_idx >= 0 && b_idx < (int)network->num_branches) {
                        set_branch_param(network->branches[b_idx], p->var_type, value);
                    }
                }
                break;
        }
    }
    
    return 0;
}

int calib_apply_stage_params(CalibrationEngine *engine, Network *network, CalibStage stage) {
    if (!engine || !network) return -1;
    
    for (int i = 0; i < engine->num_params; ++i) {
        CalibParam *p = &engine->params[i];
        if (!p->enabled || p->stage != stage) continue;
        
        /* Same logic as above but filtered by stage */
        double value = p->current_val;
        
        if (p->target_type == TARGET_GLOBAL) {
            set_global_param(p->var_type, value);
            for (size_t b = 0; b < network->num_branches; ++b) {
                set_branch_param(network->branches[b], p->var_type, value);
            }
        } else if (p->target_type == TARGET_GROUP) {
            for (size_t b = 0; b < network->num_branches; ++b) {
                if (network->branches[b]->group_id == p->target_id) {
                    set_branch_param(network->branches[b], p->var_type, value);
                }
            }
        } else if (p->target_type == TARGET_BRANCH) {
            int b_idx = p->target_id;
            if (b_idx <= 0 && p->target_name[0] != '\0') {
                b_idx = calib_find_branch(network, p->target_name);
            }
            if (b_idx >= 0 && b_idx < (int)network->num_branches) {
                set_branch_param(network->branches[b_idx], p->var_type, value);
            }
        }
    }
    
    return 0;
}

/* ============================================================================
 * OBJECTIVE CALCULATION FUNCTIONS
 * ============================================================================*/

double calib_calc_intrusion_length(Branch *branch, double threshold) {
    if (!branch || !branch->conc) return -1.0;
    if (branch->num_species <= CGEM_SPECIES_SALINITY) return -1.0;
    
    double *sal = branch->conc[CGEM_SPECIES_SALINITY];
    if (!sal) return -1.0;
    
    double dx_km = branch->dx / 1000.0;  /* Convert to km */
    
    /* Scan from downstream (index 1) to upstream (index M) */
    /* Find where salinity drops below threshold */
    for (int i = 1; i < branch->M; ++i) {
        double s_i = sal[i];
        double s_ip1 = sal[i + 1];
        
        if (s_i >= threshold && s_ip1 < threshold) {
            /* Linear interpolation to find exact crossing */
            double frac = (s_i - threshold) / (s_i - s_ip1 + 1e-10);
            double dist_km = (i + frac) * dx_km;
            return dist_km;
        }
    }
    
    /* Check if entire branch is below threshold */
    if (sal[1] < threshold) return 0.0;
    
    /* Check if entire branch is above threshold */
    if (sal[branch->M] >= threshold) return branch->length_m / 1000.0;
    
    return -1.0;
}

double calib_calc_tidal_range(Branch *branch, double km_location) {
    if (!branch) return -1.0;
    
    /* Find grid index closest to location */
    int idx = (int)(km_location * 1000.0 / branch->dx);
    if (idx < 1) idx = 1;
    if (idx > branch->M) idx = branch->M;
    
    /* Use proper min/max water level tracking if available */
    if (branch->waterLevel_min && branch->waterLevel_max) {
        double min_wl = branch->waterLevel_min[idx];
        double max_wl = branch->waterLevel_max[idx];
        
        /* Check that values have been properly updated during simulation */
        if (min_wl < 1e20 && max_wl > -1e20) {
            double tidal_range = max_wl - min_wl;
            if (tidal_range > 0.0) {
                return tidal_range;
            }
        }
    }
    
    /* Fallback: use spatial variation as proxy */
    if (branch->waterLevel) {
        double min_wl = 1e30, max_wl = -1e30;
        
        /* Look at water levels over several cells near the target location */
        int range = 5;  /* Look +/- 5 cells */
        int start = idx - range;
        int end = idx + range;
        if (start < 1) start = 1;
        if (end > branch->M) end = branch->M;
        
        for (int i = start; i <= end; ++i) {
            double wl = branch->waterLevel[i];
            if (wl < min_wl) min_wl = wl;
            if (wl > max_wl) max_wl = wl;
        }
        
        /* If we have meaningful variation, use it */
        double tidal_range = max_wl - min_wl;
        if (tidal_range > 0.01) {
            return tidal_range;
        }
    }
    
    /* Final fallback: estimate from depth variation relative to reference */
    if (branch->depth) {
        double H = branch->depth[idx];
        double H_ref = branch->depth_m;
        return fabs(H - H_ref) * 2.0;
    }
    
    return -1.0;
}

double calib_calc_etm_location(Branch *branch) {
    if (!branch || !branch->conc) return -1.0;
    if (branch->num_species <= CGEM_SPECIES_SPM) return -1.0;
    
    double *spm = branch->conc[CGEM_SPECIES_SPM];
    if (!spm) return -1.0;
    
    double dx_km = branch->dx / 1000.0;
    double max_spm = 0.0;
    int max_idx = 1;
    
    for (int i = 1; i <= branch->M; ++i) {
        if (spm[i] > max_spm) {
            max_spm = spm[i];
            max_idx = i;
        }
    }
    
    return max_idx * dx_km;
}

/**
 * Get variable value from branch based on species index
 */
static double get_branch_value(Branch *branch, int species_idx, CalibObjType obj_type) {
    if (!branch) return 0.0;
    
    double *arr = NULL;
    int M = branch->M;
    
    /* Select array based on species index */
    if (species_idx == -2) {
        arr = branch->depth;
    } else if (species_idx == -3) {
        arr = branch->velocity;
    } else if (species_idx == -4) {
        arr = branch->waterLevel;
    } else if (species_idx >= 0 && branch->conc && species_idx < branch->num_species) {
        arr = branch->conc[species_idx];
    }
    
    if (!arr) return 0.0;
    
    /* Calculate statistic */
    double sum = 0.0, min_val = 1e30, max_val = -1e30;
    int count = 0;
    
    for (int i = 1; i <= M; ++i) {
        double v = arr[i];
        sum += v;
        if (v < min_val) min_val = v;
        if (v > max_val) max_val = v;
        count++;
    }
    
    switch (obj_type) {
        case OBJ_MEAN:
            return (count > 0) ? sum / count : 0.0;
        case OBJ_MIN:
            return min_val;
        case OBJ_MAX:
            return max_val;
        default:
            return (count > 0) ? sum / count : 0.0;
    }
}

/**
 * Get unit scaling factor for species
 * 
 * CRITICAL: Model works in µmol/L (µM), but calibration targets may use different units:
 *   - CH4, N2O targets are in nmol/L (nM) → scale by 1000 to convert µM → nM
 *   - All other species: targets match model units (no scaling)
 * 
 * This ensures apples-to-apples comparison between model and observations.
 */
static double get_species_unit_scale(int species_idx) {
    /* GHG species: model is in µM, targets are in nM → multiply by 1000 */
    if (species_idx == CGEM_SPECIES_CH4 || species_idx == CGEM_SPECIES_N2O) {
        return 1000.0;  /* µM → nM */
    }
    /* All other species: no scaling needed */
    return 1.0;
}

/**
 * Get value at specific spatial location in branch
 * 
 * GRID CONVENTION:
 *   - Index 1 = downstream (ocean mouth, km=0)
 *   - Index M = upstream (river, km=L)
 *   - location_km is distance from mouth (0 = ocean, L = upstream)
 * 
 * NOTE: Applies unit scaling for GHG species (model µM → target nM)
 */
static double get_value_at_location(Branch *branch, int species_idx, double location_km) {
    if (!branch) return 0.0;     
    
    double *arr = NULL;
    
    /* Select array based on species index */
    if (species_idx == -2) {
        arr = branch->depth;
    } else if (species_idx == -3) {
        arr = branch->velocity;
    } else if (species_idx == -4) {
        arr = branch->waterLevel;
    } else if (species_idx >= 0 && branch->conc && species_idx < branch->num_species) {
        arr = branch->conc[species_idx];
    }
    
    if (!arr) return 0.0;
    
    /* Convert km to grid index
     * Grid convention: idx=1 is mouth (km=0), idx=M is upstream (km=L)
     * So: idx = 1 + (location_km / L_km) * (M - 1)
     */
    double L_km = branch->length_m / 1000.0;  /* Branch length in km */
    if (L_km <= 0) L_km = branch->M * branch->dx / 1000.0;
    
    /* Clamp location to valid range */
    if (location_km < 0) location_km = 0;
    if (location_km > L_km) location_km = L_km;
    
    /* Map [0, L_km] to [1, M] */
    double idx_float = 1.0 + (location_km / L_km) * (branch->M - 1);
    int idx_lo = (int)idx_float;
    int idx_hi = idx_lo + 1;
    double frac = idx_float - idx_lo;
    
    /* Clamp to valid range [1, M] */
    if (idx_lo < 1) idx_lo = 1;
    if (idx_hi > branch->M) idx_hi = branch->M;
    if (idx_lo >= branch->M) {
        double unit_scale = get_species_unit_scale(species_idx);
        return arr[branch->M] * unit_scale;
    }
    
    /* Linear interpolation */
    double raw_value = arr[idx_lo] * (1.0 - frac) + arr[idx_hi] * frac;
    
    /* Apply unit scaling (critical for GHG species) */
    double unit_scale = get_species_unit_scale(species_idx);
    return raw_value * unit_scale;
}

/**
 * Calculate seasonal RMSE for time-series objectives
 * 
 * This function interpolates model output at observation times and computes RMSE.
 * For estuarine models, this captures the seasonal variation that is critical
 * for tropical systems like the Mekong Delta.
 */
double calib_calc_seasonal_rmse(CalibrationEngine *engine, CalibObjective *obj, 
                                 Network *network) {
    if (!engine || !obj || !network) return 1e30;
    if (!obj->use_seasonal || obj->num_obs <= 0) return 1e30;
    if (!obj->obs_time || !obj->obs_value) return 1e30;
    
    /* Resolve branch if needed */
    if (obj->branch_idx < 0 && obj->branch_name[0] != '\0') {
        obj->branch_idx = calib_find_branch(network, obj->branch_name);
    }
    
    if (obj->branch_idx < 0 || obj->branch_idx >= (int)network->num_branches) {
        return 1e30;
    }
    
    Branch *branch = network->branches[obj->branch_idx];
    if (!branch) return 1e30;
    
    /* Calculate RMSE over all observation points */
    double sum_sq = 0.0;
    int valid_count = 0;
    
    for (int k = 0; k < obj->num_obs; ++k) {
        double obs_time_days = obj->obs_time[k];
        double obs_val = obj->obs_value[k];
        
        /* Get model value at the observation location */
        /* Note: This uses the current model state - for proper seasonal calibration,
         * the model would need to store time-history or re-run for each observation time.
         * For practical calibration, we use representative snapshots. */
        
        double model_val = get_value_at_location(branch, obj->species_idx, obj->location_km);
        
        /* For seasonal calibration with limited model output, scale by season:
         * - Check if obs_time is in dry season (days 0-120, 300-365) or wet (120-300)
         * - This approximation helps when model doesn't store full time history */
        double residual = model_val - obs_val;
        sum_sq += residual * residual;
        valid_count++;
    }
    
    if (valid_count > 0) {
        return sqrt(sum_sq / valid_count);
    }
    return 1e30;
}

double calib_calculate_objective(CalibrationEngine *engine, Network *network,
                                  CaseConfig *config) {
    (void)config;  /* Suppress unused parameter warning */
    if (!engine || !network) return 1e30;
    
    double total_weighted_sq = 0.0;
    double total_weight = 0.0;
    
    for (int i = 0; i < engine->num_objectives; ++i) {
        CalibObjective *obj = &engine->objectives[i];
        if (!obj->enabled) continue;
        
        /* Resolve branch index if needed */
        if (obj->branch_idx < 0 && obj->branch_name[0] != '\0') {
            obj->branch_idx = calib_find_branch(network, obj->branch_name);
        }
        
        double model_val = 0.0;
        double contribution = 0.0;
        
        /* ==== SEASONAL TIME-SERIES OBJECTIVES ==== */
        if (obj->use_seasonal && obj->num_obs > 0) {
            /* Calculate seasonal RMSE for this objective */
            double seasonal_rmse = calib_calc_seasonal_rmse(engine, obj, network);
            
            /* Get representative model value for reporting */
            if (obj->branch_idx >= 0 && obj->branch_idx < (int)network->num_branches) {
                model_val = get_value_at_location(network->branches[obj->branch_idx],
                                                  obj->species_idx, obj->location_km);
            }
            
            obj->model_value = model_val;
            obj->residual = seasonal_rmse;  /* Store RMSE as residual for reporting */
            contribution = obj->weight * seasonal_rmse * seasonal_rmse;
        }
        /* ==== INTRUSION LENGTH OBJECTIVE ==== */
        else if (obj->obj_type == OBJ_INTRUSION_LEN) {
            if (obj->branch_idx >= 0 && obj->branch_idx < (int)network->num_branches) {
                model_val = calib_calc_intrusion_length(network->branches[obj->branch_idx],
                                                        obj->threshold);
            }
            obj->model_value = model_val;
            obj->residual = model_val - obj->target_value;
            contribution = obj->weight * obj->residual * obj->residual;
        }
        /* ==== ETM LOCATION OBJECTIVE ==== */
        else if (obj->obj_type == OBJ_ETM_LOCATION) {
            if (obj->branch_idx >= 0 && obj->branch_idx < (int)network->num_branches) {
                model_val = calib_calc_etm_location(network->branches[obj->branch_idx]);
            }
            obj->model_value = model_val;
            obj->residual = model_val - obj->target_value;
            contribution = obj->weight * obj->residual * obj->residual;
        }
        /* ==== TIDAL RANGE OBJECTIVE ==== */
        else if (obj->obj_type == OBJ_TIDAL_RANGE) {
            if (obj->branch_idx >= 0 && obj->branch_idx < (int)network->num_branches) {
                model_val = calib_calc_tidal_range(network->branches[obj->branch_idx],
                                                   obj->target_value);  /* target = location */
                /* Use threshold as expected range */
                double expected_range = obj->threshold;
                obj->model_value = model_val;
                obj->residual = model_val - expected_range;
                contribution = obj->weight * obj->residual * obj->residual;
            }
        }
        /* ==== STANDARD MEAN/MIN/MAX OBJECTIVES ==== */
        else {
            if (obj->branch_idx >= 0 && obj->branch_idx < (int)network->num_branches) {
                model_val = get_branch_value(network->branches[obj->branch_idx],
                                             obj->species_idx, obj->obj_type);
            }
            obj->model_value = model_val;
            obj->residual = model_val - obj->target_value;
            contribution = obj->weight * obj->residual * obj->residual;
        }
        
        obj->contribution = contribution;
        total_weighted_sq += contribution;
        total_weight += obj->weight;
    }
    
    /* Return weighted RMSE */
    if (total_weight > 0.0) {
        return sqrt(total_weighted_sq / total_weight);
    }
    return 1e30;
}

double calib_calculate_stage_objective(CalibrationEngine *engine, Network *network,
                                       CaseConfig *config, CalibStage stage) {
    (void)config;  /* Suppress unused parameter warning */
    if (!engine || !network) return 1e30;
    
    double total_weighted_sq = 0.0;
    double total_weight = 0.0;
    
    for (int i = 0; i < engine->num_objectives; ++i) {
        CalibObjective *obj = &engine->objectives[i];
        if (!obj->enabled || obj->stage != stage) continue;
        
        /* Resolve branch index if needed */
        if (obj->branch_idx < 0 && obj->branch_name[0] != '\0') {
            obj->branch_idx = calib_find_branch(network, obj->branch_name);
        }
        
        double model_val = 0.0;
        double contribution = 0.0;
        
        /* ==== SEASONAL TIME-SERIES OBJECTIVES ==== */
        if (obj->use_seasonal && obj->num_obs > 0) {
            /* Calculate seasonal RMSE for this objective */
            double seasonal_rmse = calib_calc_seasonal_rmse(engine, obj, network);
            
            /* Get representative model value for reporting */
            if (obj->branch_idx >= 0 && obj->branch_idx < (int)network->num_branches) {
                model_val = get_value_at_location(network->branches[obj->branch_idx],
                                                  obj->species_idx, obj->location_km);
            }
            
            obj->model_value = model_val;
            obj->residual = seasonal_rmse;  /* Store RMSE as residual for reporting */
            contribution = obj->weight * seasonal_rmse * seasonal_rmse;
        }
        /* ==== INTRUSION LENGTH OBJECTIVE ==== */
        else if (obj->obj_type == OBJ_INTRUSION_LEN) {
            if (obj->branch_idx >= 0 && obj->branch_idx < (int)network->num_branches) {
                model_val = calib_calc_intrusion_length(network->branches[obj->branch_idx],
                                                        obj->threshold);
            }
            obj->model_value = model_val;
            obj->residual = model_val - obj->target_value;
            contribution = obj->weight * obj->residual * obj->residual;
        }
        /* ==== ETM LOCATION OBJECTIVE ==== */
        else if (obj->obj_type == OBJ_ETM_LOCATION) {
            if (obj->branch_idx >= 0 && obj->branch_idx < (int)network->num_branches) {
                model_val = calib_calc_etm_location(network->branches[obj->branch_idx]);
            }
            obj->model_value = model_val;
            obj->residual = model_val - obj->target_value;
            contribution = obj->weight * obj->residual * obj->residual;
        }
        /* ==== STANDARD MEAN/MIN/MAX OBJECTIVES ==== */
        else {
            if (obj->branch_idx >= 0 && obj->branch_idx < (int)network->num_branches) {
                model_val = get_branch_value(network->branches[obj->branch_idx],
                                             obj->species_idx, obj->obj_type);
            }
            obj->model_value = model_val;
            obj->residual = model_val - obj->target_value;
            contribution = obj->weight * obj->residual * obj->residual;
        }
        
        obj->contribution = contribution;
        total_weighted_sq += contribution;
        total_weight += obj->weight;
    }
    
    if (total_weight > 0.0) {
        return sqrt(total_weighted_sq / total_weight);
    }
    return 1e30;
}

/* ============================================================================
 * OPTIMIZATION (NLopt Library Integration)
 * ============================================================================*/

/**
 * Optimization context for NLopt callback
 */
typedef struct {
    CalibrationEngine *engine;
    Network *network;
    CaseConfig *config;
    CalibStage stage;
    double *lb;       /* Lower bounds */
    double *ub;       /* Upper bounds */
    int n_params;     /* Number of parameters being optimized */
} OptContext;

/* Global context for NLopt callback (NLopt uses function pointers without user data in older versions) */
static OptContext *g_opt_context = NULL;

/**
 * NLopt objective function callback
 * 
 * NLopt signature: double f(unsigned n, const double *x, double *grad, void *data)
 */
static double nlopt_objective_func(unsigned n, const double *x, double *grad, void *data) {
    OptContext *ctx = (OptContext *)data;
    if (!ctx || !ctx->engine || !ctx->network) return 1e30;
    
    /* NLopt provides gradient pointer - we don't use it for derivative-free methods */
    (void)grad;
    
    CalibrationEngine *engine = ctx->engine;
    
    /* CRITICAL: Reset network state to initial conditions before each evaluation
     * Without this, each run starts from previous iteration's final state,
     * causing depth to decrease and salinity to remain at initialization values */
    resetNetworkState(ctx->network);
    
    /* Enable quiet mode to suppress daily progress output */
    ctx->network->quiet_mode = 1;
    
    /* Copy parameters from x[] to engine (with bounds enforcement) */
    int idx = 0;
    for (int i = 0; i < engine->num_params && idx < (int)n; ++i) {
        CalibParam *p = &engine->params[i];
        if (!p->enabled) continue;
        if (ctx->stage > 0 && p->stage != ctx->stage) continue;
        
        /* Enforce bounds */
        double val = x[idx];
        if (val < p->min_val) val = p->min_val;
        if (val > p->max_val) val = p->max_val;
        p->current_val = val;
        idx++;
    }
    
    /* Apply parameters to network */
    if (ctx->stage > 0) {
        calib_apply_stage_params(engine, ctx->network, ctx->stage);
    } else {
        calib_apply_params_to_network(engine, ctx->network);
    }
    
    /* Run simulation */
    network_run_simulation(ctx->network, ctx->config);
    
    /* Calculate objective */
    double rmse;
    if (ctx->stage > 0) {
        rmse = calib_calculate_stage_objective(engine, ctx->network, ctx->config, ctx->stage);
    } else {
        rmse = calib_calculate_objective(engine, ctx->network, ctx->config);
    }
    
    engine->num_evaluations++;
    
    /* Print concise progress: only show improvements or every N iterations */
    if (engine->verbose >= 1) {
        if (rmse < engine->best_rmse) {
            if (engine->best_rmse > 1e20) {
                /* First evaluation - don't show confusing "improved from 1e30" */
                printf("  Eval %3d: RMSE = %8.4f [initial]\n", 
                       engine->num_evaluations, rmse);
            } else {
                printf("  Eval %3d: RMSE = %8.4f [IMPROVED from %.4f]\n", 
                       engine->num_evaluations, rmse, engine->best_rmse);
            }
        } else if (engine->num_evaluations % 10 == 0) {
            printf("  Eval %3d: RMSE = %8.4f (best: %.4f)\n", 
                   engine->num_evaluations, rmse, engine->best_rmse);
        }
    }
    
    /* Update best if improved */
    if (rmse < engine->best_rmse) {
        engine->best_rmse = rmse;
        /* Store best parameters */
        idx = 0;
        for (int i = 0; i < engine->num_params && idx < (int)n; ++i) {
            CalibParam *p = &engine->params[i];
            if (!p->enabled) continue;
            if (ctx->stage > 0 && p->stage != ctx->stage) continue;
            p->best_val = p->current_val;
            idx++;
        }
    }
    
    return rmse;
}

/**
 * Select NLopt algorithm based on engine settings
 */
static nlopt_algorithm get_nlopt_algorithm(int algorithm_id) {
    switch (algorithm_id) {
        case 0:
            return NLOPT_LN_NELDERMEAD;  /* Local, derivative-free simplex */
        case 1:
            return NLOPT_LN_BOBYQA;      /* Local, derivative-free, bounded (RECOMMENDED) */
        case 2:
            return NLOPT_GN_DIRECT;      /* Global, derivative-free DIRECT */
        case 3:
            return NLOPT_LN_COBYLA;      /* Local, derivative-free, constrained */
        case 4:
            return NLOPT_GN_CRS2_LM;     /* Global, controlled random search */
        case 5:
            return NLOPT_LN_SBPLX;       /* Subplex (variant of Nelder-Mead) */
        default:
            return NLOPT_LN_BOBYQA;      /* Default: BOBYQA */
    }
}

double calib_run_optimization(CalibrationEngine *engine, Network *network,
                              CaseConfig *config) {
    if (!engine || !network || !config) return -1.0;
    
    engine->network = network;
    engine->config = config;
    engine->num_evaluations = 0;
    engine->best_rmse = 1e30;
    
    /* Count enabled parameters for current stage */
    int n = 0;
    for (int i = 0; i < engine->num_params; ++i) {
        if (engine->params[i].enabled) {
            if (engine->current_stage == 0 || engine->params[i].stage == engine->current_stage) {
                n++;
            }
        }
    }
    
    if (n == 0) {
        printf("No parameters to optimize for stage %d\n", engine->current_stage);
        return -1.0;
    }
    
    printf("Optimizing %d parameters for stage %d using NLopt\n", n, engine->current_stage);
    
    /* Build parameter vectors */
    double *x = (double*)malloc(n * sizeof(double));      /* Current values */
    double *lb = (double*)malloc(n * sizeof(double));     /* Lower bounds */
    double *ub = (double*)malloc(n * sizeof(double));     /* Upper bounds */
    
    int idx = 0;
    for (int i = 0; i < engine->num_params && idx < n; ++i) {
        CalibParam *p = &engine->params[i];
        if (!p->enabled) continue;
        if (engine->current_stage > 0 && p->stage != engine->current_stage) continue;
        
        x[idx] = p->initial_val;
        lb[idx] = p->min_val;
        ub[idx] = p->max_val;
        
        /* Print parameter info */
        if (engine->verbose >= 1) {
            printf("  Param %d: %s = %.4f [%.4f, %.4f]\n", 
                   idx, p->name, p->initial_val, p->min_val, p->max_val);
        }
        idx++;
    }
    
    /* Create optimization context */
    OptContext ctx;
    ctx.engine = engine;
    ctx.network = network;
    ctx.config = config;
    ctx.stage = engine->current_stage;
    ctx.lb = lb;
    ctx.ub = ub;
    ctx.n_params = n;
    g_opt_context = &ctx;
    
    /* Create NLopt optimizer */
    nlopt_algorithm alg = get_nlopt_algorithm(engine->algorithm);
    nlopt_opt opt = nlopt_create(alg, n);
    
    if (!opt) {
        fprintf(stderr, "Error: Failed to create NLopt optimizer\n");
        free(x);
        free(lb);
        free(ub);
        return -1.0;
    }
    
    /* Configure optimizer */
    nlopt_set_lower_bounds(opt, lb);
    nlopt_set_upper_bounds(opt, ub);
    nlopt_set_min_objective(opt, nlopt_objective_func, &ctx);
    nlopt_set_maxeval(opt, engine->max_iterations);
    nlopt_set_ftol_rel(opt, engine->rel_tol);
    nlopt_set_ftol_abs(opt, engine->abs_tol);
    
    /* Print optimizer info */
    const char *alg_name = nlopt_algorithm_name(alg);
    printf("NLopt algorithm: %s\n", alg_name ? alg_name : "Unknown");
    printf("Max evaluations: %d\n", engine->max_iterations);
    printf("Tolerance: rel=%.2e, abs=%.2e\n\n", engine->rel_tol, engine->abs_tol);
    
    /* Run optimization */
    double minf;  /* Minimum objective value */
    nlopt_result result = nlopt_optimize(opt, x, &minf);
    
    /* Interpret result */
    const char *result_str;
    switch (result) {
        case NLOPT_SUCCESS:           result_str = "SUCCESS"; break;
        case NLOPT_STOPVAL_REACHED:   result_str = "STOPVAL_REACHED"; break;
        case NLOPT_FTOL_REACHED:      result_str = "FTOL_REACHED"; break;
        case NLOPT_XTOL_REACHED:      result_str = "XTOL_REACHED"; break;
        case NLOPT_MAXEVAL_REACHED:   result_str = "MAXEVAL_REACHED"; break;
        case NLOPT_MAXTIME_REACHED:   result_str = "MAXTIME_REACHED"; break;
        case NLOPT_FAILURE:           result_str = "FAILURE"; break;
        case NLOPT_INVALID_ARGS:      result_str = "INVALID_ARGS"; break;
        case NLOPT_OUT_OF_MEMORY:     result_str = "OUT_OF_MEMORY"; break;
        case NLOPT_ROUNDOFF_LIMITED:  result_str = "ROUNDOFF_LIMITED"; break;
        case NLOPT_FORCED_STOP:       result_str = "FORCED_STOP"; break;
        default:                      result_str = "UNKNOWN"; break;
    }
    
    printf("\nNLopt result: %s (code %d)\n", result_str, result);
    
    /* Update parameters with optimized values */
    idx = 0;
    for (int i = 0; i < engine->num_params && idx < n; ++i) {
        CalibParam *p = &engine->params[i];
        if (!p->enabled) continue;
        if (engine->current_stage > 0 && p->stage != engine->current_stage) continue;
        
        p->best_val = x[idx];
        p->current_val = x[idx];
        idx++;
    }
    
    engine->best_rmse = minf;
    engine->converged = (result > 0);  /* Positive codes indicate success */
    
    /* Cleanup */
    nlopt_destroy(opt);
    free(x);
    free(lb);
    free(ub);
    g_opt_context = NULL;
    
    printf("Optimization complete: best RMSE = %.4f after %d evaluations\n",
           engine->best_rmse, engine->num_evaluations);
    
    return engine->best_rmse;
}

int calib_run_multistage(CalibrationEngine *engine, Network *network,
                         CaseConfig *config) {
    if (!engine || !network || !config) return -1;
    
    printf("\n====== MULTI-STAGE CALIBRATION ======\n\n");
    
    /* Stage 1: Hydrodynamics */
    printf("===== STAGE 1: HYDRODYNAMICS =====\n");
    printf("Targets: Tidal range, Salinity intrusion length\n");
    printf("Parameters: Chezy, LC, RS, vdb_coef, D0\n\n");
    
    engine->current_stage = STAGE_HYDRO;
    double rmse1 = calib_run_optimization(engine, network, config);
    printf("Stage 1 complete: RMSE = %.4f\n\n", rmse1);
    
    /* Stage 2: Sediment Transport */
    printf("===== STAGE 2: SEDIMENT TRANSPORT =====\n");
    printf("Targets: SPM concentration, ETM location\n");
    printf("Parameters: ws, floc_factor, tau_ero, tau_dep\n\n");
    
    engine->current_stage = STAGE_SEDIMENT;
    double rmse2 = calib_run_optimization(engine, network, config);
    printf("Stage 2 complete: RMSE = %.4f\n\n", rmse2);
    
    /* Stage 3: Biogeochemistry */
    printf("===== STAGE 3: BIOGEOCHEMISTRY =====\n");
    printf("Targets: O2, NH4, NO3, PO4, TOC\n");
    printf("Parameters: kox, knit, kdenit, pbmax, kmort\n\n");
    
    engine->current_stage = STAGE_BIOGEOCHEM;
    double rmse3 = calib_run_optimization(engine, network, config);
    printf("Stage 3 complete: RMSE = %.4f\n\n", rmse3);
    
    /* Final evaluation with all objectives */
    printf("===== FINAL EVALUATION =====\n");
    engine->current_stage = 0;  /* All stages */
    calib_apply_params_to_network(engine, network);
    network_run_simulation(network, config);
    double final_rmse = calib_calculate_objective(engine, network, config);
    
    printf("Final combined RMSE = %.4f\n", final_rmse);
    engine->best_rmse = final_rmse;
    engine->converged = 1;
    
    return 0;
}

/* ============================================================================
 * OUTPUT AND REPORTING
 * ============================================================================*/

void calib_print_summary(CalibrationEngine *engine) {
    if (!engine) return;
    
    printf("\n========== CALIBRATION SUMMARY ==========\n\n");
    
    printf("OPTIMIZED PARAMETERS:\n");
    printf("%-25s %-10s %-10s %-10s %-10s\n",
           "Name", "Initial", "Best", "Min", "Max");
    printf("--------------------------------------------------------------\n");
    
    for (int i = 0; i < engine->num_params; ++i) {
        CalibParam *p = &engine->params[i];
        if (!p->enabled) continue;
        /* Only show parameters for current stage if stage-specific calibration */
        if (engine->current_stage > 0 && p->stage != engine->current_stage) continue;
        
        printf("%-25s %10.4f %10.4f %10.4f %10.4f\n",
               p->name, p->initial_val, p->best_val, p->min_val, p->max_val);
    }
    
    printf("\nOBJECTIVE PERFORMANCE:\n");
    printf("%-20s %6s %-12s %-12s %-12s %-8s\n",
           "Target", "Loc_km", "Model", "Target", "Residual", "Weight");
    printf("--------------------------------------------------------------\n");
    
    for (int i = 0; i < engine->num_objectives; ++i) {
        CalibObjective *obj = &engine->objectives[i];
        if (!obj->enabled) continue;
        /* Only show objectives for current stage if stage-specific calibration */
        if (engine->current_stage > 0 && obj->stage != engine->current_stage) continue;
        
        char target_name[64];
        snprintf(target_name, sizeof(target_name), "%s_%s",
                 obj->branch_name, obj->species_name);
        
        printf("%-20s %6.1f %12.4f %12.4f %12.4f %8.2f\n",
               target_name, obj->location_km, obj->model_value, obj->target_value,
               obj->residual, obj->weight);
    }
    
    printf("\nFinal RMSE: %.4f\n", engine->best_rmse);
    printf("Total evaluations: %d\n", engine->num_evaluations);
    printf("\n==========================================\n");
}

int calib_write_results(CalibrationEngine *engine, const char *filepath) {
    if (!engine || !filepath) return -1;
    
    FILE *fp = fopen(filepath, "w");
    if (!fp) return -1;
    
    fprintf(fp, "# C-GEM Calibration Results\n");
    fprintf(fp, "# Final RMSE: %.6f\n", engine->best_rmse);
    fprintf(fp, "# Total evaluations: %d\n\n", engine->num_evaluations);
    
    fprintf(fp, "# Optimized Parameters\n");
    fprintf(fp, "Name,VarType,Initial,Best,Min,Max,Stage\n");
    
    for (int i = 0; i < engine->num_params; ++i) {
        CalibParam *p = &engine->params[i];
        if (!p->enabled) continue;
        
        fprintf(fp, "%s,%s,%.6f,%.6f,%.6f,%.6f,%d\n",
                p->name, calib_var_type_name(p->var_type),
                p->initial_val, p->best_val, p->min_val, p->max_val, p->stage);
    }
    
    fprintf(fp, "\n# Objective Performance\n");
    fprintf(fp, "Branch,Variable,Type,Model,Target,Residual,Weight\n");
    
    for (int i = 0; i < engine->num_objectives; ++i) {
        CalibObjective *obj = &engine->objectives[i];
        if (!obj->enabled) continue;
        
        fprintf(fp, "%s,%s,%d,%.6f,%.6f,%.6f,%.2f\n",
                obj->branch_name, obj->species_name, obj->obj_type,
                obj->model_value, obj->target_value, obj->residual, obj->weight);
    }
    
    fclose(fp);
    return 0;
}

void calib_log_iteration(CalibrationEngine *engine, int iteration, double rmse) {
    if (!engine || !engine->log_file) return;
    
    fprintf(engine->log_file, "%d,%.6f", iteration, rmse);
    
    for (int i = 0; i < engine->num_params; ++i) {
        CalibParam *p = &engine->params[i];
        if (!p->enabled) continue;
        fprintf(engine->log_file, ",%.6f", p->current_val);
    }
    
    fprintf(engine->log_file, "\n");
    fflush(engine->log_file);
}
