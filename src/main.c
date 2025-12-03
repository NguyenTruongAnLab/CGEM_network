/**
 * @file main.c
 * @brief C-GEM Network Engine Main Entry Point
 * 
 * Initializes the network, runs simulation, and writes output.
 * Supports calibration mode with --calibrate flag.
 */

#include "define.h"
#include "network.h"
#include "optimization/calibration.h"
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef _WIN32
#include <direct.h>
#define mkdir _mkdir
#else
#include <sys/stat.h>
#endif

/* Create directories recursively (portable). Returns 0 on success, -1 on error. */
static int mkdir_recursive(const char *path) {
    if (!path || !*path) return -1;
    char tmp[CGEM_MAX_PATH];
    strncpy(tmp, path, sizeof(tmp) - 1);
    tmp[sizeof(tmp)-1] = '\0';
    size_t len = strlen(tmp);

    /* Remove trailing separators */
    while (len > 0 && (tmp[len-1] == '/' || tmp[len-1] == '\\')) { tmp[len-1] = '\0'; --len; }
    if (len == 0) return -1;

    /* Try to create directory directly */
#ifdef _WIN32
    if (_mkdir(tmp) == 0) return 0;
    if (errno == EEXIST) return 0;
#else
    if (mkdir(tmp, 0755) == 0) return 0;
    if (errno == EEXIST) return 0;
#endif

    if (errno == ENOENT) {
        /* Create parent */
        char parent[CGEM_MAX_PATH];
        strncpy(parent, tmp, sizeof(parent) - 1);
        parent[sizeof(parent)-1] = '\0';
        char *p = parent + strlen(parent) - 1;
        while (p > parent && (*p == '/' || *p == '\\')) { *p = '\0'; --p; }
        char *sep = strrchr(parent, '/');
        char *sep2 = strrchr(parent, '\\');
        char *last = sep && sep2 ? (sep > sep2 ? sep : sep2) : (sep ? sep : sep2);
        if (!last) return -1;
        *last = '\0';
        if (mkdir_recursive(parent) != 0) return -1;
#ifdef _WIN32
        if (_mkdir(tmp) == 0) return 0;
        if (errno == EEXIST) return 0;
#else
        if (mkdir(tmp, 0755) == 0) return 0;
        if (errno == EEXIST) return 0;
#endif
    }
    return -1;
}

static void free_network(Network *network) {
    if (!network) return;
    
    if (network->branches) {
        for (size_t i = 0; i < network->num_branches; ++i) {
            if (network->branches[i]) {
                free_branch(network->branches[i], network->branches[i]->num_species);
            }
        }
        free(network->branches);
        network->branches = NULL;
    }
    
    if (network->nodes) {
        for (size_t i = 0; i < network->num_nodes; ++i) {
            free(network->nodes[i].forcing_time);
            free(network->nodes[i].forcing_value);
            free(network->nodes[i].mixed_conc);
        }
        free(network->nodes);
        network->nodes = NULL;
    }
    
    network->num_branches = 0;
    network->num_nodes = 0;
}

int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr,
                "Usage: %s <path/to/case_config.txt> [options]\n"
                "       %s <path/to/case_config.txt> --calibrate [--stage N]\n\n"
                "Options:\n"
                "  --calibrate     Run calibration instead of normal simulation\n"
                "  --stage N       Run only calibration stage N (1=hydro, 2=sed, 3=biogeo)\n"
                "  --max-iter N    Maximum optimizer iterations (default 100)\n"
                "  --verbose N     Verbosity level 0-3 (default 1)\n\n"
                "Example:\n"
                "  %s INPUT/Cases/Mekong_Delta_Full/case_config.txt\n"
                "  %s INPUT/Cases/Mekong_Delta_Full/case_config.txt --calibrate\n"
                "  %s INPUT/Cases/Mekong_Delta_Full/case_config.txt --calibrate --stage 1\n",
                argv[0], argv[0], argv[0], argv[0], argv[0]);
        return EXIT_FAILURE;
    }

    const char *case_config_path = argv[1];
    
    /* Parse command line options */
    int run_calibration = 0;
    int calib_stage = 0;  /* 0 = all stages */
    int max_iter = 100;
    int verbose = 1;
    
    for (int i = 2; i < argc; ++i) {
        if (strcmp(argv[i], "--calibrate") == 0) {
            run_calibration = 1;
        } else if (strcmp(argv[i], "--stage") == 0 && i + 1 < argc) {
            calib_stage = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--max-iter") == 0 && i + 1 < argc) {
            max_iter = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--verbose") == 0 && i + 1 < argc) {
            verbose = atoi(argv[++i]);
        }
    }

    /* =================================================================
     * Load case configuration
     * ================================================================= */
    CaseConfig case_config;
    memset(&case_config, 0, sizeof(case_config));
    
    if (LoadCaseConfig(case_config_path, &case_config) != 0) {
        fprintf(stderr, "Failed to load case configuration.\n");
        return EXIT_FAILURE;
    }

    /* Set defaults if not specified */
    if (case_config.dt <= 0) case_config.dt = CGEM_DEFAULT_DT_SECONDS;
    if (case_config.dx_meters <= 0) case_config.dx_meters = CGEM_DEFAULT_DX_METERS;
    if (case_config.tidal_amplitude <= 0) case_config.tidal_amplitude = 2.0;
    if (case_config.Q_river <= 0) case_config.Q_river = 100.0;

    /* Create output directory */
    #include <errno.h>

    {
        /* Create output directory and all parents. This avoids failures when
         * the OUTPUT directory or intermediate directories haven't been
         * created yet (common when running from the INPUT case folder). */
        int rc = mkdir_recursive(case_config.output_dir);
        if (rc != 0 && errno != EEXIST) {
            fprintf(stderr, "Warning: could not create output dir '%s' (errno=%d)\n", case_config.output_dir, errno);
        }
    }

    /* =================================================================
     * Initialize network
     * ================================================================= */
    Network network = {0};
    
    /* Set species count - include all biogeochemical species */
    int num_species = CGEM_NUM_SPECIES;  /* All species */
    network.num_species = num_species;

    /* Propagate case dx/dt to network */
    network.dt = case_config.dt;
    network.dx_target = case_config.dx_meters;

    /* Load topology */
    if (LoadTopology(case_config.topology_path, &network) != 0) {
        fprintf(stderr, "Failed to load topology from %s\n", case_config.topology_path);
        free_network(&network);
        return EXIT_FAILURE;
    }

    printf("Boundary path: '%s'\n", case_config.boundary_path);

    /* Load boundaries */
    if (LoadBoundaries(case_config.boundary_path, &network) != 0) {
        fprintf(stderr, "Failed to load boundaries from '%s'\n", case_config.boundary_path);
        free_network(&network);
        return EXIT_FAILURE;
    }

    /* Initialize the entire network */
    if (initializeNetwork(&network, &case_config) != 0) {
        fprintf(stderr, "Failed to initialize network\n");
        free_network(&network);
        return EXIT_FAILURE;
    }

    /* =================================================================
     * LATERAL LOADS (Land-Use Coupling)
     * Load spatially-explicit loads from urban, agriculture, aquaculture
     * =================================================================*/
    if (LoadLateralSources(&network, case_config.base_dir) != 0) {
        fprintf(stderr, "Warning: Failed to load lateral sources\n");
        /* Continue anyway - lateral loads are optional */
    }
    
    /* Load seasonal factors for rainfall-driven lateral load variation */
    if (LoadLateralSeasonalFactors(&network, case_config.base_dir) != 0) {
        fprintf(stderr, "Warning: Failed to load lateral seasonal factors\n");
        /* Continue anyway - defaults to constant factors (1.0) */
    }

    /* Transfer config to network */
    network.dt = case_config.dt;
    network.tidal_amplitude = case_config.tidal_amplitude;
    network.Q_river = case_config.Q_river;
    network.warmup_time = case_config.warmup_days * 86400.0;
    network.total_time = case_config.duration_days * 86400.0;
    network.reaction_mode = case_config.reaction_mode;  /* Transfer reaction mode */

    /* =================================================================
     * Print configuration
     * ================================================================= */
    printf("\n");
    printf("==============================================\n");
    printf("  C-GEM Network Engine\n");
    printf("==============================================\n");
    printf("  Case        : %s\n", case_config.case_name);
    printf("  Topology    : %s\n", case_config.topology_path);
    printf("  Branches    : %zu\n", network.num_branches);
    printf("  Nodes       : %zu\n", network.num_nodes);
    printf("  Species     : %d\n", num_species);
    printf("  Time step   : %.1f s\n", network.dt);
    printf("  Duration    : %d days\n", case_config.duration_days);
    printf("  Warmup      : %d days\n", case_config.warmup_days);
    printf("  Tidal amp   : %.2f m\n", network.tidal_amplitude);
    printf("  Q_river     : %.1f m³/s\n", network.Q_river);
    printf("  Reactions   : %s\n", network.reaction_mode ? "ON" : "OFF (transport only)");
    printf("  Output dir  : %s\n", case_config.output_dir);
    printf("==============================================\n\n");

    /* Print branch info */
    printf("Branches:\n");
    for (size_t i = 0; i < network.num_branches; ++i) {
        Branch *b = network.branches[i];
        if (!b) continue;
        
        printf("  [%zu] %-20s L=%.0f m, W_down=%.0f m, W_up=%.0f m, H=%.1f m, LC=%.0f m\n",
               i + 1, b->name, b->length_m, 
               b->width_down_m, b->width_up_m, b->depth_m,
               fabs(b->lc_convergence) >= 1e8 ? -1.0 : b->lc_convergence);
    }
    printf("\n");

    /* =================================================================
     * CALIBRATION MODE
     * ================================================================= */
    if (run_calibration) {
        printf("==============================================\n");
        printf("  CALIBRATION MODE\n");
        printf("==============================================\n\n");
        
        /* Create calibration engine */
        CalibrationEngine *calib = calib_create_engine();
        if (!calib) {
            fprintf(stderr, "Failed to create calibration engine\n");
            free_network(&network);
            return EXIT_FAILURE;
        }
        
        /* Set options */
        calib->max_iterations = max_iter;
        calib->verbose = verbose;
        
        /* Build paths to calibration files */
        char params_path[CGEM_MAX_PATH];
        char targets_path[CGEM_MAX_PATH];
        char seasonal_path[CGEM_MAX_PATH];
        snprintf(params_path, sizeof(params_path), "%s/calibration_params.csv", case_config.base_dir);
        snprintf(targets_path, sizeof(targets_path), "%s/calibration_targets.csv", case_config.base_dir);
        snprintf(seasonal_path, sizeof(seasonal_path), "%s/seasonal_targets.csv", case_config.base_dir);
        
        /* Load calibration configuration */
        if (calib_load_params(calib, params_path) != 0) {
            fprintf(stderr, "Failed to load calibration parameters from %s\n", params_path);
            calib_free_engine(calib);
            free_network(&network);
            return EXIT_FAILURE;
        }
        
        if (calib_load_objectives(calib, targets_path) != 0) {
            fprintf(stderr, "Failed to load calibration targets from %s\n", targets_path);
            calib_free_engine(calib);
            free_network(&network);
            return EXIT_FAILURE;
        }
        
        /* Load seasonal targets (optional - for time-series calibration) */
        calib_load_seasonal_targets(calib, seasonal_path);
        
        /* Set output directory */
        snprintf(calib->output_dir, sizeof(calib->output_dir), "%s/calibration", case_config.output_dir);
        mkdir_recursive(calib->output_dir);
        
        /* Run calibration */
        int result;
        if (calib_stage > 0) {
            /* Single stage */
            printf("Running Stage %d calibration only...\n\n", calib_stage);
            calib->current_stage = (CalibStage)calib_stage;
            calib_run_optimization(calib, &network, &case_config);
            result = 0;
        } else {
            /* Multi-stage */
            printf("Running full 3-stage calibration...\n\n");
            result = calib_run_multistage(calib, &network, &case_config);
        }
        
        /* Print and save results */
        calib_print_summary(calib);
        
        char results_path[CGEM_MAX_PATH];
        snprintf(results_path, sizeof(results_path), "%s/calibration_results.csv", calib->output_dir);
        calib_write_results(calib, results_path);
        printf("\nResults saved to: %s\n", results_path);
        
        /* Cleanup */
        calib_free_engine(calib);
        free_network(&network);
        
        return (result == 0) ? EXIT_SUCCESS : EXIT_FAILURE;
    }

    /* =================================================================
     * NORMAL SIMULATION MODE
     * ================================================================= */
    if (network_run_simulation(&network, &case_config) != 0) {
        fprintf(stderr, "Simulation failed.\n");
        free_network(&network);
        return EXIT_FAILURE;
    }

    printf("\nSimulation completed successfully.\n");
    
    free_network(&network);
    return EXIT_SUCCESS;
}
