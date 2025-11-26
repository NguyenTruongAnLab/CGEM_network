/**
 * @file main.c
 * @brief C-GEM Network Engine Main Entry Point
 * 
 * Initializes the network, runs simulation, and writes output.
 */

#include "define.h"
#include "network.h"
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
                "Usage: %s <path/to/case_config.txt>\n"
                "Example: %s INPUT/Cases/SaigonDongNai/case_config.txt\n",
                argv[0], argv[0]);
        return EXIT_FAILURE;
    }

    const char *case_config_path = argv[1];

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

    /* Transfer config to network */
    network.dt = case_config.dt;
    network.tidal_amplitude = case_config.tidal_amplitude;
    network.Q_river = case_config.Q_river;
    network.warmup_time = case_config.warmup_days * 86400.0;
    network.total_time = case_config.duration_days * 86400.0;

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
     * Run simulation
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
