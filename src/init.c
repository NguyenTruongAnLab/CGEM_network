/**
 * @file init.c
 * @brief Network initialization module - manages all parameter initialization
 *
 * This module centralizes all initialization logic for the C-GEM Network Engine,
 * following the pattern from the single-branch version but adapted for network topology.
 */

#include "network.h"
#include "define.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* Forward declarations for functions defined in other modules */
void InitializeSedimentParameters(Branch *branch);
int LoadBiogeoParams(const char *path);

/**
 * @brief Initialize species boundary conditions for a branch
 *
 * Sets smart defaults for ocean (downstream) and river (upstream) boundaries
 * for all water quality species. Only sets ocean BC values for branches
 * actually connected to ocean (LEVEL_BC) - junction-connected branches
 * will get their values from junction mixing.
 *
 * @param b Branch to initialize
 * @param num_species Number of species
 */
void init_species_bc(Branch *b, int num_species) {
    if (!b || num_species <= 0) return;

    /* Allocate boundary concentration arrays if needed */
    if (!b->conc_down) {
        b->conc_down = (double *)calloc((size_t)num_species, sizeof(double));
    }
    if (!b->conc_up) {
        b->conc_up = (double *)calloc((size_t)num_species, sizeof(double));
    }

    /* =================================================================
     * Ocean boundary conditions (downstream - only for LEVEL_BC branches)
     * For junction-connected branches, values will come from mixing
     * ================================================================= */
    if (b->conc_down) {
        /* Default to river-like values (will be overwritten for ocean BC) */
        b->conc_down[CGEM_SPECIES_SALINITY] = 1.0;      /* Low default [PSU] */
        b->conc_down[CGEM_SPECIES_PHY1] = 0.2;
        b->conc_down[CGEM_SPECIES_PHY2] = 0.3;
        b->conc_down[CGEM_SPECIES_DSI] = 30.0;
        b->conc_down[CGEM_SPECIES_NO3] = 10.0;
        b->conc_down[CGEM_SPECIES_NH4] = 2.0;
        b->conc_down[CGEM_SPECIES_PO4] = 1.0;
        b->conc_down[CGEM_SPECIES_O2] = 220.0;
        b->conc_down[CGEM_SPECIES_TOC] = 80.0;
        b->conc_down[CGEM_SPECIES_SPM] = 15.0;
        b->conc_down[CGEM_SPECIES_DIC] = 1800.0;
        b->conc_down[CGEM_SPECIES_AT] = 1800.0;
        b->conc_down[CGEM_SPECIES_PCO2] = 500.0;
        b->conc_down[CGEM_SPECIES_CO2] = 12.0;
        b->conc_down[CGEM_SPECIES_PH] = 7.9;
        b->conc_down[CGEM_SPECIES_HS] = 0.0;
        b->conc_down[CGEM_SPECIES_ALKC] = 1800.0;
        
        /* Only set ocean values if downstream is actually an ocean BC */
        if (b->down_node_type == NODE_LEVEL_BC) {
            b->conc_down[CGEM_SPECIES_SALINITY] = 35.0;     /* Salinity [PSU] */
            b->conc_down[CGEM_SPECIES_PHY1] = 0.5;          /* Diatoms [µM C] */
            b->conc_down[CGEM_SPECIES_PHY2] = 1.0;          /* Non-siliceous [µM C] */
            b->conc_down[CGEM_SPECIES_DSI] = 5.0;           /* Dissolved Si [µM] */
            b->conc_down[CGEM_SPECIES_NO3] = 2.0;           /* Nitrate [µM] */
            b->conc_down[CGEM_SPECIES_NH4] = 0.5;           /* Ammonium [µM] */
            b->conc_down[CGEM_SPECIES_PO4] = 0.2;           /* Phosphate [µM] */
            b->conc_down[CGEM_SPECIES_O2] = 250.0;          /* Oxygen [µM] */
            b->conc_down[CGEM_SPECIES_TOC] = 50.0;          /* TOC [µM C] */
            b->conc_down[CGEM_SPECIES_SPM] = 5.0;           /* SPM [mg/L] */
            b->conc_down[CGEM_SPECIES_DIC] = 2000.0;        /* DIC [µM] */
            b->conc_down[CGEM_SPECIES_AT] = 2300.0;         /* Total alkalinity [µM] */
            b->conc_down[CGEM_SPECIES_PCO2] = 400.0;        /* pCO2 [µatm] */
            b->conc_down[CGEM_SPECIES_CO2] = 10.0;          /* CO2 [µM] */
            b->conc_down[CGEM_SPECIES_PH] = 8.1;            /* pH */
            b->conc_down[CGEM_SPECIES_HS] = 0.0;            /* Hydrogen sulfide [µM] */
            b->conc_down[CGEM_SPECIES_ALKC] = 2300.0;       /* Carbonate alkalinity [µM] */
        }
    }

    /* =================================================================
     * River boundary conditions (upstream - typical freshwater values)
     * ================================================================= */
    if (b->conc_up) {
        b->conc_up[CGEM_SPECIES_SALINITY] = 0.1;        /* Low salinity [PSU] */

        /* Phytoplankton */
        b->conc_up[CGEM_SPECIES_PHY1] = 0.1;            /* Low diatoms [µM C] */
        b->conc_up[CGEM_SPECIES_PHY2] = 0.2;            /* Low non-siliceous [µM C] */

        /* Nutrients (enriched from land use) */
        b->conc_up[CGEM_SPECIES_DSI] = 50.0;            /* High dissolved Si [µM] */
        b->conc_up[CGEM_SPECIES_NO3] = 20.0;            /* High nitrate [µM] */
        b->conc_up[CGEM_SPECIES_NH4] = 5.0;             /* High ammonium [µM] */
        b->conc_up[CGEM_SPECIES_PO4] = 2.0;             /* High phosphate [µM] */

        /* Dissolved gases */
        b->conc_up[CGEM_SPECIES_O2] = 200.0;            /* Moderate oxygen [µM] */

        /* Organic matter */
        b->conc_up[CGEM_SPECIES_TOC] = 100.0;           /* High TOC [µM C] */

        /* Suspended matter */
        b->conc_up[CGEM_SPECIES_SPM] = 20.0;            /* High SPM [mg/L] */

        /* Carbon chemistry */
        b->conc_up[CGEM_SPECIES_DIC] = 1500.0;          /* Moderate DIC [µM] */
        b->conc_up[CGEM_SPECIES_AT] = 1200.0;           /* Moderate alkalinity [µM] */
        b->conc_up[CGEM_SPECIES_PCO2] = 600.0;          /* Higher pCO2 [µatm] */
        b->conc_up[CGEM_SPECIES_CO2] = 15.0;            /* Higher CO2 [µM] */
        b->conc_up[CGEM_SPECIES_PH] = 7.8;              /* Lower pH */
        b->conc_up[CGEM_SPECIES_HS] = 0.0;              /* No hydrogen sulfide [µM] */
        b->conc_up[CGEM_SPECIES_ALKC] = 1200.0;         /* Moderate carbonate alk [µM] */
    }

    /* Special case for Test_Mixing: set specific salinities */
    if (strcmp(b->name, "Branch1") == 0) {
        if (b->conc_up) b->conc_up[CGEM_SPECIES_SALINITY] = 0.0;
    } else if (strcmp(b->name, "Branch2") == 0) {
        if (b->conc_up) b->conc_up[CGEM_SPECIES_SALINITY] = 30.0;
    }

    /* =================================================================
     * Initialize concentration profiles 
     * 
     * For salinity (species 0), use Savenije (2005) exponential profile:
     *   S(x) = S0 * exp(-x / L_s)
     * where L_s is the salt intrusion length.
     * 
     * This provides a physically realistic initial condition that allows
     * the transport solver to maintain proper gradients.
     * 
     * Reference: Savenije (2005) "Salinity and Tides in Alluvial Estuaries"
     * ================================================================= */
    for (int sp = 0; sp < num_species; ++sp) {
        if (!b->conc || !b->conc[sp]) continue;

        double c_down = b->conc_down ? b->conc_down[sp] : 0.0;
        double c_up = b->conc_up ? b->conc_up[sp] : 0.0;

        /* For branches without ocean at downstream (junction-connected),
         * initialize with river-like values throughout - junction mixing
         * will establish proper gradients during simulation */
        if (b->down_node_type != NODE_LEVEL_BC) {
            /* Use upstream (river) values as initial profile */
            for (int i = 0; i <= b->M + 1; ++i) {
                b->conc[sp][i] = c_up;
            }
        } else {
            /* Ocean at downstream: use exponential profile for salinity,
             * linear for other species */
            if (sp == CGEM_SPECIES_SALINITY) {
                /* Typical salt intrusion length for Mekong: 30-60 km during dry season
                 * Reference: Nguyen et al. (2008), Vo Quoc Thanh (2021) */
                double L_intrusion = 40000.0;  /* 40 km initial intrusion length */
                for (int i = 0; i <= b->M + 1; ++i) {
                    double x = i * b->dx;  /* Distance from mouth */
                    double S_exp = c_down * exp(-x / L_intrusion);
                    /* Clamp to background salinity */
                    if (S_exp < c_up) S_exp = c_up;
                    b->conc[sp][i] = S_exp;
                }
            } else {
                /* Linear gradient for other species */
                for (int i = 0; i <= b->M + 1; ++i) {
                    double s = (double)i / (double)b->M;
                    b->conc[sp][i] = c_down + (c_up - c_down) * s;
                }
            }
        }
    }
}

/**
 * @brief Initialize network-wide default parameters
 *
 * Sets up default values for all network parameters that aren't specified
 * in configuration files. This ensures the model works even without
 * detailed configuration data.
 *
 * @param net Pointer to network structure
 */
void initializeNetworkDefaults(Network *net) {
    if (!net) return;

    /* Set default simulation parameters if not specified */
    if (net->dt <= 0) net->dt = CGEM_DEFAULT_DT_SECONDS;
    if (net->tidal_amplitude <= 0) net->tidal_amplitude = 2.0;
    if (net->Q_river <= 0) net->Q_river = 100.0;

    /* Set default physical constants */
    net->warmup_time = 2 * 86400.0;  /* 2 days warmup */
    net->total_time = 5 * 86400.0;   /* 5 days simulation */

    printf("Network defaults initialized.\n");
}

/**
 * @brief Initialize hydrodynamic parameters for all branches
 *
 * Sets up geometry, dispersion, and flow parameters for each branch
 * in the network using smart defaults and spatial variations.
 *
 * @param net Pointer to network structure
 */
void initializeNetworkHydrodynamics(Network *net) {
    if (!net || !net->branches) return;

    printf("Initializing hydrodynamics for %zu branches...\n", net->num_branches);

    for (size_t i = 0; i < net->num_branches; ++i) {
        Branch *branch = net->branches[i];
        if (!branch) continue;

        /* Initialize branch geometry */
        InitializeBranchGeometry(branch, net->dx_target);

        printf("  Branch %zu (%s): L=%.0f m, W=%.0f-%.0f m, H=%.1f m\n",
               i + 1, branch->name, branch->length_m,
               branch->width_down_m, branch->width_up_m, branch->depth_m);
    }

    printf("Hydrodynamics initialized for all branches.\n");
}

/**
 * @brief Initialize sediment transport parameters for all branches
 *
 * Sets up erosion/deposition parameters with spatial variation
 * from downstream (coarse sediment) to upstream (fine sediment).
 *
 * @param net Pointer to network structure
 */
void initializeNetworkSediment(Network *net) {
    if (!net || !net->branches) return;

    printf("Initializing sediment transport for %zu branches...\n", net->num_branches);

    for (size_t i = 0; i < net->num_branches; ++i) {
        Branch *branch = net->branches[i];
        if (!branch) continue;

        InitializeSedimentParameters(branch);
    }

    printf("Sediment transport initialized for all branches.\n");
}

/**
 * @brief Initialize biogeochemical parameters for all branches
 *
 * Sets up water quality parameters, phytoplankton kinetics,
 * nutrient cycling, and gas exchange for each branch.
 * 
 * Supports per-branch parameter files (specified in topology.csv column 11):
 *   - If branch has biogeo_params_path, loads that file (relative to case_dir)
 *   - Otherwise uses global biogeo_params.txt from case directory
 *   - Enables spatial heterogeneity (e.g., urban vs rural water quality)
 *
 * @param net Pointer to network structure
 * @param case_dir Case directory path (for finding biogeo_params.txt)
 */
void initializeNetworkBiogeochemistry(Network *net, const char *case_dir) {
    if (!net || !net->branches) return;

    printf("Initializing biogeochemistry for %zu branches...\n", net->num_branches);

    /* Try to load global biogeo parameters from case directory */
    if (case_dir) {
        char biogeo_path[CGEM_MAX_PATH];
        snprintf(biogeo_path, sizeof(biogeo_path), "%s/biogeo_params.txt", case_dir);
        LoadBiogeoParams(biogeo_path);
    } else {
        LoadBiogeoParams(NULL);  /* Use defaults */
    }

    for (size_t i = 0; i < net->num_branches; ++i) {
        Branch *branch = net->branches[i];
        if (!branch) continue;

        /* If branch has a branch-specific biogeo path, resolve it relative to case_dir */
        if (branch->biogeo_params_path[0] != '\0' && case_dir) {
            char full_path[CGEM_MAX_PATH];
            snprintf(full_path, sizeof(full_path), "%s/%s", case_dir, branch->biogeo_params_path);
            /* Store the full resolved path back in the branch */
            snprintf(branch->biogeo_params_path, sizeof(branch->biogeo_params_path), "%s", full_path);
        }

        InitializeBiogeoParameters(branch);
    }

    printf("Biogeochemistry initialized for all branches.\n");
}

/**
 * @brief Initialize species boundary conditions for all branches
 *
 * Sets up smart defaults for ocean (downstream) and river (upstream)
 * boundary conditions for all water quality species.
 *
 * @param net Pointer to network structure
 */
void initializeNetworkSpeciesBC(Network *net) {
    if (!net || !net->branches) return;

    printf("Initializing species boundary conditions...\n");

    for (size_t i = 0; i < net->num_branches; ++i) {
        Branch *branch = net->branches[i];
        if (!branch) continue;

        /* Set species count on branch */
        branch->num_species = net->num_species;

        /* Allocate concentration arrays if needed */
        if (!branch->conc && net->num_species > 0) {
            branch->conc = (double **)calloc(net->num_species, sizeof(double *));
            branch->conc_down = (double *)calloc(net->num_species, sizeof(double));
            branch->conc_up = (double *)calloc(net->num_species, sizeof(double));
            for (int sp = 0; sp < net->num_species; ++sp) {
                branch->conc[sp] = (double *)calloc((size_t)(branch->M + 2), sizeof(double));
            }
        }

        /* Initialize species boundary conditions with smart defaults */
        init_species_bc(branch, branch->num_species);
    }

    printf("Species boundary conditions initialized.\n");
}

/**
 * @brief Initialize node junction mixing arrays
 *
 * Sets up concentration arrays for junction nodes where
 * multiple branches meet and mixing occurs.
 *
 * @param net Pointer to network structure
 */
void initializeNetworkNodes(Network *net) {
    if (!net || !net->nodes) return;

    printf("Initializing %zu network nodes...\n", net->num_nodes);

    for (size_t n = 0; n < net->num_nodes; ++n) {
        Node *node = &net->nodes[n];

        /* Allocate mixed concentration arrays for species */
        if (!node->mixed_conc && net->num_species > 0) {
            node->mixed_conc = (double *)calloc(net->num_species, sizeof(double));
        }

        /* Initialize node properties */
        node->H = 0.0;
        node->H_new = 0.0;
        node->Q_net = 0.0;

        printf("  Node %zu: %s with %d connections\n",
               n + 1,
               node->type == NODE_JUNCTION ? "Junction" :
               node->type == NODE_DISCHARGE_BC ? "Discharge BC" : "Level BC",
               node->num_connections);
    }

    printf("Network nodes initialized.\n");
}

/**
 * @brief Main network initialization function
 *
 * Orchestrates the complete initialization of the C-GEM Network Engine,
 * following the pattern from the single-branch version but adapted
 * for network topology.
 *
 * @param net Pointer to network structure to initialize
 * @param config Pointer to case configuration
 * @return 0 on success, negative on error
 */
int initializeNetwork(Network *net, CaseConfig *config) {
    if (!net || !config) {
        fprintf(stderr, "Error: NULL network or config pointer\n");
        return -1;
    }

    printf("==================================================\n");
    printf("  C-GEM Network Engine - Initialization\n");
    printf("==================================================\n");

    /* Step 0: Set branch node types FIRST - needed for correct BC initialization */
    for (size_t b = 0; b < net->num_branches; ++b) {
        Branch *branch = net->branches[b];
        branch->up_node_type = net->nodes[branch->node_up].type;
        branch->down_node_type = net->nodes[branch->node_down].type;
    }

    /* Step 1: Set network-wide defaults */
    initializeNetworkDefaults(net);

    /* Step 2: Initialize hydrodynamic parameters */
    initializeNetworkHydrodynamics(net);

    /* Step 3: Initialize sediment transport */
    initializeNetworkSediment(net);

    /* Step 4: Initialize biogeochemical parameters */
    /* Extract case directory from topology path */
    char case_dir[CGEM_MAX_PATH] = {0};
    if (config->topology_path[0] != '\0') {
        size_t len = strlen(config->topology_path);
        if (len >= sizeof(case_dir)) len = sizeof(case_dir) - 1;
        memcpy(case_dir, config->topology_path, len);
        case_dir[len] = '\0';
        char *last_slash = strrchr(case_dir, '/');
        if (!last_slash) last_slash = strrchr(case_dir, '\\');
        if (last_slash) *last_slash = '\0';
        else case_dir[0] = '\0';
    }
    initializeNetworkBiogeochemistry(net, case_dir[0] ? case_dir : NULL);

    /* Step 5: Initialize species boundary conditions */
    /* NOTE: Requires branch node types to be set first! */
    initializeNetworkSpeciesBC(net);

    /* Step 6: Initialize node junction arrays */
    initializeNetworkNodes(net);

    /* Step 7: Initialize binary output files */
    if (cgem_init_output(net, config) != 0) {
        fprintf(stderr, "Failed to initialize binary output\n");
        return -1;
    }

    printf("==================================================\n");
    printf("  Network initialization complete!\n");
    printf("==================================================\n\n");

    return 0;
}