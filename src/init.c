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
 * @param Q_river Network reference discharge [m³/s] for Savenije initialization
 */
void init_species_bc(Branch *b, int num_species, double Q_river) {
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
        
        /* ===============================================================
         * RIVE Multi-Pool Organic Matter (HD1-3, HP1-3)
         * Reference: Billen et al. (1994), Garnier et al. (2000)
         * 
         * Mekong River: High sediment/DOC load from Tonle Sap + agriculture
         * Total DOC ~2-4 mg C/L in dry season (Liu et al., 2007)
         * POC ~1-3 mg C/L (bound to SPM)
         * =============================================================== */
        /* Dissolved OC pools [µM C = mg C/L × 83.3] */
        b->conc_up[CGEM_SPECIES_HD1] = 50.0;   /* Labile DOC: 0.6 mg C/L - fast degradation */
        b->conc_up[CGEM_SPECIES_HD2] = 100.0;  /* Semi-labile DOC: 1.2 mg C/L */
        b->conc_up[CGEM_SPECIES_HD3] = 150.0;  /* Refractory DOC: 1.8 mg C/L */
        
        /* Particulate OC pools [µM C] - bound to SPM */
        b->conc_up[CGEM_SPECIES_HP1] = 30.0;   /* Labile POC: 0.36 mg C/L */
        b->conc_up[CGEM_SPECIES_HP2] = 60.0;   /* Semi-labile POC: 0.72 mg C/L */
        b->conc_up[CGEM_SPECIES_HP3] = 80.0;   /* Refractory POC: 0.96 mg C/L */
        
        /* ===============================================================
         * RIVE Heterotrophic Bacteria
         * Reference: Garnier et al. (1992), Servais & Billen (1989)
         * 
         * Typical river bacteria: 0.05-0.5 mg C/L
         * =============================================================== */
        b->conc_up[CGEM_SPECIES_BAG] = 0.1;    /* Attached bacteria [mg C/L] */
        b->conc_up[CGEM_SPECIES_BAP] = 0.05;   /* Free bacteria [mg C/L] */
        
        /* ===============================================================
         * RIVE Phosphorus and Substrates
         * =============================================================== */
        b->conc_up[CGEM_SPECIES_PIP] = 1.0;    /* Adsorbed P [µM P] - on sediment */
        b->conc_up[CGEM_SPECIES_DSS] = 0.2;    /* Simple substrates [mg C/L] - bacteria food */
        
        /* ===============================================================
         * Greenhouse Gas Species
         * Reference: Borges et al. (2015), Hu et al. (2016) - Tropical rivers
         * 
         * Rivers are typically supersaturated in CO2, CH4, N2O
         * =============================================================== */
        b->conc_up[CGEM_SPECIES_NO2] = 0.5;    /* Nitrite [µM N] - nitrification intermediate */
        b->conc_up[CGEM_SPECIES_N2O] = 20.0;   /* N2O [nmol/L] - Mekong: 10-50 nmol/L */
        b->conc_up[CGEM_SPECIES_CH4] = 200.0;  /* CH4 [nmol/L] - Mekong: 100-500 nmol/L */
    }
    
    /* ===================================================================
     * Ocean boundary conditions for RIVE species (low/depleted values)
     * =================================================================== */
    if (b->conc_down && b->down_node_type == NODE_LEVEL_BC) {
        /* Dissolved OC pools - much lower in ocean */
        b->conc_down[CGEM_SPECIES_HD1] = 5.0;   /* Labile DOC - rapidly consumed */
        b->conc_down[CGEM_SPECIES_HD2] = 20.0;  /* Semi-labile DOC */
        b->conc_down[CGEM_SPECIES_HD3] = 50.0;  /* Refractory DOC - marine CDOM */
        
        /* Particulate OC - much lower in clear ocean water */
        b->conc_down[CGEM_SPECIES_HP1] = 2.0;
        b->conc_down[CGEM_SPECIES_HP2] = 5.0;
        b->conc_down[CGEM_SPECIES_HP3] = 10.0;
        
        /* Marine bacteria - lower than freshwater */
        b->conc_down[CGEM_SPECIES_BAG] = 0.02;
        b->conc_down[CGEM_SPECIES_BAP] = 0.02;
        
        /* Phosphorus and substrates */
        b->conc_down[CGEM_SPECIES_PIP] = 0.1;   /* Less adsorbed P in clear water */
        b->conc_down[CGEM_SPECIES_DSS] = 0.05;  /* Low substrates */
        
        /* GHG species - ocean is undersaturated relative to atmosphere */
        b->conc_down[CGEM_SPECIES_NO2] = 0.1;   /* Low nitrite in oxic ocean */
        b->conc_down[CGEM_SPECIES_N2O] = 8.0;   /* Near equilibrium ~7-10 nmol/L */
        b->conc_down[CGEM_SPECIES_CH4] = 3.0;   /* Near equilibrium ~2-5 nmol/L */
    }

    /* =================================================================
     * Initialize concentration profiles 
     * 
     * SAVENIJE ANALYTICAL SOLUTION (Steady-State Salt Intrusion)
     * ===========================================================
     * Uses the exact analytical solution from Savenije (1993, 2005) to
     * initialize salinity at the THEORETICAL EQUILIBRIUM position.
     * 
     * Van den Burgh equation (steady-state):
     *   S(x)/S0 = (D(x)/D0)^(1/K)
     * 
     * Where dispersion decays as:
     *   D(x)/D0 = 1 - β * (exp(x/a) - 1)
     *   β = (K * a * Q) / (A0 * D0)
     * 
     * This approach:
     * 1. Is scientifically rigorous (exact solution of governing equations)
     * 2. Minimizes warmup time (model starts near equilibrium)
     * 3. Serves as diagnostic (if salt too far inland, check Q or geometry)
     * 
     * Reference: Savenije (1993), Savenije (2005) Ch. 9, Nguyen (2008)
     * ================================================================= */
    for (int sp = 0; sp < num_species; ++sp) {
        if (!b->conc || !b->conc[sp]) continue;

        double c_down = b->conc_down ? b->conc_down[sp] : 0.0;
        double c_up = b->conc_up ? b->conc_up[sp] : 0.0;

        /* DEFAULT: Linear profile for non-salinity species */
        for (int i = 0; i <= b->M + 1; ++i) {
            double s = (b->M > 0) ? (double)i / (double)b->M : 0.0;
            b->conc[sp][i] = c_down + (c_up - c_down) * s;
        }

        /* SAVENIJE ANALYTICAL INITIALIZATION FOR SALINITY
         * Only for branches connected to ocean (LEVEL_BC) */
        if (sp == CGEM_SPECIES_SALINITY && b->down_node_type == NODE_LEVEL_BC) {
            
            /* =========================================================
             * 1. GEOMETRY PARAMETERS
             * ========================================================= */
            double H0 = b->depth_m;                    /* Depth at mouth [m] */
            double B0 = b->width_down_m;              /* Width at mouth [m] */
            double B_up = b->width_up_m;              /* Width upstream [m] */
            double L = b->length_m;                   /* Branch length [m] */
            double dx = b->dx;                        /* Grid spacing [m] */
            
            /* Cross-sectional area at mouth */
            double A0 = H0 * B0;
            if (A0 < 100.0) A0 = 100.0;  /* Minimum area */
            
            /* Convergence length 'a' [m]
             * From exponential geometry: A(x) = A0 * exp(-x/a)
             * Therefore: a = L / ln(A0/A_up)
             * Use width convergence since depth is approximately constant */
            double convergence_a;
            if (B0 > B_up && B_up > 0) {
                convergence_a = L / log(B0 / B_up);
            } else {
                convergence_a = 1e9;  /* Prismatic channel */
            }
            /* Use stored lc_convergence if available and valid */
            if (b->lc_convergence > 0 && b->lc_convergence < 1e8) {
                convergence_a = b->lc_convergence;
            }
            
            /* =========================================================
             * 2. VAN DEN BURGH COEFFICIENT K 
             * 
             * Use the SAME K as in transport.c for consistency.
             * K is from topology.csv column 11 (vdb_coef), typically 0.25-0.50
             * ========================================================= */
            double K_analytical = b->vdb_coef;
            if (K_analytical <= 0 || K_analytical > 1.0) {
                K_analytical = 0.40;  /* Default for Mekong */
            }
            
            /* =========================================================
             * 3. FRESHWATER DISCHARGE Q
             * 
             * Use network Q_river from config (RiverDischarge in case_config.txt)
             * Scale by cross-sectional area for branch-specific discharge.
             * ========================================================= */
            
            /* Use config-driven Q_river instead of hardcoded value */
            double Q_total = (Q_river > 0) ? Q_river : 100.0;  /* Fallback if not set */
            
            /* Scale Q by branch capacity (area relative to large main stem) */
            double A_reference = 15000.0;  /* Main Tien: ~1000m × 15m */
            double Q_branch = Q_total * (A0 / A_reference);
            if (Q_branch < 200.0) Q_branch = 200.0;   /* Minimum for small branches */
            if (Q_branch > 3000.0) Q_branch = 3000.0; /* Maximum for main stems */
            
            /* =========================================================
             * 4. DISPERSION AT MOUTH D0
             * Use same formula as transport.c for consistency
             * ========================================================= */
            double alpha = b->mixing_alpha;
            if (alpha <= 0) alpha = 0.10;  /* Conservative default */
            
            /* Estimate tidal velocity */
            double tidal_amp = 1.5;  /* Typical Mekong tidal amplitude [m] */
            double U_tidal = sqrt(9.81 * H0) * (tidal_amp / H0);
            
            double D0 = alpha * U_tidal * B0;
            /* D0 bounds (consistent with transport.c):
             * - Min 30 m²/s for ocean-connected branches (tidal mixing)
             * - Max 150 m²/s to prevent over-dispersion
             * Reference: Nguyen (2008), Savenije (2005) */
            if (D0 < 30.0) D0 = 30.0;    /* Minimum tidal mixing */
            if (D0 > 150.0) D0 = 150.0;  /* Maximum realistic value */
            
            /* =========================================================
             * 5. CALCULATE SAVENIJE STEADY-STATE PROFILE
             * 
             * Van den Burgh solution:
             *   D(x)/D0 = 1 - β * (exp(x/a) - 1)
             *   S(x)/S0 = (D(x)/D0)^(1/K)
             * 
             * where β = (K * a * Q) / (A0 * D0)
             * 
             * Salt intrusion limit: x_max where D(x)/D0 = 0
             *   exp(x_max/a) = 1 + 1/β
             *   x_max = a * ln(1 + 1/β)
             * ========================================================= */
            double beta = (K_analytical * convergence_a * Q_branch) / (A0 * D0);
            
            /* Calculate intrusion limit for diagnostics */
            double x_max = 0.0;
            if (beta > 0.01) {
                x_max = convergence_a * log(1.0 + 1.0/beta);
            } else {
                x_max = L;  /* Very low discharge - salt everywhere */
            }
            
            /* Print diagnostic info */
            printf("  Savenije Init %s: K=%.2f, a=%.0f m, Q=%.0f m³/s, D0=%.0f m²/s\n",
                   b->name, K_analytical, convergence_a, Q_branch, D0);
            printf("    β=%.4f, x_max=%.1f km (branch L=%.1f km)\n",
                   beta, x_max/1000.0, L/1000.0);
            
            /* Calculate salinity at each grid point 
             * Grid convention: i=1 is the DOWNSTREAM boundary (mouth)
             * Distance from mouth: x = (i-1) * dx
             * This ensures:
             *   i=0: ghost cell (x = -dx, set to ocean value)
             *   i=1: mouth (x = 0 km)
             *   i=M: upstream boundary
             *   i=M+1: ghost cell (set to river value)
             */
            for (int i = 0; i <= b->M + 1; ++i) {
                double x = (i <= 0) ? 0.0 : (double)(i - 1) * dx;  /* Distance from mouth [m] */
                
                /* Geometric factor: exp(x/a) - 1 */
                double geom_factor = exp(x / convergence_a) - 1.0;
                
                /* Dispersion ratio: D(x)/D0 */
                double D_ratio = 1.0 - beta * geom_factor;
                
                double S_val;
                if (D_ratio <= 0.01) {
                    /* Beyond salt intrusion limit → Fresh water */
                    S_val = c_up;
                } else {
                    /* Savenije solution: S(x) = S0 * (D(x)/D0)^(1/K) */
                    S_val = c_down * pow(D_ratio, 1.0/K_analytical);
                    
                    /* Clamp to river background */
                    if (S_val < c_up) S_val = c_up;
                    if (S_val > c_down) S_val = c_down;  /* Cannot exceed ocean */
                }
                
                b->conc[sp][i] = S_val;
            }
            
            /* Ghost cells: i=0 and i=M+1 
             * These are outside the physical domain */
            b->conc[sp][0] = c_down;  /* Downstream ghost = ocean value */
            /* Upstream ghost already set by the loop (will be c_up if beyond x_max) */
        }
        /* For junction-connected branches, salinity stays at river value
         * (already set in default loop above) - salt will enter via transport */
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
        init_species_bc(branch, branch->num_species, net->Q_river);
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
               n,  /* Use actual array index, matching topology node IDs */
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

    /* Step 0a: Copy critical config values to network BEFORE defaults
     * This ensures initialization routines have access to config values */
    if (config->Q_river > 0) {
        net->Q_river = config->Q_river;
    }
    if (config->tidal_amplitude > 0) {
        net->tidal_amplitude = config->tidal_amplitude;
    }

    /* Step 0b: Set branch node types FIRST - needed for correct BC initialization */
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

/**
 * @brief Reset network state for recalibration runs
 * 
 * Reinitializes hydrodynamic variables and species concentrations to their
 * initial values. Called before each calibration iteration to ensure each
 * run starts from the same baseline state.
 * 
 * @param net Network to reset
 */
void resetNetworkState(Network *net) {
    if (!net) return;
    
    /* Reset each branch to initial conditions */
    for (size_t b = 0; b < net->num_branches; ++b) {
        Branch *branch = net->branches[b];
        if (!branch) continue;
        
        int M = branch->M;
        
        /* Reset hydrodynamic state using reference depth */
        double H_mean = branch->depth_m;  /* Use reference depth from topology */
        
        for (int i = 0; i <= M + 1; ++i) {
            branch->depth[i] = H_mean;
            branch->velocity[i] = 0.0;
            branch->waterLevel[i] = 0.0;  /* Reference level */
            
            /* Reset min/max water level tracking for tidal range calculation */
            if (branch->waterLevel_min) branch->waterLevel_min[i] = 1e30;
            if (branch->waterLevel_max) branch->waterLevel_max[i] = -1e30;
            
            /* Recalculate areas/widths from exponential geometry */
            double x = (double)i * branch->dx;
            double width = branch->width_down_m * exp(-x / branch->lc_convergence);
            if (width < 10.0) width = 10.0;
            branch->width[i] = width;
            branch->freeArea[i] = branch->depth[i] * width;
            branch->totalArea[i] = branch->freeArea[i];
            if (branch->totalArea_old) branch->totalArea_old[i] = branch->freeArea[i];
            if (branch->totalArea_old2) branch->totalArea_old2[i] = branch->freeArea[i];
        }
        
        /* Reset residual velocity filter */
        if (branch->u_residual) {
            for (int i = 0; i <= M + 1; ++i) {
                branch->u_residual[i] = 0.0;
            }
        }
        
        /* Reset species concentrations to initial/boundary values */
        /* For estuarine species like salinity, use linear gradient between boundaries */
        for (int sp = 0; sp < net->num_species; ++sp) {
            if (!branch->conc || !branch->conc[sp]) continue;
            
            double conc_up = branch->conc_up ? branch->conc_up[sp] : 0.0;
            double conc_down = branch->conc_down ? branch->conc_down[sp] : 0.0;
            
            /* Salinity (species 0) needs special handling - use gradient */
            /* For other species, uniform upstream value is appropriate */
            if (sp == CGEM_SPECIES_SALINITY && fabs(conc_down - conc_up) > 0.1) {
                /* Linear gradient from downstream (i=1) to upstream (i=M) */
                for (int i = 0; i <= M + 1; ++i) {
                    double frac = (double)i / (double)(M + 1);  /* 0 at downstream, 1 at upstream */
                    branch->conc[sp][i] = conc_down + frac * (conc_up - conc_down);
                }
            } else {
                /* For non-salinity species or uniform conditions, use upstream value */
                for (int i = 0; i <= M + 1; ++i) {
                    branch->conc[sp][i] = conc_up;
                }
            }
        }
        
        /* Clear reaction rates */
        for (int r = 0; r < branch->num_reactions; ++r) {
            if (branch->reaction_rates && branch->reaction_rates[r]) {
                for (int i = 0; i <= M + 1; ++i) {
                    branch->reaction_rates[r][i] = 0.0;
                }
            }
        }
        
        /* Reset sediment bed state */
        if (branch->bed_mass) {
            for (int i = 0; i <= M + 1; ++i) {
                branch->bed_mass[i] = 0.0;
            }
        }
        if (branch->bed_oc) {
            for (int i = 0; i <= M + 1; ++i) {
                branch->bed_oc[i] = 0.0;
            }
        }
        if (branch->benthic_zf) {
            for (int i = 0; i <= M + 1; ++i) {
                branch->benthic_zf[i] = branch->zf_init;
            }
        }
    }
    
    /* Reset node water levels */
    for (size_t n = 0; n < net->num_nodes; ++n) {
        Node *node = &net->nodes[n];
        
        /* Reset to initial forcing value if available, else 0 */
        if (node->forcing_len > 0 && node->forcing_value) {
            node->H = node->forcing_value[0];
        } else {
            node->H = 0.0;
        }
        node->H_new = node->H;
    }
}