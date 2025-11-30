/**
 * @file file.c
 * @brief C-GEM Network I/O Module - Binary and CSV output
 * 
 * Binary format (per branch):
 *   Header: M (int), num_hydro (int), num_species (int), num_reactions (int), dx (double)
 *           X grid (M doubles)
 *           Hydro names (num_hydro null-terminated strings)
 *           Species names (num_species null-terminated strings)
 *           Reaction names (num_reactions null-terminated strings)
 *   Per timestep: time (double)
 *                 hydro[num_hydro][M]
 *                 conc[num_species][M]
 *                 reaction_rates[num_reactions][M]
 * 
 * CSV format: One file per variable (same structure for hydro, species, reactions)
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "network.h"
#include "define.h"

#ifdef _WIN32
#include <direct.h>
#define mkdir _mkdir
#endif

/* --------------------------------------------------------------------------
 * Variable name arrays (used only for I/O, defined here to avoid warnings)
 * -------------------------------------------------------------------------- */

/* Hydrodynamic variable names - must match CGEM_HYDRO_* order in define.h */
static const char* HYDRO_NAMES[CGEM_NUM_HYDRO] = {
    "depth", "velocity", "waterlevel", "area", "width", "dispersion"
};

/* Species names - must match CGEM_SPECIES_* order in define.h */
static const char* SPECIES_NAMES[CGEM_NUM_SPECIES] = {
    "salinity", "phy1", "phy2", "dsi", "no3", "nh4", "po4", "o2",
    "toc", "spm", "dic", "at", "pco2", "co2", "ph", "hs", "alkc",
    /* RIVE multi-pool organic matter */
    "hd1", "hd2", "hd3", "hp1", "hp2", "hp3",
    /* RIVE bacteria */
    "bag", "bap",
    /* RIVE phosphorus */
    "pip",
    /* RIVE substrates */
    "dss"
};

/* Reaction names - must match CGEM_REACTION_* order in define.h */
static const char* REACTION_NAMES[CGEM_NUM_REACTIONS] = {
    "npp_no3", "npp_no3_1", "npp_no3_2", "npp_nh4", "npp_nh4_1", "npp_nh4_2",
    "gpp_1", "gpp_2", "npp", "phy_death", "phy_death_1", "phy_death_2",
    "si_cons", "aer_deg", "denit", "nit", "o2_ex", "o2_ex_s",
    "co2_ex", "co2_ex_s", "dic_react", "ta_react", "hs_react",
    "erosion_s", "erosion_v", "deposition_s", "deposition_v",
    /* RIVE reactions */
    "hydrolysis_hd1", "hydrolysis_hd2", "hydrolysis_hp1", "hydrolysis_hp2",
    "bac_uptake", "bac_resp", "bac_mort", "p_adsorption",
    "benthic_resp", "benthic_nh4", "benthic_po4", "benthic_o2", "benthic_dic"
};

/* --------------------------------------------------------------------------
 * Helper functions for staggered grid sampling
 * -------------------------------------------------------------------------- */

static double sample_center_value(double *arr, int idx, int M) {
    if (!arr) return 0.0;
    if ((idx % 2) != 0) return arr[idx];
    int left = (idx - 1 >= 1) ? idx - 1 : idx;
    int right = (idx + 1 <= M) ? idx + 1 : idx;
    if (left == idx && right == idx) return arr[idx];
    if (left == idx) return arr[right];
    if (right == idx) return arr[left];
    return 0.5 * (arr[left] + arr[right]);
}

static double sample_velocity_value(double *arr, int idx, int M) {
    if (!arr) return 0.0;
    if ((idx % 2) == 0) return arr[idx];
    int left = (idx - 1 >= 0) ? idx - 1 : idx;
    int right = (idx + 1 <= M + 1) ? idx + 1 : idx;
    if (left == idx && right == idx) return arr[idx];
    if (left == idx) return arr[right];
    if (right == idx) return arr[left];
    return 0.5 * (arr[left] + arr[right]);
}

/**
 * Get hydrodynamic variable value at a grid point
 * @param branch Branch pointer
 * @param hydro_idx Hydrodynamic variable index (CGEM_HYDRO_*)
 * @param grid_idx Grid index (1 to M)
 * @return Variable value
 */
static double get_hydro_value(Branch *branch, int hydro_idx, int grid_idx) {
    int M = branch->M;
    
    switch (hydro_idx) {
        case CGEM_HYDRO_DEPTH:
            return sample_center_value(branch->depth, grid_idx, M);
        case CGEM_HYDRO_VELOCITY:
            return sample_velocity_value(branch->velocity, grid_idx, M);
        case CGEM_HYDRO_WATERLEVEL:
            return sample_center_value(branch->waterLevel, grid_idx, M);
        case CGEM_HYDRO_AREA:
            return sample_center_value(branch->totalArea, grid_idx, M);
        case CGEM_HYDRO_WIDTH:
            return sample_center_value(branch->width, grid_idx, M);
        case CGEM_HYDRO_DISPERSION:
            return sample_center_value(branch->dispersion, grid_idx, M);
        default:
            return 0.0;
    }
}

/* --------------------------------------------------------------------------
 * Output initialization
 * -------------------------------------------------------------------------- */

int cgem_write_csv(const char *path) {
    (void)path;
    fprintf(stderr, "cgem_write_csv() not implemented yet.\n");
    return -1;
}

int cgem_init_output(Network *net, CaseConfig *config) {
    if (!config->write_netcdf && !config->write_csv) return 0;

    int num_hydro = CGEM_NUM_HYDRO;
    int num_species = net->num_species;
    int num_reactions = config->write_reaction_rates ? CGEM_NUM_REACTIONS : 0;

    for (size_t b = 0; b < net->num_branches; ++b) {
        Branch *branch = net->branches[b];
        if (!branch) continue;

        int M = branch->M;
        double dx = branch->dx;
        /* Output the entire physical reach (indices 1..M cover the resolved domain) */
        int max_idx = (M >= 1) ? M : 1;

        /* ============================================================
         * Binary output (for NetCDF conversion)
         * ============================================================ */
        if (config->write_netcdf) {
            char filepath[8192];
            snprintf(filepath, sizeof(filepath), "%s/%s.bin", config->output_dir, branch->name);

            branch->bin_fp = fopen(filepath, "wb");
            if (!branch->bin_fp) {
                fprintf(stderr, "Failed to open binary file: %s\n", filepath);
                return -1;
            }

            /* Write header */
            fwrite(&M, sizeof(int), 1, branch->bin_fp);
            fwrite(&num_hydro, sizeof(int), 1, branch->bin_fp);
            fwrite(&num_species, sizeof(int), 1, branch->bin_fp);
            fwrite(&num_reactions, sizeof(int), 1, branch->bin_fp);
            fwrite(&dx, sizeof(double), 1, branch->bin_fp);

            /* Write X grid */
            for (int i = 1; i <= M; ++i) {
                double x = i * dx;
                fwrite(&x, sizeof(double), 1, branch->bin_fp);
            }

            /* Write hydro names (null-terminated strings) */
            for (int h = 0; h < num_hydro; ++h) {
                const char *name = HYDRO_NAMES[h];
                fwrite(name, 1, strlen(name) + 1, branch->bin_fp);
            }

            /* Write species names (null-terminated strings) */
            for (int sp = 0; sp < num_species; ++sp) {
                const char *name = (sp < CGEM_NUM_SPECIES) ? SPECIES_NAMES[sp] : "unknown";
                fwrite(name, 1, strlen(name) + 1, branch->bin_fp);
            }

            /* Write reaction names (null-terminated strings) */
            for (int r = 0; r < num_reactions; ++r) {
                const char *name = (r < CGEM_NUM_REACTIONS) ? REACTION_NAMES[r] : "unknown";
                fwrite(name, 1, strlen(name) + 1, branch->bin_fp);
            }

            printf("Initialized binary output for %s (M=%d, hydro=%d, sp=%d, rxn=%d)\n", 
                   branch->name, M, num_hydro, num_species, num_reactions);
        }

        /* ============================================================
         * CSV output (for debugging)
         * ============================================================ */
        if (config->write_csv) {
            char csv_dir[CGEM_MAX_PATH + 16];
            snprintf(csv_dir, sizeof(csv_dir), "%s/CSV", config->output_dir);
            mkdir(csv_dir);

            /* Total: num_hydro + num_species + num_reactions */
            int total_vars = num_hydro + num_species + num_reactions;
            
            if (total_vars > 60) {
                fprintf(stderr, "Warning: Too many CSV variables (%d), truncating to 60\n", total_vars);
                total_vars = 60;
            }
            branch->num_csv_fps = total_vars;

            for (int v = 0; v < total_vars; ++v) {
                const char *var_name;
                if (v < num_hydro) {
                    var_name = HYDRO_NAMES[v];
                } else if (v < num_hydro + num_species) {
                    int sp = v - num_hydro;
                    var_name = (sp < CGEM_NUM_SPECIES) ? SPECIES_NAMES[sp] : "unknown";
                } else {
                    int r = v - num_hydro - num_species;
                    var_name = (r < CGEM_NUM_REACTIONS) ? REACTION_NAMES[r] : "unknown";
                }

                char filepath[CGEM_MAX_PATH + 128];
                snprintf(filepath, sizeof(filepath), "%s/%s_%s.csv", csv_dir, branch->name, var_name);

                branch->csv_fps[v] = fopen(filepath, "w");
                if (!branch->csv_fps[v]) {
                    fprintf(stderr, "Failed to open CSV file: %s\n", filepath);
                    return -1;
                }

                /* Write header: Time_s,x1,x2,... (km units) */
                fprintf(branch->csv_fps[v], "Time_s");
                for (int j = 1; j <= max_idx; ++j) {
                    double x = j * dx;
                    fprintf(branch->csv_fps[v], ",%.0fkm", x / 1000.0);
                }
                fprintf(branch->csv_fps[v], "\n");
            }

            printf("Initialized CSV output for %s (%d variables)\n", branch->name, total_vars);
        }
    }

    return 0;
}

/* --------------------------------------------------------------------------
 * Write timestep data
 * -------------------------------------------------------------------------- */

int cgem_write_timestep(Network *net, CaseConfig *config, double time_s) {
    int num_hydro = CGEM_NUM_HYDRO;
    int num_species = net->num_species;
    int num_reactions = (config && config->write_reaction_rates) ? CGEM_NUM_REACTIONS : 0;

    for (size_t b = 0; b < net->num_branches; ++b) {
        Branch *branch = net->branches[b];
        if (!branch) continue;

        int M = branch->M;
        /* Output the entire physical reach (indices 1..M cover the resolved domain) */
        int max_idx = (M >= 1) ? M : 1;

        /* ============================================================
         * Binary output
         * ============================================================ */
        if (branch->bin_fp) {
            /* Write time */
            fwrite(&time_s, sizeof(double), 1, branch->bin_fp);

            /* Write hydrodynamic variables (using generic function) */
            for (int h = 0; h < num_hydro; ++h) {
                for (int i = 1; i <= M; ++i) {
                    double val = get_hydro_value(branch, h, i);
                    fwrite(&val, sizeof(double), 1, branch->bin_fp);
                }
            }

            /* Write concentrations (all species) */
            for (int sp = 0; sp < num_species; ++sp) {
                for (int i = 1; i <= M; ++i) {
                    double cval = 0.0;
                    if (branch->conc && branch->conc[sp]) {
                        cval = sample_center_value(branch->conc[sp], i, M);
                    }
                    fwrite(&cval, sizeof(double), 1, branch->bin_fp);
                }
            }

            /* Write reaction rates (if reaction_rates array exists) */
            if (branch->reaction_rates) {
                for (int r = 0; r < num_reactions; ++r) {
                    for (int i = 1; i <= M; ++i) {
                        double rval = 0.0;
                        if (branch->reaction_rates[r]) {
                            rval = sample_center_value(branch->reaction_rates[r], i, M);
                        }
                        fwrite(&rval, sizeof(double), 1, branch->bin_fp);
                    }
                }
            }
        }

        /* ============================================================
         * CSV output
         * ============================================================ */
        if (branch->num_csv_fps > 0) {
            for (int v = 0; v < branch->num_csv_fps; ++v) {
                FILE *fp = branch->csv_fps[v];
                if (!fp) continue;

                /* Write time */
                fprintf(fp, "%.0f", time_s);

                /* Write values at all cell centers */
                for (int j = 1; j <= max_idx; ++j) {
                    double value = 0.0;

                    if (v < num_hydro) {
                        /* Hydrodynamic variable */
                        value = get_hydro_value(branch, v, j);
                    } else if (v < num_hydro + num_species) {
                        /* Species concentration */
                        int sp = v - num_hydro;
                        if (branch->conc && sp < branch->num_species && branch->conc[sp]) {
                            value = sample_center_value(branch->conc[sp], j, M);
                        }
                    } else {
                        /* Reaction rate */
                        int r = v - num_hydro - num_species;
                        if (branch->reaction_rates && r < CGEM_NUM_REACTIONS && branch->reaction_rates[r]) {
                            value = sample_center_value(branch->reaction_rates[r], j, M);
                        }
                    }

                    /* Use scientific notation for reaction rates (small values) */
                    if (v >= num_hydro + num_species) {
                        fprintf(fp, ",%.9e", value);
                    } else {
                        fprintf(fp, ",%.6f", value);
                    }
                }
                fprintf(fp, "\n");
            }
        }
    }

    return 0;
}

/* --------------------------------------------------------------------------
 * Close output files
 * -------------------------------------------------------------------------- */

int cgem_close_output(Network *net) {
    for (size_t b = 0; b < net->num_branches; ++b) {
        Branch *branch = net->branches[b];
        if (!branch) continue;

        if (branch->bin_fp) {
            fclose(branch->bin_fp);
            branch->bin_fp = NULL;
        }

        for (int v = 0; v < branch->num_csv_fps; ++v) {
            if (branch->csv_fps[v]) {
                fclose(branch->csv_fps[v]);
                branch->csv_fps[v] = NULL;
            }
        }
        branch->num_csv_fps = 0;
    }

    printf("Closed binary and CSV output files\n");
    return 0;
}

int cgem_write_netcdf(Network *net, CaseConfig *config) {
    /* No-op - Python script converts binary to NetCDF */
    (void)net;
    (void)config;
    return 0;
}
