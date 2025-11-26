#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "network.h"
#ifdef _WIN32
#include <direct.h>
#define mkdir _mkdir
#endif

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

int cgem_write_csv(const char *path) {
    (void)path;
    fprintf(stderr, "cgem_write_csv() not implemented yet.\n");
    return -1;
}

int cgem_init_output(Network *net, CaseConfig *config) {
    if (!config->write_netcdf && !config->write_csv) return 0;

    for (size_t b = 0; b < net->num_branches; ++b) {
        Branch *branch = net->branches[b];
        if (!branch) continue;

        if (config->write_netcdf) {
            char filepath[8192];
            snprintf(filepath, sizeof(filepath), "%s/%s.bin", config->output_dir, branch->name);

            branch->bin_fp = fopen(filepath, "wb");
            if (!branch->bin_fp) {
                fprintf(stderr, "Failed to open binary file: %s\n", filepath);
                return -1;
            }

            // Write header
            int M = branch->M;
            int num_sp = net->num_species;
            double dx = branch->dx;

            fwrite(&M, sizeof(int), 1, branch->bin_fp);
            fwrite(&num_sp, sizeof(int), 1, branch->bin_fp);
            fwrite(&dx, sizeof(double), 1, branch->bin_fp);

            // Write X grid
            for (int i = 1; i <= M; ++i) {
                double x = i * dx;
                fwrite(&x, sizeof(double), 1, branch->bin_fp);
            }

            printf("Initialized binary output for %s\n", branch->name);
        }

        if (config->write_csv) {
            char csv_dir[8192];
            snprintf(csv_dir, sizeof(csv_dir), "%s/CSV", config->output_dir);
            mkdir(csv_dir); // ignore error if exists

            const char* var_names[21] = {"depth", "velocity", "waterlevel", "dispersion", "salinity", "phy1", "phy2", "dsi", "no3", "nh4", "po4", "o2", "toc", "spm", "dic", "at", "pco2", "co2", "ph", "hs", "alkc"};

            branch->num_csv_fps = 21;

            for (int v = 0; v < 21; ++v) {
                char filepath[8192];
                snprintf(filepath, sizeof(filepath), "%s/%s_%s.csv", csv_dir, branch->name, var_names[v]);

                branch->csv_fps[v] = fopen(filepath, "w");
                if (!branch->csv_fps[v]) {
                    fprintf(stderr, "Failed to open CSV file: %s\n", filepath);
                    return -1;
                }

                // Write header: Time_s,x1,x2,... (km units) using cell centers only (odd indices), exclude ghost cells
                fprintf(branch->csv_fps[v], "Time_s");
                double dx = branch->dx;
                int M = branch->M;
                int max_idx = (M >= 3) ? (M - 3) : 0;
                for (int j = 1; j <= max_idx; ++j) {
                    double x = j * dx;
                    fprintf(branch->csv_fps[v], ",%.0fkm", x / 1000.0);
                }
                fprintf(branch->csv_fps[v], "\n");
            }

            printf("Initialized CSV output for %s\n", branch->name);
        }
    }

    return 0;
}

int cgem_write_timestep(Network *net, double time_s) {
    for (size_t b = 0; b < net->num_branches; ++b) {
        Branch *branch = net->branches[b];
        if (!branch) continue;

        if (branch->bin_fp) {
            // Write time
            fwrite(&time_s, sizeof(double), 1, branch->bin_fp);

            // Write velocity (U)
            for (int i = 1; i <= branch->M; ++i) {
                double vel = sample_velocity_value(branch->velocity, i, branch->M);
                fwrite(&vel, sizeof(double), 1, branch->bin_fp);
            }

            // Write depth (H)
            for (int i = 1; i <= branch->M; ++i) {
                double depth = sample_center_value(branch->depth, i, branch->M);
                fwrite(&depth, sizeof(double), 1, branch->bin_fp);
            }

            // Write dispersion (D)
            for (int i = 1; i <= branch->M; ++i) {
                double disp = sample_center_value(branch->dispersion, i, branch->M);
                fwrite(&disp, sizeof(double), 1, branch->bin_fp);
            }

            // Write concentrations
            for (int sp = 0; sp < net->num_species; ++sp) {
                if (branch->conc && branch->conc[sp]) {
                    for (int i = 1; i <= branch->M; ++i) {
                        double cval = sample_center_value(branch->conc[sp], i, branch->M);
                        fwrite(&cval, sizeof(double), 1, branch->bin_fp);
                    }
                }
            }
        }

        if (branch->num_csv_fps > 0) {
            // Variables: 0:depth, 1:velocity, 2:waterlevel, 3:dispersion, 4+:species
            for (int v = 0; v < branch->num_csv_fps; ++v) {
                FILE *fp = branch->csv_fps[v];
                if (!fp) continue;

                // Write time (seconds)
                fprintf(fp, "%.0f", time_s);

                // Write values for cell centers (odd indices), exclude ghost cells
                int M = branch->M;
                int max_idx = (M >= 3) ? (M - 3) : 0;
                for (int j = 1; j <= max_idx; ++j) {
                    double value = 0.0;
                    if (v == 0) {
                        value = sample_center_value(branch->depth, j, M);
                    } else if (v == 1) {
                        value = sample_velocity_value(branch->velocity, j, M);
                    } else if (v == 2) {
                        value = sample_center_value(branch->waterLevel, j, M);
                    } else if (v == 3) {
                        value = sample_center_value(branch->dispersion, j, M);
                    } else if (v >= 4) {
                        int spidx = v - 4;
                        if (branch->conc && spidx < branch->num_species && branch->conc[spidx]) {
                            value = sample_center_value(branch->conc[spidx], j, M);
                        } else {
                            value = 0.0; // safe fallback
                        }
                    }
                    fprintf(fp, ",%.6f", value);
                }
                fprintf(fp, "\n");
            }
        }
    }

    return 0;
}

int cgem_close_output(Network *net) {
    for (size_t b = 0; b < net->num_branches; ++b) {
        Branch *branch = net->branches[b];
        if (branch) {
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
                /* Matrix CSVs removed - nothing else here */
        }
    }

    printf("Closed binary and CSV output files\n");
    return 0;
}

int cgem_write_netcdf(Network *net, CaseConfig *config) {
    // This is now a no-op since we stream data
    (void)net;
    (void)config;
    return 0;
}
/* write_branch_matrix_csv removed - per-variable CSVs (1 file per variable) are used instead */
