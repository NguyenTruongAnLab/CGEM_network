/**
 * @file network_simulation.c
 * @brief C-GEM Network Simulation Driver
 * 
 * Main simulation loop with proper time stepping, output, and progress reporting.
 */

#include "network.h"
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>

/* M2 tidal period in seconds */
#define M2_PERIOD_S (12.42 * 3600.0)

/* Output intervals */
#define OUTPUT_INTERVAL_S (3600.0)   /* Write every hour */
#define PRINT_INTERVAL_S (86400.0)   /* Print progress every day */

/**
 * Write output to CSV files
 */
static void write_output(Network *net, CaseConfig *config, double time_s, FILE *fp_hydro, FILE *fp_conc) {
    if (!net) return;
    
    /* Store in time series arrays for NetCDF */
    if (config->write_netcdf || config->write_csv) {
        cgem_write_timestep(net, config, time_s);
    }
    
    /* Write hydrodynamic output for each branch */
    for (size_t b = 0; b < net->num_branches; ++b) {
        Branch *branch = net->branches[b];
        if (!branch) continue;
        
        /* Output at selected locations: head (1), middle cell center, mouth (M-1) */
        int mid = branch->M / 2;
        if (mid < 1) mid = 1;
        if (mid % 2 == 0) {
            if (mid + 1 < branch->M) mid += 1;
            else if (mid > 1) mid -= 1;
        }
        int locs[] = {1, mid, branch->M - 1};
        
        for (int k = 0; k < 3; ++k) {
            int i = locs[k];
            if (i < 1 || i > branch->M) continue;
            
                double x = i * branch->dx;
                if (fp_hydro) {
                fprintf(fp_hydro, "%.1f,%d,%s,%.1f,%.3f,%.4f,%.4f,%.2f",
                    time_s, branch->id, branch->name, x,
                    branch->depth[i],
                    branch->velocity[i],
                    branch->waterLevel[i],
                    branch->dispersion[i]);
                }
            
            /* Add salinity if available */
            if (branch->num_species > 0 && branch->conc && branch->conc[0]) {
                if (fp_hydro) fprintf(fp_hydro, ",%.4f", branch->conc[0][i]);
            }
            if (fp_hydro) fprintf(fp_hydro, "\n");
        }
    }
    
    /* Write concentration profiles (less frequent) */
    if (fp_conc && net->num_species > 0) {
        for (size_t b = 0; b < net->num_branches; ++b) {
            Branch *branch = net->branches[b];
            if (!branch || !branch->conc || !branch->conc[0]) continue;
            
            /* Full profile at odd indices */
            for (int j = 1; j <= branch->M - 1; j += 4) {
                double x = j * branch->dx;
                if (fp_conc) {
                    fprintf(fp_conc, "%.1f,%d,%s,%.1f",
                        time_s, branch->id, branch->name, x);
                }
                
                for (int sp = 0; sp < branch->num_species && sp < 5; ++sp) {
                    if (branch->conc[sp]) {
                        fprintf(fp_conc, ",%.4f", branch->conc[sp][j]);
                    }
                }
                if (fp_conc) fprintf(fp_conc, "\n");
            }
        }
    }
}

/**
 * Main simulation loop
 */
int network_run_simulation(Network *net, CaseConfig *config) {
    if (!net || !config) {
        return -1;
    }

    /* Get simulation parameters */
    double dt = net->dt;
    if (dt <= 0) dt = CGEM_DEFAULT_DT_SECONDS;
    
    double total_time = net->total_time;
    if (total_time <= 0) total_time = config->duration_days * 86400.0;
    if (total_time <= 0) total_time = 86400.0;  /* Default 1 day */
    
    double warmup_time = net->warmup_time;
    if (warmup_time <= 0) warmup_time = config->warmup_days * 86400.0;
    
    int total_steps = (int)(total_time / dt);
    int warmup_steps = (int)(warmup_time / dt);
    
    printf("Starting simulation: %d steps (%.1f days), warmup: %d steps (%.1f days)\n",
           total_steps, total_time / 86400.0,
           warmup_steps, warmup_time / 86400.0);

    FILE *fp_hydro = NULL;
    FILE *fp_conc = NULL;
    
    /* Initialize timing */
    double last_output_time = -OUTPUT_INTERVAL_S;
    double last_print_time = -PRINT_INTERVAL_S;
    clock_t start_clock = clock();

    /* Main time loop */
    for (int step = 0; step < total_steps; ++step) {
        double current_time = step * dt;
        int is_warmup = (step < warmup_steps);
        
        /* Progress output */
        if (current_time - last_print_time >= PRINT_INTERVAL_S || step == 0) {
            double elapsed_s = (double)(clock() - start_clock) / CLOCKS_PER_SEC;
            double day = current_time / 86400.0;
            double pct = 100.0 * step / total_steps;
            
            /* Get sample values from first branch */
            double sample_depth = 0.0, sample_vel = 0.0, sample_sal = 0.0;
            if (net->num_branches > 0 && net->branches[0]) {
                Branch *b = net->branches[0];
                int mid = b->M / 2;
                if (mid < 1) mid = 1;
                if (mid % 2 == 0) {
                    if (mid + 1 < b->M) mid += 1;
                    else if (mid > 1) mid -= 1;
                }
                sample_depth = b->depth[mid];
                sample_vel = b->velocity[mid];
                if (b->num_species > 0 && b->conc && b->conc[0]) {
                    sample_sal = b->conc[0][mid];
                }
            }
            
            printf("Day %6.1f (%5.1f%%) %s | H=%.2f m, U=%.3f m/s, S=%.1f | %.1f s elapsed\n",
                   day, pct, is_warmup ? "[WARMUP]" : "        ",
                   sample_depth, sample_vel, sample_sal, elapsed_s);
            
            last_print_time = current_time;
        }
        
        /* Solve hydrodynamics and transport */
        if (solve_network_step(net, current_time) != 0) {
            fprintf(stderr, "Failed at step %d (time=%.1f s)\n", step, current_time);
            fclose(fp_hydro);
            if (fp_conc) fclose(fp_conc);
            return -1;
        }
        
        /* Write output (after warmup) */
        if (!is_warmup && (current_time - last_output_time >= OUTPUT_INTERVAL_S)) {
            write_output(net, config, current_time - warmup_time, fp_hydro, fp_conc);
            last_output_time = current_time;
        }
    }

    /* Final output (per-variable CSVs) */
    double final_time = total_steps * dt;
    if (warmup_steps < total_steps) {
        /* Call cgem_write_timestep again for final per-variable snapshot */
        if (config->write_netcdf || config->write_csv) {
            cgem_write_timestep(net, config, final_time - warmup_time);
        }
        write_output(net, config, final_time - warmup_time, fp_hydro, fp_conc);
    }
    
    /* Close binary output files */
    if (config->write_netcdf) {
        cgem_close_output(net);
    }
    /* If write_netcdf enabled, attempt to convert binary files to netCDF and create plots */
    if (config->write_netcdf) {
        char cmd[4096];
        /* Check python availability */
#ifdef _WIN32
        const char *null_dev = "NUL";
#else
        const char *null_dev = "/dev/null";
#endif
        snprintf(cmd, sizeof(cmd), "python --version >%s 2>&1", null_dev);
        int rc_py = system(cmd);
        if (rc_py != 0) {
            fprintf(stderr, "Warning: Python not found on PATH, skipping .bin->.nc conversion and plotting.\n");
        } else {
            /* Convert .bin to .nc */
            char script_cmd[CGEM_MAX_PATH + 64];
            snprintf(script_cmd, sizeof(script_cmd), "python scripts/bin_to_nc.py %s", config->output_dir);
            printf("Converting .bin to .nc: %s\n", script_cmd);
            int rc_conv = system(script_cmd);
            if (rc_conv != 0) {
                fprintf(stderr, "Warning: bin_to_nc conversion returned non-zero (%d). netCDF files may be missing.\n", rc_conv);
            } else {
                /* Create plots */
                snprintf(script_cmd, sizeof(script_cmd), "python scripts/plot_netcdf.py %s", config->output_dir);
                printf("Creating plots from netCDFs: %s\n", script_cmd);
                int rc_plot = system(script_cmd);
                if (rc_plot != 0) {
                    fprintf(stderr, "Warning: plot_netcdf returned non-zero (%d).\n", rc_plot);
                }
            }
        }
    }
    
    double total_elapsed = (double)(clock() - start_clock) / CLOCKS_PER_SEC;
    printf("\nSimulation complete: %d steps in %.1f seconds (%.1f steps/s)\n",
           total_steps, total_elapsed, total_steps / total_elapsed);
    if (config->write_csv) {
        char csv_dir[CGEM_MAX_PATH + 16];
        snprintf(csv_dir, sizeof(csv_dir), "%s/CSV", config->output_dir);
        printf("CSV output written to: %s\n", csv_dir);
    }
    if (config->write_netcdf) {
        printf("Binary/netcdf output written to: %s\n", config->output_dir);
    }
    
    return 0;
}
