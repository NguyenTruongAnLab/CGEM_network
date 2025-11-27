/**
 * @file hydrodynamics.c
 * @brief C-GEM Network Hydrodynamics Module
 * 
 * Port of CGEM_Hydrodynamics.f90 to C for multi-branch networks.
 * Implements the staggered-grid hydrodynamic solver with:
 * - Tridiagonal matrix assembly and solution
 * - Iterative convergence for coupled continuity-momentum equations
 * - Tidal forcing at ocean boundary
 * - Discharge boundary at upstream
 * 
 * Grid convention (matching Fortran):
 * - Even indices (0, 2, 4, ...): Interfaces (velocity U, fluxes)
 * - Odd indices (1, 3, 5, ...): Cell centers (area AA, concentrations)
 * - M must be even
 * 
 * Reference: Savenije (2012), Gisen et al. (2015)
 */

#include "network.h"
#include "define.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* --------------------------------------------------------------------------
 * Helper functions
 * -------------------------------------------------------------------------- */

/**
 * Compute tidal elevation at the ocean boundary
 * @param t Current time [s]
 * @param amplitude Tidal amplitude [m]
 * @param period Tidal period [s] (default M2 ≈ 44712 s)
 * @return Tidal elevation [m]
 */
double TidalElevation(double t, double amplitude, double period) {
    double omega = 2.0 * CGEM_PI / period;
    return amplitude * sin(omega * t);
}

/**
 * Calculate total discharge including tributaries at a grid location
 * @param branch Pointer to branch
 * @param loc Grid location index
 * @param Q_river Base river discharge [m³/s]
 * @return Total discharge [m³/s]
 */
double TotalDischarge(Branch *branch, int loc, double Q_river) {
    (void)branch;
    (void)loc;
    /* For now, return base discharge. 
       TODO: Add tributary contributions based on location */
    return Q_river;
}

/* --------------------------------------------------------------------------
 * Geometry Initialization
 * -------------------------------------------------------------------------- */

/**
 * Initialize branch geometry arrays (width, depth, reference area, Chezy)
 * Matches CGEM_initialisation.f90 -> initialise_CGEM_arrays
 */
void InitializeBranchGeometry(Branch *branch, double target_dx) {
    if (!branch) return;
    if (!(target_dx > 0.0)) target_dx = CGEM_DEFAULT_DX_METERS; // default
    
    // 2. Calculate M based on Length and Target DX
    // M must be Even for staggered grid
    double original_length = branch->length_m;
    if (branch->length_m > 0.0) {
        /* Use nearest even integer M to minimize length adjustment */
        double steps_f = branch->length_m / target_dx;
        int Mcalc = (int)lround(steps_f);
        if (Mcalc < 2) Mcalc = 2;
        if (Mcalc % 2 != 0) {
            int M_minus = (Mcalc - 1 >= 2) ? Mcalc - 1 : 2;
            int M_plus = Mcalc + 1;
            if (fabs(steps_f - (double)M_minus) <= fabs((double)M_plus - steps_f)) {
                Mcalc = M_minus;
            } else {
                Mcalc = M_plus;
            }
        }
        branch->M = Mcalc;
        // set branch length to be divisible by target_dx; adjust length accordingly
        branch->length_m = (double)branch->M * target_dx;
        if (fabs(original_length - branch->length_m) > 1e-6) {
            fprintf(stderr, "Warning: Adjusting Branch %s length from %.2f to %.2f to match grid dx=%.2f.\n", branch->name, original_length, branch->length_m, target_dx);
        }
        branch->dx = target_dx;
    } else {
        // Fallback if length missing
        branch->M = CGEM_DEFAULT_BRANCH_CELLS; 
        branch->dx = target_dx;
        branch->length_m = branch->M * branch->dx;
    }

    printf("  > Configured %s: L=%.0f, M=%d, dx=%.2f, RS=%.1f\n", 
           branch->name, branch->length_m, branch->M, branch->dx, 
           (branch->storage_ratio > 0.0) ? branch->storage_ratio : 1.0);

    int M = branch->M;
    double length = branch->length_m;

    /* Ensure valid geometry values */
    branch->width_up_m = fmax(branch->width_up_m, CGEM_MIN_WIDTH);
    branch->width_down_m = fmax(branch->width_down_m, CGEM_MIN_WIDTH);
    
    if (branch->depth_up_m <= 0) branch->depth_up_m = branch->depth_m;
    if (branch->depth_down_m <= 0) branch->depth_down_m = branch->depth_m;
    branch->depth_up_m = fmax(branch->depth_up_m, CGEM_MIN_DEPTH);
    branch->depth_down_m = fmax(branch->depth_down_m, CGEM_MIN_DEPTH);
    branch->depth_m = fmax(branch->depth_m, CGEM_MIN_DEPTH);

    /* Compute convergence length: LC = -L / ln(W_up/W_down) 
       Note: In Fortran, W decreases from mouth (W_lb) to head (W_ub)
       Here: width_down_m = mouth (downstream), width_up_m = head (upstream) */
    double ratio = branch->width_up_m / branch->width_down_m;
    if (ratio <= 0.0 || fabs(log(ratio)) < 1e-6) {
        branch->lc_convergence = 1e9;  /* Prismatic channel */
    } else {
        branch->lc_convergence = -length / log(ratio);
    }

    /* Initialize arrays along the branch */
    for (int i = 0; i <= M + 1; ++i) {
        double x = (double)i * branch->dx;
        if (x > length) x = length;
        double s = x / length;  /* Normalized distance [0,1] */

        /* Exponential width profile: B(x) = B0 * exp(-x/LC) */
        double width;
        if (fabs(branch->lc_convergence) >= 1e8) {
            /* Prismatic channel */
            width = branch->width_down_m;
        } else {
            width = branch->width_down_m * exp(-x / branch->lc_convergence);
        }
        width = fmax(width, CGEM_MIN_WIDTH);
        branch->width[i] = width;

        /* Linear depth profile (or constant) */
        double depth_ref = branch->depth_up_m + (branch->depth_down_m - branch->depth_up_m) * (1.0 - s);
        depth_ref = fmax(depth_ref, CGEM_MIN_DEPTH);
        
        /* Reference cross-section */
        branch->refArea[i] = width * depth_ref;
        
        /* Initialize free area to zero (no tidal perturbation yet) */
        branch->freeArea[i] = 0.0;
        
        /* Total area = reference + free */
        branch->totalArea[i] = branch->refArea[i];
        branch->totalArea_old[i] = branch->totalArea[i];
        branch->totalArea_old2[i] = branch->totalArea[i];
        
        /* Water depth */
        branch->depth[i] = depth_ref;
        branch->waterLevel[i] = 0.0;
        
        /* Initialize velocity to small positive value to kickstart flow
         * This represents expected downstream flow direction */
        branch->velocity[i] = 0.1;
        
        /* Chezy coefficient (can be spatially varying) */
        /* Linear interpolation from downstream (Chezy_lb) to upstream (Chezy_ub) */
        double chezy_down = branch->chezy;
        /* Upstream Chezy ratio: configurable, default 0.8 (was hardcoded 0.6) */
        /* Higher values = less friction change upstream */
        double chezy_ratio = 0.8;  /* TODO: Load from case config */
        double chezy_up = branch->chezy * chezy_ratio;
        if (chezy_down <= 0) chezy_down = 50.0;
        if (chezy_up <= 0) chezy_up = 40.0;
        branch->chezyArray[i] = chezy_down + (chezy_up - chezy_down) * s;
        
        /* Initialize dispersion */
        branch->dispersion[i] = 200.0;
    }

    /* Initialize Van den Burgh coefficient */
    branch->vdb_coef = 1.0;

    /* Initialize concentrations if species exist */
    if (branch->num_species > 0 && branch->conc) {
        for (int sp = 0; sp < branch->num_species; ++sp) {
            if (!branch->conc[sp]) continue;
            
            double c_down = (branch->conc_down) ? branch->conc_down[sp] : 0.0;
            double c_up = (branch->conc_up) ? branch->conc_up[sp] : 0.0;
            
            /* Linear initial profile */
            for (int i = 0; i <= M + 1; ++i) {
                double s = (double)i / (double)M;
                branch->conc[sp][i] = c_down + (c_up - c_down) * s;
            }
        }
    }
}

/* --------------------------------------------------------------------------
 * Hydrodynamic Solver
 * -------------------------------------------------------------------------- */

/**
 * Set boundary conditions for the current timestep
 * Set_new_boundary_conditions
 * 
 * SCIENTIFIC APPROACH: Only set water levels here. Velocity is computed
 * by update_boundary_velocity_MOC() using Method of Characteristics.
 * This separates the Dirichlet BC (level) from the momentum calculation.
 * 
 * Reference: Stelling (1984), Abbott & Minns (1998)
 */
static void set_boundary_conditions(Branch *b, double H_down, double H_up, double Q_upstream) {
    if (!b) return;

    int M = b->M;

    /* Downstream boundary: level BC (ocean) or junction water level */
    if (b->down_node_type == NODE_LEVEL_BC || b->down_node_type == NODE_JUNCTION) {
        /* Set water level at boundary (Dirichlet condition) */
        b->waterLevel[0] = H_down;
        b->waterLevel[1] = H_down;
        
        /* Update free area consistent with water level */
        if (b->width[1] > 0.0) {
            b->freeArea[1] = H_down * b->width[1];
            b->totalArea[1] = b->refArea[1] + b->freeArea[1];
        }
        /* Note: Velocity at boundary is computed in update_boundary_velocity_MOC() */
    }

    /* Upstream boundary: discharge BC, level BC, or junction water level */
    if (b->up_node_type == NODE_DISCHARGE_BC && Q_upstream > -9000.0) {
        /* Discharge BC: set velocity directly from Q (Neumann-type condition) */
        double areaM = (b->totalArea && b->totalArea[b->M] > 0.0) ? b->totalArea[b->M] : b->refArea[b->M];
        if (fabs(areaM) < 1e-10) areaM = 1e-6;
        b->velocity[b->M] = Q_upstream / areaM;
        b->velocity[b->M+1] = b->velocity[b->M];
    } else {
        /* Level BC or Junction: set water level (Dirichlet condition) */
        b->waterLevel[M] = H_up;
        b->waterLevel[M+1] = H_up;
        
        if (b->width[M] > 0.0) {
            b->freeArea[M] = H_up * b->width[M];
            b->totalArea[M] = b->refArea[M] + b->freeArea[M];
        }
        /* Note: Velocity at boundary is computed in update_boundary_velocity_MOC() */
    }

    /* Safety: ensure non-zero areas */
    if (b->totalArea[1] < 1e-6) b->totalArea[1] = b->refArea[1];
    if (b->totalArea[M] < 1e-6) b->totalArea[M] = b->refArea[M];
}

/**
 * Assemble coefficient matrix for tridiagonal system
 * Set_coefficient_matrix
 * 
 * Matrix structure:
 *   tri_lower[i] = lower diagonal (a)
 *   tri_diag[i] = main diagonal (b)
 *   tri_upper[i] = upper diagonal (c)
 *   tri_rhs[i] = right-hand side (d)
 */
static void assemble_matrix(Branch *b, double dt) {
    int M = b->M;
    double dx = b->dx;
    double g = CGEM_GRAVITY;
    /* Use per-branch storage ratio instead of global constant */
    double RS = (b->storage_ratio > 0.0) ? b->storage_ratio : CGEM_RS;
    double inv_dt = 1.0 / dt;
    double inv_2dx = 0.5 / dx;
    
    /* Interior points */
    for (int j = 3; j <= M - 3; j += 2) {
        /* Continuity equation at odd indices (cell centers) */
        b->tri_lower[j] = -b->totalArea[j-1] * inv_2dx;
        b->tri_diag[j] = RS * inv_dt;
        b->tri_upper[j] = b->totalArea[j+1] * inv_2dx;
        b->tri_rhs[j] = RS * b->freeArea[j] * inv_dt;
        
        /* Momentum equation at even indices (interfaces) */
        int i = j + 1;
        double U_i = b->velocity[i];
        double H_i = b->depth[i];
        double C_i = b->chezyArray[i];
        
        if (H_i < CGEM_MIN_DEPTH) H_i = CGEM_MIN_DEPTH;
        if (C_i < 10.0) C_i = 50.0;
        
        double friction = (fabs(U_i) / (C_i * C_i)) / H_i;
        
        /* Convective acceleration term: (U(i+2) - U(i-2)) / (4*g*dx)
         * Matching Fortran Set_coefficient_matrix
         * This term improves momentum balance for high Froude number flows */
        double convective = 0.0;
        if (i >= 4 && i <= M - 4) {
            convective = (b->velocity[i+2] - b->velocity[i-2]) / (4.0 * g * dx);
        }
        
        b->tri_lower[i] = -inv_2dx / b->width[i-1];
        b->tri_diag[i] = 1.0/(g*dt) + friction + convective;
        b->tri_upper[i] = inv_2dx / b->width[i+1];
        b->tri_rhs[i] = U_i / (g * dt);
    }
    
    /* Last continuity equation (j = M-1) */
    {
        int j = M - 1;
        b->tri_lower[j] = -b->totalArea[j-2] * inv_2dx;
        b->tri_diag[j] = RS * inv_dt;
        b->tri_upper[j] = 0.0;
        /* Include upstream discharge in RHS */
        double Q_up = b->velocity[M] * b->totalArea[M];
        b->tri_rhs[j] = RS * b->freeArea[j] * inv_dt - Q_up * inv_2dx;
    }
    
    /* First momentum equation (i = 2) */
    {
        int i = 2;
        double U_i = b->velocity[i];
        double H_i = b->depth[i];
        double C_i = b->chezyArray[i];
        
        if (H_i < CGEM_MIN_DEPTH) H_i = CGEM_MIN_DEPTH;
        if (C_i < 10.0) C_i = 50.0;
        
        double friction = (fabs(U_i) / (C_i * C_i)) / H_i;
        double convective = (b->velocity[4] - b->velocity[2]) / (4.0 * g * dx);
        
        /* Tidal boundary forcing enters RHS */
        double tidal_term = inv_2dx * b->freeArea[1] / b->width[1];
        
        b->tri_lower[i] = 0.0;
        b->tri_diag[i] = 1.0/(g*dt) + friction + convective;
        b->tri_upper[i] = inv_2dx / b->width[3];
        b->tri_rhs[i] = U_i / (g * dt) + tidal_term;
    }
}

/**
 * Solve the tridiagonal system using Thomas algorithm
 * Solve_tridiagonal_matrix
 */
static void solve_tridiagonal(Branch *b, double *solution) {
    int M = b->M;
    double *gam = b->tri_gam;
    
    /* Forward elimination */
    double beta = b->tri_diag[2];
    if (fabs(beta) < 1e-15) beta = 1e-15;
    solution[2] = b->tri_rhs[2] / beta;
    
    for (int j = 3; j <= M - 1; ++j) {
        gam[j] = b->tri_upper[j-1] / beta;
        beta = b->tri_diag[j] - b->tri_lower[j] * gam[j];
        if (fabs(beta) < 1e-15) beta = 1e-15;
        solution[j] = (b->tri_rhs[j] - b->tri_lower[j] * solution[j-1]) / beta;
    }
    
    /* Back substitution */
    for (int j = M - 2; j >= 2; --j) {
        solution[j] -= gam[j+1] * solution[j+1];
    }
}

/**
 * Check convergence of iteration
 * @return 1 if converged, 0 otherwise
 */
static int check_convergence(Branch *b, double *old_freeArea, double *old_velocity) {
    int M = b->M;
    
    /* Check free area (odd indices) */
    for (int j = 3; j <= M - 1; j += 2) {
        double diff = fabs(b->freeArea[j] - old_freeArea[j]);
        if (diff >= CGEM_TOL) return 0;
    }
    
    /* Check velocity (even indices) */
    for (int i = 2; i <= M - 2; i += 2) {
        double diff = fabs(b->velocity[i] - old_velocity[i]);
        if (diff >= CGEM_TOL) return 0;
    }
    
    return 1;
}

/**
 * Update arrays after solving
 * Update_hydrodynamic_arrays
 * 
 * Includes explicit momentum-based boundary velocity calculation.
 * This solves the discrete momentum equation at boundaries:
 *   (U_new - U_old)/dt = -g * dH/dx - g * Sf
 * 
 * Reference: Cunge et al. (1980), Abbott & Minns (1998)
 */
static void update_arrays(Branch *b, double *solution) {
    int M = b->M;
    
    /* Extract velocity (even indices) and free area (odd indices) */
    for (int j = 2; j <= M - 2; j += 2) {
        b->velocity[j] = solution[j];
        b->freeArea[j+1] = solution[j+1];
    }
    
    /* Update derived quantities at even (interface) points */
    for (int i = 2; i <= M - 2; i += 2) {
        double AAf_im1 = b->freeArea[i-1];
        double AAf_ip1 = b->freeArea[i+1];
        b->freeArea[i] = 0.5 * (AAf_im1 + AAf_ip1);
        
        double tmp_im1 = AAf_im1 + b->refArea[i-1];
        double tmp_ip1 = AAf_ip1 + b->refArea[i+1];
        b->totalArea[i] = 0.5 * (tmp_im1 + tmp_ip1);
        b->depth[i] = 0.5 * (tmp_im1/b->width[i-1] + tmp_ip1/b->width[i+1]);
    }
    
    /* Update at odd (center) points */
    for (int j = 3; j <= M - 1; j += 2) {
        b->totalArea[j] = b->freeArea[j] + b->refArea[j];
        b->depth[j] = b->totalArea[j] / b->width[j];
        b->velocity[j] = 0.5 * (b->velocity[j-1] + b->velocity[j+1]);
    }
    
    /* Upstream boundary extrapolation (area and depth) */
    b->totalArea[M] = 1.5 * (b->freeArea[M-1] + b->refArea[M-1]) - 
                      0.5 * (b->freeArea[M-3] + b->refArea[M-3]);
    b->depth[M] = b->totalArea[M] / b->width[M];
    b->freeArea[M] = 1.5 * b->freeArea[M-1] - 0.5 * b->freeArea[M-3];
    
    /* ======================================================================
     * BOUNDARY VELOCITY: PRESERVE JUNCTION-ADJUSTED VALUES
     * 
     * For junction-connected boundaries, we DON'T extrapolate from interior.
     * The junction solver will directly adjust these velocities for mass balance.
     * 
     * Only extrapolate for ocean boundaries where tidal flow is naturally
     * established by the interior solution.
     * ====================================================================== */
    
    /* 1. Downstream Boundary (Index 1/0) */
    if (b->down_node_type == NODE_LEVEL_BC) {
        /* Ocean boundary - extrapolate from interior */
        b->velocity[1] = b->velocity[2];
        b->velocity[0] = b->velocity[2];
    }
    /* For NODE_JUNCTION: velocity[1] is set by junction solver - don't overwrite */

    /* 2. Upstream Boundary (Index M) */
    if (b->up_node_type == NODE_LEVEL_BC) {
        /* Level BC - extrapolate from interior */
        b->velocity[M] = (M >= 2) ? b->velocity[M-2] : 0.1;
        b->velocity[M+1] = b->velocity[M];
    }
    /* For NODE_JUNCTION: velocity[M] is set by junction solver - don't overwrite */
    /* For NODE_DISCHARGE_BC: velocity[M] is already set in set_boundary_conditions */
}

/**
 * Main hydrodynamic solver for a single branch
 * CGEM_Hydrodynamics
 * 
 * @param branch Branch to solve
 * @param H_down Water level at downstream (ocean) boundary [m]
 * @param H_up Water level at upstream boundary [m] (not used if Q_up > 0)
 * @param Q_up Discharge at upstream boundary [m³/s]
 * @param dt Time step [s]
 * @return 0 on success, negative on error
 */
int Hyd_Branch(Branch *branch, double H_down, double H_up, double Q_up, double dt) {
    if (!branch || branch->M <= 0) return -1;
    
    int M = branch->M;
    (void)H_up;  /* Not used when discharge BC is applied */
    
    /* Save old cross-section for transport */
    for (int i = 0; i <= M + 1; ++i) {
        branch->totalArea_old2[i] = branch->totalArea[i];
    }
    
    /* Set boundary conditions */
    set_boundary_conditions(branch, H_down, H_up, Q_up);
    
    /* Allocate temporary arrays for convergence check */
    double old_freeArea[CGEM_MAX_BRANCH_CELLS + 2];
    double old_velocity[CGEM_MAX_BRANCH_CELLS + 2];
    double solution[CGEM_MAX_BRANCH_CELLS + 2];
    
    for (int i = 0; i <= M + 1; ++i) {
        old_freeArea[i] = branch->freeArea[i];
        old_velocity[i] = branch->velocity[i];
    }
    
    /* Iterative solution */
    int converged = 0;
    for (int iter = 0; iter < CGEM_MAX_ITER; ++iter) {
        /* Assemble matrix */
        assemble_matrix(branch, dt);
        
        /* Solve tridiagonal system */
        solve_tridiagonal(branch, solution);
        
        /* Update arrays (includes boundary velocity momentum update) */
        update_arrays(branch, solution);
        
        /* Check convergence */
        converged = check_convergence(branch, old_freeArea, old_velocity);
        
        if (converged) break;
        
        /* Save for next iteration */
        for (int i = 0; i <= M + 1; ++i) {
            old_freeArea[i] = branch->freeArea[i];
            old_velocity[i] = branch->velocity[i];
        }
    }
    
    /* Update final derived quantities */
    for (int i = 0; i <= M + 1; ++i) {
        /* Ensure velocity is not exactly zero (numerical stability) */
        if (fabs(branch->velocity[i]) < 1e-10) {
            branch->velocity[i] = 1e-10;
        }
        
        /* Save for next timestep */
        branch->totalArea_old[i] = branch->totalArea[i];
        
        /* Update water level relative to reference */
        branch->waterLevel[i] = branch->depth[i] - 
            (branch->refArea[i] / branch->width[i]);
    }
    
    return converged ? 0 : -2;
}
