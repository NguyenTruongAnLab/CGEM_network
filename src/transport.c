/**
 * @file transport.c
 * @brief C-GEM Network Transport Module
 * 
 * Port of CGEM_Transport.f90 to C for multi-branch networks.
 * Implements:
 * - Van den Burgh dispersion coefficient calculation
 * - TVD (Total Variation Diminishing) advection with flux limiters
 * - Crank-Nicolson implicit dispersion scheme
 * - Open boundary condition handling
 * 
 * Grid convention (matching Fortran):
 * - Concentrations defined at odd indices (cell centers)
 * - Fluxes calculated at even indices (interfaces)
 * 
 * Reference: Savenije (2012), Volta et al. (2014), Gisen et al. (2015)
 */

#include "network.h"
#include "define.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    int up_cell;
    int down_cell;
    int up_neighbor;
    int down_neighbor;
    int up_vel_idx;
    int down_vel_idx;
    int up_ghost_primary;
    int up_ghost_secondary;
    int down_ghost_primary;
    int down_ghost_secondary;
    double up_sign;
    double down_sign;
} BranchBoundaryConfig;

static void determine_boundary_config(const Branch *b, BranchBoundaryConfig *cfg) {
    if (!b || !cfg) return;

    int M = b->M;
    double up_sign = 1.0;
    double down_sign = -1.0;
    int up_cell = M - 1;
    int down_cell = 1;

    /* When the downstream node is a level BC, the ocean is at the high-index end */
    if (b->down_node_type == NODE_LEVEL_BC) {
        down_cell = M - 1;
        up_cell = 1;
        down_sign = 1.0;
        up_sign = -1.0;
    }

    cfg->up_cell = up_cell;
    cfg->down_cell = down_cell;
    cfg->up_sign = up_sign;
    cfg->down_sign = down_sign;

    cfg->up_neighbor = up_cell + ((up_sign < 0.0) ? 2 : -2);
    if (cfg->up_neighbor < 1 || cfg->up_neighbor > M - 1) cfg->up_neighbor = up_cell;

    cfg->down_neighbor = down_cell + ((down_sign < 0.0) ? 2 : -2);
    if (cfg->down_neighbor < 1 || cfg->down_neighbor > M - 1) cfg->down_neighbor = down_cell;

    cfg->up_vel_idx = (up_cell == 1) ? 2 : M - 1;
    cfg->down_vel_idx = (down_cell == 1) ? 2 : M - 1;

    cfg->up_ghost_primary = (up_cell == 1) ? 0 : M;
    cfg->down_ghost_primary = (down_cell == 1) ? 0 : M;

    cfg->up_ghost_secondary = (up_cell == M - 1) ? M + 1 : -1;
    cfg->down_ghost_secondary = (down_cell == M - 1) ? M + 1 : -1;
}

static void set_dirichlet_value(double *arr, int idx, int ghost_a, int ghost_b, double value, int max_idx) {
    if (!arr || idx < 0 || idx > max_idx) return;
    arr[idx] = value;
    if (ghost_a >= 0 && ghost_a <= max_idx) arr[ghost_a] = value;
    if (ghost_b >= 0 && ghost_b <= max_idx) arr[ghost_b] = value;
}

/* --------------------------------------------------------------------------
 * Van den Burgh Dispersion Coefficient
 * -------------------------------------------------------------------------- */

/**
 * Calculate dispersion coefficient profile using Van den Burgh formulation
 * Matches Fortran: Calculate_Dispersion_Coefficient
 * 
 * Reference equations:
 * - Canter-Cremers number: N = -π*Q / (H0*B0)
 * - D0 = 26 * sqrt(N*g) * H0^1.5
 * - K = 4.38 * H0^0.36 * B0^-0.21 * LC^-0.14
 * - β = -(K * LCD * Q) / (D[i-1] * A[i-1])
 * - D[i] = D[i-1] * (1 - β*(exp(dx/LCD) - 1))
 * 
 * @param branch Branch to compute dispersion for
 * @param Q_total Total upstream discharge [m³/s]
 */
void ComputeDispersionCoefficient(Branch *branch, double Q_total) {
    if (!branch || branch->M <= 0) return;
    
    int M = branch->M;
    double dx = branch->dx;
    double g = CGEM_GRAVITY;
    
    /* Reference values at mouth (i=0) */
    double H0 = branch->depth[1];        /* Depth at mouth */
    double B0 = branch->width[1];        /* Width at mouth */
    double LC = branch->lc_convergence;  /* Convergence length */
    
    if (H0 < CGEM_MIN_DEPTH) H0 = CGEM_MIN_DEPTH;
    if (B0 < CGEM_MIN_WIDTH) B0 = CGEM_MIN_WIDTH;
    if (LC < 1000.0) LC = 1e9;  /* Very large for prismatic */
    
    /* Canter-Cremers estuary number */
    /* N = -π * Q / (H0 * B0) - negative Q means inflow from river */
    double N = CGEM_PI * fabs(Q_total) / (H0 * B0);
    
    /* Apply Van den Burgh tuning coefficient */
    N *= branch->vdb_coef;
    
    /* Effective dispersion at mouth: D0 = 26 * sqrt(N*g) * H0^1.5 */
    double D0 = 26.0 * sqrt(N * g) * pow(H0, 1.5);
    branch->D0 = D0;
    branch->dispersion[0] = D0;
    branch->dispersion[1] = D0;
    
    /* Calculate dispersion along the branch */
    for (int i = 2; i <= M + 1; ++i) {
        /* Local convergence length for cross-section */
        double AA_prev = branch->totalArea[i-1];
        double AA_curr = branch->totalArea[i];
        double LCD;
        
        if (AA_prev > 0 && AA_curr > 0 && fabs(AA_curr - AA_prev) > 1e-10) {
            LCD = -dx / log(AA_curr / AA_prev);
        } else {
            LCD = LC;
        }
        if (!isfinite(LCD) || fabs(LCD) < 100.0) LCD = LC;
        
        /* Van den Burgh K coefficient */
        /* K = 4.38 * H0^0.36 * B0^-0.21 * LC^-0.14 */
        double K = 4.38 * pow(H0, 0.36) * pow(B0, -0.21) * pow(fabs(LC), -0.14);
        K = CGEM_CLAMP(K, 0.0, 1.0);
        
        /* Beta coefficient */
        double D_prev = branch->dispersion[i-1];
        double A_prev = branch->totalArea[i-1];
        
        double beta = 0.0;
        if (D_prev > 0 && A_prev > 0) {
            beta = -(K * LCD * Q_total) / (D_prev * A_prev);
        }
        
        /* Dispersion at current location */
        double exp_term = exp(dx / LCD) - 1.0;
        double D_curr = D_prev * (1.0 - beta * exp_term);
        
        /* Ensure non-negative dispersion */
        if (D_curr < 0.0) D_curr = 0.0;
        if (D_curr > 10000.0) D_curr = 10000.0;  /* Cap at reasonable max */
        
        branch->dispersion[i] = D_curr;
    }
}

/* --------------------------------------------------------------------------
 * Open Boundary Conditions
 * -------------------------------------------------------------------------- */

/**
 * Apply open boundary conditions for transport
 * Matches Fortran: Open_Boundary_Conditions
 * 
 * @param branch Branch
 * @param species Species index
 * @param c_down Downstream (ocean) boundary concentration
 * @param c_up Upstream (river) boundary concentration
 * @param dt Time step [s]
 */
static void apply_open_boundaries(Branch *b, int species, double c_down, 
                                   double c_up, double dt) {
    int M = b->M;
    double dx = b->dx;
    double *conc = b->conc[species];
    BranchBoundaryConfig cfg;
    determine_boundary_config(b, &cfg);

    /* Distance to open-sea boundary condition */
    double OSBCdist = 2.0 * fmax(dx, 1.0);

    /* Handle downstream boundary */
    if (b->down_node_type == NODE_LEVEL_BC) {
        /* Dirichlet ocean boundary resides at the high-index side */
        set_dirichlet_value(conc, cfg.down_cell, cfg.down_ghost_primary,
                            cfg.down_ghost_secondary, c_down, M + 1);
    } else {
        double U_down = b->velocity[cfg.down_vel_idx];
        double inward = (U_down * cfg.down_sign <= 0.0) ? fabs(U_down) : 0.0;
        if (inward > 0.0) {
            conc[cfg.down_cell] = conc[cfg.down_cell] - 
                (conc[cfg.down_cell] - c_down) * inward * dt / OSBCdist;
        } else if (cfg.down_neighbor != cfg.down_cell) {
            conc[cfg.down_cell] = conc[cfg.down_cell] - 
                (conc[cfg.down_neighbor] - conc[cfg.down_cell]) * U_down * dt / dx;
        }
        set_dirichlet_value(conc, cfg.down_cell, cfg.down_ghost_primary,
                            cfg.down_ghost_secondary, conc[cfg.down_cell], M + 1);
    }

    /* Handle upstream boundary */
    if (b->up_node_type == NODE_DISCHARGE_BC || b->up_node_type == NODE_JUNCTION) {
        set_dirichlet_value(conc, cfg.up_cell, cfg.up_ghost_primary,
                            cfg.up_ghost_secondary, c_up, M + 1);
    } else {
        double U_up = b->velocity[cfg.up_vel_idx];
        double inward = (U_up * cfg.up_sign <= 0.0) ? fabs(U_up) : 0.0;
        if (inward > 0.0) {
            conc[cfg.up_cell] = conc[cfg.up_cell] - 
                (conc[cfg.up_cell] - c_up) * inward * dt / OSBCdist;
        } else if (cfg.up_neighbor != cfg.up_cell) {
            conc[cfg.up_cell] = conc[cfg.up_cell] - 
                (conc[cfg.up_cell] - conc[cfg.up_neighbor]) * U_up * dt / dx;
        }
        set_dirichlet_value(conc, cfg.up_cell, cfg.up_ghost_primary,
                            cfg.up_ghost_secondary, conc[cfg.up_cell], M + 1);
    }
}

/* --------------------------------------------------------------------------
 * TVD Advection Scheme
 * -------------------------------------------------------------------------- */

/**
 * Superbee flux limiter
 * φ(r) = max(0, min(2r, 1), min(r, 2))
 */
static double superbee_limiter(double r) {
    if (r <= 0.0) return 0.0;
    double a = CGEM_MIN(2.0 * r, 1.0);
    double b = CGEM_MIN(r, 2.0);
    return CGEM_MAX(0.0, CGEM_MAX(a, b));
}

/**
 * Calculate advective transport using TVD scheme
 * Matches Fortran: Calculate_Advection
 * 
 * Uses MUSCL-type reconstruction with Superbee limiter
 */
static void calculate_advection(Branch *b, int species, double dt) {
    int M = b->M;
    double dx = b->dx;
    double *conc = b->conc[species];
    
    /* Store old concentrations */
    double cold[CGEM_MAX_BRANCH_CELLS + 2];
    for (int j = 1; j <= M - 1; j += 2) {
        cold[j] = conc[j];
    }
    
    /* Calculate advective fluxes at interfaces (even indices) */
    double flux[CGEM_MAX_BRANCH_CELLS + 2];
    memset(flux, 0, sizeof(flux));
    
    for (int j = 1; j <= M - 2; j += 2) {
        int iface = j + 1;  /* Interface between j and j+2 */
        double vx = b->velocity[iface];
        double cfl = fabs(vx) * dt / (2.0 * dx);
        
        double f, rg, philen, conc_face;
        
        if (vx > 0.0) {
            /* Flow to the right (downstream to upstream in our convention) */
            f = cold[j+2] - cold[j];
            
            if (fabs(f) > 1e-35) {
                if (j > 1) {
                    /* Ratio of consecutive gradients */
                    rg = (cold[j] - cold[j-2]) / f;
                    /* Use Superbee limiter */
                    philen = superbee_limiter(rg);
                } else {
                    philen = 0.0;  /* No limiter at boundary */
                }
            } else {
                philen = 0.0;
            }
            
            /* Reconstructed face value */
            conc_face = cold[j] + 0.5 * (1.0 - cfl) * philen * f;
            flux[iface] = vx * b->totalArea[iface] * conc_face;
            
        } else if (vx < 0.0) {
            /* Flow to the left (upstream to downstream) */
            f = cold[j] - cold[j+2];
            
            if (fabs(f) > 1e-35) {
                if (j < M - 2) {
                    rg = (cold[j+2] - cold[j+4]) / f;
                    /* Use Superbee limiter */
                    philen = superbee_limiter(rg);
                } else {
                    philen = 0.0;
                }
            } else {
                philen = 0.0;
            }
            
            conc_face = cold[j+2] + 0.5 * (1.0 - cfl) * philen * f;
            flux[iface] = vx * b->totalArea[iface] * conc_face;
        }
        /* If vx == 0, flux remains 0 */
    }
    
    /* Update concentrations using flux divergence */
    for (int j = 3; j <= M - 2; j += 2) {
        double AA_old = b->totalArea_old2[j];
        double AA_new = b->totalArea[j];
        
        double rat = (AA_new > 0) ? AA_old / AA_new : 1.0;
        double rat1 = (AA_new > 0) ? dt / (2.0 * dx * AA_new) : 0.0;
        
        conc[j] = rat * cold[j] - rat1 * (flux[j+1] - flux[j-1]);
    }
}

/* --------------------------------------------------------------------------
 * Crank-Nicolson Dispersion Scheme
 * -------------------------------------------------------------------------- */

/**
 * Calculate dispersive transport using Crank-Nicolson implicit scheme
 * Matches Fortran: Calculate_Dispersion
 */
static void calculate_dispersion(Branch *b, int species, double dt) {
    int M = b->M;
    double dx = b->dx;
    double *conc = b->conc[species];
    
    /* Store old concentrations */
    double cold[CGEM_MAX_BRANCH_CELLS + 2];
    for (int j = 1; j <= M - 1; j += 2) {
        cold[j] = conc[j];
    }
    
    /* Coefficient: α = dt / (8 * dx²) */
    double alpha = dt / (8.0 * dx * dx);
    
    /* Tridiagonal coefficients */
    double a[CGEM_MAX_BRANCH_CELLS + 2];  /* Lower diagonal */
    double bb[CGEM_MAX_BRANCH_CELLS + 2];  /* Main diagonal */
    double c[CGEM_MAX_BRANCH_CELLS + 2];  /* Upper diagonal */
    double r[CGEM_MAX_BRANCH_CELLS + 2];  /* Explicit part (RHS before implicit) */
    double d[CGEM_MAX_BRANCH_CELLS + 2];  /* Full RHS */
    
    /* Boundary conditions: Dirichlet (concentration fixed) */
    a[1] = 0.0;
    bb[1] = 1.0;
    c[1] = 0.0;
    d[1] = conc[1];
    
    a[M-1] = 0.0;
    bb[M-1] = 1.0;
    c[M-1] = 0.0;
    d[M-1] = conc[M-1];
    
    /* Interior points */
    for (int i = 3; i <= M - 3; i += 2) {
        double D_im1 = b->dispersion[i-1];
        double D_ip1 = b->dispersion[i+1];
        double A_im1 = b->totalArea[i-1];
        double A_ip1 = b->totalArea[i+1];
        double A_i = b->totalArea[i];
        
        double g1 = (A_i > 0) ? D_im1 * A_im1 / A_i : 0.0;
        double g2 = (A_i > 0) ? D_ip1 * A_ip1 / A_i : 0.0;
        
        a[i] = -g1 * alpha;
        c[i] = -g2 * alpha;
        bb[i] = 1.0 + g2 * alpha + g1 * alpha;
        
        /* Explicit part coefficient */
        double r_coef = 1.0 - g2 * alpha - g1 * alpha;
        
        /* RHS: explicit half of Crank-Nicolson */
        r[i] = -c[i] * cold[i+2] + r_coef * cold[i] - a[i] * cold[i-2];
    }
    
    /* Copy r to d for the solve */
    for (int i = 1; i <= M - 1; i += 2) {
        d[i] = r[i];
    }
    d[1] = conc[1];
    d[M-1] = conc[M-1];
    
    /* Solve tridiagonal system (Thomas algorithm) */
    double gam[CGEM_MAX_BRANCH_CELLS + 2];
    
    if (fabs(bb[1]) < 1e-15) bb[1] = 1e-15;
    double bet = bb[1];
    conc[1] = d[1] / bet;
    
    for (int j = 3; j <= M - 1; j += 2) {
        gam[j] = c[j-2] / bet;
        bet = bb[j] - a[j] * gam[j];
        if (fabs(bet) < 1e-15) bet = 1e-15;
        conc[j] = (d[j] - a[j] * conc[j-2]) / bet;
    }
    
    for (int j = M - 3; j >= 1; j -= 2) {
        conc[j] = conc[j] - gam[j+2] * conc[j+2];
    }
}

/* --------------------------------------------------------------------------
 * Main Transport Function
 * -------------------------------------------------------------------------- */

/**
 * Main transport solver for a single branch
 * Matches Fortran: CGEM_Transport
 * 
 * @param branch Branch to compute transport for
 * @param dt Time step [s]
 * @return 0 on success, negative on error
 */
int Transport_Branch(Branch *branch, double dt) {
    if (!branch || branch->M <= 0) {
        return -1;
    }

    if (!branch->conc || branch->num_species <= 0) {
        return 0; /* Nothing to transport */
    }
    
    int M = branch->M;
    
    /* Compute dispersion coefficient profile */
    /* Estimate total discharge from upstream velocity */
    double Q_total = branch->velocity[M] * branch->totalArea[M];
    ComputeDispersionCoefficient(branch, Q_total);
    
    /* Process each species */
    for (int sp = 0; sp < branch->num_species; ++sp) {
        if (!branch->conc[sp]) continue;
        
        /* Get boundary concentrations */
        BranchBoundaryConfig cfg;
        determine_boundary_config(branch, &cfg);

        double default_down = branch->conc[sp][cfg.down_cell];
        double default_up = branch->conc[sp][cfg.up_cell];

        double c_down = (branch->conc_down) ? branch->conc_down[sp] : default_down;
        double c_up = (branch->conc_up) ? branch->conc_up[sp] : default_up;
        
        /* Apply open boundary conditions */
        apply_open_boundaries(branch, sp, c_down, c_up, dt);
        
        /* Calculate advection */
        calculate_advection(branch, sp, dt);
        
        /* Calculate dispersion */
        calculate_dispersion(branch, sp, dt);
        
        /* Ensure non-negative concentrations */
        for (int i = 0; i <= M + 1; ++i) {
            if (branch->conc[sp][i] < 0.0) {
                branch->conc[sp][i] = 0.0;
            }
        }
    }
    
    return 0;
}

/* --------------------------------------------------------------------------
 * Output Functions
 * -------------------------------------------------------------------------- */

/**
 * Write branch state to output file
 */
void write_branch_output(FILE *fp, Branch *branch, double time) {
    if (!fp || !branch) return;
    
    int M = branch->M;
    
    /* Write header on first call (time == 0) */
    static int header_written = 0;
    if (time == 0.0 && !header_written) {
        fprintf(fp, "Time_s,Branch_ID,Cell,Distance_m,Width_m,Depth_m,Velocity_m_s,WaterLevel_m,Dispersion_m2_s");
        if (branch->num_species > 0 && branch->conc) {
            fprintf(fp, ",Salinity");
            if (branch->num_species > 1) fprintf(fp, ",Phy1");
            if (branch->num_species > 2) fprintf(fp, ",Phy2");
        }
        fprintf(fp, "\n");
        header_written = 1;
    }
    
    /* Write data for odd indices (cell centers) */
    for (int j = 1; j <= M - 1; j += 2) {
        double x = j * branch->dx;
        fprintf(fp, "%.1f,%d,%d,%.1f,%.2f,%.3f,%.4f,%.4f,%.2f",
                time, branch->id, j, x,
                branch->width[j], branch->depth[j],
                branch->velocity[j], branch->waterLevel[j],
                branch->dispersion[j]);
        
        if (branch->num_species > 0 && branch->conc) {
            fprintf(fp, ",%.4f", branch->conc[0][j]);
            if (branch->num_species > 1 && branch->conc[1]) 
                fprintf(fp, ",%.4f", branch->conc[1][j]);
            if (branch->num_species > 2 && branch->conc[2]) 
                fprintf(fp, ",%.4f", branch->conc[2][j]);
        }
        fprintf(fp, "\n");
    }
}
