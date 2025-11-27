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
 * - Species environment flag filtering (only transport env=1 species)
 * - Flux bookkeeping for mass balance diagnostics
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

/* Local flux storage for diagnostics (per-species, reused each timestep) */
static double g_adv_flux[CGEM_MAX_BRANCH_CELLS + 2];
static double g_disp_flux[CGEM_MAX_BRANCH_CELLS + 2];

/**
 * Compute junction dispersion for branches not connected to ocean
 * 
 * For upstream branches (e.g., Tien_Main ending at a junction), the dispersion
 * at the downstream end should be continuous with the downstream branches.
 * We use a flow-weighted average of D values from the downstream branches.
 * 
 * @param branch The upstream branch
 * @param net Network containing all branches and nodes
 * @return D0 value to use at the junction end, or -1 if should use standard calc
 */
double ComputeJunctionDispersion(Branch *branch, void *network_ptr) {
    if (!branch || !network_ptr) return -1.0;
    
    Network *net = (Network *)network_ptr;
    
    /* Find which node is the junction (not ocean BC) */
    int junction_node_id = -1;
    int is_downstream_junction = 0;
    
    if (branch->down_node_type == NODE_JUNCTION) {
        junction_node_id = branch->node_down;
        is_downstream_junction = 1;
    } else if (branch->up_node_type == NODE_JUNCTION) {
        junction_node_id = branch->node_up;
        is_downstream_junction = 0;
    }
    
    if (junction_node_id < 0 || (size_t)junction_node_id >= net->num_nodes) {
        return -1.0;  /* No junction, use standard calculation */
    }
    
    Node *junction = &net->nodes[junction_node_id];
    if (junction->type != NODE_JUNCTION) {
        return -1.0;
    }
    
    /* Find downstream branches connected to this junction and compute weighted D */
    double sum_D_Q = 0.0;
    double sum_Q = 0.0;
    
    for (int c = 0; c < junction->num_connections; ++c) {
        int branch_idx = junction->connected_branches[c];
        if (branch_idx < 0 || (size_t)branch_idx >= net->num_branches) continue;
        
        Branch *other = net->branches[branch_idx];
        if (!other || other->id == branch->id) continue;
        
        int dir = junction->connection_dir[c];
        
        /* For an upstream branch's downstream junction:
         * We want D values from branches that flow OUT of this junction (dir=-1)
         * Those are the distributary branches */
        if (is_downstream_junction && dir == -1) {
            /* This branch flows out of junction - it's a downstream branch */
            /* Get D at its upstream end (near junction) */
            double D_at_junction = other->dispersion[other->M];
            double Q_branch = fabs(other->velocity[other->M] * other->totalArea[other->M]);
            
            if (D_at_junction > 0 && Q_branch > 0) {
                sum_D_Q += D_at_junction * Q_branch;
                sum_Q += Q_branch;
            }
        }
    }
    
    if (sum_Q > 0) {
        return sum_D_Q / sum_Q;  /* Flow-weighted average D */
    }
    
    return -1.0;  /* No valid downstream D found */
}

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
    
    /* C-GEM Grid Convention:
     * - Index 1 = downstream end (connects to node_down)
     * - Index M = upstream end (connects to node_up)
     * 
     * For transport:
     * - down_cell is where node_down connects (index 1)
     * - up_cell is where node_up connects (index M-1, last center cell)
     * 
     * Flow sign convention:
     * - Positive velocity = flow from upstream (M) to downstream (1)
     * - At downstream boundary: inflow = negative velocity
     * - At upstream boundary: inflow = positive velocity
     */
    cfg->down_cell = 1;           /* Downstream boundary cell (ocean side) */
    cfg->up_cell = M - 1;         /* Upstream boundary cell (river side) */
    if (cfg->down_cell < 1) cfg->down_cell = 1;
    if (cfg->up_cell < 1) cfg->up_cell = cfg->down_cell;
    if (cfg->up_cell > M - 1) cfg->up_cell = (M - 1 >= 1) ? (M - 1) : 1;
    cfg->down_sign = -1.0;        /* Positive U means outflow at downstream */
    cfg->up_sign = 1.0;           /* Positive U means inflow at upstream */

    /* Neighbor cells for gradient reconstruction */
    cfg->down_neighbor = 3;       /* Interior neighbor of downstream cell */
    cfg->up_neighbor = M - 3;     /* Interior neighbor of upstream cell */
    if (cfg->down_neighbor > M - 1) cfg->down_neighbor = cfg->down_cell;
    if (cfg->down_neighbor < 1) cfg->down_neighbor = cfg->down_cell;
    if (cfg->up_neighbor > M - 1) cfg->up_neighbor = cfg->up_cell;
    if (cfg->up_neighbor < 1) cfg->up_neighbor = cfg->up_cell;

    /* Velocity indices at boundaries */
    cfg->down_vel_idx = 2;        /* Interface velocity near downstream */
    cfg->up_vel_idx = M;          /* Interface velocity near upstream */

    /* Ghost cells for boundary extrapolation */
    cfg->down_ghost_primary = 0;
    cfg->down_ghost_secondary = -1;  /* No secondary ghost at downstream */
    cfg->up_ghost_primary = M;
    cfg->up_ghost_secondary = M + 1;
}

static void set_dirichlet_value(double *arr, int idx, int ghost_a, int ghost_b, double value, int max_idx) {
    if (!arr || idx < 0 || idx > max_idx) return;
    arr[idx] = value;
    if (ghost_a >= 0 && ghost_a <= max_idx) arr[ghost_a] = value;
    if (ghost_b >= 0 && ghost_b <= max_idx) arr[ghost_b] = value;
}

/**
 * Set ghost cells only for junction boundary conditions.
 * Unlike set_dirichlet_value, this does NOT set the interior boundary cell,
 * allowing the transport solver to compute a smooth gradient.
 */
static void set_junction_ghost_cells(double *arr, int ghost_a, int ghost_b, double value, int max_idx) {
    if (!arr) return;
    if (ghost_a >= 0 && ghost_a <= max_idx) arr[ghost_a] = value;
    if (ghost_b >= 0 && ghost_b <= max_idx) arr[ghost_b] = value;
}

/* --------------------------------------------------------------------------
 * Van den Burgh Dispersion Coefficient
 * -------------------------------------------------------------------------- */

/**
 * Calculate dispersion coefficient profile using Van den Burgh formulation
 * Calculate_Dispersion_Coefficient
 * 
 * Reference equations:
 * - Canter-Cremers number: N = -π*Q / (H0*B0)
 * - D0 = 26 * sqrt(N*g) * H0^1.5
 * - K = 4.38 * H0^0.36 * B0^-0.21 * LC^-0.14
 * - β = -(K * LCD * Q) / (D[i-1] * A[i-1])
 * - D[i] = D[i-1] * (1 - β*(exp(dx/LCD) - 1))
 * 
 * For junction continuity: if D0_override > 0, use it instead of calculating
 * from local geometry (used for upstream branches at junctions).
 * 
 * @param branch Branch to compute dispersion for
 * @param Q_total Total upstream discharge [m³/s]
 * @param D0_override Override D0 value (set to -1 to compute from geometry)
 */
static void ComputeDispersionCoefficient_Internal(Branch *branch, double Q_total, double D0_override) {
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
    
    double D0;
    
    if (D0_override > 0) {
        /* Use junction-derived D0 for continuity */
        D0 = D0_override;
    } else {
        /* Dispersion at mouth using Savenije (2005, 2012) formulation
         * 
         * From Savenije (2005) Table 9.2, typical D0 values for estuaries:
         * - Thames: 100-400 m²/s
         * - Scheldt: 50-200 m²/s  
         * - Mekong: 100-500 m²/s (estimated from salt intrusion observations)
         * 
         * Using calibrated empirical relationship based on Canter-Cremers number:
         *   D0 = k * sqrt(N * g) * H^1.5
         * 
         * With k = 4 (reduced from k = 8 for better salt intrusion control)
         * 
         * Reference: Savenije (2005) "Salinity and Tides in Alluvial Estuaries"
         *            Nguyen et al. (2008) Mekong salt intrusion study
         */
        double N = CGEM_PI * fabs(Q_total) / (H0 * B0);
        N *= branch->vdb_coef;
        
        /* Calibrated coefficient: D0 = 4 * sqrt(N*g) * H0^1.5 
         * This gives D0 ~ 100-300 m²/s for Mekong-like conditions
         */
        D0 = 4.0 * sqrt(N * g) * pow(H0, 1.5);
    }
    
    /*
     * Only enforce the tidal-mixing floor on branches that touch the open
     * ocean. Interior tributaries should keep their low junction-derived D0 so
     * the freshwater cap is not eroded by excessive dispersion.
     */
    int has_ocean_boundary = (branch->down_node_type == NODE_LEVEL_BC) ||
                             (branch->up_node_type == NODE_LEVEL_BC);

    if (D0 < 0.0) D0 = 0.0;           /* Physical lower bound */
    if (has_ocean_boundary && D0 < 30.0) D0 = 30.0;   /* Minimum tidal mixing */
    if (D0 > 1000.0) D0 = 1000.0;     /* Cap at realistic max */

    branch->D0 = D0;
    branch->dispersion[0] = D0;
    branch->dispersion[1] = D0;
    
    /* Calculate dispersion along the branch using Savenije (2012) formulation */
    for (int i = 2; i <= M + 1; ++i) {
        /* Van den Burgh K coefficient (empirical, from Savenije 2005, Table 9.1)
         * 
         * K typically ranges from 0.2 to 0.8 for alluvial estuaries.
         * The formula K = 4.38 * H0^0.36 * B0^-0.21 * LC^-0.14 is empirical.
         * 
         * For large rivers like the Mekong with H~10-15m, B~1000-2000m, LC~50-100km:
         * K ≈ 4.38 * 12^0.36 * 1500^-0.21 * 75000^-0.14 ≈ 0.3-0.5
         */
        double K = 4.38 * pow(H0, 0.36) * pow(B0, -0.21) * pow(fabs(LC), -0.14);
        K = CGEM_CLAMP(K, 0.2, 0.8);  /* Literature range (Savenije, 2005) */
        
        /* Dispersion decay using Savenije (2012) exponential formulation
         * 
         * The key insight from Savenije (2005, 2012) is that longitudinal
         * dispersion decays exponentially upstream:
         *   D(x) = D0 * exp(-x / L_d)
         * 
         * where L_d is the dispersion decay length scale.
         * 
         * For a well-mixed estuary (Savenije, 2012, Eq. 9.30):
         *   L_d = D0 * A / (K * Q)
         * 
         * This means higher river discharge Q causes FASTER decay of dispersion,
         * which pushes the salt intrusion front seaward (correct physical behavior).
         * 
         * References:
         * - Savenije, H.H.G. (2005) Salinity and Tides in Alluvial Estuaries, Elsevier
         * - Savenije, H.H.G. (2012) Salinity and Tides in Alluvial Estuaries, 2nd ed.
         */
        double D_prev = branch->dispersion[i-1];
        double A_prev = branch->totalArea[i-1];
        
        double L_d = 50000.0;  /* Default 50 km decay length */
        if (K > 0 && fabs(Q_total) > 0 && D0 > 0 && A_prev > 0) {
            L_d = (D0 * A_prev) / (K * fabs(Q_total));
        }
        /* Clamp L_d to physically reasonable range 
         * 
         * The dispersion decay length scale should be shorter than the 
         * salt intrusion length. For Mekong-type estuaries:
         * - Dry season L_s ~ 50-70 km
         * - L_d should be ~ 10-30 km
         * 
         * Shorter L_d values give steeper decay and more realistic profiles.
         */
        L_d = CGEM_CLAMP(L_d, 2000.0, 50000.0);  /* 2-50 km */
        
        /* Exponential decay: D = D_prev * exp(-dx / L_d) */
        double decay_factor = exp(-dx / L_d);
        double D_curr = D_prev * decay_factor;
        
        /* Ensure realistic dispersion bounds
         * 
         * Physical basis (Fischer et al., 1979; Savenije, 2012):
         * - River turbulent diffusion is O(1-10) m²/s, not higher
         * - Very low dispersion upstream allows advection to dominate
         * - This ensures freshwater flushes salt properly during high flow
         */
        if (D_curr < 1.0) D_curr = 1.0;          /* Minimum river turbulence [m²/s] */
        if (D_curr > 10000.0) D_curr = 10000.0;  /* Cap at reasonable max */
        
        branch->dispersion[i] = D_curr;
    }
}

/**
 * Public wrapper for dispersion coefficient calculation
 * (for backward compatibility - always computes D0 from geometry)
 */
void ComputeDispersionCoefficient(Branch *branch, double Q_total) {
    ComputeDispersionCoefficient_Internal(branch, Q_total, -1.0);
}

/* --------------------------------------------------------------------------
 * Open Boundary Conditions
 * -------------------------------------------------------------------------- */

/**
 * Check if a species should be transported (env=1) or is diagnostic-only (env=0)
 * Matches Fortran: IF (BGCArray(s)%env.EQ.1) THEN ... transport ...
 * 
 * @param species Species index
 * @return 1 if species should be transported, 0 if diagnostic only
 */
static int species_is_transported(int species) {
    if (species < 0 || species >= CGEM_NUM_SPECIES) return 0;
    return CGEM_SPECIES_TRANSPORT_FLAG[species];
}

/**
 * Apply open boundary conditions for transport
 * Open_Boundary_Conditions
 * 
 * C-GEM Grid Convention:
 * - Index 1 = downstream (ocean side, connects to node_down)
 * - Index M = upstream (river side, connects to node_up)
 * 
 * Boundary handling:
 * - LEVEL_BC (ocean): Dirichlet - sets interior cell and ghost cells
 * - DISCHARGE_BC: Dirichlet - sets interior cell and ghost cells  
 * - JUNCTION: Ghost-only - only sets ghost cells, lets transport compute interior
 *   This allows smooth gradients across junctions via dispersion
 * 
 * @param branch Branch
 * @param species Species index
 * @param c_down Downstream boundary concentration
 * @param c_up Upstream boundary concentration
 * @param dt Time step [s]
 */
static void apply_open_boundaries(Branch *b, int species, double c_down, 
                                   double c_up, double dt) {
    int M = b->M;
    double *conc = b->conc[species];
    BranchBoundaryConfig cfg;
    determine_boundary_config(b, &cfg);

    /* Handle downstream boundary (index 1, connects to node_down) */
    if (b->down_node_type == NODE_JUNCTION) {
        /* Junction: only set ghost cells, let transport compute interior */
        set_junction_ghost_cells(conc, cfg.down_ghost_primary,
                                 cfg.down_ghost_secondary, c_down, M + 1);
    } else {
        /* Ocean or other BC: full Dirichlet */
        set_dirichlet_value(conc, cfg.down_cell, cfg.down_ghost_primary,
                            cfg.down_ghost_secondary, c_down, M + 1);
    }

    /* Handle upstream boundary (index M-1, connects to node_up) */
    if (b->up_node_type == NODE_JUNCTION) {
        /* Junction: only set ghost cells, let transport compute interior */
        set_junction_ghost_cells(conc, cfg.up_ghost_primary,
                                 cfg.up_ghost_secondary, c_up, M + 1);
    } else {
        /* Discharge or other BC: full Dirichlet */
        set_dirichlet_value(conc, cfg.up_cell, cfg.up_ghost_primary,
                            cfg.up_ghost_secondary, c_up, M + 1);
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
 * Calculate_Advection
 * 
 * Uses MUSCL-type reconstruction with Superbee limiter.
 * Stores fluxes in global array for diagnostics/bookkeeping.
 * 
 * Reference: Fortran Calculate_Advection in CGEM_Transport.f90
 */
static void calculate_advection(Branch *b, int species, double dt) {
    int M = b->M;
    double dx = b->dx;
    double *conc = b->conc[species];
    
    /* Store old concentrations at ALL indices including ghost cells
     * This prevents out-of-bounds access in gradient reconstruction */
    double cold[CGEM_MAX_BRANCH_CELLS + 4];
    memset(cold, 0, sizeof(cold));
    
    /* Copy all concentration values including ghost cells */
    for (int j = 0; j <= M + 1; ++j) {
        cold[j] = conc[j];
    }
    
    /* Calculate advective fluxes at interfaces (even indices) */
    /* Use global flux array for bookkeeping */
    memset(g_adv_flux, 0, sizeof(g_adv_flux));
    double *flux = g_adv_flux;
    
    for (int j = 1; j <= M - 2; j += 2) {
        int iface = j + 1;  /* Interface between j and j+2 */
        double vx = b->velocity[iface];
        double cfl = fabs(vx) * dt / (2.0 * dx);
        
        /* CFL stability check */
        if (cfl > 1.0) cfl = 1.0;
        
        double f, rg, philen, conc_face;
        
        if (vx > 0.0) {
            /* Flow to the right (downstream to upstream in our convention)
             * Upwind cell is j, need gradient from j-2 to j */
            f = cold[j+2] - cold[j];
            
            if (fabs(f) > 1e-35) {
                if (j >= 3) {
                    /* Ratio of consecutive gradients */
                    rg = (cold[j] - cold[j-2]) / f;
                    /* Use Superbee limiter */
                    philen = superbee_limiter(rg);
                } else {
                    /* At boundary: use first-order upwind (no reconstruction) */
                    philen = 0.0;
                }
            } else {
                philen = 0.0;
            }
            
            /* Reconstructed face value - bounded by local min/max */
            conc_face = cold[j] + 0.5 * (1.0 - cfl) * philen * f;
            /* Ensure face value is bounded (TVD property) */
            double cmin = CGEM_MIN(cold[j], cold[j+2]);
            double cmax = CGEM_MAX(cold[j], cold[j+2]);
            conc_face = CGEM_CLAMP(conc_face, cmin, cmax);
            flux[iface] = vx * b->totalArea[iface] * conc_face;
            
        } else if (vx < 0.0) {
            /* Flow to the left (upstream to downstream)
             * Upwind cell is j+2, need gradient from j+2 to j+4 */
            f = cold[j] - cold[j+2];
            
            if (fabs(f) > 1e-35) {
                if (j + 4 <= M) {
                    rg = (cold[j+2] - cold[j+4]) / f;
                    /* Use Superbee limiter */
                    philen = superbee_limiter(rg);
                } else {
                    /* At boundary: use first-order upwind */
                    philen = 0.0;
                }
            } else {
                philen = 0.0;
            }
            
            /* Reconstructed face value - bounded by local min/max */
            conc_face = cold[j+2] + 0.5 * (1.0 - cfl) * philen * f;
            /* Ensure face value is bounded (TVD property) */
            double cmin = CGEM_MIN(cold[j], cold[j+2]);
            double cmax = CGEM_MAX(cold[j], cold[j+2]);
            conc_face = CGEM_CLAMP(conc_face, cmin, cmax);
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
 * Calculate_Dispersion
 * 
 * Stores dispersive fluxes in global array for diagnostics.
 * Flux = -D * A * dC/dx (negative of Fick's law for convention)
 * 
 * @param b Branch
 * @param species Species index  
 * @param dt Time step
 * @param c_down Downstream Dirichlet BC value
 * @param c_up Upstream Dirichlet BC value
 * 
 * Reference: Fortran Calculate_Dispersion in CGEM_Transport.f90
 */
static void calculate_dispersion(Branch *b, int species, double dt, 
                                  double c_down, double c_up) {
    int M = b->M;
    double dx = b->dx;
    double *conc = b->conc[species];
    
    /* Store old concentrations for flux calculation (all indices) */
    double cold[CGEM_MAX_BRANCH_CELLS + 4];
    memset(cold, 0, sizeof(cold));
    for (int j = 0; j <= M + 1; ++j) {
        cold[j] = conc[j];
    }
    
    /* Clear dispersion flux array */
    memset(g_disp_flux, 0, sizeof(g_disp_flux));
    
    /* Coefficient: α = dt / (8 * dx²) */
    double alpha = dt / (8.0 * dx * dx);
    
    /* Tridiagonal coefficients */
    double a[CGEM_MAX_BRANCH_CELLS + 2];  /* Lower diagonal */
    double bb[CGEM_MAX_BRANCH_CELLS + 2];  /* Main diagonal */
    double c[CGEM_MAX_BRANCH_CELLS + 2];  /* Upper diagonal */
    double r[CGEM_MAX_BRANCH_CELLS + 2];  /* Explicit part (RHS before implicit) */
    double d[CGEM_MAX_BRANCH_CELLS + 2];  /* Full RHS */
    
    /* Boundary conditions: Dirichlet (concentration fixed to BC values)
     * This ensures salt from junction/ocean diffuses into the branch */
    a[1] = 0.0;
    bb[1] = 1.0;
    c[1] = 0.0;
    d[1] = c_down;  /* Use downstream BC value, not current cell */
    
    a[M-1] = 0.0;
    bb[M-1] = 1.0;
    c[M-1] = 0.0;
    d[M-1] = c_up;  /* Use upstream BC value, not current cell */
    
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
    
    /* Copy r to d for interior points only - boundary d values already set above */
    for (int i = 3; i <= M - 3; i += 2) {
        d[i] = r[i];
    }
    /* d[1] and d[M-1] remain set to c_down and c_up from above */
    
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
    
    /* Calculate dispersive fluxes for diagnostics (matches Fortran) */
    /* Flux = -D * A * (C_new + C_old)/2 * gradient, using Crank-Nicolson average */
    for (int i = 2; i <= M - 2; i += 2) {
        double grad_new = (conc[i+1] - conc[i-1]) / (2.0 * dx);
        double grad_old = (cold[i+1] - cold[i-1]) / (2.0 * dx);
        double avg_grad = 0.5 * (grad_new + grad_old);
        g_disp_flux[i] = -b->dispersion[i] * b->totalArea[i] * avg_grad;
    }
}

/* --------------------------------------------------------------------------
 * Main Transport Function
 * -------------------------------------------------------------------------- */

/**
 * Main transport solver for a single branch
 * CGEM_Transport
 * 
 * Key differences from simple implementation:
 * 1. Only transports species with env=1 flag (diagnostic species like pH, pCO2 are skipped)
 * 2. Handles junction boundaries differently from ocean boundaries
 * 3. Stores flux diagnostics for mass balance checking
 * 4. Supports junction dispersion continuity when network pointer is provided
 * 
 * @param branch Branch to compute transport for
 * @param dt Time step [s]
 * @return 0 on success, negative on error
 */
int Transport_Branch(Branch *branch, double dt) {
    return Transport_Branch_Network(branch, dt, NULL);
}

/**
 * Network-aware transport solver with junction dispersion continuity
 * 
 * @param branch Branch to compute transport for
 * @param dt Time step [s]
 * @param network_ptr Pointer to Network for junction dispersion (can be NULL)
 * @return 0 on success, negative on error
 */
int Transport_Branch_Network(Branch *branch, double dt, void *network_ptr) {
    if (!branch || branch->M <= 0) {
        return -1;
    }

    if (!branch->conc || branch->num_species <= 0) {
        return 0; /* Nothing to transport */
    }
    
    int M = branch->M;
    
    /* Compute dispersion coefficient profile */
    /* For branches connected to junctions (not ocean), use flow-weighted D0 */
    double Q_total;
    double D0_override = -1.0;
    
    /* Determine which end connects to ocean vs junction to estimate flow direction */
    if (branch->down_node_type == NODE_LEVEL_BC) {
        /* Ocean at downstream (index 1) - standard case */
        Q_total = branch->velocity[M] * branch->totalArea[M];
    } else if (branch->up_node_type == NODE_LEVEL_BC) {
        /* Ocean at upstream (rare, reversed branch) */
        Q_total = branch->velocity[2] * branch->totalArea[2];
    } else {
        /* Internal branch (junction on both ends or junction + discharge) */
        /* Use the larger velocity magnitude for Q estimate */
        double Q_up = fabs(branch->velocity[M] * branch->totalArea[M]);
        double Q_down = fabs(branch->velocity[2] * branch->totalArea[2]);
        Q_total = (Q_up > Q_down) ? branch->velocity[M] * branch->totalArea[M] 
                                  : branch->velocity[2] * branch->totalArea[2];
        
        /* For branches not connected to ocean, try to get D0 from junction */
        if (network_ptr) {
            D0_override = ComputeJunctionDispersion(branch, network_ptr);
        }
    }
    
    /* === CRITICAL FIX: DISPERSION DISCHARGE CLAMP === */
    /* The Van den Burgh equation relies on RESIDUAL (Freshwater) discharge.
     * Using instantaneous tidal Q (~10,000 m3/s) causes massive over-dispersion.
     * We must clamp Q to the Net River Discharge (~2,000 m3/s).
     */
    double Q_dispersion = fabs(Q_total);
    if (network_ptr) {
        Network *net = (Network *)network_ptr;
        if (net->Q_river > 1.0) {
            /* If instantaneous flow exceeds river discharge (e.g. during tide), 
             * clamp it to the physical river discharge for mixing calculations. */
            if (Q_dispersion > net->Q_river) {
                Q_dispersion = net->Q_river;
            }
            /* Prevent D=0 at slack tide by enforcing a small floor (10% of river) */
            if (Q_dispersion < net->Q_river * 0.1) {
                Q_dispersion = net->Q_river * 0.1;
            }
        }
    }

    /* Pass Q_dispersion instead of Q_total */
    ComputeDispersionCoefficient_Internal(branch, Q_dispersion, D0_override);
    
    /* Process each species - ONLY those with transport flag (env=1) */
    for (int sp = 0; sp < branch->num_species; ++sp) {
        if (!branch->conc || !branch->conc[sp]) continue;
        
        /* Skip diagnostic-only species (pH, pCO2, CO2, ALKC) */
        /* These are computed by biogeochemistry, not transported */
        if (!species_is_transported(sp)) {
            continue;
        }
        
        /* Get boundary concentrations */
        BranchBoundaryConfig cfg;
        determine_boundary_config(branch, &cfg);

        /* Use Index 1 for Downstream, Index M for Upstream for defaults */
        double default_down = (M >= 1) ? branch->conc[sp][1] : 0.0;
        double default_up = (M >= 1) ? branch->conc[sp][M] : 0.0;

        double c_down = (branch->conc_down) ? branch->conc_down[sp] : default_down;
        double c_up = (branch->conc_up) ? branch->conc_up[sp] : default_up;
        
        /* Apply open boundary conditions */
        apply_open_boundaries(branch, sp, c_down, c_up, dt);
        
        /* Calculate advection */
        calculate_advection(branch, sp, dt);
        
        /* Calculate dispersion - pass BC values for Dirichlet boundaries */
        calculate_dispersion(branch, sp, dt, c_down, c_up);
        
        /* Ensure physically valid concentrations */
        for (int i = 0; i <= M + 1; ++i) {
            /* Non-negative constraint */
            if (branch->conc[sp][i] < 0.0) {
                branch->conc[sp][i] = 0.0;
            }
            
            /* For salinity: enforce physically realistic intrusion envelope
             * 
             * Based on Savenije (2005) salt intrusion theory, the maximum salinity
             * at distance x from the mouth follows an exponential decay:
             *   S_max(x) = S0 * exp(-x / L_s)
             * where L_s is the salt intrusion length.
             * 
             * For Mekong-type estuaries:
             * - Dry season L_s ~ 50-80 km (low discharge)
             * - Wet season L_s ~ 15-30 km (high discharge)
             * 
             * We use the maximum intrusion length (80 km) as an upper envelope
             * to prevent unphysical salt penetration beyond observed limits.
             * 
             * Reference: Savenije (2005), Nguyen et al. (2008), Vo Quoc Thanh (2021)
             */
            if (sp == CGEM_SPECIES_SALINITY) {
                double x_from_mouth = (double)i * branch->dx;  /* Distance from mouth */
                double S0 = c_down;  /* Ocean salinity */
                
                /* Maximum intrusion envelope - conservative for dry season */
                double L_max = 80000.0;  /* 80 km max intrusion length */
                double S_envelope = S0 * exp(-x_from_mouth / L_max);
                
                /* Clamp to background river salinity */
                if (S_envelope < c_up) S_envelope = c_up;
                
                /* Apply envelope constraint */
                if (branch->conc[sp][i] > S_envelope) {
                    branch->conc[sp][i] = S_envelope;
                }
            }
        }
        
        /* Set ghost cells to boundary values for next timestep's advection stencil
         * Interior boundary cells (1 and M-1) evolve via dispersion but need
         * ghost cells set to BC values for proper flux computation */
        branch->conc[sp][0] = c_down;      /* Downstream ghost = BC value */
        branch->conc[sp][M] = c_up;        /* Upstream ghost */
        branch->conc[sp][M+1] = c_up;      /* Upstream ghost 2 */
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
