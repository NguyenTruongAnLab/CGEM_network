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
 * - **LOW-PASS VELOCITY FILTER for residual discharge** (Audit fix)
 * 
 * Grid convention (matching Fortran):
 * - Concentrations defined at odd indices (cell centers)
 * - Fluxes calculated at even indices (interfaces)
 * 
 * Reference: Savenije (2012), Volta et al. (2014), Gisen et al. (2015)
 */

#include "../network.h"
#include "../define.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Local flux storage for diagnostics (per-species, reused each timestep) */
static double g_adv_flux[CGEM_MAX_BRANCH_CELLS + 2];
static double g_disp_flux[CGEM_MAX_BRANCH_CELLS + 2];

/* ============================================================================
 * RESIDUAL VELOCITY LOW-PASS FILTER
 * 
 * The Van den Burgh dispersion equation requires RESIDUAL (freshwater) discharge,
 * NOT instantaneous tidal discharge. In branching deltas like the Mekong, using
 * instantaneous Q leads to massive over-dispersion during peak tidal flow.
 * 
 * This filter extracts the tidally-averaged velocity from the instantaneous signal:
 *   u_residual[new] = (1 - alpha) * u_residual[old] + alpha * u_instantaneous
 * 
 * With alpha = dt / T_filter, where T_filter ~ 25 hours (>M2 period)
 * 
 * Reference: Savenije (2005), Audit recommendation for "Residual Discharge Problem"
 * ============================================================================*/

/**
 * Update the residual velocity array using exponential low-pass filter
 * 
 * @param branch Branch to update
 * @param dt Current timestep [s]
 */
void UpdateResidualVelocity(Branch *branch, double dt) {
    if (!branch || !branch->u_residual || branch->M <= 0) return;
    
    /* Filter time constant: ~ 25 hours (longer than M2 period of 12.42 hours) */
    /* This ensures we capture the mean flow, not tidal oscillations */
    double T_filter = 25.0 * 3600.0;  /* 25 hours in seconds */
    double alpha = dt / T_filter;
    
    /* Override with branch-specific value if set */
    if (branch->residual_alpha > 0.0 && branch->residual_alpha < 1.0) {
        alpha = branch->residual_alpha;
    }
    
    /* Clamp alpha to reasonable range */
    if (alpha > 0.1) alpha = 0.1;   /* Max 10% per step (prevents instability) */
    if (alpha < 1e-6) alpha = 1e-6;
    
    int M = branch->M;
    for (int i = 0; i <= M + 1; ++i) {
        branch->u_residual[i] = (1.0 - alpha) * branch->u_residual[i] + 
                                 alpha * branch->velocity[i];
    }
}

/**
 * Get residual discharge at a branch cross-section
 * Uses the filtered velocity if available, otherwise falls back to instantaneous
 * 
 * @param branch Branch
 * @param idx Grid index
 * @return Residual discharge [m³/s]
 */
double GetResidualDischarge(Branch *branch, int idx) {
    if (!branch || idx < 0 || idx > branch->M + 1) return 0.0;
    
    double u_res = branch->velocity[idx];  /* Default to instantaneous */
    
    if (branch->u_residual && fabs(branch->u_residual[idx]) > 1e-10) {
        u_res = branch->u_residual[idx];
    }
    
    double area = (branch->totalArea) ? branch->totalArea[idx] : 1.0;
    return u_res * area;
}

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
         * - Mekong: 150-400 m²/s (Nguyen et al., 2008)
         * 
         * Using calibrated empirical relationship based on Canter-Cremers number:
         *   D0 = k * sqrt(N * g) * H^1.5
         * 
         * LITERATURE-BASED CALIBRATION:
         * 
         * The Canter-Cremers number N relates tidal mixing to geometry:
         *   N = π * Q_f / (H0 * B0)
         * 
         * For the Mekong Delta (Nguyen et al., 2008):
         * - Q_f ~ 3000-3500 m³/s (dry season)
         * - H0 ~ 12-15 m
         * - B0 ~ 2000-3000 m
         * - Observed D0 ~ 200-400 m²/s
         * 
         * Empirical coefficient k calibration:
         * - k = 8-12: Causes D0 ~ 500-800 m²/s → excessive dispersion
         * - k = 4-6:  Gives D0 ~ 200-400 m²/s → matches observations
         * 
         * Using k = 5.0 to achieve target D0 ~ 300 m²/s
         * 
         * Reference: Savenije (2005) "Salinity and Tides in Alluvial Estuaries"
         *            Nguyen et al. (2008) Mekong salt intrusion study
         */
        /* Canter-Cremers number N = (Pi * Q) / (H0 * B0)
         * This is an estuary shape parameter, NOT scaled by vdb_coef.
         * vdb_coef (K) controls dispersion DECAY, not mouth dispersion.
         */
        double N = CGEM_PI * fabs(Q_total) / (H0 * B0);
        
        /* Savenije (2005) Canter-Cremers formulation:
         *   D0 = k * sqrt(N * g) * H0^1.5
         * 
         * SCIENTIFIC BASIS:
         * The empirical coefficient k depends on estuary type:
         * - k = 14-26: Classic shallow estuaries (Savenije, 2005)
         * - k = 6-10:  Deep alluvial estuaries (H > 10m)
         * - k = 4-6:   Large tropical deltas (Mekong, Ganges)
         * 
         * For MEKONG distributaries (H = 11-15m):
         * - Target D0 = 200-400 m²/s at mouth (from Nguyen et al., 2008)
         * - Using k = 5.0 to achieve D0 ~ 250-350 m²/s
         * - This allows proper gradient formation within branch length
         * 
         * CALIBRATION NOTE (2024-12):
         * Reduced from k=8 to k=5 based on literature values.
         * Higher D0 causes dispersion to dominate over advection, creating
         * nearly constant salinity until hitting the junction boundary.
         * 
         * Reference: Savenije (2005) Table 9.2, Nguyen et al. (2008)
         */
        D0 = 5.0 * sqrt(N * g) * pow(H0, 1.5);
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
    
    /* Calculate dispersion along the branch using Savenije (2012) formulation
     * 
     * SAVENIJE EXPONENTIAL DECAY MODEL:
     * The key insight from Savenije (2005, 2012) is that in a well-mixed 
     * alluvial estuary, both salinity AND dispersion decay exponentially
     * from the mouth upstream:
     * 
     *   D(x) = D0 * exp(-x / a)
     *   S(x) = S0 * exp(-x / L)
     * 
     * where:
     *   a = convergence length (width decay scale)
     *   L = salt intrusion length (depends on a, Q, H)
     *   D0 = dispersion at mouth
     *   x = distance from MOUTH (not from upstream!)
     * 
     * The Van den Burgh K coefficient relates D decay to S decay:
     *   D(x) = D0 * exp(-K * x / a)
     * 
     * CRITICAL: x is measured from the OCEAN BOUNDARY (mouth), not upstream.
     * In our grid: x_from_mouth = branch->length - (i * dx)
     * 
     * Reference: Savenije (2005) Chapter 4, Eq. 4.22-4.25
     */
    double a = branch->lc_convergence;  /* Convergence length [m] */
    if (a < 1000.0 || a > 1e9) {
        /* No convergence (prismatic) or invalid - use branch length */
        a = branch->length_m;
    }
    
    double K = branch->vdb_coef;
    if (K <= 0.0) {
        K = 0.35;  /* Default to alluvial estuary mid-range (Savenije, 2005) */
    }
    /* Literature range: K ≈0.2–2.0 for convergent estuaries (Savenije, 2012; Gisen et al., 2015) */
    K = CGEM_CLAMP(K, 0.1, 2.0);
    
    /* Compute dispersion decay scale: L_D = a / K
     * This is the e-folding length for dispersion decay.
     * Lower K = longer L_D = slower decay = longer salt intrusion
     */
    double L_D = a / K;
    L_D = CGEM_CLAMP(L_D, 5000.0, 500000.0);  /* 5-500 km */
    
    for (int i = 2; i <= M + 1; ++i) {
        /* Distance from MOUTH (ocean boundary) */
        double x_from_mouth;
        
        if (branch->down_node_type == NODE_LEVEL_BC) {
            /* Ocean at downstream (normal estuary orientation) */
            /* Index 1 is at mouth, index M is upstream */
            x_from_mouth = (i - 1) * dx;
        } else if (branch->up_node_type == NODE_LEVEL_BC) {
            /* Ocean at upstream (reversed) */
            x_from_mouth = (M - i) * dx;
        } else {
            /* Interior branch - use distance from downstream end */
            x_from_mouth = (i - 1) * dx;
        }
        
        /* Savenije exponential decay: D(x) = D0 * exp(-x / L_D) */
        double decay_factor = exp(-x_from_mouth / L_D);
        double D_curr = D0 * decay_factor;
        
        /* Ensure realistic dispersion bounds
         * 
         * Physical basis (Fischer et al., 1979; Savenije, 2012):
         * - Tidal mixing at mouth: D0 ~ 100-1000 m²/s
         * - River turbulent diffusion upstream: ~ 1-10 m²/s
         * - The exponential decay ensures gradual transition
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
    
    /* Determine if this is an ocean-connected branch with junction upstream
     * This requires special treatment: Neumann BC at upstream instead of Dirichlet
     * 
     * Physical basis (Savenije 2012):
     * For estuaries, the salinity profile should decay exponentially from ocean.
     * At the upstream end (junction), forcing Dirichlet BC to junction concentration
     * creates an artificial cliff. Using Neumann (zero-gradient) allows natural decay.
     */
    int use_neumann_upstream = (b->down_node_type == NODE_LEVEL_BC && 
                                 b->up_node_type != NODE_LEVEL_BC &&
                                 b->up_node_type != NODE_DISCHARGE_BC);
    
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
    
    /* =========================================================================
     * DOWNSTREAM BOUNDARY (index 1)
     * Always use Dirichlet BC - concentration fixed to ocean/junction value
     * ========================================================================= */
    a[1] = 0.0;
    bb[1] = 1.0;
    c[1] = 0.0;
    d[1] = c_down;  /* Use downstream BC value */
    
    /* =========================================================================
     * UPSTREAM BOUNDARY (index M-1)
     * Use Dirichlet for river inlets, Neumann for ocean-connected branches
     * ========================================================================= */
    if (use_neumann_upstream) {
        /* NEUMANN BC: Zero-gradient condition for ocean-connected branches
         * 
         * Physical basis (Savenije 2012, Fischer et al. 1979):
         * The steady-state salinity profile S(x) = S0 * exp(-x/L) should
         * decay smoothly toward the upstream end. Using Neumann BC allows
         * the exponential profile to extend naturally without artificial cliffs.
         * 
         * Implementation: Use a "soft" Neumann by setting the upstream cell
         * to the interior extrapolation. This is more stable than modifying
         * the tridiagonal system directly.
         */
        a[M-1] = 0.0;
        bb[M-1] = 1.0;
        c[M-1] = 0.0;
        /* Use interior value as temporary - will be overridden after solve */
        d[M-1] = (M >= 4) ? cold[M-3] : cold[1];
    } else {
        /* DIRICHLET BC: Fixed concentration for river inlets and interior branches */
        a[M-1] = 0.0;
        bb[M-1] = 1.0;
        c[M-1] = 0.0;
        d[M-1] = c_up;  /* Use upstream BC value */
    }
    
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
    
    /* POST-PROCESSING: Enforce Neumann BC for ocean-connected branches
     * 
     * After the tridiagonal solve, override the upstream cell with an
     * extrapolation from interior. This ensures smooth salinity decay.
     */
    if (use_neumann_upstream && M >= 4) {
        /* Zero-gradient: C[M-1] = C[M-3] */
        conc[M-1] = conc[M-3];
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
    
    /* =======================================================================
     * CRITICAL FIX: Use RESIDUAL DISCHARGE for Van den Burgh dispersion
     * =======================================================================
     * The Van den Burgh equation: L_d = D0 * A / (K * Q_f)
     * 
     * Q_f MUST be the RESIDUAL (freshwater) discharge, NOT instantaneous tidal Q.
     * 
     * Problem: At slack tide, instantaneous Q → 0, causing L_d → ∞
     * This makes dispersion constant throughout the branch, allowing salt
     * to "teleport" upstream during every slack tide (4x per day).
     * 
     * Solution: Use the low-pass filtered residual velocity (u_residual)
     * which extracts the net river flow component (~25-hour averaging).
     * 
     * Reference: Savenije (2005, 2012) explicitly states Q must be residual.
     * ======================================================================= */
    
    /* Method 1: Use local residual discharge from filter */
    double Q_residual = fabs(GetResidualDischarge(branch, M));
    
    /* Method 2: Fallback to network-wide Q_river if filter hasn't converged */
    double net_Q_river = 100.0;  /* Safety default */
    if (network_ptr) {
        Network *net = (Network *)network_ptr;
        if (net->Q_river > 10.0) {
            net_Q_river = net->Q_river;
        }
    }
    
    /* Select dispersion discharge: prefer local residual, fallback to global */
    double Q_dispersion;
    if (Q_residual > 10.0) {
        /* Local residual is available and significant */
        Q_dispersion = Q_residual;
    } else {
        /* Fallback to network-wide Q_river (robust for early timesteps) */
        Q_dispersion = net_Q_river;
    }
    
    /* Compute D0 override for junction-connected branches */
    double D0_override = -1.0;
    if (branch->down_node_type != NODE_LEVEL_BC && 
        branch->up_node_type != NODE_LEVEL_BC && 
        network_ptr) {
        D0_override = ComputeJunctionDispersion(branch, network_ptr);
    }
    
    /* Pass Q_dispersion (residual) to dispersion calculation */
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
            
            /* Physical bounds for conservative tracers (salinity)
             * 
             * For conservative species, concentration cannot exceed the maximum
             * of the boundary values (ocean and river). This is a physical
             * constraint, not an artificial limit.
             * 
             * Numerical schemes (especially Crank-Nicolson) can produce small
             * overshoots; this correction maintains physical consistency.
             */
            if (sp == CGEM_SPECIES_SALINITY) {
                double max_bc = (c_down > c_up) ? c_down : c_up;
                double min_bc = (c_down < c_up) ? c_down : c_up;
                if (branch->conc[sp][i] > max_bc) {
                    branch->conc[sp][i] = max_bc;
                }
                if (branch->conc[sp][i] < min_bc) {
                    branch->conc[sp][i] = min_bc;
                }
            }
        }
        
        /* Set boundary cells conditionally based on boundary type
         * 
         * CRITICAL FIX for Junction Salt Intrusion:
         * 
         * For OCEAN/DISCHARGE boundaries (Dirichlet): Force both ghost AND interior cell
         * For JUNCTION boundaries (Neumann-like): Force ONLY ghost cells
         * 
         * This allows the transport solver to compute natural diffusion gradients
         * at junctions, preventing the artificial "salt cliff" effect where
         * interior cells are forced to freshwater values.
         * 
         * Reference: Fischer et al. (1979) - Junction mixing theory
         */
        
        /* DOWNSTREAM BOUNDARY (index 0, 1) */
        if (branch->down_node_type == NODE_LEVEL_BC) {
            /* Ocean: Full Dirichlet - force ghost AND interior */
            branch->conc[sp][0] = c_down;  /* Ghost cell */
            branch->conc[sp][1] = c_down;  /* Interior boundary cell */
        } else {
            /* Junction: Only ghost cell - let interior evolve via transport */
            branch->conc[sp][0] = c_down;  /* Ghost cell only */
            /* conc[1] is left alone - computed by advection-dispersion */
        }
        
        /* UPSTREAM BOUNDARY (index M, M+1) */
        if (branch->up_node_type == NODE_DISCHARGE_BC) {
            /* River inlet: Full Dirichlet - force ghost AND interior */
            branch->conc[sp][M] = c_up;    /* Ghost cell */
            branch->conc[sp][M+1] = c_up;  /* Outer ghost */
            /* Also force the last interior cell for discharge BC */
            if (M >= 2) branch->conc[sp][M-1] = c_up;
        } else if (branch->down_node_type == NODE_LEVEL_BC) {
            /* OCEAN-CONNECTED BRANCH with JUNCTION UPSTREAM
             * 
             * This is a critical case for salt intrusion physics.
             * 
             * Physical situation (Co_Chien, Ham_Luong, My_Tho, Hau_River):
             * - Ocean at downstream (index 1) with high salinity
             * - Junction at upstream (index M) with mostly fresh water
             * 
             * PROBLEM with Dirichlet BC at junction:
             * - Setting ghost cells to fresh junction concentration creates
             *   an artificial "salt cliff" at the upstream end
             * - The exponential salinity profile cannot extend to its
             *   natural intrusion length
             * 
             * SOLUTION: Use NEUMANN (zero-gradient) boundary condition
             * - Ghost cells extrapolate from interior, not forced to junction value
             * - Allows salt profile to naturally decay toward the junction
             * - Junction still influences branch via advective flux
             * 
             * Mathematical basis (Savenije, 2012):
             * The steady-state salinity profile S(x) = S0 * exp(-x/L) should
             * extend continuously. A Neumann BC at the upstream end allows
             * dS/dx to remain finite (exponential decay) rather than forcing
             * a discontinuous jump.
             * 
             * Reference: Fischer et al. (1979) Ch. 7 - Open boundary conditions
             *            Savenije (2012) Ch. 8 - Delta network salt intrusion
             */
            if (M >= 2) {
                /* Extrapolate from interior (zero-gradient approximation) */
                double dS_dx = (branch->conc[sp][M-1] - branch->conc[sp][M-2]);
                double S_extrap = branch->conc[sp][M-1] + dS_dx;
                
                /* Clamp to physical bounds: cannot exceed ocean or go below river */
                double S_ocean = c_down;
                double S_river = 0.1;  /* Minimum freshwater */
                S_extrap = CGEM_CLAMP(S_extrap, S_river, S_ocean);
                
                branch->conc[sp][M] = S_extrap;      /* Ghost cell - extrapolated */
                branch->conc[sp][M+1] = S_extrap;    /* Outer ghost */
            } else {
                /* Short branch - use interior value */
                branch->conc[sp][M] = branch->conc[sp][M-1];
                branch->conc[sp][M+1] = branch->conc[sp][M-1];
            }
        } else {
            /* Interior branch (junction on both ends)
             * Use junction concentration but allow interior to evolve
             */
            branch->conc[sp][M] = c_up;    /* Ghost cell */
            branch->conc[sp][M+1] = c_up;  /* Outer ghost */
            /* conc[M-1] is left alone - computed by advection-dispersion */
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
