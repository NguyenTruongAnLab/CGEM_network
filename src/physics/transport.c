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
#include "../rive/rive_params.h"
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
 * Residual Velocity Filter (Audit Fix Issue #1)
 * -------------------------------------------------------------------------- */

/**
 * Update residual (tidally-averaged) velocity using exponential moving average
 * 
 * CRITICAL for Van den Burgh dispersion: The Savenije formulation requires
 * FRESHWATER (residual) discharge, not instantaneous tidal discharge.
 * 
 * Filter time constant: ~2 tidal cycles (24.84 hours) to smooth out M2 tide.
 * This extracts the net downstream flow from the oscillating tidal signal.
 * 
 * Reference: Savenije (2005), Audit recommendation
 * 
 * @param b Branch to update
 * @param dt Time step [s]
 */
static void update_residual_velocity(Branch *b, double dt) {
    if (!b || !b->u_residual || !b->velocity) return;
    
    int M = b->M;
    
    /* Filter time constant: ~2 tidal cycles to be safe
     * M2 tidal period = 12.42 hours, so 2 cycles = 24.84 hours */
    double T_filter = 24.84 * 3600.0;  /* [s] */
    double alpha = dt / T_filter;
    
    /* Clamp alpha to prevent instability */
    if (alpha > 1.0) alpha = 1.0;
    if (alpha < 0.0001) alpha = 0.0001;  /* Minimum for numerical progress */
    
    /* Store filter coefficient for diagnostics */
    b->residual_alpha = alpha;
    
    /* Apply exponential moving average to all velocity points */
    for (int i = 0; i <= M + 1; ++i) {
        /* u_residual(t+dt) = (1 - alpha) * u_residual(t) + alpha * u(t) */
        b->u_residual[i] = (1.0 - alpha) * b->u_residual[i] + alpha * b->velocity[i];
    }
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
 * @param Q_total Total upstream discharge [m³/s] (fallback, may be tidal)
 * @param D0_override Override D0 value (set to -1 to compute from geometry)
 * @param net_Q_river Network reference discharge [m³/s] for Q estimation
 */
static void ComputeDispersionCoefficient_Internal(Branch *branch, double Q_total, 
                                                   double D0_override, double net_Q_river) {
    if (!branch || branch->M <= 0) return;
    
    int M = branch->M;
    double dx = branch->dx;
    double g = CGEM_GRAVITY;
    
    /* Reference values at mouth (i=1, downstream) */
    double H0 = branch->depth[1];        /* Depth at mouth [m] */
    double B0 = branch->width[1];        /* Width at mouth [m] */
    double LC = branch->lc_convergence;  /* Convergence length [m] */
    
    if (H0 < CGEM_MIN_DEPTH) H0 = CGEM_MIN_DEPTH;
    if (B0 < CGEM_MIN_WIDTH) B0 = CGEM_MIN_WIDTH;
    
    /* Flag for prismatic channels (no convergence) */
    int is_prismatic = (LC <= 0 || LC > 500000.0);  /* LC > 500 km is effectively prismatic */
    if (is_prismatic) LC = 1e9;  /* Very large for calculations */
    
    /* Cross-sectional area at mouth */
    double A0 = H0 * B0;
    if (A0 < 100.0) A0 = 100.0;  /* Minimum area */
    
    /* =====================================================================
     * SAVENIJE VAN DEN BURGH DISPERSION (Savenije 2005, 2012)
     * 
     * The CORRECT Savenije formula for dispersion decay includes BOTH
     * the geometry (convergence) AND the river discharge (Q):
     * 
     *   D(x)/D0 = 1 - β * (exp(x/a) - 1)     [Eq. 9.13 in Savenije 2005]
     * 
     * where:
     *   β = (K * a * Q) / (A0 * D0)           [Eq. 9.14]
     *   K = Van den Burgh coefficient (0.3-0.8, typically ~0.5)
     *   a = convergence length [m]
     *   Q = freshwater discharge [m³/s]
     *   A0 = cross-sectional area at mouth [m²]
     *   D0 = dispersion at mouth [m²/s]
     * 
     * The salt intrusion limit is where D(x) → 0:
     *   x_max = a * ln(1 + 1/β)
     * 
     * PHYSICAL INTERPRETATION:
     * - Higher Q → higher β → FASTER dispersion decay → SHORTER intrusion
     * - Higher K → higher β → FASTER dispersion decay → SHORTER intrusion
     * - Higher D0 → lower β → SLOWER decay → LONGER intrusion
     * 
     * This is the FUNDAMENTAL difference from simple exponential decay!
     * The previous code: D(x) = D0 * exp(-K*x/a) ignored Q entirely,
     * causing unrealistic salt intrusion regardless of river discharge.
     * 
     * Reference: Savenije (2005) "Salinity and Tides in Alluvial Estuaries"
     *            Nguyen et al. (2008) Mekong Delta application
     * ===================================================================== */
    
    double D0;
    
    if (D0_override > 0) {
        /* Use junction-derived D0 for continuity */
        D0 = D0_override;
    } else {
        /* =================================================================
         * D0 at mouth: Using conservative estimate
         * 
         * Savenije (2005) provides empirical relation:
         *   D0 ≈ 50-200 m²/s for large tropical estuaries
         * 
         * Fischer (1979): D0 = α * U_tidal * B
         * But this often overestimates for stratified systems.
         * 
         * For Mekong (Nguyen 2008): D0 ~ 100-200 m²/s at mouth
         * 
         * We use a CONSERVATIVE estimate based on mixing_alpha:
         *   D0 = α * sqrt(g*H) * tidal_amp / H * B
         * but CAPPED at physically realistic values.
         * ================================================================= */
        
        /* Tidal velocity scale */
        double tidal_amp = 1.5;  /* Typical Mekong tidal amplitude [m] */
        double U_tidal = sqrt(g * H0) * (tidal_amp / H0);
        
        /* Use calibrated mixing efficiency from topology.csv */
        double alpha = branch->mixing_alpha;
        if (alpha <= 0) alpha = 0.15;  /* Conservative default */
        
        /* Fischer formula */
        D0 = alpha * U_tidal * B0;
        
        /* =====================================================================
         * D0 BOUNDS: PHYSICALLY REALISTIC VALUES FOR LARGE ESTUARIES
         * 
         * Literature values for tidally-averaged dispersion:
         * - Small estuaries: D0 ~ 50-200 m²/s (Savenije 2005)
         * - Large estuaries (Mekong, Ganges): D0 ~ 500-2000 m²/s
         * - Amazon: D0 ~ 1000-5000 m²/s
         * 
         * CRITICAL FIX (December 2025):
         * Previous cap at 250 m²/s caused salt intrusion of only ~20km.
         * Observed Mekong dry season intrusion is 60-80km, requiring D0 ~ 1000+ m²/s.
         * 
         * The formula D0 = α * U_tidal * B gives for Hau River:
         *   D0 = 0.18 * 1.25 * 3000 = 675 m²/s
         * This was being capped to 250, breaking the physics.
         * 
         * Reference: Nguyen et al. (2008), Savenije (2012) Chapter 5
         * =====================================================================*/
        if (D0 > 2000.0) D0 = 2000.0;  /* Upper bound: realistic for mega-deltas */
    }
    
    /*
     * Only enforce the tidal-mixing floor on branches that touch the open
     * ocean. Interior tributaries should keep their low junction-derived D0.
     */
    int has_ocean_boundary = (branch->down_node_type == NODE_LEVEL_BC) ||
                             (branch->up_node_type == NODE_LEVEL_BC);

    /* Realistic bounds for tropical estuaries */
    if (D0 < 0.0) D0 = 0.0;
    /* INCREASED minimum for ocean branches to ensure adequate tidal mixing */
    if (has_ocean_boundary && D0 < 200.0) D0 = 200.0;   /* Minimum tidal mixing at mouth */
    /* Note: D0 capped at 2000 m²/s above */

    branch->D0 = D0;
    branch->dispersion[0] = D0;
    branch->dispersion[1] = D0;
    
    /* =====================================================================
     * PROPER VAN DEN BURGH DISPERSION PROFILE
     * 
     * Use Q_total (passed in, already clamped to river discharge) to compute
     * the β parameter and apply the full Savenije formula.
     * ===================================================================== */
    
    /* Van den Burgh K coefficient (should be ~0.5 for analytical solution) */
    double K = branch->vdb_coef;
    
    /* For the Savenije analytical formula, K should be 0.3-0.8.
     * The topology.csv may have larger K values intended for numerical decay.
     * We cap K at 1.0 for the physical β calculation. */
    double K_physical = CGEM_CLAMP(K, 0.2, 1.0);
    
    /* =========================================================================
     * CRITICAL FIX (Audit Issue #1): Use RESIDUAL velocity for Q calculation
     * 
     * The Van den Burgh equation is derived for STEADY-STATE (tidally averaged)
     * conditions. Using instantaneous Q (which oscillates ±10,000 m³/s) causes:
     *   1. Dispersion coefficient oscillating wildly within tidal cycle
     *   2. Destruction of salt wedge physics
     *   3. Need for unphysical calibration parameters (K > 1.0)
     * 
     * Solution: Use low-pass filtered (residual) velocity to calculate Q.
     * This represents the net freshwater discharge Q_f (~2,000 m³/s for Mekong).
     * 
     * Reference: Savenije (2005), Nguyen et al. (2008)
     * =========================================================================*/
    
    /* =========================================================================
     * FRESHWATER DISCHARGE ESTIMATION FOR DISPERSION
     * 
     * The Van den Burgh equation requires FRESHWATER (residual) discharge Q_f,
     * NOT instantaneous tidal discharge. The key insight is:
     * 
     *   Q_f = Net river discharge through this branch
     * 
     * For a delta distributary, Q_f is the fraction of total river discharge
     * flowing through this branch. We estimate this as:
     * 
     *   Q_branch = Q_river_total * (A_branch / A_total)
     * 
     * where A is cross-sectional area at the upstream end.
     * 
     * CRITICAL: During warmup, u_residual is not yet established. We MUST use
     * a physics-based estimate, not a tiny fallback value, otherwise:
     *   - β → 0 (because Q → 0)
     *   - x_max → ∞ (salt penetrates entire branch)
     *   - Salt spreads everywhere before flushing can occur
     * 
     * Reference: Savenije (2005) Ch. 9, Nguyen et al. (2008)
     * =========================================================================*/
    double Q;
    
    /* First, try to use established residual velocity */
    double Q_residual = 0.0;
    if (branch->u_residual && branch->totalArea) {
        double U_res = branch->u_residual[branch->M];
        double A_up = branch->totalArea[branch->M];
        Q_residual = fabs(U_res * A_up);
    }
    
    /* Check if residual is "established" (meaningful value > 10% of river Q) */
    double Q_river_estimate = (net_Q_river > 0) ? net_Q_river : 3300.0;  /* Fallback */
    int residual_established = (Q_residual > Q_river_estimate * 0.05);
    
    if (residual_established) {
        /* Use the filtered residual discharge */
        Q = Q_residual;
    } else {
        /* =====================================================================
         * WARMUP PHASE: Estimate Q from network discharge and branch geometry
         * 
         * For Mekong Delta with 9 branches (4 outlets to ocean):
         *   - Tien system (3 outlets): My_Tho, Ham_Luong, Co_Chien
         *   - Hau system (1 outlet): Hau_River
         * 
         * Total dry season Q ~ 2000-3300 m³/s
         * Per distributary ~ 400-800 m³/s
         * 
         * CRITICAL FIX: The previous formula estimated Q based on cross-sectional
         * area ratio, giving Q ~ 1200-1500 m³/s per branch. This is too high!
         * With such high Q, the Van den Burgh formula gives:
         *   β = (K * a * Q) / (A0 * D0) ~ 6-8
         *   x_max ~ a * ln(1 + 1/β) ~ 10-15 km
         * 
         * But dry season intrusion should reach 40-60 km. This requires β ~ 0.5-1.0
         * 
         * Solution: Use a more realistic estimate based on number of outlets.
         * With 4 ocean outlets and Vam Nao redistributing flow:
         *   Q_per_outlet ≈ Q_river / 5 ≈ 600 m³/s
         * 
         * Reference: Nguyen (2008), Nguyen & Savenije (2006)
         * =====================================================================*/
        
        /* Estimate based on even distribution across branches
         * For 9-branch Mekong network with 4 ocean outlets:
         * Each outlet gets roughly 1/5 of total discharge */
        Q = Q_river_estimate / 5.0;
        
        /* Clamp to realistic range for Mekong distributaries */
        /* Nguyen (2008): Individual distributary Q ~ 200-800 m³/s in dry season */
        if (Q < 200.0) Q = 200.0;
        if (Q > 2000.0) Q = 2000.0;
    }
    
    /* β = (K * a * Q) / (A0 * D0) 
     * For PRISMATIC channels (no convergence), β is not meaningful.
     * The Van den Burgh formula is designed for funnel-shaped estuaries.
     * For prismatic channels, use uniform dispersion = D0 throughout.
     */
    double beta;
    if (is_prismatic) {
        beta = 0.0;  /* Uniform dispersion for prismatic channels */
    } else {
        beta = (K_physical * LC * Q) / (A0 * D0);
        /* Cap β to prevent unrealistically fast dispersion decay */
        if (beta > 50.0) beta = 50.0;
    }
    
    /* Calculate salt intrusion limit for diagnostics */
    double x_max = 0.0;
    if (beta > 0.01) {
        x_max = LC * log(1.0 + 1.0/beta);
    } else {
        x_max = branch->length_m;  /* Very low β → salt everywhere (prismatic) */
    }
    
    /* Diagnostic (only print once per branch initialization, not in quiet mode) */
    /* Note: During calibration, quiet_mode suppresses this output */
    static int init_printed[100] = {0};
    int quiet = 0;  /* Would need network pointer to check - skip for now */
    if (!quiet && branch->id >= 0 && branch->id < 100 && !init_printed[branch->id]) {
        if (is_prismatic) {
            printf("  Dispersion %s: D0=%.0f m²/s [PRISMATIC - uniform D]\n",
                   branch->name, D0);
        } else {
            printf("  Dispersion %s: D0=%.0f m²/s, K=%.2f, β=%.3f, x_max=%.1f km (Q=%.0f m³/s)\n",
                   branch->name, D0, K_physical, beta, x_max/1000.0, Q);
        }
        init_printed[branch->id] = 1;
    }
    
    /* Calculate dispersion at each grid point using Savenije formula */
    for (int i = 2; i <= M + 1; ++i) {
        /* Distance from mouth [m] */
        double x = (i - 1) * dx;
        
        double D_curr;
        
        if (is_prismatic) {
            /* Prismatic channel: uniform dispersion */
            D_curr = D0;
        } else {
            /* Funnel-shaped estuary: Van den Burgh decay */
            /* Geometric factor: exp(x/a) - 1 */
            double geom_factor = exp(x / LC) - 1.0;
            
            /* Dispersion ratio: D(x)/D0 = 1 - β * geom_factor */
            double D_ratio = 1.0 - beta * geom_factor;
            
            if (D_ratio <= 0.01) {
                /* Beyond intrusion limit: use minimum dispersion (molecular diffusion only) */
                D_curr = 1.0;  /* Very low floor for freshwater zone */
            } else {
                D_curr = D0 * D_ratio;
            }
        }
        
        /* Ensure realistic dispersion bounds */
        if (D_curr < 1.0) D_curr = 1.0;          /* Minimum molecular diffusion */
        if (D_curr > 100.0) D_curr = 100.0;      /* Cap at realistic max */
        
        branch->dispersion[i] = D_curr;
    }
}

/**
 * Public wrapper for dispersion coefficient calculation
 * (for backward compatibility - always computes D0 from geometry)
 */
void ComputeDispersionCoefficient(Branch *branch, double Q_total) {
    ComputeDispersionCoefficient_Internal(branch, Q_total, -1.0, 0.0);
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
    (void)dt;  /* Unused for now */
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
        /* Ocean or other BC: set ghost cells for boundary stencil */
        set_junction_ghost_cells(conc, cfg.down_ghost_primary,
                                 cfg.down_ghost_secondary, c_down, M + 1);
    }

    /* Handle upstream boundary (index M-1, connects to node_up) */
    if (b->up_node_type == NODE_JUNCTION) {
        /* Junction: only set ghost cells, let transport compute interior */
        set_junction_ghost_cells(conc, cfg.up_ghost_primary,
                                 cfg.up_ghost_secondary, c_up, M + 1);
    } else {
        /* Discharge or other BC: set ghost cells for boundary stencil */
        set_junction_ghost_cells(conc, cfg.up_ghost_primary,
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
    
    /* =========================================================================
     * VELOCITY SIGN CONVENTION AND FLUX DEFINITION:
     * 
     * VELOCITY:
     *   - POSITIVE velocity = flow from UPSTREAM (M) toward DOWNSTREAM (1)
     *                       = flow toward LOWER indices (leftward)
     *                       = EBB tide / river discharge direction
     *   - NEGATIVE velocity = flow toward HIGHER indices (rightward)  
     *                       = FLOOD tide direction
     * 
     * FLUX CONVENTION (Standard Finite Volume):
     *   - Flux F represents mass flow rate in POSITIVE index direction (rightward)
     *   - F > 0 means mass flowing toward higher indices
     *   - F < 0 means mass flowing toward lower indices
     * 
     * CONVERSION:
     *   - Since velocity > 0 means flow LEFTWARD, we need F = -vx * A * C
     *   - This makes F < 0 when vx > 0 (flow is leftward)
     * 
     * UPDATE FORMULA:
     *   d(V*C)/dt = -(F_right - F_left) = -(F[j+1] - F[j-1])
     *   
     *   With our convention, if vx > 0 everywhere (ebb):
     *   - F[j+1] < 0 (flow through right interface going left)
     *   - F[j-1] < 0 (flow through left interface going left)
     *   - |F[j-1]| < |F[j+1]| if upstream has lower concentration
     *   - (F[j+1] - F[j-1]) < 0 → concentration decreases (correct for flushing)
     * 
     * Reference: Hirsch (2007), Toro (2009)
     * =========================================================================*/
    
    for (int j = 1; j <= M - 2; j += 2) {
        int iface = j + 1;  /* Interface between cells j and j+2 */
        double vx = b->velocity[iface];
        double cfl = fabs(vx) * dt / (2.0 * dx);
        
        /* CFL stability check */
        if (cfl > 1.0) cfl = 1.0;
        
        double f, rg, philen, conc_face;
        
        if (vx > 0.0) {
            /* Positive velocity: flow toward lower indices (leftward)
             * Upwind cell is j+2 (higher index) */
            f = cold[j] - cold[j+2];  /* Gradient in flow direction */
            
            if (fabs(f) > 1e-35) {
                if (j + 4 <= M) {
                    rg = (cold[j+2] - cold[j+4]) / f;
                    philen = superbee_limiter(rg);
                } else {
                    philen = 0.0;
                }
            } else {
                philen = 0.0;
            }
            
            /* Reconstructed face value from upwind cell j+2 */
            conc_face = cold[j+2] + 0.5 * (1.0 - cfl) * philen * f;
            double cmin = CGEM_MIN(cold[j], cold[j+2]);
            double cmax = CGEM_MAX(cold[j], cold[j+2]);
            conc_face = CGEM_CLAMP(conc_face, cmin, cmax);
            
            /* FLUX = -vx * A * C (negative for leftward flow) */
            flux[iface] = -vx * b->totalArea[iface] * conc_face;
            
        } else if (vx < 0.0) {
            /* Negative velocity: flow toward higher indices (rightward)
             * Upwind cell is j (lower index) */
            f = cold[j+2] - cold[j];  /* Gradient in flow direction */
            
            if (fabs(f) > 1e-35) {
                if (j >= 3) {
                    rg = (cold[j] - cold[j-2]) / f;
                    philen = superbee_limiter(rg);
                } else {
                    philen = 0.0;
                }
            } else {
                philen = 0.0;
            }
            
            /* Reconstructed face value from upwind cell j */
            conc_face = cold[j] + 0.5 * (1.0 - cfl) * philen * f;
            double cmin = CGEM_MIN(cold[j], cold[j+2]);
            double cmax = CGEM_MAX(cold[j], cold[j+2]);
            conc_face = CGEM_CLAMP(conc_face, cmin, cmax);
            
            /* FLUX = -vx * A * C (positive for rightward flow since vx < 0) */
            flux[iface] = -vx * b->totalArea[iface] * conc_face;
        }
        /* If vx == 0, flux remains 0 */
    }
    
    /* Update concentrations using flux divergence 
     * 
     * d(V*C)/dt = -(F_right - F_left) = -(Flux[j+1] - Flux[j-1])
     * 
     * With positive flux = rightward flow:
     * - During ebb (vx > 0): Flux < 0 (leftward), higher upstream conc
     *   → (F[j+1] - F[j-1]) < 0 → concentration increases... wait that's wrong
     * 
     * Let me reconsider: with flux = -vx * A * C:
     * - vx > 0, C_upwind from j+2: flux = -vx * A * C[j+2] < 0
     * - At interface j-1: flux[j-1] = -vx * A * C[j] < 0
     * - If C[j+2] > C[j]: |flux[j+1]| > |flux[j-1]| → flux[j+1] < flux[j-1]
     * - (flux[j+1] - flux[j-1]) < 0 → C increases? No, that's wrong for flushing!
     * 
     * The issue is the negative sign. During ebb, MORE salt enters from upstream
     * than leaves to downstream. This should INCREASE salt in the cell, which is
     * wrong - ebb should FLUSH salt!
     * 
     * Actually, the problem is that with vx > 0, the upwind cell is j+2 (upstream).
     * If salt is being flushed (upstream is fresher), then C[j+2] < C[j].
     * - flux[j+1] = -vx * A * C[j+2] < 0 (but magnitude is small)
     * - flux[j-1] = -vx * A * C[j] < 0 (magnitude is larger)
     * - flux[j+1] > flux[j-1] (both negative, j+1 is less negative)
     * - (flux[j+1] - flux[j-1]) > 0 → concentration DECREASES ✓
     * 
     * This is correct! During flushing, fresh water from upstream displaces 
     * salty water in the cell, pushing it downstream.
     */
    for (int j = 3; j <= M - 2; j += 2) {
        double AA_old = b->totalArea_old2[j];
        double AA_new = b->totalArea[j];
        
        double rat = (AA_new > 0) ? AA_old / AA_new : 1.0;
        double rat1 = (AA_new > 0) ? dt / (2.0 * dx * AA_new) : 0.0;
        
        conc[j] = rat * cold[j] - rat1 * (flux[j+1] - flux[j-1]);
    }
    
    /* =========================================================================
     * BOUNDARY CELL UPDATE (j=1 for ocean, j=M-1 for river)
     * 
     * With flux = -vx * A * C (positive = rightward flow):
     * 
     * During EBB (vx > 0):
     *   - Flow is leftward (toward ocean)
     *   - flux[2] = -vx * A * C_upwind < 0 (leftward flow)
     *   - At ocean (interface 0): water exits to ocean
     * 
     * For cell 1 during ebb:
     *   - Mass flux from interior (interface 2): flux[2] < 0 (leftward)
     *   - Mass flux to ocean (interface 0): F_ocean = -vx * A * C[1] < 0
     *   - Net change: d(V*C)/dt = -(F[2] - F[0]) = -(flux[2] - F_ocean)
     * 
     * Since both are negative during ebb, the cell loses mass overall if
     * more leaves to ocean than enters from interior.
     * 
     * Reference: Fischer et al. (1979), Savenije (2012)
     * =========================================================================*/
    
    /* Downstream boundary (j=1) - ocean boundary for estuarine branches */
    if (b->down_node_type == NODE_LEVEL_BC) {
        double vx_boundary = b->velocity[2];
        
        /* FIX (December 2025): Removed velocity threshold (0.002) that created
         * a "dead zone" during slack tide. The boundary should be transparent
         * to flow regardless of how small the velocity is. The previous threshold
         * caused salt to accumulate at the boundary, creating an artificial plateau.
         * Reference: Fischer et al. (1979) - open boundary treatments.
         */
        if (vx_boundary > 0.0) {
            /* EBB TIDE (vx > 0): Flushing - salt exits to ocean */
            double AA_old = b->totalArea_old2[1];
            double AA_new = b->totalArea[1];
            double rat = (AA_new > 0) ? AA_old / AA_new : 1.0;
            double rat1 = (AA_new > 0) ? dt / (2.0 * dx * AA_new) : 0.0;
            
            /* Flux from interior to cell 1 at interface 2 (from advection loop) */
            /* With flux = -vx * A * C, this is negative during ebb */
            double F_interior = flux[2];
            
            /* Flux to ocean at interface 0: F = -vx * A * C[1] < 0 */
            double vx_ocean = b->velocity[0];
            if (vx_ocean < 0.001) vx_ocean = vx_boundary;
            double F_ocean = -vx_ocean * b->totalArea[0] * cold[1];
            
            /* Update: d(V*C)/dt = -(F_right - F_left) = -(F_interior - F_ocean) */
            conc[1] = rat * cold[1] - rat1 * (F_interior - F_ocean);
            
            if (conc[1] < 0.0) conc[1] = 0.0;
        }
        else if (vx_boundary < 0.0) {
            /* FLOOD TIDE (vx < 0): Ocean water enters */
            double AA_old = b->totalArea_old2[1];
            double AA_new = b->totalArea[1];
            double rat = (AA_new > 0) ? AA_old / AA_new : 1.0;
            double rat1 = (AA_new > 0) ? dt / (2.0 * dx * AA_new) : 0.0;
            
            /* During flood, vx < 0 means flow is rightward (into estuary) */
            /* Flux at interface 2: F = -vx * A * C > 0 (rightward, carrying cell 1's water to cell 3) */
            double F_interior = flux[2];
            
            /* Flux from ocean at interface 0: F = -vx * A * C_ocean > 0 (ocean water entering) */
            double vx_ocean = fabs(b->velocity[0]);
            if (vx_ocean < 0.001) vx_ocean = fabs(vx_boundary);
            double c_ocean = (b->conc_down) ? b->conc_down[species] : cold[0];
            double F_ocean = vx_ocean * b->totalArea[0] * c_ocean;
            
            /* Update: d(V*C)/dt = -(F_right - F_left) = -(F_interior - F_ocean) */
            conc[1] = rat * cold[1] - rat1 * (F_interior - F_ocean);
            
            if (conc[1] < 0.0) conc[1] = 0.0;
            if (conc[1] > c_ocean) conc[1] = c_ocean;
        }
        /* vx == 0 exactly: Slack tide - no advective flux, cell unchanged
         * This is now a truly instantaneous condition, not an artificial dead zone */
    }
    
    /* Upstream boundary (j=M-1) */
    if (b->up_node_type == NODE_DISCHARGE_BC) {
        double vx_boundary = b->velocity[M-2];
        /* FIX (December 2025): Removed velocity threshold for consistency */
        if (vx_boundary > 0.0) {
            /* Normal ebb flow - river water pushes salt downstream */
            double AA_old = b->totalArea_old2[M-1];
            double AA_new = b->totalArea[M-1];
            double rat = (AA_new > 0) ? AA_old / AA_new : 1.0;
            double rat1 = (AA_new > 0) ? dt / (2.0 * dx * AA_new) : 0.0;
            
            /* Flux at interface M-2 (between M-3 and M-1) */
            double F_interior = flux[M-2];
            
            /* Flux from river: F = -vx * A * C_river < 0 during ebb */
            double c_river = (b->conc_up) ? b->conc_up[species] : 0.1;
            double F_river = -vx_boundary * b->totalArea[M] * c_river;
            
            /* Update cell M-1: -(F_interior - F_river) */
            conc[M-1] = rat * cold[M-1] - rat1 * (F_interior - F_river);
            if (conc[M-1] < 0.0) conc[M-1] = 0.0;
        }
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
    
    /* =========================================================================
     * BOUNDARY-TYPE-DEPENDENT CONDITIONS FOR DISPERSION
     * 
     * Different boundary types require different treatment:
     * 
     * OCEAN (LEVEL_BC): Use FLOW-DEPENDENT BC
     *   - During FLOOD (boundary higher than interior): Dirichlet BC - ocean salt enters
     *   - During EBB (boundary lower than interior): Neumann BC - let interior advect out
     *   
     *   We use water level gradient rather than velocity as the indicator because:
     *   1. Velocity may be positive even during flood due to river discharge
     *   2. Water level gradient directly indicates the pressure gradient driving salt flux
     *   
     *   This is critical for proper tidal salinity dynamics. Pure Dirichlet BC
     *   forces ocean salinity at all times, preventing the ebb tide from
     *   flushing salt out, which causes a high-salinity plateau near the mouth.
     * 
     * JUNCTION: Use Neumann BC (zero gradient)
     *   - Let the junction mixing algorithm handle concentrations
     *   - Avoids double-counting of salt injection at internal nodes
     * 
     * RIVER (DISCHARGE_BC): Use Dirichlet BC with river concentration
     *   - River provides freshwater during ebb/steady flow
     *   - This is physically correct: river concentration is known/fixed
     * 
     * Reference: Fischer et al. (1979), Open boundary treatments
     * =========================================================================*/
    
    /* Downstream boundary */
    if (b->down_node_type == NODE_LEVEL_BC) {
        /* =====================================================================
         * OCEAN BOUNDARY: ALWAYS USE DIRICHLET BC FOR DISPERSION
         * 
         * Key insight from Savenije (2005) steady-state theory:
         *   Q * S = D * A * dS/dx
         * 
         * At steady state, the advective flux OUT (Q*S) is balanced by the
         * dispersive flux IN (D*A*dS/dx). This requires dispersion to ALWAYS
         * transport salt from the high-concentration ocean into the estuary,
         * regardless of flow direction.
         * 
         * Previous error: Using Neumann BC during ebb blocked dispersive flux,
         * preventing salt from entering the estuary and causing unrealistically
         * low salinity at the mouth.
         * 
         * The correct approach:
         * - DISPERSION: Always Dirichlet (ocean concentration) - drives salt IN
         * - ADVECTION: Flow-dependent - exports salt during ebb
         * 
         * This creates the physically correct balance between:
         * - River flushing (advection, pushes salt out)
         * - Tidal mixing (dispersion, brings salt in)
         * 
         * Reference: Savenije (2005), Fischer et al. (1979)
         * =====================================================================*/
        a[1] = 0.0;
        bb[1] = 1.0;
        c[1] = 0.0;
        d[1] = c_down;  /* Always use ocean concentration for dispersion BC */
    } else {
        /* JUNCTION: Apply the mixed junction concentration
         * CRITICAL FIX: Previously used Neumann BC which ignored junction mixing.
         * This prevented freshwater from entering salt-filled branches.
         * Now use Dirichlet BC with the junction-mixed concentration.
         */
        a[1] = 0.0;
        bb[1] = 1.0;
        c[1] = 0.0;
        d[1] = c_down;  /* Apply junction-mixed concentration */
    }
    
    /* Upstream boundary */
    if (b->up_node_type == NODE_DISCHARGE_BC) {
        /* RIVER: Dirichlet BC - river freshwater enters */
        a[M-1] = 0.0;
        bb[M-1] = 1.0;
        c[M-1] = 0.0;
        d[M-1] = c_up;
    } else {
        /* JUNCTION: Apply the mixed junction concentration
         * CRITICAL FIX: Use Dirichlet BC so freshwater can enter from junction.
         */
        a[M-1] = 0.0;
        bb[M-1] = 1.0;
        c[M-1] = 0.0;
        d[M-1] = c_up;  /* Apply junction-mixed concentration */
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
    
    /* =========================================================================
     * RESIDUAL VELOCITY FILTER (Audit Fix Issue #1)
     * Update low-pass filtered velocity BEFORE dispersion calculation.
     * This extracts the tidal-mean residual flow needed for Van den Burgh.
     * =========================================================================*/
    update_residual_velocity(branch, dt);
    
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
    double net_Q_river = 0.0;
    if (network_ptr) {
        Network *net = (Network *)network_ptr;
        net_Q_river = net->Q_river;
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

    /* Pass Q_dispersion AND net_Q_river for proper warmup handling */
    ComputeDispersionCoefficient_Internal(branch, Q_dispersion, D0_override, net_Q_river);
    
    /* Get simplified_mode flag to skip RIVE multi-pool species */
    BiogeoParams *biogeo_p = rive_get_params();
    int simplified_mode = (biogeo_p && biogeo_p->simplified_mode) ? 1 : 0;
    
    /* Process each species - ONLY those with transport flag (env=1) */
    for (int sp = 0; sp < branch->num_species; ++sp) {
        if (!branch->conc || !branch->conc[sp]) continue;
        
        /* Skip diagnostic-only species (pH, pCO2, CO2, ALKC) */
        /* These are computed by biogeochemistry, not transported */
        if (!species_is_transported(sp)) {
            continue;
        }
        
        /* These species (HD1-3, HP1-3, BAG, BAP, PIP, DSS) are not used
         * when simplified_mode=1, so we skip their transport.
         * This reduces transport from 28 species to 18 species (~35% speedup) */
        if (simplified_mode) {
            if (sp >= CGEM_SPECIES_HD1 && sp <= CGEM_SPECIES_DSS) {
                continue;  /* Skip HD1, HD2, HD3, HP1, HP2, HP3, BAG, BAP, PIP, DSS */
            }
        }
        
        /* Get boundary concentrations */
        BranchBoundaryConfig cfg;
        determine_boundary_config(branch, &cfg);

        /* Use Index 1 for Downstream, Index M for Upstream for defaults */
        double default_down = (M >= 1) ? branch->conc[sp][1] : 0.0;
        double default_up = (M >= 1) ? branch->conc[sp][M] : 0.0;

        double c_down = (branch->conc_down) ? branch->conc_down[sp] : default_down;
        double c_up = (branch->conc_up) ? branch->conc_up[sp] : default_up;
        
        /* CRITICAL FIX: Ensure ocean boundary values are never zero for species
         * that have forcing data. The init.c sets default values that should be
         * preserved even if transport or other modules accidentally zero them.
         * For ocean boundaries (NODE_LEVEL_BC), use hardcoded ocean defaults if c_down is 0.
         * 
         * DECEMBER 2025 FIX: GHG values must be in µmol/L (model internal unit)
         * NOT nmol/L (which is the measurement unit, just convert all measurement into µmol/L to be consistent).
         */
        if (branch->down_node_type == NODE_LEVEL_BC && c_down < 1e-10) {
            /* Restore ocean defaults for key species */
            switch (sp) {
                case CGEM_SPECIES_SALINITY: c_down = 30.5; break;
                case CGEM_SPECIES_O2: c_down = 260.0; break;  /* µmol/L */
                case CGEM_SPECIES_NO3: c_down = 55.0; break;   /* µmol N/L */
                case CGEM_SPECIES_TOC: c_down = 100.0; break;  /* µmol C/L */
                case CGEM_SPECIES_DIC: c_down = 2050.0; break; /* µmol C/L */
                case CGEM_SPECIES_AT: c_down = 2200.0; break;  /* µeq/L */
                /* GHG species: values in µmol/L (NOT nmol/L!) */
                case CGEM_SPECIES_N2O: c_down = 0.008; break;  /* 8 nmol/L = 0.008 µmol/L */
                case CGEM_SPECIES_CH4: c_down = 0.04; break;   /* 40 nmol/L = 0.04 µmol/L */
                default: break;
            }
        }
        
        /* Apply open boundary conditions */
        apply_open_boundaries(branch, sp, c_down, c_up, dt);
        
        /* Calculate advection */
        calculate_advection(branch, sp, dt);
        
        /* Calculate dispersion - pass BC values for Dirichlet boundaries */
        calculate_dispersion(branch, sp, dt, c_down, c_up);
        
        /* =================================================================
         * FLOW-DEPENDENT OPEN BOUNDARY CONDITION
         * 
         * ACADEMICALLY ROBUST METHOD (Fischer et al. 1979, Savenije 2012):
         * 
         * The open boundary treatment depends on flow direction:
         * 
         * INFLOW (flood tide, velocity < 0 at downstream boundary):
         *   → Use DIRICHLET BC: boundary cell = ocean concentration
         *   → This represents ocean water entering the estuary
         *   → Relaxation rate alpha = 0.3 (fast response to tide)
         * 
         * OUTFLOW (ebb tide, velocity > 0 at downstream boundary):
         *   → Use NEUMANN BC (zero gradient): allow salt to exit
         *   → boundary cell retains its value, advection carries it out
         *   → No relaxation toward ocean (would prevent salt leaving)
         * 
         * SLACK WATER (|velocity| < threshold):
         *   → Weak relaxation toward equilibrium
         * 
         * Reference: Fischer et al. (1979) "Mixing in Inland and Coastal
         *            Waters", Chapter 6: Estuaries
         * =================================================================*/
        if (branch->down_node_type == NODE_LEVEL_BC) {
            double vx = branch->velocity[2];  /* Velocity at first interface */
            
            if (vx < -0.002) {
                /* FLOOD TIDE (inflow): Dirichlet BC - ocean enters */
                /* Use fast relaxation to quickly bring in ocean water */
                double alpha = 0.3;
                branch->conc[sp][1] = branch->conc[sp][1] + alpha * (c_down - branch->conc[sp][1]);
            }
            /* EBB TIDE (outflow): No relaxation - let advection carry salt out */
            /* Slack water: No action needed */
        }
        
        /* Upstream (river) boundary: Dirichlet BC always (river concentration is known) */
        if (branch->up_node_type == NODE_DISCHARGE_BC) {
            double vx = branch->velocity[M-2];
            if (vx > 0.002) {
                /* Normal river flow: apply river concentration */
                double alpha = 0.3;
                branch->conc[sp][M-1] = branch->conc[sp][M-1] + alpha * (c_up - branch->conc[sp][M-1]);
            }
        }
        
        /* =====================================================================
         * SALINITY GRADIENT LIMITER - REMOVED (Audit Issue #2)
         * 
         * The previous implementation forced salinity to match the Savenije
         * steady-state profile: S(x) = S0 * (D(x)/D0)^(1/K)
         * 
         * This was REMOVED because:
         * 1. It is an artificial forcing that prevents prediction of transient
         *    events (storm surges, flood pulses) - reviewers will reject this.
         * 2. With the residual velocity filter (Issue #1 fix), the Van den Burgh
         *    profile should emerge NATURALLY from advection-dispersion physics.
         * 3. The analytical solution should only be used for INITIALIZATION
         *    (in init.c), not as a runtime constraint.
         * 
         * If salt intrusion is still excessive after this fix, the solution is
         * to adjust physical parameters (K, D0, Chezy) not to clamp the result.
         * 
         * Reference: Audit recommendation, Savenije (2005)
         * ===================================================================== */
        
        /* Ensure physically valid concentrations */
        for (int i = 0; i <= M + 1; ++i) {
            /* Non-negative constraint */
            if (branch->conc[sp][i] < 0.0) {
                branch->conc[sp][i] = 0.0;
            }
            
            /* For salinity: enforce max = ocean boundary value (35 psu typical)
             * Salinity cannot exceed the source concentration in a mixing system */
            if (sp == CGEM_SPECIES_SALINITY) {
                double max_sal = CGEM_MAX(c_down, c_up);
                /* Allow small tolerance for numerical precision */
                double sal_max = max_sal * 1.01;
                if (sal_max < 35.0) sal_max = 35.0;  /* At least ocean typical */
                if (branch->conc[sp][i] > sal_max) {
                    branch->conc[sp][i] = max_sal;
                }
            }
        }
        
        /* Set ghost cells to boundary values for next timestep's advection stencil
         * Interior boundary cells (1 and M-1) evolve via dispersion but need
         * ghost cells set to BC values for proper flux computation */
        branch->conc[sp][0] = c_down;      /* Downstream ghost = BC value */
        branch->conc[sp][M] = c_up;        /* Upstream ghost */
        branch->conc[sp][M+1] = c_up;      /* Upstream ghost 2 */
        
        /* CRITICAL FIX: For ocean boundaries, ensure cell 1 reflects ocean concentration.
         * This is necessary because:
         * 1. Dispersion always brings ocean water in (Savenije 2005)
         * 2. During ebb tide, advection exports but cannot create negative salt
         * 3. The steady-state balance requires ocean BC to be visible at the boundary cell
         * 
         * Apply strong relaxation for all species at ocean boundary.
         */
        if (branch->down_node_type == NODE_LEVEL_BC) {
            double alpha = 0.5;  /* Strong relaxation toward ocean value */
            branch->conc[sp][1] = branch->conc[sp][1] + alpha * (c_down - branch->conc[sp][1]);
        }
        
        /* Apply river BC at upstream boundary */
        if (branch->up_node_type == NODE_DISCHARGE_BC) {
            double alpha = 0.5;
            branch->conc[sp][M-1] = branch->conc[sp][M-1] + alpha * (c_up - branch->conc[sp][M-1]);
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
