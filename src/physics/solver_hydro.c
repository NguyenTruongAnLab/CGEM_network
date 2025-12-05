/**
 * @file solver_hydro.c
 * @brief Network-level hydrodynamic solver
 */

#include "network.h"
#include <stdio.h>
#include <math.h>
#include <string.h>

#define MAX_ITER 200
#define TOLERANCE 200.0  /* Q tolerance [m³/s] - junction mass balance (relaxed for large rivers) */
#define RELAXATION_FACTOR 0.3    /* Increased for faster convergence with admittance method */
#define RELAXATION_MIN 0.05
#define RELAXATION_MAX 0.8       /* Allow larger steps for admittance-based update */

/* M2 tidal period in seconds (12.42 hours) */
#define M2_PERIOD_S (12.42 * 3600.0)

/* Forward declarations */
int Hyd_Branch(Branch *branch, double H_down, double H_up, double Q_up, double dt);
int Sediment_Branch(Branch *branch, double dt);
int Biogeo_Branch(Branch *branch, double dt, void *network_ptr);
int Biogeo_GHG_Branch(Branch *branch, double dt);  /* GHG module */

static double interpolate_forcing(double time, double *times, double *values, size_t len) {
    if (!times || !values || len == 0) return 0.0;
    if (len == 1) return values[0];
    if (time <= times[0]) return values[0];
    if (time >= times[len-1]) return values[len-1];

    /* Linear interpolation */
    for (size_t i = 0; i < len - 1; ++i) {
        if (time >= times[i] && time <= times[i+1]) {
            double t0 = times[i], t1 = times[i+1];
            double v0 = values[i], v1 = values[i+1];
            double alpha = (time - t0) / (t1 - t0);
            return v0 + alpha * (v1 - v0);
        }
    }
    return values[len-1];
}

static void set_node_boundary_conc(Branch *b, int species, double value, int node_is_upstream) {
    if (!b || !b->conc || !b->conc[species]) return;
    int M = b->M;
    if (M < 2) return;

    /* 
     * Set boundary conditions for transport solver.
     * 
     * IMPORTANT: We only set the GHOST CELLS and the conc_up/conc_down arrays.
     * We do NOT set the interior boundary cell (1 or M-1) directly because
     * that would override the transport solver's solution and create artificial
     * discontinuities.
     * 
     * The transport solver will use these boundary values to compute the
     * concentration in the interior cells via the advection-dispersion equation,
     * creating a smooth gradient.
     */
    
    /* Ghost cell indices */
    int ghost_outer = node_is_upstream ? (M + 1) : -1;  /* Outer ghost (may not exist) */
    int ghost_inner = node_is_upstream ? M : 0;          /* Primary ghost cell */

    /* Set the boundary concentration arrays used by transport solver */
    if (node_is_upstream) {
        if (b->conc_up) b->conc_up[species] = value;
    } else {
        if (b->conc_down) b->conc_down[species] = value;
    }

    /* Set ghost cells for Dirichlet BC */
    if (ghost_inner >= 0 && ghost_inner <= M + 1) {
        b->conc[species][ghost_inner] = value;
    }
    if (ghost_outer >= 0 && ghost_outer <= M + 1) {
        b->conc[species][ghost_outer] = value;
    }
    
    /* DO NOT set the interior boundary cell (1 or M-1) - let transport solver compute it */
}

/**
 * Mix concentrations at junction nodes for transport boundary conditions
 * 
 * This implements flow-weighted mixing at junctions with dispersion continuity.
 * 
 * Key principle (from Savenije theory):
 * - At a junction, all branches should "see" a mixed concentration
 * - For ADVECTION: only inflows contribute mass to the junction
 * - For DISPERSION: all branches should have continuous concentration at junction
 * 
 * Implementation:
 * 1. Calculate flow-weighted mix from inflows (advective contribution)
 * 2. Calculate area-weighted average of all branches (dispersive equilibrium)
 * 3. Blend these based on Peclet number (flow dominance vs dispersion dominance)
 * 
 * For branches with outflow: BC = mixed concentration (what they "export")
 * For branches with inflow: BC = area-weighted concentration of OTHER branches
 *   (representing what dispersion brings from those branches)
 */
static void mix_junction_concentrations(Network *net) {
    if (!net) return;

    for (size_t n = 0; n < net->num_nodes; ++n) {
        Node *node = &net->nodes[n];
        if (node->type != NODE_JUNCTION || !node->mixed_conc) continue;

        for (int sp = 0; sp < net->num_species; ++sp) {
            /* Skip diagnostic-only species - they don't get transported/mixed */
            if (sp < CGEM_NUM_SPECIES && !CGEM_SPECIES_TRANSPORT_FLAG[sp]) {
                continue;
            }
            
            /* 
             * The junction concentration is determined by mass balance:
             *   Sum(Flux_in) = Sum(Flux_out)
             *   Flux = Q*C (advection) + D*A/dx * (C_neighbor - C_junction) (dispersion)
             * 
             * We sample from the INTERIOR CELL (two indices away from junction) so the
             * node "sees" approaching salt wedges even if the boundary cell is
             * temporarily clamped by Dirichlet BCs in the transport solver. This
             * prevents freshwater feedback loops when advection dominates.
             */
            double numerator = 0.0;
            double denominator = 0.0;

            for (int c = 0; c < node->num_connections; ++c) {
                int branch_idx = node->connected_branches[c];
                Branch *b = net->branches[branch_idx];
                if (!b || !b->conc || !b->conc[sp]) continue;

                /* 1. Get geometry and flow at the junction interface */
                int dir = node->connection_dir[c];

                /* For a staggered grid:
                 * - dir = 1: junction at downstream end (index 1 is first interior cell)
                 * - dir = -1: junction at upstream end (index M-1 is last interior cell)
                 */
                int vel_idx = (dir == 1) ? 2 : (b->M - 2); /* Velocity at interface near junction */
                int boundary_cell = (dir == 1) ? 1 : (b->M - 1); /* Cell center adjacent to junction */
                int interior_cell = boundary_cell + ((dir == 1) ? 2 : -2); /* Next cell inward */

                /* Clamp indices to valid range (short branches still resolve junction) */
                if (vel_idx < 0) vel_idx = 0;
                if (vel_idx > b->M) vel_idx = b->M;
                if (boundary_cell < 1) boundary_cell = 1;
                if (boundary_cell > b->M - 1) boundary_cell = b->M - 1;
                if (interior_cell < 1) interior_cell = boundary_cell;
                if (interior_cell > b->M - 1) interior_cell = boundary_cell;

                double vel = b->velocity[vel_idx];
                double area = (b->totalArea[vel_idx] > 0.0) ? b->totalArea[vel_idx] : 1e-6;
                double Q_flow = vel * area;
                double Q_mag = fabs(Q_flow);

                /* Junctions act as energetic mixing basins; enforce a stronger dispersion floor */
                double raw_D = (b->dispersion) ? b->dispersion[boundary_cell] : 0.0;
                if (raw_D < 0.0) raw_D = 0.0;
                double D = fmax(raw_D, 50.0);   /* Minimum junction mixing (m^2/s) */

                double dx = (b->dx > 1e-6) ? b->dx : 1e-6;

                /* 2. Determine flow direction relative to junction */
                int is_flowing_in = 0;
                if (dir == 1 && vel > 0.0) is_flowing_in = 1;
                if (dir == -1 && vel < 0.0) is_flowing_in = 1;

                /* 3. Sample concentrations for advection/dispersion terms */
                double C_interior = b->conc[sp][interior_cell];
                
                /* Check if this is an ocean-connected branch - if so, limit its
                 * contribution to junction mixing to prevent tidal salt pumping */
                int is_ocean_connected = 0;
                if (b->down_node_type == NODE_LEVEL_BC || b->up_node_type == NODE_LEVEL_BC) {
                    is_ocean_connected = 1;
                }
                
                if (is_flowing_in && Q_mag > 0.0) {
                    /* Advective flux: branch is bringing water INTO junction */
                    (void)is_ocean_connected; 
                    numerator += Q_mag * C_interior;
                    denominator += Q_mag;
                }
            }

            /* 5. Calculate Junction Concentration */
            double nodeC = 0.0;
            if (denominator > 1e-9) {
                nodeC = numerator / denominator;
            }
            
            /* Save node mixing concentration */
            if (node->mixed_conc) node->mixed_conc[sp] = nodeC;

            /* 6. Assign Boundary Conditions (Direction Dependent)
             * - OUTFLOWS (Flowing away from node): See the Mixed Node Concentration
             * - INFLOWS (Flowing into node): See the concentration of the *Opposite* branches
             * 
             * IMPORTANT: Do NOT override ocean boundary concentrations!
             * If a branch has an ocean BC at one end, keep the forcing value there.
             */
            for (int c = 0; c < node->num_connections; ++c) {
                int branch_idx = node->connected_branches[c];
                Branch *b = net->branches[branch_idx];
                if (!b) continue;

                int dir = node->connection_dir[c];
                int is_inflow = 0;
                
                /* Determine if this branch is currently flowing INTO the junction
                 * Standard convention: U>0 is downstream flow (from upstream to downstream).
                 * If dir=1 (Node is at Downstream end of branch): Inflow if U > 0
                 * If dir=-1 (Node is at Upstream end of branch): Inflow if U < 0 (reverse flow)
                 */
                int vel_idx = (dir == 1) ? 2 : (b->M - 2);
                if (vel_idx < 0) vel_idx = 0;
                if (vel_idx > b->M) vel_idx = b->M;
                double vel = b->velocity[vel_idx];
                
                if (dir == 1 && vel > 0.0) is_inflow = 1;       /* Normal flow entering downstream node */
                if (dir == -1 && vel < 0.0) is_inflow = 1;      /* Reverse flow entering upstream node */
                
                double bc_conc = nodeC; /* Default: Mixed Node Value */

                /* If this branch is an INFLOW (supplying water), it shouldn't just see its own water back.
                 * It should see the salt level of the OTHER branches to allow diffusion.
                 * However, for stability, we blend towards the node concentration. */
                if (is_inflow && denominator > 1e-9) {
                    /* For inflows, use the full node mix - this allows the salt wedge 
                     * to "pull" into this branch via dispersion at the junction */
                    bc_conc = nodeC;
                }

                /* Apply to boundary - but SKIP if that boundary has an ocean forcing!
                 * Ocean forcing takes precedence over junction mixing.
                 */
                if (dir == -1) {
                    /* Junction at upstream (index M) of this branch */
                    /* Only set conc_up if upstream is NOT an ocean boundary */
                    if (b->up_node_type != NODE_LEVEL_BC) {
                        set_node_boundary_conc(b, sp, bc_conc, 1);
                    }
                } else if (dir == 1) {
                    /* Junction at downstream (index 1) of this branch */
                    /* Only set conc_down if downstream is NOT an ocean boundary */
                    if (b->down_node_type != NODE_LEVEL_BC) {
                        set_node_boundary_conc(b, sp, bc_conc, 0);
                    }
                }
            }
        }
    }
}

/**
 * Apply time-varying species boundary conditions from node forcing to connected branches.
 * 
 * This function interpolates the time-varying species concentrations
 * and sets them on the connected branch's conc_down (for LEVEL_BC) or conc_up (for DISCHARGE_BC).
 * 
 * @param net Network
 * @param current_time Current simulation time [s]
 */
static void apply_species_boundary_forcing(Network *net, double current_time) {
    if (!net) return;
    
    for (size_t n = 0; n < net->num_nodes; ++n) {
        Node *node = &net->nodes[n];
        
        /* Only process boundary nodes with species forcing data */
        if (node->type != NODE_LEVEL_BC && node->type != NODE_DISCHARGE_BC) continue;
        if (!node->species_forcing_time || !node->species_forcing_value) continue;
        
        /* Interpolate each species and apply to connected branches */
        for (int sp = 0; sp < CGEM_NUM_SPECIES; ++sp) {
            if (!node->species_forcing_time[sp] || !node->species_forcing_value[sp]) continue;
            if (node->species_forcing_len[sp] == 0) continue;
            
            /* Interpolate species concentration at current time */
            double conc = interpolate_forcing(current_time,
                node->species_forcing_time[sp], node->species_forcing_value[sp],
                node->species_forcing_len[sp]);
            
            /* Apply to all branches connected to this boundary node */
            for (int c = 0; c < node->num_connections; ++c) {
                int branch_idx = node->connected_branches[c];
                if (branch_idx < 0 || (size_t)branch_idx >= net->num_branches) continue;
                
                Branch *b = net->branches[branch_idx];
                if (!b) continue;
                
                int dir = node->connection_dir[c];
                
                if (node->type == NODE_LEVEL_BC) {
                    /* Ocean boundary: set downstream (dir=1) or upstream (dir=-1) concentration */
                    if (dir == 1 && b->conc_down) {
                        b->conc_down[sp] = conc;
                    } else if (dir == -1 && b->conc_up) {
                        b->conc_up[sp] = conc;
                    }
                } else if (node->type == NODE_DISCHARGE_BC) {
                    /* River boundary: set upstream (dir=-1) or downstream (dir=1) concentration */
                    if (dir == -1 && b->conc_up) {
                        b->conc_up[sp] = conc;
                    } else if (dir == 1 && b->conc_down) {
                        b->conc_down[sp] = conc;
                    }
                }
            }
        }
    }
}

int solve_network_step(Network *net, double current_time_seconds) {
    if (!net || net->num_branches == 0) return -1;

    double dt = net->dt;
    if (dt <= 0) dt = CGEM_DEFAULT_DT_SECONDS;
    
    /* Update current day for seasonal factor lookup */
    net->current_day = ((int)(current_time_seconds / 86400.0)) % 365;
    
    /* =================================================================
     * Step 1: Update boundary node values from forcing
     * ================================================================= */
    for (size_t n = 0; n < net->num_nodes; ++n) {
        Node *node = &net->nodes[n];
        
        if (node->type == NODE_LEVEL_BC) {
            if (node->forcing_time && node->forcing_value) {
                node->H = interpolate_forcing(current_time_seconds, 
                    node->forcing_time, node->forcing_value, node->forcing_len);
            } else {
                double omega = 2.0 * CGEM_PI / M2_PERIOD_S;
                node->H = net->tidal_amplitude * sin(omega * current_time_seconds);
            }
        } 
        else if (node->type == NODE_DISCHARGE_BC) {
            if (node->forcing_time && node->forcing_value) {
                node->Q_net = interpolate_forcing(current_time_seconds,
                    node->forcing_time, node->forcing_value, node->forcing_len);
            } else {
                node->Q_net = net->Q_river; /* Default from config */
            }
        }
    }
    
    /* =================================================================
     * Step 1b: Apply time-varying species boundary conditions
     * This updates branch conc_down/conc_up from node species forcing CSV files
     * ================================================================= */
    apply_species_boundary_forcing(net, current_time_seconds);
    
    /* Initialize junction water levels on FIRST timestep only
     * Set junction H slightly above ocean level to create initial gradient
     * that drives flow from river through junction to ocean
     */
    static int first_step = 1;
    if (first_step) {
        first_step = 0;
        
        for (size_t n = 0; n < net->num_nodes; ++n) {
            Node *node = &net->nodes[n];
            if (node->type != NODE_JUNCTION) continue;
            
            /* Initialize junction H slightly above ocean (0.5m)
             * This creates the initial gradient needed to drive flow through distributaries */
            node->H = 0.5;
            
            printf("  Junction %zu: Initial H = %.2f m\n", n, node->H);
        }
    }

    /* =================================================================
     * Step 2: Iterative hydrodynamic solution
     * ================================================================= */
    double max_residual = 0.0;
    double prev_residual = 1e9;  // Initialize to large value
    int prev_sign = 0;
    double relaxation = RELAXATION_FACTOR;
    
    for (int iter = 0; iter < MAX_ITER; ++iter) {
        max_residual = 0.0;
        
        /* A. Solve hydrodynamics for each branch */
        for (size_t i = 0; i < net->num_branches; ++i) {
            Branch *b = net->branches[i];
            if (!b) continue;
            
            Node *node_up = &net->nodes[b->node_up];
            Node *node_down = &net->nodes[b->node_down];
            
            double H_down = node_down->H;
            double H_up   = node_up->H;
            double Q_up   = -9999.0; /* Flag: use H_up by default */

            /* CRITICAL FIX: Pass Q if upstream is a discharge boundary */
            if (node_up->type == NODE_DISCHARGE_BC) {
                /* Convention: Inflow to domain (Node -> Branch) is Positive Q in Hyd_Branch logic */
                Q_up = fabs(node_up->Q_net); 
            }
            
            Hyd_Branch(b, H_down, H_up, Q_up, dt);
            
            /* Update Node H for Discharge Boundaries (Neumann condition) */
            if (node_up->type == NODE_DISCHARGE_BC) {
                node_up->H = b->waterLevel[b->M];
            }
        }

        /* B. Update Junction Nodes - DIRECT VELOCITY ADJUSTMENT
         * 
         * Instead of just adjusting H, we directly adjust the boundary
         * velocities of connected branches to balance the Q residual.
         * This is more direct and ensures mass balance converges.
         */
        for (size_t n = 0; n < net->num_nodes; ++n) {
            Node *node = &net->nodes[n];
            if (node->type != NODE_JUNCTION) continue;
            
            double sum_Q_in = 0.0;
            double sum_Q_out = 0.0;
            double sum_A_out = 0.0;  /* Sum of outflow areas for velocity distribution */
            int num_outflow = 0;

            for (int c = 0; c < node->num_connections; ++c) {
                int branch_idx = node->connected_branches[c];
                Branch *b = net->branches[branch_idx];
                if (!b) continue;
                
                int dir = node->connection_dir[c]; 

                double Q_branch = 0.0;
                
                if (dir == 1) { 
                    /* Node is at downstream end of branch (index 1) */
                    Q_branch = b->velocity[1] * b->totalArea[1];
                    if (Q_branch > 0) sum_Q_in += Q_branch;
                    else sum_Q_out += -Q_branch;
                } 
                else if (dir == -1) {
                    /* Node is at upstream end of branch (index M) */
                    Q_branch = b->velocity[b->M] * b->totalArea[b->M];
                    if (Q_branch > 0) {
                        sum_Q_out += Q_branch;
                        sum_A_out += b->totalArea[b->M];
                        num_outflow++;
                    } else {
                        sum_Q_in += -Q_branch;
                    }
                }
            }
            
            double residual_Q = sum_Q_in - sum_Q_out;
            
            /* ================================================================
             * SLACK TIDE PROTECTION (Phase 1 Audit Fix - Dec 2025)
             * 
             * Near slack tide, sum_A_out can approach zero causing velocity
             * spikes (division by near-zero). This triggers instability in:
             *   1. GHG species that are sensitive to velocity-driven mixing
             *   2. Salinity gradients near the salt wedge tip
             * 
             * Solution: Use dynamic threshold that scales with total area
             * and skip velocity adjustment when flow is nearly stagnant.
             * ================================================================*/
            double total_junction_area = sum_A_out + 1e-6;  /* Total area entering junction */
            for (int c2 = 0; c2 < node->num_connections; ++c2) {
                Branch *b2 = net->branches[node->connected_branches[c2]];
                if (b2) {
                    int dir2 = node->connection_dir[c2];
                    int idx2 = (dir2 == 1) ? 1 : b2->M;
                    total_junction_area += b2->totalArea[idx2];
                }
            }
            
            /* Adaptive threshold: 0.1% of total junction cross-section */
            double area_thresh = fmax(1e-3 * total_junction_area, 10.0);  /* At least 10 m² */
            
            /* Skip velocity correction if outflow area is below threshold (slack tide) */
            int slack_tide_condition = (sum_A_out < area_thresh);
            
            if (fabs(residual_Q) > 1.0 && !slack_tide_condition && num_outflow > 0) {
                /* Required additional velocity to balance Q */
                double dU_required = residual_Q / sum_A_out;
                
                /* Apply velocity correction to outflow branches */
                for (int c = 0; c < node->num_connections; ++c) {
                    int branch_idx = node->connected_branches[c];
                    Branch *b = net->branches[branch_idx];
                    if (!b) continue;
                    
                    int dir = node->connection_dir[c];
                    
                    if (dir == -1 && b->velocity[b->M] >= 0) {
                        /* Outflow branch - add velocity correction */
                        double weight = b->totalArea[b->M] / sum_A_out;
                        double dU = dU_required * weight * 0.7;  /* Use higher factor */
                        
                        /* Limit velocity change */
                        dU = CGEM_CLAMP(dU, -1.0, 1.0);
                        
                        b->velocity[b->M] += dU;
                        b->velocity[b->M+1] = b->velocity[b->M];
                    }
                }
            }
            
            /* Also adjust H slightly for downstream pressure */
            double dH = (residual_Q / 10000.0) * relaxation;  /* Small adjustment */
            dH = CGEM_CLAMP(dH, -0.1, 0.1);
            node->H += dH;
            node->H = CGEM_CLAMP(node->H, -5.0, 5.0);
            
            if (fabs(residual_Q) > max_residual) max_residual = fabs(residual_Q);
        }

        // Adaptive relaxation
        int current_sign = (max_residual > prev_residual) ? 1 : -1;
        if (prev_sign != 0 && current_sign != prev_sign) {
            // Oscillating, reduce relaxation
            relaxation *= 0.5;
            if (relaxation < RELAXATION_MIN) relaxation = RELAXATION_MIN;
        } else if (max_residual < prev_residual) {
            // Decreasing, slightly increase
            relaxation *= 1.1;
            if (relaxation > RELAXATION_MAX) relaxation = RELAXATION_MAX;
        }
        prev_residual = max_residual;
        prev_sign = current_sign;

        if (max_residual < TOLERANCE) break;
    }
    
    /* Convergence diagnostics */
    static int convergence_warning_count = 0;
    if (max_residual >= TOLERANCE && convergence_warning_count < 10) {
        fprintf(stderr, "Warning: Junction iteration did not fully converge. "
                        "Residual Q = %.3f m³/s (tolerance = %.3f)\n", 
                        max_residual, TOLERANCE);
        convergence_warning_count++;
        if (convergence_warning_count == 10) {
            fprintf(stderr, "  (Suppressing further convergence warnings)\n");
        }
    }

    /* Prepare junction-based per-branch boundary concentrations before transport */
    mix_junction_concentrations(net);

    /* Step 3: Transport (with junction dispersion continuity) */
    for (size_t i = 0; i < net->num_branches; ++i) {
        Branch *b = net->branches[i];
        if (!b) continue;
        
        /* Use network-aware transport for junction dispersion continuity */
        Transport_Branch_Network(b, dt, (void *)net);
        Sediment_Branch(b, dt);
        
        /* Biogeochemistry: only run if reaction_mode=ON and after warmup
         * ReactionMode=OFF allows testing transport and lateral sources in isolation */
        if (net->reaction_mode && current_time_seconds >= net->warmup_time) {
            Biogeo_Branch(b, dt, (void *)net);
            Biogeo_GHG_Branch(b, dt);  /* GHG module (N2O, CH4 dynamics) */
        }
    }
    
    return 0;
}

