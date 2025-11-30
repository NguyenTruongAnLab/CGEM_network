/**
 * @file solver_hydro.c
 * @brief Network-level hydrodynamic solver
 */

#include "../network.h"
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
int Biogeo_Branch(Branch *branch, double dt);

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
 * 
 * Reference: Fischer et al. (1979), Savenije (2012) junction mixing theory
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
             * GRADIENT-PRESERVING JUNCTION MIXING
             * 
             * Physical principle (Savenije, 2012):
             * At a junction, all branches must share a common concentration at the
             * interface. The salinity profile in each branch follows the dispersion
             * equation and should be CONTINUOUS across the junction.
             * 
             * The junction concentration is determined by mass balance:
             *   Sum(Flux_in) = Sum(Flux_out)
             *   Flux = Q*C (advection) + D*A/dx * (C_neighbor - C_junction) (dispersion)
             * 
             * We sample from the INTERIOR CELL (two indices away from junction) so the
             * node "sees" approaching salt wedges even if the boundary cell is
             * temporarily clamped by Dirichlet BCs in the transport solver. This
             * prevents freshwater feedback loops when advection dominates.
             */
            double adv_numerator = 0.0;
            double adv_denominator = 0.0;
            double disp_numerator = 0.0;
            double disp_denominator = 0.0;
            double total_Q_in = 0.0;
            double total_G_disp = 0.0;

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

                /* Junction dispersion - use actual branch value
                 * No artificial floor here - let branch dispersion control mixing
                 * Reference: Fischer et al. (1979) Ch. 7 on junction mixing
                 */
                double raw_D = (b->dispersion) ? b->dispersion[boundary_cell] : 0.0;
                if (raw_D < 0.0) raw_D = 0.0;
                double D = fmax(raw_D, 1.0);    /* Minimal turbulence floor (m^2/s) */

                double dx = (b->dx > 1e-6) ? b->dx : 1e-6;

                /* 2. Determine flow direction relative to junction */
                int is_flowing_in = 0;
                if (dir == 1 && vel > 0.0) is_flowing_in = 1;
                if (dir == -1 && vel < 0.0) is_flowing_in = 1;

                /* 3. Sample concentrations */
                double C_interior = b->conc[sp][interior_cell];

                /* 4. Build advective and dispersive contributions
                 * 
                 * Physical principle (Savenije, 2012; Fischer et al., 1979):
                 * - ADVECTION: Only inflows contribute mass to the node
                 * - DISPERSION: Only INFLOWS contribute dispersive mass because:
                 *   For an INFLOW branch, dispersion opposes advection and brings
                 *   higher salinity (from ocean side) toward the junction.
                 *   For an OUTFLOW branch, dispersion would bring mass FROM junction
                 *   INTO the branch, not contribute TO the junction.
                 * 
                 * This prevents the unphysical situation where an ocean-connected
                 * branch (outflow) feeds its high salinity back into a junction
                 * that should be dominated by fresh river inflow.
                 */
                double G_disp = (area * D) / dx;  /* Dispersive conductance */

                /* Advective contribution (inflows only) */
                if (is_flowing_in && Q_mag > 0.0) {
                    adv_numerator += Q_mag * C_interior;
                    adv_denominator += Q_mag;
                    total_Q_in += Q_mag;
                }

                /* Dispersive contribution (ALL branches, not just inflows)
                 * 
                 * CRITICAL FIX: Salt intrusion at junctions
                 * 
                 * Physical basis (Savenije, 2012; Fischer et al., 1979):
                 * Dispersion is a diffusive process that acts in BOTH directions,
                 * independent of flow direction. At a junction node, salt from
                 * downstream branches (estuaries) must be able to diffuse into
                 * the junction even when those branches have net outflow.
                 * 
                 * Without this, a junction receiving fresh river water will never
                 * "see" the salt in downstream estuaries, creating an artificial
                 * freshwater cap at every junction.
                 * 
                 * The dispersive mixing uses area-weighted concentrations from
                 * ALL connected branches, allowing salt to penetrate upstream
                 * through the junction via the diffusion term.
                 */
                if (G_disp > 0.0) {
                    disp_numerator += G_disp * C_interior;
                    disp_denominator += G_disp;
                    total_G_disp += G_disp;
                }
            }

            /* 5. Calculate Junction Concentration using PECLET-WEIGHTED mixing
             * 
             * Physical basis (Savenije, 2012, Chapter 6; Fischer et al., 1979):
             * 
             * At a junction, both advection and dispersion contribute to mixing.
             * The relative importance is determined by the Peclet number:
             *   Pe = U*L / D = (Q/A) * L / D
             * 
             * For Pe >> 1: Advection dominates → use flow-weighted average
             * For Pe << 1: Dispersion dominates → use area-weighted diffusive average
             * 
             * SALT INTRUSION FIX:
             * Without dispersive contribution, junctions act as freshwater caps,
             * preventing salt from penetrating upstream. The blending ensures
             * that salt from downstream estuaries can diffuse into the junction
             * even when river flow is dominant.
             * 
             * Blending formula:
             *   C_node = (Pe * C_adv + C_disp) / (Pe + 1)
             * where Pe is normalized to [0, 10] range
             */
            double nodeC = 0.0;
            double C_adv = 0.0;
            double C_disp = 0.0;
            
            /* Calculate advective concentration */
            if (adv_denominator > 1e-9) {
                C_adv = adv_numerator / adv_denominator;
            }
            
            /* Calculate dispersive concentration (from ALL branches) */
            if (disp_denominator > 1e-9) {
                C_disp = disp_numerator / disp_denominator;
            }
            
            /* Peclet number for blending */
            double Pe = 0.0;
            if (total_G_disp > 1e-9) {
                Pe = total_Q_in / total_G_disp;  /* Dimensionless ratio */
                if (Pe > 10.0) Pe = 10.0;        /* Cap for numerical stability */
            }
            
            /* Peclet-weighted blend */
            if (adv_denominator > 1e-9 && disp_denominator > 1e-9) {
                /* Both advection and dispersion active - blend */
                nodeC = (Pe * C_adv + C_disp) / (Pe + 1.0);
            } else if (adv_denominator > 1e-9) {
                /* Only advection - use advective */
                nodeC = C_adv;
            } else if (disp_denominator > 1e-9) {
                /* Only dispersion - use dispersive */
                nodeC = C_disp;
            }
            
            /* Save node mixing concentration */
            if (node->mixed_conc) node->mixed_conc[sp] = nodeC;

            /* 6. Assign Boundary Conditions (Direction Dependent)
             * - OUTFLOWS (Flowing away from node): See the Mixed Node Concentration
             * - INFLOWS (Flowing into node): See the concentration of the *Opposite* branches
             * (This allows salt diffusion upstream against the flow)
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
                
                /* Use the Peclet-weighted node concentration for all branches */
                double bc_conc = nodeC;

                /* Apply to boundary */
                if (dir == -1) {
                    /* Junction at upstream (index M) of this branch */
                    set_node_boundary_conc(b, sp, bc_conc, 1);
                } else if (dir == 1) {
                    /* Junction at downstream (index 1) of this branch */
                    set_node_boundary_conc(b, sp, bc_conc, 0);
                }
            }
        }
    }
}

int solve_network_step(Network *net, double current_time_seconds) {
    if (!net || net->num_branches == 0) return -1;

    double dt = net->dt;
    if (dt <= 0) dt = CGEM_DEFAULT_DT_SECONDS;
    
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
     * CRITICAL FIX: Update Q_river for dispersion calculations
     * =================================================================
     * Savenije's Van den Burgh equation: L_d = D0 * A / (K * Q_f)
     * 
     * Q_f MUST be the CURRENT freshwater discharge, not the static config
     * value. During dry season, Q drops from 2000 → 1000 m³/s, which
     * DOUBLES the dispersion length scale L_d and allows salt to penetrate
     * further upstream.
     * 
     * Without this fix, the model "thinks" Q=2000 even when forcing=1000,
     * artificially suppressing salt intrusion during dry season.
     * 
     * Reference: Savenije (2005, 2012) Eq. 9.30
     */
    for (size_t n = 0; n < net->num_nodes; ++n) {
        Node *node = &net->nodes[n];
        if (node->type == NODE_DISCHARGE_BC && node->Q_net > 0) {
            /* Update network Q_river to match actual upstream forcing */
            net->Q_river = fabs(node->Q_net);
            break;  /* Use first upstream discharge node */
        }
    }
    
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
         * 
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
            
            /* If there's a Q imbalance, distribute it to outflow branches */
            if (fabs(residual_Q) > 1.0 && sum_A_out > 0.0 && num_outflow > 0) {
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
    
    /* Convergence diagnostics (similar to Fortran Update_otherhydrodynamic_arrays) */
    /* Log warning if not converged - useful for debugging */
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
        
        /* Update BCs for transport */
        /* The branch-level boundary concentrations are set centrally by
         * mix_junction_concentrations() for junction nodes. Do not override
         * branch conc_up here. For discharge/level BCs, initialization/forcing
         * is left unchanged.
         */
        
        /* Use network-aware transport for junction dispersion continuity */
        Transport_Branch_Network(b, dt, (void *)net);
        Sediment_Branch(b, dt);
        
        if (current_time_seconds >= net->warmup_time) {
            Biogeo_Branch(b, dt);
        }
    }
    
    return 0;
}

