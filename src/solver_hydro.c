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

    int cell = node_is_upstream ? (M - 1) : 1;
    if (cell < 0) cell = 0;
    if (cell > M) cell = M;

    int ghost_primary = node_is_upstream ? M : 0;
    int ghost_secondary = node_is_upstream ? (M + 1) : -1;

    if (node_is_upstream) {
        if (b->conc_up) b->conc_up[species] = value;
    } else {
        if (b->conc_down) b->conc_down[species] = value;
    }

    b->conc[species][cell] = value;
    if (ghost_primary >= 0 && ghost_primary <= M + 1) {
        b->conc[species][ghost_primary] = value;
    }
    if (ghost_secondary >= 0 && ghost_secondary <= M + 1) {
        b->conc[species][ghost_secondary] = value;
    }
}

/**
 * Mix concentrations at junction nodes for transport boundary conditions
 * 
 * This implements flow-weighted mixing at junctions:
 *   mixed_conc = Σ(Q_inflow * C_inflow) / Σ(Q_inflow)
 * 
 * For branches receiving flow FROM the junction, we assign the mixed concentration
 * (excluding that branch's own contribution to avoid feedback).
 * 
 * Key concepts:
 * - Inflow: flow entering the junction node
 * - Outflow: flow leaving the junction node
 * - dir=1: junction is at downstream end of branch (index 1)
 * - dir=-1: junction is at upstream end of branch (index M)
 * 
 * Flow sign convention:
 * - Positive velocity = flow from upstream (index M) towards downstream (index 1)
 * - At dir=1: positive Q means flow INTO junction (normal river flow)
 * - At dir=-1: positive Q means flow OUT OF junction (normal distributary flow)
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
            
            double weights[CGEM_MAX_NODE_BRANCHES] = {0};
            double conc_samples[CGEM_MAX_NODE_BRANCHES] = {0};
            double sum_weight = 0.0;
            double sum_weighted_conc = 0.0;
            int flow_dirs[CGEM_MAX_NODE_BRANCHES] = {0}; /* 1=inflow, 0=outflow */

            /* First pass: gather flow and concentrations from all connected branches */
            for (int c = 0; c < node->num_connections; ++c) {
                int branch_idx = node->connected_branches[c];
                Branch *b = net->branches[branch_idx];
                if (!b || !b->conc || !b->conc[sp]) continue;

                int dir = node->connection_dir[c];
                double Q_branch = 0.0;
                double C_branch = 0.0;
                int is_inflow = 0;

                if (dir == 1) {
                    /* Node is at downstream end of branch (index 1)
                     * Sample velocity at index 2 (first interior interface) for more stability */
                    Q_branch = b->velocity[2] * b->totalArea[2];
                    /* Use interior cell for concentration (index 3) to avoid boundary contamination */
                    int idx = (b->M >= 3) ? 3 : 1;
                    C_branch = b->conc[sp][idx];
                    /* Positive velocity means flow towards downstream (INTO junction) */
                    is_inflow = (Q_branch > 0) ? 1 : 0;
                } else if (dir == -1) {
                    /* Node is at upstream end of branch (index M)
                     * Sample velocity at M-2 for stability, or M if unavailable */
                    int vel_idx = (b->M >= 4) ? (b->M - 2) : b->M;
                    Q_branch = b->velocity[vel_idx] * b->totalArea[vel_idx];
                    /* Use interior cell for concentration */
                    int conc_idx = (b->M >= 3) ? (b->M - 2) : (b->M - 1);
                    C_branch = b->conc[sp][conc_idx];
                    /* Positive velocity means flow towards downstream (OUT OF junction) */
                    is_inflow = (Q_branch < 0) ? 1 : 0;
                }

                double weight = fabs(Q_branch);
                weights[c] = weight;
                conc_samples[c] = C_branch;
                flow_dirs[c] = is_inflow;
                
                /* Only count inflows for mixing calculation */
                if (is_inflow && weight > 1.0) {  /* Minimum 1 m³/s to count */
                    sum_weight += weight;
                    sum_weighted_conc += weight * C_branch;
                }
            }

            /* Calculate mixed concentration from all inflows */
            double nodeC = 0.0;
            if (sum_weight > 1.0) {
                nodeC = sum_weighted_conc / sum_weight;
            } else {
                /* No significant inflow detected - use a different strategy */
                /* For tidal systems, use area-weighted average of all branches */
                double sum_area = 0.0;
                double sum_area_conc = 0.0;
                for (int c = 0; c < node->num_connections; ++c) {
                    int branch_idx = node->connected_branches[c];
                    Branch *b = net->branches[branch_idx];
                    if (!b) continue;
                    
                    int dir = node->connection_dir[c];
                    double area = (dir == 1) ? b->totalArea[2] : b->totalArea[b->M - 2];
                    sum_area += area;
                    sum_area_conc += area * conc_samples[c];
                }
                if (sum_area > 0) {
                    nodeC = sum_area_conc / sum_area;
                } else if (node->mixed_conc) {
                    nodeC = node->mixed_conc[sp];  /* Keep previous value */
                }
            }

            /* Save node mixing concentration for reference and output */
            if (node->mixed_conc) node->mixed_conc[sp] = nodeC;

            /* Second pass: assign boundary concentrations to ALL branches at this junction */
            for (int c = 0; c < node->num_connections; ++c) {
                int branch_idx = node->connected_branches[c];
                Branch *b = net->branches[branch_idx];
                if (!b) continue;

                int dir = node->connection_dir[c];
                int is_inflow = flow_dirs[c];
                
                /* For outflow branches: use mixed concentration (excluding own contribution if inflow)
                 * For inflow branches: also update so dispersion can work bidirectionally */
                double nodeC_excl = nodeC;
                if (is_inflow && sum_weight > weights[c] + 1.0) {
                    /* Exclude this branch's contribution from the mix it receives */
                    double remaining = sum_weight - weights[c];
                    if (remaining > 1.0) {
                        nodeC_excl = (sum_weighted_conc - weights[c] * conc_samples[c]) / remaining;
                    }
                }
                
                /* Assign to the appropriate boundary of this branch */
                if (dir == -1) {
                    /* Junction is at upstream end of this branch - set upstream BC */
                    set_node_boundary_conc(b, sp, nodeC_excl, 1);
                } else if (dir == 1) {
                    /* Junction is at downstream end of this branch - set downstream BC */
                    set_node_boundary_conc(b, sp, nodeC_excl, 0);
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

