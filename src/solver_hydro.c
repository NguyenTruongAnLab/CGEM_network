/**
 * @file solver_hydro.c
 * @brief Network-level hydrodynamic solver
 */

#include "network.h"
#include <stdio.h>
#include <math.h>
#include <string.h>

#define MAX_ITER 100
#define TOLERANCE 1e-3
#define RELAXATION_FACTOR 0.05
#define RELAXATION_MIN 0.01
#define RELAXATION_MAX 0.2

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

static void mix_junction_concentrations(Network *net) {
    if (!net) return;

    for (size_t n = 0; n < net->num_nodes; ++n) {
        Node *node = &net->nodes[n];
        if (node->type != NODE_JUNCTION || !node->mixed_conc) continue;

        for (int sp = 0; sp < net->num_species; ++sp) {
            double weights[CGEM_MAX_NODE_BRANCHES] = {0};
            double conc_samples[CGEM_MAX_NODE_BRANCHES] = {0};
            double sum_weight = 0.0;
            double sum_weighted_conc = 0.0;

            for (int c = 0; c < node->num_connections; ++c) {
                int branch_idx = node->connected_branches[c];
                Branch *b = net->branches[branch_idx];
                if (!b || !b->conc || !b->conc[sp]) continue;

                int dir = node->connection_dir[c];
                double Q_branch = 0.0;
                double C_branch = 0.0;

                if (dir == 1) {
                    Q_branch = b->velocity[1] * b->totalArea[1];
                    int idx = (b->M >= 3) ? 2 : 1;
                    C_branch = b->conc[sp][idx];
                } else if (dir == -1) {
                    Q_branch = b->velocity[b->M] * b->totalArea[b->M];
                    int idx = (b->M >= 3) ? (b->M - 2) : b->M - 1;
                    C_branch = b->conc[sp][idx];
                }

                double weight = fabs(Q_branch);
                weights[c] = weight;
                conc_samples[c] = C_branch;
                if (weight > 0.0) {
                    sum_weight += weight;
                    sum_weighted_conc += weight * C_branch;
                }
            }

            double nodeC = 0.0;
            if (sum_weight > 0.0) {
                nodeC = sum_weighted_conc / sum_weight;
            } else {
                nodeC = node->mixed_conc ? node->mixed_conc[sp] : 0.0;
            }

            /* Save node mixing concentration for reference */
            if (node->mixed_conc) node->mixed_conc[sp] = nodeC;

            /* Second pass: assign concentrations to branch boundaries */
            for (int c = 0; c < node->num_connections; ++c) {
                int branch_idx = node->connected_branches[c];
                Branch *b = net->branches[branch_idx];
                if (!b) continue;

                int dir = node->connection_dir[c];
                double branch_weight = weights[c];
                double branch_conc = conc_samples[c];
                double nodeC_excl = nodeC;
                double remaining = sum_weight - branch_weight;
                if (remaining > 0.0) {
                    nodeC_excl = (sum_weighted_conc - branch_weight * branch_conc) / remaining;
                }
                
                if (dir == -1) {
                    set_node_boundary_conc(b, sp, nodeC_excl, 1);
                } else if (dir == 1) {
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
            /* Note: Q_net is the inflow value. */
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

        /* B. Update Junction Nodes using Mass Balance */
        for (size_t n = 0; n < net->num_nodes; ++n) {
            Node *node = &net->nodes[n];
            if (node->type != NODE_JUNCTION) continue;
            
            double sum_Q_in = 0.0;
            double sum_Q_out = 0.0;
            double storage_area = 1000.0; /* Minimum storage */

            for (int c = 0; c < node->num_connections; ++c) {
                int branch_idx = node->connected_branches[c];
                Branch *b = net->branches[branch_idx];
                if (!b) continue;
                
                int dir = node->connection_dir[c]; 
                /* dir=1: Branch flows INTO Node (Node is Downstream) */
                /* dir=-1: Branch flows OUT of Node (Node is Upstream) */

                double Q_branch = 0.0;
                
                if (dir == 1) { 
                    /* Flow at index 1 is the outlet of the branch */
                    Q_branch = b->velocity[1] * b->totalArea[1];
                    /* Standard: Positive U means flow towards 1 (Downstream) */
                    /* So Q > 0 means flow INTO node */
                     if (Q_branch > 0) sum_Q_in += Q_branch;
                     else sum_Q_out += -Q_branch;
                     
                     storage_area += b->width[1] * b->dx;
                } 
                else if (dir == -1) {
                    /* Flow at index M is the inlet of the branch */
                    Q_branch = b->velocity[b->M] * b->totalArea[b->M];
                    /* Positive U means flow towards 1 (away from Node) */
                    /* So Q > 0 means flow OUT of node */
                    if (Q_branch > 0) sum_Q_out += Q_branch;
                    else sum_Q_in += -Q_branch;
                    
                    storage_area += b->width[b->M] * b->dx;
                }
            }
            
            double residual_Q = sum_Q_in - sum_Q_out;
            
            /* Relaxation update: dH = (Qin - Qout) * dt / Area */
            double dH = (residual_Q * dt / storage_area) * relaxation;
            
            /* Limit dH to prevent explosion */
            if (dH > 0.2) dH = 0.2;
            if (dH < -0.2) dH = -0.2;
            
            node->H += dH;
            
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

    /* Prepare junction-based per-branch boundary concentrations before transport */
    mix_junction_concentrations(net);

    /* Step 3: Transport */
    for (size_t i = 0; i < net->num_branches; ++i) {
        Branch *b = net->branches[i];
        if (!b) continue;
        
        /* Update BCs for transport */
        /* The branch-level boundary concentrations are set centrally by
         * mix_junction_concentrations() for junction nodes. Do not override
         * branch conc_up here. For discharge/level BCs, initialization/forcing
         * is left unchanged.
         */
        
        Transport_Branch(b, dt);
        Sediment_Branch(b, dt);
        
        if (current_time_seconds >= net->warmup_time) {
            Biogeo_Branch(b, dt);
        }
    }
    
    return 0;
}

