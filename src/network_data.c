#include "network.h"
#include <stdlib.h>
#include <string.h>

Branch *allocate_branch(int M, int num_species) {
    Branch *b = (Branch *)calloc(1, sizeof(Branch));
    if (!b) {
        return NULL;
    }

    size_t field_len = (size_t)(M + 2);
    b->M = M;
    b->num_species = num_species;

    /* Hydrodynamic arrays */
    b->totalArea = (double *)calloc(field_len, sizeof(double));
    b->freeArea = (double *)calloc(field_len, sizeof(double));
    b->refArea = (double *)calloc(field_len, sizeof(double));
    b->width = (double *)calloc(field_len, sizeof(double));
    b->depth = (double *)calloc(field_len, sizeof(double));
    b->velocity = (double *)calloc(field_len, sizeof(double));
    b->waterLevel = (double *)calloc(field_len, sizeof(double));
    b->dispersion = (double *)calloc(field_len, sizeof(double));
    b->chezyArray = (double *)calloc(field_len, sizeof(double));
    
    /* Previous timestep values */
    b->totalArea_old = (double *)calloc(field_len, sizeof(double));
    b->totalArea_old2 = (double *)calloc(field_len, sizeof(double));
    
    /* =========================================================================
     * RESIDUAL VELOCITY FILTER (Audit recommendation)
     * Low-pass filtered velocity for correct Van den Burgh dispersion
     * =========================================================================*/
    b->u_residual = (double *)calloc(field_len, sizeof(double));
    b->residual_alpha = 0.0;  /* Will be computed from dt/T_filter */
    
    /* Manning's n friction (optional alternative to Chezy) */
    b->manning_n = (double *)calloc(field_len, sizeof(double));
    b->use_manning = 0;  /* Default: use Chezy */
    
    /* Storage width ratio dynamics (floodplain effects) */
    b->H_bank = 2.0;          /* Default bank-full level [m] */
    b->RS_channel = 1.0;      /* In-channel storage ratio */
    b->RS_floodplain = 5.0;   /* Flooded storage ratio */
    
    /* Tridiagonal solver arrays */
    b->tri_lower = (double *)calloc(field_len, sizeof(double));
    b->tri_diag = (double *)calloc(field_len, sizeof(double));
    b->tri_upper = (double *)calloc(field_len, sizeof(double));
    b->tri_rhs = (double *)calloc(field_len, sizeof(double));
    b->tri_gam = (double *)calloc(field_len, sizeof(double));

    /* Concentration arrays */
    if (num_species > 0) {
        b->conc = (double **)calloc((size_t)num_species, sizeof(double *));
        b->conc_down = (double *)calloc((size_t)num_species, sizeof(double));
        b->conc_up = (double *)calloc((size_t)num_species, sizeof(double));
        if (b->conc) {
            for (int s = 0; s < num_species; ++s) {
                b->conc[s] = (double *)calloc(field_len, sizeof(double));
            }
        }
    }
    
    /* Sediment arrays */
    b->tau_ero = (double *)calloc(field_len, sizeof(double));
    b->tau_dep = (double *)calloc(field_len, sizeof(double));
    b->mero = (double *)calloc(field_len, sizeof(double));
    
    /* =======================================================================
     * SEDIMENT BED LAYER (Fluid Mud Tracking)
     * Tracks accumulated sediment and organic carbon for benthic-pelagic coupling
     * Reference: Winterwerp (2002), Wang et al. (2018)
     * =======================================================================*/
    b->bed_mass = (double *)calloc(field_len, sizeof(double));
    b->bed_oc = (double *)calloc(field_len, sizeof(double));
    
    /* Reaction rates array */
    b->num_reactions = CGEM_NUM_REACTIONS;
    b->reaction_rates = (double **)calloc((size_t)b->num_reactions, sizeof(double *));
    for (int r = 0; r < b->num_reactions; ++r) {
        b->reaction_rates[r] = (double *)calloc(field_len, sizeof(double));
    }
    
    /* =======================================================================
     * LATERAL LOADS (Land-Use Coupling)
     * Allocate arrays for spatially-distributed lateral inflows and loads
     * These will be populated by LoadLateralSources() from CSV
     * =======================================================================*/
    b->lateral_flow = (double *)calloc(field_len, sizeof(double));
    b->lateral_conc = (double **)calloc((size_t)num_species, sizeof(double *));
    for (int s = 0; s < num_species; ++s) {
        b->lateral_conc[s] = (double *)calloc(field_len, sizeof(double));
    }
    b->has_lateral_loads = 0;  /* Will be set to 1 if loads are loaded */
    
    /* Default transport parameters */
    b->D0 = 50.0;   /* Reduced to prevent salt accumulation */
    /* Van den Burgh K coefficient (Savenije, 2005 Table 9.1):
     * 
     * SCIENTIFIC BASIS:
     * K controls the exponential decay of dispersion: D(x) = D0 * exp(-x/L_d)
     * where L_d = D0 * A / (K * Q_f)
     * 
     * - K = 0.2-0.4: Flatter profile, salt penetrates 50-70 km (deep alluvial)
     * - K = 0.5-0.7: Steeper decay, salt stays near mouth (shallow estuaries)
     * - K = 0.8-1.0: Very steep, minimal intrusion (narrow channels)
     * 
     * For MEKONG distributaries (H > 10m, funnel-shaped):
     * - Observed salt intrusion: 40-60 km in dry season
     * - Calibrated K = 0.30-0.40 (Nguyen et al., 2008; Savenije, 2005)
     * 
     * TUNING NOTE: Decrease K to allow more salt intrusion upstream.
     * 
     * Reference: Savenije (2005) "Salinity and Tides in Alluvial Estuaries"
     */
    b->vdb_coef = 0.35;

    /* Initialize CSV file pointers to NULL to avoid accidental reuse */
    for (int i = 0; i < 25; ++i) b->csv_fps[i] = NULL;
    b->num_csv_fps = 0;

    /* No matrix CSV (relying on per-variable CSVs) */

    return b;
}

void free_branch(Branch *branch, int num_species) {
    if (!branch) {
        return;
    }

    free(branch->totalArea);
    free(branch->freeArea);
    free(branch->refArea);
    free(branch->width);
    free(branch->depth);
    free(branch->velocity);
    free(branch->waterLevel);
    free(branch->dispersion);
    free(branch->chezyArray);
    free(branch->totalArea_old);
    free(branch->totalArea_old2);
    free(branch->u_residual);
    free(branch->manning_n);
    free(branch->tri_lower);
    free(branch->tri_diag);
    free(branch->tri_upper);
    free(branch->tri_rhs);
    free(branch->tri_gam);
    free(branch->tau_ero);
    free(branch->tau_dep);
    free(branch->mero);
    free(branch->bed_mass);
    free(branch->bed_oc);
    free(branch->conc_down);
    free(branch->conc_up);

    if (branch->conc) {
        for (int s = 0; s < num_species; ++s) {
            free(branch->conc[s]);
        }
        free(branch->conc);
    }

    if (branch->reaction_rates) {
        for (int r = 0; r < branch->num_reactions; ++r) {
            free(branch->reaction_rates[r]);
        }
        free(branch->reaction_rates);
    }
    
    /* Free lateral load arrays */
    free(branch->lateral_flow);
    if (branch->lateral_conc) {
        for (int s = 0; s < num_species; ++s) {
            free(branch->lateral_conc[s]);
        }
        free(branch->lateral_conc);
    }

    /* No matrix CSV to close; per-variable CSVs are closed elsewhere */

    free(branch);
}

Node *allocate_node(int num_species) {
    Node *node = (Node *)calloc(1, sizeof(Node));
    if (!node) {
        return NULL;
    }

    if (num_species > 0) {
        node->mixed_conc = (double *)calloc((size_t)num_species, sizeof(double));
    }

    return node;
}

void free_node(Node *node) {
    if (!node) {
        return;
    }

    free(node->forcing_time);
    free(node->forcing_value);
    free(node->mixed_conc);
    
    /* Free species forcing arrays */
    if (node->species_forcing_time) {
        for (int sp = 0; sp < node->num_species_forcing; ++sp) {
            free(node->species_forcing_time[sp]);
        }
        free(node->species_forcing_time);
    }
    if (node->species_forcing_value) {
        for (int sp = 0; sp < node->num_species_forcing; ++sp) {
            free(node->species_forcing_value[sp]);
        }
        free(node->species_forcing_value);
    }
    free(node->species_forcing_len);
    
    free(node);
}
