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
    
    /* Reaction rates array */
    b->num_reactions = CGEM_NUM_REACTIONS;
    b->reaction_rates = (double **)calloc((size_t)b->num_reactions, sizeof(double *));
    for (int r = 0; r < b->num_reactions; ++r) {
        b->reaction_rates[r] = (double *)calloc(field_len, sizeof(double));
    }
    
    /* Default transport parameters */
    b->D0 = 200.0;
    b->vdb_coef = 4.0;

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
    free(branch->tri_lower);
    free(branch->tri_diag);
    free(branch->tri_upper);
    free(branch->tri_rhs);
    free(branch->tri_gam);
    free(branch->tau_ero);
    free(branch->tau_dep);
    free(branch->mero);
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
    free(node);
}
