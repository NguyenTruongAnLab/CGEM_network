#include "../network.h"
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

static char *trim_in_place(char *text) {
    if (!text) {
        return text;
    }
    unsigned char *start = (unsigned char *)text;
    while (*start && isspace(*start)) {
        ++start;
    }
    if (*start == '\0') {
        *text = '\0';
        return text;
    }
    unsigned char *end = start + strlen((const char *)start) - 1;
    while (end > start && isspace(*end)) {
        *end-- = '\0';
    }
    if (start != (unsigned char *)text) {
        memmove(text, start, strlen((const char *)start) + 1);
    }
    return text;
}

static int string_ieq(const char *a, const char *b) {
    if (!a || !b) {
        return 0;
    }
    while (*a && *b) {
        char ca = (char)tolower((unsigned char)*a++);
        char cb = (char)tolower((unsigned char)*b++);
        if (ca != cb) {
            return 0;
        }
    }
    return *a == '\0' && *b == '\0';
}

/* Normalize GHG concentration units: input files use nmol/L, model uses µmol/L */


/**
 * Parse species name string to species index (case-insensitive)
 * Uses canonical names matching CGEM_SPECIES_* constants in define.h
 * @param name Species name
 * @return Species index, or -1 if not found
 */
static int parse_species_name(const char *name) {
    if (!name) return -1;
    
    /* Core species (0-16) */
    if (string_ieq(name, "SALINITY")) return CGEM_SPECIES_SALINITY;
    if (string_ieq(name, "PHY1")) return CGEM_SPECIES_PHY1;
    if (string_ieq(name, "PHY2")) return CGEM_SPECIES_PHY2;
    if (string_ieq(name, "DSI")) return CGEM_SPECIES_DSI;
    if (string_ieq(name, "NO3")) return CGEM_SPECIES_NO3;
    if (string_ieq(name, "NH4")) return CGEM_SPECIES_NH4;
    if (string_ieq(name, "PO4")) return CGEM_SPECIES_PO4;
    if (string_ieq(name, "O2")) return CGEM_SPECIES_O2;
    if (string_ieq(name, "TOC")) return CGEM_SPECIES_TOC;
    if (string_ieq(name, "SPM")) return CGEM_SPECIES_SPM;
    if (string_ieq(name, "DIC")) return CGEM_SPECIES_DIC;
    if (string_ieq(name, "AT")) return CGEM_SPECIES_AT;
    if (string_ieq(name, "PCO2")) return CGEM_SPECIES_PCO2;
    if (string_ieq(name, "CO2")) return CGEM_SPECIES_CO2;
    if (string_ieq(name, "PH")) return CGEM_SPECIES_PH;
    if (string_ieq(name, "HS")) return CGEM_SPECIES_HS;
    if (string_ieq(name, "ALKC")) return CGEM_SPECIES_ALKC;
    
    /* RIVE multi-pool organic matter (17-22) */
    if (string_ieq(name, "HD1")) return CGEM_SPECIES_HD1;
    if (string_ieq(name, "HD2")) return CGEM_SPECIES_HD2;
    if (string_ieq(name, "HD3")) return CGEM_SPECIES_HD3;
    if (string_ieq(name, "HP1")) return CGEM_SPECIES_HP1;
    if (string_ieq(name, "HP2")) return CGEM_SPECIES_HP2;
    if (string_ieq(name, "HP3")) return CGEM_SPECIES_HP3;
    
    /* RIVE bacteria and phosphorus (23-26) */
    if (string_ieq(name, "BAG")) return CGEM_SPECIES_BAG;
    if (string_ieq(name, "BAP")) return CGEM_SPECIES_BAP;
    if (string_ieq(name, "PIP")) return CGEM_SPECIES_PIP;
    if (string_ieq(name, "DSS")) return CGEM_SPECIES_DSS;
    
    /* GHG species (27-29) */
    if (string_ieq(name, "NO2")) return CGEM_SPECIES_NO2;
    if (string_ieq(name, "N2O")) return CGEM_SPECIES_N2O;
    if (string_ieq(name, "CH4")) return CGEM_SPECIES_CH4;
    if (string_ieq(name, "NITROUSOXIDE")) return CGEM_SPECIES_N2O;
    if (string_ieq(name, "METHANE")) return CGEM_SPECIES_CH4;
    
    /* 2-Pool TOC model (30-31) - SCIENTIFIC FIX December 2025 */
    if (string_ieq(name, "TOC_LABILE")) return CGEM_SPECIES_TOC_LABILE;
    if (string_ieq(name, "TOC_REFRACTORY")) return CGEM_SPECIES_TOC_REFRACTORY;
    if (string_ieq(name, "TOC_LAB")) return CGEM_SPECIES_TOC_LABILE;
    if (string_ieq(name, "TOC_REF")) return CGEM_SPECIES_TOC_REFRACTORY;
    if (string_ieq(name, "LABILE_TOC")) return CGEM_SPECIES_TOC_LABILE;
    if (string_ieq(name, "REFRACTORY_TOC")) return CGEM_SPECIES_TOC_REFRACTORY;
    
    return -1;
}

/**
 * Load ALL species forcing from a multi-column CSV file
 * Header format: time_s,salinity,phy1,phy2,...
 * Parses header to determine which columns map to which species.
 * 
 * @param full_path Path to CSV file
 * @param node Node to populate with species forcing data
 * @param dt Time step between data points [s]
 * @return 0 on success, -1 on error
 */
static int load_all_species_forcing(const char *full_path, Node *node, double dt) {
    FILE *fp = fopen(full_path, "r");
    if (!fp) {
        fprintf(stderr, "Failed to open species forcing file: %s\n", full_path);
        return -1;
    }
    
    char line[4096];
    char *text;
    
    /* Read header line to get column names */
    char header[4096] = {0};
    while (fgets(line, sizeof(line), fp)) {
        text = trim_in_place(line);
        if (*text != '\0' && *text != '#') {
            strncpy(header, text, sizeof(header) - 1);
            break;
        }
    }
    
    if (header[0] == '\0') {
        fprintf(stderr, "No header found in species file: %s\n", full_path);
        fclose(fp);
        return -1;
    }
    
    /* Parse header to get column indices for each species */
    int col_to_species[64] = {-1};  /* Map column index to species index */
    int num_cols = 0;
    int time_col = -1;
    
    char header_copy[4096];
    strncpy(header_copy, header, sizeof(header_copy) - 1);
    char *token = strtok(header_copy, ",");
    while (token && num_cols < 64) {
        char *col_name = trim_in_place(token);
        if (string_ieq(col_name, "time_s") || string_ieq(col_name, "time")) {
            time_col = num_cols;
            col_to_species[num_cols] = -1;
        } else if (string_ieq(col_name, "temperature") || string_ieq(col_name, "temp")) {
            col_to_species[num_cols] = -1;  /* Skip temperature column */
        } else {
            int sp_idx = parse_species_name(col_name);
            col_to_species[num_cols] = sp_idx;
            if (sp_idx >= 0) {
                /* printf("  Column %d (%s) -> species %d\n", num_cols, col_name, sp_idx); */
            }
        }
        num_cols++;
        token = strtok(NULL, ",");
    }
    
    if (time_col < 0) {
        fprintf(stderr, "Warning: No 'time_s' column found, using row index for time\n");
    }
    
    /* Count data rows */
    size_t count = 0;
    while (fgets(line, sizeof(line), fp)) {
        text = trim_in_place(line);
        if (*text != '\0' && *text != '#') {
            count++;
        }
    }
    rewind(fp);
    
    /* Skip header again */
    while (fgets(line, sizeof(line), fp)) {
        text = trim_in_place(line);
        if (*text != '\0' && *text != '#') break;
    }
    
    /* Allocate species forcing arrays if needed */
    if (!node->species_forcing_time) {
        node->species_forcing_time = (double **)calloc(CGEM_NUM_SPECIES, sizeof(double *));
        node->species_forcing_value = (double **)calloc(CGEM_NUM_SPECIES, sizeof(double *));
        node->species_forcing_len = (size_t *)calloc(CGEM_NUM_SPECIES, sizeof(size_t));
        node->num_species_forcing = CGEM_NUM_SPECIES;
    }
    
    /* Allocate arrays for each mapped species */
    for (int c = 0; c < num_cols; c++) {
        int sp = col_to_species[c];
        if (sp >= 0 && sp < CGEM_NUM_SPECIES) {
            node->species_forcing_time[sp] = (double *)malloc(count * sizeof(double));
            node->species_forcing_value[sp] = (double *)malloc(count * sizeof(double));
            node->species_forcing_len[sp] = count;
        }
    }
    
    /* Read data rows */
    size_t row_idx = 0;
    while (fgets(line, sizeof(line), fp) && row_idx < count) {
        text = trim_in_place(line);
        if (*text == '\0' || *text == '#') continue;
        
        /* Parse each column */
        char line_copy[4096];
        strncpy(line_copy, text, sizeof(line_copy) - 1);
        
        double time_val = (time_col >= 0) ? 0.0 : (double)row_idx * dt;
        double col_vals[64] = {0.0};
        
        token = strtok(line_copy, ",");
        int col = 0;
        while (token && col < num_cols) {
            double val = strtod(trim_in_place(token), NULL);
            if (col == time_col) {
                time_val = val;
            }
            col_vals[col] = val;
            col++;
            token = strtok(NULL, ",");
        }
        
        /* Store values for each species */
        for (int c = 0; c < num_cols; c++) {
            int sp = col_to_species[c];
            if (sp >= 0 && sp < CGEM_NUM_SPECIES && 
                node->species_forcing_time[sp] && node->species_forcing_value[sp]) {
                node->species_forcing_time[sp][row_idx] = time_val;

                /* Expect input already in µmol/L for all species (GHG included) */
                node->species_forcing_value[sp][row_idx] = col_vals[c];
            }
        }
        row_idx++;
    }
    
    /* Update actual lengths */
    for (int c = 0; c < num_cols; c++) {
        int sp = col_to_species[c];
        if (sp >= 0 && sp < CGEM_NUM_SPECIES) {
            node->species_forcing_len[sp] = row_idx;
        }
    }
    
    fclose(fp);
    
    /* Count how many species were loaded */
    int loaded_count = 0;
    for (int s = 0; s < CGEM_NUM_SPECIES; s++) {
        if (node->species_forcing_len[s] > 0) loaded_count++;
    }
    printf("Loaded %d species with %zu time points from: %s\n", loaded_count, row_idx, full_path);
    
    return 0;
}

/**
 * Load forcing data from a CSV file
 * @param full_path Path to CSV file
 * @param out_time Output array for time values
 * @param out_value Output array for data values
 * @param out_len Output length
 * @param dt Time step between data points [s]
 * @return 0 on success, -1 on error
 */
static int load_forcing_file(const char *full_path, double **out_time, double **out_value,
                             size_t *out_len, double dt) {
    FILE *ffp = fopen(full_path, "r");
    if (!ffp) {
        fprintf(stderr, "Failed to open forcing file: %s\n", full_path);
        return -1;
    }

    char line[1024];
    char *text;
    
    /* Count lines, skip header */
    size_t count = 0;
    int has_header = 0;
    while (fgets(line, sizeof(line), ffp)) {
        text = trim_in_place(line);
        if (*text != '\0' && *text != '#') {
            if (!has_header) {
                has_header = 1;
                continue;
            }
            count++;
        }
    }
    rewind(ffp);

    *out_len = count;
    *out_time = (double *)malloc(count * sizeof(double));
    *out_value = (double *)malloc(count * sizeof(double));
    if (!*out_time || !*out_value) {
        fprintf(stderr, "Failed to allocate forcing arrays\n");
        fclose(ffp);
        return -1;
    }

    size_t idx = 0;
    int is_header = 1;
    while (fgets(line, sizeof(line), ffp) && idx < count) {
        text = trim_in_place(line);
        if (*text == '\0' || *text == '#') {
            continue;
        }
        if (is_header) {
            is_header = 0;
            continue;
        }

        char *comma = strchr(text, ',');
        if (comma) {
            *comma = '\0';
            double data_val = strtod(comma + 1, NULL);
            (*out_time)[idx] = (double)idx * dt;
            (*out_value)[idx] = data_val;
            idx++;
        }
    }
    *out_len = idx;
    fclose(ffp);
    return 0;
}

static int ensure_branch_capacity(Branch ***branches, size_t *capacity, size_t required) {
    if (*capacity >= required) {
        return 0;
    }
    size_t new_capacity = (*capacity == 0) ? 8 : (*capacity * 2);
    while (new_capacity < required) {
        new_capacity *= 2;
    }
    Branch **resized = (Branch **)realloc(*branches, new_capacity * sizeof(Branch *));
    if (!resized) {
        return -1;
    }
    *branches = resized;
    *capacity = new_capacity;
    return 0;
}

int LoadTopology(const char *path, Network *net) {
    if (!path || !net) {
        return -1;
    }

    FILE *fp = fopen(path, "r");
    if (!fp) {
        fprintf(stderr, "Failed to open topology file: %s\n", path);
        return -2;
    }

    Branch **branches = NULL;
    size_t count = 0;
    size_t capacity = 0;
    int status = 0;

    char line[1024];
    while (fgets(line, sizeof(line), fp)) {
        char *text = trim_in_place(line);
        if (*text == '\0' || *text == '#') {
            continue;
        }

        if (ensure_branch_capacity(&branches, &capacity, count + 1) != 0) {
            status = -3;
            goto cleanup;
        }
        /* Parse fields first to know length to choose M for allocation */

        /* Parse topology columns (Savenije-inspired parameterisation):
         *  0: ID
         *  1: Name
         *  2: NodeUp
         *  3: NodeDown
         *  4: Length_m
         *  5: Width_Up_m
         *  6: Width_Down_m
         *  7: Depth_m
         *  8: Chezy
         *  9: Group (optional, default 0)
         * 10: RS   (optional, default 1.0) - Storage width ratio for mangroves
         * 11: VDB_K (optional, default 0.35) - Van den Burgh dispersion coefficient
         * 12: Mixing_Alpha (optional, default 0.25) - D0 mixing efficiency: D0 = α * U_tidal * B
         * 13: BiogeoParams (optional) - Path to branch-specific biogeo params file
         */
        int column = 0;
        int tmp_id = -1;
        char tmp_name[CGEM_MAX_BRANCH_NAME] = {0};
        int tmp_node_up = -1;
        int tmp_node_down = -1;
        double tmp_length = 0.0;
        double tmp_width_up = 0.0;
        double tmp_width_down = 0.0;
        double tmp_depth = 0.0;
        double tmp_chezy = 0.0;
        int tmp_group = 0;
        double tmp_storage_ratio = 1.0;      /* Default RS = 1.0 (prismatic channel) */
        double tmp_vdb_coef = 0.35;          /* Default Van den Burgh K (Savenije, 2005) */
        double tmp_mixing_alpha = 0.25;      /* Default mixing efficiency for D0 (Fischer, 1979) */
        char tmp_biogeo_path[CGEM_MAX_PATH] = {0}; /* Empty = use global defaults */
        char *token = strtok(text, ",");
        while (token && column < 14) {
            char *value = trim_in_place(token);
            switch (column) {
                case 0:
                    tmp_id = (int)strtol(value, NULL, 10);
                    break;
                case 1:
                    snprintf(tmp_name, sizeof(tmp_name), "%s", value);
                    break;
                case 2:
                    /* Keep node IDs as-is from topology.csv (no -1 offset)
                     * This matches boundary_map.csv which uses same node IDs */
                    tmp_node_up = (int)strtol(value, NULL, 10);
                    break;
                case 3:
                    /* Keep node IDs as-is from topology.csv (no -1 offset) */
                    tmp_node_down = (int)strtol(value, NULL, 10);
                    break;
                case 4:
                    tmp_length = strtod(value, NULL);
                    break;
                case 5:
                    tmp_width_up = strtod(value, NULL);
                    break;
                case 6:
                    tmp_width_down = strtod(value, NULL);
                    break;
                case 7:
                    tmp_depth = strtod(value, NULL);
                    break;
                case 8:
                    tmp_chezy = strtod(value, NULL);
                    break;
                case 9:
                    tmp_group = (int)strtol(value, NULL, 10);
                    break;
                case 10:
                    /* RS: Storage width ratio - critical for mangrove/tidal flat branches */
                    tmp_storage_ratio = strtod(value, NULL);
                    if (tmp_storage_ratio < 0.1) tmp_storage_ratio = 1.0;  /* Sanity check */
                    break;
                case 11:
                    /* Van den Burgh dispersion coefficient (controls exponential decay length) */
                    tmp_vdb_coef = strtod(value, NULL);
                    break;
                case 12:
                    /* Mixing_Alpha: D0 mixing efficiency coefficient (Fischer formula)
                     * D0 = α * U_tidal * B where α typically 0.1-0.5 for estuaries
                     * Reference: Fischer et al. (1979), Savenije (2005) */
                    tmp_mixing_alpha = strtod(value, NULL);
                    if (tmp_mixing_alpha < 0.01) tmp_mixing_alpha = 0.25;  /* Sanity check */
                    if (tmp_mixing_alpha > 2.0) tmp_mixing_alpha = 2.0;    /* Cap at max */
                    break;
                case 13:
                    /* BiogeoParams: Path to branch-specific biogeochemistry parameters */
                    snprintf(tmp_biogeo_path, sizeof(tmp_biogeo_path), "%s", value);
                    break;
                default:
                    break;
            }
            token = strtok(NULL, ",");
            column++;
        }

        if (column < 9) {
            fprintf(stderr, "Invalid topology row: expected at least 9 columns, got %d\n", column);
            status = -5;
            goto cleanup;
        }

        /* Defaults for optional columns */
        if (column < 10) {
            tmp_group = 0;
        }
        if (column < 11) {
            tmp_storage_ratio = 1.0;
        }
        if (column < 12) {
            tmp_vdb_coef = 0.35;  /* Literature mid-range if not specified */
        }
        if (column < 13) {
            tmp_mixing_alpha = 0.25;  /* Fischer (1979) default for estuaries */
        }
        /* tmp_biogeo_path defaults to empty string (use global params) */

        if (tmp_length <= 0.0) {
            tmp_length = CGEM_DEFAULT_BRANCH_CELLS * CGEM_DEFAULT_DX_METERS;
        }

        /* Determine M based on net->dx_target */
        double target_dx = net->dx_target > 0.0 ? net->dx_target : CGEM_DEFAULT_DX_METERS;
        /* Determine M as nearest even integer to tmp_length / target_dx */
        double steps_f = tmp_length / target_dx;
        int Mcalc = (int)lround(steps_f);
        if (Mcalc < 2) Mcalc = 2;
        if (Mcalc % 2 != 0) {
            int M_minus = (Mcalc - 1 >= 2) ? Mcalc - 1 : 2;
            int M_plus = Mcalc + 1;
            if (fabs(steps_f - (double)M_minus) <= fabs((double)M_plus - steps_f)) {
                Mcalc = M_minus;
            } else {
                Mcalc = M_plus;
            }
        }

        Branch *branch = allocate_branch(Mcalc, net->num_species);
        if (!branch) {
            status = -4;
            goto cleanup;
        }

        /* Copy in parsed fields into allocated branch */
        branch->id = tmp_id;
        snprintf(branch->name, sizeof(branch->name), "%s", tmp_name);
        branch->node_up = tmp_node_up;
        branch->node_down = tmp_node_down;
        branch->length_m = tmp_length; /* will be rounded in geometry init */
        branch->width_up_m = tmp_width_up;
        branch->width_down_m = tmp_width_down;
        branch->depth_m = tmp_depth > 0.0 ? tmp_depth : CGEM_MIN_DEPTH;
        branch->chezy = tmp_chezy;
        branch->group_id = tmp_group;
        branch->storage_ratio = tmp_storage_ratio;  /* RS for mangrove/tidal flat storage */
        branch->has_biogeo = 1;  /* Always enable biogeo for now */
        branch->vdb_coef = tmp_vdb_coef;
        branch->mixing_alpha = tmp_mixing_alpha;
        snprintf(branch->biogeo_params_path, sizeof(branch->biogeo_params_path), "%s", tmp_biogeo_path);
        
        /* Diagnostic output for configuration verification */
        printf("  Parsed branch %d: %s (RS=%.1f, K=%.2f, α=%.2f", tmp_id, tmp_name, tmp_storage_ratio, tmp_vdb_coef, tmp_mixing_alpha);
        if (tmp_biogeo_path[0] != '\0') {
            printf(", biogeo=%s", tmp_biogeo_path);
        }
        printf(")\n");

        branches[count++] = branch;
    }

    if (count == 0) {
        fprintf(stderr, "Topology file contains no branch definitions: %s\n", path);
        status = -6;
        goto cleanup;
    }

cleanup:
    fclose(fp);

    if (status != 0) {
        if (branches) {
            for (size_t i = 0; i < count; ++i) {
                if (branches[i]) {
                    free_branch(branches[i], branches[i]->num_species);
                }
            }
            free(branches);
        }
        return status;
    }

    net->branches = branches;
    net->num_branches = count;

    // Allocate nodes
    int max_node_id = 0;
    for (size_t i = 0; i < count; ++i) {
        if (branches[i]->node_up > max_node_id) max_node_id = branches[i]->node_up;
        if (branches[i]->node_down > max_node_id) max_node_id = branches[i]->node_down;
    }
    net->num_nodes = (size_t)max_node_id + 1;
    net->nodes = (Node *)calloc(net->num_nodes, sizeof(Node));
    if (!net->nodes) {
        // cleanup
        for (size_t i = 0; i < count; ++i) {
            free_branch(branches[i], net->num_species);
        }
        free(branches);
        return -7;
    }
    for (size_t i = 0; i < net->num_nodes; ++i) {
        net->nodes[i].id = (int)i;
        // Initialize node types and connections later
    }

    // Set node connections (types will be set properly by LoadBoundaries)
    for (size_t i = 0; i < net->num_nodes; ++i) {
        int connections = 0;
        for (size_t j = 0; j < count; ++j) {
            if (branches[j]->node_up == (int)i) {
                net->nodes[i].connected_branches[connections] = (int)j;
                net->nodes[i].connection_dir[connections] = -1;  // outflow
                connections++;
            } else if (branches[j]->node_down == (int)i) {
                net->nodes[i].connected_branches[connections] = (int)j;
                net->nodes[i].connection_dir[connections] = 1;  // inflow
                connections++;
            }
        }
        net->nodes[i].num_connections = connections;
        
        /* Default type assignment based on connections:
         * - Multi-connection nodes are always junctions
         * - Single-connection nodes: default to DISCHARGE_BC, but LoadBoundaries
         *   will overwrite with LEVEL_BC for ocean boundaries based on boundary_map.csv
         */
        if (connections > 1) {
            net->nodes[i].type = NODE_JUNCTION;
        } else if (connections == 1) {
            /* Default to DISCHARGE_BC - will be corrected by LoadBoundaries if ocean */
            net->nodes[i].type = NODE_DISCHARGE_BC;
        } else {
            /* No connections - orphan node, treat as junction placeholder */
            net->nodes[i].type = NODE_JUNCTION;
        }
    }

    return 0;
}

int LoadBoundaries(const char *path, Network *net) {
    if (!path || !net) {
        return -1;
    }

    FILE *fp = fopen(path, "r");
    if (!fp) {
        fprintf(stderr, "Failed to open boundary map file: %s\n", path);
        return -2;
    }

    /* Extract directory from path for relative file paths */
    char dir[CGEM_MAX_PATH];
    strcpy(dir, path);
    char *last_slash = strrchr(dir, '/');
    if (!last_slash) last_slash = strrchr(dir, '\\');
    if (last_slash) {
        *last_slash = '\0';
    } else {
        dir[0] = '\0';
    }

    char line[1024];
    while (fgets(line, sizeof(line), fp)) {
        char *text = trim_in_place(line);
        if (*text == '\0' || *text == '#') {
            continue;
        }

        int node_id = -1;
        char type_str[32] = {0};
        char file_path[CGEM_MAX_PATH] = {0};
        char species_str[32] = {0};

        /* Parse up to 4 columns for SPECIES type support */
        int column = 0;
        char line_copy[1024];
        strcpy(line_copy, text);
        char *token = strtok(line_copy, ",");
        while (token && column < 4) {
            char *value = trim_in_place(token);
            switch (column) {
                case 0:
                    node_id = (int)strtol(value, NULL, 10);
                    break;
                case 1:
                    snprintf(type_str, sizeof(type_str), "%s", value);
                    break;
                case 2:
                    snprintf(file_path, sizeof(file_path), "%s", value);
                    break;
                case 3:
                    snprintf(species_str, sizeof(species_str), "%s", value);
                    break;
                default:
                    break;
            }
            token = strtok(NULL, ",");
            column++;
        }

        if (column < 3 || node_id < 0 || node_id >= (int)net->num_nodes) {
            fprintf(stderr, "Invalid boundary row: node_id=%d, columns=%d, max_nodes=%zu\n", 
                    node_id, column, net->num_nodes);
            continue;
        }

        /* 
         * CRITICAL FIX: Node IDs in boundary_map.csv match topology.csv node IDs directly.
         * The internal nodes array is indexed by these IDs directly (no -1 offset).
         * This matches load_network_from_csv which uses branch->node_up/node_down directly.
         */
        Node *node = &net->nodes[node_id];

        /* Build full path */
        char full_path[CGEM_MAX_PATH + 256];
        if (dir[0] != '\0') {
            snprintf(full_path, sizeof(full_path), "%s/%s", dir, file_path);
        } else {
            strcpy(full_path, file_path);
        }

        if (string_ieq(type_str, "LEVEL")) {
            node->type = NODE_LEVEL_BC;
            if (load_forcing_file(full_path, &node->forcing_time, &node->forcing_value,
                                  &node->forcing_len, 3600.0) == 0) {
                printf("Loaded %zu forcing points for node %d (LEVEL)\n", node->forcing_len, node_id);
            }
        } 
        else if (string_ieq(type_str, "DISCHARGE")) {
            node->type = NODE_DISCHARGE_BC;
            if (load_forcing_file(full_path, &node->forcing_time, &node->forcing_value,
                                  &node->forcing_len, 86400.0) == 0) {
                printf("Loaded %zu forcing points for node %d (DISCHARGE)\n", node->forcing_len, node_id);
            }
        }
        else if (string_ieq(type_str, "SPECIES")) {
            /* Time-varying species boundary condition */
            /* Format: Node_ID, SPECIES, FilePath, SpeciesName (or ALL) */
            if (column < 4) {
                fprintf(stderr, "SPECIES boundary requires 4 columns: Node_ID, SPECIES, FilePath, SpeciesName\n");
                continue;
            }
            
            /* Check if loading ALL species from multi-column CSV */
            if (string_ieq(species_str, "ALL")) {
                /* Load all species from a single multi-column CSV file */
                if (load_all_species_forcing(full_path, node, 86400.0) == 0) {
                    printf("Loaded ALL species forcing for node %d from: %s\n", node_id, file_path);
                }
            } else {
                /* Load single species from CSV */
                int species_idx = parse_species_name(species_str);
                if (species_idx < 0) {
                    fprintf(stderr, "Unknown species name: %s for node %d\n", species_str, node_id);
                    continue;
                }
                
                /* Allocate species forcing arrays if needed */
                if (!node->species_forcing_time) {
                    node->species_forcing_time = (double **)calloc(CGEM_NUM_SPECIES, sizeof(double *));
                    node->species_forcing_value = (double **)calloc(CGEM_NUM_SPECIES, sizeof(double *));
                    node->species_forcing_len = (size_t *)calloc(CGEM_NUM_SPECIES, sizeof(size_t));
                    node->num_species_forcing = CGEM_NUM_SPECIES;
                }
                
                /* Load species forcing file */
                if (load_forcing_file(full_path, 
                                      &node->species_forcing_time[species_idx],
                                      &node->species_forcing_value[species_idx],
                                      &node->species_forcing_len[species_idx],
                                  86400.0) == 0) {
                    /* Inputs must already be in µmol/L (no automatic conversion) */
                    printf("Loaded %zu forcing points for node %d (SPECIES: %s)\n", 
                           node->species_forcing_len[species_idx], node_id, species_str);
                }
            }
        }
        else {
            fprintf(stderr, "Unknown boundary type: %s for node %d\n", type_str, node_id);
            continue;
        }
    }

    fclose(fp);
    return 0;
}

/* ===========================================================================
 * LATERAL SOURCES LOADER
 * Reads spatially-explicit lateral loads from land-use analysis
 * 
 * CSV Format:
 * Branch, Segment_Index, Distance_km, Area_km2, Q_lat_m3_s, 
 * NH4_load_g_s, NO3_load_g_s, PO4_load_g_s, TOC_load_g_s, DIC_load_g_s,
 * NH4_conc_mg_L, NO3_conc_mg_L, PO4_conc_mg_L, TOC_conc_mg_L, DIC_conc_mg_L
 * 
 * Reference: JAXA land use data fusion approach
 * ===========================================================================*/

/**
 * Find branch by name (case-insensitive)
 */
static Branch *find_branch_by_name(Network *net, const char *name) {
    if (!net || !name) return NULL;
    
    for (size_t b = 0; b < net->num_branches; ++b) {
        if (string_ieq(net->branches[b]->name, name)) {
            return net->branches[b];
        }
    }
    return NULL;
}

/**
 * Load lateral sources from CSV file
 * 
 * Converts mass loads (g/s) to concentration increments per cell
 * The loads are applied as:
 *   dC = (Q_lat * C_lat) / (Area * depth) * dt
 * where Q_lat is lateral inflow, C_lat is load concentration
 * 
 * @param net Network structure
 * @param case_dir Path to case directory
 * @return 0 on success, -1 on failure
 */
int LoadLateralSources(Network *net, const char *case_dir) {
    if (!net || !case_dir) return -1;
    
    char path[CGEM_MAX_PATH];
    snprintf(path, sizeof(path), "%s/lateral_sources.csv", case_dir);
    
    FILE *fp = fopen(path, "r");
    if (!fp) {
        printf("No lateral sources file found: %s (continuing without lateral loads)\n", path);
        return 0;  /* Not an error - lateral loads are optional */
    }
    
    printf("Loading lateral sources from: %s\n", path);
    
    char line[2048];
    int line_num = 0;
    int total_loaded = 0;
    
    /* Unit conversion factors:
     * Input: g/s of mass load, m³/s of flow
     * Target: concentration in µmol/L (for nutrients) or mg/L (for OC)
     * 
     * For NH4 (g N/s → µmol N/L):
     *   g N/s ÷ (m³/s × 14 g/mol × 1e-6 mol/µmol × 1000 L/m³)
     *   = g N/s ÷ (Q_m3_s × 0.014)  [µmol N/L]
     * 
     * For TOC (g C/s → µmol C/L):
     *   g C/s ÷ (m³/s × 12 g/mol × 1e-6 mol/µmol × 1000 L/m³)
     *   = g C/s ÷ (Q_m3_s × 0.012)  [µmol C/L]
     */
    
    /* Read header line and detect format */
    if (!fgets(line, sizeof(line), fp)) {
        fclose(fp);
        return -1;
    }
    
    /* Detect which format we're reading:
     * OLD FORMAT: Branch,Segment_Index,Distance_km,Area_km2,Q_lat_m3_s,NH4_load_g_s,...
     * NEW FORMAT: Branch,Segment_Index,Distance_km,Area_km2,Runoff_C,Is_Polder_Zone,Q_lat_base_m3_s,NH4_conc_base_mg_L,...
     */
    int is_new_format = (strstr(line, "Q_lat_base") != NULL || strstr(line, "conc_base") != NULL);
    
    while (fgets(line, sizeof(line), fp)) {
        line_num++;
        
        /* Parse CSV line */
        char branch_name[64];
        int segment_idx;
        double dist_km, area_km2, Q_lat;
        double nh4_conc = 0, no3_conc = 0, po4_conc = 0, toc_conc = 0, dic_conc = 0, spm_conc = 0;
        double ch4_conc = 0, n2o_conc = 0, at_conc = 0;  /* Units: µmol/L for all species */
        
        int parsed = 0;
        
        if (is_new_format) {
            /* NEW FORMAT (v2): Branch,Segment_Index,Distance_km,Area_km2,Runoff_C,Is_Polder_Zone,
             * Q_lat_base_m3_s,NH4_conc_base_mg_L,NO3_conc_base_mg_L,PO4_conc_base_mg_L,
             * TOC_conc_base_mg_L,DIC_conc_base_mg_L,SPM_conc_base_mg_L,
             * CH4_conc_base_umol_L,N2O_conc_base_umol_L,AT_conc_base_ueq_L (canonical units)
             * NOTE: All species concentrations must be provided in µmol/L (GHG included).
             */
            double runoff_c;
            char is_polder[16];
            
            /* Try extended format first (with GHG columns) */
            parsed = sscanf(line, "%63[^,],%d,%lf,%lf,%lf,%15[^,],%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf",
                           branch_name, &segment_idx, &dist_km, &area_km2, &runoff_c, is_polder,
                           &Q_lat, &nh4_conc, &no3_conc, &po4_conc, &toc_conc, &dic_conc, &spm_conc,
                           &ch4_conc, &n2o_conc, &at_conc);
            
            if (parsed < 7) {
                fprintf(stderr, "Warning: Failed to parse lateral_sources.csv line %d (new format)\n", line_num);
                continue;
            }
            
            /* If we only parsed 13 fields, GHG columns are missing - that's OK */
            if (parsed < 16) {
                ch4_conc = 0;
                n2o_conc = 0;
                at_conc = 0;
            }
        } else {
            /* OLD FORMAT: Branch,Segment_Index,Distance_km,Area_km2,Q_lat_m3_s,
             * NH4_load_g_s,NO3_load_g_s,PO4_load_g_s,TOC_load_g_s,DIC_load_g_s,
             * NH4_conc_mg_L,NO3_conc_mg_L,PO4_conc_mg_L,TOC_conc_mg_L,DIC_conc_mg_L */
            double nh4_g_s, no3_g_s, po4_g_s, toc_g_s, dic_g_s;
            
            parsed = sscanf(line, "%63[^,],%d,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf",
                           branch_name, &segment_idx, &dist_km, &area_km2, &Q_lat,
                           &nh4_g_s, &no3_g_s, &po4_g_s, &toc_g_s, &dic_g_s,
                           &nh4_conc, &no3_conc, &po4_conc, &toc_conc, &dic_conc);
            
            if (parsed < 10) {
                fprintf(stderr, "Warning: Failed to parse lateral_sources.csv line %d (old format)\n", line_num);
                continue;
            }
        }
        
        /* Find the branch */
        trim_in_place(branch_name);
        Branch *b = find_branch_by_name(net, branch_name);
        if (!b) {
            /* Branch not in network - skip silently */
            continue;
        }
        
        /* Map distance to grid index
         * Grid indices go from 1 (downstream/mouth) to M (upstream)
         * Distance in CSV is from upstream (0 = upstream junction)
         * So: grid_idx = M - (dist_km * 1000 / dx)
         */
        int grid_idx = b->M - (int)(dist_km * 1000.0 / b->dx);
        if (grid_idx < 1) grid_idx = 1;
        if (grid_idx > b->M) grid_idx = b->M;
        
        /* Store lateral flow */
        if (b->lateral_flow) {
            b->lateral_flow[grid_idx] += Q_lat;
        }
        
        /* Store lateral concentrations
         * Convert from g/s to µmol/L (for species tracked in µmol/L)
         * 
         * Strategy: Store the CONCENTRATION in mg/L (from CSV conc columns)
         * The concentration will be mixed with the cell concentration 
         * based on flow ratio in biogeo.c
         */
        if (b->lateral_conc && Q_lat > 1e-10) {
            /* NH4: mg N/L → µmol N/L (×1000/14 = ×71.4) */
            b->lateral_conc[CGEM_SPECIES_NH4][grid_idx] = nh4_conc * 71.4;
            
            /* NO3: mg N/L → µmol N/L */
            b->lateral_conc[CGEM_SPECIES_NO3][grid_idx] = no3_conc * 71.4;
            
            /* PO4: mg P/L → µmol P/L (×1000/31 = ×32.3) */
            b->lateral_conc[CGEM_SPECIES_PO4][grid_idx] = po4_conc * 32.3;
            
            /* TOC: mg C/L → µmol C/L (×1000/12 = ×83.3) */
            b->lateral_conc[CGEM_SPECIES_TOC][grid_idx] = toc_conc * 83.3;
            
            /* DIC: mg C/L → µmol C/L */
            b->lateral_conc[CGEM_SPECIES_DIC][grid_idx] = dic_conc * 83.3;
            
            /* === GHG species === */
            if (ch4_conc > 0) {
                b->lateral_conc[CGEM_SPECIES_CH4][grid_idx] = ch4_conc;  /* Expect µmol/L input */
            }

            if (n2o_conc > 0) {
                b->lateral_conc[CGEM_SPECIES_N2O][grid_idx] = n2o_conc;  /* Expect µmol/L input */
            }
            
            /* AT (Total Alkalinity): Already in µeq/L in CSV - store directly */
            if (at_conc > 0) {
                b->lateral_conc[CGEM_SPECIES_AT][grid_idx] = at_conc;  /* µeq/L */
            }
        }
        
        total_loaded++;
        b->has_lateral_loads = 1;
    }
    
    fclose(fp);
    
    /* Print summary */
    int branches_with_loads = 0;
    double total_Q_lat = 0.0;
    double total_NH4_load = 0.0;
    
    for (size_t b = 0; b < net->num_branches; ++b) {
        Branch *br = net->branches[b];
        if (br->has_lateral_loads) {
            branches_with_loads++;
            for (int i = 1; i <= br->M; ++i) {
                total_Q_lat += br->lateral_flow[i];
                /* Estimate NH4 load from flow × concentration */
                total_NH4_load += br->lateral_flow[i] * 
                                  br->lateral_conc[CGEM_SPECIES_NH4][i] / 71.4;  /* Back to mg/s */
            }
        }
    }
    
    printf("  Loaded %d segments for %d branches\n", total_loaded, branches_with_loads);
    printf("  Total lateral Q: %.3f m³/s (%.0f m³/day)\n", total_Q_lat, total_Q_lat * 86400);
    printf("  Est. NH4 load: %.1f mg/s (%.1f kg/day)\n", total_NH4_load, total_NH4_load * 86.4);
    
    return 0;
}

/* ===========================================================================
 * LATERAL SEASONAL FACTORS LOADER
 * Reads rainfall-driven seasonal multipliers for Q and concentrations
 * 
 * CSV Format (lateral_seasonal_factors.csv):
 * Month,Month_Name,Season,Rain_mm,Q_Factor,NH4_Factor,NO3_Factor,PO4_Factor,TOC_Factor,SPM_Factor
 * 1,Jan,dry,15,1.0,1.0,1.0,1.0,1.0,1.0
 * ...
 * 12,Dec,dry,45,1.0,1.0,1.0,1.0,1.0,1.0
 * 
 * These factors multiply the base loads during simulation:
 *   Q_actual = Q_base × Q_Factor[month]
 *   C_actual = C_base × Species_Factor[month]
 * ===========================================================================*/

/**
 * Free memory allocated for seasonal factors
 */
void FreeLateralSeasonalFactors(LateralSeasonalFactors *factors) {
    if (!factors) return;
    
    if (factors->daily_Q_factor) { free(factors->daily_Q_factor); factors->daily_Q_factor = NULL; }
    if (factors->daily_NH4_factor) { free(factors->daily_NH4_factor); factors->daily_NH4_factor = NULL; }
    if (factors->daily_NO3_factor) { free(factors->daily_NO3_factor); factors->daily_NO3_factor = NULL; }
    if (factors->daily_PO4_factor) { free(factors->daily_PO4_factor); factors->daily_PO4_factor = NULL; }
    if (factors->daily_TOC_factor) { free(factors->daily_TOC_factor); factors->daily_TOC_factor = NULL; }
    if (factors->daily_SPM_factor) { free(factors->daily_SPM_factor); factors->daily_SPM_factor = NULL; }
    /* NEW: Free GHG factor arrays */
    if (factors->daily_CH4_factor) { free(factors->daily_CH4_factor); factors->daily_CH4_factor = NULL; }
    if (factors->daily_N2O_factor) { free(factors->daily_N2O_factor); factors->daily_N2O_factor = NULL; }
    if (factors->daily_DIC_factor) { free(factors->daily_DIC_factor); factors->daily_DIC_factor = NULL; }
    if (factors->daily_AT_factor) { free(factors->daily_AT_factor); factors->daily_AT_factor = NULL; }
    
    factors->num_days = 0;
    factors->use_daily = 0;
    factors->loaded = 0;
}

/**
 * Load monthly seasonal factors from CSV
 * 
 * @param net Network structure with lateral_factors field
 * @param case_dir Path to case directory
 * @return 0 on success, -1 on failure (model can continue with default factors)
 */
int LoadLateralSeasonalFactors(Network *net, const char *case_dir) {
    if (!net || !case_dir) return -1;
    
    LateralSeasonalFactors *factors = &net->lateral_factors;
    
    /* Initialize with default factors (all 1.0 = no seasonal variation) */
    for (int m = 0; m < CGEM_LATERAL_MAX_MONTHS; ++m) {
        factors->Q_factor[m] = 1.0;
        factors->NH4_factor[m] = 1.0;
        factors->NO3_factor[m] = 1.0;
        factors->PO4_factor[m] = 1.0;
        factors->TOC_factor[m] = 1.0;
        factors->SPM_factor[m] = 1.0;
        /* NEW: Initialize GHG factors */
        factors->CH4_factor[m] = 1.0;
        factors->N2O_factor[m] = 1.0;
        factors->DIC_factor[m] = 1.0;
        factors->AT_factor[m] = 1.0;
    }
    factors->use_daily = 0;
    factors->num_days = 0;
    factors->loaded = 0;
    /* Initialize daily pointers to NULL */
    factors->daily_CH4_factor = NULL;
    factors->daily_N2O_factor = NULL;
    factors->daily_DIC_factor = NULL;
    factors->daily_AT_factor = NULL;
    strncpy(factors->climate_preset, "default", sizeof(factors->climate_preset));
    
    /* Try to load monthly factors file */
    char path[CGEM_MAX_PATH];
    snprintf(path, sizeof(path), "%s/lateral_seasonal_factors.csv", case_dir);
    
    FILE *fp = fopen(path, "r");
    if (!fp) {
        printf("No lateral seasonal factors file found: %s\n", path);
        printf("  Using constant factors (Q=1.0, C=1.0 for all months)\n");
        factors->loaded = 1;  /* Mark as loaded (with defaults) */
        return 0;
    }
    
    printf("Loading lateral seasonal factors from: %s\n", path);
    
    char line[512];
    int line_num = 0;
    int months_loaded = 0;
    
    /* Read header line */
    if (!fgets(line, sizeof(line), fp)) {
        fclose(fp);
        return -1;
    }
    
    while (fgets(line, sizeof(line), fp)) {
        line_num++;
        
        int month;
        char month_name[16], season[16];
        double rain_mm, q_factor, nh4_factor, no3_factor, po4_factor, toc_factor, spm_factor;
        double ch4_factor = 1.0, n2o_factor = 1.0, dic_factor = 1.0, at_factor = 1.0;  /* NEW: GHG factors */
        
        /* Parse CSV line - try extended format first (with GHG columns) */
        int parsed = sscanf(line, "%d,%15[^,],%15[^,],%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf",
                           &month, month_name, season, &rain_mm,
                           &q_factor, &nh4_factor, &no3_factor, &po4_factor,
                           &toc_factor, &spm_factor,
                           &ch4_factor, &n2o_factor, &dic_factor, &at_factor);
        
        if (parsed < 10) {
            /* Try basic format without GHG */
            parsed = sscanf(line, "%d,%15[^,],%15[^,],%lf,%lf,%lf,%lf,%lf,%lf,%lf",
                           &month, month_name, season, &rain_mm,
                           &q_factor, &nh4_factor, &no3_factor, &po4_factor,
                           &toc_factor, &spm_factor);
            /* Use defaults for GHG factors */
            ch4_factor = 1.0;
            n2o_factor = 1.0;
            dic_factor = 1.0;
            at_factor = 1.0;
        }
        
        if (parsed < 6) {
            /* Try simpler format without SPM */
            parsed = sscanf(line, "%d,%15[^,],%15[^,],%lf,%lf,%lf,%lf,%lf,%lf",
                           &month, month_name, season, &rain_mm,
                           &q_factor, &nh4_factor, &no3_factor, &po4_factor, &toc_factor);
            spm_factor = toc_factor;  /* Default SPM = TOC factor */
        }
        
        if (parsed < 5) {
            fprintf(stderr, "Warning: Failed to parse seasonal factors line %d\n", line_num);
            continue;
        }
        
        /* Validate month (1-12) */
        if (month < 1 || month > 12) {
            fprintf(stderr, "Warning: Invalid month %d on line %d\n", month, line_num);
            continue;
        }
        
        /* Store factors (convert to 0-indexed) */
        int m = month - 1;
        factors->Q_factor[m] = q_factor;
        factors->NH4_factor[m] = nh4_factor;
        factors->NO3_factor[m] = no3_factor;
        factors->PO4_factor[m] = po4_factor;
        factors->TOC_factor[m] = toc_factor;
        factors->SPM_factor[m] = spm_factor;
        /* NEW: Store GHG factors */
        factors->CH4_factor[m] = ch4_factor;
        factors->N2O_factor[m] = n2o_factor;
        factors->DIC_factor[m] = dic_factor;
        factors->AT_factor[m] = at_factor;
        
        months_loaded++;
    }
    
    fclose(fp);
    
    if (months_loaded < 12) {
        fprintf(stderr, "Warning: Only loaded %d months (expected 12)\n", months_loaded);
    }
    
    factors->loaded = 1;
    
    /* Print summary */
    printf("  Loaded %d months of seasonal factors\n", months_loaded);
    printf("  Monthly Q factors: ");
    for (int m = 0; m < 12; ++m) {
        printf("%.1f ", factors->Q_factor[m]);
    }
    printf("\n");
    
    /* Optionally load daily factors for higher resolution */
    snprintf(path, sizeof(path), "%s/lateral_daily_factors.csv", case_dir);
    fp = fopen(path, "r");
    if (fp) {
        printf("  Loading daily factors from: %s\n", path);
        
        /* Count lines to allocate arrays */
        int num_lines = 0;
        while (fgets(line, sizeof(line), fp)) num_lines++;
        rewind(fp);
        
        /* Skip header */
        fgets(line, sizeof(line), fp);
        num_lines--;  /* Don't count header */
        
        if (num_lines > 0 && num_lines <= CGEM_LATERAL_MAX_DAYS) {
            /* Allocate daily arrays */
            factors->daily_Q_factor = (double *)malloc(num_lines * sizeof(double));
            factors->daily_NH4_factor = (double *)malloc(num_lines * sizeof(double));
            factors->daily_NO3_factor = (double *)malloc(num_lines * sizeof(double));
            factors->daily_PO4_factor = (double *)malloc(num_lines * sizeof(double));
            factors->daily_TOC_factor = (double *)malloc(num_lines * sizeof(double));
            factors->daily_SPM_factor = (double *)malloc(num_lines * sizeof(double));
            
            if (factors->daily_Q_factor && factors->daily_NH4_factor) {
                int day_idx = 0;
                while (fgets(line, sizeof(line), fp) && day_idx < num_lines) {
                    int day, doy, month;
                    char season[16];
                    double q_f, nh4_f, no3_f, po4_f, toc_f, spm_f;
                    
                    if (sscanf(line, "%d,%d,%d,%15[^,],%lf,%lf,%lf,%lf,%lf,%lf",
                              &day, &doy, &month, season,
                              &q_f, &nh4_f, &no3_f, &po4_f, &toc_f, &spm_f) >= 6) {
                        factors->daily_Q_factor[day_idx] = q_f;
                        factors->daily_NH4_factor[day_idx] = nh4_f;
                        factors->daily_NO3_factor[day_idx] = no3_f;
                        factors->daily_PO4_factor[day_idx] = po4_f;
                        factors->daily_TOC_factor[day_idx] = toc_f;
                        factors->daily_SPM_factor[day_idx] = spm_f > 0 ? spm_f : toc_f;
                        day_idx++;
                    }
                }
                factors->num_days = day_idx;
                factors->use_daily = 1;
                printf("  Loaded %d days of daily factors\n", day_idx);
            }
        }
        
        fclose(fp);
    }
    
    return 0;
}

/**
 * Get seasonal factor for a specific day and species
 * 
 * @param factors Loaded seasonal factors
 * @param day_of_year Day of year (0-364)
 * @param species Species index
 * @return Multiplication factor (default 1.0)
 */
double GetLateralFactor(LateralSeasonalFactors *factors, int day_of_year, int species) {
    if (!factors || !factors->loaded) return 1.0;
    
    /* Use daily factors if available */
    if (factors->use_daily && factors->daily_Q_factor && day_of_year < factors->num_days) {
        switch (species) {
            case CGEM_SPECIES_NH4: return factors->daily_NH4_factor[day_of_year];
            case CGEM_SPECIES_NO3: return factors->daily_NO3_factor[day_of_year];
            case CGEM_SPECIES_PO4: return factors->daily_PO4_factor[day_of_year];
            case CGEM_SPECIES_TOC: return factors->daily_TOC_factor[day_of_year];
            case CGEM_SPECIES_SPM: return factors->daily_SPM_factor[day_of_year];
            /* NEW: GHG species (December 2025) */
            case CGEM_SPECIES_CH4: return factors->daily_CH4_factor ? factors->daily_CH4_factor[day_of_year] : 1.0;
            case CGEM_SPECIES_N2O: return factors->daily_N2O_factor ? factors->daily_N2O_factor[day_of_year] : 1.0;
            case CGEM_SPECIES_DIC: return factors->daily_DIC_factor ? factors->daily_DIC_factor[day_of_year] : 1.0;
            case CGEM_SPECIES_AT:  return factors->daily_AT_factor ? factors->daily_AT_factor[day_of_year] : 1.0;
            case -1: return factors->daily_Q_factor[day_of_year];  /* -1 = Q factor */
            default: return 1.0;
        }
    }
    
    /* Fall back to monthly factors */
    int month = (day_of_year * 12) / 365;  /* Approximate month (0-11) */
    if (month < 0) month = 0;
    if (month > 11) month = 11;
    
    switch (species) {
        case CGEM_SPECIES_NH4: return factors->NH4_factor[month];
        case CGEM_SPECIES_NO3: return factors->NO3_factor[month];
        case CGEM_SPECIES_PO4: return factors->PO4_factor[month];
        case CGEM_SPECIES_TOC: return factors->TOC_factor[month];
        case CGEM_SPECIES_SPM: return factors->SPM_factor[month];
        /* NEW: GHG species (December 2025) */
        case CGEM_SPECIES_CH4: return factors->CH4_factor[month];
        case CGEM_SPECIES_N2O: return factors->N2O_factor[month];
        case CGEM_SPECIES_DIC: return factors->DIC_factor[month];
        case CGEM_SPECIES_AT:  return factors->AT_factor[month];
        case -1: return factors->Q_factor[month];  /* -1 = Q factor */
        default: return 1.0;
    }
}


/* ===========================================================================
 * POINT SOURCE LOADER
 * Loads major urban discharge points (cities, WWTP) from CSV
 * 
 * IMPORTANT: Point sources are ADDITIVE to lateral loads, but we need to 
 * avoid double-counting urban areas. When point_sources.csv is loaded:
 * 1. The point source flow and concentrations are added at specific cells
 * 2. Lateral loads from Urban land use within exclusion_radius_km are reduced
 * 
 * CSV Format (point_sources.csv):
 * Name,Branch,Segment_Index,Distance_km,Population,Treatment,Q_m3_s,NH4_mg_L,...
 * Can_Tho,Hau_River,40,80.0,1500000,primary,2.6,36.0,...
 * ===========================================================================*/

/**
 * Load point sources from CSV file
 * 
 * Point sources represent concentrated urban discharge (cities, WWTP) that are
 * separate from diffuse lateral loads. To avoid double-counting:
 * - Lateral loads within 5km of point sources are reduced by 50%
 * - This accounts for urban runoff already captured by sewage system
 * 
 * @param net Network structure
 * @param case_dir Path to case directory
 * @return Number of point sources loaded (0 if file not found, -1 on error)
 */
int LoadPointSources(Network *net, const char *case_dir) {
    if (!net || !case_dir) return -1;
    
    char path[CGEM_MAX_PATH];
    snprintf(path, sizeof(path), "%s/point_sources.csv", case_dir);
    
    FILE *fp = fopen(path, "r");
    if (!fp) {
        printf("No point sources file found: %s (continuing without point sources)\n", path);
        return 0;  /* Not an error - point sources are optional */
    }
    
    printf("Loading point sources from: %s\n", path);
    
    char line[1024];
    int line_num = 0;
    int total_loaded = 0;
    double total_Q = 0.0;
    
    /* Exclusion radius for reducing urban lateral loads near point sources */
    const double EXCLUSION_RADIUS_KM = 5.0;
    
    /* Read header line */
    if (!fgets(line, sizeof(line), fp)) {
        fclose(fp);
        return -1;
    }
    
    while (fgets(line, sizeof(line), fp)) {
        line_num++;
        
        /* Parse CSV line */
        char name[64], branch_name[64], treatment[32];
        int segment_idx;
        double dist_km;
        int population;
        double Q_m3_s;
        double nh4_mg_L = 0, no3_mg_L = 0, po4_mg_L = 0, toc_mg_L = 0, dic_mg_L = 0, spm_mg_L = 0;
        double ch4_umol_L = 0, n2o_umol_L = 0;  /* Inputs must be µmol/L */
        
        /* Parse: Name,Branch,Segment_Index,Distance_km,Population,Treatment,Q_m3_s,NH4,NO3,PO4,TOC,DIC,SPM,CH4,N2O */
        int parsed = sscanf(line, "%63[^,],%63[^,],%d,%lf,%d,%31[^,],%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf",
                   name, branch_name, &segment_idx, &dist_km, &population, treatment,
                   &Q_m3_s, &nh4_mg_L, &no3_mg_L, &po4_mg_L, &toc_mg_L, &dic_mg_L, &spm_mg_L,
                   &ch4_umol_L, &n2o_umol_L);
        
        if (parsed < 7) {
            fprintf(stderr, "Warning: Failed to parse point_sources.csv line %d\n", line_num);
            continue;
        }
        
        /* Find the branch */
        trim_in_place(branch_name);
        Branch *b = find_branch_by_name(net, branch_name);
        if (!b) {
            fprintf(stderr, "Warning: Point source '%s' branch '%s' not found\n", name, branch_name);
            continue;
        }
        
        /* Map distance to grid index (same logic as lateral loads) */
        int grid_idx = b->M - (int)(dist_km * 1000.0 / b->dx);
        if (grid_idx < 1) grid_idx = 1;
        if (grid_idx > b->M) grid_idx = b->M;
        
        /* Add point source flow to lateral flow array */
        if (b->lateral_flow) {
            b->lateral_flow[grid_idx] += Q_m3_s;
        }
        
        /* Add point source concentrations (weighted by flow if combined with lateral) */
        if (b->lateral_conc && Q_m3_s > 1e-10) {
            /* Get existing lateral flow at this cell */
            double existing_Q = b->lateral_flow[grid_idx] - Q_m3_s;  /* Subtract point source Q we just added */
            double total_Q_cell = b->lateral_flow[grid_idx];
            
            if (existing_Q > 1e-10) {
                /* Mix point source with existing lateral load (flow-weighted average) */
                double w_point = Q_m3_s / total_Q_cell;
                double w_lateral = existing_Q / total_Q_cell;
                
                /* NH4: mg N/L -> umol N/L (x71.4) */
                b->lateral_conc[CGEM_SPECIES_NH4][grid_idx] = 
                    w_lateral * b->lateral_conc[CGEM_SPECIES_NH4][grid_idx] + 
                    w_point * (nh4_mg_L * 71.4);
                
                /* NO3 */
                b->lateral_conc[CGEM_SPECIES_NO3][grid_idx] = 
                    w_lateral * b->lateral_conc[CGEM_SPECIES_NO3][grid_idx] + 
                    w_point * (no3_mg_L * 71.4);
                
                /* PO4: mg P/L -> umol P/L (x32.3) */
                b->lateral_conc[CGEM_SPECIES_PO4][grid_idx] = 
                    w_lateral * b->lateral_conc[CGEM_SPECIES_PO4][grid_idx] + 
                    w_point * (po4_mg_L * 32.3);
                
                /* TOC: mg C/L -> umol C/L (x83.3) */
                b->lateral_conc[CGEM_SPECIES_TOC][grid_idx] = 
                    w_lateral * b->lateral_conc[CGEM_SPECIES_TOC][grid_idx] + 
                    w_point * (toc_mg_L * 83.3);
                
                /* DIC */
                b->lateral_conc[CGEM_SPECIES_DIC][grid_idx] = 
                    w_lateral * b->lateral_conc[CGEM_SPECIES_DIC][grid_idx] + 
                    w_point * (dic_mg_L * 83.3);
                
                /* CH4/N2O: expect µmol/L inputs; no automatic conversion */
                if (ch4_umol_L > 0) {
                    b->lateral_conc[CGEM_SPECIES_CH4][grid_idx] = 
                        w_lateral * b->lateral_conc[CGEM_SPECIES_CH4][grid_idx] + 
                        w_point * ch4_umol_L;
                }
                
                if (n2o_umol_L > 0) {
                    b->lateral_conc[CGEM_SPECIES_N2O][grid_idx] = 
                        w_lateral * b->lateral_conc[CGEM_SPECIES_N2O][grid_idx] + 
                        w_point * n2o_umol_L;
                }
            } else {
                /* No existing lateral load - set directly */
                b->lateral_conc[CGEM_SPECIES_NH4][grid_idx] = nh4_mg_L * 71.4;
                b->lateral_conc[CGEM_SPECIES_NO3][grid_idx] = no3_mg_L * 71.4;
                b->lateral_conc[CGEM_SPECIES_PO4][grid_idx] = po4_mg_L * 32.3;
                b->lateral_conc[CGEM_SPECIES_TOC][grid_idx] = toc_mg_L * 83.3;
                b->lateral_conc[CGEM_SPECIES_DIC][grid_idx] = dic_mg_L * 83.3;
                /* GHG: Inputs must be µmol/L */
                if (ch4_umol_L > 0) b->lateral_conc[CGEM_SPECIES_CH4][grid_idx] = ch4_umol_L;
                if (n2o_umol_L > 0) b->lateral_conc[CGEM_SPECIES_N2O][grid_idx] = n2o_umol_L;
            }
        }
        
        /* Reduce urban lateral loads near this point source to avoid double-counting
         * This reduces Q_lat by 50% within EXCLUSION_RADIUS_KM of the point source
         * Rationale: Urban runoff near cities is partly collected by sewage systems
         */
        int exclusion_cells = (int)(EXCLUSION_RADIUS_KM * 1000.0 / b->dx);
        for (int i = grid_idx - exclusion_cells; i <= grid_idx + exclusion_cells; ++i) {
            if (i >= 1 && i <= b->M && i != grid_idx && b->lateral_flow) {
                /* Reduce lateral flow by 50% in exclusion zone */
                b->lateral_flow[i] *= 0.5;
            }
        }
        
        b->has_lateral_loads = 1;
        total_loaded++;
        total_Q += Q_m3_s;
        
        printf("  + %s (pop %d): Q=%.3f m3/s at %s km %.1f (cell %d)\n", 
               name, population, Q_m3_s, branch_name, dist_km, grid_idx);
    }
    
    fclose(fp);
    
    printf("  Loaded %d point sources (total Q: %.3f m3/s)\n", total_loaded, total_Q);
    if (total_loaded > 0) {
        printf("  Urban lateral loads reduced by 50%% within %.0f km of point sources\n", EXCLUSION_RADIUS_KM);
    }
    
    return total_loaded;
}


