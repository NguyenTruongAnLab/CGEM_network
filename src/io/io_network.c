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

        /* Parse topology columns:
         * 0: ID
         * 1: Name
         * 2: NodeUp
         * 3: NodeDown
         * 4: Length_m
         * 5: Width_Up_m
         * 6: Width_Down_m
         * 7: Depth_m
         * 8: Chezy
         * 9: Group (optional, default 0)
         * 10: RS (optional, default 1.0) - Storage width ratio for mangroves
         * 11: BiogeoParams (optional) - Path to branch-specific biogeo params file
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
        char tmp_biogeo_path[CGEM_MAX_PATH] = {0}; /* Empty = use global defaults */
        char *token = strtok(text, ",");
        while (token && column < 12) {
            char *value = trim_in_place(token);
            switch (column) {
                case 0:
                    tmp_id = (int)strtol(value, NULL, 10);
                    break;
                case 1:
                    snprintf(tmp_name, sizeof(tmp_name), "%s", value);
                    break;
                case 2:
                    tmp_node_up = (int)strtol(value, NULL, 10) - 1;
                    break;
                case 3:
                    tmp_node_down = (int)strtol(value, NULL, 10) - 1;
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
        snprintf(branch->biogeo_params_path, sizeof(branch->biogeo_params_path), "%s", tmp_biogeo_path);
        
        /* Diagnostic output for configuration verification */
        printf("  Parsed branch %d: %s (RS=%.1f", tmp_id, tmp_name, tmp_storage_ratio);
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

    // Set node types and connections
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
        if (connections == 1) {
            if (i == 4) {  // Assuming node 4 is ocean
                net->nodes[i].type = NODE_LEVEL_BC;
            } else {
                net->nodes[i].type = NODE_DISCHARGE_BC;
            }
        } else if (connections > 1) {
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

        if (column < 3 || node_id < 1 || node_id > (int)net->num_nodes) {
            fprintf(stderr, "Invalid boundary row: node_id=%d, columns=%d\n", node_id, column);
            continue;
        }

        Node *node = &net->nodes[node_id - 1];

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
