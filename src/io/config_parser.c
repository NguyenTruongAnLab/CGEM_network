#include "../network.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
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

static int is_absolute_path(const char *path) {
    if (!path || !*path) {
        return 0;
    }
    if (path[0] == '/' || path[0] == '\\') {
        return 1;
    }
    if (strlen(path) > 1 && path[1] == ':' && ((path[0] >= 'A' && path[0] <= 'Z') || (path[0] >= 'a' && path[0] <= 'z'))) {
        return 1;
    }
    return 0;
}

static void join_paths(const char *base, const char *relative, char *output, size_t output_len) {
    if (!output || output_len == 0) {
        return;
    }
    if (!relative || !*relative) {
        output[0] = '\0';
        return;
    }
    if (is_absolute_path(relative) || !base || !*base) {
        snprintf(output, output_len, "%s", relative);
        return;
    }

    size_t base_len = strlen(base);
    int needs_sep = !(base_len > 0 && (base[base_len - 1] == '/' || base[base_len - 1] == '\\'));
    if (needs_sep) {
        snprintf(output, output_len, "%s/%s", base, relative);
    } else {
        snprintf(output, output_len, "%s%s", base, relative);
    }
}

// No external path utilities required - keep path building simple.

int LoadCaseConfig(const char *path, CaseConfig *config) {
    if (!path || !config) {
        return -1;
    }

    FILE *fp = fopen(path, "r");
    if (!fp) {
        fprintf(stderr, "Failed to open case configuration: %s\n", path);
        return -2;
    }

    memset(config, 0, sizeof(*config));
    config->duration_days = 365;
    config->warmup_days = 0;
    config->write_csv = 0;
    config->write_netcdf = 0;
    config->write_reaction_rates = 0;
    config->reaction_mode = 1;  /* Default: reactions ON */
    config->dx_meters = CGEM_DEFAULT_DX_METERS;
    snprintf(config->start_date, sizeof(config->start_date), "1970-01-01");

    snprintf(config->config_path, sizeof(config->config_path), "%s", path);

    const char *last_slash = strrchr(path, '/');
    const char *last_backslash = strrchr(path, '\\');
    const char *split = last_slash;
    if (!split || (last_backslash && last_backslash > split)) {
        split = last_backslash;
    }
    if (split) {
        size_t len = (size_t)(split - path);
        if (len >= sizeof(config->base_dir)) {
            len = sizeof(config->base_dir) - 1;
        }
        memcpy(config->base_dir, path, len);
        config->base_dir[len] = '\0';
    } else {
        snprintf(config->base_dir, sizeof(config->base_dir), ".");
    }

    char line[1024];
    int dt_set = 0;
    int dx_set = 0;
    while (fgets(line, sizeof(line), fp)) {
        char *text = trim_in_place(line);
        if (*text == '\0' || *text == '#' || *text == ';') {
            continue;
        }

        char *equals = strchr(text, '=');
        if (!equals) {
            continue;
        }
        *equals = '\0';
        char *key = trim_in_place(text);
        char *value = trim_in_place(equals + 1);

        if (strcmp(key, "CaseName") == 0) {
            snprintf(config->case_name, sizeof(config->case_name), "%s", value);
        } else if (strcmp(key, "Topology") == 0) {
            join_paths(config->base_dir, value, config->topology_path, sizeof(config->topology_path));
        } else if (strcmp(key, "Boundaries") == 0 || strcmp(key, "BoundaryMap") == 0) {
            join_paths(config->base_dir, value, config->boundary_path, sizeof(config->boundary_path));
        } else if (strcmp(key, "TimeStepSeconds") == 0 || strcmp(key, "TimeStep") == 0 || strcmp(key, "DELTI") == 0) {
            double new_dt = strtod(value, NULL);
            if (dt_set && fabs(new_dt - config->dt) > 1e-6) {
                fprintf(stderr, "Warning: multiple TimeStep/DELTI keys in config; using latest value %.2f (overrides %.2f)\n", new_dt, config->dt);
            }
            config->dt = new_dt;
            dt_set = 1;
        } else if (strcmp(key, "DELXI") == 0) {
            double new_dx = strtod(value, NULL);
            if (dx_set && fabs(new_dx - config->dx_meters) > 1e-6) {
                fprintf(stderr, "Warning: multiple DELXI keys in config; using latest value %.2f (overrides %.2f)\n", new_dx, config->dx_meters);
            }
            config->dx_meters = new_dx;
            dx_set = 1;
        } else if (strcmp(key, "TidalAmplitude") == 0) {
            config->tidal_amplitude = strtod(value, NULL);
        } else if (strcmp(key, "RiverDischarge") == 0 || strcmp(key, "Q_River") == 0) {
            config->Q_river = strtod(value, NULL);
        } else if (strcmp(key, "Chezy_Group1") == 0) {
            // For now, store in a way that can be used later
            // This would need to be passed to branches during initialization
        } else if (strcmp(key, "Chezy_Group2") == 0) {
            // For now, store in a way that can be used later
        } else if (strcmp(key, "C_VDB") == 0) {
            // For now, store in a way that can be used later
        } else if (strcmp(key, "D0_Correction") == 0) {
            // For now, store in a way that can be used later
        } else if (strcmp(key, "StartDate") == 0) {
            snprintf(config->start_date, sizeof(config->start_date), "%s", value);
        } else if (strcmp(key, "Duration") == 0 || strcmp(key, "Duration_Days") == 0) {
            config->duration_days = (int)strtol(value, NULL, 10);
        } else if (strcmp(key, "Warmup") == 0 || strcmp(key, "Warmup_Days") == 0) {
            config->warmup_days = (int)strtol(value, NULL, 10);
        } else if (strcmp(key, "WriteCSV") == 0) {
            config->write_csv = (int)strtol(value, NULL, 10);
        } else if (strcmp(key, "WriteNetCDF") == 0) {
            config->write_netcdf = (int)strtol(value, NULL, 10);
        } else if (strcmp(key, "WriteReactionRates") == 0) {
            config->write_reaction_rates = (int)strtol(value, NULL, 10);
        } else if (strcmp(key, "ReactionMode") == 0) {
            /* ReactionMode = ON (1) or OFF (0) - allows disabling biogeochemistry
             * for testing transport, lateral sources, and boundary conditions */
            if (strcasecmp(value, "ON") == 0 || strcmp(value, "1") == 0) {
                config->reaction_mode = 1;
            } else if (strcasecmp(value, "OFF") == 0 || strcmp(value, "0") == 0) {
                config->reaction_mode = 0;
            } else {
                fprintf(stderr, "Warning: Unknown ReactionMode '%s', defaulting to ON\n", value);
                config->reaction_mode = 1;
            }
        } else if (strcmp(key, "OutputDir") == 0) {
            /* Allow explicit output directory override */
            if (value[0] == '/' || value[0] == '\\' || (strlen(value) > 1 && value[1] == ':')) {
                /* Absolute path */
                snprintf(config->output_dir, sizeof(config->output_dir), "%s", value);
            } else {
                /* Relative path - resolve from config file directory */
                join_paths(config->base_dir, value, config->output_dir, sizeof(config->output_dir));
            }
        }
    }

    fclose(fp);

    if (config->case_name[0] == '\0') {
        const char *basename = split ? split + 1 : path;
        snprintf(config->case_name, sizeof(config->case_name), "%s", basename);
        char *dot = strrchr(config->case_name, '.');
        if (dot) {
            *dot = '\0';
        }
    }

    // Always set output directory based on case name (stored under OUTPUT/)
    {
        // Get the project root directory (go up from INPUT/Cases/CaseName/)
        char project_root[2048];
        join_paths(config->base_dir, "../../..", project_root, sizeof(project_root));
        
        
        char output_path[CGEM_MAX_PATH];
        char tmp_path[CGEM_MAX_PATH];
        join_paths(".", "OUTPUT", tmp_path, sizeof(tmp_path));
        // Append case name: OUTPUT/<case_name>
        join_paths(tmp_path, config->case_name, output_path, sizeof(output_path));
        
        // Normalize path separators for Windows
        #ifdef _WIN32
        for (char *p = output_path; *p; ++p) {
            if (*p == '/') *p = '\\';
        }
        #endif
        
        size_t len_out = strlen(output_path);
        if (len_out >= sizeof(config->output_dir)) {
            fprintf(stderr, "Warning: output path longer than %zu characters, truncating: %s\n", sizeof(config->output_dir) - 1, output_path);
            len_out = sizeof(config->output_dir) - 1;
        }
        memcpy(config->output_dir, output_path, len_out);
        config->output_dir[len_out] = '\0';
    }

    return 0;
}
