#include "network.h"
#include <stdio.h>

/*
 * Transitional loader. During Phase 2 this file will delegate to the new case
 * manager (config_parser + io_network). For now it simply advertises the
 * function signature expected by future code paths.
 */

int network_load_from_case(const char *case_config_path) {
    (void)case_config_path;
    fprintf(stderr, "network_load_from_case() is not implemented yet.\n");
    return -1;
}
