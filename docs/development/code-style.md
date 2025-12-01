# Code Style

## General Principles

1. **Readability** over cleverness
2. **Config-driven** — no hardcoded values
3. **Memory safety** — allocate at init, free at shutdown
4. **Documentation** — comment all functions

## Naming Conventions

### Functions

```c
// Module prefix + action + target
int Hyd_Branch(Branch *b, ...);           // Hydrodynamics module
int Transport_Branch(Branch *b, ...);     // Transport module
int Biogeo_Branch(Branch *b, ...);        // Biogeochemistry module

// Helper functions: lowercase with underscores
double calc_bed_shear(Branch *b, int i);
void solve_tridiagonal(double *a, double *b, double *c, double *d, int n);
```

### Variables

```c
// Local variables: lowercase with underscores
double water_depth;
int cell_index;

// Loop indices: single letter
for (int i = 0; i < n; i++) { }

// Constants: UPPER_SNAKE_CASE
#define MAX_SPECIES 50
#define GRAVITY 9.81
```

### Structs

```c
// PascalCase for struct names
typedef struct {
    int id;
    char name[64];
    double *depth;
} Branch;

typedef struct {
    Branch **branches;
    size_t num_branches;
} Network;
```

## File Organization

```c
/**
 * @file module_name.c
 * @brief Brief description
 */

#include "module_name.h"
#include <stdlib.h>
#include <math.h>

/* Constants */
#define LOCAL_CONST 1.0

/* Static helper functions */
static double helper_function(double x) {
    return x * 2.0;
}

/* Public API functions */
int Module_PublicFunction(Branch *b) {
    // Implementation
}
```

## Documentation

### Function Headers

```c
/**
 * @brief Calculate erosion rate based on shear stress
 * 
 * Uses Partheniades formula for cohesive sediment erosion.
 * 
 * @param tau Bed shear stress [Pa]
 * @param tau_crit Critical erosion stress [Pa]
 * @param M_ero Erosion coefficient [kg/m²/s]
 * @return Erosion rate [kg/m²/s], 0 if tau < tau_crit
 * 
 * @note Reference: Partheniades (1965)
 */
double calc_erosion_rate(double tau, double tau_crit, double M_ero);
```

## Memory Management

### Allocation

```c
// Allocate at initialization
Branch *allocate_branch(int M, int num_species) {
    Branch *b = malloc(sizeof(Branch));
    if (!b) return NULL;
    
    b->depth = malloc((M + 2) * sizeof(double));
    if (!b->depth) {
        free(b);
        return NULL;
    }
    
    return b;
}
```

### Deallocation

```c
// Free at shutdown
void free_branch(Branch *b) {
    if (!b) return;
    
    free(b->depth);
    free(b->velocity);
    // ... free all arrays
    
    free(b);
}
```

### Rules

- Never `malloc` in time loops
- Always check return value of `malloc`
- Set pointers to `NULL` after `free`
- Use pre-allocated scratch arrays

## Error Handling

```c
int some_function(Branch *b) {
    if (!b) {
        fprintf(stderr, "Error: NULL branch pointer\n");
        return -1;
    }
    
    if (b->M <= 0) {
        fprintf(stderr, "Error: Invalid grid size M=%d\n", b->M);
        return -1;
    }
    
    // ... computation
    
    return 0;  // Success
}
```

## Loop Patterns

### Grid Loops

```c
// Standard grid loop (1 to M, excluding ghost cells)
for (int i = 1; i <= b->M; ++i) {
    double val = b->conc[species][i];
}

// Include ghost cells
for (int i = 0; i <= b->M + 1; ++i) {
    // ...
}
```

### Species Loops

```c
for (int s = 0; s < net->num_species; ++s) {
    Transport_Branch(b, s, dt, C_up[s], C_down[s]);
}
```

## Testing

Before committing:

- [ ] Compiles without warnings
- [ ] Mekong_Delta_Full case runs successfully
- [ ] Calibration mode works
- [ ] Memory leak check (valgrind on Linux)
