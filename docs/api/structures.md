# Core Data Structures

## Network

The top-level container for all simulation data.

```c
typedef struct {
    Branch **branches;          // Array of branch pointers
    Node *nodes;                // Array of nodes
    size_t num_branches;
    size_t num_nodes;
    int num_species;
    
    // Simulation parameters
    double dt;                  // Time step [s]
    double dx_target;           // Target grid spacing [m]
    double total_time;          // Total simulation time [s]
    double warmup_time;         // Warmup period [s]
    
    // Forcing
    double tidal_amplitude;
    double Q_river;
} Network;
```

## Branch

Represents a river reach with M grid cells.

```c
typedef struct {
    // Identity
    int id;
    char name[64];
    int node_up, node_down;
    
    // Geometry
    double length_m;
    double width_down_m, width_up_m;
    double depth_m;
    double chezy;
    double lc_convergence;      // Exponential convergence length [m]
    double vdb_coef;            // Van den Burgh coefficient
    double rs;                  // Storage width ratio
    
    // Grid
    int M;                      // Number of cells
    double dx;                  // Cell size [m]
    
    // State arrays [0..M+1] (includes ghost cells)
    double *depth;
    double *velocity;
    double *waterLevel;
    double *totalArea;
    double *width;
    double *dispersion;
    
    // Species [num_species][0..M+1]
    double **conc;
    
    // Reactions [num_reactions][0..M+1]
    double **reaction_rates;
    
    // Biogeochemistry parameters
    BiogeoParams *biogeo;
    
    // Tridiagonal solver workspace
    double *tri_lower, *tri_diag, *tri_upper, *tri_rhs;
} Branch;
```

## Node

Connection point in the network.

```c
typedef enum {
    NODE_JUNCTION,
    NODE_DISCHARGE_BC,
    NODE_LEVEL_BC
} NodeType;

typedef struct {
    int id;
    NodeType type;
    
    // State
    double H, H_new;            // Water level [m]
    double Q_net;               // Net discharge [m³/s]
    
    // Connectivity
    int *connected_branches;
    int *connection_dir;        // +1 = outflow, -1 = inflow
    int num_connections;
    
    // Forcing
    double *forcing_time;
    double *forcing_value;
    int forcing_len;
    
    // Species forcing
    double **species_forcing_time;
    double **species_forcing_value;
    int *species_forcing_len;
    
    // Species mixing
    double *mixed_conc;
} Node;
```

## CaseConfig

Simulation configuration.

```c
typedef struct {
    char case_name[256];
    char output_dir[512];
    
    // Time
    double duration;            // [days]
    double warmup;              // [days]
    double dt;                  // [s]
    
    // Numerics
    double dx_target;           // [m]
    
    // Output
    int write_csv;
    int write_netcdf;
    int write_reactions;
    
    // Calibration
    int calibration_mode;
    int calibration_stage;
    int max_iterations;
} CaseConfig;
```

## BiogeoParams

Biogeochemistry parameters.

```c
typedef struct {
    // Phytoplankton
    double pbmax1, pbmax2;      // Max photosynthesis rate [d-1]
    double kmort1, kmort2;      // Mortality rate [d-1]
    double kn, kp, ksi;         // Half-saturation constants
    
    // Bacteria
    double kox;                 // Aerobic degradation [d-1]
    double knit;                // Nitrification [d-1]
    double kdenit;              // Denitrification [d-1]
    
    // Oxygen
    double ko2;                 // O2 half-saturation [µmol/L]
    
    // Sediment
    double ws_base;             // Base settling velocity [m/s]
    double tau_ero;             // Critical erosion stress [Pa]
    double tau_dep;             // Critical deposition stress [Pa]
    double floc_sal_scale;      // Flocculation salinity scale [PSU]
    double floc_factor_max;     // Max flocculation factor [-]
    
    // Salinity stress
    double sal_stress_thresh;   // Stress onset [PSU]
    double sal_stress_coef;     // Mortality enhancement [-]
} BiogeoParams;
```

## Memory Management

### Allocation

```c
Branch *allocate_branch(int M, int num_species, int num_reactions);
Node *allocate_node(int num_species);
Network *allocate_network(int num_branches, int num_nodes, int num_species);
```

### Deallocation

```c
void free_branch(Branch *b);
void free_node(Node *n);
void free_network(Network *net);
```
