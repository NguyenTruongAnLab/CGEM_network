# Project Structure

## Overview

```
CGEM_network/
├── src/                    # Source code
│   ├── *.c, *.h           # Core framework
│   ├── physics/           # Hydrodynamics + Transport
│   ├── rive/              # C-RIVE Biogeochemistry
│   ├── io/                # Input/Output
│   └── optimization/      # Calibration (NLopt)
├── INPUT/                 # Input data
│   └── Cases/             # Case configurations
├── OUTPUT/                # Simulation outputs
├── external/              # External libraries
│   └── nlopt/            # NLopt optimization library
├── scripts/               # Build and utility scripts
├── docs/                  # Documentation (MkDocs)
└── .github/               # GitHub workflows + Copilot instructions
```

## Source Code (`src/`)

### Core Framework

| File | Purpose |
|------|---------|
| `define.h` | Constants, species/reaction enums, macros |
| `network.h` | Core data structures: Branch, Node, Network, CaseConfig |
| `main.c` | Entry point, argument parsing |
| `init.c` | Network initialization orchestration |
| `network_data.c` | Memory allocation (`allocate_branch()`, `free_branch()`) |
| `network_simulation.c` | Main time-stepping loop |
| `utilities.c` | Helper functions |

### Physics Module (`src/physics/`)

| File | Purpose |
|------|---------|
| `hydrodynamics.c` | Saint-Venant solver, geometry, friction |
| `solver_hydro.c` | Network-level stepping, junction iteration |
| `transport.c` | Advection-dispersion, TVD schemes, Van den Burgh |

### C-RIVE Biogeochemistry (`src/rive/`)

| File | Purpose |
|------|---------|
| `rive.h` | Master include header |
| `rive_common.h` | Shared constants, temperature functions |
| `rive_params.h/.c` | Parameter loading from `biogeo_params.txt` |
| `biogeo.c` | Main biogeochemistry driver |
| `phytoplankton.h/.c` | PHY1 (diatoms), PHY2 (greens) |
| `nutrients.h/.c` | N/P/Si cycling |
| `oxygen.h/.c` | O₂ dynamics, reaeration |
| `carbonate_chem.h/.c` | DIC, TA, pH, pCO₂ |
| `ghg_module.h/.c` | N₂O, CH₄ emissions |
| `rk4_solver.h/.c` | Runge-Kutta 4 integrator |
| `sediment.c` | SPM erosion/deposition |

### Input/Output (`src/io/`)

| File | Purpose |
|------|---------|
| `config_parser.c` | Parse `case_config.txt` |
| `io_network.c` | Load topology.csv, boundary_map.csv |
| `io_manager.c` | Output coordination |
| `file.c` | Binary and CSV writers |
| `network_loader.c` | Loading utilities |

### Calibration (`src/optimization/`)

| File | Purpose |
|------|---------|
| `calibration.h` | Calibration structures and API |
| `calibration.c` | NLopt integration, objective functions, seasonal RMSE |

## Data Structures

### Network

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

### Branch

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
} Branch;
```

### Node

```c
typedef struct {
    int id;
    NodeType type;              // JUNCTION, DISCHARGE_BC, LEVEL_BC
    
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
    
    // Species mixing
    double *mixed_conc;
} Node;
```

## Grid Convention

C-GEM uses a **staggered grid** following Preissmann:

```
Index:    0     1     2     3     4     5    ...    M    M+1
        ghost                                      last  ghost
          │     │     │     │     │     │           │     │
          ▼     ▼     ▼     ▼     ▼     ▼           ▼     ▼
Grid:   ──┬─────┬─────┬─────┬─────┬─────┬── ... ───┬─────┬──
          │     │     │     │     │     │           │     │
          U     H     U     H     U     H           H     U
         vel  scal   vel  scal   vel  scal        scal   vel
```

- **Even indices (0, 2, 4...)**: Velocity points
- **Odd indices (1, 3, 5...)**: Scalar points (depth, concentration)
- **Index 1**: Downstream boundary (ocean)
- **Index M**: Upstream boundary (river)
- **Indices 0, M+1**: Ghost cells for boundary conditions

## Build System

### Windows (`scripts/build.bat`)

```batch
@echo off
set CC=C:\msys64\mingw64\bin\gcc.exe
set CFLAGS=-std=c11 -O2 -Wall -Wextra
set INCLUDES=-Isrc -Isrc/rive -Isrc/physics -Isrc/io -Iexternal/nlopt/src/api
set LDFLAGS=-Lexternal/nlopt/build -lnlopt -lm

# Compile all modules
gcc %CFLAGS% %INCLUDES% -c src/*.c src/physics/*.c src/rive/*.c src/io/*.c src/optimization/*.c

# Link
gcc -o bin/Debug/CGEM_Network.exe *.o %LDFLAGS%
```

### Scripts

| Script | Purpose |
|--------|---------|
| `build.bat` | Windows build script |
| `build-and-run.ps1` | Build + run helper |
| `bin_to_nc.py` | Convert binary output to NetCDF |
| `plot_netcdf.py` | Quick visualization |
| `generate_mekong_inputs.py` | Generate Mekong case inputs |

## Input/Output Flow

```
                            INPUT
                              │
        ┌─────────────────────┼─────────────────────┐
        ▼                     ▼                     ▼
   case_config.txt      topology.csv        boundary_map.csv
        │                     │                     │
        ▼                     ▼                     ▼
  ┌─────────────┐      ┌─────────────┐      ┌─────────────┐
  │ CaseConfig  │      │  Network    │      │   Nodes     │
  │             │◄─────│  Branches   │◄─────│   Forcing   │
  └─────────────┘      └─────────────┘      └─────────────┘
                              │
                              ▼
                    ┌─────────────────────┐
                    │    SIMULATION       │
                    │   network_run_      │
                    │   simulation()      │
                    └─────────────────────┘
                              │
        ┌─────────────────────┼─────────────────────┐
        ▼                     ▼                     ▼
   Binary (.bin)         CSV files            NetCDF (.nc)
   (raw output)      (per-branch)          (consolidated)
                            OUTPUT
```
