# C-GEM Network Source Code

## Overview

C-GEM Network is a **1D estuarine biogeochemical model** for multi-branch river networks.
The source code follows a modular architecture with a **single-level folder structure**:

```
src/
├── Core Framework (main loop, data structures, initialization)
├── physics/         - Hydrodynamics + Transport (Saint-Venant, advection-dispersion)
├── rive/            - C-RIVE biogeochemistry + sediment
└── io/              - Input/output, config parsing
```

---

## Architecture

### The Branch-Node System

C-GEM represents river networks as a **graph of Branches connected by Nodes**:

```
    [Node 1]                              [Node 4]
    (Discharge BC)                        (Level BC - Ocean)
        │                                      │
        ▼                                      ▼
    ═══════════════════════════════════════════════
         Branch 1: Main River                 
    ═══════════════════════════════════════════════
                      │
                      ▼
                  [Node 2]  ◄──── Junction (water/mass mixing)
                   /    \
                  ▼      ▼
    ═════════════      ═════════════
     Branch 2           Branch 3
    ═════════════      ═════════════
         │                  │
         ▼                  ▼
     [Node 5]           [Node 6]
     (Ocean)            (Ocean)
```

**Key concepts:**
- **Branch**: A river reach with M grid cells, solved as 1D domain
- **Node**: Connection point - can be junction, discharge BC, or level BC
- **Staggered Grid**: Velocity at even indices (0,2,4...), scalars at odd (1,3,5...)
- **Index 1 = downstream** (ocean), **Index M = upstream** (river source)

### Data Flow

```
┌─────────────────────────────────────────────────────────────────────┐
│                         Time Step Loop                               │
├─────────────────────────────────────────────────────────────────────┤
│  1. Hydrodynamics  │  Solve H, U for all branches with junction     │
│     (solver_hydro) │  iteration until network converges             │
├────────────────────┼────────────────────────────────────────────────┤
│  2. Transport      │  Advection-dispersion for all 17 species       │
│     (transport)    │  with TVD scheme                               │
├────────────────────┼────────────────────────────────────────────────┤
│  3. Sediment       │  SPM erosion/deposition based on shear stress  │
│     (sediment)     │                                                │
├────────────────────┼────────────────────────────────────────────────┤
│  4. Biogeochemistry│  RIVE reactions: photosynthesis, respiration,  │
│     (rive/biogeo)  │  nutrient cycling, carbonate chemistry         │
├────────────────────┼────────────────────────────────────────────────┤
│  5. Output         │  Write binary/CSV at specified intervals       │
│     (io/file)      │                                                │
└─────────────────────────────────────────────────────────────────────┘
```

---

## Directory Structure

```
src/
├── README.md               # This file
│
├── ─────────── CORE FRAMEWORK ───────────
├── define.h                # Constants, species indices (17), reaction indices (27)
├── network.h               # Branch, Node, Network, CaseConfig structures
├── main.c                  # Entry point, argument parsing
├── network_simulation.c    # Main time-stepping loop
├── init.c                  # Network initialization orchestration
├── network_data.c          # Memory allocation (allocate_branch, free_branch)
├── utilities.c             # Helper functions
│
├── ─────────── HYDRODYNAMICS ───────────
├── hydrodynamics/
│   ├── hydrodynamics.c     # Saint-Venant solver, geometry, tridiagonal
│   └── solver_hydro.c      # Network stepping, junction iteration
│
├── ─────────── TRANSPORT ───────────────
├── transport/
│   └── transport.c         # Advection-dispersion, Van den Burgh, TVD
│
├── ─────────── SEDIMENT ────────────────
├── sediment/
│   └── sediment.c          # SPM erosion/deposition (Krone/Partheniades)
│
├── ─────────── BIOGEOCHEMISTRY ─────────
├── rive/                   # C-RIVE module (see rive/README.md)
│   ├── rive.h, rive_common.h
│   ├── rive_params.h/c     # Parameter loading
│   ├── biogeo.c            # Main driver
│   ├── phytoplankton.h/c   # PHY1, PHY2
│   ├── nutrients.h/c       # N, P, Si cycling
│   ├── oxygen.h/c          # O2 dynamics
│   ├── carbonate_chem.h/c  # DIC, TA, pH, pCO2
│   ├── ghg_module.h/c      # N2O, CH4
│   └── rk4_solver.h/c      # Adaptive RK4
│
└── ─────────── INPUT/OUTPUT ────────────
    io/
    ├── config_parser.c     # Parse case_config.txt
    ├── io_network.c        # Load topology.csv, boundary_map.csv
    ├── io_manager.c        # Output orchestration
    ├── file.c              # Binary/CSV writers
    └── network_loader.c    # Network loading utilities
```

---

## Module Documentation

### 1. Core Framework (`src/`)

The core provides data structures and orchestrates the simulation:

| File | Purpose |
|------|---------|
| `define.h` | All constants, 17 species indices, 27 reaction indices, macros |
| `network.h` | `Branch`, `Node`, `Network`, `CaseConfig` structs |
| `main.c` | Entry point, creates output directories, calls simulation |
| `network_simulation.c` | Main loop: hydro → transport → sediment → biogeo → output |
| `init.c` | Orchestrates initialization of all subsystems |
| `network_data.c` | `allocate_branch()`, `free_branch()`, `allocate_node()` |

**Key Data Structures:**

```c
// Branch - a river reach with M grid cells
struct Branch {
    char name[64];
    int M;                    // Number of grid cells
    double dx;                // Grid spacing [m]
    double *depth, *velocity; // Hydrodynamic state [M+2]
    double **conc;            // Species concentrations [num_species][M+2]
    double **reaction_rates;  // Reaction rates [num_reactions][M+2]
    int node_up, node_down;   // Connected nodes
    // ... parameters for hydro, transport, biogeo
};

// Node - junction or boundary condition
struct Node {
    NodeType type;            // JUNCTION, DISCHARGE_BC, LEVEL_BC
    double H, Q_net;          // Water level, net discharge
    int connected_branches[8];
    double *mixed_conc;       // Junction mixing result
};
```

### 2. Hydrodynamics (`src/hydrodynamics/`)

Solves the **Saint-Venant equations** on a staggered grid:

| File | Purpose |
|------|---------|
| `hydrodynamics.c` | Geometry setup, single-branch solver, tridiagonal matrix |
| `solver_hydro.c` | Network-level iteration until junction mass balance converges |

**Physics:**
- Continuity: ∂A/∂t + ∂Q/∂x = 0
- Momentum: ∂U/∂t + U·∂U/∂x + g·∂H/∂x + friction = 0
- Friction: Manning/Chezy formulation
- Reference: Savenije (2012), Gisen et al. (2015)

**Junction iteration:**
```
Until converged:
  1. Solve each branch with boundary conditions from nodes
  2. Update junction water levels from mass balance
  3. Check convergence of H at all junctions
```

### 3. Transport (`src/physics/`)

Solves the **advection-dispersion equation** for all species:

| File | Purpose |
|------|---------|
| `transport.c` | TVD schemes, Van den Burgh dispersion, boundary conditions |

**Physics:**
- ∂C/∂t + U·∂C/∂x = ∂/∂x(K·∂C/∂x) + R
- Van den Burgh dispersion: K = K₀·exp(-x/L_d)
- TVD flux limiter for numerical stability
- Reference: Van den Burgh (1972), Savenije (2005)

### 4. Sediment (`src/rive/`)

**SPM erosion/deposition** based on bed shear stress:

| File | Purpose |
|------|---------|
| `sediment.c` | Krone/Partheniades formulation |

**Physics:**
- Erosion: E = M·(τ/τ_crit - 1) when τ > τ_crit
- Deposition: D = w_s·C·(1 - τ/τ_dep) when τ < τ_dep
- Shear stress: τ = ρ·g·U²/C²

### 5. Biogeochemistry (`src/rive/`)

**C-RIVE module** - comprehensive process-based biogeochemistry.
See [rive/README.md](rive/README.md) for detailed documentation.

**17 State Variables:**
- Salinity (conservative tracer)
- PHY1 (diatoms), PHY2 (green algae)
- DSi, NO3, NH4, PO4 (nutrients)
- O2 (dissolved oxygen)
- TOC (organic carbon)
- SPM (suspended matter)
- DIC, AT, pCO2, CO2, pH (carbonate system)
- HS, ALKC (diagnostic)

**27 Reaction Rates** tracked for carbon/mass budgets.

### 6. Input/Output (`src/io/`)

| File | Purpose |
|------|---------|
| `config_parser.c` | Parse `case_config.txt` (simulation settings) |
| `io_network.c` | Load `topology.csv` (branches), `boundary_map.csv` (forcing) |
| `io_manager.c` | Coordinate output timing |
| `file.c` | Write binary and CSV output files |

---

## How to Extend the Model

### Adding a New Species

1. **`define.h`**: Add to species enum and increment `CGEM_NUM_SPECIES`
   ```c
   #define CGEM_SPECIES_NEWVAR 17
   #define CGEM_NUM_SPECIES 18
   ```

2. **`init.c`**: Set boundary conditions in `init_species_bc()`

3. **`rive/biogeo.c`**: Add reactions affecting the species

4. **`io/file.c`**: Add to output variable list

### Adding a New Reaction

1. **`define.h`**: Add to reaction enum
   ```c
   #define CGEM_REACTION_NEW_PROC 27
   #define CGEM_NUM_REACTIONS 28
   ```

2. **`rive/biogeo.c`**: Calculate rate in `Biogeo_Branch()`
   ```c
   branch->reaction_rates[CGEM_REACTION_NEW_PROC][i] = rate;
   ```

### Adding a New Physics Module

1. Create new folder: `src/newmodule/`
2. Add header + implementation: `newmodule.h`, `newmodule.c`
3. Update `scripts/build.bat` to compile the new folder
4. Call from `network_simulation.c` in the time loop

---

## Build & Run

```powershell
# Build
./scripts/build.bat

# Run
./bin/Debug/CGEM_Network.exe INPUT/Cases/Tien_River/case_config.txt

# Or use VS Code tasks
Ctrl+Shift+B → Build CGEM Network
```

---

## References

- **Hydrodynamics**: Savenije (2012) "Salinity and Tides in Alluvial Estuaries"
- **Transport**: Van den Burgh (1972), Savenije (2005)
- **Biogeochemistry**: Wang et al. (2018), Hasanyar et al. (2022), Billen et al. (1994)
- **Sediment**: Krone (1962), Partheniades (1965)
