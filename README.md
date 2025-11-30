# C-GEM Network — Carbon-Generic Estuary Model for Networks

## What is C-GEM?

C-GEM Network is a **specialized 1D estuarine biogeochemical model** designed for multi-branch tidal river networks (deltas, distributaries, confluences). It combines:

- **Saint-Venant hydrodynamics** on a staggered grid with Savenije's theory
- **TVD advection-dispersion transport** for scalar species
- **C-RIVE biogeochemistry** (Unified RIVE v1.0) — the state-of-the-art carbon, nutrient, and oxygen cycle model
- **Greenhouse gas emissions** (CO2, CH4, N2O) — critical for climate impact assessment
- **Sediment erosion/deposition** based on shear stress

> **Why another model?** While powerful 2D/3D models exist (TELEMAC, Delft3D, MIKE, HEC-RAS), C-GEM fills a critical niche for **rapid, process-based biogeochemical research** in deltaic networks. It runs 100-1000× faster than 2D models, enabling sensitivity analysis, parameter optimization, and scenario testing. See [Introduction: Why C-GEM?](docs/introduction.md) for the full rationale.

## 📚 Documentation

| Document | Description |
|----------|-------------|
| **[Introduction](docs/introduction.md)** | Why C-GEM exists, comparison with TELEMAC/Delft3D/MIKE |
| **[Hydrodynamics](docs/hydrodynamics.md)** | Saint-Venant equations, staggered grid, tidal propagation |
| **[Transport](docs/transport.md)** | Advection-dispersion, TVD schemes, Savenije theory |
| **[Biogeochemistry](docs/biogeochemistry.md)** | C-RIVE module, carbonate chemistry, GHG emissions |
| **[Data Requirements](docs/data_requirements.md)** | Input preparation, forcing data, calibration |
| **[Quick Start (Windows)](docs/QUICKSTART_WINDOWS.md)** | Installation and first run |

## Key Features

✅ **Multi-branch network topology** — bifurcations, confluences, distributaries  
✅ **Computationally efficient** — 1D, runs in seconds to minutes  
✅ **Complete carbon cycle** — DIC, TA, pH, pCO2, CO2 air-water flux  
✅ **Greenhouse gas emissions** — CO2, CH4, N2O with process attribution  
✅ **C-RIVE biogeochemistry** — 6-pool organic matter, explicit bacteria, 2-step nitrification  
✅ **Process-based** — mechanistic equations, not empirical correlations  
✅ **Open source** — ANSI C, portable, transparent

## Quick Start

```powershell
# Windows: Build and run
./scripts/build-and-run.ps1 -r Tien_River

# Or step-by-step
./scripts/build.bat
./bin/Debug/CGEM_Network.exe INPUT/Cases/Tien_River/case_config.txt
```

## Project Overview

C-GEM Network simulates:
- **Hydrodynamics**: Saint-Venant equations on a staggered grid (Savenije theory)
- **Transport**: Advection-dispersion with TVD schemes for scalar species
- **Biogeochemistry**: C-RIVE water quality (phytoplankton, nutrients, oxygen, carbon, GHG)
- **Sediment**: Erosion/deposition based on shear stress

The model supports **multi-branch network topologies** (distributary deltas, confluences) with explicit junction mass balance.

## Architecture

### Data Structures (`src/network.h`)

```
Network
├── branches[]          # Array of Branch pointers
├── nodes[]             # Array of Node structs  
├── num_branches, num_nodes, num_species
└── Simulation params (dt, dx_target, total_time, etc.)

Branch                  # Per-branch state
├── Geometry: id, name, node_up, node_down, length_m, width_up/down, depth, chezy, lc_convergence
├── Grid: M (cells), dx
├── Hydro arrays [0..M+1]: depth[], velocity[], waterLevel[], totalArea[], width[], dispersion[]
├── Species arrays: conc[num_species][0..M+1]
├── Reactions: reaction_rates[num_reactions][0..M+1]
├── Tridiagonal solver: tri_lower[], tri_diag[], tri_upper[], tri_rhs[]
└── Output: bin_fp, csv_fps[]

Node                    # Junction/boundary point
├── id, type (JUNCTION, DISCHARGE_BC, LEVEL_BC)
├── H, H_new, Q_net
├── connected_branches[], connection_dir[], num_connections
├── forcing_time[], forcing_value[], forcing_len
├── species_forcing_time[][], species_forcing_value[][]
└── mixed_conc[]        # Junction mixing result
```

### Grid Convention
- **Staggered grid**: velocity at even indices (0,2,4...), scalars at odd indices (1,3,5...)
- **Index 1 = downstream** (ocean side), Index M = upstream
- **Ghost cells**: 0 and M+1 for boundary conditions

### Module Files

| File | Purpose |
|------|---------|
| `src/network.h` | Core data structures (Branch, Node, Network, CaseConfig) |
| `src/define.h` | Constants, species/reaction indices, macros |
| `src/main.c` | Entry point, argument parsing |
| `src/config_parser.c` | Parse `case_config.txt` |
| `src/io_network.c` | Load topology.csv, boundary_map.csv |
| `src/network_data.c` | Memory allocation (allocate_branch, free_branch) |
| `src/init.c` | Network initialization orchestration |
| `src/hydrodynamics.c` | Branch geometry, Saint-Venant solver |
| `src/solver_hydro.c` | Network-level stepping, junction iteration |
| `src/transport.c` | Advection-dispersion, Van den Burgh dispersion |
| `src/biogeo.c` | Biogeochemical reactions (NPP, respiration, nitrification, CO2) |
| `src/sediment.c` | SPM erosion/deposition |
| `src/file.c` | Binary and CSV output |
| `src/network_simulation.c` | Main simulation loop, output orchestration |

## Output Variables

The model outputs **50 variables** per branch (6 hydro + 17 species + 27 reactions):

### Hydrodynamic (6)
| Index | Name | Description | Unit |
|-------|------|-------------|------|
| 0 | depth | Water depth H | m |
| 1 | velocity | Current velocity U | m/s |
| 2 | waterlevel | Water surface elevation | m |
| 3 | area | Cross-sectional area A | m² |
| 4 | width | Channel width B | m |
| 5 | dispersion | Dispersion coefficient K | m²/s |

### Species (20)
| Index | Name | Description | Unit |
|-------|------|-------------|------|
| 0 | salinity | Salinity | PSU |
| 1 | phy1 | Phytoplankton (diatoms) | mgC/L |
| 2 | phy2 | Phytoplankton (non-siliceous) | mgC/L |
| 3 | dsi | Dissolved silica | µmol/L |
| 4 | no3 | Nitrate | µmol/L |
| 5 | nh4 | Ammonium | µmol/L |
| 6 | po4 | Phosphate | µmol/L |
| 7 | o2 | Dissolved oxygen | µmol/L |
| 8 | toc | Total organic carbon | µmol/L |
| 9 | spm | Suspended particulate matter | mg/L |
| 10 | dic | Dissolved inorganic carbon | µmol/L |
| 11 | at | Total alkalinity | µmol/L |
| 12 | pco2 | Partial pressure CO2 (diagnostic) | µatm |
| 13 | co2 | Aqueous CO2 (diagnostic) | µmol/L |
| 14 | ph | pH (diagnostic) | - |
| 15 | hs | Sulfide | µmol/L |
| 16 | alkc | Carbonate alkalinity iteration (diagnostic) | - |
| 27 | no2 | Nitrite (2-step nitrification) | µmol N/L |
| 28 | n2o | Nitrous oxide (GHG) | nmol N/L |
| 29 | ch4 | Methane (GHG) | µmol C/L |

### Reactions (35+)
| Index | Name | Description |
|-------|------|-------------|
| 0-2 | npp_no3, npp_no3_1, npp_no3_2 | NO3-based NPP (total, phy1, phy2) |
| 3-5 | npp_nh4, npp_nh4_1, npp_nh4_2 | NH4-based NPP (total, phy1, phy2) |
| 6-7 | gpp_1, gpp_2 | Gross primary production |
| 8 | npp | Net primary production |
| 9-11 | phy_death, phy_death_1, phy_death_2 | Phytoplankton mortality |
| 12 | si_cons | Silica consumption |
| 13 | aer_deg | Aerobic degradation |
| 14 | denit | Denitrification |
| 15 | nit | Nitrification (NH4→NO2, step 1) |
| 16-17 | o2_ex, o2_ex_s | O2 exchange (volumetric, surface) |
| 18-19 | co2_ex, co2_ex_s | CO2 exchange (volumetric, surface) |
| 20 | dic_react | DIC reactions |
| 21 | ta_react | Alkalinity reactions |
| 22 | hs_react | Sulfide reactions |
| 23-24 | erosion_s, erosion_v | Erosion (surface, volumetric) |
| 25-26 | deposition_s, deposition_v | Deposition (surface, volumetric) |
| 40 | nit2 | Nitratation (NO2→NO3, step 2) |
| 41-43 | n2o_nit, n2o_denit, n2o_ex | N2O production and exchange |
| 44-48 | ch4_prod, ch4_ox, ch4_ox_anaer, ch4_ex, ch4_ebul | CH4 dynamics |

## Case Configuration

Cases are defined in `INPUT/Cases/<CaseName>/`:
- `case_config.txt`: Main configuration (simulation period, parameters)
- `topology.csv`: Branch definitions (ID, name, nodes, geometry)
- `boundary_map.csv`: Node boundary conditions (discharge/level files)
- `biogeo_params.txt`: Biogeochemical parameters
- `forcing_data/`: Time series for boundaries and species

## Build & Run

```powershell
# Build and run (Windows/PowerShell)
./scripts/build-and-run.ps1 -r Tien_River

# Or build only
./scripts/build.bat

# Then run
./bin/Debug/CGEM_Network.exe INPUT/Cases/Tien_River/case_config.txt
```

## Citation

If you use C-GEM Network in your research, please cite the following:

**For C-RIVE biogeochemistry:**
```bibtex
@article{wang2018crive,
  title={Time-dependent global sensitivity analysis of the {C-RIVE} biogeochemical model},
  author={Wang, Shuaitao and Flipo, Nicolas and Romary, Thomas},
  journal={Water Research},
  volume={144},
  pages={341--355},
  year={2018},
  publisher={Elsevier}
}

@article{hasanyar2022unified,
  title={Unified {RIVE} v1.0: a multi-environment aquatic biogeochemical model},
  author={Hasanyar, Masihullah and Flipo, Nicolas and Vilmin, Lauriane and Wang, Shuaitao},
  journal={Biogeosciences},
  year={2022}
}
```

**For Savenije's estuarine theory:**
```bibtex
@book{savenije2005salinity,
  title={Salinity and Tides in Alluvial Estuaries},
  author={Savenije, Hubert H.G.},
  year={2005},
  publisher={Elsevier}
}
```

## License

[To be determined - suggest MIT or GPL-3.0]

## Contributing

Contributions welcome! Please read the [implementation guidelines](.github/memory-bank/implementationGuidelines.md) before submitting pull requests.