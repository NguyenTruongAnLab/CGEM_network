# Introduction

## The Challenge of Delta Biogeochemistry

Tropical deltas like the **Mekong**, **Ganges-Brahmaputra**, and **Niger** face urgent environmental challenges:

- **Eutrophication** from agricultural runoff and urbanization
- **Hypoxia** threatening fisheries and aquaculture
- **Greenhouse gas emissions** (CO₂, CH₄, N₂O) contributing to climate change
- **Salt intrusion** advancing further upstream each dry season
- **Sediment starvation** from upstream damming

Understanding these processes requires models that capture both the **hydrodynamics** of complex tidal networks and the **biogeochemical transformations** in the water column.

## Why C-GEM?

### The Problem with Existing Models

| Model | Strengths | Limitations |
|-------|-----------|-------------|
| **TELEMAC-2D/3D** | Unstructured grids, morphodynamics | Complex setup (weeks), expensive computation |
| **Delft3D-FLOW/WAQ** | Comprehensive WQ | Steep learning curve, commercial license |
| **MIKE 11/21** | User-friendly GUI | Expensive licensing (~€10k/year) |
| **HEC-RAS** | Free, flood modeling | No biogeochemistry |
| **SWAT** | Watershed hydrology | Not designed for tidal systems |
| **CE-QUAL-W2** | 2D WQ modeling | Complex setup, no network topology |

### The C-GEM Solution

C-GEM Network fills a critical niche for **rapid, process-based biogeochemical modeling** of complex tidal networks:

```
                         MODEL COMPLEXITY
                    ┌─────────────────────────────────────┐
                    │                                     │
    Empirical ◄─────┼──────────────────────────────────►  │ Process-Based
    (LOADEST)       │                                     │ (Full RIVE)
                    │              ★ C-GEM               │
                    │              (Sweet Spot)           │
                    │                   \                 │
                    │   HEC-RAS ●        \   Delft3D ●   │
                    │                     \              │
                    │         MIKE11 ●     \  TELEMAC ● │
                    │                       ↓            │
                    │                    COMPUTATIONAL   │
                    │                       COST         │
                    └─────────────────────────────────────┘
```

### What Makes C-GEM Unique?

| Feature | C-GEM | Delft3D/MIKE | SWAT | HEC-RAS |
|---------|-------|--------------|------|---------|
| **Setup time** | Hours | Days-Weeks | Days | Hours |
| **Run time (30 days)** | Minutes | Hours-Days | Hours | Minutes |
| **Network topology** | ✅ Native | ⚠️ Complex meshing | ❌ | ⚠️ Limited |
| **Tidal dispersion** | ✅ Savenije theory | ⚠️ Numerical | ❌ | ❌ |
| **Full biogeochemistry** | ✅ C-RIVE (30 species) | ⚠️ Limited WQ | ⚠️ | ❌ |
| **GHG emissions** | ✅ CO₂, CH₄, N₂O | ❌ Often missing | ❌ | ❌ |
| **Auto-calibration** | ✅ NLopt built-in | ⚠️ External tools | ✅ | ❌ |
| **Land-use coupling** | ✅ Rainfall-driven | ⚠️ Manual | ✅ | ❌ |
| **Open source** | ✅ MIT License | ❌ Commercial | ✅ | ✅ |
| **Cost** | Free | €1,000-50,000/yr | Free | Free |

---

## C-GEM's Unique Innovation: Rainfall-Driven Lateral Loads

Traditional models require users to **guess** lateral pollution factors. C-GEM automates this using **globally available data**:

```
Traditional Approach:
┌──────────────────┐    ┌───────────────────┐    ┌────────────────┐
│ Manual Factor    │ →  │ "Q_Factor = 10"   │ →  │ "Why 10?"      │
│ (Expert Guess)   │    │ (Arbitrary)       │    │ (No Validation)│
└──────────────────┘    └───────────────────┘    └────────────────┘

C-GEM Rainfall-Driven Approach:
┌──────────────────┐    ┌───────────────────┐    ┌────────────────┐
│ WorldClim/TRMM   │ →  │ Rain_Sep = 340mm  │ →  │ Q = Rain/Base  │
│ (Global Dataset) │    │ (Verifiable)      │    │ (Physics)      │
└──────────────────┘    └───────────────────┘    └────────────────┘
```

### Data-to-Model Pipeline

```
┌─────────────────┐     ┌──────────────────┐     ┌─────────────────┐
│  JAXA/Sentinel  │     │   WorldClim/     │     │   WorldPop/     │
│  Land Use Map   │     │   TRMM Rainfall  │     │   Census Data   │
└────────┬────────┘     └────────┬─────────┘     └────────┬────────┘
         │                       │                        │
         ▼                       ▼                        ▼
┌────────────────────────────────────────────────────────────────────┐
│              generate_lateral_loads_v2.py                          │
│  ┌────────────────┐  ┌────────────────┐  ┌────────────────┐       │
│  │ Emission       │  │ Runoff Physics │  │ Point Sources  │       │
│  │ Factors (EMC)  │  │ Q = f(Rain)    │  │ from Population│       │
│  └────────────────┘  └────────────────┘  └────────────────┘       │
└────────────────────────────────────────────────────────────────────┘
         │                       │                        │
         ▼                       ▼                        ▼
┌────────────────┐    ┌────────────────┐     ┌────────────────┐
│ lateral_       │    │ lateral_       │     │ point_         │
│ sources.csv    │    │ seasonal_      │     │ sources.csv    │
│ (Base loads)   │    │ factors.csv    │     │ (Cities)       │
└────────────────┘    └────────────────┘     └────────────────┘
         │                       │                        │
         └───────────────────────┼────────────────────────┘
                                 ▼
                    ┌────────────────────────┐
                    │   C-GEM Simulation     │
                    │   Q_lat = Q_base ×     │
                    │         Q_factor(day)  │
                    └────────────────────────┘
```

---

## Target Applications

C-GEM is designed for:

1. **Research questions requiring many model runs**
   - Sensitivity analysis (100s of runs)
   - Uncertainty quantification
   - Climate scenarios (RCP 4.5, 8.5)
   - Monte Carlo simulations

2. **Data-sparse tropical regions**
   - Deltas with limited monitoring
   - Developing countries without expensive model licenses
   - Regions with only satellite-derived data

3. **Integrated assessment studies**
   - Coupling with land-use change models
   - Policy scenario evaluation
   - Water-food-energy nexus analysis
   - Rapid screening before detailed 2D/3D modeling

4. **GHG emission budgets**
   - CO₂ flux from estuarine respiration
   - CH₄ ebullition and oxidation
   - N₂O from nitrification/denitrification

## Key Capabilities

### Multi-Branch Network Topology

Unlike 1D models that simulate single channels, C-GEM handles complex networks:

```
          Upstream 1                    Upstream 2
              │                             │
              ▼                             ▼
         ┌────────┐                    ┌────────┐
         │ Tien   │                    │  Hau   │
         │ Main   │                    │  Main  │
         └────┬───┘                    └────┬───┘
              │          ┌──────┐          │
              └──────────┤VamNao├──────────┘
                         └──────┘
                              │
              ┌───────────────┼───────────────┐
              ▼               ▼               ▼
         ┌────────┐      ┌────────┐      ┌────────┐
         │Co Chien│      │Ham     │      │ Hau    │
         └────┬───┘      │Luong   │      │ River  │
              │          └────┬───┘      └────┬───┘
              ▼               ▼               ▼
           Ocean           Ocean           Ocean
```

### Complete Biogeochemistry

The C-RIVE module includes 30 species and 35+ reactions:

| Process | Species | Reference |
|---------|---------|-----------|
| Primary production | PHY1, PHY2, O₂, DIC | Garnier et al. (1995) |
| Organic matter degradation | TOC → DIC, O₂ consumption | Billen et al. (1994) |
| Nitrification | NH₄ → NO₂ → NO₃ | Wang et al. (2018) |
| Denitrification | NO₃ → N₂ | Hasanyar et al. (2022) |
| Carbonate chemistry | DIC, TA, pH, pCO₂ | Millero (2007) |
| GHG emissions | CO₂, CH₄, N₂O | Process-based |
| Sediment dynamics | SPM, erosion, deposition | Winterwerp & van Kesteren (2004) |
| Flocculation | Salinity-dependent settling | Wolanski (2007) |

### Automatic Calibration

Built-in calibration with NLopt optimization:

```powershell
# Run 3-stage calibration
.\bin\Debug\CGEM_Network.exe INPUT/Cases/Mekong_Delta_Full/case_config.txt --calibrate

# Stage 1: Hydrodynamics (tidal range, salinity intrusion)
# Stage 2: Sediment (SPM, ETM location)
# Stage 3: Biogeochemistry (O2, nutrients)
```

## Getting Started

Ready to dive in? Start with:

1. [Installation Guide](installation.md) - Set up your environment
2. [Quick Start](quickstart.md) - Run your first simulation
3. [Project Structure](structure.md) - Understand the codebase
