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
| **TELEMAC-2D/3D** | Unstructured grids, morphodynamics | Complex setup, expensive computation |
| **Delft3D-FLOW/WAQ** | Comprehensive WQ | Steep learning curve, commercial |
| **MIKE 11/21** | User-friendly GUI | Expensive licensing |
| **HEC-RAS** | Free, flood modeling | No biogeochemistry |

### The C-GEM Solution

C-GEM Network fills the niche for **rapid, process-based biogeochemical modeling**:

```
                    ┌─────────────────────────────────────┐
                    │         MODEL COMPLEXITY            │
                    │                                     │
    Empirical ◄─────┼──────────────────────────────────►  │ Process-Based
    (LOADEST)       │                                     │ (RIVE)
                    │              C-GEM ●                │
                    │                   \                 │
                    │   HEC-RAS ●        \   Delft3D ●   │
                    │                     \              │
                    │         MIKE11 ●     \  TELEMAC ● │
                    │                       ↓            │
                    │                    COMPUTATIONAL   │
                    │                       COST         │
                    └─────────────────────────────────────┘
```

**C-GEM combines:**

- ✅ Process-based biogeochemistry (C-RIVE)
- ✅ Network topology (graph-based)
- ✅ Fast computation (1D)
- ✅ Open source (ANSI C)

## Target Applications

C-GEM is designed for:

1. **Research questions requiring many model runs**
   - Sensitivity analysis
   - Uncertainty quantification
   - Climate scenarios

2. **Data-sparse regions**
   - Tropical deltas with limited monitoring
   - Developing countries without expensive model licenses

3. **Integrated assessment**
   - Coupling with land-use models
   - Policy scenario evaluation
   - Rapid screening studies

## Key Capabilities

### Multi-Branch Network Topology

Unlike 1D models that simulate single channels, C-GEM handles:

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

The C-RIVE module includes:

| Process | Species | Reference |
|---------|---------|-----------|
| Primary production | PHY1, PHY2, O₂, DIC | Garnier et al. (1995) |
| Organic matter degradation | TOC → DIC, O₂ consumption | Billen et al. (1994) |
| Nitrification | NH₄ → NO₂ → NO₃ | Wang et al. (2018) |
| Denitrification | NO₃ → N₂ | Hasanyar et al. (2022) |
| Carbonate chemistry | DIC, TA, pH, pCO₂ | Millero (2007) |
| GHG emissions | CO₂, CH₄, N₂O | Process-based |

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
