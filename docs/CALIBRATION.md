# C-GEM Calibration Module

## Overview

The calibration module provides a **generic, data-driven framework** for parameter optimization in C-GEM Network. It supports multi-stage calibration using config files, eliminating the need to recompile when changing calibration targets.

## Quick Start

```powershell
# Run full 3-stage calibration
.\bin\Debug\CGEM_Network.exe INPUT/Cases/Mekong_Delta_Full/case_config.txt --calibrate

# Run only Stage 1 (hydrodynamics)
.\bin\Debug\CGEM_Network.exe INPUT/Cases/Mekong_Delta_Full/case_config.txt --calibrate --stage 1

# With options
.\bin\Debug\CGEM_Network.exe INPUT/Cases/Mekong_Delta_Full/case_config.txt --calibrate --max-iter 50 --verbose 2
```

## Three-Stage Calibration Workflow

### Stage 1: Hydrodynamics (Tidal Range & Salinity Intrusion)

**Parameters:**
| Name | Description | Typical Range | Unit |
|------|-------------|---------------|------|
| `CHEZY` | Friction coefficient | 40-80 | m^0.5/s |
| `LC_CONV` | Convergence length | 10-100 | km |
| `STORAGE_RATIO` | Storage width ratio (RS) | 1-15 | - |
| `VDB_COEF` | Van den Burgh coefficient | 0.1-0.8 | - |
| `D0` | Dispersion at mouth | 100-2000 | m²/s |
| `MANNING_N` | Manning friction (alternative) | 0.02-0.15 | s/m^1/3 |

**Objectives:**
- Tidal range at gauging stations (2-3 m for Mekong)
- Salt intrusion length (4 PSU isohaline position)
- Mean velocity in transfer channels (Vam Nao)

### Stage 2: Sediment Transport (SPM & ETM)

**Parameters:**
| Name | Description | Typical Range | Unit |
|------|-------------|---------------|------|
| `WS` | Settling velocity | 0.0001-0.005 | m/s |
| `FLOC_SAL_SCALE` | Flocculation salinity scale | 1-5 | PSU |
| `FLOC_FACTOR_MAX` | Max flocculation factor | 5-20 | - |
| `TAU_ERO` | Erosion threshold | 0.1-1.0 | Pa |
| `TAU_DEP` | Deposition threshold | 0.05-0.3 | Pa |
| `MERO` | Erosion coefficient | 1e-7 - 1e-5 | kg/m²/s |

**Objectives:**
- Estuarine Turbidity Maximum (ETM) location
- Maximum SPM concentration at ETM
- Upstream SPM concentrations

### Stage 3: Biogeochemistry (Water Quality)

**Parameters:**
| Name | Description | Typical Range | Unit |
|------|-------------|---------------|------|
| `KOX` | TOC oxidation rate | 0.02-0.5 | 1/day |
| `KNIT` | Nitrification rate | 0.05-0.3 | 1/day |
| `KDENIT` | Denitrification rate | 0.01-0.15 | 1/day |
| `PBMAX1/2` | Max photosynthesis | 1-5 | 1/day |
| `KMORT1/2` | Phytoplankton mortality | 0.02-0.2 | 1/day |
| `WIND_COEFF` | Gas exchange wind factor | 0.1-0.5 | - |
| `BENTHIC_RESP` | Benthic respiration | 30-150 | µmol/m²/day |

**Objectives:**
- Dissolved oxygen (>4 mg/L in main stems)
- Nutrient gradients (NO3, NH4, PO4)
- TOC concentrations

## Configuration Files

### calibration_params.csv

Defines parameters to optimize:

```csv
Name,TargetType,TargetID,VarType,Min,Max,Initial
Chezy_Estuary,GROUP,2,CHEZY,40.0,70.0,55.0
VDB_River,GROUP,1,VDB_COEF,0.1,0.5,0.3
Kox_Global,GLOBAL,0,KOX,0.05,0.5,0.15
LC_MyTho,BRANCH,My_Tho,LC_CONV,15000,60000,30000
```

**Columns:**
- `Name`: User-defined parameter identifier
- `TargetType`: `GLOBAL`, `GROUP`, or `BRANCH`
- `TargetID`: Group ID (1=River, 2=Estuary) or Branch name
- `VarType`: Parameter type (see tables above)
- `Min, Max`: Optimization bounds
- `Initial`: Starting value

### calibration_targets.csv

Defines observations to match:

```csv
BranchName,Variable,Statistic,TargetValue,Weight,Threshold
Ham_Luong,SALINITY,INTRUSION_LEN,55.0,50.0,4.0
Vam_Nao,VELOCITY,MEAN,0.35,100.0,0.0
Hau_River,O2,MIN,125.0,50.0,0.0
My_Tho,SPM,ETM_LOCATION,25.0,40.0,0.0
```

**Columns:**
- `BranchName`: Branch to evaluate
- `Variable`: `SALINITY`, `O2`, `SPM`, `VELOCITY`, `WATERLEVEL`, etc.
- `Statistic`: 
  - `INTRUSION_LEN`: Distance where variable crosses threshold
  - `MEAN`, `MIN`, `MAX`: Statistical measures
  - `ETM_LOCATION`: Location of SPM maximum
  - `TIDAL_RANGE`: Max-min water level
- `TargetValue`: Observed/expected value
- `Weight`: Relative importance (higher = more weight)
- `Threshold`: For intrusion length (e.g., 4 PSU)

## Tropical Estuary Adaptations

C-RIVE was developed for temperate rivers. Key adaptations for tropical systems:

### 1. Higher Reaction Rates
- Tropical temperatures (25-30°C) increase all reaction rates ~2x
- Use higher `KOX`, `KNIT`, `KDENIT` initial values

### 2. Salinity Stress on Phytoplankton
- Freshwater species stressed in brackish zones
- Consider adding `SAL_STRESS_THRESH`, `SAL_STRESS_COEF` parameters

### 3. Strong Flocculation
- Large sediment loads + salt wedge = intense flocculation
- `FLOC_FACTOR_MAX` may need 10-15 (vs 5-8 temperate)

### 4. Monsoon Variability
- Dry season: Low flow, high salinity intrusion, high residence time
- Wet season: High flow, limited intrusion, rapid flushing
- Calibrate separately for each season if needed

## Optimizer Details

The calibration uses a **Nelder-Mead simplex algorithm** (no external dependencies required). For each evaluation:

1. Update network parameters from current guess
2. Run full simulation (warmup + actual)
3. Calculate weighted RMSE objective
4. Simplex decision (reflect, expand, contract, shrink)

**Options:**
- `--max-iter N`: Maximum iterations (default 100)
- `--verbose N`: 0=silent, 1=summary, 2=per-eval, 3=debug

## Output Files

After calibration, results are saved to:

```
OUTPUT/<CaseName>/calibration/
├── calibration_results.csv    # Final parameter values
└── (log files if verbose)
```

## Adding New Parameters

1. Add to `CalibVarType` enum in `src/optimization/calibration.h`
2. Add name to `VAR_TYPE_NAMES` array in `calibration.c`
3. Add case in `set_branch_param()` or `set_global_param()` functions
4. Update this documentation

## Example: Mekong Delta Dry Season

```csv
# calibration_params.csv - Focus on convergence length for intrusion
LC_HauRiver,BRANCH,Hau_River,LC_CONV,20000,80000,40000
LC_MyTho,BRANCH,My_Tho,LC_CONV,15000,60000,30000
Chezy_Estuary,GROUP,2,CHEZY,40,70,55
VDB_Estuary,GROUP,2,VDB_COEF,0.1,0.6,0.3

# calibration_targets.csv - Match observed intrusion
Hau_River,SALINITY,INTRUSION_LEN,55.0,50,4.0
My_Tho,SALINITY,INTRUSION_LEN,50.0,50,4.0
Vam_Nao,VELOCITY,MEAN,0.35,100,0
```

## References

- Savenije (2005, 2012) - Salt intrusion in alluvial estuaries
- Wang et al. (2018) - C-RIVE sensitivity analysis
- Nguyen et al. (2008) - Mekong salinity intrusion
- Winterwerp (2002) - Flocculation dynamics
