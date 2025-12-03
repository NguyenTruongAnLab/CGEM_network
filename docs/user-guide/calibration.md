# Calibration

## Overview

C-GEM includes a built-in calibration module powered by **NLopt**, supporting:

- **Multi-stage calibration**: Hydro → Sediment → Biogeochemistry
- **Multiple algorithms**: BOBYQA, Nelder-Mead, DIRECT, COBYLA
- **Seasonal objectives**: Time-series RMSE for dry/wet season
- **Config-driven**: No recompilation needed

## Quick Start

```powershell
# Full 3-stage calibration
.\bin\Debug\CGEM_Network.exe INPUT/Cases/Mekong_Delta_Full/case_config.txt --calibrate

# Single stage with options
.\bin\Debug\CGEM_Network.exe INPUT/Cases/Mekong_Delta_Full/case_config.txt --calibrate --stage 1 --max-iter 100 --verbose 2
```

## Command-Line Options

| Option | Description | Default |
|--------|-------------|---------|
| `--calibrate` | Enable calibration mode | Off |
| `--stage N` | Run only stage N (1, 2, or 3) | All stages |
| `--max-iter N` | Maximum optimizer iterations | 100 |
| `--verbose N` | Verbosity level (0-2) | 1 |

## Three-Stage Workflow

### Stage 1: Hydrodynamics

**Goal**: Match tidal propagation and salinity intrusion

**Parameters**:

| Name | Description | Typical Range | Unit |
|------|-------------|---------------|------|
| `CHEZY` | Friction coefficient | 40-80 | m^0.5/s |
| `LC_CONV` | Convergence length | 10-100 | km |
| `VDB_COEF` | Van den Burgh coefficient | 0.1-0.8 | - |
| `D0` | Dispersion at mouth | 100-2000 | m²/s |
| `RS` | Storage width ratio | 1-15 | - |

**Objectives**:

- Tidal range at gauging stations
- Salinity intrusion length (4 PSU isohaline)
- Mean velocity at key sections

### Stage 2: Sediment Transport

**Goal**: Reproduce SPM patterns and ETM location

**Parameters**:

| Name | Description | Typical Range | Unit |
|------|-------------|---------------|------|
| `WS` | Settling velocity | 0.1-3.0 | mm/s |
| `TAU_ERO` | Critical erosion stress | 0.1-0.8 | Pa |
| `TAU_DEP` | Critical deposition stress | 0.03-0.2 | Pa |
| `FLOC_FACTOR` | Flocculation enhancement | 1-15 | - |

**Objectives**:

- Mean SPM at monitoring stations
- ETM location (maximum turbidity position)

### Stage 3: Biogeochemistry

**Goal**: Match water quality observations

**Parameters**:

| Name | Description | Typical Range | Unit |
|------|-------------|---------------|------|
| `KOX` | Aerobic degradation rate | 0.02-0.3 | d⁻¹ |
| `KNIT` | Nitrification rate | 0.03-0.25 | d⁻¹ |
| `KDENIT` | Denitrification rate | 0.01-0.1 | d⁻¹ |
| `PBMAX` | Max photosynthesis rate | 1-4 | d⁻¹ |
| `KMORT` | Mortality rate | 0.02-0.2 | d⁻¹ |

**Objectives**:

- O₂ concentration at stations
- Nutrient concentrations (NO₃, NH₄, PO₄)
- Chlorophyll-a (via PHY1 + PHY2)

## Configuration Files

### calibration_params.csv

Defines parameters to optimize:

```csv
Name,TargetType,TargetID,VarType,Min,Max,Initial,Stage
Chezy_River,BRANCH,Tien_Main,CHEZY,45,75,60,1
Chezy_Estuary,BRANCH,Ham_Luong,CHEZY,40,70,55,1
VDB_Estuary,BRANCH,Co_Chien,VDB_COEF,0.1,0.6,0.3,1
WS_Global,GLOBAL,-,WS,0.0001,0.003,0.0005,2
Kox_Global,GLOBAL,-,KOX,0.02,0.3,0.08,3
```

**Columns**:

| Column | Description |
|--------|-------------|
| `Name` | Parameter identifier |
| `TargetType` | `BRANCH` or `GLOBAL` |
| `TargetID` | Branch name (if BRANCH) or `-` |
| `VarType` | Variable type (CHEZY, VDB_COEF, WS, etc.) |
| `Min`, `Max` | Parameter bounds |
| `Initial` | Starting value |
| `Stage` | Calibration stage (1, 2, or 3) |

### calibration_targets.csv

Defines calibration objectives:

```csv
BranchName,Variable,Statistic,TargetValue,Weight,Threshold
Ham_Luong,SALINITY,INTRUSION_LEN,55.0,50.0,4.0
Co_Chien,SALINITY,MEAN,15.0,30.0,
Hau_River,O2,MEAN,180.0,40.0,
```

**Columns**:

| Column | Description |
|--------|-------------|
| `BranchName` | Target branch |
| `Variable` | Species name or WATERLEVEL/VELOCITY |
| `Statistic` | MEAN, MIN, MAX, INTRUSION_LEN, ETM_LOCATION |
| `TargetValue` | Observed value |
| `Weight` | Importance (higher = more weight) |
| `Threshold` | For intrusion: salinity threshold (PSU) |

### seasonal_targets.csv (NEW)

Defines time-series observations for seasonal calibration:

```csv
Branch,Variable,Location_km,Time_Day,Value,Weight
Ham_Luong,SALINITY,30.0,15.0,25.0,2.0
Ham_Luong,SALINITY,30.0,180.0,0.5,1.0
Co_Chien,SALINITY,25.0,15.0,28.0,2.0
Co_Chien,SALINITY,25.0,180.0,0.8,1.0
```

**Columns**:

| Column | Description |
|--------|-------------|
| `Branch` | Target branch name |
| `Variable` | Species name |
| `Location_km` | Distance from mouth [km] |
| `Time_Day` | Day from simulation start |
| `Value` | Observed value |
| `Weight` | Importance |

This allows calibration against both **dry season** (Day 15) and **wet season** (Day 180) observations simultaneously.

## NLopt Algorithms

C-GEM supports multiple optimization algorithms:

| Algorithm ID | Name | Description | Best For |
|--------------|------|-------------|----------|
| 0 | `NELDERMEAD` | Derivative-free simplex | General use |
| 1 | `BOBYQA` | Bounded optimization (default) | **Recommended** |
| 2 | `DIRECT` | Global optimization | Exploring parameter space |
| 3 | `COBYLA` | Constrained optimization | With constraints |
| 4 | `CRS2_LM` | Controlled random search | Global search |
| 5 | `SBPLX` | Subplex | Noisy objectives |

Set algorithm in code or use the default (BOBYQA).

## Output

### Console Output

```
==============================================
  CALIBRATION MODE
==============================================

Loaded 30 calibration parameters from calibration_params.csv
Loaded 32 calibration objectives from calibration_targets.csv
Loaded 14 seasonal observation groups from seasonal_targets.csv

Optimizing 11 parameters for stage 1 using NLopt
NLopt algorithm: BOBYQA bound-constrained optimization
Max evaluations: 100
Tolerance: rel=1.00e-04, abs=1.00e-06

  Eval 1: RMSE = 45.2341
  Eval 2: RMSE = 42.1234
  ...
  Eval 50: RMSE = 12.4523

NLopt result: FTOL_REACHED (code 4)
Optimization complete: best RMSE = 12.4523 after 50 evaluations

========== CALIBRATION SUMMARY ==========

OPTIMIZED PARAMETERS:
Name                      Initial    Best       Min        Max
--------------------------------------------------------------
Chezy_River                  60.0000    58.5234    45.0000    75.0000
VDB_Estuary                   0.3000     0.4521     0.1000     0.6000
...

Final RMSE: 12.4523
Total evaluations: 50
```

### Results File

Results are saved to `OUTPUT/<CaseName>/calibration/calibration_results.csv`:

```csv
# C-GEM Calibration Results
# Final RMSE: 12.4523
# Total evaluations: 50

# Optimized Parameters
Name,VarType,Initial,Best,Min,Max,Stage
Chezy_River,CHEZY,60.000000,58.523400,45.000000,75.000000,1
VDB_Estuary,VDB_COEF,0.300000,0.452100,0.100000,0.600000,1
...
```

## Best Practices

### 1. Start with Hydrodynamics

Always calibrate Stage 1 first. Biogeochemistry results depend on correct hydrodynamics.

### 2. Use Reasonable Bounds

Tight bounds speed convergence but may miss optimal values:

```csv
# Too tight - may miss optimum
Chezy_River,BRANCH,Tien_Main,CHEZY,55,65,60,1

# Too loose - slow convergence
Chezy_River,BRANCH,Tien_Main,CHEZY,20,100,60,1

# Good - informed by literature
Chezy_River,BRANCH,Tien_Main,CHEZY,45,75,60,1
```

### 3. Weight Objectives Appropriately

Primary objectives should have higher weights:

```csv
# Primary: salinity intrusion (weight 50)
Ham_Luong,SALINITY,INTRUSION_LEN,55.0,50.0,4.0

# Secondary: mean concentration (weight 20)
Ham_Luong,O2,MEAN,180.0,20.0,
```

### 4. Include Seasonal Data

For tropical systems with strong wet/dry contrast, use seasonal_targets.csv:

```csv
# Dry season peak (Day 15)
Ham_Luong,SALINITY,30.0,15.0,25.0,2.0

# Wet season minimum (Day 180)  
Ham_Luong,SALINITY,30.0,180.0,0.5,1.0
```

### 5. Validate Results

After calibration:

1. Run simulation with best parameters
2. Compare time series visually
3. Check spatial patterns
4. Test with independent data (if available)
