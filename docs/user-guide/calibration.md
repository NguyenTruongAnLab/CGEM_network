# Calibration

## Overview

C-GEM includes a built-in calibration module powered by **NLopt**, supporting:

- **Multi-stage calibration**: Hydro → Sediment → Water Quality → GHG (4 stages)
- **Multiple algorithms**: BOBYQA, Nelder-Mead, DIRECT, COBYLA
- **Seasonal objectives**: Time-series RMSE for dry/wet season
- **Config-driven**: No recompilation needed

## Quick Start

```powershell
# Full 4-stage calibration
.\bin\Debug\CGEM_Network.exe INPUT/Cases/Mekong_Delta_Full/case_config.txt --calibrate

# Single stage with options
.\bin\Debug\CGEM_Network.exe INPUT/Cases/Mekong_Delta_Full/case_config.txt --calibrate --stage 1 --max-iter 10 --verbose 2
```

## Command-Line Options

| Option | Description | Default |
|--------|-------------|---------|
| `--calibrate` | Enable calibration mode | Off |
| `--stage N` | Run only stage N (1, 2, 3, or 4) | All stages |
| `--max-iter N` | Maximum optimizer iterations | 100 |
| `--verbose N` | Verbosity level (0-2) | 1 |

## Four-Stage Calibration Workflow

### Stage 1: Hydrodynamics + Salinity

**Goal**: Match tidal propagation and salinity intrusion (60-80 km for Mekong dry season)

**Parameters**:

| Name | Description | Typical Range | Unit |
|------|-------------|---------------|------|
| `CHEZY` | Friction coefficient | 40-80 | m^0.5/s |
| `LC_CONV` | Convergence length | 10-100 | km |
| `VDB_COEF` | Van den Burgh coefficient | 0.1-0.5 | - |
| `MIXING_ALPHA` | Fischer mixing efficiency | 0.2-0.8 | - |
| `D0` | Dispersion at mouth | 500-2000 | m²/s |
| `RS` | Storage width ratio | 1-5 | - |

**Objectives**:

- Tidal range at gauging stations
- Salinity intrusion length (4 PSU isohaline)
- Salinity profile along estuary

### Stage 2: Sediment Transport (SPM)

**Goal**: Reproduce SPM patterns and ETM location

**Parameters**:

| Name | Description | Typical Range | Unit |
|------|-------------|---------------|------|
| `WS` | Settling velocity | 0.0001-0.005 | m/s |
| `TAU_ERO` | Critical erosion stress | 0.05-0.5 | Pa |
| `TAU_DEP` | Critical deposition stress | 0.01-0.2 | Pa |
| `FLOC_SAL_SCALE` | Flocculation salinity scale | 1-15 | PSU |
| `FLOC_FACTOR_MAX` | Maximum flocculation factor | 2-20 | - |

**Objectives**:

- Mean SPM at monitoring stations (8-38 mg/L)
- ETM location (maximum turbidity position)

### Stage 3: Water Quality (O₂, Nutrients, TOC, Chl-a)

**Goal**: Match water quality observations

**Parameters**:

| Name | Description | Typical Range | Unit |
|------|-------------|---------------|------|
| `KOX` | TOC degradation rate | 0.005-0.10 | d⁻¹ |
| `KNIT` | Nitrification rate | 0.02-0.20 | d⁻¹ |
| `KDENIT` | Denitrification rate | 0.01-0.15 | d⁻¹ |
| `PBMAX1` | Diatom max production | 1.0-5.0 | d⁻¹ |
| `PBMAX2` | Non-siliceous max production | 0.8-4.0 | d⁻¹ |
| `KMORT1` | Diatom mortality rate | 0.03-0.20 | d⁻¹ |
| `KMORT2` | Non-siliceous mortality | 0.02-0.15 | d⁻¹ |
| `BENTHIC_RESP` | Benthic respiration | 5-100 | µmol/m²/day |

**Objectives**:

- O₂ concentration (170-260 µmol/L)
- NO₃ concentration (7-66 µmol/L)
- NH₄ concentration (0.7-5 µmol/L)
- TOC concentration (107-220 µmol/L)
- Chlorophyll-a (0.7-7 µg/L)

### Stage 4: Greenhouse Gases (pCO₂, pH, CH₄, N₂O)

**Goal**: Match GHG concentrations and emissions

**Parameters**:

| Name | Description | Typical Range | Unit |
|------|-------------|---------------|------|
| `N2O_YIELD_NIT` | N₂O yield from nitrification | 0.001-0.02 | mol/mol |
| `N2O_YIELD_DENIT` | N₂O yield from denitrification | 0.005-0.05 | mol/mol |
| `BENTHIC_CH4_FLUX` | Benthic CH₄ flux | 50-1000 | nmol/m²/day |
| `BENTHIC_N2O_FLUX` | Benthic N₂O flux | 5-100 | nmol/m²/day |
| `WIND_SPEED` | Wind speed for gas exchange | 2-8 | m/s |
| `WIND_COEFF` | Wind gas exchange coefficient | 0.15-0.50 | - |
| `PCO2_ATM` | Atmospheric pCO₂ | 400-450 | µatm |

**Objectives**:

- pCO₂ (500-4700 µatm, gradient from ocean to upstream)
- pH (7.3-8.2, gradient from ocean to upstream)
- CH₄ (38-243 nmol/L, high upstream due to rice paddies)
- N₂O (6-52 nmol/L)

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
