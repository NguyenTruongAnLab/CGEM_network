# Data Preparation

## Input Directory Structure

Each simulation case requires a directory under `INPUT/Cases/`:

```
INPUT/Cases/YourCase/
├── case_config.txt           # Main configuration
├── topology.csv              # Network structure
├── boundary_map.csv          # Boundary conditions
├── biogeo_params.txt         # Biogeochemistry parameters (optional)
├── species_river.csv         # Upstream species concentrations
├── species_ocean.csv         # Ocean species concentrations
├── lateral_sources.csv       # Point/diffuse loads (optional)
├── calibration_params.csv    # Calibration parameters (optional)
├── calibration_targets.csv   # Calibration objectives (optional)
└── seasonal_targets.csv      # Seasonal observations (optional)
```

## 1. Main Configuration

### case_config.txt

```ini
# Case identification
CaseName = Mekong_Delta_Full

# Input file paths (relative to case directory)
Topology = topology.csv
BoundaryMap = boundary_map.csv
BiogeoParams = biogeo_params.txt

# Output
OutputDir = OUTPUT/Mekong_Delta_Full
WriteCSV = 1
WriteNetCDF = 0
WriteReactionRates = 0

# Time parameters
StartDate = 2024-03-01
Duration = 30              # days
Warmup = 10                # days before recording output

# Numerical parameters
TimeStep = 300.0           # dt [s]
DELXI = 2000.0             # grid spacing [m]
```

### Parameter Descriptions

| Parameter | Type | Description |
|-----------|------|-------------|
| `CaseName` | string | Case identifier |
| `Topology` | path | Branch definitions file |
| `BoundaryMap` | path | Boundary conditions file |
| `BiogeoParams` | path | Biogeochemistry parameters (optional) |
| `OutputDir` | path | Output directory |
| `StartDate` | date | Simulation start (YYYY-MM-DD) |
| `Duration` | int | Simulation duration [days] |
| `Warmup` | int | Warmup period [days] |
| `TimeStep` | float | Time step [seconds] |
| `DELXI` | float | Target grid spacing [meters] |

## 2. Network Topology

### topology.csv

Defines branches in the network:

```csv
BranchID,BranchName,NodeUp,NodeDown,Length_m,Width_Down_m,Width_Up_m,Depth_m,Chezy,RS
1,Tien_Main,1,3,48000,1200,1000,15,55,1.0
2,Hau_Main,2,3,44000,1100,900,14,55,1.0
3,Vam_Nao,3,4,8000,600,600,18,65,1.0
```

### Column Descriptions

| Column | Type | Unit | Description |
|--------|------|------|-------------|
| `BranchID` | int | - | Unique branch identifier |
| `BranchName` | string | - | Human-readable name |
| `NodeUp` | int | - | Upstream node ID |
| `NodeDown` | int | - | Downstream node ID |
| `Length_m` | float | m | Branch length |
| `Width_Down_m` | float | m | Width at downstream end |
| `Width_Up_m` | float | m | Width at upstream end |
| `Depth_m` | float | m | Reference depth |
| `Chezy` | float | m^0.5/s | Friction coefficient |
| `RS` | float | - | Storage width ratio (RS ≥ 1) |

### Optional Columns

| Column | Type | Unit | Description |
|--------|------|------|-------------|
| `LC_m` | float | m | Convergence length (if not computed) |
| `VDB` | float | - | Van den Burgh coefficient |
| `Manning` | float | s/m^1/3 | Manning's n (alternative to Chezy) |

## 3. Boundary Conditions

### boundary_map.csv

Assigns boundary conditions to nodes:

```csv
NodeID,BoundaryType,ForcingFile,SpeciesFile
1,DISCHARGE,upstream_Q_tien.csv,species_river.csv
2,DISCHARGE,upstream_Q_hau.csv,species_river.csv
6,LEVEL,tide_cochien.csv,species_ocean.csv
7,LEVEL,tide_mytho.csv,species_ocean.csv
```

### Boundary Types

| Type | Description | Required Forcing |
|------|-------------|-----------------|
| `DISCHARGE` | River inflow | Q [m³/s] time series |
| `LEVEL` | Tidal boundary | H [m] time series |
| `JUNCTION` | Internal node | (automatic) |

### Forcing File Format

```csv
time_hours,value
0.0,2500.0
1.0,2510.0
2.0,2495.0
...
```

For tidal boundaries:
```csv
time_hours,value
0.0,0.5
0.5,1.2
1.0,1.8
1.5,1.5
2.0,0.8
...
```

## 4. Species Boundary Conditions

### species_river.csv (upstream)

```csv
time_days,salinity,phy1,phy2,dsi,no3,nh4,po4,o2,toc,spm,dic,at
0,0.1,50,30,120,30,5,1.5,280,800,50,1800,1900
14,0.1,80,50,100,35,8,2.0,260,900,80,1850,1950
...
```

### species_ocean.csv (downstream)

```csv
time_days,salinity,phy1,phy2,dsi,no3,nh4,po4,o2,toc,spm,dic,at
0,35,10,5,5,2,0.5,0.3,220,200,10,2100,2300
14,35,15,8,6,3,0.8,0.4,210,220,15,2100,2300
...
```

### Species Units

| Species | Unit | Typical Range |
|---------|------|---------------|
| salinity | PSU | 0-35 |
| phy1, phy2 | µg C/L | 0-500 |
| dsi | µmol/L | 0-200 |
| no3, nh4 | µmol N/L | 0-100 |
| po4 | µmol P/L | 0-10 |
| o2 | µmol/L | 0-400 |
| toc | µmol C/L | 100-2000 |
| spm | mg/L | 0-500 |
| dic | µmol/L | 1000-3000 |
| at | µmol/L | 1500-3000 |

## 5. Biogeochemistry Parameters

### biogeo_params.txt

```ini
# Phytoplankton
pbmax1 = 2.5          # Max photosynthesis rate PHY1 [d-1]
pbmax2 = 2.0          # Max photosynthesis rate PHY2 [d-1]
kmort1 = 0.08         # Mortality rate PHY1 [d-1]
kmort2 = 0.06         # Mortality rate PHY2 [d-1]

# Nutrients
kox = 0.08            # Aerobic degradation rate [d-1]
knit = 0.10           # Nitrification rate [d-1]
kdenit = 0.03         # Denitrification rate [d-1]

# Oxygen
wind_k = 0.25         # Wind enhancement coefficient
current_k = 0.30      # Current enhancement coefficient

# Sediment
ws_base = 0.0005      # Base settling velocity [m/s]
tau_ero = 0.30        # Critical erosion stress [Pa]
tau_dep = 0.08        # Critical deposition stress [Pa]
floc_sal_scale = 2.0  # Flocculation salinity scale [PSU]
floc_factor_max = 8.0 # Max flocculation enhancement [-]

# Salinity stress (freshwater phytoplankton)
sal_stress_thresh = 5.0   # Onset of stress [PSU]
sal_stress_coef = 0.5     # Mortality enhancement [-]
```

## 6. Lateral Sources (Optional)

### lateral_sources.csv

Point and diffuse pollution loads:

```csv
BranchName,Location_km,Q_m3_s,NH4_umol,NO3_umol,PO4_umol,TOC_umol
Ham_Luong,25.0,0.5,100,20,5,500
Co_Chien,30.0,0.3,80,15,3,400
...
```

## 7. Calibration Files

### calibration_params.csv

```csv
Name,TargetType,TargetID,VarType,Min,Max,Initial,Stage
Chezy_River,BRANCH,Tien_Main,CHEZY,45,75,60,1
Chezy_Estuary,BRANCH,Ham_Luong,CHEZY,40,70,55,1
VDB_Estuary,BRANCH,Co_Chien,VDB_COEF,0.1,0.6,0.3,1
WS_Global,GLOBAL,-,WS,0.0001,0.003,0.0005,2
Kox_Global,GLOBAL,-,KOX,0.02,0.3,0.08,3
```

### calibration_targets.csv

```csv
BranchName,Variable,Statistic,TargetValue,Weight,Threshold
Ham_Luong,SALINITY,INTRUSION_LEN,55.0,50.0,4.0
Co_Chien,SALINITY,MEAN,15.0,30.0,
Hau_River,O2,MEAN,180.0,40.0,
```

### seasonal_targets.csv

```csv
Branch,Variable,Location_km,Time_Day,Value,Weight
Ham_Luong,SALINITY,30.0,15.0,25.0,2.0
Ham_Luong,SALINITY,30.0,180.0,0.5,1.0
Co_Chien,SALINITY,25.0,15.0,28.0,2.0
```

## Data Sources

### Typical Data Sources for Tropical Deltas

| Data Type | Sources |
|-----------|---------|
| Discharge | MRC, national hydrology agencies |
| Tidal levels | NOAA, local tide gauges |
| Water quality | EPA, monitoring programs |
| Satellite | Sentinel-2/3 (Chl-a, SPM) |
| Bathymetry | Nautical charts, surveys |

### Minimum Data Requirements

| Priority | Data | Impact |
|----------|------|--------|
| **Essential** | Discharge, tidal levels | Drives entire system |
| **Important** | Salinity, SPM | Calibration targets |
| **Useful** | O₂, nutrients | Biogeochemistry validation |
| **Optional** | CO₂, CH₄ | GHG budget |
