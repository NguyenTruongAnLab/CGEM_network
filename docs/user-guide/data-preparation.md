# Data Preparation

!!! tip "Quick Start"
    **Not sure what data you need?** See [Data Requirements](data-requirements.md) for a tiered approach:
    
    - **Tier 1 (Minimal)**: Geometry + discharge + tides = salinity modeling
    - **Tier 2 (Standard)**: + land use + rainfall = full biogeochemistry
    - **Tier 3 (Advanced)**: + calibration data = publication quality

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
├── landuse_map.csv           # Land use for lateral load generation
├── calibration_params.csv    # Calibration parameters (optional)
├── calibration_targets.csv   # Calibration objectives (optional)
├── seasonal_targets.csv      # Seasonal observations (optional)
├── forcing_data/             # Time-varying boundary conditions
│   ├── Tien_Inlet.csv        # Upstream discharge time series
│   ├── Hau_Inlet.csv         # Upstream discharge time series
│   ├── MyTho_Tide.csv        # Tidal level time series
│   └── ...
├── lateral_loads_monthly/    # Monthly lateral load files (if using --annual)
└── polder_loads/             # Virtual Polder gate factors (if using --polders)
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

# Time parameters (FULL YEAR RECOMMENDED)
StartDate = 2024-01-01
Duration = 365              # days - full year captures seasonal variation
Warmup = 30                 # days before recording output

# Numerical parameters
TimeStep = 60.0             # dt [s] - 60s for stability
DELXI = 2000.0              # grid spacing [m]

# Biogeochemistry mode (NEW)
SimplifiedMode = 1          # 1 = 80/20 parsimonious mode (recommended)
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
| `SimplifiedMode` | int | 1 = 80/20 biogeochemistry mode (recommended for tropical systems) |

### Recommended Simulation Durations

| Purpose | Duration | Warmup | Notes |
|---------|----------|--------|-------|
| Quick test | 10 days | 3 days | Verify setup works |
| Tidal analysis | 30 days | 10 days | Captures spring-neap cycle |
| Seasonal study | 365 days | 30 days | **RECOMMENDED** - full seasonal variation |
| Multi-year | 730+ days | 30 days | Climate variability analysis |

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

## 6. Lateral Sources System (NEW - Rainfall-Driven)

C-GEM features a **smart lateral load system** that automates pollution input generation from:

- **Diffuse sources**: Agricultural runoff, urban stormwater (from land use)
- **Point sources**: City sewage outfalls (population-based)
- **Seasonal factors**: Rainfall-driven Q and concentration multipliers

### Why Rainfall-Driven?

Traditional models require users to **guess** seasonal load factors (e.g., "Q_Factor = 10 in wet season"). C-GEM's approach is different:

| Aspect | Traditional | C-GEM Rainfall-Driven |
|--------|-------------|----------------------|
| **Input** | Manual Q_Factor = 10 | Rainfall = 340 mm |
| **Validation** | "Why 10?" | "WorldClim shows 340mm in Sep" |
| **Physics** | None | Q = Rain/Base × Runoff_C |
| **Transferability** | Site-specific guessing | Climate presets for any delta |

### Quick Start Workflow

```powershell
# Step 1: Generate land use map (if no GIS data available)
python scripts/generate_synthetic_landuse.py

# Step 2: Generate all lateral load files using Mekong climate
python scripts/generate_lateral_loads_v2.py --climate Mekong

# Alternative: Use custom rainfall (12 monthly values in mm)
python scripts/generate_lateral_loads_v2.py --rainfall 15,8,20,55,180,260,290,310,340,270,130,45
```

This generates **4 files** automatically:

| File | Purpose |
|------|---------|
| `lateral_sources.csv` | Base loads (dry season reference) |
| `lateral_seasonal_factors.csv` | Monthly Q and concentration multipliers |
| `lateral_daily_factors.csv` | Daily interpolated factors (365 days) |
| `point_sources.csv` | City sewage with treatment efficiency |

---

### 6.1 Land Use Map (`landuse_map.csv`)

The land use map is the **foundation** for lateral load calculation. It maps river segments to surrounding land use types.

#### Format

```csv
Branch,Distance_km,Segment_Area_km2,Pct_Urban,Pct_Rice,Pct_Aqua,Pct_Mangrove,Pct_Fruit,Pct_Forest
Tien_Main,0,2.0,5,60,10,5,15,5
Tien_Main,2,2.0,5,55,15,5,15,5
Ham_Luong,0,2.0,10,30,40,15,5,0
...
```

#### Data Sources for Land Use

| Source | Resolution | Coverage | Cost | Link |
|--------|------------|----------|------|------|
| **JAXA ALOS-2** | 25m | Global | Free | [JAXA EORC](https://www.eorc.jaxa.jp/ALOS/en/dataset/lulc_e.htm) |
| **ESA WorldCover** | 10m | Global | Free | [WorldCover](https://esa-worldcover.org/) |
| **Sentinel-2** | 10m | Global | Free | [Copernicus](https://scihub.copernicus.eu/) |
| **GlobeLand30** | 30m | Global | Free | [GlobeLand30](http://www.globallandcover.com/) |
| **National GIS** | varies | Country | varies | Contact local agencies |

#### Creating Your Own Land Use Map

**Option A: GIS Analysis (Recommended)**

1. Download satellite imagery (JAXA/Sentinel/Landsat)
2. Perform supervised classification in QGIS/ArcGIS
3. Buffer river centerline (e.g., 2km each side)
4. Calculate zonal statistics for each segment
5. Export to CSV format

**Option B: Use Provided Script (Demo/Testing)**

```powershell
python scripts/generate_synthetic_landuse.py
```

This creates a synthetic land use map based on typical Mekong Delta patterns.

#### JAXA Land Use Classes → C-GEM Mapping

| JAXA Class Code | JAXA Description | C-GEM Category |
|-----------------|------------------|----------------|
| 1 | Water bodies | (excluded) |
| 2 | Urban and built-up | Urban |
| 3 | Paddy field | Rice |
| 4 | Cropland | Fruit |
| 5 | Grassland | Forest |
| 6 | Forest | Forest |
| 7 | Mangrove | Mangrove |
| 8 | Wetland | Mangrove |
| 9 | Aquaculture pond | Aqua |

---

### 6.2 Base Loads (`lateral_sources.csv`)

Contains the **dry season base loads** that will be multiplied by seasonal factors.

#### Format (NEW v2)

```csv
Branch,Segment_Index,Distance_km,Area_km2,Runoff_C,Is_Polder_Zone,Q_lat_base_m3_s,NH4_conc_base_mg_L,NO3_conc_base_mg_L,PO4_conc_base_mg_L,TOC_conc_base_mg_L,DIC_conc_base_mg_L,SPM_conc_base_mg_L
Tien_Main,0,0,2.0,0.435,True,0.00174,3.3,5.9,0.89,23.5,15.2,85.0
Ham_Luong,15,30,2.0,0.70,False,0.0028,8.0,3.0,2.0,80.0,30.0,200.0
```

#### Column Descriptions

| Column | Type | Unit | Description |
|--------|------|------|-------------|
| `Branch` | string | - | Branch name (must match topology.csv) |
| `Segment_Index` | int | - | Grid cell index (0 = upstream end) |
| `Distance_km` | float | km | Distance from upstream |
| `Area_km2` | float | km² | Contributing drainage area |
| `Runoff_C` | float | - | Runoff coefficient (0.1-0.9) |
| `Is_Polder_Zone` | bool | - | True if rice/aqua dominated (>50%) |
| `Q_lat_base_m3_s` | float | m³/s | Base lateral inflow rate |
| `NH4_conc_base_mg_L` | float | mg N/L | Base NH4 concentration |
| `NO3_conc_base_mg_L` | float | mg N/L | Base NO3 concentration |
| `PO4_conc_base_mg_L` | float | mg P/L | Base PO4 concentration |
| `TOC_conc_base_mg_L` | float | mg C/L | Base TOC concentration |
| `DIC_conc_base_mg_L` | float | mg C/L | Base DIC concentration |
| `SPM_conc_base_mg_L` | float | mg/L | Base SPM concentration |

#### JAXA Emission Factors (Built-in)

The Python script uses these literature-based **Event Mean Concentrations (EMC)**:

| Land Use | NH4 (mg/L) | NO3 (mg/L) | PO4 (mg/L) | TOC (mg/L) | SPM (mg/L) | Runoff_C | Source |
|----------|------------|------------|------------|------------|------------|----------|--------|
| **Urban** | 10.0 | 3.0 | 1.5 | 50 | 150 | 0.85 | Burton & Pitt (2002) |
| **Rice** | 2.5 | 6.0 | 0.8 | 20 | 80 | 0.40 | Yan et al. (2003) |
| **Aqua** | 8.0 | 3.0 | 2.0 | 80 | 200 | 0.70 | Páez-Osuna (2001) |
| **Mangrove** | 0.2 | 0.1 | 0.1 | 100 | 50 | 0.90 | Alongi (2014) |
| **Fruit** | 3.0 | 8.0 | 1.0 | 25 | 60 | 0.30 | Regional data |
| **Forest** | 0.3 | 0.5 | 0.05 | 15 | 20 | 0.15 | Natural background |

---

### 6.3 Seasonal Factors (`lateral_seasonal_factors.csv`)

Contains **monthly multipliers** for Q and concentrations, derived from rainfall physics.

#### Format

```csv
Month,Month_Name,Season,Rain_mm,Q_Factor,NH4_Factor,NO3_Factor,PO4_Factor,TOC_Factor,SPM_Factor
1,Jan,dry,15,1.0,1.0,1.0,1.0,1.0,1.0
2,Feb,dry,8,1.0,1.0,1.0,1.0,1.0,1.0
...
9,Sep,wet,340,13.878,0.895,0.716,1.074,1.645,1.645
...
```

#### Physics: Rainfall → Factors

**Q Factor** (Flow multiplier):
$$Q_{factor} = \frac{Rain_{month}}{Rain_{dry\_average}}$$

**Concentration Factor** (First-flush + dilution):

- **Rising limb** (Q_factor < 5): First flush wash-off
  $$C_{factor} = Q_{factor}^{0.25}$$

- **Peak flow** (Q_factor ≥ 5): Dilution dominates
  $$C_{factor} = (5^{0.25}) \times (5/Q_{factor})^{0.4}$$

#### Climate Presets Available

| Preset | Region | Annual Rain | Dry Months | Peak Month |
|--------|--------|-------------|------------|------------|
| `Mekong` | Vietnam | 1923 mm | Jan-Apr | Sep (340mm) |
| `RedRiver` | N. Vietnam | 1668 mm | Jan-Mar, Dec | Aug (320mm) |
| `SaigonDongNai` | S. Vietnam | 1892 mm | Jan-Apr | Sep (320mm) |
| `Ganges` | Bangladesh | 1720 mm | Jan-Mar, Dec | Jul (350mm) |
| `Niger` | Nigeria | 2325 mm | Jan-Feb, Dec | Jul (380mm) |
| `Irrawaddy` | Myanmar | 2395 mm | Jan-Apr, Dec | Jul (550mm) |
| `Mediterranean` | Europe | 605 mm | Jun-Aug | Dec (95mm) |

#### Data Sources for Rainfall

| Source | Resolution | Period | Cost | Link |
|--------|------------|--------|------|------|
| **WorldClim** | 1km | 1970-2000 | Free | [WorldClim](https://www.worldclim.org/) |
| **TRMM/GPM** | 0.25° | 1998-present | Free | [NASA GES DISC](https://disc.gsfc.nasa.gov/) |
| **ERA5** | 0.25° | 1979-present | Free | [Copernicus CDS](https://cds.climate.copernicus.eu/) |
| **CHIRPS** | 0.05° | 1981-present | Free | [CHC UCSB](https://www.chc.ucsb.edu/data/chirps) |
| **Local gauges** | Point | varies | varies | National met agencies |

#### Using Custom Rainfall

```powershell
# Custom 12-value array (mm): Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec
python scripts/generate_lateral_loads_v2.py --rainfall 10,5,15,40,150,250,300,320,350,280,100,30
```

---

### 6.4 Daily Factors (`lateral_daily_factors.csv`)

Optional **daily interpolated factors** for smooth seasonal transitions.

#### Format

```csv
Day,Day_of_Year,Month,Season,Q_Factor,NH4_Factor,NO3_Factor,PO4_Factor,TOC_Factor,SPM_Factor
0,0,1,dry,1.0,1.0,1.0,1.0,1.0,1.0
1,1,1,dry,1.0001,1.0001,1.0001,1.0001,1.0001,1.0001
...
243,243,9,wet,13.5,0.91,0.73,1.09,1.64,1.64
...
```

The C model automatically uses daily factors if the file exists, otherwise falls back to monthly.

---

### 6.5 Point Sources (`point_sources.csv`)

Contains **city sewage outfalls** with population-based loads.

#### Format

```csv
Name,Branch,Segment_Index,Distance_km,Population,Treatment,Q_m3_s,NH4_mg_L,NO3_mg_L,PO4_mg_L,TOC_mg_L,DIC_mg_L,SPM_mg_L
Can_Tho,Hau_River,40,80.0,1500000,primary,2.6042,72000.0,1.0,6000.0,280000.0,30.0,50.0
My_Tho,My_Tho,20,40.0,500000,none,0.8681,80000.0,1.0,6452.0,400000.0,30.0,50.0
```

#### Per-Capita Emission Factors

| Parameter | Value | Unit | Source |
|-----------|-------|------|--------|
| Water use | 150 | L/person/day | WHO guidelines |
| NH4 | 12 | g N/person/day | Metcalf & Eddy |
| PO4 | 2 | g P/person/day | Metcalf & Eddy |
| TOC | 60 | g C/person/day | BOD equivalent |

#### Treatment Removal Efficiencies

| Treatment Level | NH4 | NO3 | PO4 | TOC | SPM |
|-----------------|-----|-----|-----|-----|-----|
| none | 0% | 0% | 0% | 0% | 0% |
| primary | 10% | 0% | 10% | 30% | 50% |
| secondary | 60% | 20% | 30% | 80% | 90% |
| tertiary | 90% | 80% | 90% | 95% | 95% |

#### Population Data Sources

| Source | Coverage | Cost | Link |
|--------|----------|------|------|
| **WorldPop** | Global | Free | [WorldPop](https://www.worldpop.org/) |
| **GPWv4** | Global | Free | [SEDAC](https://sedac.ciesin.columbia.edu/data/collection/gpw-v4) |
| **National census** | Country | varies | Local statistical offices |
| **OSM/Wikipedia** | Global | Free | City population estimates |

---

### 6.6 Complete Workflow Example

```powershell
# Navigate to project directory
cd C:\Users\nguytruo\Documents\Github\CGEM_network

# 1. Generate land use map (if needed)
python scripts/generate_synthetic_landuse.py

# 2. Generate all lateral load files with Mekong climate
python scripts/generate_lateral_loads_v2.py --climate Mekong

# 3. Verify generated files
Get-ChildItem INPUT/Cases/Mekong_Delta_Full/*.csv | Select-Object Name

# 4. Run simulation
.\bin\Debug\CGEM_Network.exe INPUT\Cases\Mekong_Delta_Full\case_config.txt
```

#### Expected Output

```
======================================================================
C-GEM LATERAL LOAD GENERATOR (Rainfall-Driven v2)
======================================================================

Using climate preset: Mekong
  Description: Mekong Delta, Vietnam (Monsoon tropical)
  Rainfall (mm): [15, 8, 20, 55, 180, 260, 290, 310, 340, 270, 130, 45]
  Annual total: 1923 mm
  Dry months: [1, 2, 3, 4]

  Monthly Factors:
Month_Name  Rain_mm  Q_Factor  NH4_Factor  TOC_Factor
       Jan       15     1.000       1.000       1.300
       Feb        8     1.000       1.000       1.300
       ...
       Sep      340    13.878       0.895       1.645
       ...

Loaded 466 segments from 9 branches
Generated 4 files successfully.
```

---

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
