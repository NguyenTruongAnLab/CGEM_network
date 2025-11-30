# Data Requirements and Preparation

## Overview

This document describes the input data required to run C-GEM Network and provides guidance on data sources, formats, and preparation steps. Proper input data is critical for realistic simulations.

---

## Input Directory Structure

```
INPUT/Cases/YourCase/
├── case_config.txt        # Main configuration file
├── topology.csv           # Network structure
├── boundary_map.csv       # Boundary condition assignments
├── biogeo_params.txt      # Biogeochemistry parameters (optional)
└── forcing/               # Time-varying boundary data
    ├── tidal_levels.csv   # Water level at ocean boundaries
    ├── river_discharge.csv # Q at river boundaries
    └── species_*.csv      # Concentrations at boundaries
```

---

## 1. Network Configuration

### case_config.txt

Main configuration file controlling the simulation:

```ini
# Case identification
CaseName = Tien_River

# Input file paths (relative to case directory)
Topology = topology.csv
Boundaries = boundary_map.csv

# Simulation period
StartDate = 2017-01-01
Duration = 130              # days
Warmup_Days = 100           # days before recording output

# Numerical parameters
TimeStepSeconds = 300       # dt [s]
DELTI = 300                 # hydrodynamic timestep [s]
DELXI = 2000                # grid spacing [m]

# Physical parameters
TidalAmplitude = 1.9        # [m] - simplified tidal forcing
RiverDischarge = 2000.0     # [m³/s] - baseline river flow

# Output settings
WriteCSV = 1                # 0=off, 1=on
WriteNetCDF = 0             # 0=off, 1=on
WriteReactionRates = 0      # 0=off, 1=on
```

### topology.csv

Defines branches and their geometry:

```csv
branch_id,name,node_up,node_down,length_m,width_down_m,width_up_m,depth_m,chezy,lc_convergence
1,Tien_Main,1,2,100000,1600,1000,15,50,212764
2,Tien_Connect,2,3,16000,1100,1200,14,50,-183884
3,My_Tho,3,4,76000,1500,800,11,50,120902
4,Ham_Luong,3,5,84000,1800,850,12,50,111954
5,Co_Chien,3,6,100000,2200,900,14,50,111880
```

**Column definitions:**

| Column | Type | Description |
|--------|------|-------------|
| branch_id | int | Unique branch identifier |
| name | string | Branch name (max 64 chars) |
| node_up | int | Upstream node ID |
| node_down | int | Downstream node ID |
| length_m | float | Branch length [m] |
| width_down_m | float | Width at downstream end [m] |
| width_up_m | float | Width at upstream end [m] |
| depth_m | float | Reference depth [m] |
| chezy | float | Chézy friction coefficient [m^0.5/s] |
| lc_convergence | float | Convergence length [m], negative = diverging |

### boundary_map.csv

Assigns boundary conditions to nodes:

```csv
node_id,type,forcing_file,species_file
1,DISCHARGE,river_discharge.csv,species_river.csv
4,LEVEL,tidal_MyTho.csv,species_ocean.csv
5,LEVEL,tidal_HamLuong.csv,species_ocean.csv
6,LEVEL,tidal_CoChien.csv,species_ocean.csv
```

**Node types:**
- `JUNCTION` (or `0`): Internal connection, no forcing
- `DISCHARGE` (or `1`): Prescribed Q(t)
- `LEVEL` (or `2`): Prescribed H(t)

---

## 2. Hydrodynamic Forcing

### Tidal Levels

Format: Time-series of water level at ocean boundaries.

```csv
time_hours,level_m
0.0,0.50
1.0,0.92
2.0,1.25
3.0,1.42
...
```

**Data sources:**
- Tide gauges (Can Tho, My Thuan, Vung Tau)
- Harmonic analysis (TPXO, FES2014)
- Regional models (MIKE, TELEMAC outputs)

**Harmonic constituents for Mekong:**

| Constituent | Period [h] | Amplitude [m] | Phase [°] |
|-------------|------------|---------------|-----------|
| M2 | 12.42 | 0.9-1.2 | 150-180 |
| S2 | 12.00 | 0.3-0.5 | 180-210 |
| K1 | 23.93 | 0.4-0.6 | 30-60 |
| O1 | 25.82 | 0.3-0.4 | 0-30 |

### River Discharge

Format: Time-series of discharge at river boundaries.

```csv
time_hours,Q_m3s
0.0,2500.0
24.0,2480.0
48.0,2450.0
...
```

**Data sources:**
- MRC stations (Tan Chau, Chau Doc, Can Tho)
- National gauging stations (SIWRP, DONRE)
- Upstream dam releases (cascade operations)

**Typical seasonal discharge (Tien River at Tan Chau):**

| Month | Q [m³/s] | Description |
|-------|----------|-------------|
| Jan-Mar | 2,000-3,000 | Dry season base flow |
| Apr-May | 2,500-4,000 | Late dry, early transition |
| Jun-Jul | 5,000-10,000 | Monsoon onset |
| Aug-Sep | 15,000-25,000 | Peak flood season |
| Oct-Nov | 10,000-20,000 | Flood recession |
| Dec | 5,000-8,000 | Post-monsoon |

---

## 3. Water Quality Boundary Conditions

### Species Forcing Files

Format: Time-series of concentrations at each boundary.

**River boundary (species_river.csv):**

```csv
time_hours,salinity,phy1,phy2,dsi,no3,nh4,po4,o2,toc,spm,dic,at
0.0,0.1,5.0,3.0,120.0,35.0,12.0,2.5,280.0,300.0,50.0,1800.0,1850.0
24.0,0.1,5.2,3.1,118.0,34.0,11.5,2.4,275.0,310.0,55.0,1820.0,1870.0
...
```

**Ocean boundary (species_ocean.csv):**

```csv
time_hours,salinity,phy1,phy2,dsi,no3,nh4,po4,o2,toc,spm,dic,at
0.0,30.0,2.0,1.0,15.0,5.0,2.0,0.5,200.0,50.0,10.0,2100.0,2350.0
...
```

### Required Species

| Species | River [typical] | Ocean [typical] | Unit |
|---------|-----------------|-----------------|------|
| salinity | 0.1-1.0 | 25-35 | PSU |
| phy1 | 5-50 | 1-5 | µg C/L |
| phy2 | 3-30 | 1-3 | µg C/L |
| dsi | 100-200 | 5-20 | µmol/L |
| no3 | 20-80 | 1-10 | µmol N/L |
| nh4 | 5-30 | 1-5 | µmol N/L |
| po4 | 1-5 | 0.2-1.0 | µmol P/L |
| o2 | 200-300 | 180-220 | µmol/L |
| toc | 200-500 | 30-80 | µmol C/L |
| spm | 30-200 | 5-20 | mg/L |
| dic | 1500-2500 | 2000-2200 | µmol/L |
| at | 1600-2600 | 2200-2400 | µeq/L |

### Optional GHG Species

| Species | River [typical] | Ocean [typical] | Unit |
|---------|-----------------|-----------------|------|
| no2 | 0.5-2.0 | 0.1-0.5 | µmol N/L |
| n2o | 10-50 | 5-15 | nmol N/L |
| ch4 | 50-500 | 2-10 | nmol C/L |

---

## 4. Biogeochemistry Parameters

### biogeo_params.txt (optional)

Override default parameters:

```ini
# Environmental
water_temp = 28.0
I0 = 400.0

# Phytoplankton
pbmax1 = 2.0
pbmax2 = 1.5
kmort1 = 0.05
kmort2 = 0.04

# Nutrients
kdsi1 = 5.0
kn1 = 2.0
kpo41 = 0.5

# Decomposition
kox = 0.1
kdenit = 0.2
knit = 0.1

# Carbonate
pco2_atm = 420.0

# GHG (C-RIVE)
N2O_yield_nit = 0.004
N2O_yield_denit = 0.01
CH4_prod_rate = 0.01
```

---

## 5. Data Sources for the Mekong Delta

### Hydrology

| Source | Data | Access |
|--------|------|--------|
| MRC Data Portal | Discharge, water level | portal.mrcmekong.org |
| SIWRP Vietnam | Hydrological stations | Request required |
| DONRE Provinces | Local gauging data | Government channels |
| TPXO/FES | Tidal constituents | Online databases |

### Water Quality

| Source | Data | Access |
|--------|------|--------|
| MRC Routine Monitoring | N, P, suspended solids | portal.mrcmekong.org |
| DONRE Monitoring | DO, nutrients | Government reports |
| Research Publications | Full biogeochemistry | Scientific literature |

### Key Publications for Mekong Biogeochemistry

1. **Nguyen et al. (2019)** - Nutrient dynamics in the Mekong
2. **Toming et al. (2020)** - DOC in tropical rivers
3. **Borges et al. (2018)** - CO2 in Southeast Asian rivers
4. **Yoshimura et al. (2015)** - pCO2 in the Mekong

---

## 6. Data Preparation Steps

### Step 1: Define Network Topology

1. Identify branches (main channels, distributaries)
2. Determine branch lengths from maps/GIS
3. Measure widths at key sections
4. Estimate depths from navigation charts
5. Calculate convergence lengths using Savenije's method:

$$L_c = -\frac{x_2 - x_1}{\ln(B_2/B_1)}$$

### Step 2: Gather Hydrodynamic Forcing

1. Obtain tidal harmonics or measured levels
2. Get river discharge records
3. Align timestamps, fill gaps
4. Convert to required format (hours from start)

### Step 3: Compile Water Quality Data

1. Collect monitoring data at boundaries
2. Calculate seasonal averages if time-series unavailable
3. Estimate missing species from literature
4. Ensure consistency (mass balance at junctions)

### Step 4: Validate Input Data

Run the `generate_mekong_rive_inputs.py` script to:
- Generate synthetic data based on literature ranges
- Visualize seasonal patterns
- Check for outliers or inconsistencies

```bash
python scripts/generate_mekong_rive_inputs.py --case Tien_River
```

---

## 7. Unit Conventions

### Concentration Units

| Species Group | Unit | Conversion |
|---------------|------|------------|
| Salinity | PSU | 1 PSU ≈ 1 g/kg |
| Phytoplankton | µg C/L | Chl-a × 40-50 |
| Nutrients | µmol/L | mg/L ÷ MW × 1000 |
| Oxygen | µmol/L | mg/L × 31.25 |
| Carbon | µmol C/L | mg C/L × 83.3 |
| SPM | mg/L | - |

### Flux Units

| Flux Type | Unit |
|-----------|------|
| River discharge | m³/s |
| Nutrient loading | kg/day |
| Gas flux | mmol/m²/day |

### Rate Units

| Rate Type | Storage | Calculation |
|-----------|---------|-------------|
| Growth rates | /day | Converted to /s internally |
| Decay rates | /day | Converted to /s internally |
| Gas transfer | m/s | Direct |

---

## 8. Quality Control Checklist

Before running a simulation:

- [ ] **Topology**: All branches connected, no orphan nodes
- [ ] **Boundaries**: Every terminal node has boundary condition
- [ ] **Time alignment**: All forcing files use same time origin
- [ ] **Units**: All concentrations in model units
- [ ] **Physical ranges**: No negative concentrations, realistic values
- [ ] **Mass balance**: River inputs + ocean inputs ≈ expected loads
- [ ] **Time step**: CFL condition satisfied (dt < dx/c)
- [ ] **Duration**: Long enough for spinup + analysis period

---

## 9. Common Issues and Solutions

### Issue: Salt intrusion too far / not far enough

**Solutions:**
- Adjust river discharge
- Check tidal amplitude
- Modify dispersion coefficient (D0, Van den Burgh β)
- Verify channel geometry

### Issue: Unrealistic oxygen levels

**Solutions:**
- Check O2 saturation at boundary temperature
- Verify reaeration rate (wind speed, k600)
- Adjust organic matter decomposition rates
- Include benthic oxygen demand

### Issue: pCO2 too high / too low

**Solutions:**
- Check DIC and TA boundary conditions
- Verify respiration rates
- Adjust benthic respiration flux
- Check gas exchange coefficient

### Issue: Numerical instability

**Solutions:**
- Reduce time step
- Check for extreme gradients in forcing
- Verify positive concentrations everywhere
- Enable RK4 solver for stiff reactions

---

## 10. Example: Tien River Setup

Complete example configuration for the Tien River case:

### Network (5 branches, 6 nodes)

```
Node 1: River discharge (Tan Chau)
Node 2: Junction (Tien-Connection)
Node 3: Junction (Bifurcation)
Node 4: Ocean level (My Tho mouth)
Node 5: Ocean level (Ham Luong mouth)
Node 6: Ocean level (Co Chien mouth)
```

### Forcing Data Summary

| Boundary | Type | Typical Values |
|----------|------|----------------|
| Node 1 | Q | 2000-5000 m³/s (dry), 15000-25000 (wet) |
| Node 4 | H | ±1.5 m (mixed semidiurnal) |
| Node 5 | H | ±1.5 m (phase lag ~30 min from Node 4) |
| Node 6 | H | ±1.5 m (phase lag ~60 min from Node 4) |

### Typical Run Time

| Configuration | Duration | Wall Time |
|---------------|----------|-----------|
| 30-day spinup | 30 days | ~30 seconds |
| Full year | 365 days | ~5 minutes |
| With full biogeochemistry | 365 days | ~10 minutes |

---

## References

1. MRC (2005). Overview of the Hydrology of the Mekong Basin.
2. Nguyen, A.D. (2008). Salt Intrusion in Multi-Channel Estuaries. PhD Thesis, Delft.
3. SIWRP (2015). Mekong Delta Water Resources Assessment.
4. Savenije, H.H.G. (2005). Salinity and Tides in Alluvial Estuaries. Elsevier.

---

## Next Steps

- [Quick Start (Windows)](QUICKSTART_WINDOWS.md) - Run your first simulation
- [Introduction](introduction.md) - Why C-GEM?
- [Biogeochemistry](biogeochemistry.md) - C-RIVE module details
