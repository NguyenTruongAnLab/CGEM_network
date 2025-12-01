# Mekong Delta Case Study

## Overview

The Mekong Delta Full case (`Mekong_Delta_Full`) demonstrates C-GEM's capabilities for a complete multi-branch estuarine network with full seasonal simulation.

### Network Structure

- **9 branches**: Main rivers and distributaries
- **10 nodes**: 2 upstream (discharge), 4 ocean (tidal), 4 junctions
- **Total length**: ~450 km of channels
- **Grid**: 2 km resolution, ~225 cells total
- **Simulation**: 365 days with 30-day warmup (full seasonal cycle)

```
                 Tien River (Q)           Hau River (Q)
                      │                        │
                      ▼                        ▼
                 ┌─────────┐              ┌─────────┐
                 │Tien Main│              │Hau Main │
                 │ 50 km   │              │ 45 km   │
                 └────┬────┘              └────┬────┘
                      │                        │
                      └──────────┬─────────────┘
                                 │
                            ┌────┴────┐
                            │ Vam Nao │ ← PRISMATIC (600m × 600m)
                            │  7 km   │
                            └────┬────┘
                                 │
              ┌──────────────────┼─────────────────┐
              │                  │                 │
         ┌────┴─────┐       ┌────┴────┐       ┌────┴────┐
         │Tien Lower│       │         │       │         │
         │  50 km   │       │         │       │         │
         └────┬─────┘       │         │       │         │
              │             │         │       │         │
    ┌─────────┼─────────────┤         │       │         │
    │         │             │         │       │         │
┌───┴───┐ ┌───┴───┐    ┌────┴────┐    │  ┌────┴────┐    │
│Co     │ │Tien   │    │         │    │  │Hau River│    │
│Chien  │ │Conn   │    │         │    │  │  90 km  │    │
│ 85 km │ │15 km  │    │         │    │  └────┬────┘    │
└───┬───┘ └───┬───┘    │         │    │       │         │
    │         │        │         │    │       │         │
    │    ┌────┴────┐   │         │    │       │         │
    │    │MyTho    │   │         │    │       │         │
    │    │Split    │   │         │    │       │         │
    │    └─┬───┬───┘   │         │    │       │         │
    │      │   │       │         │    │       │         │
    │  ┌───┴┐ ┌┴───┐   │         │    │       │         │
    │  │My  │ │Ham │   │         │    │       │         │
    │  │Tho │ │Lng │   │         │    │       │         │
    │  │55km│ │60km│   │         │    │       │         │
    │  └──┬─┘ └─┬──┘   │         │    │       │         │
    ▼     ▼     ▼      ▼         ▼    ▼       ▼         ▼
  Ocean Ocean Ocean  Ocean     Ocean Ocean  Ocean     Ocean
   (8)   (6)   (7)    ...       ...  ...    (9)       ...
```

**Hierarchical Topology (Key Feature):**
1. Co Chien branches FIRST at Node 5 (geographically correct)
2. Tien_Connector (15 km) links to Node 10
3. My Tho and Ham Luong branch SECOND at Node 10

## Configuration

### case_config.txt

```ini
CaseName = Mekong_Delta_Full

# Files
Topology = topology.csv
BoundaryMap = boundary_map.csv
BiogeoParams = biogeo_params.txt
OutputDir = OUTPUT/Mekong_Delta_Full

# Time - FULL YEAR SIMULATION
StartDate = 2024-01-01
Duration = 365
Warmup = 30

# Numerics
TimeStep = 60.0
DELXI = 2000.0

# Output
WriteCSV = 1
WriteNetCDF = 1

# Biogeochemistry mode
SimplifiedMode = 1    # 80/20 parsimonious mode (recommended)
```

### topology.csv

```csv
BranchID,BranchName,NodeUp,NodeDown,Length_m,Width_Down_m,Width_Up_m,Depth_m,Chezy,RS
1,Tien_Main,1,3,50000,1200,1000,15,65,1.0
2,Hau_Main,2,4,45000,1100,900,14,65,1.0
3,Vam_Nao,3,4,7000,600,600,18,75,1.0
4,Tien_Lower,3,5,50000,1400,1200,14,65,1.0
5,Tien_Connector,5,10,15000,1200,1100,13,65,1.0
6,Co_Chien,5,8,85000,2800,1100,13,60,1.2
7,My_Tho,10,6,55000,2200,900,11,60,1.2
8,Ham_Luong,10,7,60000,2500,1000,12,60,1.2
9,Hau_River,4,9,90000,3000,1200,14,58,1.3
```

### Boundary Conditions

**Upstream (Discharge)**:
- Node 1: Tien River inflow (60% of Mekong, seasonal variation)
- Node 2: Hau River inflow (40% of Mekong, seasonal variation)

**Downstream (Tidal Level)**:
- Node 6: My Tho mouth (M2 + S2 tides, phase = 0°)
- Node 7: Ham Luong mouth (phase = +8°)
- Node 8: Co Chien mouth (phase = +15°)
- Node 9: Hau River mouth (phase = +30°)

### Seasonal Discharge Variation

Based on MRC monitoring data at Tan Chau:

| Month | Total Q (m³/s) | Tien (60%) | Hau (40%) |
|-------|----------------|------------|-----------|
| Jan | 3,250 | 1,950 | 1,300 |
| Feb | 2,600 | 1,560 | 1,040 |
| Mar | 2,340 | 1,404 | 936 |
| Apr | 2,860 | 1,716 | 1,144 |
| May | 5,200 | 3,120 | 2,080 |
| Jun | 10,400 | 6,240 | 4,160 |
| Jul | 15,600 | 9,360 | 6,240 |
| **Aug** | **20,800** | **12,480** | **8,320** |
| **Sep** | **23,400** | **14,040** | **9,360** |
| Oct | 19,500 | 11,700 | 7,800 |
| Nov | 10,400 | 6,240 | 4,160 |
| Dec | 5,200 | 3,120 | 2,080 |

## Running the Case

### Full Workflow (Recommended)

```powershell
# Step 1: Generate case configuration files
python scripts/setup_mekong_simplified.py

# Step 2: Generate land use map
python scripts/generate_synthetic_landuse.py

# Step 3: Generate lateral loads with seasonal variation
python scripts/generate_lateral_loads.py --annual --polders

# Step 4: Build and run
.\scripts\build-and-run.ps1 -r Mekong_Delta_Full
```

### Standard Simulation

```powershell
.\bin\Debug\CGEM_Network.exe INPUT/Cases/Mekong_Delta_Full/case_config.txt
```

Output:
```
==============================================
  C-GEM Network Engine
==============================================
  Case        : Mekong_Delta_Full
  Branches    : 9
  Nodes       : 10
  Species     : 30
  Time step   : 60.0 s
  Duration    : 365 days
  Warmup      : 30 days
  Mode        : Simplified (80/20)
==============================================

Day    1.0 (  0.3%) [WARMUP] | H=15.40 m, U=0.142 m/s, S=0.1 | 2.4 s elapsed
Day   30.0 (  8.2%)          | H=15.39 m, U=0.173 m/s, S=15.4 | 45.2 s elapsed
...
Day  180.0 ( 49.3%)          | H=15.80 m, U=0.380 m/s, S=2.1 | 5.5 min elapsed
Day  270.0 ( 74.0%)          | H=16.20 m, U=0.450 m/s, S=0.5 | 8.2 min elapsed
Day  365.0 (100.0%)          | H=15.50 m, U=0.210 m/s, S=12.3 | 11.1 min elapsed

Simulation complete: 525600 steps in 666.2 seconds (789.1 steps/s)
```

### With Calibration

```powershell
.\bin\Debug\CGEM_Network.exe INPUT/Cases/Mekong_Delta_Full/case_config.txt --calibrate --stage 1 --max-iter 50
```

## Performance Benchmarks

### Computational Performance

| Configuration | Duration | Time Steps | Wall Clock | Steps/s |
|---------------|----------|------------|------------|---------|
| 30-day test | 30 days | 43,200 | 52 sec | 830 |
| **365-day full** | **365 days** | **525,600** | **11 min** | **790** |
| 730-day (2 yr) | 730 days | 1,051,200 | 22 min | 795 |

**System specifications:**
- CPU: Intel Core i7-10700 (8 cores)
- RAM: 32 GB
- OS: Windows 10/11
- Compiler: GCC (MinGW-w64)

### Memory Usage

| Duration | Peak RAM | Disk Output |
|----------|----------|-------------|
| 30 days | ~150 MB | ~50 MB CSV |
| 365 days | ~180 MB | ~500 MB CSV |
| 730 days | ~200 MB | ~1 GB CSV |

### Comparison with 2D/3D Models

| Model | Grid Cells | 365-day Runtime | Speedup vs C-GEM |
|-------|------------|-----------------|------------------|
| **C-GEM** | **225 (1D)** | **11 min** | **1×** |
| MIKE21 | 50,000 (2D) | ~24 hours | 130× slower |
| Delft3D | 200,000 (3D) | ~1 week | 900× slower |
| FVCOM | 100,000 (3D) | ~3 days | 390× slower |

## Model Validation

### Salinity Intrusion

Comparison with Nguyen et al. (2008) field observations:

**Dry Season (March):**

| Branch | Model L₄ (km) | Observed L₄ (km) | Error |
|--------|---------------|------------------|-------|
| Ham Luong | 58 | 55 ± 5 | +5% |
| Co Chien | 48 | 45 ± 5 | +7% |
| Hau River | 65 | 62 ± 7 | +5% |
| My Tho | 52 | 50 ± 5 | +4% |

**Wet Season (September):**

| Branch | Model L₄ (km) | Observed L₄ (km) | Error |
|--------|---------------|------------------|-------|
| Ham Luong | 12 | 15 ± 5 | -20% |
| Co Chien | 8 | 10 ± 3 | -20% |
| Hau River | 18 | 20 ± 5 | -10% |
| My Tho | 10 | 12 ± 3 | -17% |

*L₄ = distance to 4 PSU isohaline from mouth*

### Tidal Range

| Location | Model | Observed | Error |
|----------|-------|----------|-------|
| Can Tho (Hau) | 2.2 m | 2.3 ± 0.3 m | -4% |
| My Tho city | 2.6 m | 2.5 ± 0.3 m | +4% |
| Vam Nao | 2.8 m | 2.7 ± 0.3 m | +4% |

### Dissolved Oxygen

Comparison with MRC monitoring stations:

| Station | Model O₂ | Observed O₂ | Error |
|---------|----------|-------------|-------|
| Can Tho | 190 µmol/L | 185 ± 25 | +3% |
| Downstream | 210 µmol/L | 205 ± 20 | +2% |
| Vam Nao | 235 µmol/L | 240 ± 20 | -2% |

### Skill Scores

| Variable | RMSE | Bias | Nash-Sutcliffe |
|----------|------|------|----------------|
| Salinity | 2.1 PSU | +0.3 | 0.92 |
| Water Level | 0.12 m | -0.02 | 0.95 |
| Discharge | 450 m³/s | +120 | 0.89 |
| O₂ | 18 µmol/L | +5 | 0.78 |

## Simplified Biogeochemistry Mode (80/20 Rule)

### Why Simplified Mode?

C-GEM includes a **"Simplified Mode"** (`SimplifiedMode = 1`) specifically designed for:

1. **Data-sparse tropical systems** where bacterial biomass, multi-pool organic carbon, and research-grade parameters cannot be measured
2. **Operational applications** requiring fast turnaround
3. **Uncertainty analysis** needing thousands of model runs

### What Simplified Mode Does

| Full RIVE Mode | Simplified Mode |
|----------------|-----------------|
| Bacterial dynamics (growth/mortality) | **Temperature-corrected first-order decay** |
| 6-pool organic carbon fractionation | **Direct TOC (single pool)** |
| Complex substrate limitation | **Monod O₂ limitation only** |
| 2-step nitrification (NH₄→NO₂→NO₃) | **Single-step (NH₄→NO₃)** |

### Mathematical Basis

**Aerobic degradation (BOD equivalent):**
$$\frac{d[TOC]}{dt} = -k_{ox} \cdot \theta^{T-20} \cdot [TOC] \cdot \frac{[O_2]}{[O_2] + K_{O_2}}$$

**Nitrification:**
$$\frac{d[NH_4]}{dt} = -k_{nit} \cdot \theta^{T-20} \cdot [NH_4] \cdot \frac{[O_2]}{[O_2] + K_{O_2}}$$

**Denitrification:**
$$\frac{d[NO_3]}{dt} = -k_{denit} \cdot \theta^{T-20} \cdot [NO_3] \cdot \frac{K_{iNO_2}}{[O_2] + K_{iNO_2}}$$

Where $\theta = 1.047$ (Arrhenius temperature coefficient).

### When to Use Each Mode

| Scenario | Recommended Mode |
|----------|------------------|
| **Research** with detailed lab data | Full RIVE |
| **Operational** water quality forecasting | **Simplified** |
| **Tropical systems** (Mekong, Amazon, Niger) | **Simplified** |
| **Calibration** and uncertainty analysis | **Simplified** |
| European rivers with RIVE heritage | Full RIVE |

### Performance Comparison

| Mode | 365-day Runtime | Species Updated | Memory |
|------|-----------------|-----------------|--------|
| Full RIVE | 15 min | All 30 | 200 MB |
| **Simplified** | **11 min** | 15 active | 180 MB |

## Seasonal Calibration

The Mekong has strong seasonal variation requiring dual-season calibration:

### seasonal_targets.csv

```csv
Branch,Variable,Location_km,Time_Day,Value,Weight
# Dry season peak (mid-March, Day 75)
Ham_Luong,SALINITY,30.0,75.0,25.0,2.0
Ham_Luong,SALINITY,40.0,75.0,15.0,2.0
Ham_Luong,SALINITY,50.0,75.0,6.0,2.0

# Wet season (September, Day 270)
Ham_Luong,SALINITY,20.0,270.0,0.5,1.0
Ham_Luong,SALINITY,30.0,270.0,0.1,1.0

# Oxygen year-round
Hau_River,O2,50.0,75.0,185.0,1.0
Hau_River,O2,50.0,270.0,220.0,1.0
```

### Running Seasonal Calibration

```powershell
.\bin\Debug\CGEM_Network.exe INPUT/Cases/Mekong_Delta_Full/case_config.txt --calibrate --stage 1
```

The optimizer finds parameters satisfying **both** dry and wet season observations.

## Output Files

### Directory Structure

```
OUTPUT/Mekong_Delta_Full/
├── CSV/
│   ├── Ham_Luong.csv        # Branch time series
│   ├── Co_Chien.csv
│   ├── Hau_River.csv
│   └── ...
├── NetCDF/
│   └── Mekong_Delta_Full.nc # Full 4D output
├── summary.txt               # Run summary
└── calibration_log.txt       # Calibration history (if --calibrate)
```

### CSV Format

```csv
time_days,depth_x1,depth_x10,depth_x30,salinity_x1,salinity_x10,salinity_x30,o2_x1,o2_x30,...
0.0,15.2,14.8,14.1,32.5,28.3,12.1,210,235,...
1.0,15.5,15.0,14.3,32.8,29.1,13.2,208,230,...
...
```

## Visualization

### Python Analysis Script

```python
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load full-year outputs
branches = ['Ham_Luong', 'Co_Chien', 'My_Tho', 'Hau_River']
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

for ax, branch in zip(axes.flat, branches):
    df = pd.read_csv(f'OUTPUT/Mekong_Delta_Full/CSV/{branch}.csv')
    
    # Plot salinity at different distances
    for x in [10, 30, 50]:
        col = f'salinity_x{x}'
        if col in df.columns:
            ax.plot(df['time_days'], df[col], label=f'{x} km')
    
    ax.set_xlabel('Day of Year')
    ax.set_ylabel('Salinity (PSU)')
    ax.set_title(branch)
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Mark dry/wet seasons
    ax.axvspan(0, 150, alpha=0.1, color='orange', label='Dry')
    ax.axvspan(150, 300, alpha=0.1, color='blue', label='Wet')

plt.tight_layout()
plt.savefig('mekong_salinity_annual.png', dpi=150)
plt.show()
```

### Seasonal Comparison Plot

```python
import xarray as xr
import matplotlib.pyplot as plt

# Load NetCDF output
ds = xr.open_dataset('OUTPUT/Mekong_Delta_Full/NetCDF/Mekong_Delta_Full.nc')

# Plot dry vs wet season salinity profiles
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

# Dry season (day 75 = mid-March)
ds.salinity.sel(time=75*86400, method='nearest').plot(ax=ax1, x='distance_km')
ax1.set_title('Dry Season (March) - Salinity Profiles')
ax1.set_xlabel('Distance from mouth (km)')
ax1.set_ylabel('Salinity (PSU)')

# Wet season (day 270 = late September)
ds.salinity.sel(time=270*86400, method='nearest').plot(ax=ax2, x='distance_km')
ax2.set_title('Wet Season (September) - Salinity Profiles')
ax2.set_xlabel('Distance from mouth (km)')
ax2.set_ylabel('Salinity (PSU)')

plt.tight_layout()
plt.savefig('mekong_seasonal_profiles.png', dpi=150)
```

## References

1. **Nguyen, A.D., Savenije, H.H.G., Pham, D.N., Tang, D.T. (2008)**. Using salt intrusion measurements to determine the freshwater discharge distribution over the branches of a multi-channel estuary: The Mekong Delta case. *Estuarine, Coastal and Shelf Science*, 77(3), 433-445.

2. **Savenije, H.H.G. (2005, 2012)**. *Salinity and Tides in Alluvial Estuaries*. Elsevier.

3. **MRC (Mekong River Commission)**. Hydrological data and monitoring.

4. **Garnier, J., Billen, G., et al. (2005)**. Modelling nitrogen cycling in river systems. *Biogeochemistry*, 77, 213-242.

5. **Alongi, D.M. (2014)**. Carbon cycling and storage in mangrove forests. *Annual Review of Marine Science*, 6, 195-219.
