# Mekong Delta Case Study

## Overview

The Mekong Delta Full case (`Mekong_Delta_Full`) demonstrates C-GEM's capabilities for a complete multi-branch estuarine network.

### Network Structure

- **9 branches**: Main rivers and distributaries
- **10 nodes**: 2 upstream (discharge), 4 ocean (tidal), 4 junctions
- **Total length**: ~450 km of channels
- **Grid**: 2 km resolution, ~225 cells total

```
                 Tien River (Q)           Hau River (Q)
                      │                        │
                      ▼                        ▼
                 ┌─────────┐              ┌─────────┐
                 │Tien Main│              │Hau Main │
                 │ 48 km   │              │ 44 km   │
                 └────┬────┘              └────┬────┘
                      │                        │
                      └──────────┬─────────────┘
                                 │
                            ┌────┴────┐
                            │ Vam Nao │
                            │  8 km   │
                            └────┬────┘
                                 │
              ┌──────────────────┼──────────────────┐
              │                  │                  │
         ┌────┴─────┐       ┌────┴────┐       ┌────┴────┐
         │Tien Lower│       │Tien Conn│       │         │
         │  48 km   │       │  16 km  │       │         │
         └────┬─────┘       └────┬────┘       │         │
              │                  │            │         │
    ┌─────────┼─────────┐        │            │         │
    │         │         │        │            │         │
┌───┴───┐ ┌───┴───┐     │   ┌────┴────┐ ┌────┴────┐    │
│Co     │ │My Tho │     │   │Ham Luong│ │Hau River│    │
│Chien  │ │ 56 km │     │   │  60 km  │ │  88 km  │    │
│ 84 km │ └───┬───┘     │   └────┬────┘ └────┬────┘    │
└───┬───┘     │         │        │            │         │
    │         │         │        │            │         │
    ▼         ▼         ▼        ▼            ▼         │
  Ocean     Ocean     Ocean    Ocean        Ocean       │
```

## Configuration

### case_config.txt

```ini
CaseName = Mekong_Delta_Full

# Files
Topology = topology.csv
BoundaryMap = boundary_map.csv
BiogeoParams = biogeo_params.txt
OutputDir = OUTPUT/Mekong_Delta_Full

# Time
StartDate = 2024-03-01
Duration = 30
Warmup = 10

# Numerics
TimeStep = 300.0
DELXI = 2000.0

# Output
WriteCSV = 1
```

### topology.csv

```csv
BranchID,BranchName,NodeUp,NodeDown,Length_m,Width_Down_m,Width_Up_m,Depth_m,Chezy,RS
1,Tien_Main,1,3,48000,1200,1000,15,55,1.0
2,Hau_Main,2,3,44000,1100,900,14,55,1.0
3,Vam_Nao,3,4,8000,600,600,18,65,1.0
4,Tien_Lower,4,10,48000,1400,1200,14,55,1.0
5,Tien_Connector,10,5,16000,1200,1100,13,55,1.0
6,Co_Chien,10,6,84000,2800,1100,13,50,1.2
7,My_Tho,5,7,56000,2200,900,11,50,1.2
8,Ham_Luong,5,8,60000,2500,1000,12,50,1.2
9,Hau_River,4,9,88000,3000,1200,14,50,1.3
```

### Boundary Conditions

**Upstream (Discharge)**:
- Node 1: Tien River inflow (Q from MRC data)
- Node 2: Hau River inflow (Q from MRC data)

**Downstream (Tidal Level)**:
- Node 6: Co Chien mouth (M2 + S2 tides)
- Node 7: My Tho mouth
- Node 8: Ham Luong mouth
- Node 9: Hau River mouth

## Running the Case

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
  Time step   : 300.0 s
  Duration    : 30 days
  Warmup      : 10 days
==============================================

Day    1.0 (  3.3%) [WARMUP] | H=15.40 m, U=0.142 m/s, S=0.1 | 2.4 s elapsed
Day   10.0 ( 33.3%)          | H=15.39 m, U=0.173 m/s, S=11.4 | 22.4 s elapsed
...
Simulation complete: 8640 steps in 76.6 seconds (112.8 steps/s)
```

### With Calibration

```powershell
.\bin\Debug\CGEM_Network.exe INPUT/Cases/Mekong_Delta_Full/case_config.txt --calibrate --stage 1 --max-iter 50
```

## Results

### Salinity Intrusion

During dry season (low discharge):

| Branch | Intrusion Length (4 PSU) | Observed |
|--------|-------------------------|----------|
| Ham Luong | 55-60 km | 50-65 km |
| Co Chien | 45-50 km | 40-55 km |
| Hau River | 60-70 km | 55-75 km |

### Tidal Range

| Location | Model | Observed |
|----------|-------|----------|
| Can Tho | 2.0-2.5 m | 2.0-2.8 m |
| My Tho | 2.5-3.0 m | 2.2-3.2 m |

### Typical Output

```
time_days,depth_x1,depth_x10,salinity_x1,salinity_x10,salinity_x30,...
0.0,15.2,14.8,32.5,28.3,12.1,...
0.042,15.5,15.0,32.8,29.1,13.2,...
0.083,15.1,14.5,31.9,27.5,11.8,...
...
```

## Seasonal Calibration

The Mekong has strong seasonal variation:

- **Dry season** (Dec-Apr): Low discharge, maximum salinity intrusion
- **Wet season** (Jun-Oct): High discharge, minimal intrusion

### seasonal_targets.csv

```csv
Branch,Variable,Location_km,Time_Day,Value,Weight
# Dry season peak (mid-March)
Ham_Luong,SALINITY,30.0,15.0,25.0,2.0
Ham_Luong,SALINITY,40.0,15.0,15.0,2.0
Ham_Luong,SALINITY,50.0,15.0,6.0,2.0

# Wet season (August)
Ham_Luong,SALINITY,20.0,180.0,0.5,1.0
Ham_Luong,SALINITY,30.0,180.0,0.1,1.0
```

### Running Seasonal Calibration

```powershell
.\bin\Debug\CGEM_Network.exe INPUT/Cases/Mekong_Delta_Full/case_config.txt --calibrate
```

The optimizer will find parameters that satisfy **both** dry and wet season observations.

## Visualization

### Python Script

```python
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load branch outputs
branches = ['Ham_Luong', 'Co_Chien', 'My_Tho', 'Hau_River']
fig, axes = plt.subplots(2, 2, figsize=(14, 10))

for ax, branch in zip(axes.flat, branches):
    df = pd.read_csv(f'OUTPUT/Mekong_Delta_Full/CSV/{branch}.csv')
    
    # Plot salinity at different distances
    for x in [10, 30, 50]:
        col = f'salinity_x{x}'
        if col in df.columns:
            ax.plot(df['time_days'], df[col], label=f'{x} km')
    
    ax.set_xlabel('Time (days)')
    ax.set_ylabel('Salinity (PSU)')
    ax.set_title(branch)
    ax.legend()
    ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('mekong_salinity.png', dpi=150)
plt.show()
```

### Expected Salinity Profiles

Typical dry season longitudinal profile:

```
Salinity (PSU)
    35 ┤                                            ╭──────
       │                                        ╭───╯
    30 ┤                                    ╭───╯
       │                                ╭───╯
    25 ┤                            ╭───╯
       │                        ╭───╯
    20 ┤                    ╭───╯
       │                ╭───╯
    15 ┤            ╭───╯
       │        ╭───╯
    10 ┤    ╭───╯
       │╭───╯
     5 ┼╯
       │
     0 ┼─────────────────────────────────────────────────
       0    10    20    30    40    50    60    70    80
                        Distance from mouth (km)
```

## References

1. Nguyen, A.D., Savenije, H.H.G., Pham, D.N., Tang, D.T. (2008). Using salt intrusion measurements to determine the freshwater discharge distribution over the branches of a multi-channel estuary: The Mekong Delta case. *Estuarine, Coastal and Shelf Science*, 77(3), 433-445.

2. Renaud, F.G., Kuenzer, C. (2012). *The Mekong Delta System*. Springer.

3. MRC (Mekong River Commission). Hydrological data and monitoring.
