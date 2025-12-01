# Quick Start

This guide walks you through running your first C-GEM simulation using the Mekong Delta case.

## 1. Build the Model

=== "Windows (PowerShell)"

    ```powershell
    cd C:\path\to\CGEM_network
    .\scripts\build.bat
    ```

=== "One-liner"

    ```powershell
    .\scripts\build-and-run.ps1 -r Mekong_Delta_Full
    ```

## 2. Run a Simulation

```powershell
.\bin\Debug\CGEM_Network.exe INPUT/Cases/Mekong_Delta_Full/case_config.txt
```

You'll see progress output:
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

Day    0.0 (  0.0%) [WARMUP] | H=15.00 m, U=0.100 m/s, S=0.1 | 0.0 s elapsed
Day    1.0 (  3.3%) [WARMUP] | H=15.40 m, U=0.142 m/s, S=0.1 | 2.4 s elapsed
...
Day   30.0 (100.0%)          | H=15.41 m, U=0.154 m/s, S=26.4 | 76.6 s elapsed

Simulation complete: 8640 steps in 76.6 seconds (112.8 steps/s)
CSV output written to: .\OUTPUT\Mekong_Delta_Full/CSV
```

## 3. View Results

### CSV Output

Results are written to `OUTPUT/Mekong_Delta_Full/CSV/`:

```
OUTPUT/Mekong_Delta_Full/CSV/
├── Tien_Main.csv
├── Hau_Main.csv
├── Co_Chien.csv
├── Ham_Luong.csv
├── Hau_River.csv
└── ...
```

Each CSV contains:
- Time series at hourly intervals
- Variables: depth, velocity, salinity, O2, nutrients, etc.
- One file per branch

### Quick Plot (Python)

```python
import pandas as pd
import matplotlib.pyplot as plt

# Load Ham Luong branch output
df = pd.read_csv('OUTPUT/Mekong_Delta_Full/CSV/Ham_Luong.csv')

# Plot salinity at different locations
fig, ax = plt.subplots(figsize=(10, 6))
ax.plot(df['time_days'], df['salinity_x10'], label='10 km from mouth')
ax.plot(df['time_days'], df['salinity_x30'], label='30 km from mouth')
ax.plot(df['time_days'], df['salinity_x50'], label='50 km from mouth')
ax.set_xlabel('Time (days)')
ax.set_ylabel('Salinity (PSU)')
ax.legend()
plt.show()
```

## 4. Run Calibration

C-GEM includes built-in calibration using NLopt:

```powershell
# Run full 3-stage calibration
.\bin\Debug\CGEM_Network.exe INPUT/Cases/Mekong_Delta_Full/case_config.txt --calibrate

# Run only Stage 1 (hydrodynamics)
.\bin\Debug\CGEM_Network.exe INPUT/Cases/Mekong_Delta_Full/case_config.txt --calibrate --stage 1

# With options
.\bin\Debug\CGEM_Network.exe INPUT/Cases/Mekong_Delta_Full/case_config.txt --calibrate --max-iter 50 --verbose 2
```

Output:
```
==============================================
  CALIBRATION MODE
==============================================

Loaded 30 calibration parameters from calibration_params.csv
Loaded 32 calibration objectives from calibration_targets.csv
Loaded 14 seasonal observation groups from seasonal_targets.csv

Optimizing 11 parameters for stage 1 using NLopt
NLopt algorithm: BOBYQA bound-constrained optimization

...

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

## 5. Configure Your Own Case

### Directory Structure

Create a new case in `INPUT/Cases/YourCase/`:

```
INPUT/Cases/YourCase/
├── case_config.txt           # Main configuration
├── topology.csv              # Network structure
├── boundary_map.csv          # Boundary conditions
├── biogeo_params.txt         # Biogeochemistry parameters
├── species_river.csv         # Upstream species concentrations
├── species_ocean.csv         # Ocean species concentrations
├── calibration_params.csv    # (optional) Parameters to calibrate
├── calibration_targets.csv   # (optional) Calibration objectives
└── seasonal_targets.csv      # (optional) Seasonal observations
```

### Minimal case_config.txt

```ini
CaseName = YourCase

# Files
Topology = topology.csv
BoundaryMap = boundary_map.csv
OutputDir = OUTPUT/YourCase

# Time
StartDate = 2024-01-01
Duration = 30
Warmup = 10

# Numerics
TimeStep = 300.0
DELXI = 2000.0

# Output
WriteCSV = 1
```

### Minimal topology.csv

```csv
BranchID,BranchName,NodeUp,NodeDown,Length_m,Width_Down_m,Width_Up_m,Depth_m,Chezy
1,Main_Channel,1,2,50000,500,300,10,55
```

### Minimal boundary_map.csv

```csv
NodeID,BoundaryType,ForcingFile
1,DISCHARGE,upstream_Q.csv
2,LEVEL,downstream_tide.csv
```

## Common Issues

### "No output files"

- Check `WriteCSV = 1` in case_config.txt
- Ensure `OutputDir` exists or can be created

### "Unstable simulation"

- Reduce time step: `TimeStep = 150.0`
- Check boundary conditions for sudden jumps
- Verify geometry (width ratios, depths)

### "Calibration not converging"

- Increase iterations: `--max-iter 100`
- Check parameter bounds in calibration_params.csv
- Verify observation data in calibration_targets.csv

## Next Steps

- [Data Preparation](../user-guide/data-preparation.md) - Prepare input files
- [Calibration Guide](../user-guide/calibration.md) - Detailed calibration workflow
- [Mekong Delta Case](../cases/mekong-delta.md) - Full case study
