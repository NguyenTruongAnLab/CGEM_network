# Running Simulations

## Basic Usage

```powershell
.\bin\Debug\CGEM_Network.exe <config_file> [options]
```

### Examples

```powershell
# Standard simulation
.\bin\Debug\CGEM_Network.exe INPUT/Cases/Mekong_Delta_Full/case_config.txt

# With calibration
.\bin\Debug\CGEM_Network.exe INPUT/Cases/Mekong_Delta_Full/case_config.txt --calibrate

# Calibration with options
.\bin\Debug\CGEM_Network.exe INPUT/Cases/Mekong_Delta_Full/case_config.txt \
    --calibrate --stage 1 --max-iter 100 --verbose 2
```

## Command-Line Options

| Option | Description | Default |
|--------|-------------|---------|
| `--calibrate` | Enable calibration mode | Off |
| `--stage N` | Run only stage N (1-3) | All stages |
| `--max-iter N` | Max optimizer iterations | 100 |
| `--verbose N` | Verbosity (0-2) | 1 |

## Output

### Console Progress

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

### Output Files

```
OUTPUT/Mekong_Delta_Full/
├── CSV/
│   ├── Tien_Main.csv
│   ├── Hau_Main.csv
│   ├── Vam_Nao.csv
│   ├── Co_Chien.csv
│   ├── Ham_Luong.csv
│   └── ...
└── calibration/  (if --calibrate)
    └── calibration_results.csv
```

[Learn more](user-guide/output.md) - Detailed output file formats and post-processing

## Helper Scripts

### Build and Run

```powershell
# Build + run in one command
.\scripts\build-and-run.ps1 -r Mekong_Delta_Full

# Build + run different case
.\scripts\build-and-run.ps1 -r Tien_River
```

### Convert to NetCDF

```powershell
python scripts/bin_to_nc.py OUTPUT/Mekong_Delta_Full
```

### Quick Plot

```powershell
python scripts/plot_netcdf.py OUTPUT/Mekong_Delta_Full/Mekong_Delta_Full.nc
```

## Performance Tips

### Time Step Selection

| Condition | Recommended dt |
|-----------|---------------|
| Deep channels (>10m) | 300-600 s |
| Shallow areas (<5m) | 150-300 s |
| Strong tides (>3m) | 150-300 s |
| Calibration runs | 300 s |

### Grid Spacing

| Application | DELXI |
|-------------|-------|
| Research | 500-1000 m |
| Calibration | 2000 m |
| Screening | 5000 m |

### Expected Performance

| Grid Size | Time Step | Speed |
|-----------|-----------|-------|
| 100 cells | 300 s | ~500 steps/s |
| 225 cells | 300 s | ~110 steps/s |
| 500 cells | 300 s | ~50 steps/s |

## Troubleshooting

### "Simulation unstable"

1. Reduce time step
2. Check boundary conditions for jumps
3. Verify geometry (width ratios)

### "No output files"

1. Check `WriteCSV = 1`
2. Verify `OutputDir` exists

### "NaN values in output"

1. Check initial conditions
2. Verify species boundary files
3. Reduce biogeo time step
