# C-GEM Network

This repository contains the C implementation of the C-GEM Network Engine, including a multi-branch staggered-grid hydrodynamics solver, transport equations, and biogeochemical processes.

## Quickstart - Build and Run (Windows)

Build and run the `Tien_River` case using the included script (requires msys2/mingw)

```powershell
# Build + run Tien_River
.
cd scripts
./build-and-run.ps1 -r Tien_River
```

Outputs are written by default to `OUTPUT/Tien_River/` and include CSV and NetCDF; plotting scripts generate printable images in `OUTPUT/Tien_River/NETCDF_PLOTS/`.

## License