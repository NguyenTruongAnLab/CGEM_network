# Quick Start — Windows (PowerShell) — Minimal

Install MinGW via MSYS2 if needed:

```powershell
& "C:\msys64\usr\bin\pacman.exe" -S --needed mingw-w64-x86_64-gcc mingw-w64-x86_64-make
```
Run the helper script (build + run default case):

```powershell
.\scripts\build-and-run.ps1 -r SaigonDongNai
```
Or for Mekong case:

```powershell
.\scripts\build-and-run.ps1 -r Mekong
```
You can customize time-step and grid spacing in the case config using `DELTI` (seconds) and `DELXI` (meters). For example:

```ini
DELTI = 300
DELXI = 2000
```
Manual commands (alternate):
```powershell
Set-Location "C:\Users\nguytruo\Documents\Github\CGEM_Network"
"C:\msys64\mingw64\bin\mingw32-make.exe"
.\bin\Debug\CGEM_Network.exe INPUT/Cases/SaigonDongNai/case_config.txt
```
Tip: add `C:\msys64\mingw64\bin` to PATH for convenience.

## NetCDF Output Structure
The model writes custom binary files (.bin) for high performance, which are then converted to NetCDF using Python.

For each branch, a .bin file contains:
- Header: num_time_steps, num_x_points, num_species
- Time array (double)
- X coordinates (double)
- Depth, velocity, waterlevel, dispersion arrays (double, time × x)
- Concentration arrays for each species (double, time × x)

To convert to NetCDF:
```powershell
python scripts\bin_to_nc.py OUTPUT/CaseName
```

This creates branch_name.nc files with CF-compliant NetCDF format containing full spatiotemporal data.
