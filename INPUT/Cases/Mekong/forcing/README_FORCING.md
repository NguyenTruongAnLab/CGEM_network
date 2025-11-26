This folder contains forcing CSV files for the Mekong case used by the C-GEM Network model.

Files generated:
- TanChau_discharge.csv: Daily upstream discharge time series for 2017 (365 days). Header: Date,Q
- MyTho_tide.csv, HamLuong_tide.csv, CoChien_tide.csv, CanTho_tide.csv: Hourly tidal elevation time series for 2017 (8760 points). Header: DateTime,H

Notes:
- The model loader ignores timestamp strings and reads values sequentially with a fixed dt: 86400 s for discharge, 3600 s for tidal elevation.
- Header column names: The model expects a header row in each CSV. To avoid ambiguity we follow the repository examples exactly:
  - Discharge CSV header must be 'Date,Q' (two columns: a date string and a numeric discharge value)
  - Tide CSV header must be 'DateTime,H' (two columns: a datetime string and a numeric elevation value)
- There is a Python script `generate_mekong_forcing.py` in this directory which regenerates the files. By default it writes files for year 2017 with simple harmonic tidal components and seasonal discharge.

Usage:
  python generate_mekong_forcing.py

If you prefer to use a different year or parameters, edit the script.

The generated data is simplified (synthetic) and is intended for testing and development rather than high-fidelity simulations.

Validation:
- There is a validation script `scripts/validate_mekong_forcing.py` in the repository root that checks the header names and row counts for the generated files. Run it from the workspace root:

  python scripts\validate_mekong_forcing.py
