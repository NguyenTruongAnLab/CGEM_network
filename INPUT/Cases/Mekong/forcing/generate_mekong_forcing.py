#!/usr/bin/env python3
"""
Generate year-long forcing CSV files for the Mekong case:
- 1 year of daily discharge at TanChau (TanChau_discharge.csv)
- 1 year of hourly tidal elevation for boundary nodes (MyTho_tide.csv, HamLuong_tide.csv, CoChien_tide.csv, CanTho_tide.csv)

Files follow the CGEM network expected format where the Csv header is present
and subsequent rows are 'Date/Q' or 'DateTime/H' values. The loader ignores
the date/time strings and uses the dt = 86400 for discharge and dt=3600 for tidal (hourly).

This script writes CSV files for 2017 (non-leap year with 365 days, 8760 hours).
"""
from datetime import datetime, timedelta
import os
import math
import random

# Settings
year = 2017
start = datetime(year, 1, 1, 0, 0, 0)
end = datetime(year, 12, 31, 23, 0, 0)

# Ensure files are written in the forcing folder
forcing_dir = os.path.dirname(__file__)

# Discharge params (daily)
discharge_mean = 3500.0
discharge_amp = 1400.0  # seasonality amplitude
seed = 42

# Tidal station parameters (M2 and K1 components amplitude, mean, phase shift)
stations = {
    'MyTho_tide.csv': {'mean': 0.0, 'A_M2': 1.8, 'A_K1': 0.6, 'phi_M2': 0.2, 'phi_K1': 1.0},
    'HamLuong_tide.csv': {'mean': 0.0, 'A_M2': 1.7, 'A_K1': 0.5, 'phi_M2': -0.5, 'phi_K1': 0.4},
    'CoChien_tide.csv': {'mean': 0.0, 'A_M2': 2.2, 'A_K1': 0.8, 'phi_M2': 1.0, 'phi_K1': -0.3},
    'CanTho_tide.csv': {'mean': 0.0, 'A_M2': 1.5, 'A_K1': 0.4, 'phi_M2': 0.7, 'phi_K1': -0.1},
}

# Derived constants
M2_period_hours = 12.42
K1_period_hours = 23.93

random.seed(seed)

def write_discharge_file(path, start_dt, end_dt, mean, amp, seed=None):
    if seed is not None:
        random.seed(seed)
    d = start_dt
    with open(path, 'w', newline='') as f:
        # Use the exact header expected by the model
        f.write('Date,Q\n')
        while d <= end_dt:
            doy = (d - datetime(start_dt.year, 1, 1)).days + 1
            # seasonal sinusoid: max at the middle of the year (doy=183)
            angle = 2 * math.pi * (doy / 365.0)
            Q = mean + amp * math.sin(angle)
            # small daily variability
            Q += 50.0 * math.sin(2*math.pi*((doy*7) % 365) / 365.0)  # weekly like modulation
            # small random noise
            Q += random.gauss(0, 20)
            # ensure Q is not negative
            Q = max(Q, 0.0)
            f.write('{:%Y-%m-%d},{:.3f}\n'.format(d.date(), Q))
            d = d + timedelta(days=1)

# Create tide files (hourly)
def write_tide_file(path, start_dt, end_dt, params, seed=None):
    if seed is not None:
        random.seed(seed)
    with open(path, 'w', newline='') as f:
        # Use the exact header expected by the model
        f.write('DateTime,H\n')
        t = start_dt
        idx = 0
        while t <= end_dt:
            # t in hours since start
            th = idx
            # harmonic sum
            h = params['mean']
            h += params['A_M2'] * math.sin(2 * math.pi * th / M2_period_hours + params['phi_M2'])
            h += params['A_K1'] * math.sin(2 * math.pi * th / K1_period_hours + params['phi_K1'])
            # small long period modulation (spring-neap ~14.77 days)
            h += 0.3 * math.sin(2 * math.pi * (th / (24.0*14.77)))
            # small noise
            h += random.gauss(0, 0.05)
            # write timestamp and height
            f.write('{:%Y-%m-%d %H:%M},{:.6f}\n'.format(t, h))
            t = t + timedelta(hours=1)
            idx += 1

def validate_headers_forcing(dirpath):
    """Simple validation to ensure header names are exactly as required."""
    discharge_path = os.path.join(dirpath, 'TanChau_discharge.csv')
    tide_files = ['MyTho_tide.csv', 'HamLuong_tide.csv', 'CoChien_tide.csv', 'CanTho_tide.csv']
    # Check discharge header
    with open(discharge_path, 'r') as f:
        header = f.readline().strip()
        if header != 'Date,Q':
            raise RuntimeError(f'Unexpected discharge header: {header} (expected "Date,Q")')
    # Check tide headers
    for fname in tide_files:
        with open(os.path.join(dirpath, fname), 'r') as f:
            header = f.readline().strip()
            if header != 'DateTime,H':
                raise RuntimeError(f'Unexpected tide header in {fname}: {header} (expected "DateTime,H")')
    print('Header validation passed')


def main():
    # Write discharge
    write_discharge_file(os.path.join(forcing_dir, 'TanChau_discharge.csv'), start, end, discharge_mean, discharge_amp, seed)
    # Write tide files
    for fname, params in stations.items():
        write_tide_file(os.path.join(forcing_dir, fname), start, end, params, seed)
    # Validate
    validate_headers_forcing(forcing_dir)
    print('Generated TanChau_discharge.csv and tide files for year', year)

if __name__ == '__main__':
    main()
