# Output Analysis

## Output File Structure

### CSV Output

Each branch has a CSV file with time series at multiple locations:

```
OUTPUT/CaseName/CSV/
├── Branch1.csv
├── Branch2.csv
└── ...
```

### CSV Format

```csv
time_days,depth_x1,depth_x5,depth_x10,...,salinity_x1,salinity_x5,...
0.000,15.2,15.1,15.0,...,32.5,28.3,...
0.042,15.5,15.3,15.1,...,32.8,29.1,...
0.083,15.1,14.9,14.8,...,31.9,27.5,...
```

**Columns**:
- `time_days`: Simulation time [days]
- `variable_xN`: Variable at grid index N (where N=1 is downstream)

### Variables in Output

| Category | Variables |
|----------|-----------|
| Hydrodynamics | depth, velocity, waterlevel, area, width, dispersion |
| Salinity | salinity |
| Phytoplankton | phy1, phy2 |
| Nutrients | dsi, no3, nh4, po4, no2 |
| Oxygen | o2 |
| Carbon | toc, dic, at, pco2, co2, ph |
| Sediment | spm |
| GHG | n2o, ch4 |

## Loading Data in Python

### Basic Loading

```python
import pandas as pd
import numpy as np

# Load single branch
df = pd.read_csv('OUTPUT/Mekong_Delta_Full/CSV/Ham_Luong.csv')

# Available columns
print(df.columns.tolist())

# Time series at specific location
time = df['time_days']
sal_30km = df['salinity_x15']  # x15 = index 15 ≈ 30 km from mouth
```

### Loading All Branches

```python
import os

output_dir = 'OUTPUT/Mekong_Delta_Full/CSV'
branches = {}

for f in os.listdir(output_dir):
    if f.endswith('.csv'):
        name = f.replace('.csv', '')
        branches[name] = pd.read_csv(os.path.join(output_dir, f))

# Access specific branch
ham_luong = branches['Ham_Luong']
```

## Visualization

### Time Series Plot

```python
import matplotlib.pyplot as plt

fig, axes = plt.subplots(2, 2, figsize=(14, 10))

# Salinity at different locations
ax = axes[0, 0]
for x in [5, 15, 25]:
    col = f'salinity_x{x}'
    if col in df.columns:
        ax.plot(df['time_days'], df[col], label=f'{x*2} km')
ax.set_xlabel('Time (days)')
ax.set_ylabel('Salinity (PSU)')
ax.set_title('Salinity Time Series')
ax.legend()
ax.grid(True, alpha=0.3)

# Velocity
ax = axes[0, 1]
ax.plot(df['time_days'], df['velocity_x10'])
ax.set_xlabel('Time (days)')
ax.set_ylabel('Velocity (m/s)')
ax.set_title('Velocity at 20 km')
ax.grid(True, alpha=0.3)

# Water level
ax = axes[1, 0]
ax.plot(df['time_days'], df['waterlevel_x5'])
ax.set_xlabel('Time (days)')
ax.set_ylabel('Water Level (m)')
ax.set_title('Water Level at 10 km')
ax.grid(True, alpha=0.3)

# SPM
ax = axes[1, 1]
for x in [5, 15, 25]:
    col = f'spm_x{x}'
    if col in df.columns:
        ax.plot(df['time_days'], df[col], label=f'{x*2} km')
ax.set_xlabel('Time (days)')
ax.set_ylabel('SPM (mg/L)')
ax.set_title('Suspended Sediment')
ax.legend()
ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('results_overview.png', dpi=150)
```

### Longitudinal Profile

```python
def get_longitudinal_profile(df, variable, time_idx):
    """Extract longitudinal profile at specific time"""
    cols = [c for c in df.columns if c.startswith(f'{variable}_x')]
    values = df.iloc[time_idx][[c for c in cols]].values
    
    # Extract x indices
    x_indices = [int(c.split('_x')[1]) for c in cols]
    x_km = np.array(x_indices) * 2  # Assuming 2 km grid
    
    return x_km, values

# Plot profile at day 15
time_idx = int(15 / 0.042)  # Convert to row index
x, sal = get_longitudinal_profile(df, 'salinity', time_idx)

plt.figure(figsize=(10, 6))
plt.plot(x, sal, 'b-', linewidth=2)
plt.xlabel('Distance from mouth (km)')
plt.ylabel('Salinity (PSU)')
plt.title('Salinity Profile at Day 15')
plt.grid(True, alpha=0.3)
plt.savefig('salinity_profile.png', dpi=150)
```

### Multi-Branch Comparison

```python
fig, ax = plt.subplots(figsize=(12, 6))

branches_to_plot = ['Ham_Luong', 'Co_Chien', 'My_Tho', 'Hau_River']

for branch in branches_to_plot:
    df = pd.read_csv(f'OUTPUT/Mekong_Delta_Full/CSV/{branch}.csv')
    # Plot salinity at 30 km
    col = 'salinity_x15'
    if col in df.columns:
        ax.plot(df['time_days'], df[col], label=branch)

ax.set_xlabel('Time (days)')
ax.set_ylabel('Salinity at 30 km (PSU)')
ax.legend()
ax.grid(True, alpha=0.3)
plt.savefig('branch_comparison.png', dpi=150)
```

## NetCDF Conversion

### Convert Binary to NetCDF

```powershell
python scripts/bin_to_nc.py OUTPUT/Mekong_Delta_Full
```

### Load NetCDF

```python
import netCDF4 as nc

ds = nc.Dataset('OUTPUT/Mekong_Delta_Full/Mekong_Delta_Full.nc')

# List variables
print(ds.variables.keys())

# Extract data
time = ds.variables['time'][:]
x = ds.variables['x'][:]
salinity = ds.variables['salinity'][:, :]  # (time, x)
```

## Statistical Analysis

### Basic Statistics

```python
# Time-averaged profile
sal_mean = df[[c for c in df.columns if c.startswith('salinity_x')]].mean()

# Tidal range
depth_cols = [c for c in df.columns if c.startswith('depth_x')]
tidal_range = df[depth_cols].max() - df[depth_cols].min()

# Salt intrusion length (4 PSU)
def calc_intrusion_length(df, threshold=4.0, dx=2.0):
    """Calculate salt intrusion length"""
    sal_cols = sorted([c for c in df.columns if c.startswith('salinity_x')],
                      key=lambda x: int(x.split('_x')[1]))
    
    intrusions = []
    for _, row in df.iterrows():
        for i, col in enumerate(sal_cols):
            if row[col] < threshold:
                intrusions.append(i * dx)
                break
        else:
            intrusions.append(len(sal_cols) * dx)
    
    return intrusions

intrusion = calc_intrusion_length(df)
print(f"Mean intrusion: {np.mean(intrusion):.1f} km")
print(f"Max intrusion: {np.max(intrusion):.1f} km")
```

### Tidal Analysis

```python
from scipy import signal

# Extract tidal signal
depth = df['depth_x10'].values
time_hours = df['time_days'].values * 24

# Find peaks (high tide)
peaks, _ = signal.find_peaks(depth, distance=10)
troughs, _ = signal.find_peaks(-depth, distance=10)

# Calculate tidal range
tidal_ranges = []
for i in range(min(len(peaks), len(troughs))):
    tr = depth[peaks[i]] - depth[troughs[i]]
    tidal_ranges.append(abs(tr))

print(f"Mean tidal range: {np.mean(tidal_ranges):.2f} m")
```

## Export for Other Tools

### Export to Excel

```python
with pd.ExcelWriter('results.xlsx') as writer:
    for branch, df in branches.items():
        df.to_excel(writer, sheet_name=branch, index=False)
```

### Export Summary Statistics

```python
summary = []
for branch, df in branches.items():
    sal_cols = [c for c in df.columns if c.startswith('salinity_x')]
    summary.append({
        'Branch': branch,
        'Mean_Salinity': df[sal_cols].mean().mean(),
        'Max_Salinity': df[sal_cols].max().max(),
        'Mean_Depth': df[[c for c in df.columns if c.startswith('depth_x')]].mean().mean()
    })

pd.DataFrame(summary).to_csv('summary_stats.csv', index=False)
```
