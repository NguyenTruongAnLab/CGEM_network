#!/usr/bin/env python3

import xarray as xr
import numpy as np
import os
import sys
import struct

def _destagger_center(arr):
    if arr is None or arr.size == 0:
        return arr
    arr = arr.copy()
    M = arr.shape[-1]
    if M < 3:
        return arr
    idx = np.arange(1, M + 1)
    even_idx = np.where((idx % 2) == 0)[0]
    even_idx = even_idx[(even_idx > 0) & (even_idx < M - 1)]
    if arr.ndim == 2:
        arr[:, even_idx] = 0.5 * (arr[:, even_idx - 1] + arr[:, even_idx + 1])
    elif arr.ndim == 3:
        arr[:, :, even_idx] = 0.5 * (arr[:, :, even_idx - 1] + arr[:, :, even_idx + 1])
    return arr

def _destagger_velocity(arr):
    if arr is None or arr.size == 0:
        return arr
    arr = arr.copy()
    M = arr.shape[-1]
    if M < 3:
        return arr
    idx = np.arange(1, M + 1)
    odd_idx = np.where((idx % 2) == 1)[0]
    odd_idx = odd_idx[(odd_idx > 0) & (odd_idx < M - 1)]
    if arr.ndim == 2:
        arr[:, odd_idx] = 0.5 * (arr[:, odd_idx - 1] + arr[:, odd_idx + 1])
    elif arr.ndim == 3:
        arr[:, :, odd_idx] = 0.5 * (arr[:, :, odd_idx - 1] + arr[:, :, odd_idx + 1])
    return arr

def read_binary_file(filepath):
    """Reads a .bin file with streaming format."""
    times = []
    velocities = []
    depths = []
    dispersions = []
    concs = []

    with open(filepath, 'rb') as f:
        # Read header
        M = struct.unpack('i', f.read(4))[0]
        num_sp = struct.unpack('i', f.read(4))[0]
        dx = struct.unpack('d', f.read(8))[0]

        # Read X grid
        x_buf = f.read(8 * M)
        if len(x_buf) != 8 * M:
            raise ValueError(f'Unexpected file size while reading X grid: expected {8*M} bytes, got {len(x_buf)}')
        x = np.frombuffer(x_buf, dtype=np.float64)

        # Read time steps until EOF
        while True:
            time_bytes = f.read(8)
            if len(time_bytes) < 8:
                break
            time_val = struct.unpack('d', time_bytes)[0]
            times.append(time_val)

            # Read velocity
            vel_buf = f.read(8 * M)
            if len(vel_buf) != 8 * M:
                raise ValueError(f'Unexpected file size while reading velocity: expected {8*M} bytes, got {len(vel_buf)}')
            vel = np.frombuffer(vel_buf, dtype=np.float64)
            velocities.append(vel)

            # Read depth
            dep_buf = f.read(8 * M)
            if len(dep_buf) != 8 * M:
                raise ValueError(f'Unexpected file size while reading depth: expected {8*M} bytes, got {len(dep_buf)}')
            dep = np.frombuffer(dep_buf, dtype=np.float64)
            depths.append(dep)

            # Read dispersion
            disp_buf = f.read(8 * M)
            if len(disp_buf) != 8 * M:
                raise ValueError(f'Unexpected file size while reading dispersion: expected {8*M} bytes, got {len(disp_buf)}')
            disp = np.frombuffer(disp_buf, dtype=np.float64)
            dispersions.append(disp)

            # Read concentrations
            step_conc = {}
            species_names = ['salinity', 'phy1', 'phy2', 'dsi', 'no3', 'nh4', 'po4', 'o2', 'toc', 'spm', 'dic', 'at', 'pco2', 'co2', 'ph', 'hs', 'alkc']
            for sp in range(num_sp):
                name = species_names[sp] if sp < len(species_names) else f'sp{sp}'
                conc_buf = f.read(8 * M)
                if len(conc_buf) != 8 * M:
                    raise ValueError(f'Unexpected file size while reading species {sp}: expected {8*M} bytes, got {len(conc_buf)}')
                step_conc[name] = np.frombuffer(conc_buf, dtype=np.float64)
            concs.append(step_conc)

    # Convert to arrays
    times = np.array(times)
    velocities = np.array(velocities)
    depths = np.array(depths)
    dispersions = np.array(dispersions)

    conc_dict = {}
    for name in concs[0].keys():
        conc_dict[name] = np.array([step[name] for step in concs])

    depths = _destagger_center(depths)
    velocities = _destagger_velocity(velocities)
    dispersions = _destagger_center(dispersions)
    for name, arr in conc_dict.items():
        conc_dict[name] = _destagger_center(arr)

    trim_cols = M - 3 if M >= 3 else 0
    if trim_cols > 0:
        x = x[:trim_cols]
        depths = depths[:, :trim_cols]
        velocities = velocities[:, :trim_cols]
        dispersions = dispersions[:, :trim_cols]
        for name, arr in conc_dict.items():
            conc_dict[name] = arr[:, :trim_cols]
    else:
        x = np.array([], dtype=np.float64)
        depths = depths[:, :0]
        velocities = velocities[:, :0]
        dispersions = dispersions[:, :0]
        for name in conc_dict:
            conc_dict[name] = conc_dict[name][:, :0]

    return times, x, depths, velocities, None, dispersions, conc_dict

def binary_to_netcdf(rundir, combined=True, compress=True):
    """Convert .bin files to NetCDF. By default creates a single combined netcdf.

    If combined=False, write per-branch netcdf files (backward compatible).
    """
    bin_files = [os.path.join(rundir, f) for f in os.listdir(rundir) if f.endswith('.bin')]
    if not bin_files:
        print('No .bin files found in output dir')
        return

    # Read all bin files
    branches = []
    for bin_path in bin_files:
        branch_name = os.path.splitext(os.path.basename(bin_path))[0]
        print(f'Processing {branch_name}')
        try:
            times, x, depth, velocity, waterlevel, dispersion, conc = read_binary_file(bin_path)
            branches.append({
                'name': branch_name,
                'M': x.size,
                'x': x,
                'time': times,
                'depth': depth,
                'velocity': velocity,
                'dispersion': dispersion,
                'conc': conc
            })
        except Exception as e:
            print(f'Error processing {bin_path}: {e}')

    if not combined:
        # Per-branch write
        for b in branches:
            coords = {'time': b['time'], 'x': b['x']}
            data_vars = {
                'depth': (('time', 'x'), b['depth']),
                'velocity': (('time', 'x'), b['velocity']),
                'dispersion': (('time', 'x'), b['dispersion'])
            }
            data_vars.update({name: (('time', 'x'), data) for name, data in b['conc'].items()})
            ds = xr.Dataset(data_vars=data_vars, coords=coords)
            nc_path = os.path.join(rundir, f"{b['name']}.nc")
            ds.to_netcdf(nc_path, format='NETCDF4', encoding={var: {'zlib': compress} for var in ds.data_vars})
            print(f'Saved {b["name"]}.nc')
        return

    # Combined netcdf: check time consistency across branches
    times_ref = branches[0]['time']
    same_time = all(np.array_equal(b['time'], times_ref) for b in branches)
    if not same_time:
        print('Warning: Not all branches have the same time axis. Falling back to per-branch write.')
        binary_to_netcdf(rundir, combined=False, compress=compress)
        return

    num_branches = len(branches)
    max_M = max(b['M'] for b in branches)
    M_axis = np.arange(max_M)
    time_axis = times_ref

    # Prepare arrays filled with NaNs
    depth_arr = np.full((time_axis.size, num_branches, max_M), np.nan)
    vel_arr = np.full((time_axis.size, num_branches, max_M), np.nan)
    disp_arr = np.full((time_axis.size, num_branches, max_M), np.nan)
    species_names = list(branches[0]['conc'].keys())
    conc_arr = {name: np.full((time_axis.size, num_branches, max_M), np.nan) for name in species_names}

    branch_names = []
    branch_x = np.full((num_branches, max_M), np.nan)
    for bi, b in enumerate(branches):
        branch_names.append(b['name'])
        m = b['M']
        branch_x[bi, :m] = b['x']
        depth_arr[:, bi, :m] = b['depth']
        vel_arr[:, bi, :m] = b['velocity']
        disp_arr[:, bi, :m] = b['dispersion']
        for name in species_names:
            conc_arr[name][:, bi, :m] = b['conc'][name]

    coords = {'time': time_axis, 'branch': branch_names, 'x': (('branch', 'x'), branch_x)}
    data_vars = {
        'depth': (('time', 'branch', 'x'), depth_arr),
        'velocity': (('time', 'branch', 'x'), vel_arr),
        'dispersion': (('time', 'branch', 'x'), disp_arr),
        'branch_x': (('branch', 'x'), branch_x)
    }
    for name in species_names:
        data_vars[name] = (('time', 'branch', 'x'), conc_arr[name])

    ds = xr.Dataset(data_vars=data_vars, coords=coords)
    nc_path = os.path.join(rundir, 'combined_branches.nc')
    # Use some compression
    encoding = {var: {'zlib': compress, 'complevel': 4} for var in ds.data_vars}
    ds.to_netcdf(nc_path, format='NETCDF4', encoding=encoding)
    print(f'Saved combined NetCDF: {nc_path}')

def main():
    if len(sys.argv) != 2:
        print("Usage: python dat_to_nc.py <output_dir>")
        sys.exit(1)

    rundir = sys.argv[1]

    if not os.path.isdir(rundir):
        print(f"Error: Directory '{rundir}' not found.")
        sys.exit(1)

    print(f"Converting binary files in: {rundir}")
    binary_to_netcdf(rundir)

if __name__ == "__main__":
    main()