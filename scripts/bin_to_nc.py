#!/usr/bin/env python3
"""
Binary to NetCDF Converter for C-GEM Network
---------------------------------------------
Reads .bin files with self-describing format (variable names in header)
and converts to a single combined NetCDF file.

Binary format (per branch):
  Header: M (int), num_hydro (int), num_species (int), num_reactions (int), dx (double)
          X grid (M doubles)
          Hydro names (num_hydro null-terminated strings)
          Species names (num_species null-terminated strings)
          Reaction names (num_reactions null-terminated strings)
  Per timestep: time (double)
                hydro[num_hydro][M]
                conc[num_species][M]
                reaction_rates[num_reactions][M]
"""

import xarray as xr
import numpy as np
import os
import sys
import struct


def _destagger_center(arr):
    """Interpolate staggered center values to even indices."""
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
    """Interpolate staggered velocity values to odd indices."""
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


def read_null_terminated_string(f):
    """Read a null-terminated string from binary file."""
    chars = []
    while True:
        c = f.read(1)
        if not c or c == b'\x00':
            break
        chars.append(c.decode('ascii'))
    return ''.join(chars)


def read_binary_file(filepath):
    """
    Read a .bin file with self-describing format.
    
    Returns: dict with keys:
        'M', 'dx', 'x', 'times',
        'hydro_names', 'hydro' (dict of arrays),
        'species_names', 'species' (dict of arrays),
        'reaction_names', 'reactions' (dict of arrays)
    """
    result = {
        'hydro_names': [],
        'hydro': {},
        'species_names': [],
        'species': {},
        'reaction_names': [],
        'reactions': {}
    }
    
    times = []
    hydro_data = {}
    species_data = {}
    reaction_data = {}

    with open(filepath, 'rb') as f:
        # Read header
        M = struct.unpack('i', f.read(4))[0]
        num_hydro = struct.unpack('i', f.read(4))[0]
        num_species = struct.unpack('i', f.read(4))[0]
        num_reactions = struct.unpack('i', f.read(4))[0]
        dx = struct.unpack('d', f.read(8))[0]

        result['M'] = M
        result['dx'] = dx
        
        # Read X grid
        x_buf = f.read(8 * M)
        if len(x_buf) != 8 * M:
            raise ValueError(f'Unexpected file size reading X grid: expected {8*M}, got {len(x_buf)}')
        x = np.frombuffer(x_buf, dtype=np.float64).copy()
        result['x'] = x

        # Read hydro names
        for _ in range(num_hydro):
            name = read_null_terminated_string(f)
            result['hydro_names'].append(name)
            hydro_data[name] = []

        # Read species names
        for _ in range(num_species):
            name = read_null_terminated_string(f)
            result['species_names'].append(name)
            species_data[name] = []

        # Read reaction names
        for _ in range(num_reactions):
            name = read_null_terminated_string(f)
            result['reaction_names'].append(name)
            reaction_data[name] = []

        # Read timestep data until EOF
        while True:
            time_bytes = f.read(8)
            if len(time_bytes) < 8:
                break
            time_val = struct.unpack('d', time_bytes)[0]
            times.append(time_val)

            # Read hydro variables
            for name in result['hydro_names']:
                buf = f.read(8 * M)
                if len(buf) != 8 * M:
                    raise ValueError(f'Unexpected file size reading hydro {name}')
                hydro_data[name].append(np.frombuffer(buf, dtype=np.float64).copy())

            # Read species concentrations
            for name in result['species_names']:
                conc_buf = f.read(8 * M)
                if len(conc_buf) != 8 * M:
                    raise ValueError(f'Unexpected file size reading species {name}')
                species_data[name].append(np.frombuffer(conc_buf, dtype=np.float64).copy())

            # Read reaction rates
            for name in result['reaction_names']:
                rate_buf = f.read(8 * M)
                if len(rate_buf) != 8 * M:
                    raise ValueError(f'Unexpected file size reading reaction {name}')
                reaction_data[name].append(np.frombuffer(rate_buf, dtype=np.float64).copy())

    # Convert to arrays and apply destaggering
    result['times'] = np.array(times)
    
    # Hydro variables - apply appropriate destaggering
    # Only velocity is on staggered grid
    for name in result['hydro_names']:
        arr = np.array(hydro_data[name])
        if name == 'velocity':
            result['hydro'][name] = _destagger_velocity(arr)
        else:
            result['hydro'][name] = _destagger_center(arr)

    for name in result['species_names']:
        result['species'][name] = _destagger_center(np.array(species_data[name]))
    
    for name in result['reaction_names']:
        result['reactions'][name] = _destagger_center(np.array(reaction_data[name]))

    # Trim ghost cells
    trim_cols = M - 3 if M >= 3 else 0
    if trim_cols > 0:
        result['x'] = result['x'][:trim_cols]
        for name in result['hydro_names']:
            result['hydro'][name] = result['hydro'][name][:, :trim_cols]
        for name in result['species_names']:
            result['species'][name] = result['species'][name][:, :trim_cols]
        for name in result['reaction_names']:
            result['reactions'][name] = result['reactions'][name][:, :trim_cols]
        result['M'] = trim_cols

    return result


def binary_to_netcdf(rundir, compress=True):
    """
    Convert all .bin files in rundir to a single combined NetCDF file.
    Variable names are read from binary headers, not hardcoded.
    """
    bin_files = [os.path.join(rundir, f) for f in os.listdir(rundir) if f.endswith('.bin')]
    if not bin_files:
        print('No .bin files found in output directory')
        return

    # Read all branches
    branches = []
    for bin_path in sorted(bin_files):
        branch_name = os.path.splitext(os.path.basename(bin_path))[0]
        print(f'Processing {branch_name}...')
        try:
            data = read_binary_file(bin_path)
            data['name'] = branch_name
            branches.append(data)
            print(f'  M={data["M"]}, hydro={len(data["hydro_names"])}, species={len(data["species_names"])}, reactions={len(data["reaction_names"])}')
        except Exception as e:
            print(f'  Error: {e}')

    if not branches:
        print('No branches loaded')
        return

    # Check time consistency
    times_ref = branches[0]['times']
    same_time = all(np.array_equal(b['times'], times_ref) for b in branches)
    if not same_time:
        print('Warning: Branches have different time axes')

    # Build combined dataset
    num_branches = len(branches)
    max_M = max(b['M'] for b in branches)
    time_axis = times_ref

    # Collect all variable names from all branches
    all_hydro = set()
    all_species = set()
    all_reactions = set()
    for b in branches:
        all_hydro.update(b['hydro_names'])
        all_species.update(b['species_names'])
        all_reactions.update(b['reaction_names'])
    
    # Prepare arrays filled with NaNs
    hydro_arr = {name: np.full((len(time_axis), num_branches, max_M), np.nan) for name in all_hydro}
    species_arr = {name: np.full((len(time_axis), num_branches, max_M), np.nan) for name in all_species}
    reaction_arr = {name: np.full((len(time_axis), num_branches, max_M), np.nan) for name in all_reactions}

    branch_names = []
    branch_x = np.full((num_branches, max_M), np.nan)
    
    for bi, b in enumerate(branches):
        branch_names.append(b['name'])
        m = b['M']
        branch_x[bi, :m] = b['x']
        
        for name in b['hydro_names']:
            hydro_arr[name][:, bi, :m] = b['hydro'][name]
        
        for name in b['species_names']:
            species_arr[name][:, bi, :m] = b['species'][name]
        
        for name in b['reaction_names']:
            reaction_arr[name][:, bi, :m] = b['reactions'][name]

    # Build xarray Dataset
    coords = {'time': time_axis, 'branch': branch_names}
    data_vars = {
        'branch_x': (('branch', 'x'), branch_x)
    }
    
    # Add hydro variables (no prefix)
    for name in sorted(all_hydro):
        data_vars[name] = (('time', 'branch', 'x'), hydro_arr[name])
    
    # Add species (no prefix)
    for name in sorted(all_species):
        data_vars[name] = (('time', 'branch', 'x'), species_arr[name])
    
    # Add reactions (with 'rate_' prefix to distinguish from species)
    for name in sorted(all_reactions):
        data_vars[f'rate_{name}'] = (('time', 'branch', 'x'), reaction_arr[name])

    ds = xr.Dataset(data_vars=data_vars, coords=coords)
    
    # Write NetCDF
    nc_path = os.path.join(rundir, 'combined_branches.nc')
    encoding = {var: {'zlib': compress, 'complevel': 4} for var in ds.data_vars}
    ds.to_netcdf(nc_path, format='NETCDF4', encoding=encoding)
    
    # Print complete variable summary
    print(f'\nSaved: {nc_path}')
    print(f'  Branches ({len(branch_names)}): {", ".join(branch_names)}')
    print(f'  Hydro ({len(all_hydro)}): {", ".join(sorted(all_hydro))}')
    print(f'  Species ({len(all_species)}): {", ".join(sorted(all_species))}')
    print(f'  Reactions ({len(all_reactions)}): {", ".join(sorted(all_reactions))}')


def main():
    if len(sys.argv) != 2:
        print("Usage: python bin_to_nc.py <output_dir>")
        sys.exit(1)

    rundir = sys.argv[1]

    if not os.path.isdir(rundir):
        print(f"Error: Directory '{rundir}' not found.")
        sys.exit(1)

    print(f"Converting binary files in: {rundir}")
    binary_to_netcdf(rundir)


if __name__ == "__main__":
    main()
