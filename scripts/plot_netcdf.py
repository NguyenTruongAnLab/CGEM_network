#!/usr/bin/env python3
"""
Comprehensive NetCDF Plotter for C-GEM
--------------------------------------
1. Generates per-branch analysis (Hovmoller, Timeseries, Profiles) for ALL variables.
2. Generates "Schematic Maps" with:
   - Sankey-style non-overlapping junctions (Junction Balancing).
   - Cubic Bezier curves for natural river visualization.
   - 3-panel subplots (5th, Median, 95th percentiles).
   - "Stretch" Gap Healing: Extends grid cells to nodes without fake interpolation.
   - Enhanced Annotations: Shows Value Range [Min-Max] per branch.

Usage: 
    python plot_netcdf_comprehensive.py <output_dir> [input_case_dir]
"""

import sys
import os
import glob
import numpy as np
import pandas as pd
import xarray as xr
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.patheffects as pe
import networkx as nx
from matplotlib.path import Path
from matplotlib.collections import PatchCollection

# --- Configuration for Maps ---
WIDTH_SCALE = 5.0     # Exaggerate width for visibility (1km -> 5 units)
Y_SPACING = 40.0      # Vertical spread for tree layout
CURVATURE = 0.35      # Bezier curve strength (reduced for hierarchical splits)

# --- ACADEMIC CONFIGURATION ---
# Set to True for "Seamless/Sankey" look (forces downstream sum = upstream width).
# Set to False for "Physical/Exact" look (shows actual widths from topology.csv).
BALANCE_JUNCTIONS = True 

# Geographic Y-hints for Mekong Delta branches (North = positive Y, South = negative Y)
# This ensures proper N→S ordering in the visualization
GEOGRAPHIC_Y_HINTS = {
    # Upstream inlets
    "Tan_Chau": 30,      # Tien inlet (slightly north)
    "Chau_Doc": -50,     # Hau inlet (south)
    
    # Vam Nao junctions
    "VamNao_Tien": 20,
    "VamNao_Hau": -30,
    
    # Tien splits
    "Co_Chien_Split": 10,
    "MyTho_Split": 15,
    
    # Outlets (N→S order)
    "MyTho_Mouth": 40,        # Northernmost
    "HamLuong_Mouth": 20,     # Central
    "CoChien_Mouth": 0,       # South Tien
    "Hau_Mouth": -60,         # Far South
} 

def ensure_dir(path):
    if not os.path.isdir(path):
        os.makedirs(path, exist_ok=True)

# =============================================================================
# PART 1: GRAPH LAYOUT & TOPOLOGY
# =============================================================================

def load_topology(output_dir, input_arg=None):
    candidates = []
    if input_arg: candidates.append(os.path.join(input_arg, 'topology.csv'))
    candidates.append(os.path.join(output_dir, 'topology.csv'))
    
    norm_out = os.path.normpath(output_dir)
    case_name = os.path.basename(norm_out)
    root_dir = os.path.dirname(os.path.dirname(norm_out))
    candidates.append(os.path.join(root_dir, 'INPUT', 'Cases', case_name, 'topology.csv'))
    candidates.append(os.path.join('INPUT', 'Cases', case_name, 'topology.csv'))
    
    for p in candidates:
        if os.path.exists(p):
            print(f"Loaded topology from: {p}")
            try:
                # Read with flexible column handling - topology may have 10-14 columns
                # Columns: ID, Name, NodeUp, NodeDown, Length, W_Up, W_Down, Depth, Chezy, 
                #          Group, [RS], [VDB_K], [Mixing_Alpha], [BiogeoParams]
                df = pd.read_csv(p, comment='#', skipinitialspace=True, header=None)
                
                # Assign column names based on how many columns are present
                base_cols = ['Branch_ID', 'Name', 'NodeUp', 'NodeDown', 
                             'Length', 'W_Up', 'W_Down', 'Depth', 'Chezy', 'Group']
                optional_cols = ['RS', 'VDB_K', 'Mixing_Alpha', 'BiogeoParams']
                
                ncols = len(df.columns)
                if ncols >= 10:
                    col_names = base_cols[:min(ncols, 10)]
                    if ncols > 10:
                        col_names.extend(optional_cols[:ncols-10])
                    df.columns = col_names[:ncols]
                else:
                    # Fallback for minimal topology
                    df.columns = base_cols[:ncols]
                
                # Strip whitespace from Name column
                if 'Name' in df.columns:
                    df['Name'] = df['Name'].str.strip()
                
                return df
            except Exception as e:
                print(f"Warning: Failed to parse {p}: {e}")
                continue
    print("Warning: topology.csv not found. Schematic maps will be skipped.")
    return None

def build_schematic_layout(df):
    """Auto-generates (x,y) coordinates for nodes using NetworkX.
    
    Uses GEOGRAPHIC_Y_HINTS for proper North-South ordering when available.
    Falls back to automatic tree layout for unknown nodes.
    """
    G = nx.DiGraph()
    
    # Build node name lookup from boundary_map if available (for Y-hint matching)
    node_names = {}
    
    for _, row in df.iterrows():
        G.add_edge(row['NodeUp'], row['NodeDown'], 
                   length=row['Length']/1000.0, data=row)

    roots = [n for n, d in G.in_degree() if d == 0]
    if not roots and len(G.nodes) > 0: roots = [list(G.nodes)[0]]
    
    # X-position: Cumulative distance from roots (upstream → downstream)
    pos_x = {}
    for root in roots: pos_x[root] = 0.0
    
    try:
        nodes_ordered = list(nx.topological_sort(G))
    except:
        nodes_ordered = list(G.nodes)

    for u in nodes_ordered:
        if u not in pos_x: pos_x[u] = 0.0
        for v in G.successors(u):
            dist = G[u][v]['length']
            pos_x[v] = max(pos_x.get(v, 0), pos_x[u] + dist)

    # Y-position: Use geographic hints if available, otherwise auto-layout
    pos_y = {}
    
    # Try to load boundary_map for node names
    try:
        # Look for boundary_map in same directory as topology
        import os
        topo_dir = os.path.dirname(df.attrs.get('source_path', ''))
        if not topo_dir:
            # Try default location
            topo_dir = 'INPUT/Cases/Mekong_Delta_Full'
        
        bmap_path = os.path.join(topo_dir, 'boundary_map.csv')
        if os.path.exists(bmap_path):
            bmap = pd.read_csv(bmap_path, comment='#', skipinitialspace=True, header=None)
            for _, row in bmap.iterrows():
                if len(row) >= 2:
                    node_id = int(row[0])
                    node_name = str(row[1]).strip()
                    node_names[node_id] = node_name
    except:
        pass
    
    # Apply geographic hints where available
    for node in G.nodes():
        node_name = node_names.get(node, '')
        if node_name in GEOGRAPHIC_Y_HINTS:
            pos_y[node] = GEOGRAPHIC_Y_HINTS[node_name]
    
    # For nodes without hints, use automatic tree layout
    def assign_y(node, y_center, height):
        if node in pos_y:
            return  # Already has geographic hint
        pos_y[node] = y_center
        children = list(G.successors(node))
        if not children: return
        
        # Filter children that don't have Y yet
        children_need_y = [c for c in children if c not in pos_y]
        if not children_need_y: return
        
        step = height / len(children_need_y)
        curr = y_center + (height/2.0) - (step/2.0)
        for child in sorted(children_need_y):
            assign_y(child, curr, step)
            curr -= step

    for i, root in enumerate(roots):
        if root not in pos_y:
            assign_y(root, i * -100, 100)
    
    # Ensure all nodes have Y position
    for node in G.nodes():
        if node not in pos_y:
            # Find nearest node with Y and offset
            for parent in G.predecessors(node):
                if parent in pos_y:
                    pos_y[node] = pos_y[parent] - 10
                    break
            if node not in pos_y:
                pos_y[node] = 0

    base_pos = {n: (pos_x[n], pos_y.get(n, 0)) for n in G.nodes}
    return G, base_pos

def calculate_edge_ports(G, pos):
    """
    Calculates start/end ports.
    If BALANCE_JUNCTIONS = True: Scales child widths to match parent.
    If BALANCE_JUNCTIONS = False: Uses exact physical widths (centered).
    """
    ports = {} 
    adjusted_widths = {} 

    # 1. Outgoing Splits
    for u in G.nodes():
        children = list(G.successors(u))
        if not children: continue
        
        parents = list(G.predecessors(u))
        
        # Get Physical Widths
        child_widths = []
        for v in children:
            row = G[u][v]['data']
            child_widths.append((row['W_Up'] / 1000.0) * WIDTH_SCALE)
            
        total_child_w = sum(child_widths)
        
        # --- BALANCE LOGIC ---
        scale = 1.0
        if BALANCE_JUNCTIONS and len(parents) == 1:
            p = parents[0]
            p_row = G[p][u]['data']
            w_parent_end = (p_row['W_Down'] / 1000.0) * WIDTH_SCALE
            if total_child_w > 0:
                scale = w_parent_end / total_child_w
        
        scaled_widths = [w * scale for w in child_widths]
        
        # Stack children
        children_sorted_indices = np.argsort([-pos[v][1] for v in children])
        total_stack_h = sum(scaled_widths)
        current_y = pos[u][1] + (total_stack_h / 2.0)
        
        for i in children_sorted_indices:
            v = children[i]
            w = scaled_widths[i]
            port_y = current_y - (w / 2.0)
            
            ports[(u,v)] = ports.get((u,v), {})
            ports[(u,v)]['start'] = (pos[u][0], port_y)
            
            # Store the width we decided to use (Physical or Scaled)
            adjusted_widths[(u,v,'start')] = w
            
            current_y -= w

    # 2. Incoming Merges
    for v in G.nodes():
        parents = list(G.predecessors(v))
        if not parents: continue
        parents.sort(key=lambda u: pos[u][1], reverse=True)
        
        widths = []
        for u in parents:
            row = G[u][v]['data']
            widths.append((row['W_Down'] / 1000.0) * WIDTH_SCALE)
            
        total_w = sum(widths)
        current_y = pos[v][1] + (total_w / 2.0)
        
        for i, u in enumerate(parents):
            w = widths[i]
            port_y = current_y - (w / 2.0)
            
            ports[(u,v)] = ports.get((u,v), {})
            ports[(u,v)]['end'] = (pos[v][0], port_y)
            adjusted_widths[(u,v,'end')] = w
            current_y -= w
            
    return ports, adjusted_widths

# =============================================================================
# PART 2: GEOMETRY GENERATION
# =============================================================================

def bezier_curve(p0, p1, p2, p3, n=50):
    t = np.linspace(0, 1, n).reshape(-1, 1)
    return (1-t)**3 * p0 + 3*(1-t)**2 * t * p1 + 3*(1-t)*t**2 * p2 + t**3 * p3

def get_normals(points):
    diffs = np.diff(points, axis=0)
    diffs = np.vstack([diffs, diffs[-1]])
    normals = np.zeros_like(diffs)
    normals[:, 0] = -diffs[:, 1]
    normals[:, 1] = diffs[:, 0]
    norms = np.linalg.norm(normals, axis=1, keepdims=True)
    with np.errstate(divide='ignore', invalid='ignore'):
        return np.where(norms > 0, normals / norms, 0)

def get_segmented_patches(start_xy, end_xy, w_start, w_end, row, x_grid, data_values):
    """Generates grid-exact polygons along the curve."""
    p0, p3 = np.array(start_xy), np.array(end_xy)
    
    dist_x = np.linalg.norm(p3 - p0)
    offset = dist_x * CURVATURE
    if abs(p3[0] - p0[0]) > abs(p3[1] - p0[1]):
        p1, p2 = p0 + [offset, 0], p3 - [offset, 0]
    else:
        p1, p2 = p0 + [0, -offset], p3 - [0, offset]
        
    fine_N = 300
    curve_fine = bezier_curve(p0, p1, p2, p3, fine_N)
    dists = np.linalg.norm(np.diff(curve_fine, axis=0), axis=1)
    cum_dist = np.insert(np.cumsum(dists), 0, 0)
    total_len = cum_dist[-1]
    
    L_phys = row['Length'] if row['Length'] > 1e-6 else 1.0
    
    # Map simulation x (0=Mouth -> L=Head) to visual t (0=Head -> 1=Mouth)
    t_grid = 1.0 - (x_grid / L_phys)
    
    # Sort to draw from Upstream to Downstream (t ascending)
    sort_idx = np.argsort(t_grid)
    t_grid = t_grid[sort_idx]
    data_values = data_values[sort_idx]
    
    # --- GAP HEALING (STRETCH METHOD) ---
    # Instead of inserting fake points, we force the valid grid bounds 
    # to touch 0.0 and 1.0. This stretches the first/last cells to fill the gap.
    if len(t_grid) > 0:
        t_grid[0] = 0.0
        t_grid[-1] = 1.0
    
    t_grid = np.clip(t_grid, 0, 1)
    
    # Map to curve
    normals_fine = get_normals(curve_fine)
    grid_points, grid_normals = [], []
    
    for t in t_grid:
        idx = np.searchsorted(cum_dist, t * total_len)
        idx = min(idx, fine_N - 1)
        grid_points.append(curve_fine[idx])
        grid_normals.append(normals_fine[idx])
        
    grid_points = np.array(grid_points)
    grid_normals = np.array(grid_normals)
    
    # Interpolate width
    if w_start < 1e-9: w_start = 0.1
    widths = w_start * (w_end/w_start)**t_grid
    
    patches_list = []
    colors = []
    
    for i in range(len(t_grid) - 1):
        val = data_values[i]
        if np.isnan(val): continue
        
        p1 = grid_points[i] + grid_normals[i] * widths[i]/2
        p2 = grid_points[i+1] + grid_normals[i+1] * widths[i+1]/2
        p3 = grid_points[i+1] - grid_normals[i+1] * widths[i+1]/2
        p4 = grid_points[i] - grid_normals[i] * widths[i]/2
        
        poly = patches.Polygon([p1, p2, p3, p4], closed=True)
        patches_list.append(poly)
        colors.append(val)
        
    return patches_list, colors

def plot_schematic_maps(ds, topo, out_root):
    map_dir = os.path.join(out_root, 'NETCDF_PLOTS', 'SCHEMATIC_MAPS')
    ensure_dir(map_dir)
    print(f"Generating Schematic Maps in: {map_dir}")
    
    G, pos = build_schematic_layout(topo)
    edge_ports, adj_widths = calculate_edge_ports(G, pos)
    
    vars_to_map = [v for v in ds.data_vars if 'branch' in ds[v].dims and 'x' in ds[v].dims]
    quantiles = [0.05, 0.5, 0.95]
    q_labels = ['5th Percentile', 'Median', '95th Percentile']
    
    for var_name in vars_to_map:
        if var_name in ['branch_x', 'x', 'time']: continue
        print(f"  Processing {var_name}...")
        
        try:
            data_q = ds[var_name].quantile(quantiles, dim='time', skipna=True)
        except: continue

        valid = data_q.values.flatten()
        valid = valid[~np.isnan(valid)]
        if len(valid) == 0: continue
        vmin, vmax = valid.min(), valid.max()
        
        fig, axes = plt.subplots(3, 1, figsize=(16, 24))
        cmap_name = 'viridis'
        if 'sal' in var_name.lower(): cmap_name = 'YlGnBu'
        elif 'sed' in var_name.lower() or 'spm' in var_name.lower(): cmap_name = 'YlOrBr'
        elif 'vel' in var_name.lower(): cmap_name = 'RdBu_r'
        
        for i, ax in enumerate(axes):
            q_val = data_q.isel(quantile=i)
            
            all_patches, all_colors = [], []
            
            for u, v, data in G.edges(data=True):
                row = data['data']
                try:
                    b_idx = list(ds.branch.values).index(row['Name'])
                except: continue
                
                x_vals = ds['branch_x'][b_idx, :].values
                val_vals = q_val[b_idx, :].values
                mask = ~np.isnan(x_vals)
                x_vals, val_vals = x_vals[mask], val_vals[mask]
                if len(val_vals) < 1: continue
                
                start_pt = edge_ports.get((u,v), {}).get('start', pos[u])
                end_pt = edge_ports.get((u,v), {}).get('end', pos[v])
                
                w_start = adj_widths.get((u,v,'start'), (row['W_Up']/1000)*WIDTH_SCALE)
                w_end = adj_widths.get((u,v,'end'), (row['W_Down']/1000)*WIDTH_SCALE)
                
                segs, cols = get_segmented_patches(start_pt, end_pt, w_start, w_end, row, x_vals, val_vals)
                all_patches.extend(segs)
                all_colors.extend(cols)
                
                # Stats Annotation
                b_min = np.min(val_vals)
                b_max = np.max(val_vals)
                
                mid_x = (start_pt[0] + end_pt[0])/2
                mid_y = (start_pt[1] + end_pt[1])/2
                label_txt = f"{row['Name']}\n[{b_min:.1f} - {b_max:.1f}]"
                
                ax.text(mid_x, mid_y, label_txt, fontsize=9, ha='center', va='center', 
                        fontweight='bold', color='white',
                        path_effects=[pe.withStroke(linewidth=2, foreground="black")])

            if all_patches:
                coll = PatchCollection(all_patches, cmap=cmap_name, alpha=1.0, 
                                     edgecolor='none', linewidth=0) 
                coll.set_array(np.array(all_colors))
                coll.set_clim(vmin, vmax)
                ax.add_collection(coll)
            
            for n, (nx, ny) in pos.items():
                ax.plot(nx, ny, 'ko', markersize=2)
            
            ax.autoscale()
            ax.set_aspect('equal')
            ax.axis('off')
            ax.set_title(f"{var_name}: {q_labels[i]}", fontsize=16, weight='bold')

        cbar_ax = fig.add_axes([0.92, 0.15, 0.02, 0.7])
        sm = plt.cm.ScalarMappable(cmap=cmap_name, norm=plt.Normalize(vmin=vmin, vmax=vmax))
        sm.set_array([])
        fig.colorbar(sm, cax=cbar_ax, label=var_name)
        
        out_name = f"Schematic_{var_name}_Stats.png"
        plt.savefig(os.path.join(map_dir, out_name), dpi=150, bbox_inches='tight')
        plt.close()
        print(f"  > Saved {out_name}")

# =============================================================================
# PART 3: DETAILED PLOTTING (Branch details)
# =============================================================================

def plot_branch_detailed(ds, branch_name, out_root):
    out_dir = os.path.join(out_root, 'NETCDF_PLOTS', branch_name)
    ensure_dir(out_dir)

    if 'branch' in ds.dims:
        try:
            b_idx = list(ds.branch.values).index(branch_name)
            subset = ds.isel(branch=b_idx)
        except: return
    else:
        subset = ds

    x = subset['branch_x'].values if 'branch_x' in subset else subset['x'].values
    mask = ~np.isnan(x)
    if not np.any(mask): return
    x_km = x[mask] / 1000.0
    
    # Transform for display: Distance from Upstream (assuming L is max x)
    L = np.max(x_km)
    x_display = L - x_km # 0 at upstream
    
    time = subset.time.values
    t_hours = np.arange(len(time))
    if np.issubdtype(time.dtype, np.datetime64):
        t_hours = (time - time[0]) / np.timedelta64(1, 'h')
    
    for var_name in subset.data_vars:
        if var_name in ['branch_x', 'x', 'time']: continue
        val = subset[var_name].values
        if val.ndim != 2: continue
        val = val[:, mask]
        
        # Hovmoller
        plt.figure(figsize=(10, 4))
        cmap = 'YlGnBu' if 'sal' in var_name.lower() else 'viridis'
        plt.pcolormesh(x_display, t_hours, val, shading='auto', cmap=cmap)
        plt.xlabel('Distance from Upstream (km)')
        plt.ylabel('Time (hours)')
        plt.colorbar(label=var_name)
        plt.title(f'{branch_name}: {var_name}')
        plt.tight_layout()
        plt.savefig(os.path.join(out_dir, f'hovmoller_{var_name}.png'))
        plt.close()
        
        # Time Series
        locs = [0, len(x_km)//2, len(x_km)-1]
        plt.figure(figsize=(10, 4))
        for i in locs:
            label_dist = x_display[i]
            plt.plot(t_hours, val[:, i], label=f'{label_dist:.1f} km')
        plt.legend(title="Dist. from Upstream")
        plt.xlabel('Time (hours)')
        plt.ylabel(var_name)
        plt.title(f'{branch_name}: {var_name} Time Series')
        plt.tight_layout()
        plt.savefig(os.path.join(out_dir, f'ts_{var_name}.png'))
        plt.close()

# =============================================================================
# MAIN
# =============================================================================

def main():
    if len(sys.argv) < 2:
        print("Usage: python plot_netcdf_comprehensive.py <output_dir> [input_case_dir]")
        return 1
        
    out_dir = sys.argv[1]
    in_arg = sys.argv[2] if len(sys.argv) > 2 else None
    
    if not os.path.isdir(out_dir):
        potential = os.path.join("OUTPUT", out_dir)
        if os.path.isdir(potential): out_dir = potential
        else:
            print(f"Error: {out_dir} not found.")
            return 1

    # 1. Find NetCDF
    nc_path = os.path.join(out_dir, "combined_branches.nc")
    if not os.path.exists(nc_path):
        nc_files = glob.glob(os.path.join(out_dir, "*.nc"))
        if not nc_files:
            print("No .nc files found.")
            return 1
        nc_path = nc_files[0]
    
    print(f"Loading data: {nc_path}")
    ds = xr.open_dataset(nc_path)
    
    # 2. Load Topology
    topo = load_topology(out_dir, in_arg)
    
    # 3. Generate Schematic Maps
    if topo is not None and 'branch' in ds.dims:
        plot_schematic_maps(ds, topo, out_dir)
        
    # 4. Generate Branch Plots
    if 'branch' in ds.dims:
        for b_name in ds.branch.values:
            print(f"Plotting details: {b_name}")
            plot_branch_detailed(ds, b_name, out_dir)
    else:
        bname = os.path.splitext(os.path.basename(nc_path))[0]
        print(f"Plotting details: {bname}")
        plot_branch_detailed(ds, bname, out_dir)
        
    print("\nDone.")

if __name__ == "__main__":
    main()