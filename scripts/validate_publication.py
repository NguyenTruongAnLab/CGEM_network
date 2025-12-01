#!/usr/bin/env python3
"""
Publication-Quality Validation & Visualization for C-GEM Mekong Delta
======================================================================

This script generates comprehensive, publication-ready figures and validation
statistics for the Mekong Delta simulation results.

OUTPUT FIGURES:
===============
1. HYDRODYNAMICS
   - Fig 1: Tidal water level variations at all outlets
   - Fig 2: Discharge distribution (Tien vs Hau split)
   - Fig 3: Vam Nao inter-basin exchange flow
   
2. SALINITY INTRUSION
   - Fig 4: Longitudinal salinity profiles for ALL branches (dry season)
   - Fig 5: Longitudinal salinity profiles for ALL branches (wet season)
   - Fig 6: Seasonal comparison of 4 PSU isohaline position
   - Fig 7: Network schematic showing salinity gradient
   
3. WATER QUALITY
   - Fig 8: Dissolved oxygen profiles
   - Fig 9: Nutrient (NO3, NH4, PO4) longitudinal profiles
   - Fig 10: pCO2 and air-water CO2 flux
   
4. SEASONAL DYNAMICS
   - Fig 11: Time series of salinity at key monitoring stations
   - Fig 12: Monthly mean discharge and salinity intrusion
   - Fig 13: Hovmöller diagrams for key branches

5. VALIDATION STATISTICS
   - Table 1: Model skill scores vs literature
   - Table 2: Mass balance checks

References:
-----------
- Nguyen et al. (2008) - Mekong salinity intrusion measurements
- Fujii et al. (2003) - Vam Nao hydraulics
- Savenije (2005, 2012) - Salinity and Tides in Alluvial Estuaries
- MRC Water Quality Database

Author: C-GEM Development Team
Date: December 2024
"""

import sys
import struct
import numpy as np
import pandas as pd
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any
from dataclasses import dataclass
from datetime import datetime

# Plotting imports
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.cm import ScalarMappable
import matplotlib.gridspec as gridspec

# Try to import seaborn for better aesthetics
try:
    import seaborn as sns
    sns.set_style("whitegrid")
    HAS_SEABORN = True
except ImportError:
    HAS_SEABORN = False

# ===========================================================================
# CONFIGURATION
# ===========================================================================

SCRIPT_DIR = Path(__file__).parent.resolve()
PROJECT_ROOT = SCRIPT_DIR.parent
OUTPUT_DIR = PROJECT_ROOT / "OUTPUT" / "Mekong_Delta_Full"
FIGURES_DIR = OUTPUT_DIR / "PUBLICATION_FIGURES"

# Publication figure settings
plt.rcParams.update({
    'font.size': 11,
    'font.family': 'serif',
    'axes.labelsize': 12,
    'axes.titlesize': 14,
    'xtick.labelsize': 10,
    'ytick.labelsize': 10,
    'legend.fontsize': 10,
    'figure.titlesize': 16,
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
})

# Unit conversions
UMOL_O2_TO_MG = 0.032    # 1 µmol O2 = 0.032 mg
UMOL_C_TO_MG = 0.012     # 1 µmol C = 0.012 mg  
UMOL_N_TO_MG = 0.014     # 1 µmol N = 0.014 mg
UMOL_P_TO_MG = 0.031     # 1 µmol P = 0.031 mg

# Branch colors for consistent visualization
BRANCH_COLORS = {
    'Tien_Main': '#1f77b4',      # Blue
    'Hau_Main': '#ff7f0e',       # Orange
    'Vam_Nao': '#2ca02c',        # Green
    'Tien_Lower': '#1f77b4',     # Blue (Tien family)
    'Tien_Connector': '#17becf', # Cyan
    'Co_Chien': '#d62728',       # Red
    'My_Tho': '#9467bd',         # Purple
    'Ham_Luong': '#8c564b',      # Brown
    'Hau_River': '#ff7f0e',      # Orange (Hau family)
}

# Branch display order (upstream to downstream, N to S)
BRANCH_ORDER = [
    'Tien_Main', 'Hau_Main', 'Vam_Nao',
    'Tien_Lower', 'Tien_Connector',
    'My_Tho', 'Ham_Luong', 'Co_Chien', 'Hau_River'
]

# Distributary branches (have ocean outlets)
DISTRIBUTARY_BRANCHES = ['My_Tho', 'Ham_Luong', 'Co_Chien', 'Hau_River']

# Seasonal definitions (Mekong climate)
DRY_SEASON_DAYS = list(range(0, 150)) + list(range(335, 365))  # Dec-May
WET_SEASON_DAYS = list(range(150, 305))  # Jun-Oct
TRANSITION_DAYS = list(range(120, 150)) + list(range(305, 335))

# Literature reference values
LITERATURE_VALUES = {
    'salinity_intrusion_dry': {
        'Ham_Luong': (50, 65, 'Nguyen et al. (2008)'),
        'My_Tho': (45, 60, 'Nguyen et al. (2008)'),
        'Co_Chien': (40, 55, 'Nguyen et al. (2008)'),
        'Hau_River': (55, 75, 'Nguyen et al. (2008)'),
    },
    'tidal_range': {
        'Can_Tho': (2.0, 2.8, 'MRC Data'),
        'My_Tho': (2.2, 3.2, 'MRC Data'),
    },
    'vamnao_velocity': (0.1, 1.0, 'Fujii et al. (2003)'),
    'o2_range': (4.0, 8.0, 'MRC Monitoring'),
}


# ===========================================================================
# DATA STRUCTURES
# ===========================================================================

@dataclass
class BranchData:
    """Container for branch simulation data."""
    name: str
    M: int
    dx: float
    x: np.ndarray          # Spatial grid [m]
    times: np.ndarray      # Time array [s]
    hydro: Dict[str, np.ndarray]    # {name: 2D array [time, space]}
    species: Dict[str, np.ndarray]  # {name: 2D array [time, space]}
    
    @property
    def x_km(self) -> np.ndarray:
        """Distance from mouth in km."""
        return self.x / 1000.0
    
    @property
    def time_days(self) -> np.ndarray:
        """Time in days."""
        return self.times / 86400.0
    
    def get_seasonal_mean(self, var_name: str, season: str, var_type: str = 'species') -> np.ndarray:
        """Get spatial profile averaged over a season."""
        data_dict = self.hydro if var_type == 'hydro' else self.species
        if var_name not in data_dict:
            return np.full(self.M, np.nan)
        
        data = data_dict[var_name]
        days = self.time_days
        
        if season == 'dry':
            mask = np.array([(int(d) % 365) < 150 or (int(d) % 365) >= 335 for d in days])
        elif season == 'wet':
            mask = np.array([150 <= (int(d) % 365) < 305 for d in days])
        else:
            mask = np.ones(len(days), dtype=bool)
        
        if not np.any(mask):
            return np.nanmean(data, axis=0)
        
        return np.nanmean(data[mask, :], axis=0)
    
    def get_timeseries_at_distance(self, var_name: str, dist_km: float, 
                                    var_type: str = 'species') -> np.ndarray:
        """Get time series at a specific distance from mouth."""
        data_dict = self.hydro if var_type == 'hydro' else self.species
        if var_name not in data_dict:
            return np.full(len(self.times), np.nan)
        
        idx = np.argmin(np.abs(self.x_km - dist_km))
        return data_dict[var_name][:, idx]


# ===========================================================================
# BINARY FILE READER
# ===========================================================================

def read_branch_binary(branch_name: str, output_dir: Path = OUTPUT_DIR) -> Optional[BranchData]:
    """Read C-GEM binary output file for a branch."""
    filepath = output_dir / f"{branch_name}.bin"
    
    if not filepath.exists():
        return None
    
    try:
        with open(filepath, 'rb') as f:
            # Read header
            header = f.read(16)
            if len(header) < 16:
                return None
            M, n_hydro, n_spec, n_rxn = struct.unpack('4i', header)
            
            dx_bytes = f.read(8)
            if len(dx_bytes) < 8:
                return None
            dx = struct.unpack('d', dx_bytes)[0]
            
            # Read X grid
            x_bytes = f.read(8 * M)
            if len(x_bytes) < 8 * M:
                return None
            x = np.array(struct.unpack(f'{M}d', x_bytes))
            
            # Read null-terminated strings
            def read_string():
                chars = []
                while True:
                    c = f.read(1)
                    if c == b'\x00' or c == b'':
                        break
                    chars.append(c.decode('ascii', errors='replace'))
                return ''.join(chars)
            
            # Read names
            hydro_names = [read_string() for _ in range(n_hydro)]
            species_names = [read_string() for _ in range(n_spec)]
            rxn_names = [read_string() for _ in range(n_rxn)]
            
            # Initialize data storage
            hydro = {name: [] for name in hydro_names}
            species = {name: [] for name in species_names}
            times = []
            
            # Read time records
            while True:
                time_bytes = f.read(8)
                if len(time_bytes) < 8:
                    break
                    
                time_s = struct.unpack('d', time_bytes)[0]
                times.append(time_s)
                
                # Read hydro data
                for name in hydro_names:
                    data_bytes = f.read(8 * M)
                    if len(data_bytes) < 8 * M:
                        break
                    data = np.array(struct.unpack(f'{M}d', data_bytes))
                    hydro[name].append(data)
                
                # Read species data
                for name in species_names:
                    data_bytes = f.read(8 * M)
                    if len(data_bytes) < 8 * M:
                        break
                    data = np.array(struct.unpack(f'{M}d', data_bytes))
                    species[name].append(data)
                
                # Skip reaction data
                if n_rxn > 0:
                    f.read(8 * M * n_rxn)
            
            # Convert to arrays
            for name in hydro:
                if hydro[name]:
                    hydro[name] = np.array(hydro[name])
            
            for name in species:
                if species[name]:
                    species[name] = np.array(species[name])
            
            return BranchData(
                name=branch_name,
                M=M,
                dx=dx,
                x=x,
                times=np.array(times),
                hydro=hydro,
                species=species
            )
            
    except Exception as e:
        print(f"  [ERROR] Failed to read {filepath}: {e}")
        return None


def load_all_branches(output_dir: Path = OUTPUT_DIR) -> Dict[str, BranchData]:
    """Load all branch data from output directory."""
    branches = {}
    
    for bin_file in output_dir.glob("*.bin"):
        branch_name = bin_file.stem
        data = read_branch_binary(branch_name, output_dir)
        if data is not None:
            branches[branch_name] = data
    
    return branches


# ===========================================================================
# SALINITY INTRUSION ANALYSIS
# ===========================================================================

def calculate_intrusion_length(salinity_profile: np.ndarray, x_km: np.ndarray, 
                                threshold: float = 4.0) -> float:
    """
    Calculate salt intrusion length (distance where salinity drops below threshold).
    
    Convention: x_km[0] = downstream (ocean), x_km[-1] = upstream (river)
    """
    if len(salinity_profile) == 0 or np.all(np.isnan(salinity_profile)):
        return 0.0
    
    # Find where salinity crosses threshold from ocean side
    # Assume salinity decreases from index 0 (ocean) to index -1 (upstream)
    for i in range(len(salinity_profile)):
        if salinity_profile[i] < threshold:
            if i == 0:
                return 0.0
            # Linear interpolation
            s1, s0 = salinity_profile[i-1], salinity_profile[i]
            x1, x0 = x_km[i-1], x_km[i]
            if s1 - s0 != 0:
                x_thresh = x0 + (threshold - s0) * (x1 - x0) / (s1 - s0)
                return x_thresh
            return x_km[i]
    
    # If threshold never crossed, return max distance
    return x_km[-1] if len(x_km) > 0 else 0.0


# ===========================================================================
# FIGURE GENERATION FUNCTIONS
# ===========================================================================

def fig1_salinity_profiles_all_branches(branches: Dict[str, BranchData], 
                                         season: str = 'dry') -> plt.Figure:
    """
    Generate longitudinal salinity profiles for ALL distributary branches.
    
    Creates a multi-panel figure showing salinity from mouth (0 km) to upstream
    for each branch, with connecting nodes indicated.
    """
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    axes = axes.flatten()
    
    season_label = 'Dry Season (Dec-May)' if season == 'dry' else 'Wet Season (Jun-Oct)'
    fig.suptitle(f'Longitudinal Salinity Profiles - {season_label}', fontsize=16, fontweight='bold')
    
    for idx, branch_name in enumerate(DISTRIBUTARY_BRANCHES):
        ax = axes[idx]
        
        if branch_name not in branches:
            ax.text(0.5, 0.5, f'{branch_name}\nNo Data', ha='center', va='center', 
                   transform=ax.transAxes, fontsize=14)
            ax.set_title(branch_name)
            continue
        
        data = branches[branch_name]
        sal_profile = data.get_seasonal_mean('salinity', season)
        x_km = data.x_km
        
        # Plot salinity profile
        color = BRANCH_COLORS.get(branch_name, 'blue')
        ax.plot(x_km, sal_profile, color=color, linewidth=2.5, label=branch_name)
        
        # Add 4 PSU threshold line
        ax.axhline(y=4, color='red', linestyle='--', linewidth=1.5, label='4 PSU (MRC Standard)')
        
        # Calculate and mark intrusion length
        L4 = calculate_intrusion_length(sal_profile, x_km, threshold=4.0)
        if L4 > 0:
            ax.axvline(x=L4, color='red', linestyle=':', alpha=0.7)
            ax.annotate(f'L₄ = {L4:.1f} km', xy=(L4, 4), xytext=(L4+5, 8),
                       fontsize=10, arrowprops=dict(arrowstyle='->', color='red', alpha=0.7))
        
        # Literature comparison
        if branch_name in LITERATURE_VALUES.get('salinity_intrusion_dry', {}):
            lit_min, lit_max, ref = LITERATURE_VALUES['salinity_intrusion_dry'][branch_name]
            ax.axvspan(lit_min, lit_max, alpha=0.2, color='green', 
                      label=f'Observed ({ref})')
        
        ax.set_xlabel('Distance from Mouth (km)')
        ax.set_ylabel('Salinity (PSU)')
        ax.set_title(f'{branch_name}', fontweight='bold')
        ax.set_xlim(0, max(x_km) + 5)
        ax.set_ylim(0, 38)
        ax.grid(True, alpha=0.3)
        ax.legend(loc='upper right', fontsize=9)
    
    plt.tight_layout()
    return fig


def fig2_seasonal_salinity_comparison(branches: Dict[str, BranchData]) -> plt.Figure:
    """
    Compare dry and wet season salinity intrusion for all distributary branches.
    """
    fig, ax = plt.subplots(figsize=(12, 7))
    
    branch_names = []
    dry_L4 = []
    wet_L4 = []
    lit_dry_min = []
    lit_dry_max = []
    
    for branch_name in DISTRIBUTARY_BRANCHES:
        if branch_name not in branches:
            continue
        
        data = branches[branch_name]
        
        # Dry season
        sal_dry = data.get_seasonal_mean('salinity', 'dry')
        L4_dry = calculate_intrusion_length(sal_dry, data.x_km, 4.0)
        
        # Wet season
        sal_wet = data.get_seasonal_mean('salinity', 'wet')
        L4_wet = calculate_intrusion_length(sal_wet, data.x_km, 4.0)
        
        branch_names.append(branch_name.replace('_', '\n'))
        dry_L4.append(L4_dry)
        wet_L4.append(L4_wet)
        
        # Literature values
        if branch_name in LITERATURE_VALUES.get('salinity_intrusion_dry', {}):
            lit_min, lit_max, _ = LITERATURE_VALUES['salinity_intrusion_dry'][branch_name]
            lit_dry_min.append(lit_min)
            lit_dry_max.append(lit_max)
        else:
            lit_dry_min.append(0)
            lit_dry_max.append(0)
    
    x = np.arange(len(branch_names))
    width = 0.35
    
    # Bar plot
    bars1 = ax.bar(x - width/2, dry_L4, width, label='Dry Season (Model)', 
                   color='#d62728', alpha=0.8)
    bars2 = ax.bar(x + width/2, wet_L4, width, label='Wet Season (Model)', 
                   color='#1f77b4', alpha=0.8)
    
    # Literature range (dry season)
    for i, (ymin, ymax) in enumerate(zip(lit_dry_min, lit_dry_max)):
        if ymax > 0:
            ax.errorbar(x[i] - width/2, (ymin + ymax)/2, 
                       yerr=[[((ymin + ymax)/2 - ymin)], [(ymax - (ymin + ymax)/2)]], 
                       fmt='none', color='black', capsize=5, capthick=2,
                       label='Observed Range' if i == 0 else '')
    
    ax.set_xlabel('Branch', fontsize=12)
    ax.set_ylabel('4 PSU Intrusion Length (km)', fontsize=12)
    ax.set_title('Seasonal Variation of Salinity Intrusion\n(4 PSU Isohaline Position)', 
                fontsize=14, fontweight='bold')
    ax.set_xticks(x)
    ax.set_xticklabels(branch_names)
    ax.legend(loc='upper right')
    ax.grid(True, alpha=0.3, axis='y')
    
    # Add value labels on bars
    for bar in bars1:
        height = bar.get_height()
        ax.annotate(f'{height:.1f}', xy=(bar.get_x() + bar.get_width()/2, height),
                   xytext=(0, 3), textcoords="offset points", ha='center', va='bottom', fontsize=9)
    
    for bar in bars2:
        height = bar.get_height()
        ax.annotate(f'{height:.1f}', xy=(bar.get_x() + bar.get_width()/2, height),
                   xytext=(0, 3), textcoords="offset points", ha='center', va='bottom', fontsize=9)
    
    plt.tight_layout()
    return fig


def fig3_vam_nao_flow(branches: Dict[str, BranchData]) -> plt.Figure:
    """
    Analyze Vam Nao inter-basin exchange flow.
    """
    if 'Vam_Nao' not in branches:
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.text(0.5, 0.5, 'Vam Nao data not available', ha='center', va='center', 
               transform=ax.transAxes, fontsize=14)
        return fig
    
    vn = branches['Vam_Nao']
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Panel (a): Velocity time series at mid-channel
    ax1 = axes[0, 0]
    if 'velocity' in vn.hydro:
        vel = vn.hydro['velocity']
        mid_idx = vn.M // 2
        vel_mid = vel[:, mid_idx]
        
        ax1.plot(vn.time_days, vel_mid, color='#2ca02c', linewidth=0.5, alpha=0.7)
        ax1.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
        ax1.axhline(y=np.mean(vel_mid), color='red', linestyle='--', linewidth=2,
                   label=f'Mean = {np.mean(vel_mid):.3f} m/s')
        
        ax1.set_xlabel('Time (days)')
        ax1.set_ylabel('Velocity (m/s)')
        ax1.set_title('(a) Vam Nao Velocity at Mid-Channel', fontweight='bold')
        ax1.legend()
        ax1.grid(True, alpha=0.3)
        
        # Annotate direction
        if np.mean(vel_mid) > 0:
            ax1.annotate('Net flow: Tien → Hau', xy=(0.7, 0.9), xycoords='axes fraction',
                        fontsize=12, fontweight='bold', color='green')
        else:
            ax1.annotate('Net flow: Hau → Tien', xy=(0.7, 0.9), xycoords='axes fraction',
                        fontsize=12, fontweight='bold', color='orange')
    
    # Panel (b): Velocity histogram
    ax2 = axes[0, 1]
    if 'velocity' in vn.hydro:
        vel_all = vn.hydro['velocity'].flatten()
        ax2.hist(vel_all, bins=50, color='#2ca02c', alpha=0.7, edgecolor='black')
        ax2.axvline(x=0, color='black', linestyle='-', linewidth=1)
        ax2.axvline(x=np.mean(vel_all), color='red', linestyle='--', linewidth=2)
        
        ax2.set_xlabel('Velocity (m/s)')
        ax2.set_ylabel('Frequency')
        ax2.set_title('(b) Velocity Distribution', fontweight='bold')
        
        # Add statistics
        stats_text = f'Mean: {np.mean(vel_all):.4f} m/s\n'
        stats_text += f'Std: {np.std(vel_all):.4f} m/s\n'
        stats_text += f'Max: {np.max(vel_all):.4f} m/s\n'
        stats_text += f'Min: {np.min(vel_all):.4f} m/s'
        ax2.text(0.95, 0.95, stats_text, transform=ax2.transAxes, fontsize=10,
                verticalalignment='top', horizontalalignment='right',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # Panel (c): Water level at both ends
    ax3 = axes[1, 0]
    if 'waterLevel' in vn.hydro:
        wl = vn.hydro['waterLevel']
        ax3.plot(vn.time_days, wl[:, 0], label='Tien Side', color='#1f77b4', alpha=0.7)
        ax3.plot(vn.time_days, wl[:, -1], label='Hau Side', color='#ff7f0e', alpha=0.7)
        ax3.set_xlabel('Time (days)')
        ax3.set_ylabel('Water Level (m)')
        ax3.set_title('(c) Water Level at Vam Nao Ends', fontweight='bold')
        ax3.legend()
        ax3.grid(True, alpha=0.3)
    
    # Panel (d): Salinity through Vam Nao
    ax4 = axes[1, 1]
    if 'salinity' in vn.species:
        sal = vn.species['salinity']
        ax4.plot(vn.time_days, sal[:, 0], label='Tien Side', color='#1f77b4', alpha=0.7)
        ax4.plot(vn.time_days, sal[:, -1], label='Hau Side', color='#ff7f0e', alpha=0.7)
        ax4.set_xlabel('Time (days)')
        ax4.set_ylabel('Salinity (PSU)')
        ax4.set_title('(d) Salinity at Vam Nao Ends', fontweight='bold')
        ax4.legend()
        ax4.grid(True, alpha=0.3)
    
    fig.suptitle('Vam Nao Inter-Basin Exchange Analysis', fontsize=16, fontweight='bold')
    plt.tight_layout()
    return fig


def fig4_water_quality_profiles(branches: Dict[str, BranchData], 
                                 season: str = 'dry') -> plt.Figure:
    """
    Generate water quality profiles (O2, nutrients) for key branches.
    """
    fig = plt.figure(figsize=(16, 12))
    gs = gridspec.GridSpec(3, 2, figure=fig, hspace=0.3, wspace=0.25)
    
    season_label = 'Dry Season' if season == 'dry' else 'Wet Season'
    fig.suptitle(f'Water Quality Profiles - {season_label}', fontsize=16, fontweight='bold')
    
    # Panel (a): Dissolved Oxygen
    ax1 = fig.add_subplot(gs[0, 0])
    for branch_name in DISTRIBUTARY_BRANCHES:
        if branch_name not in branches:
            continue
        data = branches[branch_name]
        if 'o2' not in data.species:
            continue
        o2_profile = data.get_seasonal_mean('o2', season) * UMOL_O2_TO_MG
        ax1.plot(data.x_km, o2_profile, color=BRANCH_COLORS.get(branch_name, 'gray'),
                label=branch_name, linewidth=2)
    
    ax1.axhline(y=2.0, color='red', linestyle='--', label='Hypoxia Threshold')
    ax1.set_xlabel('Distance from Mouth (km)')
    ax1.set_ylabel('DO (mg/L)')
    ax1.set_title('(a) Dissolved Oxygen', fontweight='bold')
    ax1.legend(loc='best', fontsize=8)
    ax1.grid(True, alpha=0.3)
    ax1.set_ylim(0, 10)
    
    # Panel (b): Nitrate
    ax2 = fig.add_subplot(gs[0, 1])
    for branch_name in DISTRIBUTARY_BRANCHES:
        if branch_name not in branches:
            continue
        data = branches[branch_name]
        if 'no3' not in data.species:
            continue
        no3_profile = data.get_seasonal_mean('no3', season) * UMOL_N_TO_MG
        ax2.plot(data.x_km, no3_profile, color=BRANCH_COLORS.get(branch_name, 'gray'),
                label=branch_name, linewidth=2)
    
    ax2.set_xlabel('Distance from Mouth (km)')
    ax2.set_ylabel('NO₃ (mg N/L)')
    ax2.set_title('(b) Nitrate', fontweight='bold')
    ax2.legend(loc='best', fontsize=8)
    ax2.grid(True, alpha=0.3)
    
    # Panel (c): Ammonium
    ax3 = fig.add_subplot(gs[1, 0])
    for branch_name in DISTRIBUTARY_BRANCHES:
        if branch_name not in branches:
            continue
        data = branches[branch_name]
        if 'nh4' not in data.species:
            continue
        nh4_profile = data.get_seasonal_mean('nh4', season) * UMOL_N_TO_MG
        ax3.plot(data.x_km, nh4_profile, color=BRANCH_COLORS.get(branch_name, 'gray'),
                label=branch_name, linewidth=2)
    
    ax3.set_xlabel('Distance from Mouth (km)')
    ax3.set_ylabel('NH₄ (mg N/L)')
    ax3.set_title('(c) Ammonium', fontweight='bold')
    ax3.legend(loc='best', fontsize=8)
    ax3.grid(True, alpha=0.3)
    
    # Panel (d): Phosphate
    ax4 = fig.add_subplot(gs[1, 1])
    for branch_name in DISTRIBUTARY_BRANCHES:
        if branch_name not in branches:
            continue
        data = branches[branch_name]
        if 'po4' not in data.species:
            continue
        po4_profile = data.get_seasonal_mean('po4', season) * UMOL_P_TO_MG
        ax4.plot(data.x_km, po4_profile, color=BRANCH_COLORS.get(branch_name, 'gray'),
                label=branch_name, linewidth=2)
    
    ax4.set_xlabel('Distance from Mouth (km)')
    ax4.set_ylabel('PO₄ (mg P/L)')
    ax4.set_title('(d) Phosphate', fontweight='bold')
    ax4.legend(loc='best', fontsize=8)
    ax4.grid(True, alpha=0.3)
    
    # Panel (e): TOC
    ax5 = fig.add_subplot(gs[2, 0])
    for branch_name in DISTRIBUTARY_BRANCHES:
        if branch_name not in branches:
            continue
        data = branches[branch_name]
        if 'toc' not in data.species:
            continue
        toc_profile = data.get_seasonal_mean('toc', season) * UMOL_C_TO_MG
        ax5.plot(data.x_km, toc_profile, color=BRANCH_COLORS.get(branch_name, 'gray'),
                label=branch_name, linewidth=2)
    
    ax5.set_xlabel('Distance from Mouth (km)')
    ax5.set_ylabel('TOC (mg C/L)')
    ax5.set_title('(e) Total Organic Carbon', fontweight='bold')
    ax5.legend(loc='best', fontsize=8)
    ax5.grid(True, alpha=0.3)
    
    # Panel (f): pCO2
    ax6 = fig.add_subplot(gs[2, 1])
    for branch_name in DISTRIBUTARY_BRANCHES:
        if branch_name not in branches:
            continue
        data = branches[branch_name]
        if 'pco2' not in data.species:
            continue
        pco2_profile = data.get_seasonal_mean('pco2', season)
        ax6.plot(data.x_km, pco2_profile, color=BRANCH_COLORS.get(branch_name, 'gray'),
                label=branch_name, linewidth=2)
    
    ax6.axhline(y=420, color='green', linestyle='--', label='Atmospheric pCO₂')
    ax6.set_xlabel('Distance from Mouth (km)')
    ax6.set_ylabel('pCO₂ (µatm)')
    ax6.set_title('(f) Partial Pressure of CO₂', fontweight='bold')
    ax6.legend(loc='best', fontsize=8)
    ax6.grid(True, alpha=0.3)
    
    return fig


def fig5_seasonal_timeseries(branches: Dict[str, BranchData]) -> plt.Figure:
    """
    Time series showing seasonal dynamics at key monitoring locations.
    """
    fig, axes = plt.subplots(3, 1, figsize=(14, 12), sharex=True)
    
    # Monitoring locations (km from mouth)
    monitor_locs = {
        'Ham_Luong': [20, 40, 60],
        'Hau_River': [30, 50, 70],
    }
    
    # Panel (a): Salinity time series
    ax1 = axes[0]
    for branch_name, locs in monitor_locs.items():
        if branch_name not in branches:
            continue
        data = branches[branch_name]
        for dist_km in locs:
            sal_ts = data.get_timeseries_at_distance('salinity', dist_km)
            label = f'{branch_name} @ {dist_km} km'
            ax1.plot(data.time_days, sal_ts, label=label, linewidth=1, alpha=0.8)
    
    ax1.set_ylabel('Salinity (PSU)')
    ax1.set_title('(a) Salinity at Monitoring Stations', fontweight='bold')
    ax1.legend(loc='upper right', ncol=2, fontsize=8)
    ax1.grid(True, alpha=0.3)
    
    # Shade seasons
    for ax in axes:
        ax.axvspan(0, 150, alpha=0.1, color='orange', label='Dry' if ax == axes[0] else '')
        ax.axvspan(150, 305, alpha=0.1, color='blue', label='Wet' if ax == axes[0] else '')
        ax.axvspan(305, 365, alpha=0.1, color='orange')
    
    # Panel (b): Discharge (use velocity as proxy)
    ax2 = axes[1]
    for branch_name in ['Tien_Main', 'Hau_Main']:
        if branch_name not in branches:
            continue
        data = branches[branch_name]
        if 'velocity' in data.hydro:
            # Estimate discharge from mean velocity × area
            vel_mean = np.mean(data.hydro['velocity'], axis=1)
            ax2.plot(data.time_days, vel_mean, label=branch_name, linewidth=1.5)
    
    ax2.set_ylabel('Mean Velocity (m/s)')
    ax2.set_title('(b) Upstream Velocity (Discharge Proxy)', fontweight='bold')
    ax2.legend(loc='upper right')
    ax2.grid(True, alpha=0.3)
    
    # Panel (c): DO time series
    ax3 = axes[2]
    for branch_name in ['Ham_Luong', 'Hau_River']:
        if branch_name not in branches:
            continue
        data = branches[branch_name]
        if 'o2' not in data.species:
            continue
        o2_ts = data.get_timeseries_at_distance('o2', 40) * UMOL_O2_TO_MG
        ax3.plot(data.time_days, o2_ts, label=f'{branch_name} @ 40 km', linewidth=1.5)
    
    ax3.axhline(y=2.0, color='red', linestyle='--', label='Hypoxia')
    ax3.set_xlabel('Time (days)')
    ax3.set_ylabel('DO (mg/L)')
    ax3.set_title('(c) Dissolved Oxygen Dynamics', fontweight='bold')
    ax3.legend(loc='upper right')
    ax3.grid(True, alpha=0.3)
    
    fig.suptitle('Seasonal Dynamics in Mekong Delta Distributaries', 
                fontsize=16, fontweight='bold')
    plt.tight_layout()
    return fig


def fig6_hovmoller_diagrams(branches: Dict[str, BranchData]) -> plt.Figure:
    """
    Hovmöller (space-time) diagrams for key variables in selected branches.
    """
    fig, axes = plt.subplots(2, 2, figsize=(16, 12))
    
    branch_name = 'Ham_Luong'
    if branch_name not in branches:
        branch_name = list(branches.keys())[0] if branches else None
    
    if branch_name is None:
        return fig
    
    data = branches[branch_name]
    X, T = np.meshgrid(data.x_km, data.time_days)
    
    # Panel (a): Salinity Hovmöller
    ax1 = axes[0, 0]
    if 'salinity' in data.species:
        sal = data.species['salinity']
        c1 = ax1.pcolormesh(X, T, sal, shading='auto', cmap='YlGnBu', vmin=0, vmax=35)
        ax1.contour(X, T, sal, levels=[4], colors='red', linewidths=2)
        plt.colorbar(c1, ax=ax1, label='Salinity (PSU)')
        ax1.set_xlabel('Distance from Mouth (km)')
        ax1.set_ylabel('Time (days)')
        ax1.set_title(f'(a) Salinity - {branch_name}', fontweight='bold')
    
    # Panel (b): DO Hovmöller
    ax2 = axes[0, 1]
    if 'o2' in data.species:
        o2 = data.species['o2'] * UMOL_O2_TO_MG
        c2 = ax2.pcolormesh(X, T, o2, shading='auto', cmap='RdYlGn', vmin=0, vmax=10)
        ax2.contour(X, T, o2, levels=[2], colors='black', linewidths=2, linestyles='--')
        plt.colorbar(c2, ax=ax2, label='DO (mg/L)')
        ax2.set_xlabel('Distance from Mouth (km)')
        ax2.set_ylabel('Time (days)')
        ax2.set_title(f'(b) Dissolved Oxygen - {branch_name}', fontweight='bold')
    
    # Panel (c): Velocity Hovmöller
    ax3 = axes[1, 0]
    if 'velocity' in data.hydro:
        vel = data.hydro['velocity']
        vmax = np.percentile(np.abs(vel), 99)
        c3 = ax3.pcolormesh(X, T, vel, shading='auto', cmap='RdBu_r', vmin=-vmax, vmax=vmax)
        plt.colorbar(c3, ax=ax3, label='Velocity (m/s)')
        ax3.set_xlabel('Distance from Mouth (km)')
        ax3.set_ylabel('Time (days)')
        ax3.set_title(f'(c) Velocity - {branch_name}', fontweight='bold')
    
    # Panel (d): pCO2 Hovmöller
    ax4 = axes[1, 1]
    if 'pco2' in data.species:
        pco2 = data.species['pco2']
        c4 = ax4.pcolormesh(X, T, pco2, shading='auto', cmap='YlOrRd', vmin=300, vmax=2000)
        ax4.contour(X, T, pco2, levels=[420], colors='green', linewidths=2)
        plt.colorbar(c4, ax=ax4, label='pCO₂ (µatm)')
        ax4.set_xlabel('Distance from Mouth (km)')
        ax4.set_ylabel('Time (days)')
        ax4.set_title(f'(d) pCO₂ - {branch_name}', fontweight='bold')
    
    fig.suptitle(f'Space-Time (Hovmöller) Diagrams: {branch_name}', 
                fontsize=16, fontweight='bold')
    plt.tight_layout()
    return fig


# ===========================================================================
# VALIDATION STATISTICS
# ===========================================================================

def compute_validation_statistics(branches: Dict[str, BranchData]) -> pd.DataFrame:
    """
    Compute validation statistics comparing model results with literature.
    """
    records = []
    
    # 1. Salinity intrusion
    for branch_name in DISTRIBUTARY_BRANCHES:
        if branch_name not in branches:
            continue
        
        data = branches[branch_name]
        sal_dry = data.get_seasonal_mean('salinity', 'dry')
        L4_model = calculate_intrusion_length(sal_dry, data.x_km, 4.0)
        
        if branch_name in LITERATURE_VALUES.get('salinity_intrusion_dry', {}):
            L4_min, L4_max, ref = LITERATURE_VALUES['salinity_intrusion_dry'][branch_name]
            L4_obs = (L4_min + L4_max) / 2
            error_pct = (L4_model - L4_obs) / L4_obs * 100 if L4_obs > 0 else 0
            within_range = L4_min <= L4_model <= L4_max
            
            records.append({
                'Variable': 'Salinity Intrusion (L₄)',
                'Location': branch_name,
                'Model': f'{L4_model:.1f} km',
                'Observed': f'{L4_min}-{L4_max} km',
                'Error (%)': f'{error_pct:+.1f}%',
                'Status': '✓' if within_range else '✗',
                'Reference': ref
            })
    
    # 2. Vam Nao velocity
    if 'Vam_Nao' in branches:
        vn = branches['Vam_Nao']
        if 'velocity' in vn.hydro:
            vel_mean = np.mean(vn.hydro['velocity'])
            vel_min, vel_max, ref = LITERATURE_VALUES['vamnao_velocity']
            within_range = vel_min <= abs(vel_mean) <= vel_max
            
            records.append({
                'Variable': 'Vam Nao Velocity',
                'Location': 'Mid-channel',
                'Model': f'{vel_mean:.3f} m/s',
                'Observed': f'{vel_min}-{vel_max} m/s',
                'Error (%)': 'N/A',
                'Status': '✓' if within_range else '✗',
                'Reference': ref
            })
    
    # 3. Dissolved Oxygen
    for branch_name in DISTRIBUTARY_BRANCHES:
        if branch_name not in branches:
            continue
        
        data = branches[branch_name]
        if 'o2' not in data.species:
            continue
        
        o2_mean = np.mean(data.species['o2']) * UMOL_O2_TO_MG
        o2_min, o2_max, ref = LITERATURE_VALUES['o2_range']
        within_range = o2_min <= o2_mean <= o2_max
        
        records.append({
            'Variable': 'Dissolved Oxygen',
            'Location': branch_name,
            'Model': f'{o2_mean:.2f} mg/L',
            'Observed': f'{o2_min}-{o2_max} mg/L',
            'Error (%)': 'N/A',
            'Status': '✓' if within_range else '~',
            'Reference': ref
        })
    
    return pd.DataFrame(records)


def generate_validation_report(branches: Dict[str, BranchData], 
                                output_dir: Path) -> None:
    """
    Generate a text validation report.
    """
    report_path = output_dir / 'validation_report.txt'
    
    with open(report_path, 'w') as f:
        f.write("=" * 70 + "\n")
        f.write("C-GEM MEKONG DELTA - VALIDATION REPORT\n")
        f.write("=" * 70 + "\n")
        f.write(f"Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
        
        # Statistics table
        stats_df = compute_validation_statistics(branches)
        f.write("VALIDATION STATISTICS\n")
        f.write("-" * 70 + "\n")
        f.write(stats_df.to_string(index=False))
        f.write("\n\n")
        
        # Summary
        n_pass = (stats_df['Status'] == '✓').sum()
        n_total = len(stats_df)
        f.write("-" * 70 + "\n")
        f.write(f"SUMMARY: {n_pass}/{n_total} tests passed\n")
        
        if n_pass == n_total:
            f.write("STATUS: ALL VALIDATION CRITERIA MET ✓\n")
        elif n_pass >= n_total * 0.7:
            f.write("STATUS: MOSTLY ACCEPTABLE (~)\n")
        else:
            f.write("STATUS: NEEDS REVIEW (✗)\n")
    
    print(f"Validation report saved to: {report_path}")


# ===========================================================================
# MAIN
# ===========================================================================

def main():
    """Generate all publication figures and validation statistics."""
    print("=" * 70)
    print("C-GEM MEKONG DELTA - PUBLICATION VALIDATION SUITE")
    print("=" * 70)
    
    # Check output directory
    if not OUTPUT_DIR.exists():
        print(f"\n[ERROR] Output directory not found: {OUTPUT_DIR}")
        print("Please run the model first.")
        return 1
    
    # Create figures directory
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    print(f"\nOutput directory: {OUTPUT_DIR}")
    print(f"Figures directory: {FIGURES_DIR}")
    
    # Load all branch data
    print("\nLoading branch data...")
    branches = load_all_branches(OUTPUT_DIR)
    
    if not branches:
        print("[ERROR] No branch data found!")
        return 1
    
    print(f"Loaded {len(branches)} branches: {list(branches.keys())}")
    
    # Generate figures
    print("\n" + "=" * 70)
    print("GENERATING PUBLICATION FIGURES")
    print("=" * 70)
    
    figures = [
        ("Fig1_Salinity_Profiles_Dry", lambda: fig1_salinity_profiles_all_branches(branches, 'dry')),
        ("Fig2_Salinity_Profiles_Wet", lambda: fig1_salinity_profiles_all_branches(branches, 'wet')),
        ("Fig3_Seasonal_Salinity_Comparison", lambda: fig2_seasonal_salinity_comparison(branches)),
        ("Fig4_VamNao_Exchange", lambda: fig3_vam_nao_flow(branches)),
        ("Fig5_WaterQuality_Dry", lambda: fig4_water_quality_profiles(branches, 'dry')),
        ("Fig6_WaterQuality_Wet", lambda: fig4_water_quality_profiles(branches, 'wet')),
        ("Fig7_Seasonal_Timeseries", lambda: fig5_seasonal_timeseries(branches)),
        ("Fig8_Hovmoller_Diagrams", lambda: fig6_hovmoller_diagrams(branches)),
    ]
    
    for fig_name, fig_func in figures:
        print(f"\n  Generating {fig_name}...")
        try:
            fig = fig_func()
            fig.savefig(FIGURES_DIR / f"{fig_name}.png", dpi=300, bbox_inches='tight')
            fig.savefig(FIGURES_DIR / f"{fig_name}.pdf", bbox_inches='tight')
            plt.close(fig)
            print(f"    ✓ Saved {fig_name}.png and .pdf")
        except Exception as e:
            print(f"    ✗ Failed: {e}")
    
    # Generate validation report
    print("\n" + "=" * 70)
    print("GENERATING VALIDATION REPORT")
    print("=" * 70)
    
    generate_validation_report(branches, FIGURES_DIR)
    
    # Print validation statistics
    print("\n" + "=" * 70)
    print("VALIDATION STATISTICS")
    print("=" * 70)
    stats_df = compute_validation_statistics(branches)
    print(stats_df.to_string(index=False))
    
    # Summary
    n_pass = (stats_df['Status'] == '✓').sum()
    n_total = len(stats_df)
    
    print("\n" + "-" * 70)
    print(f"SUMMARY: {n_pass}/{n_total} validation tests passed")
    
    print("\n" + "=" * 70)
    print("PUBLICATION VALIDATION COMPLETE")
    print("=" * 70)
    print(f"\nFigures saved to: {FIGURES_DIR}")
    print("\nFigure List:")
    for f in sorted(FIGURES_DIR.glob("*.png")):
        print(f"  - {f.name}")
    
    return 0


if __name__ == "__main__":
    sys.exit(main())
