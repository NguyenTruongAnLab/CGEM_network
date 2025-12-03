#!/usr/bin/env python3
"""
C-GEM Network Publication Figures Generator
============================================

Creates publication-quality figures showing:
1. Network topology with node connections
2. Salinity profiles and intrusion lengths
3. SPM distributions and ETM locations  
4. Water quality (O2, nutrients, phytoplankton)
5. Greenhouse gas emissions (pCO2, CH4, N2O)
6. Calibration results and validation

Author: C-GEM Team
Date: December 2025

Usage:
    python scripts/plot_publication_figures.py [--case CASE_NAME] [--output OUTPUT_DIR]
"""

import os
import sys
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Circle
from matplotlib.lines import Line2D
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.gridspec as gridspec
from pathlib import Path

# Use publication-quality settings
plt.rcParams.update({
    'font.family': 'serif',
    'font.size': 10,
    'axes.labelsize': 11,
    'axes.titlesize': 12,
    'legend.fontsize': 9,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'axes.linewidth': 0.8,
    'lines.linewidth': 1.5,
    'axes.grid': True,
    'grid.alpha': 0.3,
})

# Color scheme for branches (colorblind-friendly)
BRANCH_COLORS = {
    'Tien_Main': '#1f77b4',       # Blue
    'Hau_Main': '#ff7f0e',        # Orange  
    'Vam_Nao': '#2ca02c',         # Green
    'Tien_Lower': '#1f77b4',      # Blue (same as Tien_Main)
    'Tien_Connector': '#9467bd',  # Purple
    'Co_Chien': '#d62728',        # Red
    'My_Tho': '#8c564b',          # Brown
    'Ham_Luong': '#e377c2',       # Pink
    'Hau_River': '#ff7f0e',       # Orange (same as Hau_Main)
}

# Branch groups for plotting
RIVER_BRANCHES = ['Tien_Main', 'Hau_Main', 'Tien_Lower', 'Tien_Connector', 'Vam_Nao']
ESTUARY_BRANCHES = ['My_Tho', 'Ham_Luong', 'Co_Chien', 'Hau_River']


class MekongNetworkPlotter:
    """Publication-quality plotter for C-GEM Mekong Delta network."""
    
    def __init__(self, case_dir, output_dir=None):
        """
        Initialize plotter with case directory.
        
        Args:
            case_dir: Path to case directory (e.g., INPUT/Cases/Mekong_Delta_Full)
            output_dir: Path for output figures (default: OUTPUT/{case_name}/figures)
        """
        self.case_dir = Path(case_dir)
        self.case_name = self.case_dir.name
        
        if output_dir:
            self.output_dir = Path(output_dir)
        else:
            self.output_dir = Path(f'OUTPUT/{self.case_name}/figures')
        
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Load data
        self.csv_dir = Path(f'OUTPUT/{self.case_name}/CSV')
        self.calib_dir = Path(f'OUTPUT/{self.case_name}/calibration')
        
        self.topology = self._load_topology()
        self.branches = list(self.topology['Name'].values) if self.topology is not None else []
        
    def _load_topology(self):
        """Load network topology from CSV."""
        topo_file = self.case_dir / 'topology.csv'
        if not topo_file.exists():
            print(f"Warning: Topology file not found: {topo_file}")
            return None
        
        # Read CSV with comment handling
        lines = []
        with open(topo_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith('#'):
                    lines.append(line)
        
        if len(lines) < 1:
            return None
        
        # The topology file doesn't have a header - columns are:
        # ID, Name, NodeUp, NodeDown, Length, W_Up, W_Down, Depth, Chezy, Group, RS, VDB_K, Mixing_Alpha
        from io import StringIO
        data = '\n'.join(lines)
        df = pd.read_csv(StringIO(data), skipinitialspace=True, header=None,
                        names=['ID', 'Name', 'NodeUp', 'NodeDown', 'Length', 'W_Up', 'W_Down', 
                               'Depth', 'Chezy', 'Group', 'RS', 'VDB_K', 'Mixing_Alpha'])
        
        # Clean up branch names (strip whitespace)
        df['Name'] = df['Name'].str.strip()
        return df
    
    def _load_branch_data(self, branch_name, variable):
        """Load time-series data for a variable from a branch."""
        filename = f'{branch_name}_{variable}.csv'
        filepath = self.csv_dir / filename
        
        if not filepath.exists():
            return None
            
        df = pd.read_csv(filepath)
        return df
    
    def _load_calibration_results(self):
        """Load calibration results if available."""
        results_file = self.calib_dir / 'calibration_results.csv'
        if not results_file.exists():
            return None
        return pd.read_csv(results_file, comment='#')
        
    def plot_network_topology(self, figsize=(14, 10), save=True):
        """
        Create a schematic network topology figure showing all branches and nodes.
        
        This is a conceptual diagram, not geographically accurate.
        """
        fig, ax = plt.subplots(figsize=figsize)
        ax.set_xlim(-2, 18)
        ax.set_ylim(-2, 14)
        ax.set_aspect('equal')
        ax.axis('off')
        
        # Node positions (conceptual layout)
        # Upstream nodes at top, downstream at bottom
        node_pos = {
            1: (3, 12),    # Tien Inlet (discharge)
            2: (11, 12),   # Hau Inlet (discharge)
            3: (4, 9),     # Junction (Tien_Main/Vam_Nao/Tien_Lower)
            4: (10, 9),    # Junction (Hau_Main/Vam_Nao/Hau_River)
            5: (3, 6),     # Junction (Tien_Lower/Connector/Co_Chien)
            10: (5, 3),    # Junction (Connector/My_Tho/Ham_Luong)
            6: (2, 0),     # My_Tho outlet (tidal)
            7: (6, 0),     # Ham_Luong outlet (tidal)
            8: (1, 3),     # Co_Chien outlet (tidal)
            9: (13, 3),    # Hau_River outlet (tidal)
        }
        
        # Branch connections (from topology)
        branch_connections = {
            'Tien_Main': (1, 3),
            'Hau_Main': (2, 4),
            'Vam_Nao': (3, 4),
            'Tien_Lower': (3, 5),
            'Tien_Connector': (5, 10),
            'Co_Chien': (5, 8),
            'My_Tho': (10, 6),
            'Ham_Luong': (10, 7),
            'Hau_River': (4, 9),
        }
        
        # Draw branches
        for branch, (n1, n2) in branch_connections.items():
            if n1 in node_pos and n2 in node_pos:
                x1, y1 = node_pos[n1]
                x2, y2 = node_pos[n2]
                
                color = BRANCH_COLORS.get(branch, '#888888')
                
                # Draw branch as thick line
                ax.plot([x1, x2], [y1, y2], 
                       color=color, linewidth=8, solid_capstyle='round',
                       alpha=0.7, zorder=1)
                
                # Add branch label at midpoint
                mx, my = (x1 + x2) / 2, (y1 + y2) / 2
                
                # Offset label to avoid overlap
                dx, dy = x2 - x1, y2 - y1
                length = np.sqrt(dx*dx + dy*dy)
                if length > 0:
                    # Perpendicular offset
                    offset = 0.4
                    ox, oy = -dy/length * offset, dx/length * offset
                else:
                    ox, oy = 0, 0
                
                ax.text(mx + ox, my + oy, branch.replace('_', ' '), 
                       fontsize=8, ha='center', va='center',
                       bbox=dict(boxstyle='round,pad=0.2', facecolor='white', 
                                alpha=0.8, edgecolor='none'))
        
        # Draw nodes
        for node_id, (x, y) in node_pos.items():
            if node_id in [1, 2]:
                # Discharge inlet nodes (blue)
                circle = Circle((x, y), 0.5, color='#2166ac', zorder=3)
                ax.add_patch(circle)
                ax.text(x, y+0.9, f'Q{node_id}\nDischarge', ha='center', va='bottom', fontsize=8)
            elif node_id in [6, 7, 8, 9]:
                # Tidal outlet nodes (orange)
                circle = Circle((x, y), 0.5, color='#d95f0e', zorder=3)
                ax.add_patch(circle)
                outlet_names = {6: 'My Tho', 7: 'Ham Luong', 8: 'Co Chien', 9: 'Hau'}
                ax.text(x, y-1.0, f'{outlet_names.get(node_id, "")}\n(Tidal)', 
                       ha='center', va='top', fontsize=8)
            else:
                # Junction nodes (gray)
                circle = Circle((x, y), 0.35, color='#666666', zorder=3)
                ax.add_patch(circle)
                ax.text(x+0.5, y+0.3, f'J{node_id}', ha='left', va='bottom', fontsize=7)
        
        # Add legend
        legend_elements = [
            Circle((0, 0), 0.1, color='#2166ac', label='Discharge Inlet'),
            Circle((0, 0), 0.1, color='#d95f0e', label='Tidal Outlet'),
            Circle((0, 0), 0.1, color='#666666', label='Junction'),
        ]
        
        # Create proxy artists for legend
        from matplotlib.patches import Patch
        legend_handles = [
            Patch(facecolor='#2166ac', label='Discharge Inlet'),
            Patch(facecolor='#d95f0e', label='Tidal Outlet'),
            Patch(facecolor='#666666', label='Junction'),
        ]
        ax.legend(handles=legend_handles, loc='upper right', framealpha=0.9)
        
        # Title
        ax.set_title('Mekong Delta Network Topology\n(C-GEM Hierarchical Representation)', 
                    fontsize=14, fontweight='bold')
        
        # Add flow direction arrows
        ax.annotate('', xy=(7, 1), xytext=(7, 11),
                   arrowprops=dict(arrowstyle='->', color='gray', lw=2, alpha=0.5))
        ax.text(7.3, 6, 'Flow\nDirection', fontsize=8, color='gray', ha='left', va='center')
        
        # Add scale info
        ax.text(0, -1.5, 'Note: Schematic diagram - not to scale', 
               fontsize=8, style='italic', color='gray')
        
        plt.tight_layout()
        
        if save:
            outpath = self.output_dir / 'network_topology.png'
            plt.savefig(outpath, dpi=300, bbox_inches='tight', facecolor='white')
            print(f"Saved: {outpath}")
        
        return fig
    
    def plot_salinity_profiles(self, time_idx=-1, figsize=(12, 8), save=True):
        """
        Plot salinity profiles for all estuarine branches.
        
        Args:
            time_idx: Time index to plot (-1 for final state)
            figsize: Figure size
            save: Whether to save the figure
        """
        fig, axes = plt.subplots(2, 2, figsize=figsize, sharex=False, sharey=True)
        axes = axes.flatten()
        
        branches = ['Hau_River', 'Ham_Luong', 'My_Tho', 'Co_Chien']
        
        for idx, branch in enumerate(branches):
            ax = axes[idx]
            
            df = self._load_branch_data(branch, 'salinity')
            if df is None:
                ax.text(0.5, 0.5, f'No data for {branch}', transform=ax.transAxes,
                       ha='center', va='center')
                continue
            
            # Get distance columns
            dist_cols = [c for c in df.columns if c.endswith('km')]
            distances = [float(c.replace('km', '')) for c in dist_cols]
            
            # Get salinity profile at specified time
            if time_idx < 0:
                time_idx = len(df) + time_idx
            
            salinity = df.iloc[time_idx][dist_cols].values
            
            # Plot profile
            color = BRANCH_COLORS.get(branch, '#333333')
            ax.plot(distances, salinity, color=color, linewidth=2, label=branch.replace('_', ' '))
            
            # Add 4 PSU threshold line
            ax.axhline(y=4.0, color='red', linestyle='--', alpha=0.7, label='4 PSU threshold')
            
            # Find intrusion length (where salinity crosses 4 PSU)
            for i in range(len(salinity)-1):
                if salinity[i] >= 4.0 and salinity[i+1] < 4.0:
                    # Linear interpolation
                    frac = (salinity[i] - 4.0) / (salinity[i] - salinity[i+1] + 1e-10)
                    intrusion_km = distances[i] + frac * (distances[i+1] - distances[i])
                    ax.axvline(x=intrusion_km, color='green', linestyle=':', alpha=0.7)
                    ax.text(intrusion_km + 2, 20, f'L = {intrusion_km:.0f} km', 
                           fontsize=9, color='green')
                    break
            
            ax.set_xlabel('Distance from mouth (km)')
            ax.set_ylabel('Salinity (PSU)')
            ax.set_title(branch.replace('_', ' '))
            ax.set_xlim(0, max(distances))
            ax.set_ylim(0, 35)
            ax.grid(True, alpha=0.3)
            ax.legend(loc='upper right', fontsize=8)
        
        # Get time info
        try:
            time_s = df.iloc[time_idx]['Time_s']
            time_days = time_s / 86400.0
            time_str = f'Day {time_days:.1f}'
        except:
            time_str = ''
        
        fig.suptitle(f'Salinity Intrusion Profiles - Mekong Delta\n{time_str}', 
                    fontsize=14, fontweight='bold')
        plt.tight_layout()
        
        if save:
            outpath = self.output_dir / 'salinity_profiles.png'
            plt.savefig(outpath, dpi=300, bbox_inches='tight', facecolor='white')
            print(f"Saved: {outpath}")
        
        return fig
    
    def plot_spm_profiles(self, time_idx=-1, figsize=(12, 8), save=True):
        """Plot SPM profiles and ETM locations."""
        fig, axes = plt.subplots(2, 2, figsize=figsize, sharex=False, sharey=False)
        axes = axes.flatten()
        
        branches = ['Hau_River', 'Ham_Luong', 'My_Tho', 'Co_Chien']
        
        for idx, branch in enumerate(branches):
            ax = axes[idx]
            
            df = self._load_branch_data(branch, 'spm')
            if df is None:
                ax.text(0.5, 0.5, f'No data for {branch}', transform=ax.transAxes,
                       ha='center', va='center')
                continue
            
            dist_cols = [c for c in df.columns if c.endswith('km')]
            distances = [float(c.replace('km', '')) for c in dist_cols]
            
            if time_idx < 0:
                time_idx = len(df) + time_idx
            
            spm = df.iloc[time_idx][dist_cols].values
            
            color = BRANCH_COLORS.get(branch, '#333333')
            ax.fill_between(distances, 0, spm, color=color, alpha=0.3)
            ax.plot(distances, spm, color=color, linewidth=2, label=branch.replace('_', ' '))
            
            # Find ETM location (max SPM)
            max_idx = np.argmax(spm)
            etm_km = distances[max_idx]
            max_spm = spm[max_idx]
            
            ax.scatter([etm_km], [max_spm], color='red', s=100, marker='*', 
                      zorder=5, label=f'ETM: {etm_km:.0f} km')
            
            ax.set_xlabel('Distance from mouth (km)')
            ax.set_ylabel('SPM (mg/L)')
            ax.set_title(branch.replace('_', ' '))
            ax.set_xlim(0, max(distances))
            ax.grid(True, alpha=0.3)
            ax.legend(loc='upper right', fontsize=8)
        
        fig.suptitle('Suspended Particulate Matter (SPM) Profiles\nEstuarine Turbidity Maximum (ETM)', 
                    fontsize=14, fontweight='bold')
        plt.tight_layout()
        
        if save:
            outpath = self.output_dir / 'spm_profiles.png'
            plt.savefig(outpath, dpi=300, bbox_inches='tight', facecolor='white')
            print(f"Saved: {outpath}")
        
        return fig
    
    def plot_water_quality_summary(self, time_idx=-1, figsize=(14, 10), save=True):
        """Plot summary of water quality variables (O2, NO3, NH4, PHY1)."""
        fig, axes = plt.subplots(2, 2, figsize=figsize)
        
        variables = [
            ('o2', 'Dissolved Oxygen', 'µM', 0, 300),
            ('no3', 'Nitrate (NO₃)', 'µM', 0, 100),
            ('nh4', 'Ammonium (NH₄)', 'µM', 0, 20),
            ('phy1', 'Phytoplankton (PHY1)', 'µM C', 0, 50)
        ]
        
        branches = ['Hau_River', 'Ham_Luong', 'My_Tho', 'Co_Chien']
        
        for ax_idx, (var, title, unit, vmin, vmax) in enumerate(variables):
            ax = axes.flatten()[ax_idx]
            
            for branch in branches:
                df = self._load_branch_data(branch, var)
                if df is None:
                    continue
                
                dist_cols = [c for c in df.columns if c.endswith('km')]
                distances = [float(c.replace('km', '')) for c in dist_cols]
                
                if time_idx < 0:
                    time_idx_use = len(df) + time_idx
                else:
                    time_idx_use = time_idx
                
                values = df.iloc[time_idx_use][dist_cols].values
                
                color = BRANCH_COLORS.get(branch, '#333333')
                ax.plot(distances, values, color=color, linewidth=1.5, 
                       label=branch.replace('_', ' '), alpha=0.8)
            
            ax.set_xlabel('Distance from mouth (km)')
            ax.set_ylabel(f'{title} ({unit})')
            ax.set_title(title)
            ax.set_ylim(vmin, vmax)
            ax.grid(True, alpha=0.3)
            ax.legend(loc='best', fontsize=8)
        
        fig.suptitle('Water Quality Profiles - Mekong Delta Estuary', 
                    fontsize=14, fontweight='bold')
        plt.tight_layout()
        
        if save:
            outpath = self.output_dir / 'water_quality_summary.png'
            plt.savefig(outpath, dpi=300, bbox_inches='tight', facecolor='white')
            print(f"Saved: {outpath}")
        
        return fig
    
    def plot_ghg_emissions(self, time_idx=-1, figsize=(14, 8), save=True):
        """Plot greenhouse gas concentrations (pCO2, CH4, N2O)."""
        fig, axes = plt.subplots(1, 3, figsize=figsize)
        
        variables = [
            ('pco2', 'pCO₂', 'µatm', 200, 2000),
            ('ch4', 'Methane (CH₄)', 'nM', 0, 500),
            ('n2o', 'Nitrous Oxide (N₂O)', 'nM', 0, 50)
        ]
        
        branches = ['Hau_River', 'Ham_Luong', 'My_Tho', 'Co_Chien']
        
        for ax_idx, (var, title, unit, vmin, vmax) in enumerate(variables):
            ax = axes[ax_idx]
            
            for branch in branches:
                df = self._load_branch_data(branch, var)
                if df is None:
                    continue
                
                dist_cols = [c for c in df.columns if c.endswith('km')]
                distances = [float(c.replace('km', '')) for c in dist_cols]
                
                if time_idx < 0:
                    time_idx_use = len(df) + time_idx
                else:
                    time_idx_use = time_idx
                
                values = df.iloc[time_idx_use][dist_cols].values
                
                color = BRANCH_COLORS.get(branch, '#333333')
                ax.plot(distances, values, color=color, linewidth=1.5, 
                       label=branch.replace('_', ' '), alpha=0.8)
            
            # Add atmospheric equilibrium line for pCO2
            if var == 'pco2':
                ax.axhline(y=420, color='green', linestyle='--', alpha=0.7, 
                          label='Atmospheric (420 µatm)')
            
            ax.set_xlabel('Distance from mouth (km)')
            ax.set_ylabel(f'{title} ({unit})')
            ax.set_title(title)
            ax.set_ylim(vmin, vmax)
            ax.grid(True, alpha=0.3)
            ax.legend(loc='best', fontsize=8)
        
        fig.suptitle('Greenhouse Gas Concentrations - Mekong Delta Estuary', 
                    fontsize=14, fontweight='bold')
        plt.tight_layout()
        
        if save:
            outpath = self.output_dir / 'ghg_emissions.png'
            plt.savefig(outpath, dpi=300, bbox_inches='tight', facecolor='white')
            print(f"Saved: {outpath}")
        
        return fig
    
    def plot_combined_network_figure(self, time_idx=-1, figsize=(16, 20), save=True):
        """
        Create a comprehensive combined figure with all network branches
        showing spatial distributions along the network.
        """
        fig = plt.figure(figsize=figsize)
        gs = gridspec.GridSpec(5, 2, height_ratios=[1.5, 1, 1, 1, 1], hspace=0.3, wspace=0.25)
        
        # Panel A: Network schematic (spans both columns)
        ax_network = fig.add_subplot(gs[0, :])
        self._draw_network_schematic(ax_network)
        ax_network.set_title('A) Mekong Delta Network Topology', fontsize=12, fontweight='bold', loc='left')
        
        # Panel B: Salinity profiles
        ax_sal = fig.add_subplot(gs[1, 0])
        self._plot_variable_all_branches(ax_sal, 'salinity', 'Salinity (PSU)', time_idx, 
                                         estuary_only=True, ylim=(0, 35))
        ax_sal.axhline(y=4.0, color='red', linestyle='--', alpha=0.7, label='4 PSU')
        ax_sal.set_title('B) Salinity Intrusion', fontsize=12, fontweight='bold', loc='left')
        
        # Panel C: SPM profiles
        ax_spm = fig.add_subplot(gs[1, 1])
        self._plot_variable_all_branches(ax_spm, 'spm', 'SPM (mg/L)', time_idx, 
                                         estuary_only=True, ylim=(0, 200))
        ax_spm.set_title('C) Suspended Particulate Matter', fontsize=12, fontweight='bold', loc='left')
        
        # Panel D: Dissolved oxygen
        ax_o2 = fig.add_subplot(gs[2, 0])
        self._plot_variable_all_branches(ax_o2, 'o2', 'O₂ (µM)', time_idx, 
                                         estuary_only=True, ylim=(0, 300))
        ax_o2.axhline(y=200, color='green', linestyle='--', alpha=0.7, label='Saturation')
        ax_o2.set_title('D) Dissolved Oxygen', fontsize=12, fontweight='bold', loc='left')
        
        # Panel E: Nitrate
        ax_no3 = fig.add_subplot(gs[2, 1])
        self._plot_variable_all_branches(ax_no3, 'no3', 'NO₃ (µM)', time_idx, 
                                         estuary_only=True, ylim=(0, 60))
        ax_no3.set_title('E) Nitrate', fontsize=12, fontweight='bold', loc='left')
        
        # Panel F: pCO2
        ax_pco2 = fig.add_subplot(gs[3, 0])
        self._plot_variable_all_branches(ax_pco2, 'pco2', 'pCO₂ (µatm)', time_idx, 
                                         estuary_only=True, ylim=(200, 2000))
        ax_pco2.axhline(y=420, color='green', linestyle='--', alpha=0.7, label='Atmospheric')
        ax_pco2.set_title('F) Partial Pressure CO₂', fontsize=12, fontweight='bold', loc='left')
        
        # Panel G: Phytoplankton
        ax_phy = fig.add_subplot(gs[3, 1])
        self._plot_variable_all_branches(ax_phy, 'phy1', 'PHY1 (µM C)', time_idx, 
                                         estuary_only=True, ylim=(0, 30))
        ax_phy.set_title('G) Diatom Phytoplankton', fontsize=12, fontweight='bold', loc='left')
        
        # Panel H: Methane
        ax_ch4 = fig.add_subplot(gs[4, 0])
        self._plot_variable_all_branches(ax_ch4, 'ch4', 'CH₄ (nM)', time_idx, 
                                         estuary_only=True, ylim=(0, 500))
        ax_ch4.set_title('H) Methane', fontsize=12, fontweight='bold', loc='left')
        
        # Panel I: N2O
        ax_n2o = fig.add_subplot(gs[4, 1])
        self._plot_variable_all_branches(ax_n2o, 'n2o', 'N₂O (nM)', time_idx, 
                                         estuary_only=True, ylim=(0, 30))
        ax_n2o.set_title('I) Nitrous Oxide', fontsize=12, fontweight='bold', loc='left')
        
        # Add common legend
        handles = [Line2D([0], [0], color=BRANCH_COLORS.get(b, '#333'), linewidth=2, 
                         label=b.replace('_', ' ')) for b in ESTUARY_BRANCHES]
        fig.legend(handles=handles, loc='lower center', ncol=4, fontsize=10,
                  bbox_to_anchor=(0.5, -0.02), frameon=True)
        
        fig.suptitle('C-GEM Mekong Delta Full Network Simulation Results', 
                    fontsize=16, fontweight='bold', y=0.98)
        
        if save:
            outpath = self.output_dir / 'combined_network_figure.png'
            plt.savefig(outpath, dpi=300, bbox_inches='tight', facecolor='white')
            print(f"Saved: {outpath}")
        
        return fig
    
    def _draw_network_schematic(self, ax):
        """Draw simplified network schematic in an axis."""
        ax.set_xlim(-1, 20)
        ax.set_ylim(-1, 8)
        ax.set_aspect('equal')
        ax.axis('off')
        
        # Simplified horizontal layout
        # Nodes: Tien inlet -> junction -> split -> outlets
        
        # Draw branches as curved/angled paths
        paths = [
            # (name, [(x1,y1), (x2,y2), ...], color)
            ('Tien_Main', [(0, 6), (4, 5)], '#1f77b4'),
            ('Hau_Main', [(0, 2), (4, 3)], '#ff7f0e'),
            ('Vam_Nao', [(4, 5), (4, 3)], '#2ca02c'),
            ('Tien_Lower', [(4, 5), (8, 5)], '#1f77b4'),
            ('Tien_Connector', [(8, 5), (10, 4)], '#9467bd'),
            ('Co_Chien', [(8, 5), (14, 6.5)], '#d62728'),
            ('My_Tho', [(10, 4), (16, 5.5)], '#8c564b'),
            ('Ham_Luong', [(10, 4), (16, 3)], '#e377c2'),
            ('Hau_River', [(4, 3), (18, 1)], '#ff7f0e'),
        ]
        
        for name, coords, color in paths:
            xs = [c[0] for c in coords]
            ys = [c[1] for c in coords]
            ax.plot(xs, ys, color=color, linewidth=6, solid_capstyle='round', alpha=0.7)
            
            # Add label
            mx, my = np.mean(xs), np.mean(ys)
            ax.text(mx, my + 0.4, name.replace('_', '\n'), fontsize=7, 
                   ha='center', va='bottom', color='black')
        
        # Add inlet/outlet markers
        ax.scatter([0, 0], [6, 2], c='blue', s=150, marker='s', zorder=10, label='Discharge Inlet')
        ax.scatter([14, 16, 16, 18], [6.5, 5.5, 3, 1], c='orange', s=150, marker='^', zorder=10, label='Tidal Outlet')
        ax.scatter([4, 4, 8, 10], [5, 3, 5, 4], c='gray', s=80, marker='o', zorder=10, label='Junction')
        
        # Add legend
        ax.legend(loc='lower right', fontsize=8)
        
        # Add arrow showing flow direction
        ax.annotate('', xy=(19, 4), xytext=(1, 4),
                   arrowprops=dict(arrowstyle='->', color='gray', lw=1.5))
        ax.text(10, 4.5, 'Flow direction →', fontsize=9, color='gray', ha='center')
    
    def _plot_variable_all_branches(self, ax, variable, ylabel, time_idx=-1, 
                                    estuary_only=False, ylim=None):
        """Plot a variable for all branches on a single axis."""
        branches = ESTUARY_BRANCHES if estuary_only else self.branches
        
        for branch in branches:
            df = self._load_branch_data(branch, variable)
            if df is None:
                continue
            
            dist_cols = [c for c in df.columns if c.endswith('km')]
            distances = [float(c.replace('km', '')) for c in dist_cols]
            
            if time_idx < 0:
                time_idx_use = len(df) + time_idx
            else:
                time_idx_use = min(time_idx, len(df) - 1)
            
            values = df.iloc[time_idx_use][dist_cols].values
            
            color = BRANCH_COLORS.get(branch, '#333333')
            ax.plot(distances, values, color=color, linewidth=1.5, 
                   label=branch.replace('_', ' '), alpha=0.8)
        
        ax.set_xlabel('Distance from mouth (km)')
        ax.set_ylabel(ylabel)
        if ylim:
            ax.set_ylim(ylim)
        ax.grid(True, alpha=0.3)
    
    def plot_calibration_summary(self, figsize=(12, 8), save=True):
        """Plot calibration results summary."""
        results = self._load_calibration_results()
        if results is None:
            print("No calibration results found")
            return None
        
        fig, axes = plt.subplots(1, 2, figsize=figsize)
        
        # Left: Parameter changes
        ax1 = axes[0]
        # This would need the parameter CSV from calibration
        ax1.text(0.5, 0.5, 'Parameter optimization results\n(requires calibration_results.csv)', 
                transform=ax1.transAxes, ha='center', va='center')
        ax1.set_title('Calibrated Parameters')
        
        # Right: Objective performance
        ax2 = axes[1]
        ax2.text(0.5, 0.5, 'Objective function performance\n(requires calibration_results.csv)', 
                transform=ax2.transAxes, ha='center', va='center')
        ax2.set_title('Calibration Fit')
        
        fig.suptitle('Calibration Summary', fontsize=14, fontweight='bold')
        plt.tight_layout()
        
        if save:
            outpath = self.output_dir / 'calibration_summary.png'
            plt.savefig(outpath, dpi=300, bbox_inches='tight', facecolor='white')
            print(f"Saved: {outpath}")
        
        return fig
    
    def generate_all_figures(self):
        """Generate all publication figures."""
        print(f"\nGenerating publication figures for {self.case_name}")
        print(f"Output directory: {self.output_dir}")
        print("-" * 60)
        
        # Network topology
        print("\n1. Network topology...")
        self.plot_network_topology()
        
        # Salinity profiles
        print("\n2. Salinity profiles...")
        self.plot_salinity_profiles()
        
        # SPM profiles
        print("\n3. SPM profiles...")
        self.plot_spm_profiles()
        
        # Water quality
        print("\n4. Water quality summary...")
        self.plot_water_quality_summary()
        
        # GHG emissions
        print("\n5. GHG emissions...")
        self.plot_ghg_emissions()
        
        # Combined figure
        print("\n6. Combined network figure...")
        self.plot_combined_network_figure()
        
        # Calibration summary
        print("\n7. Calibration summary...")
        self.plot_calibration_summary()
        
        print("\n" + "=" * 60)
        print(f"All figures saved to: {self.output_dir}")
        print("=" * 60)


def main():
    parser = argparse.ArgumentParser(
        description='Generate publication-quality figures for C-GEM Network',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument('--case', '-c', type=str, default='Mekong_Delta_Full',
                       help='Case name (default: Mekong_Delta_Full)')
    parser.add_argument('--output', '-o', type=str, default=None,
                       help='Output directory (default: OUTPUT/{case}/figures)')
    parser.add_argument('--time', '-t', type=int, default=-1,
                       help='Time index for snapshots (-1 for final state)')
    
    args = parser.parse_args()
    
    # Determine case directory
    case_dir = Path(f'INPUT/Cases/{args.case}')
    if not case_dir.exists():
        print(f"Error: Case directory not found: {case_dir}")
        sys.exit(1)
    
    # Create plotter and generate figures
    plotter = MekongNetworkPlotter(case_dir, args.output)
    plotter.generate_all_figures()
    
    plt.show()


if __name__ == '__main__':
    main()
