#!/usr/bin/env python3
"""
Update species boundary files to include 2-pool TOC model columns.

Scientific basis (December 2025 Audit Fix):
- RIVER: 15% labile, 85% refractory (dominated by terrestrial humics from Tonle Sap)
- OCEAN: 5% labile, 95% refractory (marine DOM is aged and refractory)

Also increases N2O at river boundary to reflect agricultural inputs
(upstream N2O should be ~25-30 nmol/L based on field observations).

Reference: Middelburg (1989), Hopkinson & Vallino (2005), Amon & Benner (1996)
"""
import os
from pathlib import Path

SCRIPT_DIR = Path(__file__).parent.resolve()
PROJECT_ROOT = SCRIPT_DIR.parent

def update_file_with_2pool_toc(filepath, labile_fraction, n2o_upstream_nmol=None):
    """Add toc_labile and toc_refractory columns to a species CSV file."""
    
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    # Find header and skip comment lines
    header_idx = 0
    for i, line in enumerate(lines):
        if line.strip() and not line.startswith('#'):
            header_idx = i
            break
    
    # Parse header
    header_line = lines[header_idx].strip()
    header = header_line.split(',')
    
    # Check if already has toc_labile
    if 'toc_labile' in header:
        print(f"  {filepath}: Already has toc_labile column, skipping")
        return
    
    # Find toc column index
    try:
        toc_idx = header.index('toc')
    except ValueError:
        print(f"  {filepath}: No 'toc' column found, skipping")
        return
    
    # Find n2o column index for potential update
    n2o_idx = None
    try:
        n2o_idx = header.index('n2o')
    except ValueError:
        pass
    
    # Create new header with toc_labile, toc_refractory after toc
    new_header = header[:toc_idx+1] + ['toc_labile', 'toc_refractory'] + header[toc_idx+1:]
    
    # Process data rows
    new_lines = []
    
    # Add comment lines before header
    for i in range(header_idx):
        new_lines.append(lines[i].rstrip())
    
    # Add new header
    new_lines.append(','.join(new_header))
    
    # Process data rows
    for line in lines[header_idx+1:]:
        stripped = line.strip()
        if not stripped or stripped.startswith('#'):
            new_lines.append(stripped)
            continue
        
        values = stripped.split(',')
        if len(values) < toc_idx + 1:
            new_lines.append(stripped)
            continue
        
        try:
            toc_total = float(values[toc_idx])
        except (ValueError, IndexError):
            new_lines.append(stripped)
            continue
        
        toc_labile = toc_total * labile_fraction
        toc_refractory = toc_total * (1.0 - labile_fraction)
        
        # Update N2O if specified (for river boundary)
        if n2o_upstream_nmol is not None and n2o_idx is not None:
            try:
                # Convert nmol/L to µmol/L (divide by 1000)
                values[n2o_idx] = str(round(n2o_upstream_nmol / 1000.0, 4))
            except (ValueError, IndexError):
                pass
        
        # Insert new columns
        new_values = values[:toc_idx+1] + [
            str(round(toc_labile, 1)), 
            str(round(toc_refractory, 1))
        ] + values[toc_idx+1:]
        
        new_lines.append(','.join(new_values))
    
    # Write updated file
    with open(filepath, 'w') as f:
        f.write('\n'.join(new_lines) + '\n')
    
    print(f"  Updated: {filepath}")
    print(f"    TOC split: {labile_fraction*100:.0f}% labile, {(1-labile_fraction)*100:.0f}% refractory")
    if n2o_upstream_nmol is not None:
        print(f"    N2O updated to: {n2o_upstream_nmol} nmol/L")


def main():
    case_dir = PROJECT_ROOT / "INPUT" / "Cases" / "Mekong_Delta_Full"
    
    print("=" * 60)
    print("Updating boundary files for 2-Pool TOC Model")
    print("Scientific Fix: December 2025 Audit")
    print("=" * 60)
    
    # Update river boundary: 15% labile, 85% refractory
    # Also increase N2O to 25 nmol/L (observed upstream)
    river_files = [
        case_dir / "species_river_realistic.csv",
        case_dir / "species_river.csv",
    ]
    
    for f in river_files:
        if f.exists():
            print(f"\nRiver boundary (15% labile, 85% refractory):")
            # N2O at 25 nmol/L = 0.025 µmol/L (high upstream due to agriculture)
            update_file_with_2pool_toc(f, labile_fraction=0.15, n2o_upstream_nmol=25.0)
    
    # Update ocean boundary: 5% labile, 95% refractory
    ocean_files = [
        case_dir / "species_ocean_realistic.csv",
        case_dir / "species_ocean.csv",
    ]
    
    for f in ocean_files:
        if f.exists():
            print(f"\nOcean boundary (5% labile, 95% refractory):")
            update_file_with_2pool_toc(f, labile_fraction=0.05)
    
    print("\n" + "=" * 60)
    print("Done! The 2-pool TOC model is now configured.")
    print("=" * 60)
    
    print("""
Scientific basis:
- River water has ~15% labile TOC (fresh phytoplankton, sewage)
  and ~85% refractory (Tonle Sap humics, aged terrestrial DOC)
- Ocean water has ~5% labile (minimal) and ~95% refractory
  (aged marine DOM, terrestrial humics transported to sea)

The OBSERVED near-conservative TOC behavior emerges naturally
because the labile fraction degrades quickly while the refractory
fraction is transported conservatively.

This replaces the previous "salinity switch" hack which was
scientifically incorrect (bacteria don't know about salinity!).

Reference: Middelburg (1989), Hopkinson & Vallino (2005)
""")


if __name__ == "__main__":
    main()
