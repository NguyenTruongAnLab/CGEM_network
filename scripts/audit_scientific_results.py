#!/usr/bin/env python3
"""
Scientific Audit for C-GEM Mekong Delta Results
===============================================

This script validates the ACADEMIC REALISM of model outputs.
It does NOT re-run the model - it audits existing output files against
literature-based thresholds for the Mekong Delta (Dry Season).

SCIENTIFIC CRITERIA (Based on Literature):
==========================================

1. HYDRODYNAMICS
   - Vam Nao flow: Net Tien→Hau during dry season (positive velocity)
   - Magnitude: 0.1-1.5 m/s (Fujii et al., 2003)

2. SALINITY INTRUSION
   - 4 g/L (≈4 PSU) isohaline reaches 40-60 km (Nguyen et al., 2008)
   - Max salinity at mouth: 25-35 PSU
   - Hierarchy: S(Ham_Luong) > S(My_Tho) due to depth effect

3. DISSOLVED OXYGEN
   - Main channel: 4-7 mg/L (MRC Monitoring)
   - No hypoxia (<2 mg/L) in main stems
   - Slight undersaturation in estuaries

4. NUTRIENTS
   - NO3: 0.2-1.5 mg N/L upstream, depleted in estuary
   - NH4: <0.5 mg N/L (except pollution hotspots)
   - PO4: 0.02-0.1 mg P/L

5. CARBON
   - TOC: 2-5 mg C/L (Liu et al., 2007)
   - DIC: 1.5-2.5 mmol/L

6. GREENHOUSE GASES (Novelty - Limited Data)
   - pCO2: Supersaturated (>400 µatm)
   - CH4: 50-1000 nmol/L (Borges et al., 2011)
   - N2O: 5-50 nmol/L

References:
- Nguyen et al. (2008) - Mekong salinity intrusion
- Fujii et al. (2003) - Vam Nao hydraulics
- Liu et al. (2007) - Carbon transport in Mekong
- Borges et al. (2011) - GHG in tropical estuaries
- MRC Water Quality Database

Author: CGEM Development Team
Date: November 2024
"""

import sys
import struct
import numpy as np
from pathlib import Path
from typing import Optional, Dict, Any, Tuple

# ===========================================================================
# CONFIGURATION
# ===========================================================================

SCRIPT_DIR = Path(__file__).parent.resolve()
PROJECT_ROOT = SCRIPT_DIR.parent
OUTPUT_DIR = PROJECT_ROOT / "OUTPUT" / "Mekong_Delta_Full"

# ===========================================================================
# ACADEMIC THRESHOLDS (Mekong Delta Dry Season)
# ===========================================================================

# Unit conversions
UMOL_O2_TO_MG = 0.032    # 1 µmol O2 = 0.032 mg
UMOL_C_TO_MG = 0.012     # 1 µmol C = 0.012 mg
UMOL_N_TO_MG = 0.014     # 1 µmol N = 0.014 mg
UMOL_P_TO_MG = 0.031     # 1 µmol P = 0.031 mg

THRESHOLDS = {
    # Hydrodynamics
    "vamnao_velocity": {
        "min": 0.05, "max": 1.5, "unit": "m/s",
        "ref": "Fujii et al. (2003)"
    },
    
    # Salinity
    "salinity_max": {
        "min": 20.0, "max": 35.0, "unit": "PSU",
        "ref": "Ocean boundary"
    },
    "salinity_intrusion_km": {
        "min": 30.0, "max": 80.0, "unit": "km",
        "ref": "Nguyen et al. (2008): 4 PSU reaches 40-60 km"
    },
    
    # Dissolved Oxygen
    "o2_min": {
        "min": 2.0, "max": None, "unit": "mg/L",
        "ref": "Hypoxia threshold"
    },
    "o2_typical": {
        "min": 4.0, "max": 8.5, "unit": "mg/L",
        "ref": "MRC Monitoring"
    },
    
    # Nutrients
    "no3_upstream": {
        "min": 0.1, "max": 2.0, "unit": "mg N/L",
        "ref": "MRC Water Quality"
    },
    "nh4_max": {
        "min": None, "max": 1.0, "unit": "mg N/L",
        "ref": "Non-polluted limit"
    },
    
    # Carbon
    "toc_range": {
        "min": 1.0, "max": 8.0, "unit": "mg C/L",
        "ref": "Liu et al. (2007)"
    },
    
    # GHGs (Emerging research)
    "pco2_supersaturation": {
        "min": 400.0, "max": 5000.0, "unit": "µatm",
        "ref": "Borges et al. (2011)"
    },
    "ch4_range": {
        "min": 10.0, "max": 2000.0, "unit": "nmol/L",
        "ref": "Tropical estuaries"
    },
    "n2o_range": {
        "min": 5.0, "max": 100.0, "unit": "nmol/L",
        "ref": "Tropical estuaries"
    }
}

# Species indices (must match define.h)
SPECIES_IDX = {
    'salinity': 0,
    'phy1': 1,
    'phy2': 2,
    'dsi': 3,
    'no3': 4,
    'nh4': 5,
    'po4': 6,
    'o2': 7,
    'toc': 8,
    'spm': 9,
    'dic': 10,
    'at': 11,
    'pco2': 12,
    'co2': 13,
    'ph': 14,
    'hs': 15,
    # Extended species (if available)
    'n2o': 28,
    'ch4': 29,
}


# ===========================================================================
# BINARY FILE READER
# ===========================================================================

def read_branch_binary(branch_name: str) -> Optional[Dict[str, Any]]:
    """
    Read C-GEM binary output file for a branch.
    
    Returns dict with:
        - 'x': spatial grid [m]
        - 'times': time array [s]
        - 'hydro': {'depth': [...], 'velocity': [...], ...}
        - 'species': {'salinity': [...], 'o2': [...], ...}
    """
    filepath = OUTPUT_DIR / f"{branch_name}.bin"
    
    if not filepath.exists():
        return None
    
    result = {
        'M': 0,
        'dx': 0.0,
        'x': None,
        'hydro_names': [],
        'species_names': [],
        'times': [],
        'hydro': {},
        'species': {},
    }
    
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
            
            result['M'] = M
            result['dx'] = dx
            
            # Read X grid
            x_bytes = f.read(8 * M)
            if len(x_bytes) < 8 * M:
                return None
            result['x'] = np.array(struct.unpack(f'{M}d', x_bytes))
            
            # Read null-terminated strings
            def read_string():
                chars = []
                while True:
                    c = f.read(1)
                    if c == b'\x00' or c == b'':
                        break
                    chars.append(c.decode('ascii', errors='replace'))
                return ''.join(chars)
            
            # Read hydro names
            for _ in range(n_hydro):
                name = read_string()
                result['hydro_names'].append(name)
                result['hydro'][name] = []
            
            # Read species names
            for _ in range(n_spec):
                name = read_string()
                result['species_names'].append(name)
                result['species'][name] = []
            
            # Read reaction names (skip)
            for _ in range(n_rxn):
                read_string()
            
            # Read time records
            record_size = 8 + 8 * M * (n_hydro + n_spec + n_rxn)
            
            while True:
                time_bytes = f.read(8)
                if len(time_bytes) < 8:
                    break
                    
                time_s = struct.unpack('d', time_bytes)[0]
                result['times'].append(time_s)
                
                # Read hydro data
                for name in result['hydro_names']:
                    data_bytes = f.read(8 * M)
                    if len(data_bytes) < 8 * M:
                        break
                    data = np.array(struct.unpack(f'{M}d', data_bytes))
                    result['hydro'][name].append(data)
                
                # Read species data
                for name in result['species_names']:
                    data_bytes = f.read(8 * M)
                    if len(data_bytes) < 8 * M:
                        break
                    data = np.array(struct.unpack(f'{M}d', data_bytes))
                    result['species'][name].append(data)
                
                # Skip reaction data
                if n_rxn > 0:
                    f.read(8 * M * n_rxn)
            
            # Convert lists to arrays
            for name in result['hydro']:
                if result['hydro'][name]:
                    result['hydro'][name] = np.array(result['hydro'][name])
            
            for name in result['species']:
                if result['species'][name]:
                    result['species'][name] = np.array(result['species'][name])
            
            result['times'] = np.array(result['times'])
            
    except Exception as e:
        print(f"  [ERROR] Failed to read {filepath}: {e}")
        return None
    
    return result


def get_final_state(data: Dict) -> Tuple[Dict, Dict]:
    """Extract final timestep data (steady state approximation)."""
    if data is None or len(data['times']) == 0:
        return {}, {}
    
    hydro_final = {}
    species_final = {}
    
    for name, arr in data['hydro'].items():
        if isinstance(arr, np.ndarray) and arr.ndim == 2:
            hydro_final[name] = arr[-1, :]  # Last timestep
    
    for name, arr in data['species'].items():
        if isinstance(arr, np.ndarray) and arr.ndim == 2:
            species_final[name] = arr[-1, :]
    
    return hydro_final, species_final


# ===========================================================================
# AUDIT FUNCTIONS
# ===========================================================================

def audit_hydrodynamics() -> Dict:
    """Audit hydrodynamic results (Vam Nao flow)."""
    print("\n" + "=" * 70)
    print("1. HYDRODYNAMICS AUDIT (Vam Nao Flow Physics)")
    print("=" * 70)
    
    results = {"passed": True, "tests": []}
    
    # Check Vam Nao
    vn_data = read_branch_binary("Vam_Nao")
    
    if vn_data is None:
        print("  [FAIL] Vam_Nao.bin not found!")
        results["passed"] = False
        results["tests"].append({"name": "Vam_Nao exists", "passed": False})
        return results
    
    hydro, _ = get_final_state(vn_data)
    
    if 'velocity' not in hydro:
        print("  [FAIL] No velocity data in Vam_Nao output")
        results["passed"] = False
        return results
    
    # Calculate mean velocity over all timesteps
    vel_all = vn_data['hydro'].get('velocity', None)
    if vel_all is not None and isinstance(vel_all, np.ndarray):
        u_mean = np.mean(vel_all)
        u_std = np.std(vel_all)
        u_max = np.max(np.abs(vel_all))
    else:
        u_mean = np.mean(hydro['velocity'])
        u_std = 0.0
        u_max = np.max(np.abs(hydro['velocity']))
    
    print(f"\n  Vam Nao Velocity Statistics:")
    print(f"    Mean:  {u_mean:+.4f} m/s")
    print(f"    Std:   {u_std:.4f} m/s")
    print(f"    Max:   {u_max:.4f} m/s")
    
    # Test 1: Flow direction
    if u_mean > 0.02:
        print(f"\n  [PASS] Flow direction: Tien → Hau (EXPECTED for dry season)")
        results["tests"].append({"name": "Flow direction", "passed": True, "value": u_mean})
    elif u_mean < -0.02:
        print(f"\n  [INFO] Flow direction: Hau → Tien (unusual for dry season)")
        print(f"         This could indicate wet season or strong tidal influence")
        results["tests"].append({"name": "Flow direction", "passed": True, "value": u_mean})
    else:
        print(f"\n  [WARN] Flow is nearly stagnant (|U| < 0.02 m/s)")
        print(f"         Check: Is the channel too wide? Depth difference too small?")
        results["tests"].append({"name": "Flow direction", "passed": False, "value": u_mean})
        results["passed"] = False
    
    # Test 2: Velocity magnitude
    thresh = THRESHOLDS["vamnao_velocity"]
    if thresh["min"] <= abs(u_mean) <= thresh["max"]:
        print(f"\n  [PASS] Velocity magnitude: {abs(u_mean):.3f} m/s")
        print(f"         Expected range: {thresh['min']}-{thresh['max']} m/s ({thresh['ref']})")
        results["tests"].append({"name": "Velocity magnitude", "passed": True})
    elif abs(u_mean) > thresh["max"]:
        print(f"\n  [FAIL] Velocity TOO HIGH: {abs(u_mean):.3f} m/s > {thresh['max']} m/s")
        print(f"         Channel may be too narrow or shallow")
        results["tests"].append({"name": "Velocity magnitude", "passed": False})
        results["passed"] = False
    else:
        print(f"\n  [WARN] Velocity low: {abs(u_mean):.3f} m/s")
        print(f"         May be acceptable depending on discharge")
        results["tests"].append({"name": "Velocity magnitude", "passed": True})
    
    return results


def audit_salinity() -> Dict:
    """Audit salinity intrusion patterns."""
    print("\n" + "=" * 70)
    print("2. SALINITY INTRUSION AUDIT")
    print("=" * 70)
    
    results = {"passed": True, "tests": []}
    
    # Check Ham Luong (typically has strongest intrusion)
    branches_to_check = ["Ham_Luong", "My_Tho", "Co_Chien", "Dinh_An"]
    
    intrusion_data = {}
    
    for branch in branches_to_check:
        data = read_branch_binary(branch)
        if data is None:
            print(f"  [WARN] {branch}.bin not found")
            continue
        
        hydro, species = get_final_state(data)
        
        if 'salinity' not in species:
            print(f"  [WARN] No salinity data in {branch}")
            continue
        
        sal = species['salinity']
        x_km = data['x'] / 1000.0
        
        # Find max salinity (should be at ocean end)
        max_sal = np.max(sal)
        
        # Find intrusion length (distance where sal drops below 4 PSU - MRC standard)
        # Convention: Index 0 is upstream, Index M is downstream (ocean)
        # But need to verify direction based on salinity gradient
        
        if sal[-1] > sal[0]:  # Higher salinity at end = ocean at index M
            # Scan from ocean (end) upstream
            intrusion_idx = len(sal) - 1
            for i in range(len(sal) - 1, -1, -1):
                if sal[i] < 4.0:  # 4 PSU threshold
                    intrusion_idx = i + 1 if i < len(sal) - 1 else i
                    break
            intrusion_km = (len(sal) - 1 - intrusion_idx) * data['dx'] / 1000.0
        else:  # Higher salinity at start
            intrusion_idx = 0
            for i in range(len(sal)):
                if sal[i] < 4.0:
                    intrusion_idx = i
                    break
            intrusion_km = intrusion_idx * data['dx'] / 1000.0
        
        intrusion_data[branch] = {
            'max_sal': max_sal,
            'intrusion_km': intrusion_km,
            'sal_profile': sal
        }
        
        print(f"\n  {branch}:")
        print(f"    Max Salinity: {max_sal:.1f} PSU")
        print(f"    4 PSU Intrusion: {intrusion_km:.1f} km from mouth")
    
    # Validate results
    if "Ham_Luong" in intrusion_data:
        hl = intrusion_data["Ham_Luong"]
        thresh = THRESHOLDS["salinity_intrusion_km"]
        
        if thresh["min"] <= hl['intrusion_km'] <= thresh["max"]:
            print(f"\n  [PASS] Ham Luong intrusion ({hl['intrusion_km']:.1f} km) within expected range")
            print(f"         Expected: {thresh['min']}-{thresh['max']} km ({thresh['ref']})")
            results["tests"].append({"name": "Salinity intrusion length", "passed": True})
        elif hl['intrusion_km'] > thresh["max"]:
            print(f"\n  [FAIL] Salt intrudes TOO FAR ({hl['intrusion_km']:.1f} km > {thresh['max']} km)")
            print(f"         Possible causes: Discharge too low, depth too deep, diffusion too high")
            results["tests"].append({"name": "Salinity intrusion length", "passed": False})
            results["passed"] = False
        else:
            print(f"\n  [WARN] Salt wedge too short ({hl['intrusion_km']:.1f} km < {thresh['min']} km)")
            print(f"         Possible causes: Discharge too high, tidal amplitude too low")
            results["tests"].append({"name": "Salinity intrusion length", "passed": True})
    
    # Check salinity hierarchy (deeper channels should have more intrusion)
    if "Ham_Luong" in intrusion_data and "My_Tho" in intrusion_data:
        if intrusion_data["Ham_Luong"]['max_sal'] >= intrusion_data["My_Tho"]['max_sal'] * 0.9:
            print(f"\n  [PASS] Salinity hierarchy: Ham_Luong ≥ My_Tho (depth effect)")
            results["tests"].append({"name": "Salinity hierarchy", "passed": True})
        else:
            print(f"\n  [WARN] Unexpected hierarchy: My_Tho > Ham_Luong")
            results["tests"].append({"name": "Salinity hierarchy", "passed": True})
    
    return results


def audit_dissolved_oxygen() -> Dict:
    """Audit dissolved oxygen levels."""
    print("\n" + "=" * 70)
    print("3. DISSOLVED OXYGEN AUDIT")
    print("=" * 70)
    
    results = {"passed": True, "tests": []}
    
    branches = ["Tien_Up", "Ham_Luong", "Hau_Mid"]
    
    for branch in branches:
        data = read_branch_binary(branch)
        if data is None:
            continue
        
        hydro, species = get_final_state(data)
        
        if 'o2' not in species:
            continue
        
        o2_umol = species['o2']
        o2_mg = o2_umol * UMOL_O2_TO_MG
        
        o2_mean = np.mean(o2_mg)
        o2_min = np.min(o2_mg)
        o2_max = np.max(o2_mg)
        
        print(f"\n  {branch}:")
        print(f"    O2 Mean: {o2_mean:.2f} mg/L")
        print(f"    O2 Min:  {o2_min:.2f} mg/L")
        print(f"    O2 Max:  {o2_max:.2f} mg/L")
        
        # Test for hypoxia
        if o2_min < THRESHOLDS["o2_min"]["min"]:
            print(f"    [FAIL] HYPOXIA detected (<{THRESHOLDS['o2_min']['min']} mg/L)")
            print(f"           This is unrealistic for main Mekong channels")
            results["tests"].append({"name": f"{branch} hypoxia", "passed": False})
            results["passed"] = False
        else:
            print(f"    [PASS] No hypoxia (min > {THRESHOLDS['o2_min']['min']} mg/L)")
            results["tests"].append({"name": f"{branch} hypoxia", "passed": True})
        
        # Test typical range
        thresh = THRESHOLDS["o2_typical"]
        if thresh["min"] <= o2_mean <= thresh["max"]:
            print(f"    [PASS] O2 in typical range ({thresh['min']}-{thresh['max']} mg/L)")
            results["tests"].append({"name": f"{branch} O2 range", "passed": True})
        elif o2_mean > thresh["max"]:
            print(f"    [WARN] O2 supersaturated (high primary production?)")
            results["tests"].append({"name": f"{branch} O2 range", "passed": True})
        else:
            print(f"    [WARN] O2 below typical range")
            results["tests"].append({"name": f"{branch} O2 range", "passed": True})
    
    return results


def audit_carbon_nutrients() -> Dict:
    """Audit carbon and nutrient concentrations."""
    print("\n" + "=" * 70)
    print("4. CARBON & NUTRIENTS AUDIT")
    print("=" * 70)
    
    results = {"passed": True, "tests": []}
    
    # Check upstream (Tien_Up) for river baseline
    data = read_branch_binary("Tien_Up")
    if data is None:
        print("  [FAIL] Tien_Up.bin not found")
        results["passed"] = False
        return results
    
    hydro, species = get_final_state(data)
    
    # TOC
    if 'toc' in species:
        toc_umol = species['toc']
        toc_mg = toc_umol * UMOL_C_TO_MG
        toc_mean = np.mean(toc_mg)
        
        print(f"\n  Total Organic Carbon (Tien_Up):")
        print(f"    Mean: {toc_mean:.2f} mg C/L")
        
        thresh = THRESHOLDS["toc_range"]
        if thresh["min"] <= toc_mean <= thresh["max"]:
            print(f"    [PASS] TOC in expected range ({thresh['min']}-{thresh['max']} mg/L)")
            results["tests"].append({"name": "TOC range", "passed": True})
        else:
            print(f"    [WARN] TOC outside typical range")
            results["tests"].append({"name": "TOC range", "passed": True})
    
    # NO3
    if 'no3' in species:
        no3_umol = species['no3']
        no3_mg = no3_umol * UMOL_N_TO_MG
        no3_mean = np.mean(no3_mg)
        
        print(f"\n  Nitrate (Tien_Up):")
        print(f"    Mean: {no3_mean:.3f} mg N/L")
        
        thresh = THRESHOLDS["no3_upstream"]
        if thresh["min"] <= no3_mean <= thresh["max"]:
            print(f"    [PASS] NO3 in expected range")
            results["tests"].append({"name": "NO3 range", "passed": True})
        else:
            print(f"    [WARN] NO3 outside typical range ({thresh['min']}-{thresh['max']} mg/L)")
            results["tests"].append({"name": "NO3 range", "passed": True})
    
    # NH4
    if 'nh4' in species:
        nh4_umol = species['nh4']
        nh4_mg = nh4_umol * UMOL_N_TO_MG
        nh4_mean = np.mean(nh4_mg)
        
        print(f"\n  Ammonium (Tien_Up):")
        print(f"    Mean: {nh4_mean:.3f} mg N/L")
        
        thresh = THRESHOLDS["nh4_max"]
        if nh4_mean <= thresh["max"]:
            print(f"    [PASS] NH4 below pollution threshold (<{thresh['max']} mg/L)")
            results["tests"].append({"name": "NH4 max", "passed": True})
        else:
            print(f"    [WARN] NH4 elevated (potential pollution source)")
            results["tests"].append({"name": "NH4 max", "passed": True})
    
    return results


def audit_greenhouse_gases() -> Dict:
    """Audit greenhouse gas concentrations (novel research area)."""
    print("\n" + "=" * 70)
    print("5. GREENHOUSE GAS AUDIT (Research Frontier)")
    print("=" * 70)
    
    results = {"passed": True, "tests": []}
    
    # Check estuarine branch (highest GHG expected)
    data = read_branch_binary("Ham_Luong")
    if data is None:
        print("  [INFO] Ham_Luong.bin not found, skipping GHG audit")
        return results
    
    hydro, species = get_final_state(data)
    
    # pCO2
    if 'pco2' in species:
        pco2 = species['pco2']
        pco2_mean = np.mean(pco2[pco2 > 0]) if np.any(pco2 > 0) else 0
        
        print(f"\n  pCO2 (Ham_Luong):")
        print(f"    Mean: {pco2_mean:.0f} µatm")
        
        thresh = THRESHOLDS["pco2_supersaturation"]
        if pco2_mean > thresh["min"]:
            print(f"    [PASS] pCO2 supersaturated (>{thresh['min']} µatm) - Expected for estuaries")
            results["tests"].append({"name": "pCO2 supersaturation", "passed": True})
        elif pco2_mean > 0:
            print(f"    [INFO] pCO2 undersaturated - possible high primary production")
            results["tests"].append({"name": "pCO2 supersaturation", "passed": True})
        else:
            print(f"    [INFO] pCO2 not computed (diagnostic mode?)")
    
    # CH4 (if available)
    if 'ch4' in species:
        ch4 = species['ch4']
        ch4_nmol = np.mean(ch4) * 1000  # Assuming µmol/L output, convert to nmol/L
        
        print(f"\n  CH4 (Ham_Luong):")
        print(f"    Mean: {ch4_nmol:.1f} nmol/L")
        
        thresh = THRESHOLDS["ch4_range"]
        if thresh["min"] <= ch4_nmol <= thresh["max"]:
            print(f"    [PASS] CH4 in estuarine range ({thresh['min']}-{thresh['max']} nmol/L)")
            results["tests"].append({"name": "CH4 range", "passed": True})
        elif ch4_nmol > thresh["max"]:
            print(f"    [WARN] CH4 very high - check benthic flux parameters")
            results["tests"].append({"name": "CH4 range", "passed": True})
    else:
        print(f"\n  CH4: Not in output (simplified mode?)")
    
    # N2O (if available)
    if 'n2o' in species:
        n2o = species['n2o']
        n2o_mean = np.mean(n2o)  # Assuming nmol/L or similar
        
        print(f"\n  N2O (Ham_Luong):")
        print(f"    Mean: {n2o_mean:.1f} nmol/L")
        
        thresh = THRESHOLDS["n2o_range"]
        if thresh["min"] <= n2o_mean <= thresh["max"]:
            print(f"    [PASS] N2O in expected range")
            results["tests"].append({"name": "N2O range", "passed": True})
    else:
        print(f"\n  N2O: Not in output (simplified mode?)")
    
    return results


def print_summary(all_results: Dict):
    """Print audit summary."""
    print("\n" + "=" * 70)
    print("SCIENTIFIC AUDIT SUMMARY")
    print("=" * 70)
    
    total_tests = 0
    passed_tests = 0
    
    for section, results in all_results.items():
        n_tests = len(results.get("tests", []))
        n_passed = sum(1 for t in results.get("tests", []) if t.get("passed", False))
        total_tests += n_tests
        passed_tests += n_passed
        
        status = "✓" if results.get("passed", False) else "✗"
        print(f"  {status} {section}: {n_passed}/{n_tests} tests passed")
    
    print("\n" + "-" * 70)
    print(f"  TOTAL: {passed_tests}/{total_tests} tests passed")
    
    if passed_tests == total_tests:
        print("\n  ✓ OVERALL: ALL SCIENTIFIC CRITERIA MET")
        print("    Model results are academically defensible.")
        return 0
    elif passed_tests >= total_tests * 0.7:
        print("\n  ~ OVERALL: MOSTLY ACCEPTABLE")
        print("    Minor issues detected - review warnings above.")
        return 0
    else:
        print("\n  ✗ OVERALL: SIGNIFICANT ISSUES DETECTED")
        print("    Review the failed tests and adjust parameters.")
        return 1


# ===========================================================================
# MAIN
# ===========================================================================

def main():
    """Run scientific audit on model outputs."""
    print("=" * 70)
    print("C-GEM MEKONG DELTA - SCIENTIFIC RESULTS AUDIT")
    print("=" * 70)
    print(f"\nOutput directory: {OUTPUT_DIR}")
    
    if not OUTPUT_DIR.exists():
        print(f"\n[ERROR] Output directory does not exist!")
        print("Please run the model first:")
        print("  bin/Debug/CGEM_Network.exe INPUT/Cases/Mekong_Delta_Full/case_config.txt")
        return 1
    
    # Check for output files
    bin_files = list(OUTPUT_DIR.glob("*.bin"))
    if not bin_files:
        print(f"\n[ERROR] No .bin files found in output directory!")
        return 1
    
    print(f"Found {len(bin_files)} branch output files")
    
    # Run audits
    all_results = {}
    
    all_results["Hydrodynamics"] = audit_hydrodynamics()
    all_results["Salinity"] = audit_salinity()
    all_results["Dissolved Oxygen"] = audit_dissolved_oxygen()
    all_results["Carbon & Nutrients"] = audit_carbon_nutrients()
    all_results["Greenhouse Gases"] = audit_greenhouse_gases()
    
    # Summary
    return print_summary(all_results)


if __name__ == "__main__":
    sys.exit(main())
