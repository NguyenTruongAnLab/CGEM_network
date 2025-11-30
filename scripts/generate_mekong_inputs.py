#!/usr/bin/env python3
"""
CGEM Mekong Delta Input Generator - Realistic Data-Driven Version

This script generates model inputs using ACTUAL available data sources:
1. Monitoring center data (BOD5, COD, NH4, NO3, pH, O2, TSS)
2. Field campaign snapshots (POC, DOC, TOC, Chl-a, salinity, alkalinity, CO2, CH4, N2O)
3. Sentinel satellite products (Chl-a, SPM, DOC, CDOM)
4. MRC hydrology data (discharge, tidal levels)

The script implements the 80/20 simplification approach for Mekong Delta,
where monitoring is less extensive than European rivers.

Usage:
    python generate_mekong_inputs.py --case Tien_River --start 2024-01-01 --days 365

Author: CGEM Development Team
Date: 2024-2025
"""

import numpy as np
import pandas as pd
from datetime import datetime, timedelta
from pathlib import Path
import argparse
import json
import warnings

# ============================================================================
# MEKONG DELTA OBSERVED DATA RANGES
# Based on: Monitoring center + 2024-2025 field campaigns
# ============================================================================

# Monitoring center typical values (converted to model units)
MONITORING_DATA = {
    # From routine monitoring (monthly/quarterly)
    'BOD5_mg_L': {'dry': (2, 5), 'wet': (3, 8)},           # mg O2/L
    'COD_mg_L': {'dry': (8, 15), 'wet': (12, 25)},         # mg O2/L
    'NH4_mg_N_L': {'dry': (0.1, 0.5), 'wet': (0.2, 0.8)},  # mg N/L
    'NO3_mg_N_L': {'dry': (0.3, 1.2), 'wet': (0.5, 2.0)},  # mg N/L
    'O2_mg_L': {'dry': (5.0, 7.5), 'wet': (4.5, 7.0)},     # mg/L
    'pH': {'dry': (7.0, 7.8), 'wet': (6.8, 7.5)},
    'TSS_mg_L': {'dry': (30, 80), 'wet': (80, 250)},       # mg/L
}

# Field campaign 2024-2025 data (dry + wet season snapshots)
FIELD_CAMPAIGN = {
    'dry_season': {  # Jan-Apr
        'salinity_PSU': {'upstream': (0, 0.5), 'midstream': (2, 15), 'mouth': (20, 32)},
        'TOC_mg_L': (2.5, 4.5),       # mg C/L
        'DOC_mg_L': (1.8, 3.5),       # mg C/L
        'POC_mg_L': (0.5, 1.5),       # mg C/L
        'Chla_ug_L': (3, 15),         # µg/L
        'alkalinity_meq_L': (1.5, 2.5),  # meq/L
        'CO2_uatm': (1000, 4000),     # µatm (supersaturated)
        'CH4_nmol_L': (50, 500),      # nmol/L
        'N2O_nmol_L': (8, 25),        # nmol/L
        'DSi_mg_L': (3, 8),           # mg Si/L
    },
    'wet_season': {  # Jul-Oct
        'salinity_PSU': {'upstream': (0, 0.2), 'midstream': (0, 3), 'mouth': (5, 20)},
        'TOC_mg_L': (3.5, 6.0),
        'DOC_mg_L': (2.5, 4.5),
        'POC_mg_L': (1.0, 2.5),
        'Chla_ug_L': (2, 8),          # Lower due to turbidity
        'alkalinity_meq_L': (1.2, 2.0),
        'CO2_uatm': (1500, 5000),
        'CH4_nmol_L': (100, 800),
        'N2O_nmol_L': (10, 35),
        'DSi_mg_L': (4, 10),
    }
}

# Satellite-derived ranges (Sentinel-3 OLCI)
SATELLITE_DATA = {
    'Chla_ug_L': (1, 30),         # Weekly composites
    'SPM_mg_L': (10, 300),        # TSM product
    'CDOM_1_m': (0.5, 3.0),       # aCDOM(443)
}

# MRC hydrological data
HYDROLOGY = {
    'Tien_River': {
        'Q_m3s': {'dry': (1500, 3500), 'wet': (6000, 14000), 'transition': (3000, 6000)},
        'tide_amplitude_m': (1.5, 2.5),
    },
    'branches': {
        'My_Tho': {'Q_fraction': 0.35, 'tide_damping': 1.0},
        'Ham_Luong': {'Q_fraction': 0.35, 'tide_damping': 0.95},
        'Co_Chien': {'Q_fraction': 0.30, 'tide_damping': 0.90},
    }
}


# ============================================================================
# UNIT CONVERSION FUNCTIONS
# ============================================================================

def bod5_to_toc(bod5_mg_L: float) -> float:
    """Convert BOD5 (mg O2/L) to TOC (µmol C/L).
    Assumes BOD5:TOC ratio ≈ 1.5 for tropical waters."""
    toc_mg_L = bod5_mg_L / 1.5
    return toc_mg_L * 1000 / 12  # mg/L to µmol/L

def cod_to_toc(cod_mg_L: float) -> float:
    """Convert COD (mg O2/L) to TOC (µmol C/L).
    Assumes COD:TOC ratio ≈ 2.67."""
    toc_mg_L = cod_mg_L / 2.67
    return toc_mg_L * 1000 / 12


def bod_cod_to_crive_fractions(bod5: float, cod: float, tss: float) -> dict:
    """
    Maps Engineering Units (BOD/COD/TSS) to C-RIVE Biogeochemical Fractions.
    
    This is the critical function that bridges monitoring data to model inputs.
    Based on typical ratios for tropical rivers (Servais et al., 1999).
    
    Args:
        bod5: BOD5 in mg O2/L
        cod: COD in mg O2/L  
        tss: Total suspended solids in mg/L
        
    Returns:
        dict with HD1, HD2, HD3, HP1, HP2, HP3 in µmol C/L
        
    The Logic:
    1. COD ≈ Ultimate Oxygen Demand (UOD)
    2. TOC (mass) ≈ UOD / 2.67 (stoichiometry of O2:C oxidation)
    3. BOD5 ≈ 60-70% of Biodegradable COD (BCOD) in warm waters
    4. Biodegradable fraction → Labile pools (HD1/HP1)
    5. Refractory fraction → HD3/HP3
    6. Particulate/dissolved split based on TSS (sorption model)
    """
    # 1. Total Organic Carbon (Estimate from COD)
    # COD is roughly equivalent to Ultimate Oxygen Demand (UOD)
    # C (mass) ≈ UOD / 2.67 (stoichiometry of O2:C oxidation)
    toc_total_umol = (cod / 2.67) * (1000 / 12.0)
    
    # 2. Biodegradable Fraction (Labile)
    # BOD5 is approx 60-70% of Biodegradable COD (BCOD) in warm waters
    # In tropical conditions, use 65% (faster kinetics)
    bcod = bod5 / 0.65
    biodegradable_C_umol = (bcod / 2.67) * (1000 / 12.0)
    
    # Ensure biodegradable doesn't exceed total
    biodegradable_C_umol = min(biodegradable_C_umol, toc_total_umol * 0.9)
    
    # 3. Partitioning (The "Input Mapper")
    # Particulate fraction based on TSS (generic sorption model)
    # Higher TSS = more organic matter sorbed to particles
    # Typical range: 20-80% particulate
    f_particulate = min(0.8, max(0.2, tss * 0.004))  # 0.4% per mg/L TSS
    
    # 4. LABILE POOLS (HD1/HP1) - Biodegradable fraction
    # These are rapidly consumed (half-life ~1 day)
    hp1 = biodegradable_C_umol * f_particulate
    hd1 = biodegradable_C_umol * (1.0 - f_particulate)
    
    # 5. REFRACTORY POOLS (HD3/HP3) - The rest of the carbon
    refractory_C = max(0.0, toc_total_umol - biodegradable_C_umol)
    hp3 = refractory_C * f_particulate
    hd3 = refractory_C * (1.0 - f_particulate)
    
    # 6. SEMI-LABILE (HD2/HP2) - Transitional pools
    # Usually small - steal 15% from refractory (these have ~3-day half-life)
    hp2 = hp3 * 0.15
    hd2 = hd3 * 0.15
    hp3 *= 0.85
    hd3 *= 0.85
    
    return {
        'hd1': hd1,  # Labile dissolved OC [µmol C/L]
        'hd2': hd2,  # Semi-labile dissolved OC [µmol C/L]
        'hd3': hd3,  # Refractory dissolved OC [µmol C/L]
        'hp1': hp1,  # Labile particulate OC [µmol C/L]
        'hp2': hp2,  # Semi-labile particulate OC [µmol C/L]
        'hp3': hp3,  # Refractory particulate OC [µmol C/L]
        'toc_total': toc_total_umol,
        'biodegradable_fraction': biodegradable_C_umol / max(toc_total_umol, 1.0),
        'particulate_fraction': f_particulate,
    }


def mg_N_to_umol(mg_N_L: float) -> float:
    """Convert mg N/L to µmol N/L."""
    return mg_N_L * 1000 / 14

def mg_P_to_umol(mg_P_L: float) -> float:
    """Convert mg P/L to µmol P/L."""
    return mg_P_L * 1000 / 31

def mg_Si_to_umol(mg_Si_L: float) -> float:
    """Convert mg Si/L to µmol Si/L."""
    return mg_Si_L * 1000 / 28

def mg_O2_to_umol(mg_O2_L: float) -> float:
    """Convert mg O2/L to µmol O2/L."""
    return mg_O2_L * 1000 / 32

def meq_to_umol_alk(meq_L: float) -> float:
    """Convert meq/L alkalinity to µmol/L."""
    return meq_L * 1000

def chla_to_phy(chla_ug_L: float, C_Chl_ratio: float = 50) -> float:
    """Convert Chl-a (µg/L) to phytoplankton carbon (µgC/L)."""
    return chla_ug_L * C_Chl_ratio


# ============================================================================
# SEASONAL FUNCTIONS (MEKONG-SPECIFIC)
# ============================================================================

def get_season(date: datetime) -> str:
    """Determine Mekong season from date.
    Dry: Dec-Apr, Wet: Jul-Oct, Transition: May-Jun, Nov"""
    month = date.month
    if month in [12, 1, 2, 3, 4]:
        return 'dry'
    elif month in [7, 8, 9, 10]:
        return 'wet'
    else:
        return 'transition'

def seasonal_factor(date: datetime) -> dict:
    """Calculate seasonal modulation factors for Mekong Delta."""
    doy = date.timetuple().tm_yday
    
    # Wet season peaks around day 245 (early September)
    wet_phase = 2 * np.pi * (doy - 245) / 365
    wet_factor = 0.5 * (1 + np.cos(wet_phase))
    
    # Temperature peaks in April-May (day 120)
    temp_phase = 2 * np.pi * (doy - 120) / 365
    temp_factor = 0.5 * (1 + np.cos(temp_phase))
    
    return {
        'discharge': 0.2 + 0.8 * wet_factor,
        'spm': 0.15 + 0.85 * wet_factor,
        'nutrients': 0.3 + 0.7 * wet_factor,
        'temperature': 0.4 + 0.6 * temp_factor,
        'salinity_intrusion': 1.0 - 0.8 * wet_factor,  # Less intrusion in wet
    }


def interpolate_seasonal(dry_range: tuple, wet_range: tuple, 
                         date: datetime, noise_pct: float = 0.1) -> float:
    """Interpolate between dry and wet season values."""
    sf = seasonal_factor(date)
    
    if sf['discharge'] > 0.6:  # Wet season dominant
        low, high = wet_range
    elif sf['discharge'] < 0.4:  # Dry season dominant
        low, high = dry_range
    else:  # Transition
        # Blend the ranges
        low = dry_range[0] * (1 - sf['discharge']) + wet_range[0] * sf['discharge']
        high = dry_range[1] * (1 - sf['discharge']) + wet_range[1] * sf['discharge']
    
    value = np.random.uniform(low, high)
    noise = 1.0 + noise_pct * np.random.randn()
    return value * max(0.5, min(1.5, noise))


# ============================================================================
# BOUNDARY CONDITION GENERATORS
# ============================================================================

def generate_river_boundary(start_date: str, duration_days: int, 
                            dt_hours: float = 6.0,
                            include_crive_fractions: bool = True) -> pd.DataFrame:
    """
    Generate river (upstream) boundary conditions from monitoring data.
    
    Uses: BOD5/COD → TOC, NH4, NO3, O2, TSS from monitoring center
    Plus: Field campaign data for alkalinity, CO2, DOC fractionation
    
    Generates ALL 30 species required by CGEM_NUM_SPECIES in define.h.
    
    Args:
        start_date: Start date string (YYYY-MM-DD)
        duration_days: Number of days to generate
        dt_hours: Time step in hours
        include_crive_fractions: If True, include HD1-3, HP1-3 columns for full RIVE mode
    """
    start = datetime.strptime(start_date, '%Y-%m-%d')
    n_steps = int(duration_days * 24 / dt_hours) + 1
    
    data = {
        'time_s': np.arange(n_steps) * dt_hours * 3600,
    }
    
    for i in range(n_steps):
        t = start + timedelta(hours=i * dt_hours)
        season = get_season(t)
        sf = seasonal_factor(t)
        
        # Salinity (freshwater upstream)
        data.setdefault('salinity', []).append(np.random.uniform(0.0, 0.3))
        
        # Temperature (tropical, 26-31°C)
        temp = 28.5 + 2.5 * sf['temperature'] + np.random.randn() * 0.5
        data.setdefault('temperature', []).append(temp)
        
        # From monitoring data
        mon = MONITORING_DATA
        
        # NH4 and NO3
        nh4_mg = interpolate_seasonal(mon['NH4_mg_N_L']['dry'], mon['NH4_mg_N_L']['wet'], t)
        no3_mg = interpolate_seasonal(mon['NO3_mg_N_L']['dry'], mon['NO3_mg_N_L']['wet'], t)
        data.setdefault('nh4', []).append(mg_N_to_umol(nh4_mg))
        data.setdefault('no3', []).append(mg_N_to_umol(no3_mg))
        
        # PO4 (estimate: PO4 ≈ 0.1 × NH4 in Mekong)
        po4_umol = data['nh4'][-1] * 0.1 + np.random.uniform(0.5, 2.0)
        data.setdefault('po4', []).append(po4_umol)
        
        # O2
        o2_mg = interpolate_seasonal(mon['O2_mg_L']['dry'], mon['O2_mg_L']['wet'], t)
        data.setdefault('o2', []).append(mg_O2_to_umol(o2_mg))
        
        # =======================================================================
        # BOD5/COD → C-RIVE Organic Fractions (The Input Mapper)
        # =======================================================================
        bod5 = interpolate_seasonal(mon['BOD5_mg_L']['dry'], mon['BOD5_mg_L']['wet'], t)
        cod = interpolate_seasonal(mon['COD_mg_L']['dry'], mon['COD_mg_L']['wet'], t)
        spm = interpolate_seasonal(mon['TSS_mg_L']['dry'], mon['TSS_mg_L']['wet'], t)
        
        # Map to C-RIVE fractions
        oc_fractions = bod_cod_to_crive_fractions(bod5, cod, spm)
        
        # Total OC (for simplified mode - diagnostic only, will be sum of pools)
        data.setdefault('toc', []).append(oc_fractions['toc_total'])
        
        # SPM (TSS)
        data.setdefault('spm', []).append(spm)
        
        # From field campaign
        fc = FIELD_CAMPAIGN['dry_season'] if season == 'dry' else FIELD_CAMPAIGN['wet_season']
        
        # DSi (dissolved silica)
        dsi_mg = np.random.uniform(*fc['DSi_mg_L'])
        data.setdefault('dsi', []).append(mg_Si_to_umol(dsi_mg))
        
        # Alkalinity
        alk_meq = np.random.uniform(*fc['alkalinity_meq_L'])
        data.setdefault('at', []).append(meq_to_umol_alk(alk_meq))
        
        # DIC (estimate: DIC ≈ AT for low-salinity freshwater)
        dic = data['at'][-1] * (0.95 + 0.1 * np.random.rand())
        data.setdefault('dic', []).append(dic)
        
        # Phytoplankton (from Chl-a)
        chla = np.random.uniform(*fc['Chla_ug_L'])
        phy = chla_to_phy(chla)
        data.setdefault('phy1', []).append(phy * 0.7)  # Diatoms dominant
        data.setdefault('phy2', []).append(phy * 0.3)
        
        # ===================================================================
        # Diagnostic species (computed by model, but need initial/BC values)
        # ===================================================================
        data.setdefault('pco2', []).append(0.0)   # Computed from carbonate system
        data.setdefault('co2', []).append(0.0)    # Computed from carbonate system
        data.setdefault('ph', []).append(7.2)     # Initial pH estimate
        data.setdefault('hs', []).append(0.0)     # Hydrogen sulfide (anoxic only)
        data.setdefault('alkc', []).append(0.0)   # Iteration counter - diagnostic
        
        # ===================================================================
        # RIVE Multi-pool Organic Matter (from BOD/COD mapping)
        # ===================================================================
        data.setdefault('hd1', []).append(oc_fractions['hd1'])  # Labile dissolved OC
        data.setdefault('hd2', []).append(oc_fractions['hd2'])  # Semi-labile dissolved OC
        data.setdefault('hd3', []).append(oc_fractions['hd3'])  # Refractory dissolved OC
        data.setdefault('hp1', []).append(oc_fractions['hp1'])  # Labile particulate OC
        data.setdefault('hp2', []).append(oc_fractions['hp2'])  # Semi-labile particulate OC
        data.setdefault('hp3', []).append(oc_fractions['hp3'])  # Refractory particulate OC
        
        # ===================================================================
        # RIVE Bacteria (literature values for tropical rivers)
        # Reference: Servais et al. (1999), Garnier et al. (2002)
        # ===================================================================
        data.setdefault('bag', []).append(2.0 + np.random.randn() * 0.3)   # Attached bacteria [mg C/L]
        data.setdefault('bap', []).append(1.0 + np.random.randn() * 0.2)   # Free bacteria [mg C/L]
        
        # ===================================================================
        # RIVE Phosphorus and Substrates
        # ===================================================================
        # PIP: Particulate inorganic P (estimated from SPM and PO4)
        pip = spm * 0.001 * 0.5  # ~0.05% P content in SPM, 50% inorganic
        data.setdefault('pip', []).append(pip)
        
        # DSS: Dissolved simple substrates (labile DOC fraction for bacteria)
        dss = oc_fractions['hd1'] * 0.1 * 12 / 1000  # 10% of HD1 as mg C/L
        data.setdefault('dss', []).append(dss)
        
        # ===================================================================
        # GHG Species (from field campaign data)
        # ===================================================================
        data.setdefault('no2', []).append(0.5 + np.random.rand() * 0.5)  # Nitrite [µmol N/L]
        n2o_nmol = np.random.uniform(*fc['N2O_nmol_L'])
        data.setdefault('n2o', []).append(n2o_nmol)  # N2O [nmol/L]
        ch4_nmol = np.random.uniform(*fc['CH4_nmol_L'])
        data.setdefault('ch4', []).append(ch4_nmol / 1000.0)  # CH4 [µmol/L] (convert from nmol)
    
    return pd.DataFrame(data)


def generate_ocean_boundary(start_date: str, duration_days: int,
                            dt_hours: float = 6.0) -> pd.DataFrame:
    """
    Generate ocean (downstream) boundary conditions.
    
    South China Sea typical values with tidal modulation.
    Generates ALL 30 species required by CGEM_NUM_SPECIES.
    """
    start = datetime.strptime(start_date, '%Y-%m-%d')
    n_steps = int(duration_days * 24 / dt_hours) + 1
    
    data = {
        'time_s': np.arange(n_steps) * dt_hours * 3600,
    }
    
    for i in range(n_steps):
        t = start + timedelta(hours=i * dt_hours)
        sf = seasonal_factor(t)
        
        # Ocean salinity (varies with intrusion)
        sal_base = 32.0 - 5.0 * sf['discharge']  # Lower when discharge high
        data.setdefault('salinity', []).append(sal_base + np.random.randn() * 1.0)
        
        # Temperature (warmer than river due to mixing)
        temp = 29.0 + 1.5 * sf['temperature'] + np.random.randn() * 0.3
        data.setdefault('temperature', []).append(temp)
        
        # Ocean nutrients (low, oligotrophic)
        data.setdefault('no3', []).append(np.random.uniform(1.0, 5.0))
        data.setdefault('nh4', []).append(np.random.uniform(0.5, 2.0))
        data.setdefault('po4', []).append(np.random.uniform(0.1, 0.5))
        data.setdefault('dsi', []).append(np.random.uniform(5.0, 15.0))
        
        # Ocean O2 (well-mixed, near saturation)
        data.setdefault('o2', []).append(np.random.uniform(200.0, 230.0))
        
        # Ocean TOC (lower than river)
        toc_ocean = np.random.uniform(50.0, 100.0)
        data.setdefault('toc', []).append(toc_ocean)
        
        # Ocean SPM (low)
        spm_ocean = np.random.uniform(5.0, 20.0)
        data.setdefault('spm', []).append(spm_ocean)
        
        # Ocean carbonate system
        data.setdefault('dic', []).append(np.random.uniform(2000.0, 2150.0))
        data.setdefault('at', []).append(np.random.uniform(2250.0, 2350.0))
        
        # Ocean phytoplankton (higher near coast)
        chla_ocean = np.random.uniform(1.0, 5.0)
        phy = chla_to_phy(chla_ocean)
        data.setdefault('phy1', []).append(phy * 0.5)
        data.setdefault('phy2', []).append(phy * 0.5)
        
        # ===================================================================
        # Diagnostic species (computed by model)
        # ===================================================================
        data.setdefault('pco2', []).append(0.0)
        data.setdefault('co2', []).append(0.0)
        data.setdefault('ph', []).append(8.1)     # Ocean pH
        data.setdefault('hs', []).append(0.0)
        data.setdefault('alkc', []).append(0.0)
        
        # ===================================================================
        # RIVE Multi-pool Organic Matter (lower in ocean)
        # ===================================================================
        # In ocean, most OC is refractory
        data.setdefault('hd1', []).append(toc_ocean * 0.05)   # 5% labile dissolved
        data.setdefault('hd2', []).append(toc_ocean * 0.15)   # 15% semi-labile dissolved
        data.setdefault('hd3', []).append(toc_ocean * 0.50)   # 50% refractory dissolved
        data.setdefault('hp1', []).append(toc_ocean * 0.02)   # 2% labile particulate
        data.setdefault('hp2', []).append(toc_ocean * 0.08)   # 8% semi-labile particulate
        data.setdefault('hp3', []).append(toc_ocean * 0.20)   # 20% refractory particulate
        
        # ===================================================================
        # RIVE Bacteria (lower in ocean - oligotrophic)
        # ===================================================================
        data.setdefault('bag', []).append(0.5 + np.random.randn() * 0.1)
        data.setdefault('bap', []).append(0.5 + np.random.randn() * 0.1)
        
        # ===================================================================
        # RIVE Phosphorus and Substrates
        # ===================================================================
        data.setdefault('pip', []).append(spm_ocean * 0.0005)  # Lower P content in ocean
        data.setdefault('dss', []).append(0.5 + np.random.rand() * 0.5)
        
        # ===================================================================
        # GHG Species (ocean - near equilibrium with atmosphere)
        # ===================================================================
        data.setdefault('no2', []).append(0.1 + np.random.rand() * 0.1)
        data.setdefault('n2o', []).append(8.0 + np.random.randn() * 1.0)  # Near saturation
        data.setdefault('ch4', []).append(0.005 + np.random.rand() * 0.005)  # Very low
    
    return pd.DataFrame(data)


def generate_discharge_timeseries(start_date: str, duration_days: int,
                                  dt_hours: float = 1.0) -> pd.DataFrame:
    """Generate Mekong discharge using MRC seasonal patterns."""
    start = datetime.strptime(start_date, '%Y-%m-%d')
    n_steps = int(duration_days * 24 / dt_hours) + 1
    
    data = {
        'time_s': np.arange(n_steps) * dt_hours * 3600,
        'Q_m3s': []
    }
    
    Q_ranges = HYDROLOGY['Tien_River']['Q_m3s']
    
    for i in range(n_steps):
        t = start + timedelta(hours=i * dt_hours)
        season = get_season(t)
        sf = seasonal_factor(t)
        
        # Base discharge from season
        if season == 'dry':
            Q_base = np.random.uniform(*Q_ranges['dry'])
        elif season == 'wet':
            Q_base = np.random.uniform(*Q_ranges['wet'])
        else:
            Q_base = np.random.uniform(*Q_ranges['transition'])
        
        # Smooth seasonal transition
        Q_smooth = Q_ranges['dry'][0] + sf['discharge'] * (Q_ranges['wet'][1] - Q_ranges['dry'][0])
        Q = 0.7 * Q_smooth + 0.3 * Q_base
        
        # Add daily variation (dam releases, etc.)
        daily_var = 1.0 + 0.05 * np.sin(2 * np.pi * t.hour / 24)
        
        # Add noise
        noise = 1.0 + 0.08 * np.random.randn()
        
        data['Q_m3s'].append(Q * daily_var * max(0.7, min(1.3, noise)))
    
    return pd.DataFrame(data)


def generate_tide_timeseries(start_date: str, duration_days: int,
                             dt_hours: float = 0.5,
                             branch: str = 'My_Tho') -> pd.DataFrame:
    """Generate M2 tidal elevation for a branch mouth."""
    n_steps = int(duration_days * 24 / dt_hours) + 1
    time_s = np.arange(n_steps) * dt_hours * 3600
    
    # M2 tidal period
    M2_period_s = 12.42 * 3600
    M2_freq = 2 * np.pi / M2_period_s
    
    # Spring-neap modulation
    Msf_period_s = 14.77 * 86400
    Msf_freq = 2 * np.pi / Msf_period_s
    
    # Get branch-specific parameters
    amp_range = HYDROLOGY['Tien_River']['tide_amplitude_m']
    amp_base = np.mean(amp_range)
    damping = HYDROLOGY['branches'].get(branch, {'tide_damping': 1.0})['tide_damping']
    
    H = []
    for ts in time_s:
        # Spring-neap modulation (±20%)
        spring_neap = 1.0 + 0.2 * np.cos(Msf_freq * ts)
        
        # M2 tide with phase shift per branch
        phase_shift = hash(branch) % 360 * np.pi / 180  # Pseudo-random phase
        H_val = amp_base * damping * spring_neap * np.cos(M2_freq * ts + phase_shift)
        H.append(H_val)
    
    return pd.DataFrame({'time_s': time_s, 'H_m': H})


# ============================================================================
# BIOGEO PARAMETERS (MEKONG-CALIBRATED)
# ============================================================================

def generate_biogeo_params_mekong(output_path: str):
    """
    Generate biogeo_params.txt optimized for Mekong Delta.
    
    Uses simplified 80/20 approach with parameters calibrated
    from available monitoring and field campaign data.
    
    CRITICAL: ghg_passive_mode = 1 by default (SAFE mode)
    This prevents mass balance errors from GHG feedback.
    """
    params = """# ============================================================================
# CGEM BIOGEOCHEMISTRY PARAMETERS: MEKONG DELTA (80/20 Mode)
# ============================================================================
# Simplified configuration for data-limited tropical deltas.
# Based on: Monitoring center data + 2024-2025 field campaigns
# 
# Key simplifications:
# - Single TOC pool (no HD1-3, HP1-3 unless skip_multipool_oc=0)
# - No explicit bacteria (bulk kox instead)
# - Single-step nitrification (in core module)
# - GHG in PASSIVE MODE (diagnostic only, no feedback to O2/NH4/NO3)
# - Literature-based benthic fluxes
# ============================================================================

# --- Environmental Parameters ---
water_temp = 28.0           # Mean temperature [°C] - tropical
ws = 0.0008                 # SPM settling velocity [m/s]

# --- Light Parameters (from SPM/Chl-a relationships) ---
I0 = 350.0                  # Surface irradiance [W/m²] - tropical, partly cloudy
kd1 = 0.20                  # Base attenuation [1/m] - higher for turbid delta
kd2_spm = 0.018             # SPM contribution [m²/mg] - calibrated to field Kd
kd2_phy1 = 0.008            # Phytoplankton self-shading
kd2_phy2 = 0.006

# --- Phytoplankton (calibrated to satellite Chl-a) ---
# Using single PHY for simplicity (PHY1=diatoms, PHY2 minimal)
alpha1 = 0.02               # P-I curve slope
pbmax1 = 2.8                # Max photosynthesis [1/day] - tropical, high
kexc1 = 0.10                # Excretion fraction
kgrowth1 = 0.10             # Growth respiration
kmaint1 = 0.03              # Maintenance [1/day]
kmort1 = 0.08               # Mortality [1/day] - higher in tropics

alpha2 = 0.015
pbmax2 = 2.0
kexc2 = 0.10
kgrowth2 = 0.10
kmaint2 = 0.02
kmort2 = 0.06

# --- Nutrient Half-Saturation (literature + field adjustment) ---
kdsi1 = 3.0                 # DSi [µmol/L] - Mekong has high Si
kn1 = 1.5                   # N [µmol/L]
kpo41 = 0.3                 # PO4 [µmol/L]
kn2 = 2.0
kpo42 = 0.5

# --- Decomposition Rates (from BOD5/COD data) ---
# kox calibrated to: kox ≈ BOD5 / (TOC × 5 days)
# Mekong BOD5 ~5 mg/L, TOC ~250 µmol/L → kox ~ 0.15/day
kox = 0.15                  # Aerobic TOC decay [1/day]
kdenit = 0.04               # Denitrification [1/day] - limited by O2>30 µM
knit = 0.12                 # Nitrification [1/day] - calibrated to NH4 decay

# --- Half-Saturation Constants ---
ktox = 40.0                 # TOC half-sat [µmol/L]
ko2 = 8.0                   # O2 for respiration [µmol/L]
ko2_nit = 4.0               # O2 for nitrification [µmol/L]
kno3 = 4.0                  # NO3 for denitrification [µmol/L]
knh4 = 2.0                  # NH4 for nitrification [µmol/L]
kino2 = 6.0                 # O2 inhibition for denitrification [µmol/L]

# --- Stoichiometry (Redfield) ---
redn = 0.151                # N:C ratio
redp = 0.0094               # P:C ratio
redsi = 0.12                # Si:C for diatoms (Mekong high Si)

# --- Gas Exchange (from field pCO2 + literature k600) ---
pco2_atm = 420.0            # Atmospheric pCO2 [µatm]
wind_speed = 3.0            # Mean wind [m/s] - sheltered delta
wind_coeff = 0.251          # Wanninkhof coefficient
schmidt_exp = -0.5
current_k_factor = 0.30     # Current contribution (important in tidal channels)

# --- Benthic Fluxes (CRITICAL - literature for tropical deltas) ---
# Reference: Borges et al. (2015), Lara et al. (2017)
# Higher than temperate due to temperature and organic loading
benthic_resp_20C = 70.0     # Benthic CO2 flux [mmol C/m²/day]
benthic_Q10 = 2.0           # Temperature coefficient
benthic_NH4_flux = 2.5      # Benthic NH4 flux [mmol N/m²/day]
benthic_PO4_flux = 0.15     # Benthic PO4 flux [mmol P/m²/day]

# ============================================================================
# GREENHOUSE GAS PARAMETERS (Passive Mode - Diagnostic Only)
# ============================================================================
# Reference: Garnier et al. (2007), Marescaux et al. (2019)
# 
# In PASSIVE MODE (ghg_passive_mode=1), the GHG module:
# - Reads nitrification/denitrification rates from core module
# - Calculates N2O as YIELD FRACTION of these rates (no recalculation)
# - CH4 is from benthic flux only (does not feedback to O2)
# - Does NOT modify O2, NH4, NO3 (prevents double-counting)
#
# CALIBRATION TIP: If modeled N2O is too low vs field data:
#   → Increase N2O_yield_nit (typical range: 0.003-0.01)
#   → DO NOT change the code, change the yield parameter!
# ============================================================================

# N2O yield fractions (fraction of N processed that becomes N2O)
N2O_yield_nit = 0.005       # N2O yield from nitrification (0.3-1%)
N2O_yield_denit = 0.015     # N2O yield from denitrification (0.5-3%)

# Benthic GHG fluxes (for passive mode - independent of water column)
benthic_CH4_flux = 150.0    # CH4 sediment flux [µmol/m²/day] - tropical, high
benthic_N2O_flux = 8.0      # N2O sediment flux [nmol/m²/day]

# ============================================================================
# SAFETY MODE FLAGS (80/20 Simplification)
# ============================================================================
# These flags control which processes are active. 
# Default settings are SAFE for data-limited applications.
#
# CRITICAL: ghg_passive_mode = 1 prevents:
#   - Double counting of nitrification (core vs GHG module)
#   - Mass balance errors from uncalibrated CH4 oxidation
#   - O2 crashing to negative values
# ============================================================================

ghg_passive_mode = 1        # 1 = SAFE (GHG is diagnostic only)
                            # 0 = ACTIVE (GHG feeds back to O2/NH4/NO3)
                            # WARNING: Only set to 0 if you have calibrated GHG params!

simplified_mode = 1         # 1 = Enable 80/20 simplified mode
skip_bacteria = 1           # 1 = Use bulk kox instead of explicit BAG/BAP bacteria
skip_multipool_oc = 1       # 1 = Use single TOC instead of HD1-3, HP1-3
skip_ghg_dynamics = 0       # 1 = Skip GHG entirely (fastest, least output)
                            # 0 = Calculate GHG in passive diagnostic mode
skip_p_adsorption = 1       # 1 = Skip dynamic PO4-PIP equilibrium

# ============================================================================
# VALIDATION TARGETS (from field campaign 2024-2025)
# Use these to check model performance
# ============================================================================
# Salinity intrusion (dry season): 40-60 km from mouth
# pCO2 range: 1000-4000 µatm (supersaturated)
# O2 minimum: 150-200 µmol/L in turbidity maximum
# Chl-a maximum: 5-15 µg/L at salinity 5-15 PSU
# N2O: 8-35 nmol/L (field range), saturation 100-300%
# CH4: 50-800 nmol/L (field range), highest in wet season
"""
    
    with open(output_path, 'w') as f:
        f.write(params)
    print(f"Generated: {output_path}")


# ============================================================================
# MAIN CASE GENERATOR
# ============================================================================

def generate_tien_river_case(case_dir: str, start_date: str = '2024-01-01',
                              duration_days: int = 365):
    """
    Generate complete input files for Tien River case.
    
    Creates:
    - case_config.txt
    - topology.csv
    - boundary_map.csv
    - biogeo_params.txt
    - forcing_data/*.csv
    """
    case_path = Path(case_dir)
    case_path.mkdir(parents=True, exist_ok=True)
    forcing_path = case_path / 'forcing_data'
    forcing_path.mkdir(exist_ok=True)
    
    print(f"\n{'='*60}")
    print(f"Generating Tien River Case (Mekong 80/20 Mode)")
    print(f"{'='*60}")
    print(f"Start: {start_date}, Duration: {duration_days} days")
    print(f"Output: {case_dir}")
    print()
    
    # 1. Discharge (MRC-based)
    print("1. Generating discharge time series...")
    Q_df = generate_discharge_timeseries(start_date, duration_days)
    Q_df.to_csv(forcing_path / 'Tien_Input.csv', index=False)
    
    # 2. Tidal levels for each branch
    print("2. Generating tidal boundary conditions...")
    for branch in ['My_Tho', 'Ham_Luong', 'Co_Chien']:
        tide_df = generate_tide_timeseries(start_date, duration_days, branch=branch)
        tide_df.to_csv(forcing_path / f'{branch}_Tide.csv', index=False)
    
    # 3. River (upstream) species boundary
    print("3. Generating river boundary conditions...")
    river_bc = generate_river_boundary(start_date, duration_days)
    river_bc.to_csv(forcing_path / 'species_river.csv', index=False)
    
    # 4. Ocean (downstream) species boundary
    print("4. Generating ocean boundary conditions...")
    ocean_bc = generate_ocean_boundary(start_date, duration_days)
    ocean_bc.to_csv(forcing_path / 'species_ocean.csv', index=False)
    
    # 5. Biogeochemistry parameters
    print("5. Generating biogeochemistry parameters...")
    generate_biogeo_params_mekong(case_path / 'biogeo_params.txt')
    
    print(f"\n{'='*60}")
    print("Input generation complete!")
    print(f"{'='*60}")
    print(f"\nFiles created:")
    print(f"  - forcing_data/Tien_Input.csv (discharge)")
    print(f"  - forcing_data/My_Tho_Tide.csv")
    print(f"  - forcing_data/Ham_Luong_Tide.csv")
    print(f"  - forcing_data/Co_Chien_Tide.csv")
    print(f"  - forcing_data/species_river.csv")
    print(f"  - forcing_data/species_ocean.csv")
    print(f"  - biogeo_params.txt")
    print(f"\nNote: Uses 80/20 simplified mode for Mekong Delta.")
    print(f"See docs/MEKONG_80_20_ANALYSIS.md for calibration guidance.")


def main():
    parser = argparse.ArgumentParser(
        description='Generate C-GEM inputs for Mekong Delta (80/20 simplified mode)'
    )
    parser.add_argument('--case', '-c', default='Tien_River',
                        help='Case name (default: Tien_River)')
    parser.add_argument('--output', '-o', default=None,
                        help='Output directory')
    parser.add_argument('--start', '-s', default='2017-01-01',
                        help='Start date YYYY-MM-DD')
    parser.add_argument('--days', '-d', type=int, default=365,
                        help='Duration in days')
    parser.add_argument('--seed', type=int, default=42,
                        help='Random seed')
    
    args = parser.parse_args()
    np.random.seed(args.seed)
    
    # Determine output path
    if args.output:
        case_dir = args.output
    else:
        script_dir = Path(__file__).parent
        project_root = script_dir.parent
        case_dir = project_root / 'INPUT' / 'Cases' / args.case
    
    generate_tien_river_case(
        case_dir=str(case_dir),
        start_date=args.start,
        duration_days=args.days
    )


if __name__ == '__main__':
    main()
