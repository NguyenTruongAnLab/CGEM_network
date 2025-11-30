# C-GEM Mekong Delta: 80/20 Simplification Analysis

## Executive Summary

This document analyzes the C-GEM Network model's data requirements for **realistic application in the Mekong Delta**, where monitoring data is limited compared to European rivers like the Seine or Loire. We apply the **80/20 rule** (Pareto principle) to identify which parameters and processes contribute most to model accuracy.

---

## 1. Available Data Assessment

### What You Have

| Data Source | Variables | Coverage | Resolution | Quality |
|-------------|-----------|----------|------------|---------|
| **Monitoring Center** | BOD5, COD, NH4, NO3, pH, O2, TSS | Multi-year | Monthly/quarterly | ⭐⭐⭐ |
| **Field Campaign 2024-25** | POC, DOC, TOC, Chl-a, salinity, alkalinity, TSS, nutrients, CO2, CH4, N2O | 2 snapshots | All branches, upstream→downstream | ⭐⭐⭐⭐ |
| **Sentinel Satellite** | Chl-a, SPM, DOC, POC, CDOM | 2016-2025 | Weekly (>300m), Monthly (<300m) | ⭐⭐⭐ |
| **Hydrology** | Discharge (MRC), tidal levels | Multi-year | Daily/hourly | ⭐⭐⭐⭐ |

### What You're Missing

| Critical Gap | Impact | Workaround |
|--------------|--------|------------|
| **Benthic fluxes** | HIGH - affects 30-50% of nutrient/CO2 budget | Use literature values from similar tropical deltas |
| **Continuous species BC** | MEDIUM - boundaries drive gradients | Interpolate from snapshots + satellite |
| **Light attenuation (Kd)** | MEDIUM - controls phytoplankton | Calculate from SPM/Chl-a |
| **Bacterial biomass** | LOW - internal pools equilibrate | Remove explicit bacteria, use bulk TOC decay |
| **Organic matter pools (HD1-3, HP1-3)** | LOW - over-parameterized for available data | Aggregate to TOC only |

---

## 2. Process Importance Ranking

### HIGH Priority (80% of accuracy) - KEEP

| Process | Impact on Key Outputs | Data Required | You Have It? |
|---------|----------------------|---------------|--------------|
| **Salt intrusion** | Controls ecosystem gradients | Salinity BC | ✅ Field data |
| **O2 dynamics** | Hypoxia prediction | O2 BC, reaeration | ✅ Monitoring |
| **SPM transport** | Light, adsorption | TSS BC | ✅ Monitoring + satellite |
| **Phytoplankton** | Chl-a patterns, NPP | Chl-a, nutrients | ✅ Satellite + field |
| **Basic N cycle** | NH4→NO3 | NH4, NO3 | ✅ Monitoring |
| **CO2 exchange** | GHG flux | pCO2, alkalinity | ✅ Field campaign |

### MEDIUM Priority (15% of accuracy) - SIMPLIFY

| Process | Current Complexity | Simplified Approach |
|---------|-------------------|---------------------|
| **6-pool organic matter** | HD1, HD2, HD3, HP1, HP2, HP3 | → Single TOC pool |
| **Explicit bacteria** | BAG, BAP with growth | → Bulk decay rate (kox) |
| **P adsorption** | Dynamic PO4 ↔ PIP | → Constant PO4 ratio |
| **2-step nitrification** | NH4→NO2→NO3 | → Direct NH4→NO3 |

### LOW Priority (5% of accuracy) - REMOVE/HARDCODE

| Process | Why Remove | Default Value |
|---------|-----------|---------------|
| **N2O dynamics** | Needs continuous monitoring | Calculate post-hoc from rates |
| **CH4 ebullition** | Needs benthic measurements | Literature flux estimate |
| **Sulfide (HS)** | Rare in Mekong main channel | Set to 0 |
| **DSS (dissolved substrates)** | Internal bacteria food | Remove |

---

## 3. Recommended "Mekong Mode" Configuration

### Simplified Species List (12 instead of 30)

```
ESSENTIAL (transport + reactions):
  1. Salinity      - Conservative tracer, mixing indicator
  2. PHY (total)   - Chlorophyll proxy (combine PHY1+PHY2)
  3. DSi           - Diatom limitation (important in Mekong)
  4. NO3           - Nitrate
  5. NH4           - Ammonium
  6. PO4           - Phosphate
  7. O2            - Dissolved oxygen
  8. TOC           - Total organic carbon (aggregate)
  9. SPM           - Suspended sediment
  10. DIC          - Dissolved inorganic carbon
  11. AT           - Total alkalinity

DIAGNOSTIC (calculated, not transported):
  12. pCO2, pH     - From DIC/AT equilibrium
```

### Simplified Reaction List (8 instead of 53)

```
ESSENTIAL REACTIONS:
  1. NPP           - Net primary production (from Chl-a/light)
  2. AER_DEG       - Aerobic degradation (TOC mineralization)
  3. NIT           - Nitrification (NH4→NO3, single step)
  4. DENIT         - Denitrification (NO3→N2)
  5. O2_EX         - O2 air-water exchange
  6. CO2_EX        - CO2 air-water exchange
  7. EROSION       - SPM erosion
  8. DEPOSITION    - SPM deposition
```

### Key Parameter Sources

| Parameter | Source | Approach |
|-----------|--------|----------|
| **kox** (TOC decay) | BOD5/COD data | kox ≈ BOD5/(TOC × 5 days) |
| **knit** | NH4 monitoring | Calibrate to NH4 profiles |
| **pbmax** | Chl-a + light | Calibrate to satellite Chl-a |
| **Kd** (light atten.) | SPM data | Kd = 0.15 + 0.015 × SPM |
| **k600** | Wind + current | Raymond or Abril formula |
| **benthic_resp** | Literature | 40-80 mmol C/m²/day (tropical) |

---

## 4. Data Conversion Formulas

### From Monitoring Data to Model Units

```python
# BOD5 → TOC approximation
TOC_umol = BOD5_mg_L * 1000 / 12 / 1.5  # Assume BOD:TOC ≈ 1.5

# COD → TOC approximation  
TOC_umol = COD_mg_L * 1000 / 12 / 2.67  # COD:TOC ≈ 2.67

# Chl-a → Phytoplankton carbon
PHY_ugC_L = Chla_ug_L * 50  # C:Chl ratio = 50 (tropical)

# TSS → SPM
SPM_mg_L = TSS_mg_L  # Direct conversion

# mg/L → µmol/L conversions
NO3_umol = NO3_mg_N_L * 1000 / 14
NH4_umol = NH4_mg_N_L * 1000 / 14
PO4_umol = PO4_mg_P_L * 1000 / 31
O2_umol = O2_mg_L * 1000 / 32
```

### Satellite Data Integration

```python
# Sentinel-3 OLCI products → Model inputs
Chl_a_model = Chl_a_satellite  # Direct use [µg/L]
SPM_model = TSM_satellite      # Total suspended matter [mg/L]

# DOC from CDOM (regional relationship for Mekong)
# Based on tropical river studies: DOC ≈ 50 + 0.8 × CDOM_aCDOM(443)
DOC_umol = (50 + 0.8 * aCDOM_443) * 1000 / 12

# Weekly satellite → Daily interpolation
# Use cubic spline for gap filling
```

---

## 5. Validation Strategy

### Using Your 2024-2025 Field Campaign

**Dry Season Snapshot:**
- Compare model salinity intrusion length
- Validate pCO2 spatial pattern
- Check O2 minimum zones

**Wet Season Snapshot:**
- Validate nutrient flushing
- Check SPM transport
- Compare Chl-a patterns

### Using Satellite Time Series

**Monthly Validation:**
- Chl-a mean per branch
- SPM seasonal pattern
- DOC/CDOM relationships

### Performance Targets

| Variable | Metric | Target |
|----------|--------|--------|
| Salinity | RMSE | < 3 PSU |
| Chl-a | Bias | < 50% |
| O2 | RMSE | < 30 µmol/L |
| pCO2 | Bias | < 200 µatm |
| SPM | Log-RMSE | < 0.5 |

---

## 6. Implementation Recommendations

### Phase 1: Core Model (Now)

1. **Enable "Mekong Mode"** in config:
   ```
   SimplifiedMode = 1
   NumSpecies = 12
   SkipBacteria = 1
   SkipGHG = 1
   ```

2. **Use monitoring data for BCs**:
   - Monthly interpolation for NH4, NO3, O2
   - Seasonal mean for TOC (from BOD5/COD)

3. **Calibrate key parameters**:
   - kox, knit, pbmax using dry season snapshot

### Phase 2: Enhanced Model (6 months)

1. **Integrate satellite data**:
   - Weekly Chl-a, SPM boundary conditions
   - Spatial validation maps

2. **Add benthic flux**:
   - Literature-based SOD estimate
   - Scaled by depth and organic content

### Phase 3: Full Model (1 year+)

1. **Full RIVE activation** when/if:
   - Benthic chamber measurements available
   - Bacterial biomass data collected
   - Long-term monitoring established

---

## 7. Literature Values for Missing Data

### Benthic Fluxes (Tropical Deltas)

| Flux | Range | Recommended | Reference |
|------|-------|-------------|-----------|
| SOD | 30-100 mmol O2/m²/d | 50 | Lara et al. (2017) Amazon |
| Benthic CO2 | 40-120 mmol C/m²/d | 60 | Borges et al. (2015) tropics |
| NH4 flux | 0.5-5 mmol N/m²/d | 2 | Billen et al. (1994) |
| PO4 flux | 0.05-0.5 mmol P/m²/d | 0.15 | Various |

### Rate Constants (Tropical Estuaries)

| Parameter | Range | Mekong | Reference |
|-----------|-------|--------|-----------|
| kox | 0.05-0.3 /day | 0.15 | Volta et al. (2016) |
| knit | 0.05-0.2 /day | 0.10 | Billen (1994) |
| kdenit | 0.02-0.1 /day | 0.05 | Seitzinger (1988) |
| pbmax | 1.5-4.0 /day | 2.5 | Garnier (2002) |

---

## 8. Summary: 80/20 Action Items

### DO (High Impact, Feasible)

✅ Use salinity, O2, NH4, NO3 from monitoring  
✅ Use Chl-a, SPM from satellite  
✅ Use CO2, alkalinity from field campaign  
✅ Validate with spatial snapshots  
✅ Literature-based benthic flux  

### SIMPLIFY (Medium Impact)

⚡ Aggregate organic matter to single TOC  
⚡ Remove explicit bacteria  
⚡ Single-step nitrification  
⚡ Skip GHG species for now  

### SKIP (Low Impact, High Data Need)

❌ Multi-pool organic matter (needs DOC fractionation)  
❌ 2-step nitrification (needs NO2 data)  
❌ CH4 ebullition (needs benthic probes)  
❌ N2O dynamics (needs continuous monitoring)  

---

## References

- Billen et al. (1994) RIVE model - Hydrobiologia
- Borges et al. (2015) CO2 in tropical rivers - Biogeosciences
- Garnier et al. (2002) Nutrient dynamics - Hydrobiologia
- Lara et al. (2017) Amazon carbon budget - Nature Geoscience
- Volta et al. (2016) C-GEM Scheldt - GMD
- Wang et al. (2018) C-RIVE sensitivity - Water Research
