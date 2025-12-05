# Comprehensive C-GEM Model Audit Report v1.2

**Date:** December 5, 2025 (Final Audit)  
**Auditor:** AI-assisted technical review  
**Model Version:** CGEM_RIVE branch  
**Test Case:** Mekong Delta, Vietnam (March 2025 dry season field data)  
**Status:** Model limitations identified - requires fundamental restructuring for some variables

---

## Executive Summary

After extensive debugging of the salinity-based benthic scaling and lateral load corrections, the model shows **mixed performance**:

- **Excellent**: Salinity (R² = 0.65-0.99), pH (R² = 0.61-0.84)
- **Good**: NO₃ (R² = 0.95 for Co_Chien), O₂ (significantly improved), Alkalinity
- **Poor**: TOC (systematic underestimation), pCO₂ (underestimated), CH₄/N₂O (missing sources)

The **fundamental limitation** is that TOC and O₂ cannot be simultaneously fit with the current model structure because:
1. O₂ requires minimal benthic consumption to match the conservative mixing gradient
2. pCO₂ requires significant respiration (benthic + water column) to generate the supersaturation
3. TOC requires minimal consumption to maintain observed concentrations

These competing requirements cannot be satisfied with a single set of parameters.

---

## 1. Validation Results by Zone (Co_Chien Branch)

### 1.1 Zonal Analysis

| Zone | Distance | Salinity | O₂ Obs | O₂ Model | TOC Obs | TOC Model | Status |
|------|----------|----------|--------|----------|---------|-----------|--------|
| **Downstream** | 0-10 km | 25-30 PSU | 250-260 | 222-258 | 160-170 | 99-148 | O₂ improved, TOC too low |
| **Transition** | 10-30 km | 5-25 PSU | 220-250 | 213-223 | 135-160 | 90-100 | O₂ good, TOC too low |
| **Upstream** | 60-85 km | 0-2 PSU | 170-190 | 180-182 | 130-145 | 105-109 | O₂ slightly high, TOC low |

### 1.2 Key Metrics (December 5, 2025 - Final)

| Variable | Branch | RMSE | Bias | R² | Assessment |
|----------|--------|------|------|-----|------------|
| **Salinity** | Co_Chien | 1.38 | -0.77 | **0.987** | ✅ Excellent |
| **Salinity** | My_Tho | 4.98 | +1.25 | **0.912** | ✅ Excellent |
| **O₂** | Co_Chien | 29.03 | -10.61 | **0.597** | ✅ Good (improved from R²=0.30) |
| **O₂** | My_Tho | 18.31 | -13.51 | **0.804** | ✅ Very Good |
| **pH** | Co_Chien | 0.15 | +0.07 | **0.810** | ✅ Excellent |
| **pH** | My_Tho | 0.12 | +0.06 | **0.841** | ✅ Excellent |
| **NO₃** | Co_Chien | 4.75 | -1.34 | **0.950** | ✅ Excellent |
| **TOC** | Co_Chien | 50.82 | **-46.69** | 0.230 | 🔴 Model too low everywhere |
| **TOC** | Hau_River | 27.41 | -5.98 | **0.802** | ✅ Best branch for TOC |
| **pCO₂** | Co_Chien | 1461 | **-841** | 0.549 | 🔴 Model too low upstream |
| **pCO₂** | Hau_River | 1052 | -92.77 | **0.746** | ✅ Best branch for pCO₂ |
| **CH₄** | All | 24-55 | -10 to -22 | 0.20-0.38 | ⚠️ Missing lateral sources |
| **N₂O** | All | 8-20 | -6 to -14 | 0.27-0.46 | ⚠️ Missing agricultural inputs |

---

## 2. Root Causes of Remaining Biases

### 2.1 TOC Underestimation (-47 µmol/L bias in Co_Chien)

**Problem**: Model TOC is systematically lower than observed across all zones.

**Observed Pattern**:
- Downstream (ocean): ~160-170 µmol/L
- Upstream (river): ~130-145 µmol/L
- Gradient: DECREASING from ocean to river

**Model Pattern** (after fixes):
- Downstream: ~100-148 µmol/L (correct trend)
- Upstream: ~105-109 µmol/L
- Gap: 50-60 µmol/L too low

**Root Causes**:
1. **Ocean boundary TOC too low**: Set at 148 µmol/L but observed is 160-170
2. **Water column degradation (kox)**: Even at 0.012/day, still consuming too much
3. **Lateral TOC inputs were over-corrected**: Fixed from 30,000 → 400 µmol/L, may be too low now

**Potential Fixes** (not implemented):
- Increase ocean boundary TOC to 165 µmol/L
- Further reduce kox to 0.005/day (but this will break pCO₃)
- Increase lateral TOC to ~800-1000 µmol/L (balanced addition)

### 2.2 pCO₂ Underestimation (-841 µatm bias in Co_Chien)

**Problem**: Model pCO₂ is too low, especially in the transition and upstream zones.

**Observed Pattern**:
- Downstream: ~800-1500 µatm
- Mid-estuary: ~2000-3000 µatm
- Upstream: ~3500-4500 µatm

**Model Pattern**:
- Downstream: ~500-900 µatm (reasonable)
- Mid-estuary: ~1200-1500 µatm (too low)
- Upstream: ~1800-2200 µatm (too low)

**Root Causes**:
1. **Insufficient benthic respiration**: Reduced to fix O₂, but this also reduced CO₂ production
2. **Low water column respiration**: kox reduced means less CO₂ from TOC degradation
3. **Trade-off with O₂**: Cannot increase respiration without over-consuming O₂

**The Fundamental Conflict**:
```
O₂ balance:    O₂ loss = benthic SOD + water column respiration
CO₂ balance:   CO₂ gain = benthic respiration + water column respiration

If SOD is high:  O₂ drops too fast, pCO₂ correct
If SOD is low:   O₂ correct, pCO₂ too low
```

### 2.3 CH₄ and N₂O Underestimation

**Problem**: Greenhouse gas concentrations are systematically too low upstream.

**Observed**:
- CH₄: 50-200 nmol/L upstream, 5-30 nmol/L at mouth
- N₂O: 10-25 nmol/L upstream, 5-10 nmol/L at mouth

**Model**:
- CH₄: 20-50 nmol/L (missing major sources)
- N₂O: 5-10 nmol/L (missing agricultural inputs)

**Root Causes**:
1. **Missing rice paddy CH₄ flux**: Rice paddies are THE major CH₄ source in Mekong
2. **Missing aquaculture pond inputs**: Significant CH₄ and NH₄ sources
3. **Lateral load system inadequate**: Current system distributes loads evenly, not spatially-targeted

---

## 3. Model Architecture Limitations

### 3.1 Single Benthic Flux Rate

**Current Design**: One benthic_resp_20C parameter for entire network.

**Limitation**: Cannot have:
- Low SOD in downstream (for O₂ preservation)
- High SOD in upstream (for pCO₂ generation)

**Proposed Fix**: Branch-specific or segment-specific benthic parameters in topology.csv

### 3.2 No Decoupling of O₂ and CO₂ Benthic Fluxes

**Current Design**: Benthic respiration consumes O₂ AND produces CO₂ with 1:1 stoichiometry.

**Reality**: Estuarine sediments can have:
- Anaerobic CO₂ production (no O₂ consumption)
- CO₂ from carbonate dissolution (no O₂ involved)
- Sulfate reduction producing CO₂ without O₂

**Proposed Fix**: Separate benthic_o2_flux and benthic_co2_flux parameters.

### 3.3 Lateral Load Spatial Resolution

**Current Design**: Lateral loads per segment, but TOC/nutrient concentrations are uniform.

**Limitation**: Cannot represent:
- Point sources (cities) vs diffuse sources (agriculture)
- Tidal-influenced lateral inputs
- Sub-daily variation

**Proposed Fix**: Implement sub-segment lateral load targeting.

---

## 4. Recommended Priority Fixes

### Priority 1: Decouple O₂ and CO₂ Benthic Fluxes (Impact: High)

```c
// New parameters in biogeo_params.txt:
benthic_o2_rate = 30.0    // mmol O₂/m²/day (lower for O₂ fit)
benthic_co2_rate = 80.0   // mmol CO₂/m²/day (higher for pCO₂ fit)
```

This allows fitting O₂ and pCO₂ independently.

### Priority 2: Increase Ocean TOC Boundary (Impact: Medium)

```csv
# In species_ocean_realistic.csv:
# Change TOC from 140-160 to 165-180 µmol/L
```

### Priority 3: Branch-Specific Benthic Parameters (Impact: Medium)

```csv
# In topology.csv, add column:
# Benthic_Scale
# Co_Chien: 0.5 (lower SOD, more marine)
# Tien_Main: 1.5 (higher SOD, more organic sediment)
```

### Priority 4: Targeted CH₄ Sources (Impact: High for GHG)

Create spatial CH₄ source map based on:
- Rice paddy area (JAXA land use)
- Aquaculture pond density
- Distance from major tributaries

---

## 5. What the Model DOES Capture Well

Despite limitations, the model successfully reproduces:

1. **Salinity intrusion length and gradient** (R² > 0.9)
2. **Tidal variations** in water level and velocity
3. **NO₃ gradient** (conservative mixing + some reactions)
4. **pH gradient** (controlled by carbonate equilibrium)
5. **Alkalinity spatial pattern** (mixing-dominated)
6. **O₂ TREND** (decreasing upstream, even if absolute values off)
7. **pCO₂ TREND** (increasing upstream, gradient correct)

---

## 6. Conclusion

The C-GEM model in its current form is **suitable for**:
- Salinity intrusion studies
- Nutrient transport modeling (NO₃, PO₄)
- pH and carbonate chemistry (first-order)
- Relative spatial comparisons

The model is **NOT YET suitable for**:
- Quantitative O₂/pCO₂ predictions (requires decoupled benthic fluxes)
- TOC budget calculations (requires better boundary conditions)
- CH₄/N₂O emission inventories (requires spatially-explicit sources)

**Recommended Next Steps**:
1. Implement decoupled benthic O₂/CO₂ fluxes
2. Add branch-specific benthic parameters
3. Develop spatially-explicit CH₄ source module for rice paddies
4. Acquire better boundary condition data (field measurements)

---

## Appendix: Parameter Values Used in This Audit

```
# Core degradation rates
kox = 0.012 /day (water column TOC oxidation)
knit = 0.05 /day (nitrification)
kdenit = 0.03 /day (denitrification)

# Benthic fluxes
benthic_resp_20C = 45.0 mmol C/m²/day
benthic_Q10 = 1.8
benthic_ocean_scale = 0.2 (S > 2 PSU)
benthic_upstream_scale = 2.0 (S ≤ 2 PSU)

# Gas exchange
wind_speed = 5.0 m/s
current_k_factor = 1.0
pco2_atm = 420 ppm

# Lateral loads (after correction)
TOC_lateral = 4-6 mg C/L (≈330-500 µmol/L)
CH4_lateral = 0.03-0.25 µmol/L (30-250 nmol/L)
N2O_lateral = 0.01-0.08 µmol/L (10-80 nmol/L)
```

---

*End of Audit Report*
