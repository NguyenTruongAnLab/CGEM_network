# Model Validation and Known Issues

This document summarizes the validation results against Mekong Delta field data (March 2025) and identifies known model limitations requiring future development.

---

## Validation Dataset

**Source**: Mekong Delta field campaign, March 2025 (dry season)

**Branches sampled**: Hau River, Ham Luong, Co Chien, My Tho

**Variables**: Salinity, O₂, pCO₂, pH, NO₃, NH₄, PO₄, TOC, SPM, AT, CH₄, N₂O

---

## Validation Summary (v1.2.0)

### Reactions OFF (Transport Only)

| Variable | Best Branch | R² | Worst Branch | R² | Assessment |
|----------|-------------|-----|--------------|-----|------------|
| **Salinity** | Co_Chien | 0.972 | Ham_Luong | 0.613 | ✅ Good |
| **O₂** | My_Tho | 0.499 | Hau_River | 0.186 | ⚠️ Bias -19 to -36 µM |
| **pCO₂** | Hau_River | 0.968 | Ham_Luong | 0.806 | ✅ Good |
| **pH** | My_Tho | 0.974 | Co_Chien | 0.888 | ✅ Good |
| **NO₃** | Co_Chien | 0.767 | My_Tho | 0.005 | ⚠️ Mixed |
| **NH₄** | - | <0.2 | - | - | ❌ Poor |
| **TOC** | Hau_River | 0.568 | Co_Chien | 0.006 | ⚠️ Bias +13 to -60 µM |
| **SPM** | Hau_River | 0.575 | Ham_Luong | 0.015 | ⚠️ Mixed |
| **AT** | Co_Chien | 0.938 | Ham_Luong | 0.533 | ✅ Reasonable |
| **CH₄** | Ham_Luong | 0.372 | Co_Chien | 0.220 | ❌ Bias -52 to -88 nmol/L |
| **N₂O** | My_Tho | 0.574 | Hau_River | 0.030 | ⚠️ Mixed |

### Reactions ON (Full Biogeochemistry)

| Variable | Before → After | Assessment |
|----------|----------------|------------|
| **O₂** | Bias -19→0 µM | ✅ Much improved |
| **pCO₂** | R² 0.97→0.001 | ❌ Much worse (too high) |
| **pH** | R² 0.94→0.08 | ❌ Much worse (too low) |
| **NO₃** | Similar | - |
| **NH₄** | Bias +0.25→-0.73 | Changed direction |
| **TOC** | Bias +14→+3 | ✅ Slightly improved |

---

## Known Issues by Variable

### 1. O₂ (Dissolved Oxygen)

**Observed Pattern**: Decreases from 260 µM at mouth to 150-180 µM upstream over 80 km

**Model Behavior**: 
- Reactions OFF: Stays ~220 µM (transport only preserves ocean value)
- Reactions ON: Drops appropriately due to respiration

**Root Cause**: With reactions OFF, O₂ behaves conservatively. Real O₂ consumption requires:
1. TOC lateral inputs from agriculture
2. Active respiration (reactions ON)

**Recommendation**: Run with reactions ON for O₂ validation

---

### 2. pCO₂ / pH

**Observed Pattern**: 
- pCO₂ increases from 400-800 µatm at mouth to 3000-6000 µatm upstream
- pH decreases from 8.0-8.1 at mouth to 7.2-7.5 upstream

**Model Behavior**:
- Reactions OFF: Good match (R² > 0.8, transport-dominated)
- Reactions ON: pCO₂ too high by 1000-2000 µatm, pH too low by 0.3

**Root Cause**: Respiration rate overestimated, producing too much CO₂

**Potential Fixes**:
1. Reduce `k_resp` parameter in biogeo_params.txt
2. Increase photosynthesis to offset respiration
3. Add alkalinity inputs from weathering

---

### 3. TOC (Total Organic Carbon)

**Observed Pattern**: 140-200 µM throughout, slight increase upstream

**Model Behavior**: 80-120 µM (too low by 40-80 µM)

**Root Cause**: Missing lateral TOC inputs from:
- Agricultural runoff (rice paddies release organic matter)
- Urban/industrial effluent
- Aquaculture pond discharge
- Mangrove leaf litter

**Required Fix**: Add TOC to lateral sources with appropriate concentrations:
```csv
# In lateral_sources.csv, add TOC column
TOC_conc_base_uM = 200  # for agricultural zones
TOC_conc_base_uM = 300  # for urban zones
```

---

### 4. CH₄ (Methane)

**Observed Pattern**: 
- 40-80 nmol/L at mouth (ocean value)
- 100-200 nmol/L upstream (production zones)

**Model Behavior**: 
- v1.2.0: 40 nmol/L at mouth (fixed!)
- 1-10 nmol/L mid-estuary (too low)
- Near 0 upstream (much too low)

**Root Cause**: 
1. ✅ Ocean BC bug fixed in v1.2.0
2. ❌ Missing lateral CH₄ from rice paddies (major source)
3. ❌ Missing aquaculture pond CH₄
4. ❌ Benthic flux may be too low

**Required Fixes**:
1. Add CH₄ to lateral sources:
   ```csv
   CH4_conc_base_nM = 500  # for rice paddy zones
   CH4_conc_base_nM = 300  # for aquaculture zones
   ```
2. Increase `benthic_CH4_flux` to 300-500 µmol/m²/day
3. Add rice paddy drainage as point sources during wet season

---

### 5. N₂O (Nitrous Oxide)

**Observed Pattern**:
- 8-15 nmol/L at mouth
- 10-25 nmol/L mid-estuary
- Up to 175 nmol/L at some upstream locations

**Model Behavior**:
- 8 nmol/L at mouth (correct after v1.2.0)
- 5-15 nmol/L mid-estuary (reasonable)
- <5 nmol/L upstream (too low)

**Root Cause**:
1. Missing NH₄ inputs → less nitrification → less N₂O
2. Missing agricultural N₂O sources

**Required Fixes**:
1. Add NH₄ to lateral sources from fertilizer runoff
2. Consider direct N₂O inputs from agricultural drainage

---

### 6. SPM (Suspended Particulate Matter)

**Observed Pattern**: Highly dynamic, 8-35 mg/L with spatial variability

**Model Behavior**: Relatively stable 10-15 mg/L

**Root Cause**:
1. Erosion/deposition too weak
2. Missing tidal pumping effects
3. Constant settling velocity (no flocculation dynamics)

**Required Fixes**:
1. Tune erosion rate (increase `M_ero`)
2. Add spatial variation in bed composition
3. Implement dynamic flocculation with salinity

---

## Priority Development Items

### High Priority (Critical for Science)

1. **Lateral TOC/CH₄/N₂O inputs**: Add these species to the lateral sources system
2. **Tune respiration rate**: Reduce pCO₂ overestimation when reactions ON
3. **Rice paddy drainage**: Major CH₄ source not captured

### Medium Priority (Improve Validation)

4. **Dynamic SPM**: Better erosion/deposition with tidal variation
5. **NH₄ lateral inputs**: Drive N₂O production upstream
6. **Alkalinity weathering**: Improve pH in freshwater zones

### Lower Priority (Refinements)

7. **Ebullition**: Direct CH₄ bubble flux in shallow areas
8. **Sediment-water N cycling**: Coupled nitrification-denitrification
9. **Phytoplankton seasonal dynamics**: Currently too static

---

## Recommended Validation Workflow

1. **Run with ReactionMode = OFF first** to validate transport and boundary conditions
2. **Check all species at ocean boundary** (should match forcing file values)
3. **Compare conservative species** (salinity, alkalinity if no reactions)
4. **Enable reactions** and compare O₂, pCO₂, nutrients
5. **Identify gaps** → likely missing lateral inputs
6. **Add lateral sources** for species with upstream bias
7. **Re-validate** with updated lateral sources

---

## References

- Field data: Mekong Delta Campaign, March 2025
- Validation script: `scripts/compare_with_validation.py`
- Validation data: `INPUT/Cases/Mekong_Delta_Full/Validation/Mekong_Mar2025_model_units.csv`
- Output figures: `OUTPUT/Mekong_Delta_Full/figures/model_validation_comparison.png`
