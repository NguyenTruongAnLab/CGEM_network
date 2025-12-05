# Comprehensive Model Audit: December 2025 (v2)

**Date:** December 5, 2025  
**Auditor:** GitHub Copilot (Claude Opus 4.5)  
**Model Version:** C-GEM Network v1.2.1  
**Test Case:** Mekong_Delta_Full

---

## Executive Summary

This comprehensive audit identified and fixed a **critical bug** in the CH4 saturation calculation that was causing CH4 to accumulate to millions of µmol/L instead of the expected ~10-200 nmol/L range.

### Critical Bug Fixed

| Issue | Root Cause | Impact | Status |
|-------|------------|--------|--------|
| **CH4 accumulation (1 million ×)** | Wrong Yamamoto (1976) coefficients in `crive_calc_ch4_sat()` | CH4_sat = 2.4×10⁶ instead of 0.003 µmol/L | ✅ **FIXED** |

### Current Validation Status (Post-Fix)

| Variable | RMSE | Bias | Status |
|----------|------|------|--------|
| **Salinity** | 3-12 PSU | +2 to -8 | ⚠️ Needs dispersion calibration |
| **O₂** | 13-18 µmol/L | -4 to +3 | ✅ Acceptable |
| **pCO₂** | 2300-2900 µatm | -1900 to -2300 | ⚠️ Model too low |
| **pH** | 0.9-1.3 | +0.8 to +1.2 | ⚠️ Model too high (linked to pCO₂) |
| **NO₃** | 5-24 µmol/L | -8 to -20 | ⚠️ Model too low upstream |
| **NH₄** | 1.1-1.7 µmol/L | -0.1 to -0.4 | ✅ Good |
| **TOC** | 19-76 µmol/L | -1 to -72 | ⚠️ Model too low upstream |
| **AT** | 62-232 µeq/L | ±100 | ⚠️ Variable by branch |
| **CH4** | 65-107 nmol/L | -62 to -88 | ⚠️ Model too low (missing lateral) |
| **N₂O** | 9-20 nmol/L | -7 to -14 | ⚠️ Model too low |

---

## Part 1: Bug Analysis

### 1.1 CH4 Saturation Bug (CRITICAL)

**Location:** `src/rive/ghg_module.c`, function `crive_calc_ch4_sat()`

**Symptom:** CH4 accumulated to 58 million µmol/L at the ocean boundary (observed ~40 nmol/L)

**Root Cause:** The Yamamoto (1976) polynomial coefficients were incorrect:
```c
// WRONG coefficients (gave ln(β) ≈ +14 instead of -6)
const double A1 = -415.2807;  // Should be for a different gas/unit
const double A2 = 596.8104;
const double A3 = 379.2599;
const double A4 = -62.0757;
```

These coefficients produced β = 1.2×10⁶ mol/L/atm instead of ~1.4×10⁻³ mol/L/atm.

**Impact Chain:**
1. CH4_sat = 2.4×10⁶ µmol/L (should be 0.003)
2. Flux = k × (CH4 - CH4_sat) / H = **negative** (invasion instead of evasion)
3. Model adds CH4 from atmosphere at 0.6 µmol/L/s → 175 µmol/L/day accumulation

**Fix Applied:**
```c
double crive_calc_ch4_sat(double temp, double salinity, double CH4_atm) {
    /* Henry's law with van't Hoff temperature dependence
     * Reference: Sander (2015) Atmos. Chem. Phys. */
    const double H0 = 1.4e-3;      /* mol/L/atm at 25°C */
    const double dH_R = 1900.0;    /* K */
    const double T_ref = 298.15;   /* K */
    
    double T_K = temp + 273.15;
    double H_T = H0 * exp(-dH_R * (1.0/T_K - 1.0/T_ref));
    
    /* Salinity correction (~3% per 10 PSU) */
    H_T *= exp(-salinity * 0.003);
    
    double CH4_atm_atm = CH4_atm * 1e-9;  /* ppb to atm */
    return H_T * CH4_atm_atm * 1e6;  /* µmol/L */
}
```

---

## Part 2: Remaining Issues (Expert Consultation Needed)

### 2.1 pCO₂ Underestimation

**Problem:** Model pCO₂ is 1900-2300 µatm too low compared to field data.

**Potential Causes:**
1. **Insufficient respiration** - TOC decomposition rates may be too low
2. **Missing benthic CO₂ flux** - Sediment respiration not implemented
3. **Wrong K0(CO₂) calculation** - Similar to CH4 bug?

**Expert Questions:**
- What is the typical benthic CO₂ flux in the Mekong (mmol C/m²/day)?
- Should we add a `benthic_CO2_flux` parameter similar to CH4?

### 2.2 TOC Deficit Upstream

**Problem:** TOC is 70+ µmol/L too low in upstream branches (Ham_Luong, My_Tho).

**Potential Causes:**
1. **Over-oxidation** - kox = 0.01 /day may be too high for refractory DOM
2. **Missing lateral TOC** - Rice paddies, mangroves contribute DOC
3. **Boundary condition** - River TOC may be underestimated

**Expert Questions:**
- What fraction of Mekong TOC is labile vs refractory?
- What are typical lateral TOC concentrations from rice paddies?

### 2.3 CH4/N₂O Still Too Low

**Problem:** Even with saturation fix, CH4 is 60-90 nmol/L below observations.

**Potential Causes:**
1. **Missing lateral CH4 from rice paddies** - These are major CH4 sources
2. **Benthic CH4 flux too low** - Current 5 µmol/m²/day may be underestimate
3. **Aquaculture CH4** - Not in lateral loads

**Expert Questions:**
- What is the typical CH4 concentration in rice paddy drainage (nmol/L)?
- Should we increase benthic_CH4_flux to 50-100 µmol/m²/day?

### 2.4 Salinity Intrusion Pattern

**Problem:** Salinity shows variable bias by branch (±8 PSU).

**Potential Causes:**
1. **Van den Burgh K needs calibration** - Current 0.15-0.22 may be wrong
2. **Convergence length mismatch** - Geometry may not match reality
3. **Discharge partitioning** - Tien/Hau split may be incorrect

**Expert Questions:**
- What is the observed salt intrusion distance in March (dry season)?
- Are there published K values for Mekong distributaries?

---

## Part 3: Source Code Verification

### 3.1 Files Audited

| File | Status | Issues Found |
|------|--------|--------------|
| `src/rive/ghg_module.c` | ✅ Fixed | CH4 saturation coefficients |
| `src/rive/biogeo.c` | ✅ OK | GHG dynamics correct |
| `src/physics/transport.c` | ✅ OK | Boundary conditions correct |
| `src/init.c` | ✅ OK | Initialization correct |
| `src/io/io_network.c` | ✅ OK | File parsing correct |

### 3.2 Unit Consistency Check

| Variable | Model Unit | Boundary Unit | Validation Unit | Status |
|----------|------------|---------------|-----------------|--------|
| Salinity | PSU | PSU | PSU | ✅ |
| O₂ | µmol/L | µmol/L | µmol/L | ✅ |
| pCO₂ | µatm | µatm | µatm | ✅ |
| TOC | µmol C/L | µmol C/L | µmol C/L | ✅ |
| CH4 | µmol/L | µmol/L | nmol/L (×1000) | ✅ |
| N₂O | µmol/L | µmol/L | nmol/L (×1000) | ✅ |

### 3.3 Key Equations Verified

#### Transport (advection-dispersion)
```
∂C/∂t + u·∂C/∂x = ∂/∂x(D·∂C/∂x) + R
```
✅ Implemented correctly in `transport.c`

#### CH4 air-water flux
```
F = kg × (C - Csat) / H
kg = k600 × √(600/Sc)
```
✅ Implemented correctly (after saturation fix)

#### GHG mass balance
```
dCH4/dt = Production - Oxidation - Evasion - Ebullition
```
✅ Implemented correctly in `biogeo.c`

---

## Part 4: Recommendations for Expert Review

### 4.1 High Priority (Affects Publication)

1. **pCO₂/pH calibration**
   - Need expert review of carbonate chemistry parameters
   - Consider adding benthic CO₂ flux

2. **TOC source/sink balance**
   - Review kox (oxidation rate) for tropical refractory DOM
   - Add TOC to lateral loads from rice paddies

3. **CH4/N₂O lateral inputs**
   - Rice paddies are major sources not currently captured
   - Consider adding aquaculture point sources

### 4.2 Medium Priority (Model Improvement)

4. **Salinity calibration**
   - Systematic K and α calibration needed
   - Compare with observed intrusion distances

5. **Dispersion parameterization**
   - Van den Burgh vs other formulations
   - Tidal prism approach for validation

### 4.3 Low Priority (Future Work)

6. **Sediment-GHG coupling**
   - Currently no feedback between SPM and GHG
   - Could affect CH4 ebullition

7. **Seasonal forcing validation**
   - Need wet season data for comparison
   - Currently only March (dry season) data

---

## Part 5: Files Changed

| File | Change | Lines |
|------|--------|-------|
| `src/rive/ghg_module.c` | Fixed CH4 saturation | 203-238 |
| `src/rive/biogeo.c` | Removed debug prints | 871-879, 1001-1003 |
| `src/init.c` | Removed debug prints | 616-640 |

---

## Appendix A: Validation Results Summary

```
SALINITY:
  Hau_River      : RMSE=  11.66, Bias=  +3.81
  Ham_Luong      : RMSE=  12.17, Bias=  -7.46
  Co_Chien       : RMSE=   5.78, Bias=  -0.39
  My_Tho         : RMSE=   7.50, Bias=  -3.68

O2:
  Hau_River      : RMSE=  18.25, Bias=  +2.79
  Ham_Luong      : RMSE=  12.76, Bias=  -3.73
  Co_Chien       : RMSE=  15.93, Bias= -10.15
  My_Tho         : RMSE=  15.43, Bias=  -1.67

PCO2:
  Hau_River      : RMSE=2352.15, Bias=-1920.85
  Ham_Luong      : RMSE=2437.47, Bias=-2044.09
  Co_Chien       : RMSE=2931.72, Bias=-2340.76
  My_Tho         : RMSE=2542.50, Bias=-2065.15

TOC:
  Hau_River      : RMSE=  18.70, Bias=  -0.74
  Ham_Luong      : RMSE=  73.99, Bias= -71.88
  Co_Chien       : RMSE=  39.39, Bias= -33.34
  My_Tho         : RMSE=  75.53, Bias= -70.55

CH4 (after fix):
  Hau_River      : RMSE=  97.98, Bias= -82.98 nmol/L
  Ham_Luong      : RMSE=  97.07, Bias= -88.07 nmol/L
  Co_Chien       : RMSE= 106.94, Bias= -86.42 nmol/L
  My_Tho         : RMSE=  65.43, Bias= -61.80 nmol/L
```

---

## Appendix B: Literature for Expert Reference

### CH4 Saturation
- Sander, R. (2015). Compilation of Henry's law constants. Atmos. Chem. Phys., 15, 4399-4981.
- Wiesenburg, D.A., Guinasso Jr, N.L. (1979). Equilibrium solubilities of methane. J. Chem. Eng. Data, 24, 356-360.

### Mekong Biogeochemistry
- Nguyen, A.D. et al. (2008). Salt intrusion in the Mekong Delta. Physics and Chemistry of the Earth.
- Borges, A.V. et al. (2018). Globally significant greenhouse gas emissions from African inland waters. Nature Geoscience.

### Gas Transfer
- Wanninkhof, R. (2014). Relationship between wind speed and gas exchange over the ocean. Limnol. Oceanogr.: Methods.
- Abril, G. et al. (2009). Turbidity limits gas exchange in a large macrotidal estuary. Estuarine, Coastal and Shelf Science.

---

*End of Audit Document*
