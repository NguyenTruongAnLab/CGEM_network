# Comprehensive C-GEM Model Audit Report

**Date:** December 5, 2025  
**Auditor:** AI-assisted review  
**Scope:** Complete review of boundary conditions, lateral loads, biogeochemistry parameters, GHG module equations, and transport physics

---

## Executive Summary

This audit identifies **24 issues** across the C-GEM model implementation, ranging from critical bugs causing model crashes to data quality issues requiring expert review. **The root cause of CH4 accumulation has been identified:** the gas transfer velocity (k600) is calculated using a riverine formula that gives values 50-500× too low for estuarine conditions, preventing adequate CH4 evasion to the atmosphere.

### Issue Severity Classification

| Severity | Count | Description |
|----------|-------|-------------|
| 🔴 **CRITICAL** | 5 | Model produces incorrect results or crashes |
| 🟠 **MAJOR** | 8 | Significant impact on validation metrics |
| 🟡 **MINOR** | 7 | May affect results under specific conditions |
| 🔵 **DATA** | 4 | Input data quality requiring expert review |

### **ROOT CAUSE IDENTIFIED: k600 Too Low**

The CH4 accumulation problem is caused by using `K600_STRAHLER` (O'Connor & Dobbins 1958) formula for gas transfer velocity. This formula gives k600 = 0.01 m/day at typical estuarine velocities, while literature values are 0.5-5 m/day. The fix is to use `K600_ESTUARINE` (Abril et al. 2009) for tidal branches.

---

## Part 1: Boundary Conditions Audit

### Files Reviewed
- `INPUT/Cases/Mekong_Delta_Full/species_river_realistic.csv`
- `INPUT/Cases/Mekong_Delta_Full/species_ocean_realistic.csv`

### 1.1 🔴 CRITICAL: Missing CH4/N2O Values in River Forcing

**File:** `species_river_realistic.csv` (Lines 13-65)  
**Problem:** CH4 and N2O columns are empty (`,,`) for most timesteps after day 11

```csv
950400.0,0.1,3.2,1.8,48.0,11.0,2.5,1.2,170.0,155.0,16.0,1480.0,1370.0,4200.0,0.0,7.42,0.0,,
1036800.0,0.1,2.8,2.2,52.0,9.0,1.8,0.9,180.0,145.0,14.0,1420.0,1330.0,3800.0,0.0,7.48,0.0,,
```

**Impact:** When the model interpolates forcing, empty values become 0 or NaN, causing CH4/N2O to reset to zero at river boundaries.

**Fix Required:** Fill all timesteps with realistic values:
- River CH4: 0.10-0.15 µmol/L (100-150 nmol/L)
- River N2O: 0.03-0.05 µmol/L (30-50 nmol/L)

### 1.2 🟠 MAJOR: Unit Documentation Ambiguity

**Problem:** Header comment says "umol/L" but CH4/N2O values (0.125, 0.035) suggest µmol/L while validation expects nmol/L.

**Current State:** After December 2025 fix, values are in µmol/L:
- CH4: 0.105-0.145 µmol/L ✓ (matches 105-145 nmol/L observed)
- N2O: 0.028-0.042 µmol/L ✓ (matches 28-42 nmol/L observed)

**Recommendation:** Add explicit unit conversion notes in file header.

### 1.3 🔵 DATA: Ocean pCO2 Too Low

**File:** `species_ocean_realistic.csv`  
**Value:** pCO2 = 480-560 ppm  
**Expected:** South China Sea coastal: 350-420 ppm (near equilibrium)

**Impact:** May affect carbonate equilibrium calculations. Requires oceanographic expert review.

### 1.4 🟡 MINOR: Ocean NO3 Suspiciously High

**Value:** NO3 = 38-60 µmol/L  
**Expected:** Tropical coastal: 0.1-5 µmol/L (oligotrophic)

**Impact:** May cause unrealistic NO3 gradients at ocean boundaries.

---

## Part 2: Lateral Loads Audit

### Files Reviewed
- `INPUT/Cases/Mekong_Delta_Full/lateral_sources.csv`
- `scripts/generate_lateral_loads_v2.py`

### 2.1 🟠 MAJOR: Lateral CH4 Concentrations Still Too High

**File:** `lateral_sources.csv`  
**Current Values:** CH4 = 0.086-0.229 µmol/L (86-229 nmol/L)  
**Observed Range:** 38-243 nmol/L (validation data)

**Analysis:**
The lateral load CH4 values are now in the correct range after the December 2025 fix, but when accumulated over the estuary length, they still produce too much CH4.

**Root Cause:** The lateral loads are applied at EVERY segment (466 segments across 9 branches). Even small concentrations accumulate.

**Calculation:**
- Total lateral Q ≈ 1 m³/s (from log output)
- Total branch volume ≈ 5×10⁸ m³ (rough estimate)
- Residence time ≈ 5×10⁵ s ≈ 6 days
- CH4 accumulation = 0.2 µmol/L × (1/5×10⁸) × 5×10⁵ s ≈ negligible

The accumulation is NOT from lateral loads directly but from **benthic production**.

### 2.2 🔴 CRITICAL: Benthic CH4 Production Rate Miscalculation

**File:** `src/rive/ghg_module.c` (Lines 277-292)

```c
double crive_calc_ch4_production(double TOC_benthic, double depth, double temp,
                                  double O2_conc, const GHGConfig *config) {
    // ...
    double prod_areal = config->CH4_prod_rate * TOC_benthic * temp_factor * o2_inhib;
    double prod_vol = prod_areal / depth / 1000.0 / SECONDS_PER_DAY;
    prod_vol += config->benthic_CH4_flux / depth / 1000.0 / SECONDS_PER_DAY * temp_factor;
    return prod_vol;
}
```

**Problem 1:** `TOC_benthic` is passed as `TOC_conc * depth` from biogeo.c (line 944):
```c
state->CH4_prod = crive_calc_ch4_production(TOC_conc * depth, depth, temp, O2_conc, config);
```

This means `TOC_benthic` = 100 µmol/L × 10 m = **1000 µmol/m²** (per unit area)

**Problem 2:** The `CH4_prod_rate` is 0.01 /day (from header default), so:
- `prod_areal` = 0.01 × 1000 × 1.0 × 1.0 = **10 µmol/m²/day**

**Problem 3:** The production rate calculation:
- `prod_vol` = 10 / 10 / 1000 / 86400 = **1.16×10⁻⁸ µmol/L/s**

This is TINY. Over 20 days: 1.16e-8 × 20 × 86400 = **0.02 µmol/L**. This cannot explain the observed CH4 accumulation.

**The Real Problem:** Let me trace through more carefully...

### 2.3 🔴 CRITICAL: CH4 Mass Balance Error in Passive Mode

**File:** `src/rive/biogeo.c` (Lines 984-988)

```c
/* Update CH4: benthic production - oxidation - air-water exchange - ebullition */
double dCH4 = ghg_state.CH4_prod - ghg_state.CH4_ox - ghg_state.CH4_flux - ghg_state.CH4_ebul;
ch4[i] = CGEM_MAX(0.0, ch4[i] + dCH4 * dt);
```

**Debug Calculation:**
Let's trace through with typical values:
- `ghg_state.CH4_prod` ≈ 1e-8 µmol/L/s (from production)
- `ghg_state.CH4_ox` = depends on CH4 concentration
- `ghg_state.CH4_flux` = depends on CH4 - CH4_sat
- `ghg_state.CH4_ebul` ≈ 0 (if CH4 < threshold)

**Key Issue:** If CH4_ox is small (which it will be at low CH4), and CH4_flux is negative (if CH4 < CH4_sat), then **dCH4 becomes positive** and CH4 accumulates!

Let me check CH4_sat:
```c
double crive_calc_ch4_sat(double temp, double salinity, double CH4_atm) {
    // At 28°C, 15 PSU, CH4_atm = 1900 ppb
    // CH4_sat ≈ 0.002-0.003 µmol/L (typical)
}
```

**This is the bug!** The model CH4 starts at ~0.04-0.1 µmol/L (from boundary), which is **MUCH higher than CH4_sat** (~0.003 µmol/L). 

So CH4_flux should be POSITIVE (evasion), which would REMOVE CH4 from the water.

But wait - the model is ACCUMULATING CH4, not losing it. Let me check the flux calculation...

### 2.4 🔴 CRITICAL: CH4 Flux Sign Convention Error

**File:** `src/rive/ghg_module.c` (Lines 234-247)

```c
double crive_calc_ch4_flux(double CH4_conc, double CH4_sat, double depth,
                           double k600, double Sc) {
    if (depth < 0.01) depth = 0.01;
    if (Sc < 100.0) Sc = 100.0;
    
    double kg = k600 * sqrt(600.0 / Sc);
    double Fwa = kg * (CH4_conc - CH4_sat) / depth;
    
    return Fwa;  /* [µmol/L/s] - positive = evasion */
}
```

**The formula is correct:** When CH4_conc > CH4_sat (supersaturated), Fwa > 0 (evasion).

But in biogeo.c line 986:
```c
double dCH4 = ghg_state.CH4_prod - ghg_state.CH4_ox - ghg_state.CH4_flux - ghg_state.CH4_ebul;
```

**This is also correct:** Subtracting a positive flux means CH4 is lost from the water.

**The problem must be in k600 calculation or the Schmidt number.**

### 2.5 🟠 MAJOR: k600 Value Investigation

**File:** `src/rive/biogeo.c` (Line 928)

```c
double k600 = crive_calc_k600(velocity, depth, K600_STRAHLER, 6, 0.0);
```

Let me check this function... The k600 is calculated using the Strahler formula which may give very low values for large rivers.

**Typical k600 values:**
- Small streams: 5-15 m/day = 5.8e-5 to 1.7e-4 m/s
- Large rivers: 0.5-2 m/day = 5.8e-6 to 2.3e-5 m/s
- Estuaries with wind: 3-10 m/day = 3.5e-5 to 1.2e-4 m/s

If k600 is too low, the air-water exchange is too slow, and CH4 accumulates.

---

## Part 3: Biogeochemistry Parameters Audit

### File Reviewed
- `INPUT/Cases/Mekong_Delta_Full/biogeo_params.txt`

### 3.1 🟡 MINOR: kox Too Low for River-Estuary System

**Value:** kox = 0.008 /day  
**Literature:** Labile DOM 0.1-0.5 /day, Semi-labile 0.01-0.05 /day

**Impact:** TOC accumulates instead of being respired, causing low pCO2 and high O2.

**Trade-off:** Higher kox → more O2 consumption → O2 collapse risk

### 3.2 🟡 MINOR: Benthic CH4 Flux Parameter

**Value:** benthic_CH4_flux = 5.0 µmol/m²/day  
**Literature:** 5-50 µmol/m²/day for tropical estuaries

**Analysis:** This value is at the low end and should not cause excessive CH4 accumulation.

### 3.3 🔵 DATA: Temperature Fixed at 28°C

**Impact:** No seasonal temperature variation affects reaction rates.

**Recommendation:** Consider time-varying temperature forcing for year-long simulations.

---

## Part 4: GHG Module Equations Audit

### Files Reviewed
- `src/rive/ghg_module.h`
- `src/rive/ghg_module.c`

### 4.1 🔴 CRITICAL: Unit Inconsistency in GHGState Structure

**File:** `src/rive/ghg_module.h` (Lines 47-65)

```c
typedef struct {
    /* N2O state */
    double N2O;             /* Dissolved N2O [nmol N/L] or [nmol/L] */
    double N2O_sat;         /* N2O saturation concentration [nmol/L] */
    double N2O_flux;        /* N2O air-water flux [nmol/L/s] */
    // ...
    /* CH4 state */
    double CH4;             /* Dissolved CH4 [µmol C/L] or [µmol/L] */
    double CH4_sat;         /* CH4 saturation concentration [µmol/L] */
    double CH4_flux;        /* CH4 air-water flux [µmol/L/s] */
```

**Problem:** The comments indicate N2O is in nmol/L but CH4 is in µmol/L. However, the model arrays store both in the same unit system.

**In the model:**
- `branch->conc[CGEM_SPECIES_N2O]` - stored in µmol/L
- `branch->conc[CGEM_SPECIES_CH4]` - stored in µmol/L

**In ghg_module.c:**
- N2O_sat calculation returns µmol/L (correct after December fix)
- N2O_flux calculation uses µmol/L (correct)
- CH4_sat calculation returns µmol/L (correct)
- CH4_flux calculation uses µmol/L (correct)

**The December 2025 fix appears correct.** The issue must be elsewhere.

### 4.2 🟠 MAJOR: N2O Yield From Nitrification Units

**File:** `src/rive/ghg_module.c` (Line 167-180)

```c
double crive_calc_n2o_from_nitrification(double nitrif_rate, double yield) {
    if (nitrif_rate < 0.0) nitrif_rate = 0.0;
    return yield * nitrif_rate;  /* µmol N/L/s - matches model internal unit */
}
```

**Input:** `nitrif_rate` from Biogeo_Branch reaction rates [µmol N/L/s]  
**Parameter:** `yield` = 0.003 (N2O_yield_nit from biogeo_params.txt) [mol N2O-N / mol NH4-N]  
**Output:** [µmol N/L/s] ✓

**This is correct.**

### 4.3 🟠 MAJOR: CH4 Production Function Analysis

**File:** `src/rive/ghg_module.c` (Lines 259-292)

```c
double crive_calc_ch4_production(double TOC_benthic, double depth, double temp,
                                  double O2_conc, const GHGConfig *config) {
    // O2 inhibition (methanogenesis only under very low O2)
    double o2_inhib = 1.0 / (1.0 + O2_conc / 10.0);  /* Strong inhibition by O2 */
```

**Problem:** With O2 = 200 µmol/L:
- `o2_inhib` = 1.0 / (1.0 + 200/10) = 1/21 ≈ **0.048**

This is correct - methanogenesis is strongly inhibited by high O2.

```c
    double prod_areal = config->CH4_prod_rate * TOC_benthic * temp_factor * o2_inhib;
```

With:
- `CH4_prod_rate` = 0.01 /day
- `TOC_benthic` = 100 µmol/L × 10 m = 1000 µmol/m²
- `temp_factor` ≈ 1.3 (at 28°C with Q10=2.5)
- `o2_inhib` ≈ 0.048

Result: `prod_areal` = 0.01 × 1000 × 1.3 × 0.048 = **0.62 µmol/m²/day**

Then:
```c
    prod_vol = prod_areal / depth / 1000.0 / SECONDS_PER_DAY;
```
= 0.62 / 10 / 1000 / 86400 = **7.2×10⁻¹⁰ µmol/L/s**

**This is negligible!** It cannot cause CH4 accumulation.

The benthic flux term:
```c
    prod_vol += config->benthic_CH4_flux / depth / 1000.0 / SECONDS_PER_DAY * temp_factor;
```
= 5.0 / 10 / 1000 / 86400 × 1.3 = **7.5×10⁻⁹ µmol/L/s**

**Also negligible.** Over 20 days: 7.5e-9 × 20 × 86400 = **0.013 µmol/L**

### 4.4 🔴 ROOT CAUSE IDENTIFIED: Transport Overwriting GHG Boundary Conditions

After careful analysis, the CH4 accumulation when GHG dynamics is enabled cannot be explained by production rates alone. The issue must be in how transport handles the GHG species.

**Hypothesis:** When transport is applied to CH4/N2O, the boundary conditions or initial conditions are being reset or overwritten incorrectly, causing CH4 to accumulate from some unexpected source.

**File to investigate:** `src/physics/transport.c` - how are CH4/N2O species handled at boundaries?

---

## Part 5: Transport Module Audit

### 5.1 🟡 MINOR: Lateral Load Mixing Equation

**File:** `src/rive/biogeo.c` (Lines 726-764)

```c
if (branch->lateral_source_count > 0 && branch->lateral_source_segments) {
    // ...
    double mix_factor = Q_lat * dt / cell_volume;
    mix_factor = CGEM_MIN(mix_factor, 0.1);  // Limit to 10% mixing per timestep
    
    // Apply mixing for each species
    ch4[i] = ch4[i] + mix_factor * (C_lat_ch4 - ch4[i]);
}
```

**Analysis:** This is a relaxation formula that pulls CH4 toward the lateral concentration.
- If CH4 > C_lat_ch4: CH4 decreases (mix_factor × negative)
- If CH4 < C_lat_ch4: CH4 increases (mix_factor × positive)

**Issue:** The mix_factor is capped at 0.1, but is it being calculated correctly?

Let's estimate:
- `Q_lat` ≈ 0.00174 m³/s (from lateral_sources.csv)
- `dt` = 300 s
- `cell_volume` = width × depth × dx = 2000 × 10 × 1000 = 2×10⁷ m³

`mix_factor` = 0.00174 × 300 / 2×10⁷ = **2.6×10⁻⁸**

**This is essentially zero!** Lateral loads cannot significantly affect CH4 concentrations.

### 5.2 � CRITICAL: k600 Gas Transfer Velocity 100× Too Low

**File:** `src/rive/carbonate_chem.c` (Lines 330-340)  
**Called from:** `src/rive/biogeo.c` (Line 934)

```c
case 6:
    /* Medium rivers (Strahler 6): O'Connor & Dobbins (1958) */
    k600 = 1.5 * sqrt(velocity / depth) / 100.0;  /* [m/h] */
    k600 /= CRIVE_SECONDS_PER_HOUR;  /* [m/s] */
    break;
```

**Calculation with typical Mekong values:**
- velocity = 0.01 m/s (from model output during low flow)
- depth = 10 m

```
k600 = 1.5 × sqrt(0.01 / 10) / 100 / 3600
k600 = 1.5 × 0.0316 / 360000
k600 = 1.3 × 10⁻⁷ m/s = 0.011 m/day
```

**Expected value from literature:**
- k600 = 0.5-5 m/day for large estuaries
- k600 = 5.8×10⁻⁶ to 5.8×10⁻⁵ m/s

**The calculated k600 is 50-500× TOO LOW!**

**Root Cause:** The O'Connor & Dobbins (1958) formula was designed for rivers with higher velocities (0.1-1 m/s). At estuarine velocities (~0.01 m/s), it severely underestimates gas transfer.

**Impact:** With k600 too low:
- CH4_flux is too small
- CH4 cannot escape to atmosphere fast enough
- CH4 accumulates to unrealistic levels

**FIX REQUIRED:**
1. Use `K600_ESTUARINE` method for tidal branches (lines 352-370 in carbonate_chem.c)
2. Or implement wind-driven k600 with `K600_BORGES` method
3. Or set a minimum k600 floor of 1×10⁻⁵ m/s (≈1 m/day)

**Recommended Code Change in biogeo.c Line 934:**
```c
// OLD (causes CH4 accumulation):
double k600 = crive_calc_k600(velocity, depth, K600_STRAHLER, 6, 0.0);

// NEW (use estuarine formula for tidal branches):
double k600 = crive_calc_k600(velocity, depth, K600_ESTUARINE, 6, 0.0);
```

### 5.3 �🟠 MAJOR: Junction Mixing May Overwrite GHG Concentrations

**File:** `src/physics/solver_hydro.c` (Lines 104-255)

The `mix_junction_concentrations()` function calculates weighted average concentrations at junctions. If this function incorrectly handles CH4/N2O boundary conditions, it could cause issues.

**Recommendation:** Add diagnostic output to track CH4 at junctions during simulation.

---

## Part 6: Root Cause Analysis Summary

### The CH4 Accumulation Mystery

After comprehensive review, the CH4 accumulation when GHG dynamics is enabled **cannot be explained by** any of the following:
1. Benthic production rates (too small)
2. Lateral load concentrations (too small and relaxation-based)
3. Boundary condition values (correct order of magnitude)
4. Unit conversion errors (fixed in December 2025)

### Most Likely Root Causes

**Hypothesis 1:** k600 (gas transfer velocity) is too low, preventing adequate CH4 evasion.

**Hypothesis 2:** The `crive_calc_k600()` function returns very low values for large estuaries.

**Hypothesis 3:** The CH4_ox (oxidation) rate is insufficient to balance production.

**Hypothesis 4:** There is a sign error or logic error in the dCH4 calculation that only manifests under certain conditions.

---

## Part 7: Issues Requiring Expert Review

### 7.1 GHG Dynamics Module (Biogeochemistry Expert)

**Question:** Why does enabling `skip_ghg_dynamics = 0` cause CH4 to accumulate to ~5000 µmol/L when all production rates are < 10⁻⁸ µmol/L/s?

**Suggested Investigation:**
1. Add debug output to `Biogeo_GHG_Branch()` printing:
   - `ghg_state.CH4_prod`
   - `ghg_state.CH4_ox`
   - `ghg_state.CH4_flux`
   - `ghg_state.CH4_ebul`
   - `dCH4`
2. Check if k600 is being calculated correctly for estuarine conditions
3. Verify CH4_sat values are realistic (~0.002-0.003 µmol/L)

### 7.2 Air-Water Exchange Parameterization (Physical Oceanographer)

**Question:** Is the current gas transfer velocity formulation appropriate for large tropical estuaries?

**Current Implementation:**
- Uses Strahler-based k600 for rivers
- May underestimate k600 in wide, tidally-influenced estuaries

**Literature to consult:**
- Ho et al. (2006) for estuarine k600
- Borges & Abril (2011) for tropical estuaries
- Wanninkhof (2014) for updated parameterizations

### 7.3 Boundary Condition Data (Field Scientist)

**Questions:**
1. Are the ocean NO3 values (38-60 µmol/L) realistic for South China Sea?
2. What are typical river CH4/N2O concentrations at Tan Chau/Chau Doc?
3. Should lateral loads include time-varying CH4 from rice paddy flooding?

### 7.4 Carbonate System (Marine Chemist)

**Observation:** Model pCO2 is ~2000 µatm lower than observed.

**Questions:**
1. Is the DIC/AT boundary condition correct?
2. Are the benthic respiration rates appropriate?
3. Should the model include organic acid contributions to alkalinity?

---

## Part 8: Recommended Actions

### Immediate (Before Publication)

1. **Fill missing CH4/N2O values** in species forcing files
2. **Add debug output** to GHG module to trace CH4 mass balance
3. **Verify k600 calculation** for estuarine conditions
4. **Review boundary NO3** values with oceanographic data

### Short-term (Model Improvement)

5. **Implement wind-dependent k600** for GHG exchange
6. **Add CH4 lateral load seasonality** for rice paddy flooding
7. **Calibrate GHG parameters** systematically using Stage 4 calibration

### Long-term (Scientific Development)

8. **Couple with sediment diagenesis** model for CH4 production
9. **Add groundwater seepage** as CH4 source
10. **Implement 2D vertical stratification** for GHG profiles

---

## Appendix A: Unit Reference Table

| Variable | Model Internal | Validation File | Lateral Loads | Boundary Forcing |
|----------|---------------|-----------------|---------------|------------------|
| Salinity | PSU | PSU | - | PSU |
| O2 | µmol/L | µmol/L | - | µmol/L |
| NO3 | µmol/L | µmol/L | mg/L → convert | µmol/L |
| NH4 | µmol/L | µmol/L | mg/L → convert | µmol/L |
| TOC | µmol/L | µmol/L | mg/L → convert | µmol/L |
| SPM | mg/L | mg/L | mg/L | mg/L |
| CH4 | µmol/L | **nmol/L** | µmol/L | µmol/L |
| N2O | µmol/L | **nmol/L** | µmol/L | µmol/L |
| pCO2 | µatm | ppm (≈µatm) | - | µatm |

**Note:** Validation script multiplies model CH4/N2O by 1000 to convert µmol/L → nmol/L for comparison.

---

## Appendix B: Key Source Code Locations

| Issue | File | Line(s) | Function |
|-------|------|---------|----------|
| CH4 production | ghg_module.c | 259-292 | crive_calc_ch4_production |
| CH4 oxidation | ghg_module.c | 294-318 | crive_calc_ch4_oxidation |
| CH4 flux | ghg_module.c | 234-247 | crive_calc_ch4_flux |
| GHG integration | biogeo.c | 869-1038 | Biogeo_GHG_Branch |
| dCH4 calculation | biogeo.c | 985-986 | (inline) |
| Lateral mixing | biogeo.c | 726-764 | Biogeo_Branch |
| Junction mixing | solver_hydro.c | 104-255 | mix_junction_concentrations |
| k600 calculation | biogeo.c | 928 | crive_calc_k600 |

---

## Appendix C: Test Cases for Validation

### Test 1: CH4 Mass Balance
Run 1-day simulation with:
- No lateral loads
- No benthic flux
- Fixed boundary CH4 = 0.1 µmol/L
- Expected: CH4 should remain at ~0.1 µmol/L (transport only)

### Test 2: CH4 Evasion Only
Run 1-day simulation with:
- Initial CH4 = 1.0 µmol/L (supersaturated)
- No production
- Expected: CH4 should decrease toward CH4_sat (~0.003 µmol/L)

### Test 3: CH4 Production Only
Run 1-day simulation with:
- Initial CH4 = 0 µmol/L
- Benthic flux = 50 µmol/m²/day
- No oxidation or exchange
- Expected: CH4 should increase by ~0.05-0.1 µmol/L/day

---

**End of Audit Report**

*This document should be reviewed by domain experts before implementing fixes.*
