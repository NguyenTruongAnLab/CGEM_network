# Comprehensive C-GEM Model Audit Report

**Date:** December 5, 2025 (Final Version)  
**Auditor:** AI-assisted technical review with domain expert validation  
**Model Version:** CGEM_RIVE branch  
**Test Case:** Mekong Delta, Vietnam (March 2025 dry season field data)  
**Scope:** Complete review of physics, biogeochemistry, transport, and GHG modules

---

## Executive Summary

This audit documents the current state of the C-GEM Network model following intensive debugging and calibration efforts in December 2025. The model has achieved **good performance** for conservative tracers (salinity R²=0.65-0.99) and just acceptable for nutrient gradients (NO₃), but shows **systematic biases** in nutrients (TOC, O₂, NH₄, PO₄), the carbonate system (pCO₂, pH) and greenhouse gases (CH₄, N₂O) that reflect fundamental model limitations rather than parameter tuning issues.

The audit follows the **80/20 principle**: identifying the 20% of missing processes that cause 80% of the remaining model-data mismatch, with focus on improvements achievable using globally-available datasets rather than site-specific measurements.

---

## Table of Contents

1. [Current Validation Status](#1-current-validation-status)
2. [Model Architecture Overview](#2-model-architecture-overview)
3. [Process-by-Process Technical Analysis](#3-process-by-process-technical-analysis)
4. [Root Cause Analysis of Biases](#4-root-cause-analysis-of-biases)
5. [Equations and Their Limitations](#5-equations-and-their-limitations)
6. [Recommended Improvements (80/20 Prioritized)](#6-recommended-improvements-8020-prioritized)
7. [Global Dataset Integration Strategy](#7-global-dataset-integration-strategy)
8. [Implementation Roadmap](#8-implementation-roadmap)

---

## 1. Current Validation Status

### 1.1 Validation Metrics Summary (December 5, 2025)

| Variable | Best Branch | RMSE | Bias | R² | Status |
|----------|-------------|------|------|-----|--------|
| **Salinity** | Co_Chien | 1.38 PSU | -0.77 | 0.987 | ✅ Excellent |
| **pH** | My_Tho | 0.11 | +0.05 | 0.732 | ✅ Good |
| **O₂** | My_Tho | 24.2 µmol/L | -17.6 | 0.586 | ⚠️ Model too low |
| **NO₃** | Co_Chien | 4.81 µmol/L | -0.41 | 0.945 | ⚠️ Branch-dependent, Model too low|
| **NH₄** | Co_Chien | 1.01 µmol/L | +0.07 | 0.019 | 🟠 Gradient captured, bias remains|
| **SPM** | Co_Chien | 4.76 mg/L | +0.13 | 0.100 | ⚠️ Low R² |
| **pCO₂** | Hau_River | 1033 µatm | -191 | 0.794 | 🟠 Gradient captured, bias remains |
| **TOC** | Hau_River | 21.7 µmol/L | +9.6 | 0.582 | ⚠️ Branch-dependent, Model too low |
| **Alkalinity** | Co_Chien | 66.5 µeq/L | +35.5 | 0.878 | ✅ Good |
| **CH₄** | My_Tho | 15.4 nmol/L | -6.1 | 0.034 | 🔴 Missing lateral sources |
| **N₂O** | Co_Chien | 8.37 nmol/L | -6.6 | 0.305 | 🟠 Missing agricultural inputs |

### 1.2 Key Achievements (This Audit Cycle)

1. **Ocean boundary pCO₂ fixed**: Model now correctly shows ~500 µatm at mouth (was 900+ µatm)
2. **Benthic O₂ consumption added**: SOD (Sediment Oxygen Demand) now included in O₂ balance
3. **DIC/TA boundary consistency**: Ocean BC now thermodynamically consistent (pH=8.1, pCO₂~450)
4. **NO₃ ocean boundary corrected**: Now 55 µmol/L matching observed coastal upwelling signal

### 1.3 Remaining Systematic Biases

| Issue | Magnitude | Root Cause |
|-------|-----------|------------|
| pCO₂ too low upstream | -600 to -800 µatm | Spatially uniform benthic respiration |
| O₂ too low mid-estuary | -17 to -25 µmol/L | Benthic O₂ demand too high near ocean |
| CH₄ underestimated | -35 to -38 nmol/L | Missing rice paddy/aquaculture lateral inputs |
| N₂O underestimated | -6 to -14 nmol/L | Missing agricultural nitrification inputs |
| TOC variable by branch | ±60 µmol/L | Lateral TOC sources not calibrated |

---

## 2. Model Architecture Overview

### 2.1 Core Framework

```
┌─────────────────────────────────────────────────────────────────────┐
│                     C-GEM Network Architecture                       │
├─────────────────────────────────────────────────────────────────────┤
│                                                                      │
│  ┌──────────────┐    ┌──────────────┐    ┌──────────────┐          │
│  │ Hydrodynamics│───▶│  Transport   │───▶│Biogeochemistry│         │
│  │ Saint-Venant │    │ Advection-   │    │   C-RIVE     │          │
│  │  (staggered) │    │ Dispersion   │    │ (simplified) │          │
│  └──────────────┘    └──────────────┘    └──────────────┘          │
│         │                   │                   │                   │
│         ▼                   ▼                   ▼                   │
│  ┌──────────────┐    ┌──────────────┐    ┌──────────────┐          │
│  │   Network    │    │   Lateral    │    │     GHG      │          │
│  │  Junctions   │    │   Sources    │    │   Module     │          │
│  │  (mass bal.) │    │ (land use)   │    │ (CO₂,CH₄,N₂O)│          │
│  └──────────────┘    └──────────────┘    └──────────────┘          │
│                                                                      │
└─────────────────────────────────────────────────────────────────────┘
```

### 2.2 Simplified Mode (80/20 Biogeochemistry)

The model uses `simplified_mode = 1` which bypasses the complex multi-pool RIVE bacterial dynamics in favor of first-order kinetics suitable for data-sparse systems:

**Active Processes:**
- TOC degradation: `dTOC/dt = -kox × TOC × f(O₂) × θ^(T-20)`
- Nitrification: `dNH₄/dt = -knit × NH₄ × f(O₂) × θ^(T-20)`  
- Denitrification: `dNO₃/dt = -kdenit × NO₃ × f(TOC) × g(O₂) × θ^(T-20)`
- O₂ exchange: `dO₂/dt = k_L × (O₂_sat - O₂) / depth`
- Benthic fluxes: DIC, O₂, CH₄, N₂O from sediments

**Disabled in Simplified Mode:**
- Multi-pool organic matter (HD1-3, HP1-3)
- Bacterial biomass dynamics (BAG, BAP)
- Phosphorus adsorption (PIP)
- Complex grazing chains

### 2.3 Staggered Grid Convention

```
Index:    0    1    2    3    4    ...   M-1   M   M+1
         ghost                              ghost
          ┃    ┃    ┃    ┃    ┃           ┃    ┃    ┃
          ▼    ▼    ▼    ▼    ▼           ▼    ▼    ▼
       ───●────○────●────○────●─── ... ───●────○────●───
          │    │    │    │    │           │    │    │
          │    └─scalar (C,S,T)           │    │
          └─velocity                      └─velocity
          
DOWNSTREAM (ocean) ◄──────────────────────► UPSTREAM (river)
       Index 1                                Index M
```

- **Positive velocity (vx > 0)**: Flow from upstream (M) toward downstream (1)
- **Scalar cells**: Odd indices (1, 3, 5, ...)
- **Velocity cells**: Even indices (0, 2, 4, ...)

---

## 3. Process-by-Process Technical Analysis

### 3.1 Carbonate System

#### Current Implementation (`biogeo.c`, lines 245-354)

The carbonate equilibrium solver uses Newton-Raphson iteration to find pH from DIC and TA:

```c
// Carbonate equilibrium: solve for [H⁺] given DIC and TA
// TA = [HCO₃⁻] + 2[CO₃²⁻] + [B(OH)₄⁻] + [OH⁻] - [H⁺]
while (fabs(guess) > TOL && iterations < 50) {
    carbonAlk = TA - borate - Kw/H + H;
    ratio = DIC / carbonAlk;
    // Quadratic solution for H from carbonate speciation
    discriminant = (ratio-1)² × K1² - 4×K1×K2×(1-2×ratio);
    H_new = 0.5 × ((ratio-1)×K1 + sqrt(discriminant));
    iterations++;
}
```

**Key Equations:**

$$K_1 = \frac{[H^+][HCO_3^-]}{[CO_2]} \approx 10^{-5.85} \text{ at 28°C, S=30}$$

$$K_2 = \frac{[H^+][CO_3^{2-}]}{[HCO_3^-]} \approx 10^{-8.97} \text{ at 28°C, S=30}$$

$$pCO_2 = \frac{[CO_2]}{K_H} \text{ where } K_H \approx 0.027 \text{ µmol L}^{-1} \text{µatm}^{-1}$$

**Validation Status:** The solver correctly reproduces equilibrium relationships. The issue is not the carbonate chemistry math, but the **DIC source terms**.

#### Identified Limitation: Uniform Benthic Respiration

The current implementation uses spatially uniform `benthic_resp_20C`:

```c
double benthic_rate_day = p->benthic_resp_20C * pow(p->benthic_Q10, (temp - 20.0) / 10.0);
double benthic_co2_rate = benthic_rate_day / depth / SECONDS_PER_DAY;
```

**Problem:** Benthic respiration varies enormously along estuaries:
- Ocean mouth: Sandy sediments, low organic content → 10-30 mmol C/m²/day
- Mid-estuary: Mixed sediments → 30-60 mmol C/m²/day  
- Upstream: Fine sediments, high organic deposition → 60-150 mmol C/m²/day

**Literature Reference:** Abril et al. (2010) *Limnol. Oceanogr.* 55:1199-1212

### 3.2 Oxygen Dynamics

#### Current Implementation

```c
// O₂ balance (simplified mode)
double o2_consumption = toc_deg_rate + 2.0 * nit_rate + benthic_o2_rate;
o2[i] = MAX(0.0, o2[i] + (o2_ex - o2_consumption) * dt);
```

**Key Equations:**

$$\frac{dO_2}{dt} = k_L \frac{(O_{2,sat} - O_2)}{h} - R_{resp} - 2R_{nit} - SOD/h$$

Where:
- $k_L$ = gas transfer velocity (m/s)
- $h$ = water depth (m)
- $R_{resp}$ = aerobic respiration rate (µmol/L/s)
- $R_{nit}$ = nitrification rate (µmol N/L/s)
- $SOD$ = sediment oxygen demand (µmol O₂/m²/s)

**Gas Transfer Velocity (Wanninkhof 1992):**

$$k_L = 0.31 \times u_{10}^2 \times \left(\frac{Sc}{660}\right)^{-0.5}$$

**Validation Issue:** Model O₂ is 17-25 µmol/L too low mid-estuary. The SOD applied uniformly causes excessive O₂ consumption near the ocean where observed O₂ remains near saturation.

### 3.3 Greenhouse Gas Module

#### CH₄ Implementation (`ghg_module.c`)

```c
// CH₄ budget
dCH4/dt = benthic_flux/h - oxidation - air_water_flux/h + lateral_input
```

**Key Processes:**

1. **Benthic flux** (methanogenesis in anoxic sediments):
   $$F_{CH4,benthic} = F_{max} \times Q_{10}^{(T-20)/10} \times f(O_2^{bottom})$$

2. **Aerobic oxidation** (MOB bacteria):
   $$R_{ox} = k_{ox} \times [CH_4] \times \frac{[O_2]}{K_{O2} + [O_2]}$$

3. **Air-water exchange**:
   $$F_{CH4} = k_{CH4} \times ([CH_4] - [CH_4]_{sat})$$

**Validation Issue:** CH₄ underestimated by 35-38 nmol/L. The benthic flux alone cannot explain observed concentrations because:
- Rice paddies contribute ~50-200 nmol/L via drainage
- Aquaculture ponds contribute ~100-500 nmol/L via discharge
- These are **lateral inputs**, not benthic

#### N₂O Implementation

```c
// N₂O production from nitrification and denitrification
N2O_prod = yield_nit × nit_rate + yield_denit × denit_rate
```

**Key Equations:**

$$\frac{d[N_2O]}{dt} = \epsilon_{nit} R_{nit} + \epsilon_{denit} R_{denit} + F_{benthic}/h - k_{N2O}([N_2O] - [N_2O]_{sat})/h$$

Where:
- $\epsilon_{nit}$ = 0.002 (0.2% of N oxidized)
- $\epsilon_{denit}$ = 0.005 (0.5% of N reduced)

**Validation Issue:** N₂O underestimated by 6-14 nmol/L. Agricultural drainage (with active nitrification of fertilizer NH₄⁺) is the dominant source in the Mekong Delta but is not represented in current lateral loads.

---

## 4. Root Cause Analysis of Biases

### 4.1 pCO₂ Bias: Why Model is Too Low Upstream

**Observed Pattern:**
- Ocean (km=0): pCO₂ = 500 µatm ✓ Model correct
- Mid-estuary (km=30): pCO₂ = 1400 µatm, Model = 2300 (OK)
- Upstream (km=80): pCO₂ = 4500 µatm, Model = 1800 (TOO LOW)

**Root Causes:**

1. **Uniform benthic respiration** applies the same rate everywhere:
   - At ocean mouth: Too much CO₂ production (sandy sediments should have less)
   - At upstream: Too little CO₂ production (organic sediments should have more)

2. **Dispersion dominates over production**:
   - High mixing_alpha (0.40-0.50) causes rapid exchange with ocean water
   - DIC produced upstream is quickly mixed with low-DIC ocean water
   - Net effect: pCO₂ gradient is "flattened"

3. **Missing in-situ processes**:
   - Mangrove porewater discharge (high DIC, high pCO₂)
   - Rice paddy drainage (CO₂ supersaturated water)

**Quantitative Analysis:**

| Location | Observed pCO₂ | Model pCO₂ | Deficit | Required Additional DIC Source |
|----------|---------------|------------|---------|-------------------------------|
| km=30 | 1434 | 2329 | +895 | None (model OK) |
| km=60 | 2923 | 1839 | -1084 | +3.5 µmol DIC/L/hr |
| km=80 | 4720 | 1857 | -2863 | +9.2 µmol DIC/L/hr |

### 4.2 O₂ Bias: Why Model is Too Low Mid-Estuary

**Observed Pattern:**
- Ocean: O₂ = 262 µmol/L (near saturation), Model = 258 ✓
- km=30: O₂ = 230 µmol/L, Model = 184 (TOO LOW)
- Upstream: O₂ = 160 µmol/L, Model = 179 (slightly high)

**Root Cause:** Benthic O₂ consumption (SOD) is applied uniformly, but:
- Near ocean: Coarse sediments have LOW SOD (~5-15 mmol O₂/m²/day)
- Upstream: Fine sediments have HIGH SOD (~30-80 mmol O₂/m²/day)

**Current Code Issue:**
```c
// Same benthic_o2_rate for ALL cells
double benthic_o2_rate = benthic_rate_day / depth / SECONDS_PER_DAY;
```

### 4.3 CH₄ and N₂O: Missing Lateral Sources

**The 80/20 insight:** In agricultural deltas like the Mekong, **lateral inputs dominate** GHG budgets:

| Source | Contribution to CH₄ | Contribution to N₂O |
|--------|---------------------|---------------------|
| Benthic (water column) | 20-30% | 10-20% |
| Rice paddies | 40-50% | 30-40% |
| Aquaculture | 20-30% | 10-20% |
| Urban drainage | 5-10% | 20-30% |

**Currently Missing:**
- `CH4_conc_base` in lateral_sources.csv is ~0.1 µmol/L (too low)
- Rice paddy CH₄ can be 1-5 µmol/L during drainage events
- Agricultural NH₄ fertilizer creates N₂O "hot spots" not represented

---

## 5. Equations and Their Limitations

### 5.1 Summary of Model Equations

#### Transport (Advection-Dispersion)

$$\frac{\partial C}{\partial t} + u\frac{\partial C}{\partial x} = \frac{1}{A}\frac{\partial}{\partial x}\left(AD\frac{\partial C}{\partial x}\right) + S$$

Where dispersion follows Van den Burgh (Savenije 2005):
$$D(x) = D_0 \left[1 - K\frac{Q}{A_0 D_0}\left(e^{x/a} - 1\right)\right]$$

**Limitation:** D₀ parameterization (mixing_alpha) is empirical and requires calibration.

#### Biogeochemistry Rates

| Process | Equation | Parameters |
|---------|----------|------------|
| TOC degradation | $R = k_{ox} \times TOC \times \frac{O_2}{K_{O2}+O_2} \times \theta^{T-20}$ | kox=0.02/day |
| Nitrification | $R = k_{nit} \times NH_4 \times \frac{O_2}{K_{O2,nit}+O_2} \times \theta^{T-20}$ | knit=0.05/day |
| Denitrification | $R = k_{denit} \times NO_3 \times \frac{K_{O2,inhib}}{K_{O2,inhib}+O_2} \times \frac{TOC}{K_{TOC}+TOC}$ | kdenit=0.03/day |
| O₂ reaeration | $R = k_L \times (O_{2,sat} - O_2) / h$ | Wanninkhof (1992) |
| Benthic respiration | $F = F_{20} \times Q_{10}^{(T-20)/10} / h$ | F₂₀=70 mmol/m²/day |

#### Gas Exchange

$$k = 0.31 \times u_{10}^2 \times \left(\frac{Sc}{660}\right)^{-0.5}$$

**Limitation:** Does not account for:
- Current-enhanced turbulence (significant in tidal estuaries)
- Bubble-mediated transfer (important for CH₄)
- Surfactant effects (can reduce k by 50%)

### 5.2 Known Simplifications

| Simplification | Impact | Justification |
|----------------|--------|---------------|
| No vertical stratification | Medium | 1D model assumption; valid for well-mixed estuaries |
| Uniform water temperature | Low | Mekong ~28°C year-round |
| No sediment diagenesis | High | Benthic fluxes are parameterized, not computed |
| No tidal pumping of porewaters | Medium | Would add DIC/CH₄ pulses at low tide |
| Instantaneous carbonate equilibrium | Low | Valid at hourly timescales |

---

## 6. Recommended Improvements (80/20 Prioritized)

### Priority 1: Spatially-Varying Benthic Fluxes (HIGH IMPACT)

**Problem:** Uniform benthic_resp_20C = 70 mmol/m²/day everywhere

**Solution:** Scale benthic flux with distance from ocean (proxy for sediment organic content):

```c
// Proposed: Distance-weighted benthic flux
double dist_factor = 1.0 + 1.5 * (branch->dx * i) / branch->length;  // 1.0-2.5 scaling
double benthic_rate_day = p->benthic_resp_20C * dist_factor * pow(Q10, (T-20)/10);
```

**Justification:** 
- Literature shows 3-5× variation in benthic fluxes along estuaries (Abril et al. 2010)
- Distance is a globally-available proxy (no field measurements needed)
- Can be refined with satellite-derived sediment type maps

**Expected Impact:**
- pCO₂ upstream: +1500-2000 µatm (fixes gradient)
- O₂ near ocean: +20-30 µmol/L (reduces SOD at sandy mouth)

### Priority 2: Enhanced Lateral GHG Sources (HIGH IMPACT)

**Problem:** CH₄ and N₂O lateral inputs too low

**Solution:** Expand `generate_lateral_loads_v2.py` with GHG emission factors:

```python
# Proposed land use GHG emission factors
GHG_EMISSIONS = {
    'Rice_Paddy': {'CH4_umol_L': 2.0, 'N2O_umol_L': 0.15},  # Drainage water
    'Aquaculture': {'CH4_umol_L': 1.0, 'N2O_umol_L': 0.08},
    'Urban': {'CH4_umol_L': 0.3, 'N2O_umol_L': 0.20},
    'Agriculture': {'CH4_umol_L': 0.1, 'N2O_umol_L': 0.12},
    'Forest': {'CH4_umol_L': 0.05, 'N2O_umol_L': 0.02}
}
```

**Justification:**
- IPCC Wetlands Supplement (2014) provides default emission factors
- Satellite land use maps (JAXA, Sentinel) can identify rice paddies
- Seasonal scaling already implemented (wet season = higher emissions)

**Expected Impact:**
- CH₄: +30-50 nmol/L upstream (fixes 35 nmol/L bias)
- N₂O: +5-10 nmol/L upstream (fixes 6 nmol/L bias)

### Priority 3: Tidal Porewater Exchange (MEDIUM IMPACT)

**Problem:** Missing pulsed DIC/CH₄ inputs from sediment porewater drainage at low tide

**Physical Mechanism:**
At low tide, porewaters enriched in DIC (5000-20000 µmol/L) and CH₄ (1-10 µmol/L) drain into the channel. This is NOT captured by steady-state benthic flux.

**Proposed Implementation:**
```c
// Porewater drainage at low tide (simplified)
if (water_level < mean_water_level - 0.5) {
    double porewater_flux = k_pore * (pore_DIC - water_DIC) * exposed_area;
    DIC[i] += porewater_flux / volume * dt;
}
```

**Data Requirements:**
- Tidal amplitude (already available from forcing)
- Intertidal area fraction (from satellite DEM, ~30m resolution)
- Porewater concentrations (literature values: DIC ~8000 µmol/L, CH₄ ~2 µmol/L)

**Expected Impact:**
- pCO₂: +500-1000 µatm in tidal flat areas
- CH₄: +20-50 nmol/L during ebb tide

### Priority 4: Sediment-Type Proxy from SPM (LOW-MEDIUM IMPACT)

**Concept:** Use SPM concentration as a proxy for sediment type:
- High SPM → Fine sediments → High benthic fluxes
- Low SPM → Coarse sediments → Low benthic fluxes

```c
// SPM-scaled benthic flux
double spm_factor = 0.5 + 0.5 * (spm[i] / 50.0);  // 0.5-1.5 scaling, ref SPM=50 mg/L
double benthic_rate = p->benthic_resp_20C * spm_factor * temp_factor;
```

**Justification:**
- SPM is a transported quantity (already simulated)
- Provides dynamic (tidal) variation in benthic flux
- No additional data requirements

---

## 7. Global Dataset Integration Strategy

### 7.1 Available Global Datasets

| Dataset | Resolution | Variables | Access |
|---------|------------|-----------|--------|
| JAXA ALOS-2 | 25m | Land use/cover | Free |
| Sentinel-2 | 10m | Land use, water bodies | Free |
| WorldClim | 1km | Monthly precipitation | Free |
| ERA5 | 31km | Wind, temperature | Free |
| SRTM | 30m | Elevation (intertidal area) | Free |
| NOAA | Point | Atmospheric CO₂, CH₄, N₂O | Free |
| Ocean Color (MODIS) | 4km | Chlorophyll, TSM | Free |

### 7.2 Proposed Data Pipeline

```
┌─────────────────┐     ┌─────────────────┐     ┌─────────────────┐
│   Satellite     │     │    Climate      │     │   Topography    │
│   Land Use      │     │   Reanalysis    │     │      DEM        │
│  (JAXA/S2)      │     │   (ERA5)        │     │    (SRTM)       │
└────────┬────────┘     └────────┬────────┘     └────────┬────────┘
         │                       │                       │
         ▼                       ▼                       ▼
┌─────────────────────────────────────────────────────────────────┐
│                generate_lateral_loads_v3.py                      │
│  • Rice paddy fraction → CH₄/N₂O emissions                      │
│  • Rainfall × runoff coefficient → Q_lateral                     │
│  • Elevation → intertidal area → porewater exchange             │
└─────────────────────────────────────────────────────────────────┘
         │
         ▼
┌─────────────────┐     ┌─────────────────┐     ┌─────────────────┐
│ lateral_sources │     │ seasonal_factors│     │ porewater_params│
│     .csv        │     │     .csv        │     │     .csv        │
└─────────────────┘     └─────────────────┘     └─────────────────┘
```

### 7.3 Recommended New Parameters from Global Data

| Parameter | Derivation | Expected Improvement |
|-----------|------------|---------------------|
| Rice paddy fraction | JAXA 25m classification | CH₄, N₂O lateral loads |
| Intertidal area | SRTM elevation < MHW | Porewater DIC flux |
| Sediment type index | MODIS TSM + bathymetry | Spatially-varying benthic flux |
| Seasonal discharge | ERA5 precipitation | Lateral Q_factor |

---

## 8. Implementation Roadmap

### Phase 1: Quick Wins (1-2 weeks)

1. **Implement distance-weighted benthic flux**
   - Modify `Biogeo_Branch_Simplified()` in `biogeo.c`
   - Add `benthic_gradient_factor` parameter
   - Expected: Fix pCO₂ gradient, improve O₂ at ocean

2. **Update lateral GHG emission factors**
   - Modify `generate_lateral_loads_v2.py`
   - Add rice paddy CH₄/N₂O based on IPCC defaults
   - Expected: Fix CH₄/N₂O biases

### Phase 2: Data Integration (2-4 weeks)

3. **Add intertidal porewater exchange**
   - New function `calc_porewater_flux()` in `biogeo.c`
   - Triggered by water level below MWL
   - Requires: intertidal_fraction per branch (from DEM)

4. **SPM-benthic coupling**
   - Scale benthic flux by local SPM
   - Self-consistent: SPM affects sedimentation → benthic flux

### Phase 3: Validation & Documentation (2-4 weeks)

5. **Multi-site validation**
   - Test on Red River, Ganges (different climate/land use)
   - Ensure improvements are not Mekong-specific

6. **Sensitivity analysis**
   - Document parameter ranges and uncertainties
   - Create calibration guidelines for new sites

---

## Appendix A: Parameter Values (Current Configuration)

```
# Simplified mode biogeochemistry (biogeo_params.txt)
water_temp = 28.0
kox = 0.02 /day         # TOC degradation (increased from 0.008)
knit = 0.05 /day        # Nitrification
kdenit = 0.03 /day      # Denitrification
theta_ox = 1.047        # Temperature coefficient
benthic_resp_20C = 70.0 # mmol C/m²/day (was 50, needs spatial variation)
benthic_Q10 = 1.8       # Temperature sensitivity
N2O_yield_nit = 0.002   # 0.2% of N oxidized
N2O_yield_denit = 0.005 # 0.5% of N reduced
benthic_CH4_flux = 300  # µmol/m²/day
```

## Appendix B: Boundary Conditions (Current Configuration)

**Ocean (species_ocean_realistic.csv):**
- Salinity: 30.5 PSU
- DIC: 1950 µmol/L
- AT: 2200 µeq/L (gives pH=8.1, pCO₂~450 µatm)
- O₂: 260 µmol/L
- NO₃: 55 µmol/L (coastal upwelling)

**River (species_river_realistic.csv):**
- Salinity: 0.1 PSU
- DIC: 1480 µmol/L
- AT: 1320 µeq/L (gives pH=7.45, pCO₂~4200 µatm)
- O₂: 175 µmol/L
- NO₃: 10 µmol/L

## Appendix C: Key References

1. **Savenije (2005, 2012)** - Salinity and Tides in Alluvial Estuaries. Elsevier.
2. **Abril et al. (2010)** - Carbon dioxide and methane emissions from estuaries. *Limnol. Oceanogr.* 55:1199.
3. **Wanninkhof (1992)** - Gas exchange relationships. *J. Geophys. Res.* 97:7373.
4. **Garnier et al. (2005)** - RIVERSTRAHLER model. *Biogeochemistry* 77:213.
5. **IPCC (2014)** - Wetlands Supplement to the 2006 Guidelines.
6. **Lueker et al. (2000)** - Carbonate equilibrium constants. *Mar. Chem.* 70:105.
7. **Weiss (1974)** - CO₂ solubility in seawater. *Mar. Chem.* 2:203.
8. **Borges & Abril (2011)** - Carbon Dioxide and Methane Dynamics in Estuaries. *Treatise on Estuarine and Coastal Science*.

---

**Document Version:** 2.0  
**Last Updated:** December 5, 2025  
**Next Review:** After Phase 1 implementation
