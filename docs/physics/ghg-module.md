# Greenhouse Gas (GHG) Module

## Overview

The C-GEM GHG module simulates the dynamics of dissolved greenhouse gases in estuarine waters:

- **CO₂** (Carbon dioxide) - via carbonate chemistry
- **CH₄** (Methane) - production, oxidation, gas exchange
- **N₂O** (Nitrous oxide) - from nitrification and denitrification

This document explains the underlying physics, parameterization, and known limitations.

---

## Species Indices

| Species | Index | Unit | Transport | Description |
|---------|-------|------|-----------|-------------|
| N2O | 28 | nmol/L | Yes | Nitrous oxide |
| CH4 | 29 | nmol/L | Yes | Methane |
| CO2 | 13 | µmol/L | No (diagnostic) | Computed from DIC/TA |
| pCO2 | 12 | µatm | No (diagnostic) | Computed from DIC/TA |

---

## CH₄ (Methane) Dynamics

### ⚠️ Important: Water Column vs Benthic Sources

In oxic estuaries like the Mekong Delta (O₂ > 30% saturation), **CH₄ does NOT form in the water column**. Methanogenesis requires strict anoxia (O₂ < 1 µM).

**Primary CH₄ Sources in Tropical Estuaries:**

| Source | Mechanism | Typical Flux |
|--------|-----------|--------------|
| **Benthic sediment** | Diffusion from anoxic mud | 50-500 µmol/m²/day |
| **Rice paddy drainage** | Flooded anoxic soils | High during wet season |
| **Aquaculture ponds** | Anoxic pond bottom sediments | Variable |
| **Groundwater seepage** | CH₄-rich groundwater | Localized hotspots |

### Current Model Implementation

The model includes:

1. **Benthic CH₄ flux** (primary source):
   ```c
   F_benthic = benthic_CH4_flux / (H * 1000.0);  // µmol/m²/day → µM/day
   ```
   Default: 150 µmol/m²/day (adjustable in `biogeo_params.txt`)

2. **CH₄ oxidation** (sink):
   ```c
   R_oxid = k_CH4_oxid * CH4 * O2 / (K_O2_CH4 + O2);
   ```
   
3. **Air-water gas exchange** (typically emission):
   ```c
   F_gas = k_CH4 * (CH4 - CH4_sat);  // Usually CH4 > CH4_sat → emission
   ```

### Known Limitations

1. **No lateral CH₄ inputs**: Rice paddies, aquaculture, and urban drainage are major CH₄ sources not yet in the lateral loads system

2. **Constant benthic flux**: Should vary with:
   - Organic matter deposition
   - Temperature (Q10 ≈ 2-3)
   - Sediment type (higher in fine muds)

3. **No ebullition**: Direct bubble flux can dominate in shallow, organic-rich areas

### Validation Comparison (Mekong March 2025)

| Location | Observed CH₄ | Model (v1.2.0) | Notes |
|----------|--------------|----------------|-------|
| Mouth (0-10 km) | 40-80 nmol/L | 40 nmol/L | Now correct (was 0) |
| Mid-estuary (30-50 km) | 50-150 nmol/L | 1-10 nmol/L | Underestimated |
| Upstream (60-80 km) | 100-200 nmol/L | Near 0 | Major gap |

**Key gap**: The model underestimates upstream CH₄ because it lacks:
- Rice paddy drainage inputs
- Urban sewage CH₄
- Enhanced benthic flux in freshwater zones

### Recommendations for Improvement

1. Add CH₄ to lateral sources system:
   ```csv
   Branch,Segment_Index,...,CH4_conc_base_nM
   Ham_Luong,45,...,500
   ```

2. Make benthic flux spatially variable:
   ```c
   F_benthic = base_flux * (1 + organic_deposition_factor) * f(temperature);
   ```

3. Add rice paddy drainage as point sources during wet season

---

## N₂O (Nitrous Oxide) Dynamics

### Production Pathways

N₂O is produced through two distinct microbial pathways:

#### 1. Nitrification Pathway (oxic)

During NH₄ → NO₃ oxidation, a small fraction escapes as N₂O:

```c
R_N2O_nit = f_N2O * k_nit * NH4 * O2 / (K_O2_nit + O2);
```

Where:
- `f_N2O` = yield factor (0.001-0.01, typically 0.002)
- `k_nit` = nitrification rate constant

This pathway is **enhanced at low O₂** (2-6 mg/L) because incomplete oxidation increases N₂O yield.

#### 2. Denitrification Pathway (micro-anoxic)

During NO₃ → N₂ reduction, N₂O is an intermediate:

```c
R_N2O_denit = g_N2O * k_denit * NO3 * TOC * K_O2_denit / (K_O2_denit + O2);
```

Where:
- `g_N2O` = yield factor (0.01-0.05)
- Suppressed when O₂ is high, maximized in transition zone

### N₂O Consumption

N₂O is reduced to N₂ under anoxia:

```c
R_N2O_red = k_N2O_red * N2O * K_O2_N2O / (K_O2_N2O + O2);
```

### Gas Exchange

N₂O is typically supersaturated in estuaries:

```c
N2O_sat = 8.0 * f(T, S);  // nmol/L at equilibrium with atmosphere
F_N2O = k_N2O * (N2O - N2O_sat);
```

### Validation Comparison (Mekong March 2025)

| Location | Observed N₂O | Model (v1.2.0) | Notes |
|----------|--------------|----------------|-------|
| Mouth | 8-15 nmol/L | 8 nmol/L | Correct (ocean value) |
| Mid-estuary | 10-25 nmol/L | 5-15 nmol/L | Reasonable |
| Upstream | 15-175 nmol/L | <5 nmol/L | Underestimated |

**Key issues**:
- Model underestimates N₂O in upstream freshwater zones
- Nitrification-sourced N₂O requires more NH₄ input (lateral loads)
- Urban sewage is a major N₂O source not captured

---

## CO₂ / pCO₂ Dynamics

### Carbonate Chemistry

CO₂ is computed diagnostically from DIC and Total Alkalinity (TA) using the carbonate equilibrium system:

```c
// Carbonic acid equilibria
H2CO3* ⇌ H⁺ + HCO3⁻     (K1)
HCO3⁻  ⇌ H⁺ + CO3²⁻     (K2)

// CO2 concentration
CO2 = DIC / (1 + K1/H + K1*K2/H²)

// pCO2 from Henry's law
pCO2 = CO2 / K_H(T,S)
```

### Sources and Sinks of DIC

| Process | Effect on DIC | Effect on TA |
|---------|---------------|--------------|
| Respiration | +1 | 0 |
| Photosynthesis | -1 | 0 |
| Nitrification | 0 | -2 |
| Denitrification | -0.8 | +0.8 |
| CaCO₃ dissolution | +1 | +2 |
| Air-water CO₂ flux | ±1 | 0 |

### Validation Comparison (Mekong March 2025)

| Location | Observed pCO₂ | Model (Reactions ON) | Notes |
|----------|---------------|---------------------|-------|
| Mouth | 400-600 µatm | 600-800 µatm | Slight overestimate |
| Mid-estuary | 1000-3000 µatm | 2000-3500 µatm | Too high |
| Upstream | 3000-6000 µatm | 4000-5000 µatm | Reasonable |

**Issues**: 
- pCO₂ too high because respiration rates may be overestimated
- pH consequently too low (7.3-7.5 vs observed 7.5-8.1)

---

## Common Issues and Solutions

### Issue 1: GHG Species = 0 at Ocean Mouth

**Symptom**: CH₄, N₂O show 0 nmol/L at km 0-5

**Cause**: Ocean boundary concentration bug (fixed in v1.2.0)

**Solution**: Update to v1.2.0 or apply the boundary condition fix

### Issue 2: CH₄/N₂O Too Low Upstream

**Symptom**: Model < 50 nmol/L while observed > 100 nmol/L

**Cause**: Missing lateral inputs from:
- Rice paddy drainage
- Urban/industrial discharge
- Aquaculture effluent

**Solution**: 
1. Add CH₄/N₂O columns to `lateral_sources.csv`
2. Increase benthic CH₄ flux parameter
3. Add agricultural point sources

### Issue 3: pCO₂ Too High / pH Too Low

**Symptom**: pCO₂ > 5000 µatm, pH < 7.4

**Cause**: Excessive respiration rates relative to photosynthesis

**Solution**:
1. Reduce `k_resp_base` in biogeo_params.txt
2. Increase light availability for photosynthesis
3. Check TOC boundary values (may be driving too much respiration)

### Issue 4: O₂ Too Stable (Not Dropping Upstream)

**Symptom**: O₂ stays ~220 µM from ocean to upstream

**Cause**: 
- Reaeration too fast
- Respiration too slow
- TOC too low to drive oxygen demand

**Solution**:
1. Check TOC boundary values (should be 150-200 µM upstream)
2. Reduce reaeration coefficient
3. Increase respiration rate or TOC degradation rate

---

## Key Parameters (`biogeo_params.txt`)

| Parameter | Default | Unit | Description |
|-----------|---------|------|-------------|
| `benthic_CH4_flux` | 150.0 | µmol/m²/day | CH₄ diffusion from sediment |
| `k_CH4_oxid` | 0.05 | 1/day | CH₄ oxidation rate |
| `f_N2O_nit` | 0.002 | - | N₂O yield from nitrification |
| `g_N2O_denit` | 0.02 | - | N₂O yield from denitrification |
| `k_nit` | 0.1 | 1/day | Nitrification rate |
| `k_denit` | 0.05 | 1/day | Denitrification rate |
| `k_resp` | 0.03 | 1/day | TOC respiration rate |

---

## References

- **CH₄ in estuaries**: Borges & Abril (2011), Biogeosciences
- **N₂O from nitrification**: Goreau et al. (1980), Applied & Environmental Microbiology
- **N₂O from denitrification**: Seitzinger & Kroeze (1998), Nutrient Cycling
- **Carbonate chemistry**: Zeebe & Wolf-Gladrow (2001), CO2 in Seawater
- **Gas exchange**: Wanninkhof (1992), Journal of Geophysical Research
