# Biogeochemistry Module

## Overview

The C-GEM biogeochemistry module implements **C-RIVE** (Unified RIVE v1.0), a comprehensive process-based model for carbon, nutrient, and oxygen cycling. Key features include:

- **Complete carbonate chemistry** (DIC, TA, pH, pCO₂, CO₂ flux)
- **Greenhouse gas dynamics** (CO₂, CH₄, N₂O)
- **2-step nitrification** (NH₄ → NO₂ → NO₃)
- **RK4 adaptive solver** for numerical stability

---

## Species and State Variables

### Core Biogeochemical Species

| Index | Symbol | Name | Unit |
|-------|--------|------|------|
| 0 | SAL | Salinity | PSU |
| 1 | PHY1 | Diatoms | µg C/L |
| 2 | PHY2 | Green algae | µg C/L |
| 3 | DSi | Dissolved silica | µmol/L |
| 4 | NO₃ | Nitrate | µmol N/L |
| 5 | NH₄ | Ammonium | µmol N/L |
| 6 | PO₄ | Phosphate | µmol P/L |
| 7 | O₂ | Dissolved oxygen | µmol/L |
| 8 | TOC | Total organic carbon | µmol C/L |
| 9 | SPM | Suspended matter | mg/L |

### Carbonate System

| Index | Symbol | Name | Unit |
|-------|--------|------|------|
| 10 | DIC | Dissolved inorganic C | µmol/L |
| 11 | TA | Total alkalinity | µeq/L |
| 12 | pCO₂ | CO₂ partial pressure | µatm |
| 13 | CO₂ | Dissolved CO₂ | µmol/L |
| 14 | pH | Acidity | - |

### Greenhouse Gases

| Index | Symbol | Name | Unit |
|-------|--------|------|------|
| 27 | NO₂ | Nitrite | µmol N/L |
| 28 | N₂O | Nitrous oxide | nmol N/L |
| 29 | CH₄ | Methane | µmol C/L |

---

## Phytoplankton Module

### Growth Model

$$\mu = \mu_{max} \cdot f_T \cdot f_L \cdot f_N$$

### Light Limitation

Depth-integrated photosynthesis:

$$f_L = \frac{1}{k_d H} \left[e^{-k_d H \cdot e^{-I_0/I_k}} - e^{-I_0/I_k}\right]$$

### Nutrient Limitation

**Diatoms (PHY1)** - Liebig's law:

$$f_N = \min\left(\frac{DSi}{DSi + K_{Si}}, \frac{DIN}{DIN + K_N}, \frac{PO_4}{PO_4 + K_P}\right)$$

### Salinity Stress

Freshwater diatoms in tropical estuaries experience salinity stress:

$$f_{sal} = 1 + k_{sal} \cdot \max(0, S - S_{thresh})$$

Applied to mortality:
$$M = k_M \cdot f_{sal} \cdot PHY$$

---

## Nitrogen Cycling

### Two-Step Nitrification

**Step 1**: NH₄ → NO₂ (Ammonia oxidation)
$$R_{nit1} = k_{nit1} \cdot f_T \cdot \frac{O_2}{K_{O2} + O_2} \cdot NH_4$$

**Step 2**: NO₂ → NO₃ (Nitrite oxidation)
$$R_{nit2} = k_{nit2} \cdot f_T \cdot \frac{O_2}{K_{O2} + O_2} \cdot NO_2$$

### Denitrification

$$R_{denit} = k_{denit} \cdot f_T \cdot \frac{K_{O2,denit}}{K_{O2,denit} + O_2} \cdot NO_3$$

---

## Oxygen Dynamics

### Production/Consumption

$$\frac{dO_2}{dt} = GPP - R_{resp} - R_{nit} - R_{deg} + F_{ex}$$

### Air-Water Exchange

$$F_{ex} = k_w \cdot (O_{2,sat} - O_2)$$

Gas transfer velocity:
$$k_w = k_{600} \cdot \left(\frac{Sc}{600}\right)^{-0.5}$$

### O₂ Saturation

$$O_{2,sat} = f(T, S)$$

Using Weiss (1970) equation.

---

## Carbonate Chemistry

### Equilibrium System

$$DIC = [CO_2] + [HCO_3^-] + [CO_3^{2-}]$$
$$TA = [HCO_3^-] + 2[CO_3^{2-}] + [OH^-] - [H^+]$$

### pH Calculation

Newton-Raphson iteration to solve carbonate equilibria.

### CO₂ Flux

$$F_{CO2} = k_w \cdot K_0 \cdot (pCO_{2,water} - pCO_{2,atm})$$

---

## Greenhouse Gas Module

### N₂O Production

From nitrification:
$$R_{N2O,nit} = f_{N2O,nit} \cdot R_{nit}$$

From denitrification:
$$R_{N2O,denit} = f_{N2O,denit} \cdot R_{denit}$$

Typical yields:
- Nitrification: 0.3-0.5%
- Denitrification: 0.1-2.0%

### CH₄ Dynamics

Production in sediments:
$$R_{CH4,prod} = k_{CH4} \cdot f_T \cdot TOC$$

Oxidation in water column:
$$R_{CH4,ox} = k_{ox,CH4} \cdot \frac{O_2}{K_{O2} + O_2} \cdot CH_4$$

Ebullition:
$$R_{CH4,eb} = k_{eb} \cdot \max(0, CH_4 - CH_{4,sat})$$

---

## Parameters

### Phytoplankton

| Parameter | Symbol | Default | Unit |
|-----------|--------|---------|------|
| Max growth rate | $\mu_{max}$ | 2.0 | d⁻¹ |
| Mortality rate | $k_M$ | 0.08 | d⁻¹ |
| Light saturation | $I_k$ | 100 | W/m² |

### Bacteria

| Parameter | Symbol | Default | Unit |
|-----------|--------|---------|------|
| Aerobic degradation | $k_{ox}$ | 0.08 | d⁻¹ |
| Nitrification step 1 | $k_{nit1}$ | 0.10 | d⁻¹ |
| Nitrification step 2 | $k_{nit2}$ | 0.15 | d⁻¹ |
| Denitrification | $k_{denit}$ | 0.03 | d⁻¹ |

---

## References

1. Wang, S. et al. (2018). Time-dependent global sensitivity analysis of C-RIVE. *Water Research*, 144, 341-355.
2. Hasanyar, M. et al. (2022). Unified RIVE v1.0. *Biogeosciences*.
3. Weiss, R.F. (1970). Solubility of gases in seawater. *Deep-Sea Research*.
4. Wanninkhof, R. (1992). Gas exchange. *J. Geophys. Res.*
