# Biogeochemistry Module

## Overview

The C-GEM biogeochemistry module implements **C-RIVE** (Unified RIVE v1.0), a comprehensive process-based model for carbon, nutrient, and oxygen cycling in aquatic systems. Key features include:

- **Complete carbonate chemistry** (DIC, TA, pH, pCO2, CO2 flux)
- **Greenhouse gas dynamics** (CO2, CH4, N2O)
- **Multi-pool organic matter** (6 pools with different reactivity)
- **Explicit microbial communities** (heterotrophic bacteria, nitrifiers)
- **2-step nitrification** (NH4 → NO2 → NO3)
- **RK4 adaptive solver** for numerical stability

---

## Species and State Variables

### Core Biogeochemical Species

| Index | Symbol | Name | Unit | Description |
|-------|--------|------|------|-------------|
| 0 | SAL | Salinity | PSU | Conservative tracer |
| 1 | PHY1 | Diatoms | µg C/L | Siliceous phytoplankton |
| 2 | PHY2 | Green algae | µg C/L | Non-siliceous phytoplankton |
| 3 | DSi | Dissolved silica | µmol/L | Nutrient for diatoms |
| 4 | NO3 | Nitrate | µmol N/L | Oxidized nitrogen |
| 5 | NH4 | Ammonium | µmol N/L | Reduced nitrogen |
| 6 | PO4 | Phosphate | µmol P/L | Dissolved reactive P |
| 7 | O2 | Dissolved oxygen | µmol/L | ~312 µmol/L at saturation |
| 9 | SPM | Suspended matter | mg/L | Affects light, adsorption |

### Carbonate System

| Index | Symbol | Name | Unit | Description |
|-------|--------|------|------|-------------|
| 10 | DIC | Dissolved inorganic C | µmol/L | CO2 + HCO3⁻ + CO3²⁻ |
| 11 | TA | Total alkalinity | µeq/L | Acid-neutralizing capacity |
| 12 | pCO2 | CO2 partial pressure | µatm | *Diagnostic* |
| 13 | CO2 | Dissolved CO2 | µmol/L | *Diagnostic* |
| 14 | pH | Acidity | - | *Diagnostic* |

### RIVE Organic Matter Pools

| Index | Symbol | Name | Reactivity |
|-------|--------|------|------------|
| 17 | HD1 | Labile dissolved OC | Fast (k = 0.5-1.0 /day) |
| 18 | HD2 | Semi-labile dissolved OC | Medium (k = 0.1-0.3 /day) |
| 19 | HD3 | Refractory dissolved OC | Slow (k = 0.01-0.05 /day) |
| 20 | HP1 | Labile particulate OC | Fast |
| 21 | HP2 | Semi-labile particulate OC | Medium |
| 22 | HP3 | Refractory particulate OC | Slow |

### RIVE Microbial Pools

| Index | Symbol | Name | Description |
|-------|--------|------|-------------|
| 23 | BAG | Attached bacteria | Large, particle-associated |
| 24 | BAP | Free bacteria | Small, planktonic |
| 25 | PIP | Particulate inorganic P | Adsorbed on SPM |
| 26 | DSS | Dissolved simple substrates | Bacterial food source |

### Greenhouse Gas Species

| Index | Symbol | Name | Unit | Source |
|-------|--------|------|------|--------|
| 27 | NO2 | Nitrite | µmol N/L | Nitrification intermediate |
| 28 | N2O | Nitrous oxide | nmol N/L | GHG from nitrif/denit |
| 29 | CH4 | Methane | µmol C/L | GHG from methanogenesis |

---

## Phytoplankton Module

### Growth Model

Phytoplankton growth follows the **Eilers-Peeters photosynthesis model**:

$$\mu = \mu_{max} \cdot f_T \cdot f_L \cdot f_N$$

Where:
- $\mu_{max}$ = Maximum growth rate [/day]
- $f_T$ = Temperature function
- $f_L$ = Light limitation
- $f_N$ = Nutrient limitation

### Light Limitation

Depth-integrated photosynthesis:

$$f_L = \frac{1}{k_d H} \left[e^{-k_d H \cdot e^{-I_0/I_k}} - e^{-I_0/I_k}\right]$$

Where:
- $k_d$ = Light extinction [/m]
- $H$ = Water depth [m]
- $I_0$ = Surface irradiance [W/m²]
- $I_k$ = Light saturation parameter [W/m²]

### Light Extinction

$$k_d = k_{d,water} + k_{d,SPM} \cdot SPM + k_{d,PHY} \cdot (PHY1 + PHY2)$$

Typical values:
- $k_{d,water}$ = 0.1 /m (clear water)
- $k_{d,SPM}$ = 0.01 /m per mg/L
- $k_{d,PHY}$ = 0.005 /m per µg C/L

### Nutrient Limitation

**Diatoms (PHY1)** - Liebig's law:

$$f_N = \min\left(\frac{DSi}{DSi + K_{Si}}, \frac{DIN}{DIN + K_N}, \frac{PO4}{PO4 + K_P}\right)$$

**Non-siliceous (PHY2)**:

$$f_N = \min\left(\frac{DIN}{DIN + K_N}, \frac{PO4}{PO4 + K_P}\right)$$

Where $DIN = NO3 + NH4$.

### NH4/NO3 Preference

Phytoplankton preferentially uptake NH4 (energetically cheaper):

$$f_{NH4} = \frac{NH4}{K_{NH4,pref} + NH4}$$

NH4 uptake fraction: $f_{NH4}$  
NO3 uptake fraction: $(1 - f_{NH4})$

### Loss Processes

- **Respiration**: $R = k_R \cdot f_T \cdot PHY$
- **Mortality**: $M = k_M \cdot PHY$ (feeds detritus)
- **Settling**: $S = w_s / H \cdot PHY$

---

## Organic Matter Dynamics (RIVE)

### Multi-Pool Conceptual Model

```
         River Input
              ↓
    ┌─────────┴─────────┐
    │   Phytoplankton   │
    │   Mortality       │
    └─────────┬─────────┘
              ↓
    ┌─────────┴─────────┐
    │  Particulate OC   │
    │  HP1 → HP2 → HP3  │
    └─────────┬─────────┘
              ↓ Hydrolysis
    ┌─────────┴─────────┐
    │   Dissolved OC    │
    │  HD1 → HD2 → HD3  │
    └─────────┬─────────┘
              ↓ Uptake
    ┌─────────┴─────────┐
    │     Bacteria      │
    │    BAG + BAP      │
    └─────────┬─────────┘
              ↓ Respiration
            CO2 + NH4 + PO4
```

### Hydrolysis Rates

$$\frac{dHD1}{dt} = k_{hydr,1} \cdot HD1 + k_{hydr,P1} \cdot HP1$$

Where hydrolysis converts:
- HP (particulate) → HD (dissolved)
- HD (high reactivity) → DSS (simple substrates)

### Temperature Dependence

All rates follow Arrhenius:

$$k(T) = k_{20} \cdot \theta^{(T-20)}$$

Where:
- $k_{20}$ = Rate at 20°C
- $\theta$ = Temperature coefficient (1.04-1.08)

---

## Bacterial Dynamics (RIVE)

### Two-Pool Bacteria Model

RIVE distinguishes:
- **BAG**: Large (gros), attached, slower growth, higher yield
- **BAP**: Small (petit), free, faster growth, lower yield

### Growth on Substrates

$$\mu_B = \mu_{max,B} \cdot f_T \cdot \frac{DSS}{DSS + K_S} \cdot \frac{O2}{O2 + K_{O2}}$$

Where DSS = dissolved simple substrates (from hydrolysis).

### Bacterial Respiration

$$R_B = (1 - Y) \cdot \mu_B \cdot B + k_{maint} \cdot B$$

Where:
- $Y$ = Growth yield (0.2-0.4)
- $k_{maint}$ = Maintenance respiration [/day]

### Mortality and Recycling

Dead bacteria become particulate OC:
$$\text{BAG mortality} → HP1$$
$$\text{BAP mortality} → HD1$$

---

## Nitrogen Cycle

### 2-Step Nitrification (from C-RIVE)

**Step 1 - Nitrosation (NH4 → NO2):**

$$R_{nit1} = k_{nit1} \cdot f_T \cdot \frac{NH4}{NH4 + K_{NH4}} \cdot \frac{O2}{O2 + K_{O2,nit}}$$

**Step 2 - Nitratation (NO2 → NO3):**

$$R_{nit2} = k_{nit2} \cdot f_T \cdot \frac{NO2}{NO2 + K_{NO2}} \cdot \frac{O2}{O2 + K_{O2,nit}}$$

Typical rates at 20°C:
- $k_{nit1}$ = 0.3 /day
- $k_{nit2}$ = 0.6 /day (faster, prevents NO2 accumulation)

### Denitrification

Under low O2 conditions:

$$R_{denit} = k_{denit} \cdot f_T \cdot \frac{OC}{OC + K_{OC}} \cdot \frac{NO3}{NO3 + K_{NO3}} \cdot \frac{K_{O2,inhib}}{O2 + K_{O2,inhib}}$$

Product: N2 (to atmosphere) + some N2O (GHG)

### N2O Production

From C-RIVE:

$$R_{N2O,nit} = \epsilon_{nit} \cdot R_{nit1}$$
$$R_{N2O,denit} = \epsilon_{denit} \cdot R_{denit}$$

Where:
- $\epsilon_{nit}$ ≈ 0.004 (0.4% of NH4 oxidized)
- $\epsilon_{denit}$ ≈ 0.01 (1% of NO3 reduced)

---

## Carbonate Chemistry (C-RIVE)

### DIC-TA-pH System

DIC (Dissolved Inorganic Carbon):
$$DIC = [CO2] + [HCO3^-] + [CO3^{2-}]$$

Total Alkalinity:
$$TA = [HCO3^-] + 2[CO3^{2-}] + [OH^-] - [H^+] + ...$$

### pH Calculation

The Cardano-Vieta analytical solution for [H⁺]:

$$[H^+]^3 + X[H^+]^2 + Y[H^+] + Z = 0$$

Where X, Y, Z depend on DIC, TA, K1, K2.

### CO2 Speciation

$$[CO2] = DIC \cdot \frac{[H^+]^2}{[H^+]^2 + K_1[H^+] + K_1 K_2}$$

$$[HCO3^-] = DIC \cdot \frac{K_1[H^+]}{[H^+]^2 + K_1[H^+] + K_1 K_2}$$

### Dissociation Constants

Temperature-dependent (Harned & Davis):

$$pK_1 = -126.34 + \frac{6320.8}{T} + 19.57 \ln(T)$$

### CO2 Air-Water Exchange

$$F_{CO2} = k_{CO2} \cdot (C_{sat} - [CO2])$$

Where:
- $k_{CO2}$ = Gas transfer velocity [m/s]
- $C_{sat}$ = Saturation concentration [µmol/L]

### k600 Parameterization

For tidal estuaries, C-GEM uses the **Abril et al. (2009)** formulation, which is specifically developed for macrotidal estuaries like the Scheldt, Gironde, Mekong, and Saigon rivers:

$$k_{600} = 1.0 + 1.719 \sqrt{v}$$ [cm/h]

Where $v$ is the current velocity [m/s].

This formulation accounts for:
- **Tidal currents** generating turbulence at the water surface
- **Typical estuarine mixing conditions**
- **Valid for large tidal estuaries** (width > 100 m, tidal range > 1 m)

**Alternative methods available:**

| Method | Formula | Application |
|--------|---------|-------------|
| Abril et al. (2009) | $k_{600} = 1.0 + 1.719\sqrt{v}$ | Tidal estuaries (default) |
| Borges et al. (2004) | $k_{600} = 1.0 + U_{10} + 0.25 U_{10}^2$ | Wind-driven, open estuaries |
| Raymond et al. (2012) | Strahler-dependent | Rivers and streams |

**Note:** The Raymond et al. (2012) Strahler-based formulation is designed for rivers, not estuaries. For large tidal systems like the Mekong Delta or Saigon River, the Abril et al. (2009) method is more appropriate as it was calibrated using measurements in similar environments.

### Schmidt Number Correction

$$k_{CO2} = k_{600} \cdot \sqrt{\frac{600}{Sc}}$$

Where Sc is the Schmidt number:
$$Sc = 1911 - 118T + 3.45T^2 - 0.041T^3$$

---

## Methane (CH4) Module

### Sources

**Methanogenesis** (anoxic sediments):
$$R_{CH4,prod} = k_{CH4,prod} \cdot f_T \cdot OC_{sed} \cdot f_{anoxic}$$

**Benthic flux**:
$$F_{CH4,sed} = F_{CH4,0} \cdot f_T / H$$

### Sinks

**Aerobic oxidation**:
$$R_{CH4,ox} = k_{CH4,ox} \cdot f_T \cdot \frac{CH4}{CH4 + K_{CH4}} \cdot \frac{O2}{O2 + K_{O2,CH4}}$$

**Air-water evasion**:
$$F_{CH4,air} = k_{CH4} \cdot (CH4 - CH4_{sat})$$

**Ebullition** (bubble release):
$$R_{CH4,ebul} = k_{ebul} \cdot \max(0, CH4 - CH4_{thresh})$$

Typical threshold: 50 µmol/L in shallow water.

---

## Oxygen Dynamics

### Sources

$$\frac{dO2}{dt}_{sources} = R_{NPP} + F_{O2,air}$$

Where NPP includes photoautotrophic oxygen production.

### Sinks

$$\frac{dO2}{dt}_{sinks} = -R_{OC,ox} - 1.5 R_{nit1} - 0.5 R_{nit2} - R_{CH4,ox}$$

Stoichiometry:
- OC oxidation: 1 mol O2 per mol C
- Nitrosation: 1.5 mol O2 per mol NH4
- Nitratation: 0.5 mol O2 per mol NO2
- CH4 oxidation: 2 mol O2 per mol CH4

### Oxygen Saturation

Weiss (1970):
$$\ln(O2_{sat}) = A_1 + \frac{A_2}{T} + A_3 \ln(T) + A_4 T + S[B_1 + B_2 T + B_3 T^2]$$

At 25°C, freshwater: O2_{sat} ≈ 258 µmol/L (8.3 mg/L)

---

## Reaction Network Summary

```
                    ┌──────────────────────────────────────┐
                    │           ATMOSPHERE                 │
                    │    CO2    CH4    N2O    O2          │
                    └────↑↓─────↑↓─────↑↓─────↑↓──────────┘
                         │      │      │      │
    ┌────────────────────┴──────┴──────┴──────┴────────────┐
    │                   WATER COLUMN                        │
    │                                                       │
    │  Phyto ─────→ OC ─────→ CO2 + NH4 + PO4              │
    │    ↑           ↑                                      │
    │    │           │       NH4 ──→ NO2 ──→ NO3           │
    │  Light    Bacteria        ↓ N2O    ↓                 │
    │  + N,P,Si                          ↓ N2O             │
    │                              NO3 ──→ N2              │
    │                                                       │
    │  DIC ←→ CO2 ←→ HCO3⁻ ←→ CO3²⁻                       │
    │                                                       │
    └──────────────────────────┬───────────────────────────┘
                               ↓ settling, diffusion
    ┌──────────────────────────┴───────────────────────────┐
    │                   SEDIMENT                            │
    │                                                       │
    │  OC ──→ CO2 (aerobic)                                │
    │  OC ──→ CH4 (anaerobic)                              │
    │                                                       │
    └──────────────────────────────────────────────────────┘
```

---

## Numerical Solver

### RK4 Adaptive Scheme

C-GEM uses a 4th-order Runge-Kutta solver from C-RIVE:

$$k_1 = f(t, y)$$
$$k_2 = f(t + \frac{\Delta t}{2}, y + \frac{\Delta t}{2} k_1)$$
$$k_3 = f(t + \frac{\Delta t}{2}, y + \frac{\Delta t}{2} k_2)$$
$$k_4 = f(t + \Delta t, y + \Delta t \cdot k_3)$$
$$y_{n+1} = y_n + \frac{\Delta t}{6}(k_1 + 2k_2 + 2k_3 + k_4)$$

### Adaptive Time-Stepping

For stiff reactions, internal subdivision:
$$\Delta t_{bio} = \min\left(\Delta t, \frac{C}{|dC/dt|} \cdot \alpha\right)$$

Where $\alpha$ ≈ 0.5 is a safety factor.

### Mass Balance Checking

Each species tracks:
- Production from each reaction
- Consumption by each reaction
- Net change per timestep

---

## Parameters

### Phytoplankton

| Parameter | Symbol | PHY1 | PHY2 | Unit |
|-----------|--------|------|------|------|
| Max growth | µmax | 2.0 | 1.5 | /day |
| Light slope | α | 0.02 | 0.015 | m²/W |
| Half-sat Si | KSi | 5 | - | µM |
| Half-sat N | KN | 2 | 3 | µM |
| Half-sat P | KP | 0.5 | 0.7 | µM |
| Mortality | kM | 0.05 | 0.04 | /day |
| Settling | ws | 0.5 | 0.3 | m/day |

### Organic Matter (RIVE)

| Pool | k_hydr [/day] | Fraction |
|------|---------------|----------|
| HD1 | 0.75 | 0.15 |
| HD2 | 0.25 | 0.35 |
| HD3 | 0.01 | 0.50 |
| HP1 | 0.50 | 0.10 |
| HP2 | 0.20 | 0.30 |
| HP3 | 0.01 | 0.60 |

### GHG Parameters

| Parameter | Value | Unit | Reference |
|-----------|-------|------|-----------|
| N2O yield (nit) | 0.004 | mol/mol | Cébron 2005 |
| N2O yield (denit) | 0.01 | mol/mol | Garnier 2007 |
| CH4 prod rate | 0.01 | /day | Borges 2011 |
| CH4 ox rate | 0.5 | /day | Borges 2011 |
| CH4 ebul thresh | 50 | µmol/L | DelSontro 2018 |

---

## References

1. Billen, G. et al. (1994). A biogeochemical model of the Scheldt estuary.
2. Garnier, J. et al. (2002). Modelling the nitrogen cycle in the Seine River.
3. Wang, S. et al. (2018). C-RIVE sensitivity analysis. Water Research 144, 341-355.
4. Hasanyar, M. et al. (2022). Unified RIVE v1.0. Biogeosciences.
5. Weiss, R.F. (1974). CO2 in seawater. Marine Chemistry 2, 203-215.
6. **Abril, G. et al. (2009). Turbidity limits gas exchange in a large macrotidal estuary. Estuarine, Coastal and Shelf Science 83, 342-348.** *(k600 for tidal estuaries)*
7. Borges, A.V. et al. (2004). Variability of the gas transfer velocity of CO2 in a macrotidal estuary. Estuaries 27, 593-603.
8. Raymond, P.A. et al. (2012). Scaling gas exchange. Limnol. Oceanogr. Fluids Environ. 2, 1-14. *(for rivers only)*
9. Borges, A.V. & Abril, G. (2011). Carbon dioxide and methane dynamics in estuaries.

---

## Next: [Data Requirements](data_requirements.md)
