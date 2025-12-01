# Transport Module

## Overview

The C-GEM transport module solves the **advection-dispersion equation** for scalar concentrations, including salinity, nutrients, organic matter, and dissolved gases. It employs TVD (Total Variation Diminishing) schemes for numerical stability and Savenije's theory for tidal dispersion.

---

## Governing Equation

### 1D Advection-Dispersion

$$\frac{\partial (AC)}{\partial t} + \frac{\partial (QC)}{\partial x} = \frac{\partial}{\partial x}\left(A K \frac{\partial C}{\partial x}\right) + A \cdot S$$

Where:
- $C$ = Concentration [unit/m³]
- $A$ = Cross-sectional area [m²]
- $Q$ = Discharge [m³/s]
- $K$ = Longitudinal dispersion coefficient [m²/s]
- $S$ = Source/sink term [unit/m³/s]

---

## Tidal Dispersion: Savenije's Theory

### Physical Mechanism

In tidal estuaries, longitudinal dispersion arises from:

1. **Tidal pumping**: Correlation between tidal velocity and concentration
2. **Trapping**: Side embayments fill on flood, release on ebb
3. **Residual circulation**: Gravitational circulation in stratified areas
4. **Shear dispersion**: Vertical/lateral velocity gradients

### Van den Burgh's Coefficient

Savenije (2005) showed that dispersion in alluvial estuaries follows:

$$\frac{dK}{dx} = -K \cdot \beta$$

Where $\beta$ is the **Van den Burgh coefficient** (dimensionless, typically 0.2-0.8).

This leads to:

$$K(x) = K_0 \exp(-\beta x / L_c)$$

### Typical Values for Mekong Delta

| Parameter | Dry Season | Wet Season |
|-----------|------------|------------|
| $K_0$ | 500-1000 m²/s | 200-500 m²/s |
| $\beta$ | 0.3-0.5 | 0.4-0.6 |
| Salt intrusion | 50-70 km | 10-30 km |

---

## Numerical Scheme

### Operator Splitting

The transport equation is split into:

1. **Advection step**: $\frac{\partial C}{\partial t} + U\frac{\partial C}{\partial x} = 0$
2. **Dispersion step**: $\frac{\partial C}{\partial t} = \frac{1}{A}\frac{\partial}{\partial x}\left(AK\frac{\partial C}{\partial x}\right)$
3. **Reaction step**: $\frac{\partial C}{\partial t} = S$

### TVD Advection Schemes

To prevent numerical oscillations near sharp fronts, C-GEM uses **Total Variation Diminishing** schemes with flux limiters.

### Flux Limiters

| Limiter | Formula | Properties |
|---------|---------|------------|
| Minmod | $\max(0, \min(1, r))$ | Most diffusive TVD |
| Van Leer | $\frac{r + \|r\|}{1 + \|r\|}$ | Smooth |
| Superbee | $\max(0, \min(2r, 1), \min(r, 2))$ | Least diffusive |

C-GEM default: **Superbee**

### Dispersion: Implicit Crank-Nicolson

$$\frac{C^{n+1}_i - C^n_i}{\Delta t} = \frac{\theta}{A_i \Delta x^2}\left[A_{i+1/2}K_{i+1/2}(C^{n+1}_{i+1} - C^{n+1}_i) - A_{i-1/2}K_{i-1/2}(C^{n+1}_i - C^{n+1}_{i-1})\right]$$

---

## Boundary Conditions

### Downstream (Ocean)

**Flood tide** (flow into estuary):
$$C_{boundary} = C_{ocean}$$

**Ebb tide** (flow out of estuary):
$$\frac{\partial C}{\partial x} = 0$$

### Upstream (River)

$$C_{boundary} = C_{river}(t)$$

### Junction Mixing

$$C_{mixed} = \frac{\sum_i Q_i C_i}{\sum_i Q_i}$$

---

## Species Transport

### Transported Species

| Index | Species | Notes |
|-------|---------|-------|
| 0 | Salinity | Conservative tracer |
| 1-2 | PHY1, PHY2 | Phytoplankton |
| 3 | DSi | Dissolved silica |
| 4-5 | NO₃, NH₄ | Nitrogen |
| 6 | PO₄ | Phosphate |
| 7 | O₂ | Dissolved oxygen |
| 9 | SPM | Suspended sediment |
| 10-11 | DIC, TA | Carbonate system |
| 27-29 | NO₂, N₂O, CH₄ | GHG species |

### Diagnostic Species (Not Transported)

| Index | Species | Computed From |
|-------|---------|---------------|
| 8 | TOC | Sum of OC pools |
| 12 | pCO₂ | Carbonate equilibrium |
| 13 | CO₂ | Carbonate equilibrium |
| 14 | pH | DIC and TA |

---

## Numerical Stability

### Peclet Number

$$Pe = \frac{U \Delta x}{K} < 2$$

### Courant Number

$$Cr = \frac{U \Delta t}{\Delta x} < 1$$

For typical conditions ($U$ = 1 m/s, $\Delta x$ = 2000 m):
$$\Delta t < 2000 \text{ s}$$

---

## Salt Intrusion Validation

### Savenije's Analytical Solution

$$S(x) = S_0 \exp\left(-\frac{x}{L_s}\right)$$

### Mekong Observations

| Season | Intrusion (4 PSU) |
|--------|-------------------|
| Dry | 50-70 km |
| Wet | 10-20 km |

---

## References

1. Savenije, H.H.G. (2005). *Salinity and Tides in Alluvial Estuaries*. Elsevier.
2. Fischer, H.B. et al. (1979). *Mixing in Inland and Coastal Waters*. Academic Press.
3. Van Leer, B. (1979). Towards the ultimate conservative difference scheme. *J. Comput. Phys.*
4. Nguyen, A.D. et al. (2008). Salt intrusion in multi-channel estuaries. *Water Resour. Res.*
