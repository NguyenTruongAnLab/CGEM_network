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

### Conservative Form

For numerical stability, rewritten as:

$$\frac{\partial C}{\partial t} + U\frac{\partial C}{\partial x} = \frac{1}{A}\frac{\partial}{\partial x}\left(A K \frac{\partial C}{\partial x}\right) + S$$

Where $U = Q/A$ is the cross-sectional average velocity.

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

Where:
- $K_0$ = Dispersion at the mouth [m²/s]
- $L_c$ = Convergence length [m]
- $\beta$ = Van den Burgh coefficient [-]

### Dispersion at the Mouth

$$K_0 = \frac{1400 \cdot H_0 \cdot \sqrt{g H_0}}{f \cdot N_r}$$

Where:
- $H_0$ = Depth at mouth [m]
- $f$ = Darcy-Weisbach friction factor ≈ $8g/C^2$
- $N_r$ = Canter-Cremers estuary number = $Q_r T / P_t$
- $Q_r$ = River discharge [m³/s]
- $T$ = Tidal period [s]
- $P_t$ = Tidal prism [m³]

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

To prevent numerical oscillations near sharp fronts (e.g., salinity intrusion), C-GEM uses **Total Variation Diminishing** schemes.

!!! warning "Grid Convention: Velocity Sign and Upwind Direction"
    In C-GEM's staggered grid:
    
    - **Index 1 = downstream (ocean)**, **Index M = upstream (river)**
    - **Positive velocity (vx > 0)** = flow from upstream toward downstream = **leftward in array**
    - **Upwind cell for vx > 0** = higher index cell (j+2), NOT lower index (j)
    
    This is opposite to many textbook conventions where positive velocity means rightward flow!

**First-order upwind** (diffusive but stable):

$$F_{i+1/2} = \begin{cases} U_{i+1/2} C_{i+1} & \text{if } U > 0 \text{ (flow toward lower indices)} \\ U_{i+1/2} C_i & \text{if } U < 0 \text{ (flow toward higher indices)} \end{cases}$$

!!! note "Flux Sign Convention"
    C-GEM uses `flux = -vx * A * C` so that positive flux represents flow in the positive x direction (toward higher indices). This ensures consistency with the standard finite volume update formula.

**Second-order with flux limiter** (less diffusive):

$$F_{i+1/2} = F^{low}_{i+1/2} + \psi(r) (F^{high}_{i+1/2} - F^{low}_{i+1/2})$$

Where $\psi(r)$ is a limiter function and $r$ is the ratio of consecutive gradients.

### Flux Limiters

| Limiter | Formula | Properties |
|---------|---------|------------|
| Minmod | $\max(0, \min(1, r))$ | Most diffusive TVD |
| Van Leer | $\frac{r + |r|}{1 + |r|}$ | Smooth |
| Superbee | $\max(0, \min(2r, 1), \min(r, 2))$ | Least diffusive TVD |
| MC | $\max(0, \min(2r, (1+r)/2, 2))$ | Monotonized central |

C-GEM default: **Van Leer** (good balance of accuracy and stability)

### Dispersion: Implicit Crank-Nicolson

$$\frac{C^{n+1}_i - C^n_i}{\Delta t} = \frac{\theta}{A_i \Delta x^2}\left[A_{i+1/2}K_{i+1/2}(C^{n+1}_{i+1} - C^{n+1}_i) - A_{i-1/2}K_{i-1/2}(C^{n+1}_i - C^{n+1}_{i-1})\right] + (1-\theta)[\text{same at } n]$$

With $\theta = 0.5$ (Crank-Nicolson) for second-order time accuracy.

This yields a tridiagonal system solved by Thomas algorithm.

---

## Boundary Conditions

### Downstream (Ocean) Boundary

!!! warning "Critical: Separate Treatment for Advection vs Dispersion"
    The advection and dispersion operators require **different** boundary conditions at the ocean. Mixing these up causes either salinity plateaus (blocked dispersion) or unrealistic salt accumulation.

#### Advection Boundary (Flow-Dependent)

**Flood tide** (vx < 0, flow into estuary):
$$C_{boundary} = C_{ocean}$$

**Ebb tide** (vx > 0, flow out of estuary):
$$\frac{\partial C}{\partial x} = 0 \quad \text{(zero gradient - let interior water exit)}$$

#### Dispersion Boundary (Always Dirichlet)

$$C_{boundary} = C_{ocean} \quad \text{(always, regardless of flow direction)}$$

This is the key insight from **Savenije (2005)**: at steady state, the advective flux OUT (river flushing) is balanced by the dispersive flux IN (tidal mixing):

$$Q \cdot S = D \cdot A \cdot \frac{dS}{dx}$$

Dispersion must **always** transport salt up the concentration gradient (from ocean into estuary), even during ebb when advection is exporting salt. Using Neumann BC for dispersion during ebb would block this mechanism and cause unrealistically low salinity at the mouth.

!!! danger "Common Bug: Salinity Plateau or Zero Salt at Mouth"
    - **Plateau at 25-30 psu extending 50+ km**: Usually caused by inverted advection sign convention (upwinding from wrong cell).
    - **Near-zero salinity at mouth**: Usually caused by using Neumann BC for dispersion during ebb, blocking salt entry.

### Upstream (River) Boundary

$$C_{boundary} = C_{river}(t)$$

Time-varying boundary conditions can represent:
- Seasonal nutrient loads
- Episodic pollution events
- Diurnal patterns (e.g., photosynthesis)

### Junction Mixing

At network junctions, concentrations mix according to flow:

$$C_{mixed} = \frac{\sum_i Q_i C_i}{\sum_i Q_i}$$

Branches receiving flow from the junction use $C_{mixed}$ as their boundary condition.

---

## Species-Specific Transport

### Transported Species

| Index | Species | Notes |
|-------|---------|-------|
| 0 | Salinity | Conservative tracer |
| 1-2 | PHY1, PHY2 | Phytoplankton (settling) |
| 3 | DSi | Dissolved silica |
| 4-5 | NO3, NH4 | Nitrogen species |
| 6 | PO4 | Phosphate |
| 7 | O2 | Dissolved oxygen |
| 9 | SPM | Suspended sediment (settling) |
| 10-11 | DIC, TA | Carbonate system |
| 17-26 | RIVE pools | Organic matter, bacteria |
| 27-29 | NO2, N2O, CH4 | GHG species |

### Diagnostic Species (Not Transported)

| Index | Species | Computed From |
|-------|---------|---------------|
| 8 | TOC | Sum of OC pools |
| 12 | pCO2 | Carbonate equilibrium |
| 13 | CO2 | Carbonate equilibrium |
| 14 | pH | DIC and TA |

### Settling

Particulate species settle at velocity $w_s$:

$$S_{settling} = -\frac{w_s}{H} C$$

Typical settling velocities:
- PHY: 0.1-1.0 m/day
- SPM: 1-10 m/day (size-dependent)
- HP (particulate OC): 0.5-2.0 m/day

---

## Numerical Stability

### Peclet Number

For advection-dispersion stability:

$$Pe = \frac{U \Delta x}{K} < 2$$

If $Pe > 2$, numerical oscillations may occur. Solutions:
- Reduce $\Delta x$
- Use TVD scheme
- Increase $K$ (artificial dispersion)

### Courant Number

For advection:

$$Cr = \frac{U \Delta t}{\Delta x} < 1$$

For typical estuarine conditions ($U$ = 1 m/s, $\Delta x$ = 2000 m):

$$\Delta t < 2000 \text{ s}$$

The default $\Delta t$ = 300 s provides a safety factor of ~7.

### Mass Conservation

C-GEM monitors total mass in each branch:

$$M = \sum_i A_i C_i \Delta x$$

Mass should be conserved to < 0.1% over a tidal cycle (excluding reactions).

---

## Salt Intrusion Validation

### Savenije's Analytical Solution

For a converging estuary with exponential width:

$$S(x) = S_0 \exp\left(-\frac{x}{L_s}\right)$$

Where the intrusion length $L_s$ depends on:
- River discharge $Q_r$
- Tidal range
- Geometry ($B_0$, $H_0$, $L_c$)

### Salt Intrusion Length

Empirical relation (Savenije 2005):

$$L_s = \frac{H_0 \cdot C^2}{2g} \cdot \frac{1}{N_r}$$

For the Mekong in dry season:
- $H_0$ ≈ 10 m
- $C$ ≈ 50 m^{1/2}/s
- $N_r$ ≈ 0.1
- **$L_s$ ≈ 60-70 km** (matches observations)

### Calibration Procedure

1. Run with observed river discharge
2. Compare simulated vs. observed salinity at stations
3. Adjust $K_0$ and $\beta$ to match intrusion length
4. Verify tidal salinity variation

---

## Implementation

### Main Transport Function

```c
int Transport_Branch(Branch *branch, double dt);
int Transport_Branch_Network(Branch *branch, double dt, void *network_ptr);
```

**Steps:**
1. Compute dispersion coefficients
2. Apply boundary conditions (flow-dependent)
3. Build tridiagonal matrix for dispersion
4. Solve for each transported species
5. Apply flux limiters for advection

### Dispersion Configuration

```c
void ComputeDispersionCoefficient(Branch *branch, double Q_total);
```

Uses Savenije's theory:
```c
// K = K0 * exp(-beta * x / Lc)
double K0 = branch->D0;
double beta = branch->vdb_coef;
double Lc = branch->lc_convergence;
```

### Species Transport Flags

In `define.h`:
```c
static const int CGEM_SPECIES_TRANSPORT_FLAG[CGEM_NUM_SPECIES] = {
    1,  // SALINITY - transport
    1,  // PHY1 - transport
    ...
    0,  // TOC - diagnostic (computed)
    0,  // pCO2 - diagnostic (computed)
    ...
};
```

---

## Example: Salinity Profile

### Initial Condition

Linear gradient from ocean to river:
```
x=0 (ocean):   S = 30 PSU
x=L (river):   S = 0 PSU
```

### After 1 Tidal Cycle

Tidal mixing creates characteristic profile:
```
          30 |____
             |    \___
Salinity     |        \___
(PSU)     15 |            \___
             |                \___
           0 |____________________\____
             0     20    40    60    80 km
```

### Seasonal Variation

```
Dry season:  Intrusion 60-70 km, S > 4 at Can Tho
Wet season:  Intrusion 10-20 km, S < 1 at Can Tho
```

---

## Troubleshooting Salinity Transport

### Common Bugs and Solutions

| Symptom | Root Cause | Solution |
|---------|------------|----------|
| **Salinity plateau 25-30 psu for 50+ km** | Inverted upwind direction in advection | Use j+2 as upwind when vx > 0 (flow toward lower indices) |
| **Near-zero salinity at mouth (2-8 psu)** | Neumann BC for dispersion during ebb | Always use Dirichlet BC for dispersion at ocean |
| **Salt accumulating unrealistically** | Dirichlet BC for advection during ebb | Use Neumann BC for advection during ebb |
| **No salt intrusion at all** | Dispersion disabled | Set `DISABLE_DISPERSION_FOR_TEST = 0` |
| **Checkerboard oscillations** | Wrong flux convention sign | Use `flux = -vx * A * C` |

### Verification Procedure

After running a simulation, check salinity at the mouth:

```powershell
$csv = Import-Csv "OUTPUT/YourCase/CSV/Branch_salinity.csv"
$last = $csv[-1]
"2km: $($last.'2km')"   # Should be ~32-33 psu
"10km: $($last.'10km')" # Should be 8-15 psu for dry season
```

Expected values for dry season (Q ~ 3000 m³/s):

| Distance | Expected Salinity |
|----------|-------------------|
| 2 km     | 32-33 psu         |
| 6 km     | 15-25 psu         |
| 10 km    | 8-15 psu          |
| 20 km    | < 2 psu           |
| 40+ km   | < 0.5 psu         |

---

## References

1. Savenije, H.H.G. (2005). Salinity and Tides in Alluvial Estuaries. Elsevier.
2. Fischer, H.B. et al. (1979). Mixing in Inland and Coastal Waters. Academic Press.
3. Van Leer, B. (1979). Towards the ultimate conservative difference scheme. J. Comput. Phys.
4. Sweby, P.K. (1984). High resolution schemes using flux limiters. SIAM J. Numer. Anal.
5. Nguyen, A.D. et al. (2008). Salt intrusion in multi-channel estuaries: Mekong Delta. Water Resour. Res.

---

## Next: [Biogeochemistry Module](biogeochemistry.md)

## Next: [Sediment Module](sediment.md)