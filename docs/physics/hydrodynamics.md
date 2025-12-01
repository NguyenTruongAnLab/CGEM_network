# Hydrodynamics Module

## Overview

The C-GEM hydrodynamics module solves the **1D Saint-Venant equations** on a staggered grid, capturing tidal propagation, river discharge, and water level variations in multi-branch estuarine networks.

---

## Governing Equations

### Saint-Venant System

The shallow water equations in 1D, assuming hydrostatic pressure and well-mixed cross-sections:

**Continuity (Mass Conservation):**

$$\frac{\partial A}{\partial t} + \frac{\partial Q}{\partial x} = 0$$

**Momentum:**

$$\frac{\partial Q}{\partial t} + \frac{\partial}{\partial x}\left(\frac{Q^2}{A}\right) + gA\frac{\partial H}{\partial x} + gA S_f = 0$$

Where:
- $A$ = Cross-sectional area [m²]
- $Q$ = Discharge [m³/s]
- $H$ = Water surface elevation [m]
- $g$ = Gravitational acceleration (9.81 m/s²)
- $S_f$ = Friction slope [-]

### Friction Formulation

Using the Chézy-Manning approach:

$$S_f = \frac{|Q|Q}{C^2 A^2 R}$$

Where:
- $C$ = Chézy coefficient [m^{1/2}/s]
- $R$ = Hydraulic radius ≈ $A/B$ for wide channels [m]
- $B$ = Channel width [m]

The Chézy coefficient relates to Manning's $n$:

$$C = \frac{R^{1/6}}{n}$$

Typical values for tropical deltas:
- $n$ = 0.02-0.03 for main channels
- $n$ = 0.03-0.05 for vegetated areas

---

## Numerical Scheme

### Staggered Grid (Preissmann Box)

C-GEM uses a staggered arrangement following Abbott (1979):

```
     U[0]    H[1]    U[2]    H[3]    U[4]    H[5]
      |       o       |       o       |       o
      |       |       |       |       |       |
    --|-------|-------|-------|-------|-------|--
   x=0      dx/2     dx     3dx/2    2dx    5dx/2

   Even indices: Velocity (U), Discharge (Q)
   Odd indices:  Scalars (H, A, B, C)
```

**Advantages:**
- Natural treatment of boundary conditions
- Stable for subcritical flow
- Straightforward extension to networks

### Grid Indexing Convention

| Index | Variable Type | Location |
|-------|--------------|----------|
| 0 | Velocity | Downstream boundary |
| 1 | Scalar | First internal point |
| 2 | Velocity | First internal flux |
| ... | ... | ... |
| M | Scalar | Last internal point |
| M+1 | Velocity | Upstream boundary |

**Direction convention:**
- Index 1 = downstream (ocean side)
- Index M = upstream (river side)
- Positive velocity = upstream flow (flood tide)

### Implicit Time Integration

The momentum equation is discretized implicitly to allow larger time steps:

$$\frac{Q^{n+1} - Q^n}{\Delta t} + \text{advection} + g A \frac{\partial H^{n+1}}{\partial x} + \text{friction}^{n+1} = 0$$

This yields a tridiagonal system solved by the Thomas algorithm.

### CFL Condition

Stability requires:

$$\Delta t < \frac{\Delta x}{\sqrt{gH} + |U|}$$

For typical delta parameters ($H$ = 10 m, $U$ = 1 m/s, $\Delta x$ = 2000 m):

$$\Delta t < \frac{2000}{\sqrt{9.81 \times 10} + 1} \approx 180 \text{ s}$$

The default time step of 300 s is stable for most conditions with implicit friction.

---

## Channel Geometry

### Exponentially Varying Width

Following Savenije's theory, channel width varies exponentially:

$$B(x) = B_0 \exp\left(-\frac{x}{L_c}\right)$$

Where:
- $B_0$ = Width at mouth [m]
- $L_c$ = Convergence length [m] (negative = diverging)

This is the characteristic shape of alluvial estuaries formed by tidal-fluvial interaction.

### Cross-Sectional Area

For a trapezoidal or approximately rectangular section:

$$A = B \cdot H$$

The model computes:
- `totalArea`: Total cross-sectional area
- `freeArea`: Area above reference level
- `refArea`: Reference area at datum

### Storage Width Ratio

For channels with significant floodplain storage:

$$R_s = \frac{B_{storage}}{B_{conveyance}}$$

- $R_s$ = 1.0 for prismatic channels
- $R_s$ = 2-10 for mangrove-fringed channels

This affects tidal damping and phase.

---

## Boundary Conditions

### Downstream (Ocean) Boundary

**Tidal level forcing:**

$$H(t) = H_{MSL} + A_0 \sin\left(\frac{2\pi t}{T}\right)$$

For realistic multi-constituent tides:

$$H(t) = H_{MSL} + \sum_i A_i \cos(\omega_i t - \phi_i)$$

Where the dominant constituents in the Mekong are:
- M2: $T$ = 12.42 hours, $A$ ≈ 0.8-1.2 m
- S2: $T$ = 12.00 hours, $A$ ≈ 0.3-0.5 m
- K1: $T$ = 23.93 hours, $A$ ≈ 0.4-0.6 m
- O1: $T$ = 25.82 hours, $A$ ≈ 0.3-0.4 m

### Upstream (River) Boundary

**Discharge forcing:**

$$Q(t) = Q_{river}(t)$$

Seasonal patterns in the Mekong:
- Dry season (Jan-May): 2,000-5,000 m³/s
- Wet season (Jul-Oct): 20,000-40,000 m³/s
- Transition (Jun, Nov-Dec): 5,000-15,000 m³/s

### Junction Conditions

At network junctions, C-GEM enforces:

1. **Mass conservation:**
$$\sum_i Q_i = 0$$

2. **Common water level:**
$$H_1 = H_2 = ... = H_n$$

Where subscripts denote connected branches.

---

## Network Topology

### Branch Structure

Each branch is defined by:

```c
typedef struct {
    int id;
    char name[64];
    int node_up;      // Upstream node ID
    int node_down;    // Downstream node ID
    double length_m;  // Channel length [m]
    double width_up_m, width_down_m;  // Width at ends [m]
    double depth_m;   // Reference depth [m]
    double chezy;     // Friction coefficient
    double lc_convergence;  // Convergence length [m]
    // ... arrays for H, U, A, B
} Branch;
```

### Node Types

| Type | Description | Boundary Condition |
|------|-------------|-------------------|
| `NODE_JUNCTION` | Internal connection | Mass balance, common H |
| `NODE_LEVEL_BC` | Tidal boundary | Prescribed H(t) |
| `NODE_DISCHARGE_BC` | River boundary | Prescribed Q(t) |

### Example: Tien River Network

```
              [2] MyTho ← (ocean)
             ↗
[river] → [1] Tien_Main → Junction
             ↘
              [3] HamLuong ← (ocean)
             ↘
              [4] CoChien ← (ocean)
```

---

## Tidal Propagation

### Tidal Amplification/Damping

In a converging estuary, tidal amplitude can increase upstream due to funnel effect, or decrease due to friction:

$$\frac{dA}{dx} = A \left(\frac{1}{L_c} - \delta\right)$$

Where $\delta$ is the damping coefficient depending on:
- Friction (Manning's n)
- Geometry (width, depth)
- River discharge

### Phase Lag

The tidal wave propagates at the celerity:

$$c = \sqrt{gH}$$

Phase lag between mouth and upstream location:

$$\Delta\phi = \int_0^x \frac{dx'}{c(x')}$$

---

## Code Implementation

### Main Hydrodynamics Function

```c
int Hyd_Branch(Branch *branch, double H_down, double H_up, double Q_up, double dt);
```

**Algorithm:**
1. Set boundary conditions
2. Compute friction terms (implicit)
3. Build tridiagonal system
4. Solve for H, U at new time
5. Update area, velocity arrays

### Solver Configuration

```c
// In case_config.txt
TimeStepSeconds = 300   // dt [s]
DELXI = 2000           // dx [m]
TidalAmplitude = 1.9   // A0 [m]
RiverDischarge = 2000  // Q_river [m³/s]
```

---

## Validation

### Analytical Solutions

For idealized cases, compare with:
- Tidal wave in prismatic channel
- Steady flow with friction
- Savenije's salt intrusion curves

### Field Data

For the Mekong Delta:
- Water level gauges (Can Tho, My Thuan, Tan Chau)
- Discharge measurements (MRC stations)
- ADCP transects (MRC, SIWRP)

### Performance Metrics

- RMSE of water level < 0.1 m
- Phase error < 30 minutes
- Amplitude error < 10%

---

## References

1. Abbott, M.B. (1979). Computational Hydraulics. Pitman.
2. Savenije, H.H.G. (2005). Salinity and Tides in Alluvial Estuaries. Elsevier.
3. Cunge, J.A., Holly, F.M., Verwey, A. (1980). Practical Aspects of Computational River Hydraulics. Pitman.
4. Nguyen, A.D. et al. (2008). Water Resources Research. Mekong Delta salt intrusion.

---

## Next: [Transport Module](transport.md)
