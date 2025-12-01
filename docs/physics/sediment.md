# Sediment Transport

## Overview

The sediment module simulates suspended particulate matter (SPM) dynamics including:

- Erosion and deposition based on shear stress
- Salinity-dependent flocculation (ETM formation)
- Light attenuation effects on phytoplankton

## Governing Equation

SPM transport with settling:

$$\frac{\partial (A \cdot SPM)}{\partial t} + \frac{\partial (Q \cdot SPM)}{\partial x} = \frac{\partial}{\partial x}\left(A K \frac{\partial SPM}{\partial x}\right) + A(E - D)$$

Where:
- $E$ = Erosion rate [mg/L/s]
- $D$ = Deposition rate [mg/L/s]

## Erosion

### Partheniades Formula

$$E = M_e \left(\frac{\tau_b}{\tau_{crit,e}} - 1\right)^n \quad \text{when } \tau_b > \tau_{crit,e}$$

Where:
- $M_e$ = Erosion coefficient [mg/m²/s]
- $\tau_b$ = Bed shear stress [Pa]
- $\tau_{crit,e}$ = Critical erosion stress [Pa]
- $n$ = Exponent (typically 1)

### Bed Shear Stress

From hydrodynamics:

$$\tau_b = \rho \frac{g U |U|}{C^2}$$

## Deposition

### Krone Formula

$$D = w_s \cdot SPM \cdot \left(1 - \frac{\tau_b}{\tau_{crit,d}}\right) \quad \text{when } \tau_b < \tau_{crit,d}$$

Where:
- $w_s$ = Settling velocity [m/s]
- $\tau_{crit,d}$ = Critical deposition stress [Pa]

## Flocculation Model

### Salinity-Dependent Settling

In estuaries, salinity enhances flocculation:

$$w_s(S) = w_{s,0} \cdot f_{floc}(S)$$

$$f_{floc}(S) = 1 + (f_{max} - 1) \cdot \tanh\left(\frac{S}{S_{scale}}\right)$$

Where:
- $w_{s,0}$ = Base settling velocity (freshwater)
- $f_{max}$ = Maximum flocculation factor
- $S_{scale}$ = Salinity scale for flocculation

### Typical Values

| Parameter | Mekong | Loire | Unit |
|-----------|--------|-------|------|
| $w_{s,0}$ | 0.5 | 0.3 | mm/s |
| $f_{max}$ | 8 | 5 | - |
| $S_{scale}$ | 2 | 3 | PSU |

## Estuarine Turbidity Maximum (ETM)

The flocculation model creates the ETM:

```
SPM (mg/L)
    200 ┤                  ╭────╮
        │                 ╱      ╲
    150 ┤               ╱         ╲
        │              ╱           ╲
    100 ┤            ╱              ╲
        │          ╱                 ╲
     50 ┤       ╱                     ╲────────
        │    ╱                             
      0 ┼───╯──────────────────────────────────
        0    10    20    30    40    50    60
                 Distance from mouth (km)
                    ↑
                   ETM
```

## Parameters

| Parameter | Symbol | Default | Range | Unit |
|-----------|--------|---------|-------|------|
| Base settling velocity | $w_{s,0}$ | 0.0005 | 0.0001-0.003 | m/s |
| Critical erosion stress | $\tau_{crit,e}$ | 0.3 | 0.1-0.8 | Pa |
| Critical deposition stress | $\tau_{crit,d}$ | 0.08 | 0.03-0.2 | Pa |
| Erosion coefficient | $M_e$ | 0.0001 | 0.00001-0.001 | kg/m²/s |
| Flocculation max | $f_{max}$ | 8 | 3-15 | - |
| Salinity scale | $S_{scale}$ | 2 | 1-5 | PSU |

## Implementation

### Key Functions

```c
// src/rive/sediment.c

int Sediment_Branch(Branch *b, double dt) {
    // For each grid cell
    for (int i = 1; i <= b->M; i++) {
        double tau = calc_bed_shear(b, i);
        double sal = b->conc[SPECIES_SALINITY][i];
        
        // Flocculation-enhanced settling
        double ws = calc_settling_velocity(sal, b->biogeo);
        
        // Erosion
        double E = 0.0;
        if (tau > tau_crit_ero) {
            E = M_ero * (tau / tau_crit_ero - 1.0);
        }
        
        // Deposition
        double D = 0.0;
        if (tau < tau_crit_dep) {
            D = ws * b->conc[SPECIES_SPM][i] * (1.0 - tau / tau_crit_dep);
        }
        
        // Update SPM
        b->conc[SPECIES_SPM][i] += (E - D) * dt / b->depth[i];
    }
}
```

## Light Attenuation

SPM affects phytoplankton through light attenuation:

$$I(z) = I_0 \exp(-k_d z)$$

$$k_d = k_{bg} + k_{SPM} \cdot SPM$$

Where:
- $k_{bg}$ = Background attenuation [m⁻¹]
- $k_{SPM}$ = SPM-specific attenuation [m²/g]

Depth-averaged light:

$$\bar{I} = \frac{I_0}{k_d H}\left(1 - e^{-k_d H}\right)$$

## References

1. Partheniades, E. (1965). Erosion and deposition of cohesive soils. *J. Hydraul. Div.*, 91(1), 105-139.

2. Krone, R.B. (1962). *Flume studies of the transport of sediment in estuarial shoaling processes*. University of California.

3. Winterwerp, J.C. (2002). On the flocculation and settling velocity of estuarine mud. *Continental Shelf Research*, 22(9), 1339-1360.
