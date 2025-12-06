# Sediment Transport

## Overview

The sediment module simulates suspended particulate matter (SPM) dynamics including:

- Erosion and deposition based on shear stress
- Salinity-dependent flocculation (ETM formation)
- Light attenuation effects on phytoplankton
- **Active Sediment Layer (SOC)** - Dynamic benthic fluxes from sediment organic carbon (v1.3.0)

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

---

## Active Sediment Layer (SOC Pool) - v1.3.0

### Why SOC Matters for Scenario Analysis

**Without SOC** (`enable_soc=0`): Benthic fluxes are fixed parameters. If you simulate a 50% reduction in upstream pollution, the sediment still releases the same fluxes forever. **This is WRONG for scenario analysis.**

**With SOC** (`enable_soc=1`): Benthic fluxes are proportional to the accumulated sediment organic carbon pool. Pollution reduction → less POC deposition → SOC pool depletes → benthic fluxes decrease over time. **This captures the "legacy load" effect correctly.**

### SOC Pool Dynamics

The SOC pool [g C/m²] evolves according to:

$$\frac{dSOC}{dt} = F_{dep} - R_{min} - R_{burial}$$

Where:
- $F_{dep}$ = POC deposition flux [g C/m²/day]
- $R_{min}$ = Mineralization rate [g C/m²/day]
- $R_{burial}$ = Permanent burial rate [g C/m²/day]

### POC Deposition Sources

Total deposition is the sum of three sources:

1. **Dead phytoplankton settling**:
$$F_{phy} = \text{phy\_death} \times f_{settle} \times H \times 0.012$$

2. **SPM-bound POC settling**:
$$F_{SPM} = SPM \times w_s \times \text{poc\_ratio} \times 86400 \times 0.001$$

3. **Labile TOC settling** (particulate fraction):
$$F_{TOC} = TOC \times f_{labile} \times f_{part} \times 0.012 \times w_{s,poc} \times 86400$$

### Mineralization

Temperature-dependent first-order decay:

$$R_{min} = k_{soc}(T) \times SOC$$

$$k_{soc}(T) = k_{soc,20} \times Q_{10}^{(T-20)/10}$$

### Benthic Fluxes from SOC

Mineralization produces benthic fluxes proportional to SOC pool:

| Flux | Formula | Unit |
|------|---------|------|
| **SOD** | $R_{min} \times (1 - f_{anaerobic}) / 12$ | mmol O₂/m²/day |
| **NH₄** | $SOD \times N:C_{Redfield}$ | mmol N/m²/day |
| **PO₄** | $SOD \times P:C_{Redfield}$ | mmol P/m²/day |
| **DIC** | $R_{min} \times RQ / 12$ | mmol C/m²/day |
| **CH₄** | $R_{min} \times f_{anaerobic} \times y_{CH4} / 12 \times 1000$ | µmol/m²/day |
| **N₂O** | $F_{NH4} \times y_{N2O} \times 10^6$ | nmol/m²/day |

### SOC Parameters

Configure in `biogeo_params.txt`:

| Parameter | Symbol | Default | Range | Unit | Description |
|-----------|--------|---------|-------|------|-------------|
| `enable_soc` | - | 0 | 0-1 | - | Enable dynamic SOC |
| `k_soc_20C` | $k_{soc,20}$ | 0.003 | 0.001-0.01 | 1/day | Decay rate at 20°C |
| `soc_Q10` | $Q_{10}$ | 2.0 | 1.5-3.0 | - | Temperature coefficient |
| `soc_init` | $SOC_0$ | 500 | 100-2000 | g C/m² | Initial pool |
| `soc_max` | $SOC_{max}$ | 5000 | 1000-10000 | g C/m² | Maximum capacity |
| `k_burial` | $k_b$ | 0.0001 | 0-0.001 | 1/day | Burial rate |
| `soc_f_anaerobic` | $f_{an}$ | 0.3 | 0.1-0.5 | - | Anaerobic fraction |
| `soc_ch4_yield` | $y_{CH4}$ | 0.5 | 0.3-0.6 | mol/mol | CH₄ yield |
| `soc_n2o_yield` | $y_{N2O}$ | 0.02 | 0.01-0.05 | mol/mol | N₂O yield |

### When to Use SOC

| Scenario | Recommendation |
|----------|----------------|
| **Baseline calibration** | `enable_soc=0` (fixed fluxes, fewer parameters) |
| **Pollution reduction scenarios** | `enable_soc=1` (captures legacy effects) |
| **Climate change scenarios** | `enable_soc=1` (temperature affects decay) |
| **Long-term projections** | `enable_soc=1` (sediment pool evolves) |

### References

- Chapra, S.C. (2008). *Surface Water-Quality Modeling*, Chapter 24. McGraw-Hill.
- DiToro, D.M. (2001). *Sediment Flux Modeling*. Wiley.
- Soetaert, K., et al. (1996). *Est. Coast. Shelf Sci.*, 43, 371-403.
- Middelburg, J.J. (1989). *Geochim. Cosmochim. Acta*, 53, 1577-1581.

---

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
