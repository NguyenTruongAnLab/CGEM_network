# Boundary Conditions in C-GEM Network

## Overview

This document explains the boundary condition system in C-GEM Network, including common bugs and their solutions. This is **critical reading** for anyone modifying the transport or solver code.

---

## Types of Boundary Nodes

C-GEM supports three types of boundary nodes:

| Type | Code | Description | Typical Use |
|------|------|-------------|-------------|
| `NODE_LEVEL_BC` | 0 | Water level (Dirichlet) | Ocean/tidal boundaries |
| `NODE_DISCHARGE_BC` | 1 | Discharge (Neumann) | River inlets |
| `NODE_JUNCTION` | 2 | Internal junction | Branch confluences |

---

## Ocean Boundary Conditions (NODE_LEVEL_BC)

### Hydrodynamics

At ocean boundaries, water level $H$ is prescribed from tidal forcing files:

```c
node->H = interpolate_forcing(time, node->forcing_time, node->forcing_value, node->forcing_len);
```

### Species Transport

For species concentrations, ocean boundaries require **different treatment for advection vs dispersion**:

#### Advection: Flow-Dependent BC

```c
if (vx > 0) {
    // EBB TIDE (outflow): Neumann BC
    // Let water exit freely - use zero gradient
    conc_ghost = conc_interior;
}
if (vx < 0) {
    // FLOOD TIDE (inflow): Dirichlet BC  
    // Ocean water enters - use forcing value
    conc_ghost = c_ocean;
}
```

#### Dispersion: ALWAYS Dirichlet BC

```c
// Dispersion ALWAYS uses Dirichlet BC regardless of flow direction
// This is critical! Dispersion brings ocean water IN even during ebb tide.
a[1] = 0.0;  // No coupling to ghost cell
bb[1] = 1.0; // Diagonal = 1
c[1] = 0.0;  // No coupling to interior  
d[1] = c_ocean;  // Set to ocean concentration
```

### Physical Basis: Savenije's Steady-State Theory

At steady state, the advection-dispersion equation becomes:

$$Q \times S = D \times A \times \frac{dS}{dx}$$

Where:
- Left side: Advective salt flux OUT (ebb-dominated estuaries)
- Right side: Dispersive salt flux IN (always toward lower concentration)

**Key insight**: Dispersion ALWAYS transports salt INTO the estuary from the high-concentration ocean, regardless of flow direction. Using Neumann BC for dispersion during ebb blocks this mechanism and causes unrealistically low salinity at the mouth.

---

## ⚠️ Critical Bug: Ocean BC Zeroing

### The Problem

In version 1.1.0 and earlier, a critical bug caused all species (except salinity) to have incorrect ocean boundary values. The symptoms were:

- **CH4**: 0 nmol/L at mouth instead of 40 nmol/L
- **N2O**: 0 nmol/L at mouth instead of 8 nmol/L  
- **O2**: Drops to ~200 µM instead of staying at 260 µM
- **TOC**: Incorrect boundary mixing during flood tide

### Root Cause

The junction mixing algorithm `mix_junction_concentrations()` was designed to compute mixed concentrations for species at junctions. However, it had a flaw:

```c
// PROBLEMATIC CODE (before fix):
for (int c = 0; c < node->num_connections; ++c) {
    Branch *b = net->branches[node->connected_branches[c]];
    int dir = node->connection_dir[c];
    
    // This set boundary concentrations for ALL connected branches
    if (dir == 1) {
        set_node_boundary_conc(b, sp, bc_conc, 0);  // Sets conc_down!
    }
}
```

The problem: When no branches flow INTO a junction during ebb tide, the mixed concentration `bc_conc` becomes 0. This zero value was then applied to `conc_down` for all connected branches - including those whose downstream end was an **ocean boundary**, not a junction!

### The Fix (v1.2.0)

Three-part solution:

#### Part 1: Skip Ocean-Connected Boundaries in Junction Mixing

```c
if (dir == 1) {
    // Only set conc_down if downstream is NOT an ocean boundary
    if (b->down_node_type != NODE_LEVEL_BC) {
        set_node_boundary_conc(b, sp, bc_conc, 0);
    }
}
```

#### Part 2: Fallback Ocean Defaults in Transport

```c
// In Transport_Branch_Network():
if (branch->down_node_type == NODE_LEVEL_BC && c_down < 1e-10) {
    // Restore ocean defaults for key species
    switch (sp) {
        case CGEM_SPECIES_SALINITY: c_down = 30.5; break;
        case CGEM_SPECIES_O2: c_down = 260.0; break;
        case CGEM_SPECIES_CH4: c_down = 40.0; break;
        case CGEM_SPECIES_N2O: c_down = 8.0; break;
        // ... other species
    }
}
```

#### Part 3: Strong Boundary Relaxation

```c
// Apply relaxation to ensure cell 1 reflects ocean BC
if (branch->down_node_type == NODE_LEVEL_BC) {
    double alpha = 0.5;  // Strong relaxation
    branch->conc[sp][1] = branch->conc[sp][1] + alpha * (c_down - branch->conc[sp][1]);
}
```

---

## River Boundary Conditions (NODE_DISCHARGE_BC)

### Hydrodynamics

At river inlets, discharge $Q$ is prescribed:

```c
Q_up = interpolate_forcing(time, node->forcing_time, node->forcing_value, node->forcing_len);
```

The water level at the inlet is computed from the momentum equation (Neumann condition).

### Species Transport

River boundaries typically use Dirichlet BC for species:

```c
conc_up[sp] = river_concentration[sp];
```

Species forcing files (`species_river_realistic.csv`) provide time-varying river concentrations.

---

## Junction Conditions (NODE_JUNCTION)

### Hydrodynamics

Junctions enforce:
1. **Mass conservation**: Sum of inflows = sum of outflows
2. **Water level continuity**: All connected branches have same water level at junction

Solved iteratively until residual discharge < tolerance (0.1 m³/s).

### Species Transport

Junction mixing uses flow-weighted averaging:

$$C_{junction} = \frac{\sum_{in} Q_i \times C_i}{\sum_{in} Q_i}$$

Where sums are over branches flowing INTO the junction.

**Important**: Branches flowing OUT of the junction see this mixed concentration as their boundary condition.

---

## Debugging Boundary Issues

### Checklist for Species Boundary Problems

1. **Check forcing files exist and are loaded**:
   ```
   Loaded ALL species forcing for node 6 from: species_ocean_realistic.csv
   ```

2. **Verify species forcing values**:
   - Check that the species column exists in the CSV
   - Verify values are non-zero and reasonable

3. **Check node type assignment**:
   - Ocean nodes should be `NODE_LEVEL_BC`
   - River nodes should be `NODE_DISCHARGE_BC`

4. **Trace `conc_down`/`conc_up` values**:
   - After initialization: should match forcing file values
   - After `apply_species_boundary_forcing()`: should be updated each timestep
   - After `mix_junction_concentrations()`: ocean boundaries should be unchanged

### Common Symptoms and Causes

| Symptom | Likely Cause | Fix |
|---------|--------------|-----|
| Species = 0 at ocean mouth | Junction mixing overwrote ocean BC | Check `down_node_type` before setting BC |
| Species too low at mouth | Weak or missing Dirichlet BC for dispersion | Ensure dispersion always uses Dirichlet |
| Species values oscillating | Timestep too large for dispersion | Reduce dt or increase dispersion |
| No gradient from ocean to river | Transport not enabled for species | Check `CGEM_SPECIES_TRANSPORT_FLAG` |

---

## Best Practices

1. **Always test transport-only first** (`ReactionMode = OFF`) to isolate boundary issues from biogeochemistry

2. **Check all species at boundaries**, not just salinity - the physics is the same but bugs may only manifest for some species

3. **Use validation data at the mouth** to verify ocean BCs are being applied correctly

4. **Monitor `conc_down`/`conc_up` arrays** - these should always match forcing values for boundary nodes

5. **Junction mixing should never modify ocean/river boundaries** - only internal junction concentrations

---

## References

- Savenije, H.H.G. (2005). Salinity and Tides in Alluvial Estuaries. Elsevier.
- Van den Burgh, P. (1972). Ontwikkeling van een methode voor het voorspellen van zoutverdelingen in estuaria, kanalen en zeeen. Rijkswaterstaat.
