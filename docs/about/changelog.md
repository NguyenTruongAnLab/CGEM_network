# Changelog

All notable changes to C-GEM Network will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.3.1] - 2025-12-06

### Fixed

- **Critical N2O Unit Mismatch Bug** - N2O evasion was 1000× too weak:
  - `crive_calc_n2o_flux()` returned nmol/L/s but was used in µmol/L/s context
  - Fixed by converting both input (N2O_conc µmol/L → nmol/L) and output (nmol/L/s → µmol/L/s)
  - Benthic N2O flux conversion was dividing by 1e6 instead of 1000
  - **Result**: N2O PBIAS went from +30,000% to -50% (now realistic!)

### Changed

- **`crive_calc_n2o_flux()` in ghg_module.c**:
  - Now accepts N2O_conc in µmol/L (model internal units)
  - Internally converts to nmol/L for flux calculation
  - Returns rate in µmol/L/s (consistent with biogeo.c)

- **Benthic N2O flux in biogeo.c**:
  - Fixed unit conversion: `benthic_N2O_flux / depth / RIVE_SECONDS_PER_DAY / 1000.0`
  - (Was incorrectly dividing by 1e6, making benthic N2O flux 1000× too small)

## [1.3.0] - 2025-12-06

### Added

- **Active Sediment Layer (SOC Pool)** - Dynamic benthic fluxes for scenario analysis:
  - New `sediment_diagenesis.c/h` module implementing single-layer sediment organic carbon pool
  - POC deposition from dead phytoplankton, SPM-bound organic carbon, and settling TOC
  - Temperature-dependent mineralization with Q10 correction
  - Dynamic benthic fluxes proportional to SOC pool (SOD, NH4, PO4, DIC, CH4, N2O)
  - Burial to permanent storage
  - **Key benefit**: When `enable_soc=1`, pollution reduction scenarios correctly show decreasing benthic fluxes over time (vs. fixed fluxes that never change)

- **New parameters in `biogeo_params.txt`**:
  - `enable_soc` - Toggle dynamic SOC (0=fixed fluxes, 1=dynamic)
  - `k_soc_20C` - SOC mineralization rate [1/day]
  - `soc_Q10` - Temperature coefficient
  - `soc_init` - Initial SOC pool [g C/m²]
  - `soc_max` - Maximum SOC capacity [g C/m²]
  - `k_burial` - Permanent burial rate [1/day]
  - `soc_f_anaerobic` - Anaerobic fraction
  - `soc_ch4_yield` - CH4 yield from anaerobic decay
  - `soc_n2o_yield` - N2O yield from sediment N cycling

- **New reaction indices** for SOC-related processes:
  - `CGEM_REACTION_SOC_DEPOSITION` - POC settling to sediment
  - `CGEM_REACTION_SOC_MINERAL` - SOC mineralization
  - `CGEM_REACTION_SOC_BURIAL` - Permanent burial
  - `CGEM_REACTION_SOC_SOD` - Dynamic sediment oxygen demand
  - `CGEM_REACTION_SOC_NH4/PO4/DIC/CH4/N2O` - Dynamic benthic nutrient/GHG fluxes

### Changed

- **Conditional benthic flux calculation** in `biogeo.c`:
  - When `enable_soc=0`: Uses fixed benthic fluxes with spatial scaling (backward compatible)
  - When `enable_soc=1`: Fixed benthic rates set to zero, fluxes calculated by SOC module
  - Prevents double-counting of benthic processes

### Technical Details

The SOC module implements simplified early diagenesis:

```
dSOC/dt = Deposition - k_soc(T) × SOC - k_burial × SOC

Deposition = phy_death × f_settle × depth × 0.012 
           + SPM × ws × poc_ratio × 86400 × 0.001
           + TOC_labile × f_particulate × ws_poc × 86400 × 0.012

Fluxes (scaled by SOC × k_soc):
  - SOD = mineralization × (1 - f_anaerobic) [mmol O2/m²/day]
  - NH4 = SOD × N:C_ratio [mmol N/m²/day]  
  - PO4 = SOD × P:C_ratio [mmol P/m²/day]
  - DIC = mineralization × RQ [mmol C/m²/day]
  - CH4 = mineralization × f_anaerobic × ch4_yield [µmol/m²/day]
  - N2O = NH4_flux × n2o_yield [nmol/m²/day]
```

**Validation Status**: The SOC module is implemented but `enable_soc=0` by default for backward compatibility. Enable for scenario analysis after baseline calibration.

**References**: 
- Chapra (2008) Surface Water-Quality Modeling, Chapter 24
- DiToro (2001) Sediment Flux Modeling
- Soetaert et al. (1996) Est. Coast. Shelf Science

---

## [1.2.0] - 2025-12-03

### Added

- **ReactionMode Configuration** - Toggle biogeochemistry ON/OFF in `case_config.txt`:
  - `ReactionMode = ON` (default): Full biogeochemistry with nutrient cycling, respiration, GHG emissions
  - `ReactionMode = OFF`: Transport-only mode for testing boundary conditions and advection-dispersion
  - Useful for debugging, isolating transport vs reaction effects, and performance optimization

### Fixed

- **Critical: Ocean Boundary Concentration Bug** - All species (not just salinity) now correctly maintain ocean boundary values:
  
  **Root Cause**: Junction mixing algorithm (`mix_junction_concentrations()`) was inadvertently zeroing `conc_down[]` for species when no flow entered junctions during ebb tide. This affected CH4, N2O, O2, TOC, and all other transported species.
  
  **Symptoms Fixed**:
  - CH4/N2O showing 0 at estuary mouth instead of ocean values (40/8 nmol/L)
  - O2 dropping to unrealistic values at downstream boundary
  - TOC not matching ocean boundary during flood tide
  - All species showing incorrect boundary behavior at ocean nodes
  
  **Solution** (three-part fix):
  1. Modified `mix_junction_concentrations()` to skip ocean-connected boundaries (`NODE_LEVEL_BC`)
  2. Added fallback ocean defaults in transport.c for critical species when `c_down` becomes zero
  3. Added strong boundary relaxation (α=0.5) at ocean boundaries to ensure cell 1 reflects forcing

- **Species Boundary Forcing Application** - Ensured `apply_species_boundary_forcing()` correctly sets all species from CSV files every timestep, not just during initialization

### Changed

- Ocean boundary treatment now robust for all species (was previously only tested for salinity)
- Transport solver applies stronger relaxation at open boundaries for scalar stability
- Ghost cell values now properly maintained for all species throughout simulation

### Technical Details

The ocean boundary fix addresses a fundamental issue with the junction-transport interaction:

```
Problem: During ebb tide (net outflow), junction mixing computed nodeC = 0 for GHG species
         (no inflow → denominator = 0 → nodeC = 0), then propagated zero to all connected branches.

Solution: Check branch endpoint type before applying junction-derived BC:
```

```c
// In mix_junction_concentrations():
if (dir == 1) {
    // Only set conc_down if downstream is NOT an ocean boundary
    if (b->down_node_type != NODE_LEVEL_BC) {
        set_node_boundary_conc(b, sp, bc_conc, 0);
    }
}
```

**Validation Results** (Mekong Delta, Reactions OFF):
| Species | Before Fix (1km) | After Fix (1km) | Expected |
|---------|-----------------|-----------------|----------|
| CH4 | 0 nmol/L | 40 nmol/L | ~40 nmol/L |
| N2O | 0 nmol/L | 8 nmol/L | ~8 nmol/L |
| O2 | ~220 µM | 260 µM | 260 µM |
| Salinity | 30.5 PSU | 30.5 PSU | 30.5 PSU |

---

## [1.1.0] - 2025-12-02

### Added

- **Rainfall-Driven Lateral Loads System** - A major new feature for data-sparse regions:
  - `LateralSeasonalFactors` structure for monthly/daily load multipliers
  - Physics-based factor calculation from rainfall data (first-flush wash-off + dilution)
  - 7 climate presets (Mekong, RedRiver, Ganges, Niger, Irrawaddy, SaigonDongNai, Mediterranean)
  - JAXA/Sentinel land use emission factor mapping (Urban, Rice, Aquaculture, Mangrove, Fruit, Forest)
  - Automatic point source generation for cities with treatment efficiency
  - Daily interpolation of monthly factors for smooth seasonal transitions

- **New Python Generator Scripts**:
  - `generate_lateral_loads_v2.py` - Smart rainfall-driven lateral load generator
  - Climate preset system using WorldClim/TRMM rainfall data
  - Automatic conversion from rainfall (mm) to Q and concentration factors
  - Support for custom rainfall arrays (12 monthly values)

- **Enhanced Input File Formats**:
  - `lateral_sources.csv` - Base loads with JAXA-derived emission concentrations
  - `lateral_seasonal_factors.csv` - Monthly Q and species multipliers
  - `lateral_daily_factors.csv` - Daily interpolated factors (365 days)
  - `point_sources.csv` - City sewage with population and treatment level

### Changed

- `Biogeo_Branch()` now accepts network pointer for seasonal factor access
- Lateral loads section applies seasonal factors: `Q_actual = Q_base × Q_factor`
- CSV parser auto-detects old vs new lateral_sources.csv format

### Technical Details

The seasonal factor physics follows:

$$Q_{factor} = \frac{Rain_{month}}{Rain_{dry\_base}}$$

$$C_{factor} = \begin{cases} Q_{factor}^{0.25} & Q_{factor} < 5 \text{ (first flush)} \\ (5^{0.25}) \times (5/Q)^{0.4} & Q_{factor} \geq 5 \text{ (dilution)} \end{cases}$$

---

## [1.0.1] - 2025-12-02

### Fixed

- **Critical: Salinity Transport Bug** - Two fundamental issues in `transport.c`:
  1. **Advection sign convention**: Corrected upwind direction for staggered grid. With positive velocity meaning flow from upstream (M) to downstream (1), the upwind cell must be the higher index (j+2), not j.
  2. **Dispersion boundary condition**: Ocean boundary now always uses Dirichlet BC for dispersion, allowing salt to enter via tidal mixing even during ebb-dominated flow. Previous Neumann BC during ebb blocked dispersive salt flux.

- **Dispersion coefficient initialization**: Removed conflicting caps that limited D0 to unrealistically low values.

### Changed

- Ocean boundary treatment follows Savenije (2005) steady-state theory: dispersion constantly transports salt IN while advection exports during ebb.
- Flux convention changed to `flux = -vx * A * C` for proper mass balance with the grid convention.

---

## [1.0.0] - 2025-12-01

### Added

- **Calibration Module** with NLopt integration
  - BOBYQA, Nelder-Mead, DIRECT, COBYLA, CRS2, SBPLX algorithms
  - Multi-stage calibration workflow (hydro → sediment → biogeo)
  - Seasonal objective functions with time-series RMSE
  - `seasonal_targets.csv` for dry/wet season calibration

- **Greenhouse Gas Module**
  - N₂O production from nitrification and denitrification
  - CH₄ production, oxidation, and ebullition
  - Air-water gas exchange

- **Flocculation Model**
  - Salinity-dependent settling velocity
  - Estuarine Turbidity Maximum (ETM) formation

- **Tropical Physics Enhancements**
  - Salinity toxicity for freshwater phytoplankton
  - Dynamic storage ratio RS(H)
  - Manning's n friction option
  - Residual discharge filter

- **Documentation**
  - Full MkDocs Material theme setup
  - Comprehensive user guide
  - API reference
  - Case studies (Mekong Delta)

### Changed

- Default test case changed from Tien_River to Mekong_Delta_Full
- Species count increased from 17 to 30
- Reaction count increased from 27 to 35+
- Build system updated to include NLopt library

### Fixed

- Junction mass balance conservation
- Boundary condition interpolation edge cases
- Memory leaks in calibration module

## [0.9.0] - 2025-11-15

### Added

- C-RIVE biogeochemistry module
- 2-step nitrification (NH₄ → NO₂ → NO₃)
- Carbonate chemistry (DIC, TA, pH, pCO₂)
- RK4 solver for biogeochemistry

## [0.8.0] - 2025-10-01

### Added

- Sediment erosion/deposition module
- SPM transport with settling

## [0.7.0] - 2025-09-01

### Added

- Multi-branch network topology
- Junction iteration solver
- CSV output for each branch

## [0.6.0] - 2025-08-01

### Added

- TVD advection schemes (Superbee limiter)
- Van den Burgh tidal dispersion

## [0.5.0] - 2025-07-01

### Added

- Saint-Venant hydrodynamics
- Staggered grid implementation
- Implicit time stepping

## [0.1.0] - 2025-06-01

### Added

- Initial project structure
- Basic data structures (Branch, Node, Network)
- Configuration file parsing
