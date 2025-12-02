# Changelog

All notable changes to C-GEM Network will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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
