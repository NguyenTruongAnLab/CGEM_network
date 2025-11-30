# C-GEM Network Documentation

## Model Overview

**C-GEM Network** (Carbon-Generic Estuary Model for Networks) is a specialized 1D estuarine biogeochemical model designed for multi-branch tidal river networks. It couples:

- **Saint-Venant hydrodynamics** on a staggered grid
- **TVD advection-dispersion transport** with Savenije's theory
- **C-RIVE biogeochemistry** (Unified RIVE v1.0) for carbon, nutrients, and greenhouse gases

## Why C-GEM?

While powerful 2D/3D models exist (TELEMAC, Delft3D, MIKE), C-GEM fills a critical niche for **rapid, process-based biogeochemical modeling** of complex deltaic networks. See [Introduction](introduction.md) for the full rationale.

## Documentation Structure

| Document | Description |
|----------|-------------|
| [Introduction](introduction.md) | Why C-GEM exists, comparison with other models |
| [Hydrodynamics](hydrodynamics.md) | Saint-Venant equations, staggered grid, tidal propagation |
| [Transport](transport.md) | Advection-dispersion, TVD schemes, Savenije theory |
| [Biogeochemistry](biogeochemistry.md) | C-RIVE module, carbonate chemistry, GHG emissions |
| [Data Requirements](data_requirements.md) | Input preparation, forcing data, calibration |
| [Quick Start (Windows)](QUICKSTART_WINDOWS.md) | Installation and first run |

## Key Features

- ✅ **Multi-branch network topology** (bifurcations, confluences)
- ✅ **Computationally efficient** (1D, runs in seconds-minutes)
- ✅ **Complete carbon cycle** (DIC, TA, pH, pCO2, CO2 flux)
- ✅ **Greenhouse gas emissions** (CO2, CH4, N2O)
- ✅ **Process-based** (mechanistic, not empirical)
- ✅ **Open source** (ANSI C, portable)

## Quick Links

- **Build**: `scripts/build.bat` or `scripts/build-and-run.ps1`
- **Test case**: `INPUT/Cases/Tien_River/`
- **Output**: `OUTPUT/Tien_River/`

## Citation

If you use C-GEM Network, please cite:

## License

[To be determined - suggest MIT or GPL-3.0]
