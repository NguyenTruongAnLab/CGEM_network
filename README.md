# C-GEM Network

<div align="center">
  <h3>Carbon-Generic Estuary Model for Multi-Branch Networks</h3>
  <p><em>A specialized 1D biogeochemical model for tidal river deltas</em></p>
  
  [![Documentation](https://img.shields.io/badge/docs-mkdocs-blue)](https://nguytruonganlab.github.io/CGEM_network/)
  [![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
</div>

---

## What is C-GEM?

C-GEM Network is a **1D estuarine biogeochemical model** for multi-branch tidal river networks. It combines:

- **Saint-Venant hydrodynamics** on a staggered grid
- **TVD advection-dispersion transport** with Savenije's tidal dispersion theory
- **C-RIVE biogeochemistry** — carbon, nutrients, oxygen, and greenhouse gases
- **NLopt-powered calibration** with seasonal objective functions

> **Full documentation**: [nguytruonganlab.github.io/CGEM_network](https://nguytruonganlab.github.io/CGEM_network/)

## Key Features

| Feature | Description |
|---------|-------------|
| **Multi-branch topology** | Bifurcations, confluences, distributary networks |
| **Fast computation** | 100-1000× faster than 2D/3D models |
| **Complete biogeochemistry** | 30 species, 35+ reactions, carbonate chemistry |
| **GHG emissions** | CO₂, CH₄, N₂O with process attribution |
| **Built-in calibration** | NLopt optimization with seasonal targets |

## Quick Start

```powershell
# Clone and build
git clone https://github.com/nguytruo/CGEM_network.git
cd CGEM_network

# Build and run Mekong Delta case
.\scripts\build-and-run.ps1 -r Mekong_Delta_Full

# Or run with calibration
.\bin\Debug\CGEM_Network.exe INPUT/Cases/Mekong_Delta_Full/case_config.txt --calibrate
```

## Project Structure

```
CGEM_network/
├── src/                    # Source code
│   ├── physics/           # Hydrodynamics + Transport
│   ├── rive/              # C-RIVE Biogeochemistry
│   ├── io/                # Input/Output
│   └── optimization/      # Calibration (NLopt)
├── INPUT/Cases/           # Case configurations
│   └── Mekong_Delta_Full/ # Main test case
├── docs/                  # MkDocs documentation
└── scripts/               # Build and utility scripts
```

## Documentation

| Section | Content |
|---------|---------|
| [Getting Started](https://nguytruonganlab.github.io/CGEM_network/getting-started/quickstart/) | Installation, quick start |
| [User Guide](https://nguytruonganlab.github.io/CGEM_network/user-guide/overview/) | Data preparation, running, output |
| [Physics](https://nguytruonganlab.github.io/CGEM_network/physics/hydrodynamics/) | Hydrodynamics, transport, biogeochemistry |
| [API Reference](https://nguytruonganlab.github.io/CGEM_network/api/structures/) | Data structures, functions |
| [Case Studies](https://nguytruonganlab.github.io/CGEM_network/cases/mekong-delta/) | Mekong Delta application |
## Citation

```bibtex
@software{cgem_network_2025,
  title = {{C-GEM Network}: Carbon-Generic Estuary Model for Multi-Branch Networks},
  author = {Nguyen, Truong-An and {Nguyen Truong An}},
  year = {2025},
  url = {https://github.com/nguytruo/CGEM_network}
}
```

For C-RIVE biogeochemistry, also cite:
- Wang et al. (2018). *Water Research*, 144, 341-355.
- Wang et al. (2024). *GMD*.

## License

MIT License — see [LICENSE](LICENSE) for details.

## Contributing

See [Contributing Guide](https://nguytruonganlab.github.io/CGEM_network/development/contributing/) for guidelines.