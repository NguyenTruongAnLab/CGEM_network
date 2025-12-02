# C-GEM Network

<div align="center">
  <h2>Carbon-Generic Estuary Model for Multi-Branch Networks</h2>
  <p><em>A specialized 1D biogeochemical model for tropical river deltas</em></p>
</div>

---

<div class="grid cards" markdown>

-   :material-book-open-variant:{ .lg .middle } __Introduction__

    ---

    Why C-GEM? Comparison with other models and unique features

    [:octicons-arrow-right-24: Learn more](getting-started/introduction.md)

-   :material-water:{ .lg .middle } __Hydrodynamics__

    ---

    Saint-Venant equations on staggered grid with Savenije's estuarine theory

    [:octicons-arrow-right-24: Learn more](physics/hydrodynamics.md)

-   :material-waves:{ .lg .middle } __Transport__

    ---

    TVD advection-dispersion with tidal dispersion and salinity intrusion

    [:octicons-arrow-right-24: Learn more](physics/transport.md)

-   :material-molecule:{ .lg .middle } __Biogeochemistry__

    ---

    C-RIVE module: carbon, nutrients, oxygen, and greenhouse gases (30 species)

    [:octicons-arrow-right-24: Learn more](physics/biogeochemistry.md)

-   :material-tune:{ .lg .middle } __Calibration__

    ---

    NLopt-powered optimization with seasonal data support

    [:octicons-arrow-right-24: Learn more](user-guide/calibration.md)

-   :material-map-marker-path:{ .lg .middle } __Lateral Loads (NEW)__

    ---

    **Rainfall-driven** land-use coupling with global data support (JAXA, WorldClim)

    [:octicons-arrow-right-24: Learn more](user-guide/data-preparation.md#6-lateral-sources-system-new---rainfall-driven)

</div>

## Why C-GEM?

C-GEM Network fills a critical niche for **rapid, process-based biogeochemical modeling** of complex deltaic networks:

| Feature | C-GEM | Delft3D/MIKE | SWAT | HEC-RAS |
|---------|-------|--------------|------|---------|
| **Setup time** | Hours | Days-Weeks | Days | Hours |
| **Run time (30 days)** | Minutes | Hours-Days | Hours | Minutes |
| **Network topology** | ✅ Native | ⚠️ Complex meshing | ❌ | ⚠️ |
| **Tidal dispersion** | ✅ Savenije theory | ⚠️ Numerical | ❌ | ❌ |
| **Full biogeochemistry** | ✅ C-RIVE (30 species) | ⚠️ Limited WQ | ⚠️ | ❌ |
| **GHG emissions** | ✅ CO₂, CH₄, N₂O | ❌ | ❌ | ❌ |
| **Calibration** | ✅ NLopt built-in | ⚠️ External | ✅ | ❌ |
| **Land-use coupling** | ✅ Rainfall-driven | ⚠️ Manual | ✅ | ❌ |
| **Open source** | ✅ MIT License | ❌ Commercial | ✅ | ✅ |

[Learn more](getting-started/introduction.md) about what makes C-GEM unique.

## Quick Start

=== "Windows (PowerShell)"

    ```powershell
    # Clone and build
    git clone https://github.com/nguytruonganlab/CGEM_network.git
    cd CGEM_network
    
    # Generate lateral loads (optional - uses built-in climate presets)
    python scripts/generate_lateral_loads_v2.py --climate Mekong
    
    # Build and run
    .\scripts\build-and-run.ps1 -r Mekong_Delta_Full
    ```

=== "Linux/macOS"

    ```bash
    # Clone and build
    git clone https://github.com/nguytruonganlab/CGEM_network.git
    cd CGEM_network
    
    # Generate lateral loads (optional)
    python scripts/generate_lateral_loads_v2.py --climate Mekong
    
    # Build and run
    make
    ./bin/CGEM_Network INPUT/Cases/Mekong_Delta_Full/case_config.txt
    ```

[Quick start guide](getting-started/quickstart.md) • [Installation](getting-started/installation.md) • [Data preparation](user-guide/data-preparation.md)

## Key Features

-  **Multi-branch network topology** — bifurcations, confluences, distributaries
-  **Computationally efficient** — 1D, runs in seconds to minutes
-  **Complete carbon cycle** — DIC, TA, pH, pCO₂, CO₂ air-water flux
-  **Greenhouse gas emissions** — CO₂, CH₄, N₂O with process attribution
-  **C-RIVE biogeochemistry** — 6-pool organic matter, 2-step nitrification
-  **Automatic calibration** — NLopt integration with seasonal objectives
-  **Open source** — ANSI C, portable, transparent

## Architecture Overview

```mermaid
flowchart TB
    subgraph cgem[C-GEM Network]
        hydro[Hydrodynamics<br/>Saint-Venant<br/>Staggered grid]
        transport[Transport<br/>Advective-dispersive<br/>TVD schemes]
        biogeo[Biogeochemistry<br/>C-RIVE + GHG module]
        solver[Network Solver<br/>Junction mass balance<br/>Multi-branch iteration<br/>Implicit time stepping<br/>Boundary conditions]
        hydro --> transport --> biogeo --> solver
        solver --> calibration[Calibration<br/>NLopt + seasonal objectives]
        solver --> ioManager[I/O Manager<br/>CSV & NetCDF<br/>Config]
        solver --> output[Output<br/>Binary/CSV<br/>Time series]
    end
```

[Learn more](getting-started/structure.md) - Detailed architecture and source code structure

## Case Studies

<div class="grid cards" markdown>

-   :material-map:{ .lg .middle } __Mekong Delta__

    ---

    4-branch network, seasonal calibration, full biogeochemistry

    [:octicons-arrow-right-24: View case study](cases/mekong-delta.md)

-   :material-map:{ .lg .middle } __Tien River__

    ---

    Single distributary test case for validation

    [:octicons-arrow-right-24: View case study](cases/tien-river.md)

</div>

[Learn more](user-guide/overview.md) - Full user guide and case study of Mekong Delta

## Citation

If you use C-GEM Network, please cite:

```bibtex
@software{cgem_network_2025,
  title = {C-GEM Network: Carbon-Generic Estuary Model for Multi-Branch Networks},
  author = {Nguyen, Truong and Nguyen Truong An},
  year = {2025},
  url = {https://github.com/nguytruonganlab/CGEM_network}
}
```

## License

C-GEM Network is released under the [MIT License](about/license.md).
