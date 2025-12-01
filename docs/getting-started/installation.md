# Installation

## Prerequisites

### Windows

1. **MinGW-w64 compiler** (via MSYS2):
   ```powershell
   # Install MSYS2 from https://www.msys2.org/
   # Then in MSYS2 terminal:
   pacman -S mingw-w64-x86_64-gcc mingw-w64-x86_64-make
   ```

2. **Add to PATH**:
   ```powershell
   # Add C:\msys64\mingw64\bin to your system PATH
   ```

3. **Python 3.8+** (for output processing):
   ```powershell
   # Install from https://www.python.org/
   pip install numpy pandas netCDF4 matplotlib
   ```

### Linux / macOS

1. **GCC compiler**:
   ```bash
   # Ubuntu/Debian
   sudo apt-get install build-essential

   # macOS (Xcode command line tools)
   xcode-select --install
   ```

2. **Python 3.8+**:
   ```bash
   pip install numpy pandas netCDF4 matplotlib
   ```

## Download

### Option 1: Clone Repository

```bash
git clone https://github.com/nguytruonganlab/CGEM_network.git
cd CGEM_network
```

### Option 2: Download Release

Download the latest release from [GitHub Releases](https://github.com/nguytruonganlab/CGEM_network/releases).

## Build

=== "Windows (PowerShell)"

    ```powershell
    cd CGEM_network
    .\scripts\build.bat
    ```

=== "Linux/macOS"

    ```bash
    cd CGEM_network
    make
    ```

The executable will be created at:

- Windows: `bin\Debug\CGEM_Network.exe`
- Linux/macOS: `bin/CGEM_Network`

## Verify Installation

Run the test case:

=== "Windows"

    ```powershell
    .\bin\Debug\CGEM_Network.exe INPUT/Cases/Tien_River/case_config.txt
    ```

=== "Linux/macOS"

    ```bash
    ./bin/CGEM_Network INPUT/Cases/Tien_River/case_config.txt
    ```

Expected output:
```
==============================================
  C-GEM Network Engine
==============================================
  Case        : Tien_River
  Branches    : 5
  Nodes       : 6
  Duration    : 30 days
==============================================

Day    1.0 (  3.3%)          | H=12.50 m, U=0.15 m/s, S=5.2 | 1.2 s elapsed
...
Simulation complete: 8640 steps in 45.2 seconds
```

## Optional: NLopt Library (for Calibration)

The calibration module requires NLopt for optimization:

=== "Windows"

    NLopt is pre-built and included in `external/nlopt/build/`.
    The build script automatically links it.

=== "Linux"

    ```bash
    # Ubuntu/Debian
    sudo apt-get install libnlopt-dev

    # Or build from source
    git clone https://github.com/stevengj/nlopt.git
    cd nlopt
    mkdir build && cd build
    cmake ..
    make
    sudo make install
    ```

=== "macOS"

    ```bash
    brew install nlopt
    ```

## Directory Structure After Installation

```
CGEM_network/
├── bin/
│   └── Debug/
│       ├── CGEM_Network.exe    # Main executable
│       └── libnlopt.dll        # NLopt library (Windows)
├── INPUT/
│   └── Cases/
│       ├── Mekong_Delta_Full/  # Main test case
│       ├── Tien_River/         # Simple test case
│       └── ...
├── OUTPUT/
│   └── (simulation outputs)
├── scripts/
│   ├── build.bat               # Windows build
│   ├── build-and-run.ps1       # Build + run helper
│   └── bin_to_nc.py            # Output conversion
└── src/
    └── (source code)
```

## Troubleshooting

### "gcc not found"

Ensure MinGW bin directory is in PATH:
```powershell
$env:PATH += ";C:\msys64\mingw64\bin"
```

### "libnlopt.dll not found"

Copy the DLL to the bin directory:
```powershell
copy external\nlopt\build\libnlopt.dll bin\Debug\
```

### Build errors with math functions

Ensure `-lm` is included in the linker flags (check `scripts/build.bat`).

## Next Steps

- [Quick Start](quickstart.md) - Run your first simulation
- [Project Structure](structure.md) - Understand the codebase
