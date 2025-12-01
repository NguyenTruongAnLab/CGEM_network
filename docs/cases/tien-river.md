# Tien River Case Study

## Overview

The Tien River case (`Tien_River`) is a simplified test case for validating C-GEM's multi-branch functionality.

### Network Structure

- **5 branches**
- **6 nodes**
- Single bifurcation point

```
          Upstream (Q)
              │
              ▼
         ┌─────────┐
         │  Tien   │
         │  Main   │
         └────┬────┘
              │
         ┌────┴────┐
         │   Bi-   │
         │furcation│
         └────┬────┘
              │
    ┌─────────┼─────────┐
    │         │         │
    ▼         ▼         ▼
┌───────┐ ┌───────┐ ┌───────┐
│Branch │ │Branch │ │Branch │
│  A    │ │   B   │ │   C   │
└───┬───┘ └───┬───┘ └───┬───┘
    │         │         │
    ▼         ▼         ▼
  Ocean     Ocean     Ocean
```

## Running the Case

```powershell
.\bin\Debug\CGEM_Network.exe INPUT/Cases/Tien_River/case_config.txt
```

## Use Case

The Tien River case is ideal for:

1. **Testing new features** - Simpler than full Mekong
2. **Debugging** - Faster runtime
3. **Learning** - Understand model behavior
4. **Validation** - Compare with analytical solutions

## Expected Results

- Tidal propagation upstream
- Salt intrusion during low discharge
- Junction mass balance conservation
