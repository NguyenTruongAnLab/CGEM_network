# C-GEM Network Data Requirements

## Overview

C-GEM Network is designed for **data-sparse tropical estuaries** while supporting advanced applications when more data is available. This guide defines three tiers of data requirements:

| Tier | Purpose | Data Effort | Typical Application |
|------|---------|-------------|---------------------|
| **MINIMAL** | Quick assessment | 1-2 days | Screening studies, rapid assessment |
| **STANDARD** | Seasonal modeling | 1-2 weeks | Research, policy support |
| **ADVANCED** | Full calibration | 1+ months | Publication, regulatory compliance |

---

## Quick Reference: What Data Do You REALLY Need?

### Absolute Minimum (Tier 1)
```
✓ Network topology (branch lengths, widths, depths)
✓ Upstream river discharge (single value or monthly means)
✓ Ocean tidal range (single amplitude)
✓ Ocean salinity (single value, e.g., 30-35 PSU)
```
**That's it!** Everything else has sensible defaults.

### Recommended Standard (Tier 2)
```
+ Upstream species concentrations (O₂, nutrients)
+ Monthly rainfall (for lateral loads)
+ Land use map (for NPS estimation)
```

### Full Advanced (Tier 3)
```
+ Hourly discharge time series
+ Multi-station validation data
+ Custom biogeochemistry parameters
+ Calibration targets
```

---

## Comparison with Other Models

| Data Type | SWAT | WASP | Delft3D-WAQ | CE-QUAL-W2 | **C-GEM** |
|-----------|------|------|-------------|------------|-----------|
| DEM | Required | - | - | - | **Optional** |
| Soil map | Required | - | - | - | **Not needed** |
| Hydrodynamic model | - | Recommended | Required | Built-in | **Built-in** |
| Watershed model for NPS | Built-in | Required* | Required* | Required* | **Built-in** |
| Tidal boundary | - | Manual | Required | Required | **Auto-generated** |
| Meteorology | Required | Optional | Required | Required | **Optional** |

*WASP, Delft3D, CE-QUAL-W2 require external watershed models (SWAT, HSPF, LSPC) to generate nonpoint source loads. C-GEM generates them internally from land use + rainfall.

---

## Tier 1: MINIMAL Data Requirements

### Purpose
Rapid assessment of salinity intrusion and basic water quality for planning or screening.

### Required Files

#### 1. `case_config.txt` (Required)
```ini
CaseName = MyEstuary
Topology = topology.csv
Duration = 30
TimeStep = 300.0
```

#### 2. `topology.csv` (Required)
Minimum 4 columns:
```csv
BranchID,BranchName,NodeUp,NodeDown,Length_m,Width_m,Depth_m
1,Main_Channel,1,2,50000,500,8
```

**How to estimate geometry WITHOUT surveys:**
- **Length**: Google Earth/Maps measurement
- **Width**: Satellite imagery at mouth and upstream
- **Depth**: Navigation charts, literature, or estimate 5-15m for tropical deltas

#### 3. Discharge and Tidal Forcing
Two options:

**Option A: Constant values in `case_config.txt`**
```ini
RiverDischarge = 2000.0    # m³/s (annual mean)
TidalAmplitude = 1.5       # meters
```

**Option B: Simple time series (one file)**
```csv
time_s,Q_m3s
0,2000
86400,2100
```

### What Happens with Minimal Data
- Salinity intrusion: ✓ Works (geometry-driven)
- Tidal dynamics: ✓ Works (sinusoidal tide)
- O₂/CO₂: Uses defaults (may have 30-50% uncertainty)
- Nutrients: Conservative transport only
- GHGs: Rough estimates from defaults

### Default Values Applied

| Parameter | Default | Source |
|-----------|---------|--------|
| Ocean salinity | 32 PSU | Typical tropical ocean |
| River O₂ | 200 µmol/L | ~6.4 mg/L, typical tropical |
| River pCO₂ | 3000 ppm | Elevated, typical tropical river |
| Ocean pCO₂ | 420 ppm | Atmospheric equilibrium |
| Chezy coefficient | 55 m^0.5/s | Sand/silt bed |

---

## Tier 2: STANDARD Data Requirements

### Purpose
Seasonal water quality modeling with realistic lateral inputs.

### Additional Files Needed

#### 4. `boundary_map.csv` (Recommended)
```csv
NodeID,Type,FilePath
1,DISCHARGE,forcing_data/river_discharge.csv
2,LEVEL,forcing_data/ocean_tide.csv
```

#### 5. Species Boundary Conditions (Recommended)

**`species_river.csv`** - Upstream endmember:
```csv
time_s,salinity,o2,toc,dic,at,no3,nh4,po4,pco2
0,0.1,180,150,1400,1300,10,2,1,3500
```

**`species_ocean.csv`** - Ocean endmember:
```csv
time_s,salinity,o2,toc,dic,at,no3,nh4,po4,pco2
0,32,250,100,2000,2200,5,0.5,0.3,420
```

**Where to get endmember values:**
- Literature from similar systems
- Regional databases (e.g., MRC for Mekong, GEMS/Water)
- Single field measurement (dry season snapshot)

#### 6. Monthly Rainfall (For Lateral Loads)
The script `generate_lateral_loads_v2.py` accepts:

**Option A: Climate preset**
```bash
python scripts/generate_lateral_loads_v2.py --climate Mekong
```
Available presets: `Mekong`, `RedRiver`, `Ganges`, `Niger`, `Irrawaddy`, `SaigonDongNai`, `Mediterranean`

**Option B: Custom rainfall**
```bash
python scripts/generate_lateral_loads_v2.py --rainfall 15,10,25,60,200,280,300,320,350,280,140,50
```
Format: 12 monthly values in mm (Jan-Dec)

#### 7. `landuse_map.csv` (For Lateral Loads)
Simple format - can be uniform:
```csv
Branch,Distance_km,Pct_Urban,Pct_Rice,Pct_Aqua,Pct_Mangrove,Pct_Fruit
Main_Channel,0,5,70,10,5,10
Main_Channel,10,5,70,10,5,10
```

**How to create WITHOUT detailed GIS:**
1. Use Google Earth to visually estimate percentages per 10km
2. Apply regional averages from literature
3. Use `scripts/generate_synthetic_landuse.py` for demo

### Generated Files
Running the lateral loads script creates:
- `lateral_sources.csv` - Base loads (dry season)
- `lateral_seasonal_factors.csv` - Monthly multipliers
- `point_sources.csv` - Urban sewage estimates

---

## Tier 3: ADVANCED Data Requirements

### Purpose
Publication-quality modeling with calibration and validation.

### Additional Files

#### 8. `biogeo_params.txt` (Customization)
Biogeochemistry parameters. See [Parameter Reference](../api/parameters.md).

```ini
# Decomposition rates (adjust for your system)
kox = 0.005           # TOC oxidation [/day]
knit = 0.12           # Nitrification [/day]

# Sediment-water exchange
benthic_resp_20C = 50.0   # Benthic O₂ demand [mmol/m²/day]
benthic_NH4_flux = 2.0    # NH₄ flux from sediment [mmol/m²/day]
```

#### 9. `calibration_params.csv` (For Optimization)
```csv
Name,TargetType,TargetID,VarType,Min,Max,Initial,Stage
Chezy_Main,BRANCH,Main_Channel,CHEZY,40,70,55,1
Van_den_Burgh_K,BRANCH,Main_Channel,VDB_K,0.2,0.8,0.5,1
```

#### 10. `calibration_targets.csv` (Observations)
```csv
Branch,Variable,Distance_km,Observed_Value,Weight
Main_Channel,SALINITY,20,15.5,1.0
Main_Channel,O2,20,190,0.5
```

### Forcing Data Resolution

| Data Type | Minimum | Recommended | High-Resolution |
|-----------|---------|-------------|-----------------|
| Discharge | Daily | 3-hourly | Hourly |
| Tidal level | 3-hourly | Hourly | 10-min |
| Species BC | Weekly | Daily | Hourly |
| Meteorology | - | Daily | Hourly |

---

## Data Sources by Region

### Global Free Datasets

| Data | Source | Resolution | Link |
|------|--------|------------|------|
| **Discharge** | GRDC | Monthly | grdc.bafg.de |
| **Tides** | TPXO, FES2014 | Harmonic | aviso.altimetry.fr |
| **Rainfall** | CHIRPS, GPM | Daily-monthly | chc.ucsb.edu/data/chirps |
| **Land Use** | ESA CCI, MODIS | 300m-1km | esa-landcover-cci.org |
| **Climate** | ERA5 | Hourly | cds.climate.copernicus.eu |
| **Ocean WQ** | CMEMS | Monthly | marine.copernicus.eu |

### Regional Data

| Region | Hydrology | Water Quality | Contact |
|--------|-----------|---------------|---------|
| Mekong | MRC Data Portal | MRC WQMN | mrcmekong.org |
| Ganges | India-WRIS | CPCB | india-wris.gov.in |
| Red River | MARD Vietnam | MONRE | dwrm.gov.vn |
| Niger | NBA | - | abn.ne |

---

## Input Data Complexity Analysis

### Current C-GEM Input Structure

| File | Columns/Params | Essential? | Can Default? |
|------|----------------|------------|--------------|
| `topology.csv` | 13 | 7 essential | 6 can default |
| `boundary_map.csv` | 4 | All essential | - |
| `species_river.csv` | 19 | 2-3 essential | Rest can default |
| `species_ocean.csv` | 19 | 2-3 essential | Rest can default |
| `lateral_sources.csv` | 16 | 3 essential | Auto-generated |
| `lateral_seasonal_factors.csv` | 14 | 4 essential | Auto-generated |
| `biogeo_params.txt` | 50+ | 0 essential | All have defaults |

### Truly ESSENTIAL Data (Cross-Model Consensus)

Based on comparison with SWAT, WASP, Delft3D-WAQ, and CE-QUAL-W2:

1. **Network geometry** - Required by all models
2. **River discharge** - Primary transport driver
3. **Tidal boundary** - Mixing energy (estuaries)
4. **Ocean salinity** - Calibration anchor

Everything else can be:
- Estimated from literature
- Generated from land use + rainfall
- Set to regional default values

---

## Simplification Recommendations

### If You Have LIMITED Data

1. **Start with salinity-only mode**
   ```ini
   ReactionMode = OFF
   ```
   This tests hydrodynamics before adding biogeochemistry.

2. **Use climate presets**
   ```bash
   python scripts/generate_lateral_loads_v2.py --climate Mekong
   ```

3. **Use default biogeochemistry**
   Delete `biogeo_params.txt` to use built-in tropical defaults.

4. **Focus calibration on Chezy + K**
   These two parameters control 80% of salinity intrusion.

### If You Have ABUNDANT Data

1. Enable full biogeochemistry with custom parameters
2. Use hourly forcing data
3. Multi-stage calibration (hydro → transport → biogeochem)
4. Validate against multiple stations

---

## Quick Start Checklist

### Minimum Viable Setup (30 minutes)
- [ ] Create case directory: `INPUT/Cases/MyCaseName/`
- [ ] Copy and edit `case_config.txt`
- [ ] Create `topology.csv` from Google Earth measurements
- [ ] Set `RiverDischarge` and `TidalAmplitude` in config
- [ ] Run: `./bin/Debug/CGEM_Network.exe INPUT/Cases/MyCaseName/case_config.txt`

### Standard Setup (1 day)
- [ ] All minimum steps
- [ ] Create `boundary_map.csv` with forcing files
- [ ] Generate `species_river.csv` and `species_ocean.csv`
- [ ] Run: `python scripts/generate_lateral_loads_v2.py --climate YourRegion`
- [ ] Run simulation

### Advanced Setup (1 week+)
- [ ] All standard steps
- [ ] Collect field validation data
- [ ] Create `calibration_params.csv` and `calibration_targets.csv`
- [ ] Run calibration: `--calibrate --stage 1`
- [ ] Iterate on parameters

---

## Comparison: C-GEM vs. "SWAT-lite"

Your concern about C-GEM becoming a "SWAT-lite" is valid. Here's the key difference:

| Aspect | SWAT | C-GEM Network |
|--------|------|---------------|
| **Primary focus** | Watershed processes | Estuarine processes |
| **Spatial resolution** | HRU-based (variable) | 1D along-channel |
| **Hydrodynamics** | Simple routing | Saint-Venant equations |
| **Tidal processes** | Not included | Core feature |
| **Salinity** | Not included | Core feature |
| **Lateral loads** | Process-based (soil, crop) | Emission coefficient |
| **Data requirement** | High (DEM, soil, land use) | Low-Medium |
| **Run time** | Hours-days | Minutes-hours |

**C-GEM's Philosophy:**
> "Emit loads from land use, not compute them from soil physics"

This is the **NEWS (Nutrient Export from Watersheds)** approach, validated globally and appropriate for data-sparse regions where detailed soil and crop data don't exist.

---

## References

1. Garnier, J., et al. (2005). N, P, Si nutrient export from Seine watershed. Biogeochemistry 77, 213-242.
2. Mayorga, E., et al. (2010). Global Nutrient Export from WaterSheds 2 (NEWS 2). Global Biogeochemical Cycles 24.
3. Savenije, H.H.G. (2005). Salinity and Tides in Alluvial Estuaries. Elsevier.
4. Nguyen, A.D., et al. (2008). Salt intrusion in the Mekong Delta. J. Hydrology 356, 356-368.
