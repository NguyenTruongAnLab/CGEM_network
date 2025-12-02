#ifndef CGEM_DEFINE_H
#define CGEM_DEFINE_H

/* ===========================================================================
 * CGEM-RIVE: Enhanced biogeochemistry based on C-RIVE (Unified RIVE v1.0)
 * 
 * Features:
 * - Full C-RIVE carbonate chemistry (DIC, TA, pH, pCO2)
 * - Complete GHG module (CO2, CH4, N2O)
 * - 2-step nitrification (NH4 → NO2 → NO3)
 * - RK4 adaptive solver from C-RIVE
 * - Multi-pool organic matter and RIVE bacteria
 * 
 * Reference: Wang et al. (2018), Hasanyar et al. (2022), Billen et al. (1994)
 * ===========================================================================*/

#define CGEM_MAX_SPECIES 35  /* Increased for GHG species (NO2, N2O, CH4) */
#define CGEM_MAX_BRANCH_NAME 64
#define CGEM_MAX_PATH 4096

#define CGEM_DEFAULT_BRANCH_CELLS 80
#define CGEM_MAX_BRANCH_CELLS 1000
#define CGEM_MAX_BRANCHES 20
#define CGEM_MAX_TS_STEPS 10000

/* Physical constants (match legacy single-branch solver) */
#define CGEM_GRAVITY 9.81
#define CGEM_PI 3.14159265358979323846
#define CGEM_RHO_WATER 1000.0
#define CGEM_T_KELVIN 273.15    /* Kelvin offset */

/* Simulation defaults (overridden by global/case parameters) */
#define CGEM_DEFAULT_DT_SECONDS 300.0
#define CGEM_DEFAULT_DX_METERS  2000.0

/* Numerical thresholds */
#define CGEM_MIN_DEPTH 0.1
#define CGEM_MIN_WIDTH 1.0
#define CGEM_TOL 1e-5
#define CGEM_MAX_ITER 1000
#define CGEM_EPS 1e-10         /* Small number for numerical stability */

/* RK4 solver constants (from C-RIVE) */
#define CGEM_RK4_STAGES 4
#define CGEM_RK4_MAX_DIV_DT 100.0  /* Maximum time step subdivision */

/* Storage width ratio (RS in Fortran) */
#define CGEM_RS 1.0

/* Tidal parameters */
#define CGEM_M2_FREQ 0.080536912751677847  /* M2 tidal frequency [cycle/hr] */

/* ===========================================================================
 * CONCENTRATION-FLOW MODES FOR LATERAL LOADS
 * 
 * Different land use types have different concentration-flow relationships.
 * These modes are applied per-segment based on dominant land use.
 * 
 * Reference: 
 *   - Dilution: Haygarth & Jarvis (1999), Journal of Env. Quality
 *   - First-flush: Lee et al. (2002), Water Research
 *   - Tidal export: Bouillon et al. (2008), Global Biogeochem. Cycles
 * ===========================================================================*/
typedef enum {
    CONC_FLOW_DILUTION = 0,   /* Agricultural: conc = base / sqrt(flow_mult) */
    CONC_FLOW_CONSTANT = 1,   /* Wastewater: conc = base (independent of flow) */
    CONC_FLOW_INVERSE  = 2,   /* Strong dilution: conc = base / flow_mult */
    CONC_FLOW_TIDAL    = 3    /* Mangrove: conc = base * sqrt(flow_mult) */
} ConcFlowMode;

/* Hydrodynamic variable indices */
#define CGEM_HYDRO_DEPTH       0   /* Water depth H [m] */
#define CGEM_HYDRO_VELOCITY    1   /* Current velocity U [m/s] */
#define CGEM_HYDRO_WATERLEVEL  2   /* Water surface elevation [m] */
#define CGEM_HYDRO_AREA        3   /* Cross-sectional area A [m²] */
#define CGEM_HYDRO_WIDTH       4   /* Channel width B [m] */
#define CGEM_HYDRO_DISPERSION  5   /* Dispersion coefficient K [m²/s] */

#define CGEM_NUM_HYDRO 6

/* ===========================================================================
 * SPECIES INDICES
 * =========================================================================== */

/* Original CGEM species (0-16) */
#define CGEM_SPECIES_SALINITY 0
#define CGEM_SPECIES_PHY1     1
#define CGEM_SPECIES_PHY2     2
#define CGEM_SPECIES_DSI      3
#define CGEM_SPECIES_NO3      4
#define CGEM_SPECIES_NH4      5
#define CGEM_SPECIES_PO4      6
#define CGEM_SPECIES_O2       7
#define CGEM_SPECIES_TOC      8   /* Total OC (sum of HD1+HD2+HD3+HP1+HP2+HP3) - diagnostic */
#define CGEM_SPECIES_SPM      9
#define CGEM_SPECIES_DIC      10
#define CGEM_SPECIES_AT       11  /* Total Alkalinity [µeq/L] */
#define CGEM_SPECIES_PCO2     12  /* pCO2 [µatm] - diagnostic */
#define CGEM_SPECIES_CO2      13  /* Dissolved CO2 [µmol/L] - diagnostic */
#define CGEM_SPECIES_PH       14  /* pH [-] - diagnostic */
#define CGEM_SPECIES_HS       15  /* Hydrogen sulfide [µmol/L] */
#define CGEM_SPECIES_ALKC     16  /* Carbonate alkalinity iteration count - diagnostic */

/* RIVE multi-pool organic matter (17-22) */
#define CGEM_SPECIES_HD1      17  /* Labile dissolved OC [µmol C/L] */
#define CGEM_SPECIES_HD2      18  /* Semi-labile dissolved OC [µmol C/L] */
#define CGEM_SPECIES_HD3      19  /* Refractory dissolved OC [µmol C/L] */
#define CGEM_SPECIES_HP1      20  /* Labile particulate OC [µmol C/L] */
#define CGEM_SPECIES_HP2      21  /* Semi-labile particulate OC [µmol C/L] */
#define CGEM_SPECIES_HP3      22  /* Refractory particulate OC [µmol C/L] */

/* RIVE heterotrophic bacteria (23-24) */
#define CGEM_SPECIES_BAG      23  /* Attached (gros) bacteria [mg C/L] */
#define CGEM_SPECIES_BAP      24  /* Free (petit) bacteria [mg C/L] */

/* RIVE phosphorus (25) */
#define CGEM_SPECIES_PIP      25  /* Particulate inorganic P (adsorbed) [µmol P/L] */

/* RIVE dissolved substrates (26) */
#define CGEM_SPECIES_DSS      26  /* Dissolved simple substrates [mg C/L] - bacteria food */

/* ===========================================================================
 * GREENHOUSE GAS SPECIES (27-29) - NEW FROM C-RIVE
 * Reference: Garnier et al. (2007), Wang et al. (2018)
 * =========================================================================== */
#define CGEM_SPECIES_NO2      27  /* Nitrite [µmol N/L] - intermediate in nitrification */
#define CGEM_SPECIES_N2O      28  /* Nitrous oxide [nmol N/L] - GHG from nitrification/denit */
#define CGEM_SPECIES_CH4      29  /* Methane [µmol C/L] - GHG from methanogenesis */

/* ===========================================================================
 * REACTION INDICES
 * =========================================================================== */

#define CGEM_REACTION_NPP_NO3     0
#define CGEM_REACTION_NPP_NO3_1   1
#define CGEM_REACTION_NPP_NO3_2   2
#define CGEM_REACTION_NPP_NH4     3
#define CGEM_REACTION_NPP_NH4_1   4
#define CGEM_REACTION_NPP_NH4_2   5
#define CGEM_REACTION_GPP_1       6
#define CGEM_REACTION_GPP_2       7
#define CGEM_REACTION_NPP         8
#define CGEM_REACTION_PHY_DEATH   9
#define CGEM_REACTION_PHY_DEATH_1 10
#define CGEM_REACTION_PHY_DEATH_2 11
#define CGEM_REACTION_SI_CONS     12
#define CGEM_REACTION_AER_DEG     13
#define CGEM_REACTION_DENIT       14
#define CGEM_REACTION_NIT         15  /* NH4 → NO2 (nitrosation) */
#define CGEM_REACTION_O2_EX       16
#define CGEM_REACTION_O2_EX_S     17
#define CGEM_REACTION_CO2_EX      18
#define CGEM_REACTION_CO2_EX_S    19
#define CGEM_REACTION_DIC_REACT   20
#define CGEM_REACTION_TA_REACT    21
#define CGEM_REACTION_HS_REACT    22
#define CGEM_REACTION_EROSION_S   23
#define CGEM_REACTION_EROSION_V   24
#define CGEM_REACTION_DEPOSITION_S 25
#define CGEM_REACTION_DEPOSITION_V 26

/* RIVE additional reactions (27-39) */
#define CGEM_REACTION_HYDROLYSIS_HD1  27  /* HD1 → DSS hydrolysis */
#define CGEM_REACTION_HYDROLYSIS_HD2  28  /* HD2 → DSS hydrolysis */
#define CGEM_REACTION_HYDROLYSIS_HP1  29  /* HP1 → HD1 hydrolysis (particulate to dissolved) */
#define CGEM_REACTION_HYDROLYSIS_HP2  30  /* HP2 → HD2 hydrolysis */
#define CGEM_REACTION_BAC_UPTAKE      31  /* DSS uptake by bacteria */
#define CGEM_REACTION_BAC_RESP        32  /* Bacterial respiration */
#define CGEM_REACTION_BAC_MORT        33  /* Bacterial mortality */
#define CGEM_REACTION_P_ADSORPTION    34  /* PO4 ↔ PIP adsorption */
#define CGEM_REACTION_BENTHIC_RESP    35  /* Benthic OC respiration */
#define CGEM_REACTION_BENTHIC_NH4     36  /* Benthic NH4 flux */
#define CGEM_REACTION_BENTHIC_PO4     37  /* Benthic PO4 flux */
#define CGEM_REACTION_BENTHIC_O2      38  /* Benthic O2 demand (SOD) */
#define CGEM_REACTION_BENTHIC_DIC     39  /* Benthic DIC flux */

/* ===========================================================================
 * GHG REACTIONS (40-52) - NEW FROM C-RIVE
 * Reference: Garnier et al. (2007), Marescaux et al. (2019)
 * =========================================================================== */
#define CGEM_REACTION_NIT2            40  /* NO2 → NO3 (nitratation) */
#define CGEM_REACTION_N2O_NIT         41  /* N2O production from nitrification */
#define CGEM_REACTION_N2O_DENIT       42  /* N2O production from denitrification */
#define CGEM_REACTION_N2O_EX          43  /* N2O air-water exchange */
#define CGEM_REACTION_CH4_PROD        44  /* Methanogenesis */
#define CGEM_REACTION_CH4_OX          45  /* Methane oxidation (aerobic) */
#define CGEM_REACTION_CH4_OX_SO4      46  /* Methane oxidation (anaerobic - sulfate) */
#define CGEM_REACTION_CH4_EX          47  /* CH4 air-water exchange */
#define CGEM_REACTION_CH4_EBUL        48  /* CH4 ebullition */
#define CGEM_REACTION_CO2_BENTHIC     49  /* CO2 benthic flux */
#define CGEM_REACTION_CH4_BENTHIC     50  /* CH4 benthic flux */
#define CGEM_REACTION_N2O_BENTHIC     51  /* N2O benthic flux */
#define CGEM_REACTION_ANAMMOX         52  /* Anaerobic ammonium oxidation */

#define CGEM_NUM_SPECIES 30       /* Updated for GHG species */
#define CGEM_NUM_REACTIONS 53     /* Updated for GHG reactions */

/* ===========================================================================
 * SPECIES TRANSPORT FLAGS
 * env=1 means transport, env=0 means diagnostic only
 * Matching Fortran BGCArray(s)%env convention
 * =========================================================================== */
static const int CGEM_SPECIES_TRANSPORT_FLAG[CGEM_NUM_SPECIES] = {
    1,  /* SALINITY - transport */
    1,  /* PHY1 - transport */
    1,  /* PHY2 - transport */
    1,  /* DSI - transport */
    1,  /* NO3 - transport */
    1,  /* NH4 - transport */
    1,  /* PO4 - transport */
    1,  /* O2 - transport */
    0,  /* TOC - diagnostic only (computed from OC pools) */
    1,  /* SPM - transport */
    1,  /* DIC - transport */
    1,  /* AT (alkalinity) - transport */
    0,  /* PCO2 - diagnostic only (computed from carbonate system) */
    0,  /* CO2 - diagnostic only (computed from carbonate system) */
    0,  /* PH - diagnostic only (computed from carbonate system) */
    1,  /* HS - transport */
    0,  /* ALKC - diagnostic only (iteration counter) */
    /* RIVE multi-pool organic matter */
    1,  /* HD1 - labile dissolved OC - transport */
    1,  /* HD2 - semi-labile dissolved OC - transport */
    1,  /* HD3 - refractory dissolved OC - transport */
    1,  /* HP1 - labile particulate OC - transport */
    1,  /* HP2 - semi-labile particulate OC - transport */
    1,  /* HP3 - refractory particulate OC - transport */
    /* RIVE bacteria */
    1,  /* BAG - attached bacteria - transport */
    1,  /* BAP - free bacteria - transport */
    /* RIVE phosphorus */
    1,  /* PIP - particulate inorganic P - transport */
    /* RIVE substrates */
    1,  /* DSS - dissolved simple substrates - transport */
    /* GHG species */
    1,  /* NO2 - nitrite - transport */
    1,  /* N2O - nitrous oxide - transport */
    1   /* CH4 - methane - transport */
};

/* ===========================================================================
 * C-RIVE CARBONATE CHEMISTRY CONSTANTS
 * Reference: Zeebe & Wolf-Gladrow (2001), Marescaux et al. (2019)
 * =========================================================================== */

/* Atmospheric pCO2 [µatm] (current ~420 ppm) */
#define CGEM_PCO2_ATM_DEFAULT 420.0

/* Atmospheric N2O [ppb] (current ~335 ppb) */
#define CGEM_N2O_ATM_DEFAULT 335.0

/* Atmospheric CH4 [ppb] (current ~1900 ppb) */
#define CGEM_CH4_ATM_DEFAULT 1900.0

/* Default open-sea boundary condition distance [m] */
#define CGEM_DEFAULT_OSBC_DIST 10000.0

/* Utility macros */
#define CGEM_CLAMP(x, lo, hi) (((x) < (lo)) ? (lo) : (((x) > (hi)) ? (hi) : (x)))
#define CGEM_MAX(a, b) (((a) > (b)) ? (a) : (b))
#define CGEM_MIN(a, b) (((a) < (b)) ? (a) : (b))
#define CGEM_SQUARE(x) ((x) * (x))

#endif /* CGEM_DEFINE_H */
