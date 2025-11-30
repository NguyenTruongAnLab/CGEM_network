/**
 * @file carbonate_chem.h
 * @brief C-RIVE Carbonate Chemistry Module for CGEM
 * 
 * Full implementation of the C-RIVE unified inorganic carbon module including:
 * - Carbonate equilibrium (DIC ↔ CO2 ↔ HCO3- ↔ CO3²-)
 * - pH calculation from DIC and Total Alkalinity
 * - CO2 solubility (Henry's law with temperature and salinity dependence)
 * - Air-water CO2 exchange with k600 parameterization
 * 
 * Ported from: C-RIVE calc_co2.c, calc_param.c (Unified RIVE v1.0)
 * 
 * References:
 * - Wang et al. (2018) Water Research 144, 341-355
 * - Hasanyar et al. (2022) Biogeosciences
 * - Marescaux et al. (2019) Carbonate speciation
 * - Weiss (1974) Marine Chemistry - CO2 solubility
 * - Zeebe & Wolf-Gladrow (2001) CO2 in Seawater
 * 
 * CITATION for C-RIVE:
 * Wang, S., Flipo, N., Romary, T., 2018. Time-dependent global sensitivity analysis
 * of the C-RIVE biogeochemical model. Water Research 144, 341-355.
 */

#ifndef CGEM_CARBONATE_CHEM_H
#define CGEM_CARBONATE_CHEM_H

#include "../define.h"
#include <math.h>

/* ===========================================================================
 * C-RIVE Constants
 * ===========================================================================*/

#define CRIVE_T_KELVIN 273.15       /* Celsius to Kelvin offset */
#define CRIVE_GC_MMOL 0.012         /* g C per mmol (12/1000) */
#define CRIVE_RHO_WATER 1000000.0   /* Water density [g/m³] */

/* Default dissociation constant for boric acid (freshwater) */
#define CRIVE_KB_DEFAULT 5.81e-10   /* [mol/L] */

/* Default total boron (freshwater = 0, seawater = 0.41 mM) */
#define CRIVE_BOR0_DEFAULT 0.0      /* [mM] */

/* ===========================================================================
 * Carbonate Chemistry Parameters Structure
 * ===========================================================================*/

/**
 * Carbonate system state at a single point
 * All concentrations in [mmol/m³] = [µmol/L]
 */
typedef struct {
    double DIC;         /* Dissolved Inorganic Carbon [mmol/m³] */
    double TA;          /* Total Alkalinity [mmol/m³] = [µeq/L] */
    double pH;          /* pH [-] */
    double CO2;         /* Dissolved CO2 [mmol/m³] */
    double HCO3;        /* Bicarbonate [mmol/m³] */
    double CO3;         /* Carbonate [mmol/m³] */
    double pCO2;        /* Partial pressure of CO2 [µatm] */
    double Csat;        /* Saturation concentration [mmol/m³] */
    double k0;          /* CO2 solubility [mmol/m³/µatm] */
    double k600;        /* Gas transfer velocity for Sc=600 [m/s] */
    double Sc;          /* Schmidt number [-] */
} CarbonateState;

/**
 * k600 calculation method (from C-RIVE, extended for estuaries)
 * 
 * For tidal estuaries (Mekong, Saigon, Scheldt), use K600_ESTUARINE (Abril et al. 2009)
 * For rivers, use K600_STRAHLER (Raymond et al. 2012)
 */
typedef enum {
    K600_RESERVOIRS = 0,    /* Alin et al. (2011) for reservoirs */
    K600_STRAHLER = 1,      /* Varies with Strahler order (rivers) */
    K600_USER = 2,          /* User-defined constant */
    K600_ESTUARINE = 3,     /* Abril et al. (2009) for tidal estuaries */
    K600_BORGES = 4         /* Borges et al. (2004) wind-driven estuarine */
} K600Method;

/**
 * Carbonate chemistry configuration parameters
 */
typedef struct {
    double kb;              /* Boric acid dissociation constant [mol/L] */
    double bor0;            /* Total boron concentration [mM] */
    double pCO2_atm;        /* Atmospheric pCO2 [µatm] */
    K600Method k600_method; /* Method for k600 calculation */
    double k600_user;       /* User-defined k600 [m/s] if method = K600_USER */
    int strahler_order;     /* Strahler stream order (1-8) */
    double wind_speed;      /* Wind speed at 10m [m/s] for estuarine methods */
    double fetch;           /* Fetch length [m] for wind-driven exchange */
} CarbonateConfig;

/* ===========================================================================
 * Function Prototypes - C-RIVE Carbonate Chemistry
 * ===========================================================================*/

/**
 * Initialize default carbonate configuration
 * @param config Pointer to configuration structure to initialize
 */
void crive_carbonate_init_config(CarbonateConfig *config);

/**
 * Calculate carbonate system equilibrium from DIC and TA
 * Ported directly from C-RIVE calc_CO2() and calc_pH()
 * 
 * @param state Pointer to carbonate state (DIC and TA must be set, others computed)
 * @param temp Water temperature [°C]
 * @param salinity Salinity [PSU]
 * @param config Carbonate configuration parameters
 * @return 0 on success, -1 on error
 */
int crive_calc_carbonate_system(CarbonateState *state, double temp, 
                                 double salinity, const CarbonateConfig *config);

/**
 * Calculate pH from DIC and TA using analytical Cardano solution
 * Direct port from C-RIVE calc_pH()
 * 
 * Solves the cubic equation for [H+] using the Cardano-Vieta formula:
 *   X = (1 - bor0/TA) * kb + (1 - DIC/TA) * K1
 *   Y = (1 - (bor0 + DIC)/TA) * K1*kb + (1 - 2*DIC/TA) * K1*K2
 *   Z = (1 - (bor0 + 2*DIC)/TA) * K1*K2*kb
 * 
 * @param DIC_mmol Dissolved inorganic carbon [mmol/m³]
 * @param TA_mmol Total alkalinity [mmol/m³]
 * @param temp Water temperature [°C]
 * @param kb Boric acid dissociation constant [mol/L]
 * @param bor0 Total boron [mM]
 * @return Hydrogen ion activity [mol/L]
 */
double crive_calc_pH(double DIC_mmol, double TA_mmol, double temp,
                     double kb, double bor0);

/**
 * Calculate CO2 solubility (Henry's law constant k0)
 * Direct port from C-RIVE calc_co2_solubility()
 * 
 * Uses Weiss (1974) formulation:
 *   ln(k0) = A1 + A2*(100/T) + A3*ln(T/100) + S*(B1 + B2*(T/100) + B3*(T/100)²)
 * 
 * @param temp_K Water temperature [K]
 * @param salinity Salinity [PSU]
 * @return CO2 solubility [mol/kg/atm] = [mmol/g/atm]
 */
double crive_calc_co2_solubility(double temp_K, double salinity);

/**
 * Calculate CO2 saturation concentration in water
 * Direct port from C-RIVE calc_co2_sat()
 * 
 * Csat = k0 * xCO2 * 1e-6 * rho_water
 * where xCO2 is dry air CO2 converted from wet air pCO2
 * 
 * @param temp_water Water temperature [°C]
 * @param temp_air Air temperature [°C]
 * @param pCO2_atm Atmospheric pCO2 [µatm]
 * @param salinity Salinity [PSU]
 * @return CO2 saturation concentration [mmol/m³]
 */
double crive_calc_co2_sat(double temp_water, double temp_air, 
                          double pCO2_atm, double salinity);

/**
 * Calculate Schmidt number for CO2
 * From Wanninkhof (1992) / Marescaux et al. (2019)
 * 
 * Sc = A - B*T + C*T² - D*T³
 * 
 * @param temp Water temperature [°C]
 * @return Schmidt number [-]
 */
double crive_calc_co2_schmidt(double temp);

/**
 * Calculate gas transfer velocity k600 for CO2
 * Direct port from C-RIVE calc_k600_co2()
 * 
 * Methods:
 * - RESERVOIRS: k600 = (13.82 + 0.35*v*100) / 100 [m/h] (Alin et al. 2011)
 * - STRAHLER: Varies with stream order (O'Connor-Dobbins for large rivers)
 * - USER: Constant user-defined value
 * 
 * @param velocity Water velocity [m/s]
 * @param depth Water depth [m]
 * @param method k600 calculation method
 * @param strahler Strahler stream order (used if method = K600_STRAHLER)
 * @param k600_user User-defined k600 [m/s] (used if method = K600_USER)
 * @return Gas transfer velocity k600 [m/s]
 */
double crive_calc_k600(double velocity, double depth, K600Method method,
                       int strahler, double k600_user);

/**
 * Calculate CO2 air-water flux
 * Direct port from C-RIVE rea_degassing_CO2()
 * 
 * Flux = k * (Csat - C) / depth [mmol/m³/s]
 * where k = k600 * sqrt(600/Sc) [m/s]
 * 
 * @param CO2_conc Current dissolved CO2 [mmol/m³]
 * @param CO2_sat CO2 saturation concentration [mmol/m³]
 * @param depth Water depth [m]
 * @param k600 Gas transfer velocity at Sc=600 [m/s]
 * @param Sc Schmidt number [-]
 * @return CO2 flux [mmol/m³/s] (positive = into water)
 */
double crive_calc_co2_flux(double CO2_conc, double CO2_sat, double depth,
                           double k600, double Sc);

/**
 * Calculate water density
 * Direct port from C-RIVE calc_water_density()
 * Uses International One Atmosphere Equation (Millero & Poisson 1981)
 * 
 * @param temp_water Water temperature [°C]
 * @param salinity Salinity [PSU]
 * @param pressure Applied pressure [bar] (0 for surface)
 * @return Water density [g/m³]
 */
double crive_calc_water_density(double temp_water, double salinity, double pressure);

/**
 * Dry-wet air conversion for CO2
 * Direct port from C-RIVE calc_co2_dry_wet_conversion()
 * Weiss & Price (1980) formulation
 * 
 * @param temp_air Air temperature [°C]
 * @param pCO2_atm Atmospheric pCO2 in wet air [µatm]
 * @param salinity Salinity [PSU]
 * @return xCO2 in dry air [µatm]
 */
double crive_calc_co2_dry_wet(double temp_air, double pCO2_atm, double salinity);

/* ===========================================================================
 * Carbonate Dissociation Constants
 * ===========================================================================*/

/**
 * First carbonate dissociation constant K1
 * H2CO3 ↔ H+ + HCO3-
 * Uses Harned & Davis (1943) for freshwater
 * 
 * pK1 = -126.34048 + 6320.813/T + 19.568224*ln(T)
 * 
 * @param temp_K Temperature [K]
 * @return K1 [mol/L]
 */
double crive_calc_K1(double temp_K);

/**
 * Second carbonate dissociation constant K2
 * HCO3- ↔ H+ + CO3²-
 * Uses Harned & Scholes (1941) for freshwater
 * 
 * pK2 = -90.18333 + 5143.692/T + 14.613358*ln(T)
 * 
 * @param temp_K Temperature [K]
 * @return K2 [mol/L]
 */
double crive_calc_K2(double temp_K);

#endif /* CGEM_CARBONATE_CHEM_H */
