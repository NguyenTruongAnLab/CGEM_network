/**
 * @file ghg_module.h
 * @brief C-RIVE Greenhouse Gas Module for CGEM
 * 
 * Complete GHG module implementing CO2, CH4, and N2O dynamics including:
 * - CO2: Air-water exchange, carbonate chemistry (see carbonate_chem.h)
 * - N2O: Production from nitrification and denitrification, air-water exchange
 * - CH4: Methanogenesis, aerobic/anaerobic oxidation, ebullition, air-water exchange
 * 
 * Ported from: C-RIVE rea_degassing.c, growth.c (Unified RIVE v1.0)
 * 
 * References:
 * - Garnier et al. (2007) N2O production from nitrification and denitrification
 * - Cébron et al. (2005) N2O from nitrification in rivers
 * - Wang et al. (2018) C-RIVE sensitivity analysis
 * - Borges & Abril (2011) CH4 in estuaries
 * - DelSontro et al. (2018) CH4 ebullition in lakes
 * 
 * CITATION for C-RIVE:
 * Wang, S., Flipo, N., Romary, T., 2018. Time-dependent global sensitivity analysis
 * of the C-RIVE biogeochemical model. Water Research 144, 341-355.
 */

#ifndef CGEM_GHG_MODULE_H
#define CGEM_GHG_MODULE_H

#include "../define.h"
#include <math.h>

/* ===========================================================================
 * GHG Constants from C-RIVE
 * ===========================================================================*/

/* Atmospheric concentrations (defaults) */
#define CRIVE_N2O_ATM_PPB     335.0    /* Atmospheric N2O [ppb] (current ~335) */
#define CRIVE_CH4_ATM_PPB     1900.0   /* Atmospheric CH4 [ppb] (current ~1900) */

/* Kelvin offset */
#define CRIVE_T_KELVIN_GHG    273.15

/* N2O production factors from C-RIVE param_RIVE.h */
/* PROD_N2O: specific rate of N2O production
 * ~ 4 gN/gC for nitrification (Cébron et al. 2005)
 * ~ 0.01 molN2O/molNO3 for denitrification (Billen, Modifications à Sénèque)
 */
#define CRIVE_N2O_YIELD_NIT     0.004   /* gN-N2O / gC-bacteria during nitrification */
#define CRIVE_N2O_YIELD_DENIT   0.01    /* mol N2O / mol NO3 reduced during denitrification */

/* CH4 production/oxidation parameters (estuarine defaults) */
#define CRIVE_CH4_PROD_RATE     0.01    /* Base methanogenesis rate [/day] */
#define CRIVE_CH4_OX_RATE       0.5     /* CH4 aerobic oxidation rate [/day] */
#define CRIVE_CH4_KS_O2         10.0    /* O2 half-sat for CH4 oxidation [µmol/L] */
#define CRIVE_CH4_KS_SO4        500.0   /* SO4 half-sat for anaerobic CH4 ox [µmol/L] */
#define CRIVE_CH4_EBUL_THRESH   50.0    /* CH4 threshold for ebullition [µmol/L] */

/* ===========================================================================
 * GHG State Structure
 * ===========================================================================*/

/**
 * Complete GHG state at a single point
 * All concentrations in model units (see comments)
 */
typedef struct {
    /* N2O state */
    double N2O;             /* Dissolved N2O [nmol N/L] or [nmol/L] */
    double N2O_sat;         /* N2O saturation concentration [nmol/L] */
    double N2O_flux;        /* N2O air-water flux [nmol/L/s] (positive = evasion) */
    double N2O_from_nit;    /* N2O production from nitrification [nmol/L/s] */
    double N2O_from_denit;  /* N2O production from denitrification [nmol/L/s] */
    double Sc_N2O;          /* Schmidt number for N2O */
    
    /* CH4 state */
    double CH4;             /* Dissolved CH4 [µmol C/L] or [µmol/L] */
    double CH4_sat;         /* CH4 saturation concentration [µmol/L] */
    double CH4_flux;        /* CH4 air-water flux [µmol/L/s] (positive = evasion) */
    double CH4_prod;        /* CH4 production from methanogenesis [µmol/L/s] */
    double CH4_ox;          /* CH4 aerobic oxidation [µmol/L/s] */
    double CH4_ox_anaer;    /* CH4 anaerobic oxidation (sulfate) [µmol/L/s] */
    double CH4_ebul;        /* CH4 ebullition flux [µmol/L/s] */
    double Sc_CH4;          /* Schmidt number for CH4 */
    
    /* NO2 (intermediate in 2-step nitrification) */
    double NO2;             /* Nitrite [µmol N/L] */
    
} GHGState;

/**
 * GHG configuration parameters
 */
typedef struct {
    /* N2O parameters */
    double N2O_atm_ppb;         /* Atmospheric N2O [ppb] */
    double N2O_yield_nit;       /* N2O yield from nitrification [mol N2O / mol NH4] */
    double N2O_yield_denit;     /* N2O yield from denitrification [mol N2O / mol NO3] */
    double kmich_no2_n2o;       /* Michaelis constant for NO2 in N2O production */
    
    /* CH4 parameters */
    double CH4_atm_ppb;         /* Atmospheric CH4 [ppb] */
    double CH4_prod_rate;       /* Methanogenesis rate at 20°C [/day] */
    double CH4_ox_rate;         /* Aerobic oxidation rate at 20°C [/day] */
    double CH4_ks_o2;           /* O2 half-saturation for CH4 oxidation [µmol/L] */
    double CH4_ks_so4;          /* SO4 half-saturation for anaerobic CH4 ox [µmol/L] */
    double CH4_ebul_thresh;     /* CH4 threshold for ebullition [µmol/L] */
    double CH4_ebul_rate;       /* Ebullition rate constant [/day] */
    
    /* 2-step nitrification parameters (from C-RIVE growth_bactn.c) */
    double k_nit1;              /* NH4 → NO2 rate at 20°C [/day] */
    double k_nit2;              /* NO2 → NO3 rate at 20°C [/day] */
    double ks_nh4;              /* NH4 half-saturation for nitrification [µmol/L] */
    double ks_no2;              /* NO2 half-saturation for nitratation [µmol/L] */
    double ks_o2_nit;           /* O2 half-saturation for nitrification [µmol/L] */
    
    /* Temperature coefficients */
    double Q10_nit;             /* Q10 for nitrification (typical 1.08) */
    double Q10_ch4_prod;        /* Q10 for methanogenesis (typical 2.0-3.0) */
    double Q10_ch4_ox;          /* Q10 for CH4 oxidation (typical 2.0) */
    
    /* Benthic fluxes */
    double benthic_CH4_flux;    /* Sediment CH4 flux [µmol/m²/day] */
    double benthic_N2O_flux;    /* Sediment N2O flux [nmol/m²/day] */
    
} GHGConfig;

/* ===========================================================================
 * Function Prototypes - N2O Module (from C-RIVE rea_degassing.c, growth_bactn.c)
 * ===========================================================================*/

/**
 * Initialize default GHG configuration
 * @param config Pointer to configuration structure to initialize
 */
void crive_ghg_init_config(GHGConfig *config);

/**
 * Calculate N2O saturation concentration
 * Weiss & Price (1980) formulation
 * 
 * Uses Bunsen solubility β [mol/(L·atm)] with atmospheric N2O = 335 ppb.
 * C_sat = β × p_N2O × 10⁹ [nmol/L]
 * 
 * Typical values:
 *   - At 25°C, S=0:  ~7 nmol/L
 *   - At 25°C, S=35: ~5 nmol/L
 * 
 * Reference: Weiss, R.F., Price, B.A. (1980) Marine Chemistry, 8, 347-359.
 * 
 * @param temp Water temperature [°C]
 * @param salinity Salinity [PSU]
 * @return N2O saturation concentration [nmol/L] (matches GHGState.N2O units)
 */
double crive_calc_n2o_sat(double temp, double salinity);

/**
 * Calculate N2O Schmidt number
 * From C-RIVE param_RIVE.h / rea_degassing.c
 * 
 * Sc = poly3(a, b, c, d, T)
 * 
 * @param temp Water temperature [°C]
 * @return Schmidt number for N2O [-]
 */
double crive_calc_n2o_schmidt(double temp);

/**
 * Calculate N2O air-water flux
 * Direct port from C-RIVE rea_degassing_N2O()
 * 
 * Flux = kg * (C - Csat) / depth [µmol/L/s]
 * where kg = 1.719 * sqrt(600*v / (Sc*depth))
 * 
 * DECEMBER 2025 AUDIT FIX: Now returns µmol/L/s (was nmol/L/s)
 * to match biogeo.c internal units. The model stores N2O in µmol/L.
 * 
 * @param N2O_conc Current dissolved N2O [µmol/L] (model internal units)
 * @param N2O_sat N2O saturation concentration [nmol/L] (from crive_calc_n2o_sat)
 * @param depth Water depth [m]
 * @param velocity Water velocity [m/s]
 * @param Sc Schmidt number [-]
 * @return N2O flux [µmol/L/s] (positive = evasion to atmosphere)
 */
double crive_calc_n2o_flux(double N2O_conc, double N2O_sat, double depth,
                           double velocity, double Sc);

/**
 * Calculate N2O production from nitrification
 * Based on Cébron et al. (2005) and C-RIVE growth_bactn.c
 * 
 * N2O_nit = yield_nit * nitrification_rate
 * 
 * @param nitrif_rate Nitrification rate [µmol N/L/s]
 * @param yield N2O yield factor [mol N2O-N / mol NH4-N]
 * @return N2O production rate [nmol N/L/s]
 */
double crive_calc_n2o_from_nitrification(double nitrif_rate, double yield);

/**
 * Calculate N2O production from denitrification
 * Based on Billen modifications to Sénèque and C-RIVE
 * 
 * N2O_denit = yield_denit * denitrification_rate
 * 
 * @param denit_rate Denitrification rate [µmol N/L/s]
 * @param yield N2O yield factor [mol N2O / mol NO3]
 * @return N2O production rate [nmol N/L/s]
 */
double crive_calc_n2o_from_denitrification(double denit_rate, double yield);

/* ===========================================================================
 * Function Prototypes - CH4 Module
 * Based on Borges & Abril (2011), DelSontro et al. (2018)
 * ===========================================================================*/

/**
 * Calculate CH4 saturation concentration
 * Yamamoto et al. (1976) formulation
 * 
 * @param temp Water temperature [°C]
 * @param salinity Salinity [PSU]
 * @param CH4_atm Atmospheric CH4 [ppb]
 * @return CH4 saturation concentration [µmol/L]
 */
double crive_calc_ch4_sat(double temp, double salinity, double CH4_atm);

/**
 * Calculate CH4 Schmidt number
 * Wanninkhof (2014) coefficients
 * 
 * @param temp Water temperature [°C]
 * @return Schmidt number for CH4 [-]
 */
double crive_calc_ch4_schmidt(double temp);

/**
 * Calculate CH4 air-water flux
 * Similar formulation to CO2 and N2O
 * 
 * @param CH4_conc Current dissolved CH4 [µmol/L]
 * @param CH4_sat CH4 saturation concentration [µmol/L]
 * @param depth Water depth [m]
 * @param k600 Gas transfer velocity at Sc=600 [m/s]
 * @param Sc Schmidt number [-]
 * @return CH4 flux [µmol/L/s] (positive = evasion to atmosphere)
 */
double crive_calc_ch4_flux(double CH4_conc, double CH4_sat, double depth,
                           double k600, double Sc);

/**
 * Calculate CH4 production from methanogenesis
 * Temperature-dependent production in anoxic sediments
 * 
 * @param TOC_benthic Benthic organic carbon [mg C/m²]
 * @param depth Water depth [m]
 * @param temp Water temperature [°C]
 * @param O2_conc Oxygen concentration [µmol/L]
 * @param config GHG configuration
 * @return CH4 production rate [µmol/L/s]
 */
double crive_calc_ch4_production(double TOC_benthic, double depth, double temp,
                                  double O2_conc, const GHGConfig *config);

/**
 * Calculate aerobic CH4 oxidation
 * Michaelis-Menten kinetics with O2 limitation
 * 
 * @param CH4_conc CH4 concentration [µmol/L]
 * @param O2_conc O2 concentration [µmol/L]
 * @param temp Water temperature [°C]
 * @param config GHG configuration
 * @return CH4 oxidation rate [µmol/L/s]
 */
double crive_calc_ch4_oxidation(double CH4_conc, double O2_conc, double temp,
                                 const GHGConfig *config);

/**
 * Calculate anaerobic CH4 oxidation (sulfate-mediated)
 * Important in marine/estuarine sediments
 * 
 * @param CH4_conc CH4 concentration [µmol/L]
 * @param SO4_conc Sulfate concentration [µmol/L]
 * @param temp Water temperature [°C]
 * @param config GHG configuration
 * @return Anaerobic CH4 oxidation rate [µmol/L/s]
 */
double crive_calc_ch4_oxidation_anaerobic(double CH4_conc, double SO4_conc, 
                                           double temp, const GHGConfig *config);

/**
 * Calculate CH4 ebullition (bubble release)
 * Triggered when CH4 exceeds solubility threshold
 * DelSontro et al. (2018) parameterization
 * 
 * @param CH4_conc CH4 concentration [µmol/L]
 * @param CH4_sat CH4 saturation [µmol/L]
 * @param depth Water depth [m]
 * @param config GHG configuration
 * @return CH4 ebullition rate [µmol/L/s]
 */
double crive_calc_ch4_ebullition(double CH4_conc, double CH4_sat, double depth,
                                  const GHGConfig *config);

/* ===========================================================================
 * Function Prototypes - 2-Step Nitrification (from C-RIVE growth_bactn.c)
 * NH4 → NO2 (nitrosation) and NO2 → NO3 (nitratation)
 * ===========================================================================*/

/**
 * Calculate nitrosation rate (NH4 → NO2)
 * First step of 2-step nitrification
 * 
 * @param NH4_conc Ammonium concentration [µmol/L]
 * @param O2_conc Oxygen concentration [µmol/L]
 * @param temp Water temperature [°C]
 * @param config GHG configuration
 * @return Nitrosation rate [µmol N/L/s]
 */
double crive_calc_nitrosation(double NH4_conc, double O2_conc, double temp,
                               const GHGConfig *config);

/**
 * Calculate nitratation rate (NO2 → NO3)
 * Second step of 2-step nitrification
 * 
 * @param NO2_conc Nitrite concentration [µmol/L]
 * @param O2_conc Oxygen concentration [µmol/L]
 * @param temp Water temperature [°C]
 * @param config GHG configuration
 * @return Nitratation rate [µmol N/L/s]
 */
double crive_calc_nitratation(double NO2_conc, double O2_conc, double temp,
                               const GHGConfig *config);

/* ===========================================================================
 * Integrated GHG Calculation
 * ===========================================================================*/

/**
 * Calculate all GHG fluxes and transformations
 * Main entry point for GHG module
 * 
 * @param state GHG state structure (input: current concentrations, output: rates)
 * @param temp Water temperature [°C]
 * @param salinity Salinity [PSU]
 * @param depth Water depth [m]
 * @param velocity Water velocity [m/s]
 * @param O2_conc Oxygen concentration [µmol/L]
 * @param NO3_conc Nitrate concentration [µmol/L]
 * @param NH4_conc Ammonium concentration [µmol/L]
 * @param TOC_conc Total organic carbon [µmol C/L]
 * @param denit_rate Denitrification rate [µmol N/L/s]
 * @param k600 Gas transfer velocity at Sc=600 [m/s]
 * @param config GHG configuration
 * @return 0 on success, -1 on error
 */
int crive_calc_ghg_system(GHGState *state, double temp, double salinity,
                          double depth, double velocity, double O2_conc,
                          double NO3_conc, double NH4_conc, double TOC_conc,
                          double denit_rate, double k600, const GHGConfig *config);

#endif /* CGEM_GHG_MODULE_H */
