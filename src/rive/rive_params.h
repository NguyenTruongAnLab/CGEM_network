/**
 * @file rive_params.h
 * @brief RIVE Biogeochemistry Parameters
 * 
 * Defines BiogeoParams structure and parameter loading functions.
 */

#ifndef RIVE_PARAMS_H
#define RIVE_PARAMS_H

#include "../network.h"

/* Seconds per day - for unit conversion */
#define RIVE_SECONDS_PER_DAY 86400.0

/**
 * Global biogeo parameters structure (loaded from config file)
 */
typedef struct {
    int loaded;
    
    /* Water and environmental */
    double water_temp;
    double ws;
    
    /* Light */
    double I0;
    double kd1;
    double kd2_spm;
    double kd2_phy1;
    double kd2_phy2;
    
    /* Phytoplankton */
    double alpha1, alpha2;
    double pbmax1, pbmax2;
    double kexc1, kexc2;
    double kgrowth1, kgrowth2;
    double kmaint1, kmaint2;
    double kmort1, kmort2;
    
    /* Nutrients */
    double kdsi1, kn1, kpo41;
    double kn2, kpo42;
    
    /* Decomposition - rates stored as [1/day] */
    double kox, kdenit, knit;
    double ktox, ko2, ko2_nit, kno3, knh4, kino2;
    
    /* ===========================================================================
     * 2-POOL TOC MODEL (SCIENTIFIC FIX - December 2025 Audit)
     * 
     * Replaces the salinity-based kox scaling hack with proper mechanistic model.
     * 
     * TOC_LABILE: Fast-degrading fraction (fresh phytoplankton, sewage)
     *   - kox_labile ~ 0.1-0.3 /day at 20°C (literature: Hopkinson 2005)
     *   - River boundary: ~20 µmol/L (10-15% of total)
     *   - Ocean boundary: ~5 µmol/L
     *   
     * TOC_REFRACTORY: Slow-degrading fraction (terrestrial humics, aged DOC)
     *   - kox_refractory ~ 0.005-0.02 /day at 20°C (Middelburg 1989)
     *   - River boundary: ~130 µmol/L (85-90% of total)
     *   - Ocean boundary: ~160 µmol/L
     *
     * Reference: Middelburg (1989), Hopkinson & Vallino (2005), Amon & Benner (1996)
     * ===========================================================================*/
    double kox_labile;            /* Labile TOC decay rate [1/day] (default 0.15) */
    double kox_refractory;        /* Refractory TOC decay rate [1/day] (default 0.008) */
    
    /* ===========================================================================
     * DECOUPLED BENTHIC FLUXES (SCIENTIFIC FIX - December 2025 Audit)
     * 
     * PROBLEM: Previous code coupled benthic O2 and CO2 1:1 (RQ=1).
     * In tropical estuaries with anaerobic sediments, RQ > 1 because:
     *   - Sulfate reduction produces CO2 without consuming O2
     *   - Denitrification produces CO2 with less O2 cost
     *   - Methanogenesis produces CO2 and CH4 without O2
     * 
     * SOLUTION: Separate SOD and DIC flux parameters
     *   benthic_SOD: Sediment Oxygen Demand [mmol O2/m²/day]
     *   benthic_DIC_flux: DIC (CO2) release from sediments [mmol C/m²/day]
     *   RQ_benthic = benthic_DIC_flux / benthic_SOD (typically 1.0-3.0)
     *
     * For aerobic-dominated sediments: RQ ~ 1.0
     * For anaerobic sediments (sulfate reduction): RQ ~ 1.5-2.0
     * For highly reduced sediments: RQ ~ 2.0-3.0
     *
     * Reference: Cai (2011), Borges & Abril (2011), Middelburg et al. (2005)
     * ===========================================================================*/
    double benthic_SOD;           /* Sediment O2 demand [mmol O2/m²/day] (default = benthic_resp_20C) */
    double benthic_DIC_flux;      /* Sediment DIC release [mmol C/m²/day] (default = SOD * RQ) */
    double RQ_benthic;            /* Respiratory quotient for benthic [mol CO2/mol O2] (default 1.5) */
    
    /* Stoichiometry */
    double redn, redp, redsi;
    
    /* Gas exchange */
    double pco2_atm;
    
    /* Wind-driven gas exchange */
    double wind_speed;
    double wind_coeff;
    double schmidt_exp;
    double current_k_factor;
    
    /* Benthic respiration */
    double benthic_resp_20C;
    double benthic_Q10;
    
    /* Spatially-varying benthic flux (December 2025 Audit Fix)
     * Benthic fluxes vary along estuaries due to sediment type:
     * - Ocean mouth: Sandy sediments → LOW benthic flux
     * - Upstream: Fine organic sediments → HIGH benthic flux
     * 
     * The scaling factors multiply benthic_resp_20C:
     *   benthic_rate(x) = benthic_resp_20C * scale_factor(x)
     * where scale_factor = ocean_scale at x=0, upstream_scale at x=L
     * 
     * Reference: Abril et al. (2010), Cai (2011)
     */
    double benthic_ocean_scale;     /* Scaling at ocean mouth (default 0.3) */
    double benthic_upstream_scale;  /* Scaling at upstream end (default 2.0) */
    double benthic_S_high;          /* Salinity threshold for ocean scale (default 15) */
    double benthic_S_low;           /* Salinity threshold for upstream scale (default 2) */
    
    /* RIVE organic matter */
    double khydr1, khydr2, khydr3;
    double frac_hd1, frac_hd2;
    double frac_hp1, frac_hp2;
    
    /* RIVE bacteria */
    double bag_bmax20, bap_bmax20;
    double bag_kdb20, bap_kdb20;
    double bac_ks, bac_yield;
    double bag_topt, bap_topt;
    double bag_dti, bap_dti;
    double bag_vs;
    
    /* Phosphorus adsorption (Langmuir isotherm) */
    double kpads, pac;
    
    /* ===========================================================================
     * SEDIMENT-BIOGEOCHEMISTRY COUPLING (From C-RIVE)
     * Reference: Wang et al. (2018) Water Research 144, 341-355
     * ===========================================================================*/
    
    /* POC/SPM ratio: Fraction of TOC that settles/erodes with mud */
    double poc_spm_ratio;         /* Typically 0.01-0.03 (1-3%) for Mekong */
    
    /* Benthic */
    double zf_init;
    double benthic_porosity;
    double benthic_density;
    
    /* GHG parameters */
    double N2O_atm_ppb;
    double N2O_yield_nit, N2O_yield_denit;
    double CH4_atm_ppb;
    double CH4_prod_rate, CH4_ox_rate;
    double CH4_ks_o2, CH4_ebul_thresh;
    double k_nit1, k_nit2, ks_no2;
    
    /* Benthic GHG fluxes (for passive mode) */
    double benthic_CH4_flux;      /* CH4 sediment flux [µmol/m²/day] */
    double benthic_N2O_flux;      /* N2O sediment flux [nmol/m²/day] */
    
    /* Solver settings */
    int use_rk4;
    
    /* ===========================================================================
     * TROPICAL ESTUARY PHYSICS (Critical for Mekong Delta)
     * Reference: Audit recommendations for "Academic Solidity"
     * ===========================================================================*/
    
    /* Salinity stress on freshwater phytoplankton */
    double sal_stress_thresh;     /* Salinity threshold for stress onset [PSU] (default 5) */
    double sal_stress_coef;       /* Mortality increase coefficient (default 0.5) */
    double sal_marine_opt;        /* Optimal salinity for marine species [PSU] (default 30) */
    
    /* Flocculation - salinity-dependent settling velocity */
    double floc_sal_scale;        /* Salinity scale for flocculation [PSU] (default 2) */
    double floc_factor_max;       /* Maximum flocculation factor (default 10) */
    int enable_flocculation;      /* 1 = Enable flocculation physics */
    
    /* Tidal pumping for benthic fluxes */
    double tidal_pump_coef;       /* Tidal pumping enhancement coefficient */
    double tidal_pump_vel_ref;    /* Reference velocity for tidal pumping [m/s] */
    
    /* ===========================================================================
     * SAFETY MODE FLAGS (80/20 Simplification for Data-Limited Systems)
     * ===========================================================================
     * When ghg_passive_mode = 1 (default, SAFE):
     *   - GHG module is a "listener" that reads core rates but does NOT modify
     *     core state variables (O2, NH4, NO3)
     *   - N2O is calculated as yield fraction of core nitrification/denitrification
     *   - CH4 dynamics are independent (don't feed back to O2)
     *   - Prevents double-counting and mass balance errors
     * 
     * When ghg_passive_mode = 0 (ACTIVE, RISKY):
     *   - GHG module applies feedback to core variables
     *   - Only use if you have calibrated GHG parameters with field data
     * ===========================================================================*/
    int ghg_passive_mode;         /* 1 = Diagnostic only (Safe), 0 = Active feedback (Risky) */
    int simplified_mode;          /* 1 = 80/20 mode (skip bacteria, multi-pool OC) */
    int skip_bacteria;            /* 1 = Use bulk kox instead of explicit bacteria */
    int skip_multipool_oc;        /* 1 = Use single TOC instead of HD1-3, HP1-3 */
    int skip_ghg_dynamics;        /* 1 = Only calculate GHG as post-processing */
    int skip_p_adsorption;        /* 1 = Skip dynamic PO4-PIP equilibrium */
    int skip_carbonate_reactions; /* 1 = Skip DIC/TA/pCO2 reactions (transport only) - Dec 2025 */
    
    /* ===========================================================================
     * SIMPLIFIED MODE KINETICS (BOD/COD-based for data-sparse regions)
     * When simplified_mode = 1:
     *   - Bypasses bacterial dynamics (BAG, BAP)
     *   - Uses temperature-corrected first-order decay from BOD/COD
     *   - Maps monitoring data directly to model state
     * ===========================================================================*/
    double theta_ox;              /* Arrhenius temp coefficient for oxidation (default 1.047) */
    double theta_nit;             /* Arrhenius temp coefficient for nitrification */
    double theta_denit;           /* Arrhenius temp coefficient for denitrification */
    
    /* ===========================================================================
     * ACTIVE SEDIMENT LAYER (SOC Pool) - December 2025 Audit Priority #1
     * 
     * This feature makes benthic fluxes DYNAMIC instead of fixed!
     * SOC pool receives POC deposition and releases benthic fluxes proportional
     * to accumulated organic carbon. This captures "legacy load" effects.
     * 
     * Without this: If you reduce pollution by 50%, sediment still releases
     * the same fluxes forever - WRONG for scenario analysis!
     * 
     * With SOC: Benthic fluxes decrease as sediment pool depletes - CORRECT!
     * 
     * Reference: Chapra (2008), DiToro (2001), Soetaert et al. (1996)
     * ===========================================================================*/
    int enable_soc;               /* 1 = Enable dynamic SOC, 0 = fixed benthic fluxes */
    double k_soc_20C;             /* SOC mineralization rate at 20°C [1/day] (default 0.003) */
    double soc_Q10;               /* Q10 for SOC decay (default 2.0) */
    double soc_init;              /* Initial SOC pool [g C/m²] (default 500) */
    double soc_max;               /* Maximum SOC capacity [g C/m²] (default 5000) */
    double k_burial;              /* SOC burial rate [1/day] (default 0.0001) */
    double soc_f_anaerobic;       /* Anaerobic fraction of mineralization (default 0.3) */
    double soc_ch4_yield;         /* CH4 yield from anaerobic decay [mol/mol C] (default 0.5) */
    double soc_n2o_yield;         /* N2O yield from sediment N cycling (default 0.02) */
    
} BiogeoParams;

/* Global params accessor */
BiogeoParams* rive_get_params(void);

/* Load parameters from file */
int LoadBiogeoParams(const char *path);

/* Initialize branch parameters from global */
void InitializeBiogeoParameters(Branch *branch);

#endif /* RIVE_PARAMS_H */
