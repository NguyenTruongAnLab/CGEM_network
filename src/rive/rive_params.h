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
    
    /* Phosphorus adsorption */
    double kpads, pac;
    
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
    
} BiogeoParams;

/* Global params accessor */
BiogeoParams* rive_get_params(void);

/* Load parameters from file */
int LoadBiogeoParams(const char *path);

/* Initialize branch parameters from global */
void InitializeBiogeoParameters(Branch *branch);

#endif /* RIVE_PARAMS_H */
