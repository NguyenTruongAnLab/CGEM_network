/**
 * @file rive_params.c
 * @brief RIVE Biogeochemistry Parameter Loading
 * 
 * Handles loading and initialization of biogeochemical parameters.
 */

#include "rive_params.h"
#include "../define.h"
#include "ghg_module.h"
#include "carbonate_chem.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

/* Global biogeo parameters */
static BiogeoParams g_biogeo_params = {0};

/* Global configurations */
static GHGConfig g_ghg_config = {0};
static int g_ghg_config_initialized = 0;
static CarbonateConfig g_carb_config = {0};
static int g_carb_config_initialized = 0;

BiogeoParams* rive_get_params(void) {
    return &g_biogeo_params;
}

/**
 * Helper: Trim whitespace in place
 */
static char *biogeo_trim(char *s) {
    if (!s) return s;
    while (*s && isspace((unsigned char)*s)) s++;
    if (*s == '\0') return s;
    char *end = s + strlen(s) - 1;
    while (end > s && isspace((unsigned char)*end)) *end-- = '\0';
    return s;
}

/**
 * Load biogeochemistry parameters from a config file
 */
int LoadBiogeoParams(const char *path) {
    /* Set defaults first */
    g_biogeo_params.water_temp = 25.0;
    g_biogeo_params.ws = 1e-3;
    g_biogeo_params.I0 = 200.0;
    g_biogeo_params.kd1 = 0.1;
    g_biogeo_params.kd2_spm = 0.01;
    g_biogeo_params.kd2_phy1 = 0.005;
    g_biogeo_params.kd2_phy2 = 0.005;
    g_biogeo_params.alpha1 = 0.02;
    g_biogeo_params.alpha2 = 0.015;
    g_biogeo_params.pbmax1 = 2.0;
    g_biogeo_params.pbmax2 = 1.5;
    g_biogeo_params.kexc1 = 0.1;
    g_biogeo_params.kexc2 = 0.1;
    g_biogeo_params.kgrowth1 = 0.1;
    g_biogeo_params.kgrowth2 = 0.1;
    g_biogeo_params.kmaint1 = 0.02;
    g_biogeo_params.kmaint2 = 0.015;
    g_biogeo_params.kmort1 = 0.05;
    g_biogeo_params.kmort2 = 0.04;
    g_biogeo_params.kdsi1 = 5.0;
    g_biogeo_params.kn1 = 2.0;
    g_biogeo_params.kpo41 = 0.5;
    g_biogeo_params.kn2 = 3.0;
    g_biogeo_params.kpo42 = 0.7;
    g_biogeo_params.kox = 0.1;
    g_biogeo_params.kdenit = 0.05;
    g_biogeo_params.knit = 0.1;
    g_biogeo_params.ktox = 50.0;
    g_biogeo_params.ko2 = 10.0;
    g_biogeo_params.ko2_nit = 5.0;
    g_biogeo_params.kno3 = 5.0;
    g_biogeo_params.knh4 = 2.0;
    g_biogeo_params.kino2 = 5.0;
    g_biogeo_params.redn = 0.151;
    g_biogeo_params.redp = 0.01;
    g_biogeo_params.redsi = 0.1;
    g_biogeo_params.pco2_atm = 420.0;
    
    /* Enhanced pCO2 parameters */
    g_biogeo_params.wind_speed = 3.5;
    g_biogeo_params.wind_coeff = 0.251;
    g_biogeo_params.schmidt_exp = -0.5;
    g_biogeo_params.current_k_factor = 0.25;
    g_biogeo_params.benthic_resp_20C = 60.0;
    g_biogeo_params.benthic_Q10 = 2.0;
    
    /* RIVE organic matter */
    g_biogeo_params.khydr1 = 0.75;
    g_biogeo_params.khydr2 = 0.25;
    g_biogeo_params.khydr3 = 0.01;
    g_biogeo_params.frac_hd1 = 0.15;
    g_biogeo_params.frac_hd2 = 0.35;
    g_biogeo_params.frac_hp1 = 0.10;
    g_biogeo_params.frac_hp2 = 0.30;
    
    /* RIVE bacteria */
    g_biogeo_params.bag_bmax20 = 0.6;
    g_biogeo_params.bap_bmax20 = 0.16;
    g_biogeo_params.bag_kdb20 = 0.05;
    g_biogeo_params.bap_kdb20 = 0.02;
    g_biogeo_params.bac_ks = 0.1;
    g_biogeo_params.bac_yield = 0.25;
    g_biogeo_params.bag_topt = 22.0;
    g_biogeo_params.bap_topt = 20.0;
    g_biogeo_params.bag_dti = 12.0;
    g_biogeo_params.bap_dti = 17.0;
    g_biogeo_params.bag_vs = 0.02;
    
    /* Phosphorus adsorption (Langmuir) */
    g_biogeo_params.kpads = 0.05;   /* Equilibrium constant [L/µmol] - adjusted for delta */
    g_biogeo_params.pac = 0.005;    /* P capacity [µmol P / mg SPM] */
    
    /* Sediment-TOC coupling (From C-RIVE) */
    g_biogeo_params.poc_spm_ratio = 0.02;  /* POC/SPM ratio: 2% (typical for Mekong) */
    
    /* Benthic */
    g_biogeo_params.zf_init = 0.001;
    g_biogeo_params.benthic_porosity = 0.88;
    g_biogeo_params.benthic_density = 2300000.0;
    
    /* GHG */
    g_biogeo_params.N2O_atm_ppb = CGEM_N2O_ATM_DEFAULT;
    g_biogeo_params.N2O_yield_nit = 0.004;    /* N2O yield from nitrification (0.3-1%) */
    g_biogeo_params.N2O_yield_denit = 0.01;   /* N2O yield from denitrification (0.5-2%) */
    g_biogeo_params.CH4_atm_ppb = CGEM_CH4_ATM_DEFAULT;
    g_biogeo_params.CH4_prod_rate = 0.01;
    g_biogeo_params.CH4_ox_rate = 0.5;
    g_biogeo_params.CH4_ks_o2 = 10.0;
    g_biogeo_params.CH4_ebul_thresh = 50.0;
    g_biogeo_params.k_nit1 = 0.3;
    g_biogeo_params.k_nit2 = 0.6;
    g_biogeo_params.ks_no2 = 1.0;
    
    /* Benthic GHG fluxes (for passive mode) */
    g_biogeo_params.benthic_CH4_flux = 100.0;  /* µmol/m²/day - typical for tropical */
    g_biogeo_params.benthic_N2O_flux = 5.0;    /* nmol/m²/day */
    
    /* Solver */
    g_biogeo_params.use_rk4 = 1;
    
    /* ===========================================================================
     * SAFETY MODE FLAGS - Default to SAFE (passive) mode
     * This prevents mass balance errors in data-limited applications
     * ===========================================================================*/
    g_biogeo_params.ghg_passive_mode = 1;     /* Default: Safe diagnostic mode */
    g_biogeo_params.simplified_mode = 0;      /* Default: Full RIVE mode */
    g_biogeo_params.skip_bacteria = 0;
    g_biogeo_params.skip_multipool_oc = 0;
    g_biogeo_params.skip_ghg_dynamics = 0;
    g_biogeo_params.skip_p_adsorption = 0;
    g_biogeo_params.skip_carbonate_reactions = 0;  /* Default: full carbonate chemistry */
    
    /* ===========================================================================
     * TROPICAL ESTUARY PHYSICS (Critical for Mekong Delta)
     * Reference: Audit recommendations for "Academic Solidity"
     * ===========================================================================*/
    
    /* Salinity stress on freshwater phytoplankton */
    g_biogeo_params.sal_stress_thresh = 5.0;   /* PSU - onset of stress */
    g_biogeo_params.sal_stress_coef = 0.5;     /* Mortality increase rate */
    g_biogeo_params.sal_marine_opt = 30.0;     /* Optimal salinity for marine spp */
    
    /* Flocculation - salinity-dependent settling velocity */
    g_biogeo_params.floc_sal_scale = 2.0;      /* PSU - hyperbolic scale */
    g_biogeo_params.floc_factor_max = 10.0;    /* Maximum ws enhancement */
    g_biogeo_params.enable_flocculation = 1;   /* ON by default */
    
    /* Tidal pumping for benthic fluxes */
    g_biogeo_params.tidal_pump_coef = 1.5;     /* Enhancement factor */
    g_biogeo_params.tidal_pump_vel_ref = 0.5;  /* Reference velocity [m/s] */
    
    /* Simplified mode Arrhenius coefficients */
    g_biogeo_params.theta_ox = 1.047;          /* Oxidation temp coefficient */
    g_biogeo_params.theta_nit = 1.047;         /* Nitrification temp coefficient */
    g_biogeo_params.theta_denit = 1.047;       /* Denitrification temp coefficient */
    
    /* Initialize GHG config */
    if (!g_ghg_config_initialized) {
        crive_ghg_init_config(&g_ghg_config);
        g_ghg_config_initialized = 1;
    }
    
    /* Initialize carbonate config */
    if (!g_carb_config_initialized) {
        crive_carbonate_init_config(&g_carb_config);
        g_carb_config_initialized = 1;
    }
    
    if (!path) {
        g_biogeo_params.loaded = 1;
        return 0;
    }
    
    FILE *fp = fopen(path, "r");
    if (!fp) {
        g_biogeo_params.loaded = 1;
        return -1;
    }
    
    printf("Loading biogeochemistry parameters from: %s\n", path);
    
    char line[256];
    while (fgets(line, sizeof(line), fp)) {
        char *text = biogeo_trim(line);
        if (*text == '\0' || *text == '#') continue;
        
        char *eq = strchr(text, '=');
        if (!eq) continue;
        
        *eq = '\0';
        char *key = biogeo_trim(text);
        char *val = biogeo_trim(eq + 1);
        double v = strtod(val, NULL);
        
        /* Match parameter names */
        if (strcmp(key, "water_temp") == 0) g_biogeo_params.water_temp = v;
        else if (strcmp(key, "ws") == 0) g_biogeo_params.ws = v;
        else if (strcmp(key, "I0") == 0) g_biogeo_params.I0 = v;
        else if (strcmp(key, "kd1") == 0) g_biogeo_params.kd1 = v;
        else if (strcmp(key, "kd2_spm") == 0) g_biogeo_params.kd2_spm = v;
        else if (strcmp(key, "kd2_phy1") == 0) g_biogeo_params.kd2_phy1 = v;
        else if (strcmp(key, "kd2_phy2") == 0) g_biogeo_params.kd2_phy2 = v;
        else if (strcmp(key, "alpha1") == 0) g_biogeo_params.alpha1 = v;
        else if (strcmp(key, "alpha2") == 0) g_biogeo_params.alpha2 = v;
        else if (strcmp(key, "pbmax1") == 0) g_biogeo_params.pbmax1 = v;
        else if (strcmp(key, "pbmax2") == 0) g_biogeo_params.pbmax2 = v;
        else if (strcmp(key, "kexc1") == 0) g_biogeo_params.kexc1 = v;
        else if (strcmp(key, "kexc2") == 0) g_biogeo_params.kexc2 = v;
        else if (strcmp(key, "kgrowth1") == 0) g_biogeo_params.kgrowth1 = v;
        else if (strcmp(key, "kgrowth2") == 0) g_biogeo_params.kgrowth2 = v;
        else if (strcmp(key, "kmaint1") == 0) g_biogeo_params.kmaint1 = v;
        else if (strcmp(key, "kmaint2") == 0) g_biogeo_params.kmaint2 = v;
        else if (strcmp(key, "kmort1") == 0) g_biogeo_params.kmort1 = v;
        else if (strcmp(key, "kmort2") == 0) g_biogeo_params.kmort2 = v;
        else if (strcmp(key, "kdsi1") == 0) g_biogeo_params.kdsi1 = v;
        else if (strcmp(key, "kn1") == 0) g_biogeo_params.kn1 = v;
        else if (strcmp(key, "kpo41") == 0) g_biogeo_params.kpo41 = v;
        else if (strcmp(key, "kn2") == 0) g_biogeo_params.kn2 = v;
        else if (strcmp(key, "kpo42") == 0) g_biogeo_params.kpo42 = v;
        else if (strcmp(key, "kox") == 0) g_biogeo_params.kox = v;
        else if (strcmp(key, "kdenit") == 0) g_biogeo_params.kdenit = v;
        else if (strcmp(key, "knit") == 0) g_biogeo_params.knit = v;
        else if (strcmp(key, "ktox") == 0) g_biogeo_params.ktox = v;
        else if (strcmp(key, "ko2") == 0) g_biogeo_params.ko2 = v;
        else if (strcmp(key, "ko2_nit") == 0) g_biogeo_params.ko2_nit = v;
        else if (strcmp(key, "kno3") == 0) g_biogeo_params.kno3 = v;
        else if (strcmp(key, "knh4") == 0) g_biogeo_params.knh4 = v;
        else if (strcmp(key, "kino2") == 0) g_biogeo_params.kino2 = v;
        else if (strcmp(key, "redn") == 0) g_biogeo_params.redn = v;
        else if (strcmp(key, "redp") == 0) g_biogeo_params.redp = v;
        else if (strcmp(key, "redsi") == 0) g_biogeo_params.redsi = v;
        else if (strcmp(key, "pco2_atm") == 0) g_biogeo_params.pco2_atm = v;
        else if (strcmp(key, "wind_speed") == 0) g_biogeo_params.wind_speed = v;
        else if (strcmp(key, "wind_coeff") == 0) g_biogeo_params.wind_coeff = v;
        else if (strcmp(key, "schmidt_exp") == 0) g_biogeo_params.schmidt_exp = v;
        else if (strcmp(key, "current_k_factor") == 0) g_biogeo_params.current_k_factor = v;
        else if (strcmp(key, "benthic_resp_20C") == 0) g_biogeo_params.benthic_resp_20C = v;
        else if (strcmp(key, "benthic_Q10") == 0) g_biogeo_params.benthic_Q10 = v;
        else if (strcmp(key, "khydr1") == 0) g_biogeo_params.khydr1 = v;
        else if (strcmp(key, "khydr2") == 0) g_biogeo_params.khydr2 = v;
        else if (strcmp(key, "khydr3") == 0) g_biogeo_params.khydr3 = v;
        else if (strcmp(key, "frac_hd1") == 0) g_biogeo_params.frac_hd1 = v;
        else if (strcmp(key, "frac_hd2") == 0) g_biogeo_params.frac_hd2 = v;
        else if (strcmp(key, "frac_hp1") == 0) g_biogeo_params.frac_hp1 = v;
        else if (strcmp(key, "frac_hp2") == 0) g_biogeo_params.frac_hp2 = v;
        else if (strcmp(key, "bag_bmax20") == 0) g_biogeo_params.bag_bmax20 = v;
        else if (strcmp(key, "bap_bmax20") == 0) g_biogeo_params.bap_bmax20 = v;
        else if (strcmp(key, "bag_kdb20") == 0) g_biogeo_params.bag_kdb20 = v;
        else if (strcmp(key, "bap_kdb20") == 0) g_biogeo_params.bap_kdb20 = v;
        else if (strcmp(key, "bac_ks") == 0) g_biogeo_params.bac_ks = v;
        else if (strcmp(key, "bac_yield") == 0) g_biogeo_params.bac_yield = v;
        else if (strcmp(key, "bag_topt") == 0) g_biogeo_params.bag_topt = v;
        else if (strcmp(key, "bap_topt") == 0) g_biogeo_params.bap_topt = v;
        else if (strcmp(key, "bag_dti") == 0) g_biogeo_params.bag_dti = v;
        else if (strcmp(key, "bap_dti") == 0) g_biogeo_params.bap_dti = v;
        else if (strcmp(key, "bag_vs") == 0) g_biogeo_params.bag_vs = v;
        /* Phosphorus adsorption & Sediment-TOC coupling */
        else if (strcmp(key, "kpads") == 0) g_biogeo_params.kpads = v;
        else if (strcmp(key, "pac") == 0) g_biogeo_params.pac = v;
        else if (strcmp(key, "poc_spm_ratio") == 0) g_biogeo_params.poc_spm_ratio = v;
        else if (strcmp(key, "zf_init") == 0) g_biogeo_params.zf_init = v;
        else if (strcmp(key, "benthic_porosity") == 0) g_biogeo_params.benthic_porosity = v;
        else if (strcmp(key, "benthic_density") == 0) g_biogeo_params.benthic_density = v;
        /* GHG parameters */
        else if (strcmp(key, "N2O_yield_nit") == 0) g_biogeo_params.N2O_yield_nit = v;
        else if (strcmp(key, "N2O_yield_denit") == 0) g_biogeo_params.N2O_yield_denit = v;
        else if (strcmp(key, "benthic_CH4_flux") == 0) g_biogeo_params.benthic_CH4_flux = v;
        else if (strcmp(key, "benthic_N2O_flux") == 0) g_biogeo_params.benthic_N2O_flux = v;
        /* Safety mode flags */
        else if (strcmp(key, "ghg_passive_mode") == 0) g_biogeo_params.ghg_passive_mode = (int)v;
        else if (strcmp(key, "simplified_mode") == 0) g_biogeo_params.simplified_mode = (int)v;
        else if (strcmp(key, "skip_bacteria") == 0) g_biogeo_params.skip_bacteria = (int)v;
        else if (strcmp(key, "skip_multipool_oc") == 0) g_biogeo_params.skip_multipool_oc = (int)v;
        else if (strcmp(key, "skip_ghg_dynamics") == 0) g_biogeo_params.skip_ghg_dynamics = (int)v;
        else if (strcmp(key, "skip_p_adsorption") == 0) g_biogeo_params.skip_p_adsorption = (int)v;
        else if (strcmp(key, "skip_carbonate_reactions") == 0) g_biogeo_params.skip_carbonate_reactions = (int)v;
        /* Tropical physics parameters */
        else if (strcmp(key, "sal_stress_thresh") == 0) g_biogeo_params.sal_stress_thresh = v;
        else if (strcmp(key, "sal_stress_coef") == 0) g_biogeo_params.sal_stress_coef = v;
        else if (strcmp(key, "sal_marine_opt") == 0) g_biogeo_params.sal_marine_opt = v;
        else if (strcmp(key, "floc_sal_scale") == 0) g_biogeo_params.floc_sal_scale = v;
        else if (strcmp(key, "floc_factor_max") == 0) g_biogeo_params.floc_factor_max = v;
        else if (strcmp(key, "enable_flocculation") == 0) g_biogeo_params.enable_flocculation = (int)v;
        else if (strcmp(key, "tidal_pump_coef") == 0) g_biogeo_params.tidal_pump_coef = v;
        else if (strcmp(key, "tidal_pump_vel_ref") == 0) g_biogeo_params.tidal_pump_vel_ref = v;
        else if (strcmp(key, "theta_ox") == 0) g_biogeo_params.theta_ox = v;
        else if (strcmp(key, "theta_nit") == 0) g_biogeo_params.theta_nit = v;
        else if (strcmp(key, "theta_denit") == 0) g_biogeo_params.theta_denit = v;
    }
    
    fclose(fp);
    g_biogeo_params.loaded = 1;
    printf("Biogeochemistry parameters loaded.\n");
    return 0;
}

/**
 * Initialize biogeochemical parameters for a branch
 */
void InitializeBiogeoParameters(Branch *branch) {
    if (!branch) return;
    
    /* Ensure global params are initialized */
    if (!g_biogeo_params.loaded) {
        LoadBiogeoParams(NULL);
    }
    
    BiogeoParams *p = &g_biogeo_params;
    
    /* Water and environmental */
    branch->water_temp = p->water_temp;
    branch->ws = p->ws;
    
    /* Light */
    branch->I0 = p->I0;
    branch->kd1 = p->kd1;
    branch->kd2_spm = p->kd2_spm;
    branch->kd2_phy1 = p->kd2_phy1;
    branch->kd2_phy2 = p->kd2_phy2;
    
    /* Phytoplankton */
    branch->alpha1 = p->alpha1;
    branch->alpha2 = p->alpha2;
    branch->pbmax1 = p->pbmax1;
    branch->pbmax2 = p->pbmax2;
    branch->kexc1 = p->kexc1;
    branch->kexc2 = p->kexc2;
    branch->kgrowth1 = p->kgrowth1;
    branch->kgrowth2 = p->kgrowth2;
    branch->kmaint1 = p->kmaint1;
    branch->kmaint2 = p->kmaint2;
    branch->kmort1 = p->kmort1;
    branch->kmort2 = p->kmort2;
    
    /* Nutrients */
    branch->kdsi1 = p->kdsi1;
    branch->kn1 = p->kn1;
    branch->kpo41 = p->kpo41;
    branch->kn2 = p->kn2;
    branch->kpo42 = p->kpo42;
    
    /* Decomposition */
    branch->kox = p->kox;
    branch->kdenit = p->kdenit;
    branch->knit = p->knit;
    branch->ktox = p->ktox;
    branch->ko2 = p->ko2;
    branch->ko2_nit = p->ko2_nit;
    branch->kno3 = p->kno3;
    branch->knh4 = p->knh4;
    branch->kino2 = p->kino2;
    
    /* Stoichiometry */
    branch->redn = p->redn;
    branch->redp = p->redp;
    branch->redsi = p->redsi;
    
    /* Gas exchange */
    branch->pco2_atm = p->pco2_atm;
    
    /* RIVE organic matter */
    branch->khydr1 = p->khydr1;
    branch->khydr2 = p->khydr2;
    branch->khydr3 = p->khydr3;
    branch->frac_hd1 = p->frac_hd1;
    branch->frac_hd2 = p->frac_hd2;
    branch->frac_hp1 = p->frac_hp1;
    branch->frac_hp2 = p->frac_hp2;
    
    /* RIVE bacteria */
    branch->bag_bmax20 = p->bag_bmax20;
    branch->bap_bmax20 = p->bap_bmax20;
    branch->bag_kdb20 = p->bag_kdb20;
    branch->bap_kdb20 = p->bap_kdb20;
    branch->bac_ks = p->bac_ks;
    branch->bac_yield = p->bac_yield;
    branch->bag_topt = p->bag_topt;
    branch->bap_topt = p->bap_topt;
    branch->bag_dti = p->bag_dti;
    branch->bap_dti = p->bap_dti;
    branch->bag_vs = p->bag_vs;
    
    /* Phosphorus */
    branch->kpads = p->kpads;
    branch->pac = p->pac;
    
    /* Benthic */
    branch->zf_init = p->zf_init;
    branch->benthic_porosity = p->benthic_porosity;
    branch->benthic_density = p->benthic_density;
}
