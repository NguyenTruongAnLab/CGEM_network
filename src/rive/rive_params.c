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

/* ===========================================================================
 * REGIONAL DEFAULT PARAMETER PRESETS
 * 
 * Pre-calibrated parameter sets for common tropical delta systems.
 * These provide reasonable starting points for data-sparse applications.
 * 
 * Usage: Call SetRegionalDefaults("Mekong") AFTER LoadBiogeoParams()
 *        Only modifies parameters that weren't explicitly set in config file.
 *
 * References:
 * - Mekong: Nguyen et al. (2008), MRC State of Basin (2018), DRAGON (2015)
 * - Ganges: Mukhopadhyay et al. (2006), Jahan & Strezov (2017)
 * - Niger: Marty et al. (2002), Borges & Frankignoulle (2002)
 * ===========================================================================*/

typedef struct {
    const char *name;
    const char *description;
    double water_temp;      /* Mean annual temperature [°C] */
    double kox;             /* TOC oxidation rate [/day] */
    double knit;            /* Nitrification rate [/day] */
    double kdenit;          /* Denitrification rate [/day] */
    double benthic_resp_20C;/* Benthic O2 demand [mmol/m²/day] */
    double benthic_NH4_flux;/* Benthic NH4 flux [mmol/m²/day] */
    double N2O_yield_nit;   /* N2O from nitrification */
    double N2O_yield_denit; /* N2O from denitrification */
    double CH4_prod_rate;   /* CH4 production rate [/day] */
    double poc_spm_ratio;   /* POC/SPM ratio */
    double sal_stress_thresh;/* Salinity stress threshold [PSU] */
} RegionalPreset;

static const RegionalPreset REGIONAL_PRESETS[] = {
    /* Mekong Delta: Highly productive, warm, high DOC from rice paddies
     * Validated against March 2025 field data */
    {
        .name = "Mekong",
        .description = "Mekong Delta, Vietnam - Monsoon tropical with rice agriculture",
        .water_temp = 28.0,
        .kox = 0.005,           /* Low - refractory DOC from paddies */
        .knit = 0.12,           /* Moderate nitrification */
        .kdenit = 0.04,         /* Active denitrification in paddies */
        .benthic_resp_20C = 70.0,/* High organic sediments */
        .benthic_NH4_flux = 2.5,/* NH4 from anaerobic sediments */
        .N2O_yield_nit = 0.004, /* 0.4% yield */
        .N2O_yield_denit = 0.01,/* 1% yield */
        .CH4_prod_rate = 0.02,  /* Active methanogenesis */
        .poc_spm_ratio = 0.02,  /* 2% POC in SPM */
        .sal_stress_thresh = 5.0
    },
    /* Red River Delta: Similar to Mekong but cooler, less aquaculture */
    {
        .name = "RedRiver",
        .description = "Red River Delta, Vietnam - Subtropical monsoon",
        .water_temp = 24.0,
        .kox = 0.008,           /* Slightly higher degradation (less refractory) */
        .knit = 0.10,
        .kdenit = 0.05,
        .benthic_resp_20C = 60.0,
        .benthic_NH4_flux = 2.0,
        .N2O_yield_nit = 0.005,
        .N2O_yield_denit = 0.012,
        .CH4_prod_rate = 0.015,
        .poc_spm_ratio = 0.025, /* Higher POC due to less dilution */
        .sal_stress_thresh = 4.0
    },
    /* Ganges-Brahmaputra: Very high sediment load, low DOC, high nutrients */
    {
        .name = "Ganges",
        .description = "Ganges-Brahmaputra Delta, Bangladesh - High sediment monsoon",
        .water_temp = 27.0,
        .kox = 0.01,            /* Faster degradation (more labile) */
        .knit = 0.15,           /* High N inputs from agriculture */
        .kdenit = 0.06,         /* High denitrification */
        .benthic_resp_20C = 80.0,/* Very organic sediments */
        .benthic_NH4_flux = 3.5,/* High NH4 from intense farming */
        .N2O_yield_nit = 0.006,
        .N2O_yield_denit = 0.015,
        .CH4_prod_rate = 0.025, /* High rice paddies */
        .poc_spm_ratio = 0.01,  /* Low POC due to high mineral SPM */
        .sal_stress_thresh = 5.0
    },
    /* Niger Delta: Tropical equatorial, mangrove-dominated, oil pollution */
    {
        .name = "Niger",
        .description = "Niger Delta, Nigeria - Equatorial mangrove system",
        .water_temp = 28.5,
        .kox = 0.003,           /* Very slow - refractory mangrove DOC */
        .knit = 0.08,           /* Lower N inputs */
        .kdenit = 0.03,
        .benthic_resp_20C = 50.0,/* Lower due to sandy sediments */
        .benthic_NH4_flux = 1.5,
        .N2O_yield_nit = 0.003,
        .N2O_yield_denit = 0.008,
        .CH4_prod_rate = 0.03,  /* High from mangrove sediments */
        .poc_spm_ratio = 0.03,  /* Higher POC from mangroves */
        .sal_stress_thresh = 3.0 /* Mangrove species more tolerant */
    },
    /* Irrawaddy Delta: Myanmar - Similar to Mekong but less developed */
    {
        .name = "Irrawaddy",
        .description = "Irrawaddy Delta, Myanmar - Less impacted monsoon system",
        .water_temp = 27.0,
        .kox = 0.006,
        .knit = 0.10,
        .kdenit = 0.04,
        .benthic_resp_20C = 55.0,
        .benthic_NH4_flux = 2.0,
        .N2O_yield_nit = 0.004,
        .N2O_yield_denit = 0.01,
        .CH4_prod_rate = 0.015,
        .poc_spm_ratio = 0.02,
        .sal_stress_thresh = 5.0
    },
    /* Saigon-Dong Nai: Urban-impacted tropical estuary */
    {
        .name = "SaigonDongNai",
        .description = "Saigon-Dong Nai basin, Vietnam - Urban tropical estuary",
        .water_temp = 28.5,
        .kox = 0.015,           /* Higher - urban organic matter more labile */
        .knit = 0.18,           /* High - lots of NH4 from sewage */
        .kdenit = 0.08,         /* High - hypoxic conditions */
        .benthic_resp_20C = 100.0,/* Very high - organic pollution */
        .benthic_NH4_flux = 5.0,/* Very high from sewage sludge */
        .N2O_yield_nit = 0.008, /* High - WWTP contribution */
        .N2O_yield_denit = 0.02,/* High - low O2 conditions */
        .CH4_prod_rate = 0.04,  /* High - anaerobic sediments */
        .poc_spm_ratio = 0.015,
        .sal_stress_thresh = 4.0
    },
    /* Mediterranean type - for comparison */
    {
        .name = "Mediterranean",
        .description = "Mediterranean climate estuary - Dry summer, wet winter",
        .water_temp = 18.0,     /* Cooler */
        .kox = 0.02,            /* Faster at lower temp reference */
        .knit = 0.15,
        .kdenit = 0.03,
        .benthic_resp_20C = 40.0,/* Lower in cooler water */
        .benthic_NH4_flux = 1.0,
        .N2O_yield_nit = 0.005,
        .N2O_yield_denit = 0.01,
        .CH4_prod_rate = 0.005, /* Low - cooler, less anaerobic */
        .poc_spm_ratio = 0.02,
        .sal_stress_thresh = 8.0 /* Higher tolerance */
    },
    /* Sentinel - marks end of array */
    { .name = NULL }
};

/**
 * Apply regional default parameters
 * 
 * @param region_name Name of the regional preset (case-insensitive)
 * @return 0 on success, -1 if region not found
 * 
 * Note: This should be called AFTER loading biogeo_params.txt to apply
 * regional defaults only for parameters not explicitly specified.
 */
int SetRegionalDefaults(const char *region_name) {
    if (!region_name) return -1;
    
    /* Find matching preset */
    const RegionalPreset *preset = NULL;
    for (int i = 0; REGIONAL_PRESETS[i].name != NULL; ++i) {
        if (strcasecmp(REGIONAL_PRESETS[i].name, region_name) == 0) {
            preset = &REGIONAL_PRESETS[i];
            break;
        }
    }
    
    if (!preset) {
        fprintf(stderr, "Warning: Regional preset '%s' not found. Available:\n", region_name);
        for (int i = 0; REGIONAL_PRESETS[i].name != NULL; ++i) {
            fprintf(stderr, "  - %s: %s\n", REGIONAL_PRESETS[i].name, REGIONAL_PRESETS[i].description);
        }
        return -1;
    }
    
    printf("Applying regional defaults: %s\n", preset->description);
    
    /* Apply preset values (only if not already loaded from file) */
    /* Note: LoadBiogeoParams sets defaults first, so we override here */
    g_biogeo_params.water_temp = preset->water_temp;
    g_biogeo_params.kox = preset->kox;
    g_biogeo_params.knit = preset->knit;
    g_biogeo_params.kdenit = preset->kdenit;
    g_biogeo_params.benthic_resp_20C = preset->benthic_resp_20C;
    g_biogeo_params.N2O_yield_nit = preset->N2O_yield_nit;
    g_biogeo_params.N2O_yield_denit = preset->N2O_yield_denit;
    g_biogeo_params.CH4_prod_rate = preset->CH4_prod_rate;
    g_biogeo_params.poc_spm_ratio = preset->poc_spm_ratio;
    g_biogeo_params.sal_stress_thresh = preset->sal_stress_thresh;
    
    return 0;
}

/**
 * List available regional presets
 */
void ListRegionalPresets(void) {
    printf("\nAvailable regional parameter presets:\n");
    printf("--------------------------------------\n");
    for (int i = 0; REGIONAL_PRESETS[i].name != NULL; ++i) {
        printf("  %-15s : %s\n", REGIONAL_PRESETS[i].name, REGIONAL_PRESETS[i].description);
        printf("                    T=%.0f°C, kox=%.3f/d, knit=%.2f/d\n",
               REGIONAL_PRESETS[i].water_temp, REGIONAL_PRESETS[i].kox, REGIONAL_PRESETS[i].knit);
    }
    printf("\nUsage: Set 'RegionalPreset = Mekong' in case_config.txt\n\n");
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
