/**
 * @file biogeo.c
 * @brief C-GEM Network Biogeochemical Module (RIVE-Enhanced)
 * 
 * Implements water quality reactions including:
 * - RIVE multi-pool organic matter (HD1/HD2/HD3, HP1/HP2/HP3)
 * - RIVE heterotrophic bacteria (BAG/BAP)
 * - Phytoplankton growth with nutrient limitation
 * - Nutrient cycling (N, P, Si)
 * - Oxygen dynamics with benthic-pelagic coupling
 * - Carbon chemistry (pH, pCO2, CO2 flux)
 * - Phosphorus adsorption equilibrium on SPM
 * 
 * Reference: Billen et al. (1994), Garnier et al. (2002), Volta et al. (2014)
 */

#include "network.h"
#include "define.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

/* Seconds per day - for unit conversion */
#define SECONDS_PER_DAY 86400.0

/* RIVE constants */
#define RIVE_SQUARE(x) ((x) * (x))

/* Global biogeo parameters (loaded from config file) */
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
    
    /* Decomposition - NOTE: All rates stored as [1/day], converted to [1/s] in calculations */
    double kox, kdenit, knit;
    double ktox, ko2, ko2_nit, kno3, knh4, kino2;
    
    /* Stoichiometry */
    double redn, redp, redsi;
    
    /* Gas exchange */
    double pco2_atm;
    
    /* ========== ENHANCED pCO2 PARAMETERS ========== */
    
    /* Wind-driven gas exchange (Wanninkhof 2014) */
    double wind_speed;          /* [m/s] - mean wind at 10m height */
    double wind_coeff;          /* Wanninkhof coefficient (0.251 for long-term avg) */
    double schmidt_exp;         /* Schmidt number exponent (-0.5 typical) */
    double current_k_factor;    /* Current contribution factor (0.1-0.5) */
    
    /* Benthic respiration (major CO2 source in shallow estuaries) */
    /* NOTE: Stored as [mmol C/m²/day], converted in calculations */
    double benthic_resp_20C;    /* [mmol C/m²/day] at 20°C - sediment respiration */
    double benthic_Q10;         /* Temperature coefficient (typically 2.0) */
    
    /* ========== RIVE ORGANIC MATTER PARAMETERS ========== */
    double khydr1;              /* Labile hydrolysis rate [1/day] (0.75) */
    double khydr2;              /* Semi-labile hydrolysis rate [1/day] (0.25) */
    double khydr3;              /* Refractory hydrolysis rate [1/day] (0.01) */
    double frac_hd1;            /* Fraction of DOC as HD1 (0.15) */
    double frac_hd2;            /* Fraction of DOC as HD2 (0.35) */
    double frac_hp1;            /* Fraction of POC as HP1 (0.10) */
    double frac_hp2;            /* Fraction of POC as HP2 (0.30) */
    
    /* ========== RIVE BACTERIA PARAMETERS ========== */
    double bag_bmax20;          /* BAG max growth rate at 20°C [1/h] (0.6) */
    double bap_bmax20;          /* BAP max growth rate at 20°C [1/h] (0.16) */
    double bag_kdb20;           /* BAG mortality rate at 20°C [1/h] (0.05) */
    double bap_kdb20;           /* BAP mortality rate at 20°C [1/h] (0.02) */
    double bac_ks;              /* Bacteria half-saturation for DSS [mg C/L] (0.1) */
    double bac_yield;           /* Bacteria growth yield [-] (0.25) */
    double bag_topt;            /* BAG optimal temperature [°C] (22) */
    double bap_topt;            /* BAP optimal temperature [°C] (20) */
    double bag_dti;             /* BAG temperature spread [°C] (12) */
    double bap_dti;             /* BAP temperature spread [°C] (17) */
    double bag_vs;              /* BAG settling velocity [m/h] (0.02) */
    
    /* ========== RIVE PHOSPHORUS ADSORPTION ========== */
    double kpads;               /* P adsorption equilibrium constant (3.43) */
    double pac;                 /* P adsorption capacity [µmol P/mg SPM] (0.37) */
    
    /* ========== RIVE BENTHIC PARAMETERS ========== */
    double zf_init;             /* Initial fluid sediment depth [m] (0.001) */
    double benthic_porosity;    /* Sediment porosity [-] (0.88) */
    double benthic_density;     /* Sediment density [g/m³] (2300000) */
    
} BiogeoParams;

static BiogeoParams g_biogeo_params = {0};

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
 * 
 * File format (key = value pairs):
 *   # Comment
 *   water_temp = 25.0
 *   I0 = 200.0
 *   ...
 * 
 * @param path Path to biogeo_params.txt
 * @return 0 on success, -1 on file not found (uses defaults)
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
    
    /* Enhanced pCO2 parameters - Mekong Delta defaults */
    g_biogeo_params.wind_speed = 3.5;           /* Mean wind [m/s] - Mekong typical */
    g_biogeo_params.wind_coeff = 0.251;         /* Wanninkhof (2014) for long-term avg */
    g_biogeo_params.schmidt_exp = -0.5;         /* Standard Schmidt exponent */
    g_biogeo_params.current_k_factor = 0.25;    /* Current contribution (25%) */
    g_biogeo_params.benthic_resp_20C = 60.0;    /* [mmol C/m²/day] - tropical deltaic mud */
    g_biogeo_params.benthic_Q10 = 2.0;          /* Temperature coefficient */
    
    /* RIVE organic matter defaults (from RIVE organic.cpp) */
    g_biogeo_params.khydr1 = 0.75;              /* Labile hydrolysis [1/day] */
    g_biogeo_params.khydr2 = 0.25;              /* Semi-labile hydrolysis [1/day] */
    g_biogeo_params.khydr3 = 0.01;              /* Refractory hydrolysis [1/day] */
    g_biogeo_params.frac_hd1 = 0.15;            /* Fraction DOC as labile */
    g_biogeo_params.frac_hd2 = 0.35;            /* Fraction DOC as semi-labile */
    g_biogeo_params.frac_hp1 = 0.10;            /* Fraction POC as labile */
    g_biogeo_params.frac_hp2 = 0.30;            /* Fraction POC as semi-labile */
    
    /* RIVE bacteria defaults (from RIVE bacteria.cpp) */
    g_biogeo_params.bag_bmax20 = 0.6;           /* BAG max growth [1/h] */
    g_biogeo_params.bap_bmax20 = 0.16;          /* BAP max growth [1/h] */
    g_biogeo_params.bag_kdb20 = 0.05;           /* BAG mortality [1/h] */
    g_biogeo_params.bap_kdb20 = 0.02;           /* BAP mortality [1/h] */
    g_biogeo_params.bac_ks = 0.1;               /* DSS half-saturation [mg C/L] */
    g_biogeo_params.bac_yield = 0.25;           /* Growth yield [-] */
    g_biogeo_params.bag_topt = 22.0;            /* BAG T optimum [°C] */
    g_biogeo_params.bap_topt = 20.0;            /* BAP T optimum [°C] */
    g_biogeo_params.bag_dti = 12.0;             /* BAG T spread [°C] */
    g_biogeo_params.bap_dti = 17.0;             /* BAP T spread [°C] */
    g_biogeo_params.bag_vs = 0.02;              /* BAG settling [m/h] */
    
    /* RIVE phosphorus adsorption defaults (from phosphorus.cpp) */
    g_biogeo_params.kpads = 3.43;               /* P adsorption constant */
    g_biogeo_params.pac = 0.37;                 /* P adsorption capacity */
    
    /* RIVE benthic defaults */
    g_biogeo_params.zf_init = 0.001;            /* Initial fluid sediment [m] */
    g_biogeo_params.benthic_porosity = 0.88;    /* Porosity [-] */
    g_biogeo_params.benthic_density = 2300000.0; /* Density [g/m³] */
    
    if (!path) {
        g_biogeo_params.loaded = 1;
        return 0;
    }
    
    FILE *fp = fopen(path, "r");
    if (!fp) {
        /* File not found - use defaults silently */
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
        /* Enhanced pCO2 parameters */
        else if (strcmp(key, "wind_speed") == 0) g_biogeo_params.wind_speed = v;
        else if (strcmp(key, "wind_coeff") == 0) g_biogeo_params.wind_coeff = v;
        else if (strcmp(key, "schmidt_exp") == 0) g_biogeo_params.schmidt_exp = v;
        else if (strcmp(key, "current_k_factor") == 0) g_biogeo_params.current_k_factor = v;
        else if (strcmp(key, "benthic_resp_20C") == 0) g_biogeo_params.benthic_resp_20C = v;
        else if (strcmp(key, "benthic_Q10") == 0) g_biogeo_params.benthic_Q10 = v;
        /* RIVE organic matter parameters */
        else if (strcmp(key, "khydr1") == 0) g_biogeo_params.khydr1 = v;
        else if (strcmp(key, "khydr2") == 0) g_biogeo_params.khydr2 = v;
        else if (strcmp(key, "khydr3") == 0) g_biogeo_params.khydr3 = v;
        else if (strcmp(key, "frac_hd1") == 0) g_biogeo_params.frac_hd1 = v;
        else if (strcmp(key, "frac_hd2") == 0) g_biogeo_params.frac_hd2 = v;
        else if (strcmp(key, "frac_hp1") == 0) g_biogeo_params.frac_hp1 = v;
        else if (strcmp(key, "frac_hp2") == 0) g_biogeo_params.frac_hp2 = v;
        /* RIVE bacteria parameters */
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
        /* RIVE phosphorus parameters */
        else if (strcmp(key, "kpads") == 0) g_biogeo_params.kpads = v;
        else if (strcmp(key, "pac") == 0) g_biogeo_params.pac = v;
        /* RIVE benthic parameters */
        else if (strcmp(key, "zf_init") == 0) g_biogeo_params.zf_init = v;
        else if (strcmp(key, "benthic_porosity") == 0) g_biogeo_params.benthic_porosity = v;
        else if (strcmp(key, "benthic_density") == 0) g_biogeo_params.benthic_density = v;
    }
    
    fclose(fp);
    g_biogeo_params.loaded = 1;
    printf("Biogeochemistry parameters loaded.\n");
    return 0;
}

/**
 * Load biogeochemistry parameters from file into a specific BiogeoParams struct
 * Used for per-branch parameter loading
 * 
 * @param path Path to biogeo_params.txt
 * @param params Pointer to BiogeoParams struct to populate
 * @return 0 on success, -1 on file not found (params unchanged)
 */
static int LoadBiogeoParamsToStruct(const char *path, BiogeoParams *params) {
    if (!path || !params) return -1;
    
    FILE *fp = fopen(path, "r");
    if (!fp) {
        return -1;
    }
    
    printf("  Loading branch-specific biogeo params from: %s\n", path);
    
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
        
        /* Match parameter names - only override non-zero values */
        if (strcmp(key, "water_temp") == 0) params->water_temp = v;
        else if (strcmp(key, "ws") == 0) params->ws = v;
        else if (strcmp(key, "I0") == 0) params->I0 = v;
        else if (strcmp(key, "kd1") == 0) params->kd1 = v;
        else if (strcmp(key, "kd2_spm") == 0) params->kd2_spm = v;
        else if (strcmp(key, "kd2_phy1") == 0) params->kd2_phy1 = v;
        else if (strcmp(key, "kd2_phy2") == 0) params->kd2_phy2 = v;
        else if (strcmp(key, "alpha1") == 0) params->alpha1 = v;
        else if (strcmp(key, "alpha2") == 0) params->alpha2 = v;
        else if (strcmp(key, "pbmax1") == 0) params->pbmax1 = v;
        else if (strcmp(key, "pbmax2") == 0) params->pbmax2 = v;
        else if (strcmp(key, "kexc1") == 0) params->kexc1 = v;
        else if (strcmp(key, "kexc2") == 0) params->kexc2 = v;
        else if (strcmp(key, "kgrowth1") == 0) params->kgrowth1 = v;
        else if (strcmp(key, "kgrowth2") == 0) params->kgrowth2 = v;
        else if (strcmp(key, "kmaint1") == 0) params->kmaint1 = v;
        else if (strcmp(key, "kmaint2") == 0) params->kmaint2 = v;
        else if (strcmp(key, "kmort1") == 0) params->kmort1 = v;
        else if (strcmp(key, "kmort2") == 0) params->kmort2 = v;
        else if (strcmp(key, "kdsi1") == 0) params->kdsi1 = v;
        else if (strcmp(key, "kn1") == 0) params->kn1 = v;
        else if (strcmp(key, "kpo41") == 0) params->kpo41 = v;
        else if (strcmp(key, "kn2") == 0) params->kn2 = v;
        else if (strcmp(key, "kpo42") == 0) params->kpo42 = v;
        else if (strcmp(key, "kox") == 0) params->kox = v;
        else if (strcmp(key, "kdenit") == 0) params->kdenit = v;
        else if (strcmp(key, "knit") == 0) params->knit = v;
        else if (strcmp(key, "ktox") == 0) params->ktox = v;
        else if (strcmp(key, "ko2") == 0) params->ko2 = v;
        else if (strcmp(key, "ko2_nit") == 0) params->ko2_nit = v;
        else if (strcmp(key, "kno3") == 0) params->kno3 = v;
        else if (strcmp(key, "knh4") == 0) params->knh4 = v;
        else if (strcmp(key, "kino2") == 0) params->kino2 = v;
        else if (strcmp(key, "redn") == 0) params->redn = v;
        else if (strcmp(key, "redp") == 0) params->redp = v;
        else if (strcmp(key, "redsi") == 0) params->redsi = v;
        else if (strcmp(key, "pco2_atm") == 0) params->pco2_atm = v;
        /* Enhanced pCO2 parameters */
        else if (strcmp(key, "wind_speed") == 0) params->wind_speed = v;
        else if (strcmp(key, "wind_coeff") == 0) params->wind_coeff = v;
        else if (strcmp(key, "schmidt_exp") == 0) params->schmidt_exp = v;
        else if (strcmp(key, "current_k_factor") == 0) params->current_k_factor = v;
        else if (strcmp(key, "benthic_resp_20C") == 0) params->benthic_resp_20C = v;
        else if (strcmp(key, "benthic_Q10") == 0) params->benthic_Q10 = v;
        /* RIVE organic matter parameters */
        else if (strcmp(key, "khydr1") == 0) params->khydr1 = v;
        else if (strcmp(key, "khydr2") == 0) params->khydr2 = v;
        else if (strcmp(key, "khydr3") == 0) params->khydr3 = v;
        else if (strcmp(key, "frac_hd1") == 0) params->frac_hd1 = v;
        else if (strcmp(key, "frac_hd2") == 0) params->frac_hd2 = v;
        else if (strcmp(key, "frac_hp1") == 0) params->frac_hp1 = v;
        else if (strcmp(key, "frac_hp2") == 0) params->frac_hp2 = v;
        /* RIVE bacteria parameters */
        else if (strcmp(key, "bag_bmax20") == 0) params->bag_bmax20 = v;
        else if (strcmp(key, "bap_bmax20") == 0) params->bap_bmax20 = v;
        else if (strcmp(key, "bag_kdb20") == 0) params->bag_kdb20 = v;
        else if (strcmp(key, "bap_kdb20") == 0) params->bap_kdb20 = v;
        else if (strcmp(key, "bac_ks") == 0) params->bac_ks = v;
        else if (strcmp(key, "bac_yield") == 0) params->bac_yield = v;
        else if (strcmp(key, "bag_topt") == 0) params->bag_topt = v;
        else if (strcmp(key, "bap_topt") == 0) params->bap_topt = v;
        else if (strcmp(key, "bag_dti") == 0) params->bag_dti = v;
        else if (strcmp(key, "bap_dti") == 0) params->bap_dti = v;
        else if (strcmp(key, "bag_vs") == 0) params->bag_vs = v;
        /* RIVE phosphorus parameters */
        else if (strcmp(key, "kpads") == 0) params->kpads = v;
        else if (strcmp(key, "pac") == 0) params->pac = v;
    }
    
    fclose(fp);
    return 0;
}

/* -------------------------------------------------------------------------- */

typedef enum {
    TEMP_KIN_PBMAX = 0,
    TEMP_KIN_KMAINT,
    TEMP_KIN_KMORT,
    TEMP_KIN_KHET,
    TEMP_KIN_KDENIT,
    TEMP_KIN_KNIT
} TempKinID;

/**
 * Temperature-dependent rate function (matches Fortran Temp_kin)
 */
static double temp_kin(TempKinID id, double base_rate, double temp) {
    double delta = temp - 20.0;
    double factor = 1.0;

    switch (id) {
        case TEMP_KIN_KNIT:
            factor = pow(1.08, delta);
            break;
        case TEMP_KIN_KDENIT:
            factor = pow(1.07, delta);
            break;
        case TEMP_KIN_KHET:
            factor = pow(2.5, delta / 10.0);
            break;
        case TEMP_KIN_KMORT:
            factor = exp(0.07 * (delta + 20.0));
            break;
        case TEMP_KIN_KMAINT:
            factor = exp(0.0322 * delta);
            break;
        case TEMP_KIN_PBMAX:
        default:
            factor = pow(1.067, delta);
            break;
    }

    return base_rate * factor;
}

/* ===========================================================================
 * RIVE Temperature and Kinetic Functions
 * Reference: QualNET chem_utils.cpp
 * ===========================================================================*/

/**
 * RIVE bell-curve temperature function (ftp)
 * More realistic than Q10 for biological processes with optimal temperatures
 * 
 * ftp(T) = exp(-((T - Topt) / dti)²)
 * 
 * @param temp Current temperature [°C]
 * @param topt Optimal temperature [°C]
 * @param dti Temperature spread parameter [°C]
 * @return Temperature factor [0-1]
 * 
 * Reference: Billen et al. (1994), RIVE chem_utils.cpp
 */
static double rive_ftp(double temp, double topt, double dti) {
    if (dti <= 0.0) return 1.0;
    double delta = (temp - topt) / dti;
    return exp(-delta * delta);
}

/**
 * Michaelis-Menten kinetics
 * 
 * mich(S, K) = S / (S + K)
 * 
 * @param substrate Substrate concentration
 * @param k_half Half-saturation constant
 * @return Limitation factor [0-1]
 */
static double rive_mich(double substrate, double k_half) {
    if (substrate <= 0.0) return 0.0;
    if (k_half <= 0.0) return 1.0;
    return substrate / (substrate + k_half);
}

/**
 * RIVE Phosphorus adsorption equilibrium
 * Calculates dissolved PO4 from total inorganic P (PIT = PO4 + PIP)
 * Based on Langmuir-type adsorption on SPM
 * 
 * Reference: RIVE phosphorus.cpp
 * 
 * @param pit Total inorganic P [µmol/L]
 * @param spm Suspended particulate matter [mg/L]
 * @param kpads Adsorption constant
 * @param pac P adsorption capacity
 * @return Dissolved PO4 [µmol/L]
 */
static double rive_phosphorus_equilibrium(double pit, double spm, double kpads, double pac) {
    if (pit <= 0.0) return 0.0;
    if (spm <= 0.0) return pit;  /* No SPM = all dissolved */
    
    /* From RIVE phosphorus.cpp:
     * rp = kpads - pit + mes * pac
     * po4 = (sqrt(rp² + 4*pit*kpads) - rp) / 2
     */
    double rp = kpads - pit + spm * pac;
    double po4 = (sqrt(rp * rp + 4.0 * pit * kpads) - rp) / 2.0;
    
    /* Ensure physical bounds */
    if (po4 < 0.0) po4 = 0.0;
    if (po4 > pit) po4 = pit;
    
    return po4;
}

/**
 * Light limitation function
 * Port of the Fortran LightLim implementation
 */
static double light_limitation(double alpha, double pbmax, double I0, double kd, double depth) {
    const double EULER_GAMMA = 0.5772156649015328606;
    if (depth < CGEM_MIN_DEPTH) depth = CGEM_MIN_DEPTH;
    if (alpha <= 0.0 || pbmax <= 0.0 || kd <= 0.0 || I0 <= 0.0) {
        return 0.0;
    }

    double Ebottom = I0 * exp(-kd * depth);
    if (Ebottom < 1e-30) {
        return 0.0;
    }

    double psurf = (I0 * alpha) / pbmax;
    double pbot = (Ebottom * alpha) / pbmax;
    if (psurf <= 0.0 || pbot <= 0.0) {
        return 0.0;
    }

    double ps2 = psurf * psurf;
    double ps3 = ps2 * psurf;
    double ps4 = ps3 * psurf;
    double ps5 = ps4 * psurf;

    double pb2 = pbot * pbot;
    double pb3 = pb2 * pbot;
    double pb4 = pb3 * pbot;
    double pb5 = pb4 * pbot;

    double appGAMMAsurf;
    double appGAMMAbot;

    if (psurf <= 1.0) {
        appGAMMAsurf = -(log(psurf) + EULER_GAMMA - psurf + (ps2 / 4.0) - (ps3 / 18.0) + (ps4 / 96.0) - (ps5 / 600.0));
    } else {
        double frac = psurf + 9.0;
        frac = psurf + 7.0 - (16.0 / frac);
        frac = psurf + 5.0 - (9.0 / frac);
        frac = psurf + 3.0 - (4.0 / frac);
        frac = psurf + 1.0 - (1.0 / frac);
        appGAMMAsurf = exp(-psurf) * (1.0 / frac);
    }

    if (pbot <= 1.0) {
        appGAMMAbot = -(log(pbot) + EULER_GAMMA + pbot - (pb2 / 4.0) + (pb3 / 18.0) - (pb4 / 96.0) + (pb5 / 600.0));
    } else {
        double frac = pbot + 9.0;
        frac = pbot + 7.0 - (16.0 / frac);
        frac = pbot + 5.0 - (9.0 / frac);
        frac = pbot + 3.0 - (4.0 / frac);
        frac = pbot + 1.0 - (1.0 / frac);
        appGAMMAbot = exp(-pbot) * (1.0 / frac);
    }

    double integral = (appGAMMAsurf - appGAMMAbot + log(I0 / Ebottom)) / kd;
    if (!isfinite(integral) || integral < 0.0) {
        return 0.0;
    }
    return integral;
}

/**
 * Oxygen saturation concentration [µM]
 */
static double oxygen_saturation(double temp, double salinity) {
    /* Simplified Weiss (1970) equation */
    double A1 = -173.4292, A2 = 249.6339, A3 = 143.3483, A4 = -21.8492;
    double B1 = -0.033096, B2 = 0.014259, B3 = -0.0017;
    
    double Ts = log((298.15 - temp) / (273.15 + temp));
    double O2_sat = exp(A1 + A2*Ts + A3*Ts*Ts + A4*Ts*Ts*Ts + 
                       salinity*(B1 + B2*Ts + B3*Ts*Ts));
    return O2_sat; /* [µmol/kg] */
}

/**
 * Enhanced piston velocity for CO2 gas exchange [m/s]
 * Combines wind-driven (Wanninkhof 2014) and current-driven components
 * Critical for accurate pCO2 flux calculation in estuaries
 * 
 * @param velocity Water velocity [m/s]
 * @param depth Water depth [m]
 * @param salinity Salinity [PSU]
 * @param temp Temperature [°C]
 * @return Gas transfer velocity k [m/s]
 */
static double piston_velocity(double velocity, double depth, double salinity, double temp) {
    if (depth < CGEM_MIN_DEPTH) depth = CGEM_MIN_DEPTH;
    (void)salinity; /* Salinity affects Schmidt number minimally for CO2 */
    
    /* Schmidt number for CO2 in seawater (Wanninkhof 2014, Table 1) */
    /* Sc = A - B*T + C*T² - D*T³ */
    double Sc = 2116.8 - 136.25*temp + 4.7353*temp*temp - 0.092307*temp*temp*temp 
                + 0.0007555*temp*temp*temp*temp;
    if (Sc < 100.0) Sc = 100.0;  /* Prevent unrealistic values at high T */
    
    /* Wind-driven component: k_wind = a * U10² * (Sc/660)^n */
    /* Wanninkhof (2014): a = 0.251 for long-term average winds */
    double U10 = g_biogeo_params.wind_speed;  /* Wind speed at 10m [m/s] */
    double a = g_biogeo_params.wind_coeff;    /* 0.251 typical */
    double n = g_biogeo_params.schmidt_exp;   /* -0.5 typical */
    
    /* k in cm/hr, then convert to m/s */
    double k_wind = a * U10 * U10 * pow(Sc / 660.0, n);  /* [cm/hr] */
    
    /* Current-driven component (O'Connor & Dobbins type, important in rivers) */
    /* k_current ~ factor * sqrt(U/H) */
    double abs_vel = fabs(velocity);
    double k_current = g_biogeo_params.current_k_factor * sqrt(abs_vel / depth) * 100.0 * 3600.0;  /* [cm/hr] */
    
    /* Combined transfer velocity */
    double k_total = sqrt(k_wind * k_wind + k_current * k_current);  /* Quadratic combination */
    
    /* Convert from cm/hr to m/s */
    return k_total / (100.0 * 3600.0);
}

/**
 * Henry's constant for CO2 (µmol m⁻³ µatm⁻¹)
 */
static double henry_co2(double temp, double salinity) {
    double Tabs = temp + 273.15;
    if (Tabs <= 0.0) Tabs = 273.15;
    double lnK0 = -574.70126 + (21541.52 / Tabs) - 0.000147759 * Tabs * Tabs + 89.892 * log(Tabs);
    double f = 0.029941 - 0.00027455 * Tabs + 0.00000053407 * Tabs * Tabs;
    return exp(lnK0 + f * salinity) * 1e6;
}

static void carbonate_constants(double temp, double salinity,
                                double *K1, double *K2, double *Kw,
                                double *KH2S, double *KBOH4) {
    double Tabs = temp + 273.15;
    if (Tabs <= 0.0) Tabs = 273.15;
    double sqrtS = sqrt(CGEM_MAX(salinity, 0.0));
    double sal = CGEM_MAX(salinity, 0.0);

    double pK = -14.8425 + (3404.71 / Tabs) + (0.032786 * Tabs);
    double f1 = -0.0230848 - (14.3456 / Tabs);
    double f2 = 0.000691881 + (0.429955 / Tabs);
    *K1 = pow(10.0, -(pK + (f1 * sqrtS) + (f2 * sal)));

    pK = -6.4980 + (2902.39 / Tabs) + (0.02379 * Tabs);
    f1 = -0.458898 + (41.24048 / Tabs);
    f2 = 0.0284743 - (2.55895 / Tabs);
    *K2 = pow(10.0, -(pK + (f1 * sqrtS) + (f2 * sal)));

    pK = (-13847.26 / Tabs) + 148.9652 - (23.6521 * log(Tabs));
    f1 = (118.67 / Tabs) - 5.977 + 1.0495 * log(Tabs);
    f2 = -0.01615;
    *Kw = exp(pK + (f1 * sqrtS) + (f2 * sal));

    pK = (-4276.1 / Tabs) + 141.328 - (23.093 * log(Tabs));
    double ionic = (19.919 * sal) / CGEM_MAX(1000.0 - 1.00198 * sal, 1.0);
    double lnKH2S = (((-13856.0 / Tabs) + 324.57 - 47.986 * log(Tabs)) * sqrt(ionic)) +
                    (((35474.0 / Tabs) - 771.54 + 114.723 * log(Tabs)) * ionic) -
                    ((2698.0 / Tabs) * pow(ionic, 1.5)) +
                    ((1776.0 / Tabs) * ionic * ionic) + pK;
    *KH2S = exp(lnKH2S);

    double term1 = (-8966.90 - 2890.53 * sqrtS - 77.942 * sal + 1.728 * pow(sal, 1.5) - 0.0996 * sal * sal) / Tabs;
    double term2 = (148.0248 + 137.1942 * sqrtS + 1.62142 * sal);
    double term3 = (-24.4344 - 25.085 * sqrtS - 0.2474 * sal) * log(Tabs);
    double term4 = 0.053105 * sqrtS * Tabs;
    *KBOH4 = exp(term1 + term2 + term3 + term4);
}

static double carbonate_residual(double H, double dic, double alk, double sulf, double TotB,
                                 double K1, double K2, double Kw, double KH2S, double KBOH4,
                                 double *xi1, double *xi2, double *gamma1, double *beta1) {
    double Hsafe = CGEM_MAX(H, 1e-12);
    double denom = (Hsafe * Hsafe) + (K1 * Hsafe) + (K1 * K2);
    if (denom <= 0.0) denom = 1e-12;

    *xi1 = (K1 * Hsafe) / denom;
    *xi2 = (*xi1) * K2 / Hsafe;
    *gamma1 = KH2S / (KH2S + Hsafe);
    *beta1 = KBOH4 / (KBOH4 + Hsafe);

    double residual = alk - (((*xi1 + 2.0 * (*xi2)) * dic) + ((*gamma1) * sulf) + ((*beta1) * TotB) +
                              (Kw / Hsafe) - Hsafe);
    return residual;
}

static double calculate_carbonate_system(Branch *branch, int idx, double temp,
                                         double depth, double salinity, double piston_vel) {
    double *dic = branch->conc[CGEM_SPECIES_DIC];
    double *alk = branch->conc[CGEM_SPECIES_AT];
    double *co2 = branch->conc[CGEM_SPECIES_CO2];
    double *pco2 = branch->conc[CGEM_SPECIES_PCO2];
    double *ph = branch->conc[CGEM_SPECIES_PH];
    double *alkc = branch->conc[CGEM_SPECIES_ALKC];
    double *hs = branch->conc[CGEM_SPECIES_HS];

    double henry = henry_co2(temp, salinity);
    if (!dic || !alk || !co2 || !pco2 || !ph || !alkc || !branch->reaction_rates) {
        return henry;
    }

    double sal = CGEM_MAX(salinity, 0.0);
    double Hprot = (ph[idx] > 0.0) ? pow(10.0, -ph[idx]) : 1e-8;
    double sulf = (hs) ? CGEM_MAX(hs[idx], 0.0) : 0.0;
    double TotB = 416.0 * (sal / 35.0);

    double K1, K2, Kw, KH2S, KBOH4;
    carbonate_constants(temp, sal, &K1, &K2, &Kw, &KH2S, &KBOH4);

    double xi1 = 0.0, xi2 = 0.0, gamma1 = 0.0, beta1 = 0.0;
    double guess = carbonate_residual(Hprot, dic[idx], alk[idx], sulf, TotB, K1, K2, Kw, KH2S, KBOH4,
                                      &xi1, &xi2, &gamma1, &beta1);

    double newH = CGEM_MAX(Hprot, 1e-9);
    int iterations = 0;
    while (fabs(guess) > CGEM_TOL && iterations < 50) {
        double carbonAlk = alk[idx] - (gamma1 * sulf) - (beta1 * TotB) - (Kw / CGEM_MAX(newH, 1e-12)) + newH;
        double ratio = dic[idx] / CGEM_MAX(carbonAlk, 1e-12);
        double inner = 1.0 - ratio;
        double discriminant = inner * inner * (K1 * K1) - 4.0 * K1 * K2 * (1.0 - 2.0 * ratio);
        if (discriminant < 0.0) discriminant = 0.0;
        double update = 0.5 * (((ratio) - 1.0) * K1 + sqrt(discriminant));
        if (!isfinite(update) || update <= 0.0) {
            break;
        }
        newH = update;
        guess = carbonate_residual(newH, dic[idx], alk[idx], sulf, TotB, K1, K2, Kw, KH2S, KBOH4,
                                    &xi1, &xi2, &gamma1, &beta1);
        iterations++;
    }

    Hprot = CGEM_MAX(newH, 1e-9);
    ph[idx] = -log10(Hprot);
    alkc[idx] = (double)iterations;

    double denom = 1.0 + (K1 / Hprot) + ((K1 * K2) / (Hprot * Hprot));
    double co2_spec = (denom > 0.0) ? dic[idx] / denom : dic[idx];
    if (!isfinite(co2_spec) || co2_spec < 0.0) {
        co2_spec = CGEM_MAX(dic[idx], 0.0);
    }
    co2[idx] = co2_spec;

    double depth_eff = (depth < CGEM_MIN_DEPTH) ? CGEM_MIN_DEPTH : depth;
    double co2_flux = 0.0;
    if (depth_eff > 0.0) {
        co2_flux = -(piston_vel * 0.913 / depth_eff) * (co2_spec - henry * branch->pco2_atm);
    }
    branch->reaction_rates[CGEM_REACTION_CO2_EX][idx] = co2_flux;
    branch->reaction_rates[CGEM_REACTION_CO2_EX_S][idx] = co2_flux * depth_eff;
    branch->reaction_rates[CGEM_REACTION_HS_REACT][idx] = 0.0;

    double aer_deg = branch->reaction_rates[CGEM_REACTION_AER_DEG][idx];
    double denit = branch->reaction_rates[CGEM_REACTION_DENIT][idx];
    double npp_total = branch->reaction_rates[CGEM_REACTION_NPP][idx];
    double npp_no3 = branch->reaction_rates[CGEM_REACTION_NPP_NO3][idx];
    double npp_nh4 = branch->reaction_rates[CGEM_REACTION_NPP_NH4][idx];
    double nitrif = branch->reaction_rates[CGEM_REACTION_NIT][idx];

    /* =========================================================================
     * BENTHIC RESPIRATION (major pCO2 driver in shallow estuaries)
     * CRITICAL: Proper unit conversion from [mmol/m²/day] to [µM/s]
     * ========================================================================= */
    
    /* Temperature-corrected benthic respiration [mmol C/m²/day] */
    double benthic_rate_day = g_biogeo_params.benthic_resp_20C * 
                              pow(g_biogeo_params.benthic_Q10, (temp - 20.0) / 10.0);
    
    /* Convert [mmol/m²/day] → [µM/s] (mmol/m³/s):
     * 1. Divide by depth (m) to get volumetric rate [mmol/m³/day]
     * 2. Divide by 86400 to convert day to seconds */
    double benthic_co2_source = benthic_rate_day / depth_eff / SECONDS_PER_DAY;

    /* DIC source = CO2 flux + water column respiration + benthic + denitrification - NPP */
    /* Note: aer_deg already converted to [1/s] in loop 1, so units are consistent */
    branch->reaction_rates[CGEM_REACTION_DIC_REACT][idx] = co2_flux + aer_deg + denit 
                                                           + benthic_co2_source - npp_total;
    
    branch->reaction_rates[CGEM_REACTION_TA_REACT][idx] = (15.0 / 106.0) * aer_deg +
                                                          (93.4 / 106.0) * denit -
                                                          2.0 * nitrif +
                                                          (-15.0 / 106.0) * npp_nh4 +
                                                          (17.0 / 106.0) * npp_no3;

    (void)pco2; /* Updated later after flux integration */
    return henry;
}

/**
 * Initialize biogeochemical parameters for a branch
 * 
 * If the branch has a biogeo_params_path specified, load those values
 * (overriding global defaults). This enables spatial heterogeneity in
 * water quality processes (e.g., high kox in urban Saigon vs low kox
 * in rural Dong Nai).
 * 
 * Scientific motivation:
 *   - Urban rivers typically have higher organic matter decay rates
 *     (kox ~ 0.2-0.3 /day) due to high inputs and warmer temperatures
 *   - Rural/natural rivers have lower decay rates (kox ~ 0.05-0.1 /day)
 *   - Mangrove estuaries have high benthic respiration (60-100 mmol C/m²/day)
 * 
 * Reference: Volta et al. (2016), Regnier et al. (2013)
 */
void InitializeBiogeoParameters(Branch *branch) {
    if (!branch) return;
    
    /* Ensure global params are initialized */
    if (!g_biogeo_params.loaded) {
        LoadBiogeoParams(NULL);  /* Load defaults */
    }
    
    /* Start with global parameters */
    BiogeoParams params = g_biogeo_params;
    
    /* Override with branch-specific parameters if specified */
    if (branch->biogeo_params_path[0] != '\0') {
        /* Load branch-specific params (overlay on global) */
        if (LoadBiogeoParamsToStruct(branch->biogeo_params_path, &params) == 0) {
            printf("  Branch '%s' using custom biogeo parameters\n", branch->name);
        } else {
            fprintf(stderr, "Warning: Could not load branch biogeo params from '%s', using global defaults\n",
                    branch->biogeo_params_path);
        }
    }
    
    /* Water and environmental parameters */
    branch->water_temp = params.water_temp;
    branch->ws = params.ws;
    
    /* Light parameters */
    branch->I0 = params.I0;
    branch->kd1 = params.kd1;
    branch->kd2_spm = params.kd2_spm;
    branch->kd2_phy1 = params.kd2_phy1;
    branch->kd2_phy2 = params.kd2_phy2;
    
    /* Phytoplankton parameters */
    branch->alpha1 = params.alpha1;
    branch->alpha2 = params.alpha2;
    branch->pbmax1 = params.pbmax1;
    branch->pbmax2 = params.pbmax2;
    branch->kexc1 = params.kexc1;
    branch->kexc2 = params.kexc2;
    branch->kgrowth1 = params.kgrowth1;
    branch->kgrowth2 = params.kgrowth2;
    branch->kmaint1 = params.kmaint1;
    branch->kmaint2 = params.kmaint2;
    branch->kmort1 = params.kmort1;
    branch->kmort2 = params.kmort2;
    
    /* Nutrient limitation parameters */
    branch->kdsi1 = params.kdsi1;
    branch->kn1 = params.kn1;
    branch->kpo41 = params.kpo41;
    branch->kn2 = params.kn2;
    branch->kpo42 = params.kpo42;
    
    /* Decomposition parameters */
    branch->kox = params.kox;
    branch->kdenit = params.kdenit;
    branch->knit = params.knit;
    branch->ktox = params.ktox;
    branch->ko2 = params.ko2;
    branch->ko2_nit = params.ko2_nit;
    branch->kno3 = params.kno3;
    branch->knh4 = params.knh4;
    branch->kino2 = params.kino2;
    
    /* Stoichiometric ratios */
    branch->redn = params.redn;
    branch->redp = params.redp;
    branch->redsi = params.redsi;
    
    /* Gas exchange */
    branch->pco2_atm = params.pco2_atm;
    
    /* RIVE organic matter parameters */
    branch->khydr1 = g_biogeo_params.khydr1;
    branch->khydr2 = g_biogeo_params.khydr2;
    branch->khydr3 = g_biogeo_params.khydr3;
    branch->frac_hd1 = g_biogeo_params.frac_hd1;
    branch->frac_hd2 = g_biogeo_params.frac_hd2;
    branch->frac_hp1 = g_biogeo_params.frac_hp1;
    branch->frac_hp2 = g_biogeo_params.frac_hp2;
    
    /* RIVE bacteria parameters */
    branch->bag_bmax20 = g_biogeo_params.bag_bmax20;
    branch->bap_bmax20 = g_biogeo_params.bap_bmax20;
    branch->bag_kdb20 = g_biogeo_params.bag_kdb20;
    branch->bap_kdb20 = g_biogeo_params.bap_kdb20;
    branch->bac_ks = g_biogeo_params.bac_ks;
    branch->bac_yield = g_biogeo_params.bac_yield;
    branch->bag_topt = g_biogeo_params.bag_topt;
    branch->bap_topt = g_biogeo_params.bap_topt;
    branch->bag_dti = g_biogeo_params.bag_dti;
    branch->bap_dti = g_biogeo_params.bap_dti;
    branch->bag_vs = g_biogeo_params.bag_vs;
    
    /* RIVE phosphorus adsorption parameters */
    branch->kpads = g_biogeo_params.kpads;
    branch->pac = g_biogeo_params.pac;
    
    /* RIVE benthic parameters */
    branch->zf_init = g_biogeo_params.zf_init;
    branch->benthic_porosity = g_biogeo_params.benthic_porosity;
    branch->benthic_density = g_biogeo_params.benthic_density;
}
int Biogeo_Branch(Branch *branch, double dt) {
    if (!branch || branch->M <= 0 || branch->num_species < CGEM_NUM_SPECIES || !branch->reaction_rates) {
        return 0;
    }

    int M = branch->M;
    double temp = branch->water_temp;

    double *salinity = branch->conc[CGEM_SPECIES_SALINITY];
    double *phy1 = branch->conc[CGEM_SPECIES_PHY1];
    double *phy2 = branch->conc[CGEM_SPECIES_PHY2];
    double *dsi = branch->conc[CGEM_SPECIES_DSI];
    double *no3 = branch->conc[CGEM_SPECIES_NO3];
    double *nh4 = branch->conc[CGEM_SPECIES_NH4];
    double *po4 = branch->conc[CGEM_SPECIES_PO4];
    double *o2 = branch->conc[CGEM_SPECIES_O2];
    double *toc = branch->conc[CGEM_SPECIES_TOC];
    double *spm = branch->conc[CGEM_SPECIES_SPM];
    double *dic = branch->conc[CGEM_SPECIES_DIC];
    double *at = branch->conc[CGEM_SPECIES_AT];
    double *co2 = branch->conc[CGEM_SPECIES_CO2];
    double *pco2 = branch->conc[CGEM_SPECIES_PCO2];
    double *ph = branch->conc[CGEM_SPECIES_PH];
    double *alkc = branch->conc[CGEM_SPECIES_ALKC];

    if (!salinity || !phy1 || !phy2 || !dsi || !no3 || !nh4 || !po4 || !o2 || !toc ||
        !dic || !at || !co2 || !pco2 || !ph || !alkc) {
        return 0;
    }

    double *pis_vel = (double *)calloc((size_t)(M + 2), sizeof(double));
    if (!pis_vel) {
        return -1;
    }
    for (int i = 1; i <= M; ++i) {
        pis_vel[i] = piston_velocity(branch->velocity[i], branch->depth[i], salinity[i], temp);
    }

    /* Primary loop: compute reaction rates */
    for (int i = 1; i <= M - 1; i += 2) {
        double depth = CGEM_MAX(branch->depth[i], CGEM_MIN_DEPTH);

        double sal = salinity[i];
        double no3nh4 = no3[i] + nh4[i];

        double kd = branch->kd1 +
                    branch->kd2_spm * (spm ? spm[i] * 1e-3 : 0.0) +
                    branch->kd2_phy1 * phy1[i] +
                    branch->kd2_phy2 * phy2[i];

        double pbmax1 = temp_kin(TEMP_KIN_PBMAX, branch->pbmax1, temp);
        double pbmax2 = temp_kin(TEMP_KIN_PBMAX, branch->pbmax2, temp);
        double lightlim1 = light_limitation(branch->alpha1, pbmax1, branch->I0, kd, depth);
        double lightlim2 = light_limitation(branch->alpha2, pbmax2, branch->I0, kd, depth);

        double dsi_term = dsi[i] + branch->kdsi1;
        double n_term1 = no3nh4 + branch->kn1;
        double n_term2 = no3nh4 + branch->kn2;
        double po4_term1 = po4[i] + branch->kpo41;
        double po4_term2 = po4[i] + branch->kpo42;

        double nlim1 = 0.0;
        if (dsi_term > 0.0 && n_term1 > 0.0 && po4_term1 > 0.0) {
            nlim1 = (dsi[i] / dsi_term) * ((no3nh4) / n_term1) * (po4[i] / po4_term1);
        }
        double nlim2 = 0.0;
        if (n_term2 > 0.0 && po4_term2 > 0.0) {
            nlim2 = (no3nh4 / n_term2) * (po4[i] / po4_term2);
        }

        double nswitch = nh4[i] / (5.0 + CGEM_MAX(nh4[i], 0.0));

        double gpp1 = pbmax1 * phy1[i] * nlim1 * lightlim1;
        double gpp2 = pbmax2 * phy2[i] * nlim2 * lightlim2;

        double maint1 = temp_kin(TEMP_KIN_KMAINT, branch->kmaint1, temp) * phy1[i];
        double maint2 = temp_kin(TEMP_KIN_KMAINT, branch->kmaint2, temp) * phy2[i];

        double npp_no31 = (1.0 - nswitch) * (gpp1 / depth) * (1.0 - branch->kexc1) * (1.0 - branch->kgrowth1) - maint1;
        double npp_no32 = (1.0 - nswitch) * (gpp2 / depth) * (1.0 - branch->kexc2) * (1.0 - branch->kgrowth2) - maint2;
        double npp_nh41 = nswitch * (gpp1 / depth) * (1.0 - branch->kexc1) * (1.0 - branch->kgrowth1) - maint1;
        double npp_nh42 = nswitch * (gpp2 / depth) * (1.0 - branch->kexc2) * (1.0 - branch->kgrowth2) - maint2;

        double mort1 = temp_kin(TEMP_KIN_KMORT, branch->kmort1, temp) * phy1[i];
        double mort2 = temp_kin(TEMP_KIN_KMORT, branch->kmort2, temp) * phy2[i];

        double si_cons = fmin(0.0, -branch->redsi * (npp_nh41 + npp_nh42));

        /* CRITICAL: Convert rates from [1/day] to [1/s] */
        double kox_s = branch->kox / SECONDS_PER_DAY;
        double kdenit_s = branch->kdenit / SECONDS_PER_DAY;
        double knit_s = branch->knit / SECONDS_PER_DAY;
        
        double aer_deg = temp_kin(TEMP_KIN_KHET, kox_s, temp) *
                         (toc[i] / (toc[i] + branch->ktox)) *
                         (o2[i] / (o2[i] + branch->ko2));

        double denit = temp_kin(TEMP_KIN_KDENIT, kdenit_s, temp) *
                       (toc[i] / (toc[i] + branch->ktox)) *
                       (branch->kino2 / (o2[i] + branch->kino2)) *
                       (no3[i] / (no3[i] + branch->kno3));

        double nitrif = temp_kin(TEMP_KIN_KNIT, knit_s, temp) *
                        (o2[i] / (o2[i] + branch->ko2_nit)) *
                        (nh4[i] / (nh4[i] + branch->knh4));

        double o2_sat = oxygen_saturation(temp, sal);
        double o2_ex = (pis_vel[i] / depth) * (o2_sat - o2[i]);

        double npp_no3_sum = npp_no31 + npp_no32;
        double npp_nh4_sum = npp_nh41 + npp_nh42;
        double npp_total = npp_no3_sum + npp_nh4_sum;
        double phy_death = mort1 + mort2;

        branch->reaction_rates[CGEM_REACTION_GPP_1][i] = gpp1;
        branch->reaction_rates[CGEM_REACTION_GPP_2][i] = gpp2;
        branch->reaction_rates[CGEM_REACTION_NPP_NO3_1][i] = npp_no31;
        branch->reaction_rates[CGEM_REACTION_NPP_NO3_2][i] = npp_no32;
        branch->reaction_rates[CGEM_REACTION_NPP_NH4_1][i] = npp_nh41;
        branch->reaction_rates[CGEM_REACTION_NPP_NH4_2][i] = npp_nh42;
        branch->reaction_rates[CGEM_REACTION_NPP_NO3][i] = npp_no3_sum;
        branch->reaction_rates[CGEM_REACTION_NPP_NH4][i] = npp_nh4_sum;
        branch->reaction_rates[CGEM_REACTION_NPP][i] = npp_total;
        branch->reaction_rates[CGEM_REACTION_PHY_DEATH_1][i] = mort1;
        branch->reaction_rates[CGEM_REACTION_PHY_DEATH_2][i] = mort2;
        branch->reaction_rates[CGEM_REACTION_PHY_DEATH][i] = phy_death;
        branch->reaction_rates[CGEM_REACTION_SI_CONS][i] = si_cons;
        branch->reaction_rates[CGEM_REACTION_AER_DEG][i] = aer_deg;
        branch->reaction_rates[CGEM_REACTION_DENIT][i] = denit;
        branch->reaction_rates[CGEM_REACTION_NIT][i] = nitrif;
        branch->reaction_rates[CGEM_REACTION_O2_EX][i] = o2_ex;
        branch->reaction_rates[CGEM_REACTION_O2_EX_S][i] = o2_ex * depth;
    }

    /* Second loop: carbonate system and state updates */
    for (int i = 1; i <= M - 1; i += 2) {
        double depth = CGEM_MAX(branch->depth[i], CGEM_MIN_DEPTH);
        double sal = salinity[i];

        double henry = calculate_carbonate_system(branch, i, temp, depth, sal, pis_vel[i]);

        double npp_no31 = branch->reaction_rates[CGEM_REACTION_NPP_NO3_1][i];
        double npp_no32 = branch->reaction_rates[CGEM_REACTION_NPP_NO3_2][i];
        double npp_nh41 = branch->reaction_rates[CGEM_REACTION_NPP_NH4_1][i];
        double npp_nh42 = branch->reaction_rates[CGEM_REACTION_NPP_NH4_2][i];
        double npp_no3 = branch->reaction_rates[CGEM_REACTION_NPP_NO3][i];
        double npp_nh4 = branch->reaction_rates[CGEM_REACTION_NPP_NH4][i];
        double npp_total = branch->reaction_rates[CGEM_REACTION_NPP][i];
        double mort1 = branch->reaction_rates[CGEM_REACTION_PHY_DEATH_1][i];
        double mort2 = branch->reaction_rates[CGEM_REACTION_PHY_DEATH_2][i];
        double phy_death = branch->reaction_rates[CGEM_REACTION_PHY_DEATH][i];
        double si_cons = branch->reaction_rates[CGEM_REACTION_SI_CONS][i];
        double aer_deg = branch->reaction_rates[CGEM_REACTION_AER_DEG][i];
        double denit = branch->reaction_rates[CGEM_REACTION_DENIT][i];
        double nitrif = branch->reaction_rates[CGEM_REACTION_NIT][i];
        double o2_ex = branch->reaction_rates[CGEM_REACTION_O2_EX][i];

        phy1[i] = CGEM_MAX(0.0, phy1[i] + (npp_no31 + npp_nh41 - mort1) * dt);
        phy2[i] = CGEM_MAX(0.0, phy2[i] + (npp_no32 + npp_nh42 - mort2) * dt);
        dsi[i] = CGEM_MAX(0.0, dsi[i] + si_cons * dt);
        po4[i] = CGEM_MAX(0.0, po4[i] + branch->redp * (aer_deg + denit - npp_total) * dt);
        o2[i] = CGEM_MAX(0.0, o2[i] + (-aer_deg + npp_nh4 + (138.0 / 106.0) * npp_no3 - 2.0 * nitrif + o2_ex) * dt);
        toc[i] = CGEM_MAX(0.0, toc[i] + (-aer_deg - denit + phy_death) * dt);
        nh4[i] = CGEM_MAX(0.0, nh4[i] + (branch->redn * (aer_deg - npp_nh4) - nitrif) * dt);
        no3[i] = CGEM_MAX(0.0, no3[i] + (-94.4 / 106.0 * denit + nitrif - branch->redn * npp_no3) * dt);

        dic[i] = CGEM_MAX(0.0, dic[i] + branch->reaction_rates[CGEM_REACTION_DIC_REACT][i] * dt);
        at[i] = CGEM_MAX(0.0, at[i] + branch->reaction_rates[CGEM_REACTION_TA_REACT][i] * dt);
        co2[i] = CGEM_MAX(0.0, co2[i] + (branch->reaction_rates[CGEM_REACTION_CO2_EX][i] + denit - npp_total) * dt);
        if (henry > 0.0) {
            pco2[i] = CGEM_MAX(0.0, co2[i] / henry);
        } else {
            pco2[i] = 0.0;
        }
    }

    free(pis_vel);
    return 0;
}
