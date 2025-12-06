/**
 * @file sediment_diagenesis.c
 * @brief Active Sediment Layer (Simplified Diagenesis) Implementation
 * 
 * CRITICAL MISSING FEATURE IMPLEMENTATION (December 2025 Audit)
 * 
 * This module implements a single active sediment layer that:
 * 1. Receives POC deposition from the water column
 * 2. Mineralizes SOC with temperature-dependent kinetics  
 * 3. Produces DYNAMIC benthic fluxes proportional to SOC pool
 * 
 * This allows the model to capture "legacy loads" - historical pollution
 * accumulated in sediments that continues to affect water quality even
 * after pollution is reduced.
 * 
 * Reference:
 * - Chapra (2008) Surface Water-Quality Modeling, Chapter 24
 * - DiToro (2001) Sediment Flux Modeling
 * - Soetaert et al. (1996) Est. Coast. Shelf Sci.
 */

#include "sediment_diagenesis.h"
#include "rive_params.h"
#include "../define.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* ============================================================================
 * GLOBAL SOC CONFIGURATION
 * ============================================================================*/
static SOCConfig g_soc_config = {0};
static int g_soc_initialized = 0;

/* Unit conversion constants */
#define SECONDS_PER_DAY   86400.0
#define UMOL_TO_MMOL      0.001
#define MG_TO_G           0.001
#define MOL_C_TO_G        12.0    /* g C per mol C */

/* ============================================================================
 * CONFIGURATION FUNCTIONS
 * ============================================================================*/

void soc_init_config(SOCConfig *config) {
    if (!config) return;
    
    /* SOC decay kinetics
     * k_soc ~ 0.001-0.01 /day for estuarine sediments
     * Reference: Middelburg et al. (1996), Boudreau (1996)
     * 
     * Mekong sediments are organic-rich → use moderate rate
     */
    config->k_soc_20C = 0.003;    /* [1/day] - gives ~1 year turnover */
    config->soc_Q10 = 2.0;        /* Standard Q10 for microbial processes */
    
    /* Initial conditions
     * SOC pool sizes from literature:
     * - Oligotrophic systems: 100-300 g C/m²
     * - Mesotrophic estuaries: 300-1000 g C/m²
     * - Eutrophic deltas: 500-2000 g C/m²
     * 
     * Mekong is mesotrophic-eutrophic
     */
    config->soc_init = 500.0;     /* [g C/m²] */
    config->soc_max = 5000.0;     /* [g C/m²] - prevents runaway accumulation */
    
    /* Burial rate
     * k_burial removes SOC to deep sediment (permanent burial)
     * Typical: 0.0001-0.001 /day (1-10% per year)
     */
    config->k_burial = 0.0001;    /* [1/day] - ~3.5% per year */
    
    /* Stoichiometry (from rive_params defaults) */
    config->soc_redn = 0.151;     /* mol N / mol C (Redfield ~ 16/106) */
    config->soc_redp = 0.0094;    /* mol P / mol C (Redfield ~ 1/106) */
    config->soc_RQ = 1.2;         /* mol CO2 / mol O2 (>1 for anaerobic contribution) */
    
    /* Anaerobic fraction
     * In well-oxygenated waters: 20-40% of benthic mineralization is anaerobic
     * In hypoxic waters: 50-80%
     * This affects CH4 production
     */
    config->f_anaerobic = 0.30;   /* 30% anaerobic */
    
    /* CH4 production
     * Only from anaerobic fraction of SOC decay
     * Typical yield: 0.3-0.7 mol CH4 / mol C (methanogenesis)
     */
    config->ch4_yield = 0.5;
    config->ch4_o2_inhibit = 50.0;  /* µmol/L - CH4 production inhibited above this */
    
    /* N2O production from sediment N cycling
     * Coupled nitrification-denitrification in sediments
     * Typical: 1-5% of N mineralized ends up as N2O
     */
    config->n2o_yield_sediment = 0.02;  /* 2% */
    
    /* Control flags */
    config->enable_soc = 1;       /* Enabled by default */
    config->loaded = 1;
    
    g_soc_initialized = 1;
}

SOCConfig* soc_get_config(void) {
    if (!g_soc_initialized) {
        soc_init_config(&g_soc_config);
    }
    return &g_soc_config;
}

/* ============================================================================
 * BRANCH INITIALIZATION
 * ============================================================================*/

void soc_init_branch(Branch *branch, SOCConfig *config) {
    if (!branch || !config) return;
    
    /* Initialize SOC pool (using bed_oc array from network_data.c)
     * 
     * SPATIAL VARIATION:
     * - Upstream (fine mud): Higher initial SOC
     * - Downstream (sandy): Lower initial SOC
     * 
     * We use a linear gradient from 50% at downstream to 150% at upstream
     */
    int M = branch->M;
    
    for (int i = 0; i <= M + 1; ++i) {
        /* Spatial factor: 0.5 at i=1 (downstream), 1.5 at i=M (upstream) */
        double s = (M > 1) ? (double)(i - 1) / (double)(M - 1) : 0.5;
        double spatial_factor = 0.5 + 1.0 * s;  /* Range: 0.5-1.5 */
        
        /* Initialize SOC pool [g C/m²] */
        branch->bed_oc[i] = config->soc_init * spatial_factor;
        
        /* Clamp to max */
        if (branch->bed_oc[i] > config->soc_max) {
            branch->bed_oc[i] = config->soc_max;
        }
    }
    
    printf("  SOC pools initialized: %.0f-%.0f g C/m² (downstream to upstream)\n",
           branch->bed_oc[1], branch->bed_oc[M]);
}

/* ============================================================================
 * DEPOSITION CALCULATION
 * ============================================================================*/

double soc_calc_deposition(Branch *branch, int i, double phy_death, 
                           double toc, double spm, double ws, double depth,
                           SOCConfig *config) {
    if (!branch || !config || depth < CGEM_MIN_DEPTH) return 0.0;
    
    BiogeoParams *p = rive_get_params();
    double poc_spm_ratio = (p && p->poc_spm_ratio > 0.0) ? p->poc_spm_ratio : 0.02;
    
    /* =========================================================================
     * POC DEPOSITION SOURCES
     * 
     * 1. Dead phytoplankton sinking
     *    - phy_death [µmol C/L/day] * f_settling
     *    - Not all dead phytoplankton settles (some stays suspended)
     * 
     * 2. SPM-bound POC settling
     *    - SPM [mg/L] * ws [m/s] * poc_spm_ratio [g C/g SPM]
     *    - Converts to [g C/m²/day]
     * 
     * 3. Direct POC settling (from TOC labile fraction)
     *    - Assumes 30% of labile TOC is particulate
     * =========================================================================*/
    
    /* Fraction of phy_death that settles (rest is dissolved or stays suspended) */
    double f_phy_settle = 0.5;
    
    /* Convert phy_death: [µmol C/L/day] → [g C/m²/day]
     * µmol/L/day * depth[m] * 12 g/mol * 1e-6 mol/µmol * 1000 L/m³ = g/m²/day
     * = µmol/L/day * depth * 0.012
     */
    double dep_phy = phy_death * f_phy_settle * depth * 0.012;  /* [g C/m²/day] */
    
    /* SPM-bound POC settling
     * SPM [mg/L] * ws [m/s] * poc_ratio * 86400 s/day * 0.001 g/mg = g/m²/day
     */
    double dep_spm = spm * ws * SECONDS_PER_DAY * MG_TO_G * poc_spm_ratio;  /* [g C/m²/day] */
    
    /* Direct labile TOC settling (30% of labile is particulate)
     * Assumes labile TOC ~ 15% of total TOC (from 2-pool model)
     */
    double toc_labile_frac = 0.15;
    double f_particulate = 0.30;
    double ws_poc = ws * 0.3;  /* POC settles slower than SPM */
    
    /* toc [µmol/L] * 0.012 → g/L, then * ws [m/s] * 86400 [s/day] = g/m²/day */
    double dep_toc = toc * toc_labile_frac * f_particulate * 0.012 * ws_poc * SECONDS_PER_DAY;
    
    /* Total deposition */
    double total_deposition = dep_phy + dep_spm + dep_toc;
    
    /* Sanity limits: max 50 g C/m²/day (very high but possible in blooms) */
    if (total_deposition > 50.0) total_deposition = 50.0;
    if (total_deposition < 0.0) total_deposition = 0.0;
    
    return total_deposition;
}

/* ============================================================================
 * FLUX CALCULATION (Core of simplified diagenesis)
 * ============================================================================*/

void soc_calc_fluxes(SOCState *state, double soc_current, double deposition,
                     double temp, double o2, SOCConfig *config) {
    if (!state || !config) return;
    
    /* Initialize state */
    state->soc_pool = soc_current;
    state->deposition_flux = deposition;
    
    /* =========================================================================
     * SOC MINERALIZATION RATE
     * 
     * First-order decay with Q10 temperature correction:
     *   R_min = k_soc(T) * SOC
     *   k_soc(T) = k_soc_20C * Q10^((T-20)/10)
     * 
     * Units: [g C/m²/day]
     * =========================================================================*/
    double k_T = config->k_soc_20C * pow(config->soc_Q10, (temp - 20.0) / 10.0);
    double mineralization = k_T * soc_current;  /* [g C/m²/day] */
    state->mineralization_rate = mineralization;
    
    /* =========================================================================
     * BURIAL RATE
     * 
     * Permanent removal to deep sediments
     * R_burial = k_burial * SOC
     * =========================================================================*/
    double burial = config->k_burial * soc_current;  /* [g C/m²/day] */
    state->burial_rate = burial;
    
    /* =========================================================================
     * AEROBIC vs ANAEROBIC PARTITIONING
     * 
     * Anaerobic fraction increases at low O2:
     *   f_anaer = f_anaerobic_base + (1 - f_anaerobic_base) * K_o2 / (O2 + K_o2)
     * 
     * At O2 = 0:     f_anaer → 1.0 (fully anaerobic)
     * At O2 = 200:   f_anaer → f_anaerobic_base (well oxygenated)
     * =========================================================================*/
    double K_o2 = config->ch4_o2_inhibit;  /* Half-saturation for O2 inhibition */
    double f_anaer = config->f_anaerobic + 
                     (1.0 - config->f_anaerobic) * K_o2 / (o2 + K_o2 + 1e-10);
    
    double aerobic_mineral = mineralization * (1.0 - f_anaer);
    double anaerobic_mineral = mineralization * f_anaer;
    
    /* =========================================================================
     * OUTPUT FLUX CALCULATIONS
     * 
     * All fluxes are proportional to mineralization rate, making them DYNAMIC!
     * =========================================================================*/
    
    /* Convert g C/m²/day → mmol C/m²/day (÷12 g/mol × 1000 mmol/mol = ×83.33) */
    double mmol_C_per_day = mineralization * 83.33;
    
    /* SOD: Sediment Oxygen Demand [mmol O2/m²/day]
     * Only aerobic mineralization consumes O2
     * Stoichiometry: 1 mol O2 per mol C oxidized
     */
    state->sod_flux = aerobic_mineral * 83.33;  /* [mmol O2/m²/day] */
    
    /* NH4 flux [mmol N/m²/day]
     * N:C ratio from stoichiometry
     */
    state->nh4_flux = mmol_C_per_day * config->soc_redn;
    
    /* PO4 flux [mmol P/m²/day]
     * P:C ratio from stoichiometry
     * Note: Actual P release depends on redox conditions (not modeled here)
     */
    state->po4_flux = mmol_C_per_day * config->soc_redp;
    
    /* DIC flux [mmol C/m²/day]
     * RQ > 1 accounts for anaerobic CO2 production
     * Aerobic: 1 mol CO2 per mol O2
     * Anaerobic: produces CO2 without O2 consumption
     */
    state->dic_flux = state->sod_flux * config->soc_RQ + 
                      anaerobic_mineral * 83.33;  /* Anaerobic CO2 */
    
    /* CH4 flux [µmol/m²/day]
     * Only from anaerobic mineralization
     * CH4 yield from methanogenesis ~ 50%
     */
    double ch4_production = anaerobic_mineral * 83.33 * config->ch4_yield * 1000.0;  /* µmol/m²/day */
    
    /* CH4 oxidation in sediment surface (50% of production)
     * This is a simplification - actual oxidation depends on O2 penetration
     */
    double f_ch4_oxidized = 0.5 * o2 / (o2 + 50.0 + 1e-10);
    state->ch4_flux = ch4_production * (1.0 - f_ch4_oxidized);  /* [µmol/m²/day] */
    
    /* N2O flux [nmol/m²/day]
     * From coupled nitrification-denitrification in sediments
     * Scales with NH4 flux (which drives sediment N cycling)
     */
    state->n2o_flux = state->nh4_flux * config->n2o_yield_sediment * 1e6;  /* [nmol/m²/day] */
}

/* ============================================================================
 * BRANCH UPDATE (Main integration function)
 * ============================================================================*/

int soc_update_branch(Branch *branch, double dt) {
    if (!branch || branch->M <= 0) return -1;
    
    SOCConfig *config = soc_get_config();
    if (!config || !config->enable_soc) {
        /* SOC disabled - use fixed benthic fluxes (backward compatible) */
        return 0;
    }
    
    BiogeoParams *p = rive_get_params();
    int M = branch->M;
    double temp = branch->water_temp;
    double dt_day = dt / SECONDS_PER_DAY;
    
    /* Get species arrays */
    double *o2 = branch->conc[CGEM_SPECIES_O2];
    double *nh4 = branch->conc[CGEM_SPECIES_NH4];
    double *po4 = branch->conc[CGEM_SPECIES_PO4];
    double *dic = branch->conc[CGEM_SPECIES_DIC];
    double *ch4 = branch->conc[CGEM_SPECIES_CH4];
    double *n2o = branch->conc[CGEM_SPECIES_N2O];
    double *toc = branch->conc[CGEM_SPECIES_TOC];
    double *spm = branch->conc[CGEM_SPECIES_SPM];
    
    if (!o2 || !toc) return -1;  /* Minimum required species */
    
    /* Process each grid cell */
    for (int i = 1; i <= M - 1; i += 2) {  /* Only odd indices (scalar points) */
        double depth = CGEM_MAX(branch->depth[i], CGEM_MIN_DEPTH);
        
        /* Get phytoplankton death rate from reaction_rates (if available) */
        double phy_death = 0.0;
        if (branch->reaction_rates) {
            phy_death = branch->reaction_rates[CGEM_REACTION_PHY_DEATH][i];
        }
        
        /* Calculate POC deposition */
        double deposition = soc_calc_deposition(branch, i, phy_death,
                                                 toc[i], 
                                                 spm ? spm[i] : 50.0,
                                                 branch->ws,
                                                 depth, config);
        
        /* Calculate benthic fluxes from SOC pool */
        SOCState state;
        soc_calc_fluxes(&state, branch->bed_oc[i], deposition, temp, o2[i], config);
        
        /* =====================================================================
         * UPDATE SOC POOL
         * 
         * dSOC/dt = Deposition - Mineralization - Burial
         * 
         * Explicit Euler integration (stable for small dt/k < 1)
         * =====================================================================*/
        double dSOC = (deposition - state.mineralization_rate - state.burial_rate) * dt_day;
        branch->bed_oc[i] += dSOC;
        
        /* Enforce bounds */
        if (branch->bed_oc[i] < 0.0) branch->bed_oc[i] = 0.0;
        if (branch->bed_oc[i] > config->soc_max) branch->bed_oc[i] = config->soc_max;
        
        /* =====================================================================
         * APPLY BENTHIC FLUXES TO WATER COLUMN
         * 
         * Convert from [mmol/m²/day] to [µmol/L/s]:
         *   flux [mmol/m²/day] / depth [m] / 86400 [s/day] × 1000 [µmol/mmol]
         *   = flux / depth / 86.4 [µmol/L/s]
         * =====================================================================*/
        double convert = 1000.0 / (depth * SECONDS_PER_DAY);  /* mmol/m²/day → µmol/L/s */
        
        /* O2: Consume due to SOD (negative flux from sediment perspective) */
        if (o2) {
            o2[i] -= state.sod_flux * convert * dt;
            if (o2[i] < 0.0) o2[i] = 0.0;
        }
        
        /* NH4: Release from sediment */
        if (nh4) {
            nh4[i] += state.nh4_flux * convert * dt;
        }
        
        /* PO4: Release from sediment */
        if (po4) {
            po4[i] += state.po4_flux * convert * dt;
        }
        
        /* DIC: Release from sediment */
        if (dic) {
            dic[i] += state.dic_flux * convert * dt;
        }
        
        /* CH4: Release from sediment [µmol/m²/day → µmol/L/s] */
        if (ch4) {
            double ch4_convert = 1.0 / (depth * SECONDS_PER_DAY);
            ch4[i] += state.ch4_flux * ch4_convert * dt;
        }
        
        /* N2O: Release from sediment [nmol/m²/day → nmol/L/s] */
        if (n2o) {
            double n2o_convert = 1.0 / (depth * SECONDS_PER_DAY);
            n2o[i] += state.n2o_flux * n2o_convert * dt;
        }
        
        /* =====================================================================
         * STORE RATES FOR OUTPUT
         * =====================================================================*/
        if (branch->reaction_rates) {
            /* Use the SOC-specific reaction indices if available */
            /* For now, store in benthic flux slots */
            branch->reaction_rates[CGEM_REACTION_BENTHIC_O2][i] = state.sod_flux / depth / SECONDS_PER_DAY;
            branch->reaction_rates[CGEM_REACTION_BENTHIC_DIC][i] = state.dic_flux / depth / SECONDS_PER_DAY;
            branch->reaction_rates[CGEM_REACTION_BENTHIC_NH4][i] = state.nh4_flux / depth / SECONDS_PER_DAY;
            branch->reaction_rates[CGEM_REACTION_BENTHIC_PO4][i] = state.po4_flux / depth / SECONDS_PER_DAY;
        }
    }
    
    return 0;
}

/* ============================================================================
 * UTILITY FUNCTIONS
 * ============================================================================*/

double soc_get_inventory(Branch *branch) {
    if (!branch || branch->M <= 0) return 0.0;
    
    double total = 0.0;
    int M = branch->M;
    
    for (int i = 1; i <= M; ++i) {
        /* Area = width × dx [m²] */
        double area = branch->width[i] * branch->dx;
        /* Inventory += SOC [g C/m²] × area [m²] / 1000 [kg/g] */
        total += branch->bed_oc[i] * area / 1000.0;
    }
    
    return total;  /* [kg C] */
}
