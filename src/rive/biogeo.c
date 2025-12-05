/**
 * @file biogeo.c
 * @brief RIVE Biogeochemistry Main Driver
 * 
 * Main entry point for biogeochemical calculations.
 * Coordinates phytoplankton, nutrients, oxygen, and carbonate modules.
 * 
 * Reference: Wang et al. (2018), Hasanyar et al. (2022)
 */

#include "../network.h"
#include "../define.h"
#include "rive_params.h"
#include "phytoplankton.h"
#include "nutrients.h"
#include "oxygen.h"
#include "carbonate_chem.h"
#include "ghg_module.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* Forward declaration for seasonal factor helper (defined in io_network.c) */
double GetLateralFactor(LateralSeasonalFactors *factors, int day_of_year, int species);

/* ============================================================================
 * DATA CONVERSION UTILITIES FOR SIMPLIFIED MODE
 * ============================================================================
 * These conversion factors allow C-GEM to be applied in data-sparse regions
 * where only standard water quality monitoring data (BOD5, COD) is available.
 * 
 * SIMPLIFIED MODE (simplified_mode = 1 in biogeo_params.txt):
 * - Uses single TOC pool instead of multi-pool RIVE HD1/HD2/HD3/HP1/HP2/HP3
 * - Suitable for data-sparse tropical systems with only bulk measurements
 * - Reduces calibration parameters while maintaining key biogeochemical fluxes
 * 
 * ============================================================================
 * CONVERSION FACTORS (Reference: Chapra, 2008; Tchobanoglous et al., 2014)
 * ============================================================================
 * 
 * 1. BOD5 → TOC Conversion:
 *    -----------------------
 *    BOD5 represents O2 consumed in 5 days of biodegradation.
 *    For first-order decay with rate k_bod (typically 0.23/day for domestic):
 *    
 *      BOD_ultimate = BOD5 / (1 - exp(-5 * k_bod))
 *      
 *    Using stoichiometry C:O2 = 12:32 (mineralization of glucose):
 *      TOC = BOD_ultimate * (12/32) = BOD_ultimate * 0.375
 *      
 *    For k_bod = 0.23/day:
 *      TOC ≈ BOD5 / (1 - exp(-1.15)) * 0.375 ≈ BOD5 * 0.53
 *      
 *    USAGE: TOC [mg C/L] = BOD5 [mg O2/L] * 0.53
 *    
 *    Note: This assumes mostly biodegradable organic matter. For industrial
 *    wastewater or refractory DOM, use site-specific calibration.
 * 
 * 2. COD → TOC Conversion:
 *    ----------------------
 *    COD represents total O2 demand for complete oxidation.
 *    Theoretical stoichiometry for organic carbon:
 *      CH2O + O2 → CO2 + H2O
 *      12 g C requires 32 g O2
 *      
 *    Therefore:
 *      TOC = COD / 2.67
 *      
 *    USAGE: TOC [mg C/L] = COD [mg O2/L] / 2.67
 *    
 *    Note: This assumes organic matter only. If significant inorganic
 *    reductants (sulfide, Fe2+) are present, adjust accordingly.
 * 
 * 3. Chlorophyll-a → Phytoplankton Biomass:
 *    --------------------------------------
 *    Chl-a is commonly measured; convert to carbon biomass:
 *      
 *      PHY [µg C/L] = Chl-a [µg/L] × C:Chl ratio
 *      
 *    C:Chl ratios (from literature):
 *      - Diatoms (PHY1): 30-60, use 40 (typical eutrophic)
 *      - Non-siliceous (PHY2): 40-80, use 50 (flagellates)
 *      
 *    USAGE: 
 *      PHY1 = Chl_a_diatoms * 40  (or total Chl-a * 0.6 * 40 for diatom fraction)
 *      PHY2 = Chl_a_green * 50    (or total Chl-a * 0.4 * 50 for green fraction)
 * 
 * ============================================================================
 * EXAMPLE: Converting monitoring data for Mekong Delta
 * ============================================================================
 * Monitoring station reports: BOD5 = 8 mg/L, COD = 25 mg/L, Chl-a = 15 µg/L
 * 
 * → TOC (from BOD5): 8 * 0.53 = 4.2 mg C/L  (biodegradable fraction)
 * → TOC (from COD):  25 / 2.67 = 9.4 mg C/L (total oxidizable)
 * → Use average or COD-based for tropical rivers with refractory DOM
 * → PHY1 (60% diatoms): 15 * 0.6 * 40 / 1000 = 0.36 mg C/L = 30 µM C
 * → PHY2 (40% green):   15 * 0.4 * 50 / 1000 = 0.30 mg C/L = 25 µM C
 * 
 * Reference: 
 *   - Chapra (2008) Surface Water-Quality Modeling, Ch. 22-25
 *   - Tchobanoglous et al. (2014) Wastewater Engineering, 5th Ed.
 *   - Garnier et al. (2000) RIVE model documentation
 * ============================================================================
 */

/* Conversion factor: BOD5 → TOC (assumes k_bod = 0.23/day) */
#define BOD5_TO_TOC_FACTOR  0.53

/* Conversion factor: COD → TOC (stoichiometric, C:O2 = 12:32) */
#define COD_TO_TOC_FACTOR   (1.0 / 2.67)

/* C:Chl-a ratios for phytoplankton conversion */
#define CHLA_TO_PHY1_C_RATIO  40.0  /* Diatoms */
#define CHLA_TO_PHY2_C_RATIO  50.0  /* Non-siliceous algae */

/**
 * Convert BOD5 to TOC for simplified mode input
 * @param bod5 Five-day biochemical oxygen demand [mg O2/L]
 * @return Total organic carbon [mg C/L]
 */
static inline double bod5_to_toc(double bod5) {
    return bod5 * BOD5_TO_TOC_FACTOR;
}

/**
 * Convert COD to TOC for simplified mode input
 * @param cod Chemical oxygen demand [mg O2/L]
 * @return Total organic carbon [mg C/L]
 */
static inline double cod_to_toc(double cod) {
    return cod * COD_TO_TOC_FACTOR;
}

/**
 * Convert Chlorophyll-a to phytoplankton carbon biomass
 * @param chla Chlorophyll-a concentration [µg/L]
 * @param c_to_chl C:Chl ratio (use CHLA_TO_PHY1_C_RATIO or CHLA_TO_PHY2_C_RATIO)
 * @return Phytoplankton biomass [µg C/L]
 */
static inline double chla_to_phy_carbon(double chla, double c_to_chl) {
    return chla * c_to_chl;
}

/* Global configurations */
static GHGConfig g_ghg_config = {0};
static int g_ghg_config_initialized = 0;
static CarbonateConfig g_carb_config = {0};
static int g_carb_config_initialized = 0;

/* ============================================================================
 * CARBONATE CHEMISTRY HELPERS
 * ============================================================================*/

/**
 * Henry's constant for CO2
 * 
 * UNIT CHAIN DOCUMENTATION (December 2025 Audit):
 * ================================================
 * 1. The formula (modified Weiss 1974) gives ln(K0) where K0 is in mol/(kg·atm)
 * 2. For freshwater/seawater, mol/kg ≈ mol/L (ρ ≈ 1 kg/L)
 * 3. Multiply by 1e6 to convert: mol/(L·atm) → µmol/(m³·µatm)
 *    - Numerator: 1 mol = 1e6 µmol
 *    - Denominator: 1 atm = 1e6 µatm
 *    - Volume: 1 L = 0.001 m³, so mol/L = 1000 mol/m³
 *    - Net: mol/(L·atm) × 1e6 × 1000/1e6 = 1000 × mol/(m³·µatm)
 *    - Actually: the ×1e6 here encodes both µmol conversion AND m³ scaling
 * 
 * 4. When used in calculate_carbonate_system():
 *    - Divide by 1e6 to convert m³ → L for concentration comparison
 *    - pCO2 [µatm] = CO2 [µmol/L] / (henry / 1e6) [µmol/(L·µatm)]
 * 
 * Typical values at 25°C:
 *   - Freshwater (S=0): K0 ≈ 0.034 mol/(L·atm) = 34000 µmol/(m³·µatm) after ×1e6
 *   - Seawater (S=35): K0 ≈ 0.029 mol/(L·atm) = 29000 µmol/(m³·µatm) after ×1e6
 * 
 * Reference: Weiss (1974) Marine Chemistry 2:203-215
 * 
 * @param temp Water temperature [°C]
 * @param salinity Salinity [PSU]
 * @return Henry's constant [µmol/(m³·µatm)] - divide by 1e6 for [µmol/(L·µatm)]
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

    return alk - (((*xi1 + 2.0 * (*xi2)) * dic) + ((*gamma1) * sulf) + ((*beta1) * TotB) +
                  (Kw / Hsafe) - Hsafe);
}

/**
 * Calculate carbonate system for a single grid cell
 */
static double calculate_carbonate_system(Branch *branch, int idx, double temp,
                                         double depth, double salinity, double piston_vel) {
    BiogeoParams *p = rive_get_params();
    
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
        if (!isfinite(update) || update <= 0.0) break;
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

    /* Benthic respiration */
    double benthic_rate_day = p->benthic_resp_20C * 
                              pow(p->benthic_Q10, (temp - 20.0) / 10.0);
    double benthic_co2_source = benthic_rate_day / depth_eff / RIVE_SECONDS_PER_DAY;

    /* DIC and TA reactions */
    branch->reaction_rates[CGEM_REACTION_DIC_REACT][idx] = co2_flux + aer_deg + denit 
                                                           + benthic_co2_source - npp_total;
    
    branch->reaction_rates[CGEM_REACTION_TA_REACT][idx] = (15.0 / 106.0) * aer_deg +
                                                          (93.4 / 106.0) * denit -
                                                          2.0 * nitrif +
                                                          (-15.0 / 106.0) * npp_nh4 +
                                                          (17.0 / 106.0) * npp_no3;

    /* Calculate and update pCO2 [µatm] from dissolved CO2 and Henry's constant
     * pCO2 = CO2 / K_H where K_H is Henry's constant [µmol/L/µatm]
     * 
     * Typical values for tropical estuaries:
     * - Atmosphere: ~415 µatm
     * - Freshwater (supersaturated): 1000-5000 µatm
     * - Marine (near equilibrium): 350-450 µatm
     */
    double henry_coef = henry / 1e6;  /* Convert from µmol m⁻³ µatm⁻¹ to µmol L⁻¹ µatm⁻¹ */
    if (henry_coef > 1e-12) {
        pco2[idx] = co2_spec / henry_coef;
    } else {
        pco2[idx] = 400.0;  /* Default atmospheric value */
    }
    
    /* Clamp to physical range */
    if (pco2[idx] < 100.0) pco2[idx] = 100.0;
    if (pco2[idx] > 20000.0) pco2[idx] = 20000.0;  /* Max for highly supersaturated waters */
    
    return henry;
}

/* ============================================================================
 * MAIN BIOGEOCHEMISTRY DRIVER
 * ============================================================================*/

/* ============================================================================
 * SIMPLIFIED "MEKONG MODE" (Parsimonious Biogeochemistry)
 * 
 * For data-sparse tropical systems where bacterial biomass, multi-pool OC,
 * and other research-grade parameters cannot be measured, this mode uses
 * robust engineering kinetics based on standard monitoring data (BOD, COD).
 * 
 * Key simplifications:
 * 1. Temperature-corrected first-order decay (Arrhenius) instead of bacterial dynamics
 * 2. Direct use of TOC instead of 6-pool OC fractionation
 * 3. Monod limitation on O2 only (not bacteria substrate)
 * 
 * ============================================================================*/
static int Biogeo_Branch_Simplified(Branch *branch, double dt, void *network_ptr) {
    if (!branch || branch->M <= 0) return -1;
    
    /* Get network for seasonal factors (can be NULL) */
    Network *net = (Network *)network_ptr;
    int day_of_year = net ? net->current_day : 0;
    LateralSeasonalFactors *factors = net ? &net->lateral_factors : NULL;
    
    BiogeoParams *p = rive_get_params();
    int M = branch->M;
    double temp = branch->water_temp;
    
    /* Arrhenius temperature coefficients */
    double theta_ox = (p && p->theta_ox > 0.0) ? p->theta_ox : 1.047;
    double theta_nit = (p && p->theta_nit > 0.0) ? p->theta_nit : 1.047;
    double theta_denit = (p && p->theta_denit > 0.0) ? p->theta_denit : 1.047;
    
    /* Get species arrays */
    double *salinity = branch->conc[CGEM_SPECIES_SALINITY];
    double *toc = branch->conc[CGEM_SPECIES_TOC];
    double *nh4 = branch->conc[CGEM_SPECIES_NH4];
    double *no3 = branch->conc[CGEM_SPECIES_NO3];
    double *o2 = branch->conc[CGEM_SPECIES_O2];
    double *dic = branch->conc[CGEM_SPECIES_DIC];
    double *at = branch->conc[CGEM_SPECIES_AT];
    
    if (!toc || !nh4 || !no3 || !o2) return -1;
    
    /* Pre-calculate piston velocity for gas exchange */
    double *pis_vel = (double *)calloc((size_t)(M + 2), sizeof(double));
    if (!pis_vel) return -1;
    
    for (int i = 1; i <= M; ++i) {
        pis_vel[i] = rive_piston_velocity(branch->velocity[i], branch->depth[i], temp);
    }
    
    for (int i = 1; i <= M - 1; i += 2) {
        double depth = CGEM_MAX(branch->depth[i], CGEM_MIN_DEPTH);
        double sal = salinity ? salinity[i] : 0.0;
        
        /* =======================================================================
         * TEMPERATURE CORRECTION (Arrhenius) - CRITICAL FOR TROPICS
         * Tropical waters (30°C) decay waste 2x faster than temperate (20°C)
         * =======================================================================*/
        double k_ox_T = (branch->kox / RIVE_SECONDS_PER_DAY) * pow(theta_ox, temp - 20.0);
        double k_nit_T = (branch->knit / RIVE_SECONDS_PER_DAY) * pow(theta_nit, temp - 20.0);
        double k_denit_T = (branch->kdenit / RIVE_SECONDS_PER_DAY) * pow(theta_denit, temp - 20.0);
        
        /* =======================================================================
         * AEROBIC DEGRADATION (The "BOD" equivalent)
         * Single Monod limitation on O2 - suitable for monitoring data
         * =======================================================================*/
        double o2_lim = o2[i] / (o2[i] + branch->ko2);
        double toc_deg_rate = k_ox_T * toc[i] * o2_lim;
        
        /* =======================================================================
         * NITRIFICATION (NH4 → NO3)
         * Single-step for simplicity - 2-step is in GHG module
         * =======================================================================*/
        double nh4_lim = nh4[i] / (nh4[i] + branch->knh4);
        double nit_rate = k_nit_T * nh4[i] * o2_lim * nh4_lim;
        
        /* =======================================================================
         * DENITRIFICATION (NO3 → N2) - CRITICAL FOR ALKALINITY BALANCE
         * 
         * Denitrification is the KEY pH buffer mechanism in estuaries!
         * It restores alkalinity consumed by nitrification.
         * 
         * Requirements:
         *   - Low O2 (inhibited by O2 > ~30 µM)
         *   - Organic carbon substrate (TOC)
         *   - Nitrate availability
         * 
         * Stoichiometry: 5CH2O + 4NO3- + 4H+ → 2N2 + 5CO2 + 7H2O
         * Produces 0.93 mol TA per mol N reduced (critical buffer!)
         * =======================================================================*/
        double denit_o2_inhib = branch->kino2 / (o2[i] + branch->kino2);
        double no3_lim = no3[i] / (no3[i] + branch->kno3);
        double toc_lim_denit = toc[i] / (toc[i] + branch->ktox);  /* TOC limitation */
        double denit_rate = k_denit_T * no3[i] * denit_o2_inhib * no3_lim * toc_lim_denit;
        
        /* =======================================================================
         * OXYGEN EXCHANGE
         * =======================================================================*/
        double o2_sat = rive_oxygen_saturation(temp, sal);
        double o2_ex = rive_calc_o2_exchange(o2[i], o2_sat, pis_vel[i], depth);
        
        /* =======================================================================
         * OXYGEN BALANCE
         * - Degradation (C oxidation): 1 mol O2 per mol C
         * - Nitrification: 2 mol O2 per mol N (NH4 → NO3)
         * + Reaeration
         * =======================================================================*/
        double o2_consumption = toc_deg_rate + 2.0 * nit_rate;
        
        /* Store reaction rates for output */
        branch->reaction_rates[CGEM_REACTION_AER_DEG][i] = toc_deg_rate;
        branch->reaction_rates[CGEM_REACTION_NIT][i] = nit_rate;
        branch->reaction_rates[CGEM_REACTION_DENIT][i] = denit_rate;
        branch->reaction_rates[CGEM_REACTION_O2_EX][i] = o2_ex;
        branch->reaction_rates[CGEM_REACTION_O2_EX_S][i] = o2_ex * depth;
        
        /* Zero out complex rates not used in simplified mode */
        branch->reaction_rates[CGEM_REACTION_GPP_1][i] = 0.0;
        branch->reaction_rates[CGEM_REACTION_GPP_2][i] = 0.0;
        branch->reaction_rates[CGEM_REACTION_NPP][i] = 0.0;
        
        /* =======================================================================
         * UPDATE CONCENTRATIONS
         * =======================================================================*/
        toc[i] = CGEM_MAX(0.0, toc[i] - toc_deg_rate * dt);
        nh4[i] = CGEM_MAX(0.0, nh4[i] - nit_rate * dt);
        no3[i] = CGEM_MAX(0.0, no3[i] + (nit_rate - denit_rate) * dt);
        o2[i] = CGEM_MAX(0.0, o2[i] + (o2_ex - o2_consumption) * dt);
        
        /* =======================================================================
         * CARBONATE CHEMISTRY (DIC/TA) - SCIENTIFIC MODE
         * 
         * STOICHIOMETRY (Zeebe & Wolf-Gladrow 2001, Frankignoulle et al. 1996):
         * 
         * DIC Budget:
         *   + Aerobic respiration: C6H12O6 + 6O2 → 6CO2 + 6H2O
         *   + Denitrification: also produces CO2
         *   - CO2 air-water exchange (positive = CO2 escaping to atmosphere)
         * 
         * TA Budget (H+ equivalent):
         *   + Denitrification: +0.93 TA per mol N (consumes H+ = adds TA)
         *   - Nitrification: -2.0 TA per mol N (produces 2 H+)
         *   + Aerobic degradation: +0.14 TA per mol C (minor)
         *   
         * CRITICAL: CO2 exchange does NOT directly change TA (CO2 is uncharged)
         * but indirectly affects pH through carbonate equilibrium.
         * 
         * The key to stability is BALANCING nitrification (acid) with 
         * denitrification (base). This requires:
         *   1. Sufficient TOC for denitrification substrate
         *   2. Sufficient NO3 from nitrification
         *   3. Low O2 zones where denitrification can occur
         * =======================================================================*/
        if (!p->skip_carbonate_reactions) {
            /* DIC update: respiration adds CO2, gas exchange removes it */
            if (dic) {
                /* CO2 flux: positive = entering water, negative = leaving to atm */
                double co2_flux_rate = (pis_vel[i] / depth) * 
                    (dic[i] / (1.0 + 10.0) - 420.0 * 0.035);  /* Approximate CO2 flux */
                
                /* DIC balance: respiration - CO2 evasion */
                double dic_change = toc_deg_rate + denit_rate - co2_flux_rate;
                dic[i] = CGEM_MAX(0.0, dic[i] + dic_change * dt);
            }
            
            /* TA update: scientifically correct stoichiometry
             * Reference: Wolf-Gladrow et al. (2007) Mar. Chem. 106:287-300
             * 
             * Key insight: The balance of nitrification vs denitrification
             * determines whether the system acidifies or alkalinizes.
             * 
             * With proper kox (0.05), enough TOC is degraded to sustain
             * denitrification, which restores TA and prevents pH crash.
             */
            if (at) {
                double ta_change = 0.0;
                
                /* Aerobic respiration: minor TA production
                 * (from organic N,P release) ~0.14 per mol C */
                ta_change += (15.0/106.0) * toc_deg_rate;
                
                /* Nitrification: STRONG ACID PRODUCTION
                 * NH4+ + 2O2 → NO3- + 2H+ + H2O
                 * Consumes 2 mol TA per mol N oxidized */
                ta_change -= 2.0 * nit_rate;
                
                /* Denitrification: CRITICAL ALKALINITY RESTORATION
                 * 5CH2O + 4NO3- + 4H+ → 2N2 + 5CO2 + 7H2O
                 * Produces ~0.93 mol TA per mol N reduced
                 * This is the KEY buffer mechanism in estuaries! */
                ta_change += 0.93 * denit_rate;
                
                at[i] = CGEM_MAX(0.0, at[i] + ta_change * dt);
            }
        }
        /* else: DIC and TA are transported conservatively (no reaction term) */
        
        /* =======================================================================
         * LATERAL LOADS (Land-Use Coupling) with SEASONAL FACTORS
         * 
         * Add mass flux from urban, agriculture, aquaculture, and mangroves.
         * The loads are modulated by rainfall-driven seasonal factors.
         * 
         * The source term is:
         *   dC/dt = (Q_base × Q_factor) / V * (C_base × C_factor - C)
         * 
         * where:
         *   Q_base = base lateral inflow [m³/s] (from lateral_sources.csv)
         *   Q_factor = seasonal Q multiplier (from lateral_seasonal_factors.csv)
         *   V = cell volume = Area × depth × dx
         *   C_base = base concentration [µmol/L]
         *   C_factor = seasonal concentration multiplier
         *   C = current concentration [µmol/L]
         * 
         * Seasonal factors capture:
         * - Wet season: Higher Q (runoff), but diluted concentrations
         * - Dry season: Lower Q, but concentrated pollutants
         * - First flush: Initial wet season spike in concentrations
         * 
         * Reference: Garnier et al. (2005), Billen et al. (2007), MRC (2018)
         * 
         * === UPDATED December 2025 ===
         * Added CH4, N2O, and AT lateral inputs for GHG validation
         * =======================================================================*/
        if (branch->has_lateral_loads && branch->lateral_flow && branch->lateral_conc) {
            double Q_lat_base = branch->lateral_flow[i];  /* m³/s (base/dry) */
            
            if (Q_lat_base > 1e-10) {
                /* Apply seasonal Q factor */
                double Q_factor = 1.0;
                double NH4_factor = 1.0, NO3_factor = 1.0, PO4_factor = 1.0;
                double TOC_factor = 1.0, DIC_factor = 1.0;
                /* NEW: GHG factors (December 2025) */
                double CH4_factor = 1.0, N2O_factor = 1.0, AT_factor = 1.0;
                
                if (factors && factors->loaded) {
                    Q_factor = GetLateralFactor(factors, day_of_year, -1);  /* -1 = Q */
                    NH4_factor = GetLateralFactor(factors, day_of_year, CGEM_SPECIES_NH4);
                    NO3_factor = GetLateralFactor(factors, day_of_year, CGEM_SPECIES_NO3);
                    PO4_factor = GetLateralFactor(factors, day_of_year, CGEM_SPECIES_PO4);
                    TOC_factor = GetLateralFactor(factors, day_of_year, CGEM_SPECIES_TOC);
                    DIC_factor = GetLateralFactor(factors, day_of_year, CGEM_SPECIES_DIC);
                    /* NEW: GHG factors */
                    CH4_factor = GetLateralFactor(factors, day_of_year, CGEM_SPECIES_CH4);
                    N2O_factor = GetLateralFactor(factors, day_of_year, CGEM_SPECIES_N2O);
                    AT_factor = GetLateralFactor(factors, day_of_year, CGEM_SPECIES_AT);
                }
                
                double Q_lat = Q_lat_base * Q_factor;
                
                /* Cell volume [m³] = width × depth × dx */
                double cell_volume = branch->width[i] * depth * branch->dx;
                
                if (cell_volume > 1e-6) {
                    /* Mixing factor: fraction of cell replaced per dt */
                    double mix_factor = (Q_lat * dt) / cell_volume;
                    
                    /* Clamp to prevent numerical instability */
                    if (mix_factor > 0.1) mix_factor = 0.1;  /* Max 10% replacement per step */
                    
                    /* Apply lateral concentrations with seasonal factors */
                    double C_lat_nh4 = branch->lateral_conc[CGEM_SPECIES_NH4][i] * NH4_factor;
                    double C_lat_no3 = branch->lateral_conc[CGEM_SPECIES_NO3][i] * NO3_factor;
                    double C_lat_po4 = branch->lateral_conc[CGEM_SPECIES_PO4][i] * PO4_factor;
                    double C_lat_toc = branch->lateral_conc[CGEM_SPECIES_TOC][i] * TOC_factor;
                    double C_lat_dic = branch->lateral_conc[CGEM_SPECIES_DIC][i] * DIC_factor;
                    /* NEW: GHG lateral concentrations (December 2025) */
                    double C_lat_ch4 = branch->lateral_conc[CGEM_SPECIES_CH4][i] * CH4_factor;
                    double C_lat_n2o = branch->lateral_conc[CGEM_SPECIES_N2O][i] * N2O_factor;
                    double C_lat_at = branch->lateral_conc[CGEM_SPECIES_AT][i] * AT_factor;
                    
                    /* Mix with current concentrations */
                    nh4[i] = CGEM_MAX(0.0, nh4[i] + mix_factor * (C_lat_nh4 - nh4[i]));
                    no3[i] = CGEM_MAX(0.0, no3[i] + mix_factor * (C_lat_no3 - no3[i]));
                    
                    if (branch->conc[CGEM_SPECIES_PO4]) {
                        branch->conc[CGEM_SPECIES_PO4][i] = CGEM_MAX(0.0, 
                            branch->conc[CGEM_SPECIES_PO4][i] + mix_factor * (C_lat_po4 - branch->conc[CGEM_SPECIES_PO4][i]));
                    }
                    
                    toc[i] = CGEM_MAX(0.0, toc[i] + mix_factor * (C_lat_toc - toc[i]));
                    
                    /* DIC and AT lateral inputs - skip if skip_carbonate_reactions=1
                     * This ensures pCO2/pH remain transport-dominated when this flag is set */
                    if (!p->skip_carbonate_reactions) {
                        if (dic) {
                            dic[i] = CGEM_MAX(0.0, dic[i] + mix_factor * (C_lat_dic - dic[i]));
                        }
                        
                        /* AT (Total Alkalinity): Mix from weathering and fertilizer inputs */
                        if (branch->conc[CGEM_SPECIES_AT] && C_lat_at > 0) {
                            branch->conc[CGEM_SPECIES_AT][i] = CGEM_MAX(0.0,
                                branch->conc[CGEM_SPECIES_AT][i] + mix_factor * (C_lat_at - branch->conc[CGEM_SPECIES_AT][i]));
                        }
                    }
                    
                    /* === NEW: Mix GHG species (December 2025) === */
                    /* CH4: Mix from rice paddies/aquaculture drainage */
                    if (branch->conc[CGEM_SPECIES_CH4] && C_lat_ch4 > 0) {
                        branch->conc[CGEM_SPECIES_CH4][i] = CGEM_MAX(0.0,
                            branch->conc[CGEM_SPECIES_CH4][i] + mix_factor * (C_lat_ch4 - branch->conc[CGEM_SPECIES_CH4][i]));
                    }
                    
                    /* N2O: Mix from agricultural nitrification/denitrification */
                    if (branch->conc[CGEM_SPECIES_N2O] && C_lat_n2o > 0) {
                        branch->conc[CGEM_SPECIES_N2O][i] = CGEM_MAX(0.0,
                            branch->conc[CGEM_SPECIES_N2O][i] + mix_factor * (C_lat_n2o - branch->conc[CGEM_SPECIES_N2O][i]));
                    }
                }
            }
        }
        
        /* =======================================================================
         * CARBONATE CHEMISTRY (compute pH, pCO2 from DIC/TA)
         * Skip if skip_carbonate_reactions=1 - let transport handle pCO2/pH
         * =======================================================================*/
        if (!p->skip_carbonate_reactions) {
            calculate_carbonate_system(branch, i, temp, depth, sal, pis_vel[i]);
        }
    }
    
    free(pis_vel);
    return 0;
}

/**
 * Main biogeochemistry function for a branch
 * Calculates all reaction rates and updates species concentrations.
 * 
 * @param branch Branch to process
 * @param dt Time step [s]
 * @param network_ptr Pointer to Network (for seasonal factors)
 */
int Biogeo_Branch(Branch *branch, double dt, void *network_ptr) {
    if (!branch || branch->M <= 0 || branch->num_species < CGEM_NUM_SPECIES || !branch->reaction_rates) {
        return 0;
    }
    
    /* Get network for seasonal factors (can be NULL for backward compatibility) */
    Network *net = (Network *)network_ptr;
    
    /* =======================================================================
     * SIMPLIFIED MODE CHECK
     * 
     * If simplified_mode = 1, bypass complex RIVE bacterial/multi-pool OC
     * dynamics and use robust engineering kinetics suitable for data-sparse
     * tropical systems like the Mekong Delta.
     * 
     * Reference: Audit recommendation for "80/20 Mode"
     * =======================================================================*/
    BiogeoParams *p = rive_get_params();
    if (p && p->simplified_mode) {
        return Biogeo_Branch_Simplified(branch, dt, network_ptr);
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
    double *dic = branch->conc[CGEM_SPECIES_DIC];
    double *at = branch->conc[CGEM_SPECIES_AT];
    double *co2 = branch->conc[CGEM_SPECIES_CO2];
    double *pco2 = branch->conc[CGEM_SPECIES_PCO2];

    if (!salinity || !phy1 || !phy2 || !dsi || !no3 || !nh4 || !po4 || !o2 || !toc ||
        !dic || !at || !co2 || !pco2) {
        return 0;
    }

    /* Allocate piston velocity array */
    double *pis_vel = (double *)calloc((size_t)(M + 2), sizeof(double));
    if (!pis_vel) return -1;
    
    for (int i = 1; i <= M; ++i) {
        pis_vel[i] = rive_piston_velocity(branch->velocity[i], branch->depth[i], temp);
    }

    /* ========================================================================
     * LOOP 1: Calculate all reaction rates
     * ========================================================================*/
    for (int i = 1; i <= M - 1; i += 2) {
        double depth = CGEM_MAX(branch->depth[i], CGEM_MIN_DEPTH);
        double sal = salinity[i];

        /* Phytoplankton */
        double gpp1, gpp2, npp_no31, npp_no32, npp_nh41, npp_nh42, mort1, mort2;
        rive_calc_phytoplankton(branch, i, temp, depth,
                                &gpp1, &gpp2,
                                &npp_no31, &npp_no32,
                                &npp_nh41, &npp_nh42,
                                &mort1, &mort2);

        /* Silica consumption */
        double si_cons = rive_calc_silica_consumption(npp_nh41 + npp_nh42, branch->redsi);

        /* Nutrient transformations */
        double aer_deg = rive_calc_aerobic_degradation(toc[i], o2[i], temp,
                                                       branch->kox, branch->ktox, branch->ko2);
        
        double denit = rive_calc_denitrification(no3[i], o2[i], toc[i], temp,
                                                 branch->kdenit, branch->kino2, 
                                                 branch->kno3, branch->ktox);
        
        double nitrif = rive_calc_nitrification(nh4[i], o2[i], temp,
                                                branch->knit, branch->ko2_nit, branch->knh4);

        /* Oxygen exchange */
        double o2_sat = rive_oxygen_saturation(temp, sal);
        double o2_ex = rive_calc_o2_exchange(o2[i], o2_sat, pis_vel[i], depth);

        /* Aggregate rates */
        double npp_no3_sum = npp_no31 + npp_no32;
        double npp_nh4_sum = npp_nh41 + npp_nh42;
        double npp_total = npp_no3_sum + npp_nh4_sum;
        double phy_death = mort1 + mort2;

        /* Store reaction rates */
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

    /* ========================================================================
     * LOOP 2: Carbonate system and state updates
     * ========================================================================*/
    for (int i = 1; i <= M - 1; i += 2) {
        double depth = CGEM_MAX(branch->depth[i], CGEM_MIN_DEPTH);
        double sal = salinity[i];

        double henry = calculate_carbonate_system(branch, i, temp, depth, sal, pis_vel[i]);

        /* Get rates */
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

        /* Update state variables */
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

/**
 * Calculate GHG dynamics for a branch (N2O, CH4)
 * 
 * PASSIVE MODE (ghg_passive_mode = 1, DEFAULT):
 *   - Reads nitrification/denitrification rates from core module (Biogeo_Branch)
 *   - Calculates N2O as yield fraction of these rates
 *   - CH4 dynamics are independent (benthic flux + oxidation + exchange)
 *   - Does NOT modify O2, NH4, NO3 (prevents double-counting)
 * 
 * ACTIVE MODE (ghg_passive_mode = 0, USE WITH CAUTION):
 *   - Recalculates 2-step nitrification (NH4 → NO2 → NO3)
 *   - Applies feedback to O2, NH4, NO3
 *   - Only use if you have calibrated GHG parameters
 */
int Biogeo_GHG_Branch(Branch *branch, double dt) {
    if (!branch || branch->M <= 0) return -1;
    
    BiogeoParams *p = rive_get_params();
    
    /* Initialize GHG config if needed */
    if (!g_ghg_config_initialized) {
        crive_ghg_init_config(&g_ghg_config);
        /* Override defaults with values from biogeo_params.txt */
        g_ghg_config.benthic_CH4_flux = p->benthic_CH4_flux;
        g_ghg_config.benthic_N2O_flux = p->benthic_N2O_flux;
        g_ghg_config_initialized = 1;
    }
    
    /* Skip GHG if disabled */
    if (p->skip_ghg_dynamics) return 0;
    
    /* Check if GHG species arrays exist */
    if (branch->num_species < CGEM_SPECIES_CH4 + 1) return 0;
    
    double *no2 = branch->conc[CGEM_SPECIES_NO2];
    double *n2o = branch->conc[CGEM_SPECIES_N2O];
    double *ch4 = branch->conc[CGEM_SPECIES_CH4];
    
    if (!no2 || !n2o || !ch4) return 0;
    
    int M = branch->M;
    double temp = branch->water_temp;
    
    double *salinity = branch->conc[CGEM_SPECIES_SALINITY];
    double *o2 = branch->conc[CGEM_SPECIES_O2];
    double *nh4 = branch->conc[CGEM_SPECIES_NH4];
    double *no3 = branch->conc[CGEM_SPECIES_NO3];
    double *toc = branch->conc[CGEM_SPECIES_TOC];
    
    for (int i = 1; i <= M - 1; i += 2) {
        double depth = CGEM_MAX(branch->depth[i], CGEM_MIN_DEPTH);
        double velocity = fabs(branch->velocity[i]);
        double sal = salinity ? salinity[i] : 0.0;
        
        /* =====================================================================
         * GET CORE REACTION RATES (The "Listener" Approach)
         * These were calculated by Biogeo_Branch - we use them, not recalculate
         * ===================================================================== */
        double core_nit_rate = 0.0;
        double core_denit_rate = 0.0;
        if (branch->reaction_rates) {
            core_nit_rate = branch->reaction_rates[CGEM_REACTION_NIT][i];
            core_denit_rate = branch->reaction_rates[CGEM_REACTION_DENIT][i];
        }
        
        /* =====================================================================
         * N2O CALCULATION (Yield-based from core rates)
         * N2O is just a fraction of the N processed by nitrification/denitrification
         * Reference: Garnier et al. (2007), Marescaux et al. (2019)
         * 
         * DECEMBER 2025 FIX: All rates now in µmol/L/s to match model internal units.
         * The ×1000 factor was removed from ghg_module.c to fix unit mismatch.
         * ===================================================================== */
        /* N2O from nitrification: yield × nit_rate [µmol N/L/s] */
        double n2o_from_nit = core_nit_rate * p->N2O_yield_nit;
        
        /* N2O from denitrification: yield × denit_rate [µmol N/L/s] */
        double n2o_from_denit = core_denit_rate * p->N2O_yield_denit;
        
        /* N2O air-water exchange */
        double k600 = crive_calc_k600(velocity, depth, K600_STRAHLER, 6, 0.0);
        double n2o_sat = crive_calc_n2o_sat(temp, sal);
        double n2o_flux = crive_calc_n2o_flux(n2o[i], n2o_sat, depth, k600, temp);
        
        /* =====================================================================
         * CH4 CALCULATION (Independent, does not affect O2 budget)
         * CH4 is primarily from benthic sources in estuaries
         * ===================================================================== */
        GHGState ghg_state;
        ghg_state.NO2 = no2[i];
        ghg_state.N2O = n2o[i];
        ghg_state.CH4 = ch4[i];
        
        crive_calc_ghg_system(&ghg_state, temp, sal, depth, velocity,
                              o2 ? o2[i] : 200.0,
                              no3 ? no3[i] : 10.0,
                              nh4 ? nh4[i] : 5.0,
                              toc ? toc[i] : 100.0,
                              core_denit_rate, k600, &g_ghg_config);
        
        /* NOTE: benthic_CH4_flux is already included in crive_calc_ghg_system()
         * via crive_calc_ch4_production(). Do NOT add it again here!
         */
        
        /* =====================================================================
         * STORE REACTION RATES (for output/diagnostics)
         * ===================================================================== */
        if (branch->reaction_rates) {
            branch->reaction_rates[CGEM_REACTION_N2O_NIT][i] = n2o_from_nit;
            branch->reaction_rates[CGEM_REACTION_N2O_DENIT][i] = n2o_from_denit;
            branch->reaction_rates[CGEM_REACTION_N2O_EX][i] = n2o_flux;
            branch->reaction_rates[CGEM_REACTION_CH4_PROD][i] = ghg_state.CH4_prod;
            branch->reaction_rates[CGEM_REACTION_CH4_OX][i] = ghg_state.CH4_ox;
            branch->reaction_rates[CGEM_REACTION_CH4_EX][i] = ghg_state.CH4_flux;
            branch->reaction_rates[CGEM_REACTION_CH4_EBUL][i] = ghg_state.CH4_ebul;
        }
        
        /* =====================================================================
         * UPDATE GHG CONCENTRATIONS
         * ===================================================================== */
        if (p->ghg_passive_mode == 1) {
            /* =================================================================
             * PASSIVE MODE (SAFE): Only update N2O/CH4, DO NOT touch O2/NH4/NO3
             * This prevents double-counting with the core biogeochemistry module
             * ================================================================= */
            
            /* Update N2O: production from yields - air-water exchange */
            double dN2O = n2o_from_nit + n2o_from_denit - n2o_flux;
            n2o[i] = CGEM_MAX(0.0, n2o[i] + dN2O * dt);
            
            /* Update CH4: benthic production - oxidation - air-water exchange - ebullition */
            double dCH4 = ghg_state.CH4_prod - ghg_state.CH4_ox - ghg_state.CH4_flux - ghg_state.CH4_ebul;
            ch4[i] = CGEM_MAX(0.0, ch4[i] + dCH4 * dt);
            
            /* NO2 in passive mode: just track intermediate (no feedback) */
            /* NO2 dynamics are diagnostic only - they don't affect NO3 budget */
            
        } else {
            /* =================================================================
             * ACTIVE MODE (RISKY): Apply full 2-step nitrification feedback
             * WARNING: Only use if you have calibrated GHG parameters!
             * ================================================================= */
            
            /* 2-step nitrification */
            double nitrosation = crive_calc_nitrosation(nh4 ? nh4[i] : 0.0, 
                                                         o2 ? o2[i] : 200.0, 
                                                         temp, &g_ghg_config);
            double nitratation = crive_calc_nitratation(no2[i], 
                                                         o2 ? o2[i] : 200.0, 
                                                         temp, &g_ghg_config);
            
            if (branch->reaction_rates) {
                branch->reaction_rates[CGEM_REACTION_NIT][i] = nitrosation;
                branch->reaction_rates[CGEM_REACTION_NIT2][i] = nitratation;
            }
            
            /* Update NO2 (intermediate) */
            double dNO2 = nitrosation - nitratation - n2o_from_nit * 0.001;
            no2[i] = CGEM_MAX(0.0, no2[i] + dNO2 * dt);
            
            /* Update N2O */
            double dN2O = n2o_from_nit + n2o_from_denit - n2o_flux;
            n2o[i] = CGEM_MAX(0.0, n2o[i] + dN2O * dt);
            
            /* Update CH4 */
            double dCH4 = ghg_state.CH4_prod - ghg_state.CH4_ox - ghg_state.CH4_ox_anaer 
                        - ghg_state.CH4_flux - ghg_state.CH4_ebul;
            ch4[i] = CGEM_MAX(0.0, ch4[i] + dCH4 * dt);
            
            /* ACTIVE MODE: Apply feedback to core variables (RISKY!) */
            if (nh4) nh4[i] = CGEM_MAX(0.0, nh4[i] - nitrosation * dt);
            if (no3) no3[i] = CGEM_MAX(0.0, no3[i] + nitratation * dt);
            if (o2) {
                /* O2 consumption: 1.5 mol O2 per mol NH4 → NO2, 0.5 mol O2 per mol NO2 → NO3 */
                double o2_nit_consumption = 1.5 * nitrosation + 0.5 * nitratation;
                /* CH4 oxidation also consumes O2: 2 mol O2 per mol CH4 */
                double o2_ch4_consumption = 2.0 * ghg_state.CH4_ox;
                o2[i] = CGEM_MAX(0.0, o2[i] - (o2_nit_consumption + o2_ch4_consumption) * dt);
            }
        }
    }
    
    return 0;
}
