/**
 * @file ghg_module.c
 * @brief C-RIVE Greenhouse Gas Module Implementation
 * 
 * Complete implementation of CO2, N2O, and CH4 dynamics for CGEM.
 * Ported from C-RIVE rea_degassing.c, growth_bactn.c (Unified RIVE v1.0)
 * 
 * CITATION:
 * Wang, S., Flipo, N., Romary, T., 2018. Time-dependent global sensitivity analysis
 * of the C-RIVE biogeochemical model. Water Research 144, 341-355.
 * 
 * COPYRIGHT: Original C-RIVE code (c) 2023 Contributors to the librive library.
 * Eclipse Public License v2.0
 */

#include "ghg_module.h"
#include <stdio.h>
#include <stdlib.h>

/* ===========================================================================
 * Local Helper Macros
 * ===========================================================================*/

#define CRIVE_POLY2(a, b, c, x) ((a) + (b)*(x) + (c)*(x)*(x))
#define CRIVE_POLY3(a, b, c, d, x) ((a) + (b)*(x) + (c)*(x)*(x) + (d)*(x)*(x)*(x))

/* Time conversions */
#define SECONDS_PER_DAY   86400.0
#define SECONDS_PER_HOUR  3600.0

/* Small value to prevent division by zero */
#define GHG_EPS 1e-10

/* ===========================================================================
 * Configuration Initialization
 * ===========================================================================*/

void crive_ghg_init_config(GHGConfig *config) {
    if (!config) return;
    
    /* N2O parameters (from C-RIVE and Garnier et al. 2007) */
    config->N2O_atm_ppb = CRIVE_N2O_ATM_PPB;
    config->N2O_yield_nit = CRIVE_N2O_YIELD_NIT;
    config->N2O_yield_denit = CRIVE_N2O_YIELD_DENIT;
    config->kmich_no2_n2o = 5.0;  /* µmol/L */
    
    /* CH4 parameters (from Borges & Abril 2011, tropical estuaries) */
    config->CH4_atm_ppb = CRIVE_CH4_ATM_PPB;
    config->CH4_prod_rate = CRIVE_CH4_PROD_RATE;
    config->CH4_ox_rate = CRIVE_CH4_OX_RATE;
    config->CH4_ks_o2 = CRIVE_CH4_KS_O2;
    config->CH4_ks_so4 = CRIVE_CH4_KS_SO4;
    config->CH4_ebul_thresh = CRIVE_CH4_EBUL_THRESH;
    config->CH4_ebul_rate = 1.0;  /* /day */
    
    /* 2-step nitrification parameters (from C-RIVE growth_bactn.c) */
    config->k_nit1 = 0.3;   /* /day at 20°C - nitrosation */
    config->k_nit2 = 0.6;   /* /day at 20°C - nitratation (faster than nitrosation) */
    config->ks_nh4 = 2.0;   /* µmol/L - NH4 half-saturation */
    config->ks_no2 = 1.0;   /* µmol/L - NO2 half-saturation */
    config->ks_o2_nit = 5.0; /* µmol/L - O2 half-saturation for nitrification */
    
    /* Temperature coefficients */
    config->Q10_nit = 1.08;      /* Q10 for nitrification */
    config->Q10_ch4_prod = 2.5;  /* Q10 for methanogenesis (high for anaerobic) */
    config->Q10_ch4_ox = 2.0;    /* Q10 for CH4 oxidation */
    
    /* Benthic fluxes (tropical delta defaults) */
    config->benthic_CH4_flux = 100.0;  /* µmol/m²/day */
    config->benthic_N2O_flux = 5.0;    /* nmol/m²/day */
}

/* ===========================================================================
 * N2O Functions - Direct port from C-RIVE
 * ===========================================================================*/

double crive_calc_n2o_sat(double temp, double salinity) {
    /* Weiss & Price (1980) formulation for N2O saturation
     * From C-RIVE calc_N2O_sat() and init_N2O()
     * 
     * Coefficients for freshwater (salinity effect minimal):
     * Csat = (a0 + a1*T + a2*T²) * 14 * 1e-6 [ngN/L]
     * 
     * C-RIVE uses: sat[0]=9.55, sat[1]=-0.38, sat[2]=0.00435
     */
    const double a0 = 9.55;
    const double a1 = -0.38;
    const double a2 = 0.00435;
    
    double Csat = CRIVE_POLY2(a0, a1, a2, temp) * 14.0 * 1e-6;  /* ngN/L */
    
    /* Convert ngN/L to nmol/L
     * 1 nmol N = 14 ng N
     * So ngN/L / 14 = nmol/L
     * But C-RIVE already multiplies by 14, so result is in ~nmol equivalent
     * Actually: Csat in nmol/L = CRIVE value (which is in scaled units)
     */
    
    /* For a cleaner implementation, use standard Weiss & Price:
     * ln(Csat) = A1 + A2*(100/T) + A3*ln(T/100) + S*(B1 + B2*(T/100))
     * with atmospheric N2O = 335 ppb
     */
    double T_K = temp + CRIVE_T_KELVIN_GHG;
    double T100 = T_K / 100.0;
    
    /* Weiss & Price (1980) coefficients for N2O */
    const double A1_wp = -165.8806;
    const double A2_wp = 222.8743;
    const double A3_wp = 92.0792;
    const double A4_wp = -1.48425;
    const double B1_wp = -0.056235;
    const double B2_wp = 0.031619;
    const double B3_wp = -0.0048472;
    
    double ln_beta = A1_wp + A2_wp * (100.0 / T_K) + A3_wp * log(T100) + A4_wp * T100 * T100;
    ln_beta += salinity * (B1_wp + B2_wp * T100 + B3_wp * T100 * T100);
    
    double beta = exp(ln_beta);  /* mol/L/atm */
    
    /* N2O saturation at atmospheric equilibrium [µmol/L]
     * DECEMBER 2025 FIX: Return µmol/L to match model internal units.
     * Previously returned nmol/L which caused unit mismatch in flux calculation.
     */
    double N2O_atm_atm = CRIVE_N2O_ATM_PPB * 1e-9;  /* ppb to atm */
    Csat = beta * N2O_atm_atm * 1e6;  /* mol/L/atm * atm * 1e6 = µmol/L */
    
    return Csat;
}

double crive_calc_n2o_schmidt(double temp) {
    /* Schmidt number for N2O from C-RIVE
     * Sc = a0 + a1*T + a2*T² + a3*T³
     * 
     * From Wanninkhof (1992) / C-RIVE init_N2O():
     * sc[0]=2055.4, sc[1]=-137.11, sc[2]=4.3173, sc[3]=-0.054350
     */
    const double a0 = 2055.4;
    const double a1 = -137.11;
    const double a2 = 4.3173;
    const double a3 = -0.054350;
    
    return CRIVE_POLY3(a0, a1, a2, a3, temp);
}

double crive_calc_n2o_flux(double N2O_conc, double N2O_sat, double depth,
                           double velocity, double Sc) {
    /* Direct port from C-RIVE rea_degassing_N2O()
     * 
     * kg = 1.719 * sqrt((600 * v) / (Sc * h))
     * Fwa = kg * (C - Csat)
     * 
     * Note: Flux is positive when N2O evades (C > Csat)
     */
    if (depth < 0.01) depth = 0.01;
    if (Sc < 100.0) Sc = 100.0;
    
    double kg = 1.719 * sqrt((600.0 * fabs(velocity)) / (Sc * depth));
    double Fwa = kg * (N2O_conc - N2O_sat);
    
    /* DECEMBER 2025 FIX: Return µmol/L/s to match model internal units.
     * N2O_conc and N2O_sat are now both in µmol/L.
     */
    return Fwa;  /* [µmol/L/s] - positive = evasion */
}

double crive_calc_n2o_from_nitrification(double nitrif_rate, double yield) {
    /* N2O production from nitrification
     * From C-RIVE growth_bactn.c and Cébron et al. (2005)
     * 
     * N2O_nit = yield * nitrification_rate
     * 
     * Typical yield: 0.003-0.01 mol N2O-N / mol NH4-N oxidized
     * 
     * DECEMBER 2025 FIX: Keep in µmol/L/s to match model's internal units.
     * The ×1000 conversion was causing unit mismatch with n2o[] array.
     */
    if (nitrif_rate < 0.0) nitrif_rate = 0.0;
    
    return yield * nitrif_rate;  /* µmol N/L/s - matches model internal unit */
}

double crive_calc_n2o_from_denitrification(double denit_rate, double yield) {
    /* N2O production from incomplete denitrification
     * From C-RIVE and Garnier et al. (2007)
     * 
     * N2O_denit = yield * denitrification_rate
     * 
     * Typical yield: 0.005-0.02 mol N2O / mol NO3 reduced
     * 
     * DECEMBER 2025 FIX: Keep in µmol/L/s to match model's internal units.
     */
    if (denit_rate < 0.0) denit_rate = 0.0;
    
    return yield * denit_rate;  /* µmol N/L/s - matches model internal unit */
}

/* ===========================================================================
 * CH4 Functions
 * Based on Borges & Abril (2011), adapted for CGEM
 * ===========================================================================*/

double crive_calc_ch4_sat(double temp, double salinity, double CH4_atm) {
    /* CH4 saturation using Yamamoto et al. (1976) solubility
     * 
     * ln(beta) = A1 + A2*(100/T) + A3*ln(T/100) + S*(B1 + B2*(T/100) + B3*(T/100)²)
     * 
     * beta in mol/L/atm
     */
    double T_K = temp + CRIVE_T_KELVIN_GHG;
    double T100 = T_K / 100.0;
    
    /* Yamamoto et al. (1976) coefficients for CH4 */
    const double A1 = -415.2807;
    const double A2 = 596.8104;
    const double A3 = 379.2599;
    const double A4 = -62.0757;
    const double B1 = -0.059160;
    const double B2 = 0.032174;
    const double B3 = -0.0048198;
    
    double ln_beta = A1 + A2 * (100.0 / T_K) + A3 * log(T100) + A4 * T100;
    ln_beta += salinity * (B1 + B2 * T100 + B3 * T100 * T100);
    
    double beta = exp(ln_beta);  /* mol/L/atm */
    
    /* CH4 saturation at atmospheric equilibrium [µmol/L] */
    double CH4_atm_atm = CH4_atm * 1e-9;  /* ppb to atm */
    double Csat = beta * CH4_atm_atm * 1e6;  /* mol/L/atm * atm * 1e6 = µmol/L */
    
    return Csat;
}

double crive_calc_ch4_schmidt(double temp) {
    /* Schmidt number for CH4
     * Wanninkhof (2014) coefficients
     * Sc = a0 - a1*T + a2*T² - a3*T³
     */
    const double a0 = 1909.4;
    const double a1 = 120.78;
    const double a2 = 4.1555;
    const double a3 = 0.052569;
    
    return a0 - a1*temp + a2*temp*temp - a3*temp*temp*temp;
}

double crive_calc_ch4_flux(double CH4_conc, double CH4_sat, double depth,
                           double k600, double Sc) {
    /* CH4 air-water flux similar to CO2
     * 
     * Flux = kg * (C - Csat) / depth
     * where kg = k600 * sqrt(600/Sc)
     * 
     * Positive = evasion
     */
    if (depth < 0.01) depth = 0.01;
    if (Sc < 100.0) Sc = 100.0;
    
    double kg = k600 * sqrt(600.0 / Sc);
    double Fwa = kg * (CH4_conc - CH4_sat) / depth;
    
    return Fwa;  /* [µmol/L/s] - positive = evasion */
}

double crive_calc_ch4_production(double TOC_benthic, double depth, double temp,
                                  double O2_conc, const GHGConfig *config) {
    /* Methanogenesis in anoxic sediments
     * 
     * Production occurs when O2 is depleted, scaling with:
     * - Benthic organic carbon availability
     * - Temperature (high Q10 for anaerobic processes)
     * - Inverse of O2 concentration (inhibition)
     * 
     * Returns production as flux into water column [µmol/L/s]
     */
    if (!config) return 0.0;
    if (depth < 0.01) depth = 0.01;
    
    /* O2 inhibition (methanogenesis only under very low O2) */
    double o2_inhib = 1.0 / (1.0 + O2_conc / 10.0);  /* Strong inhibition by O2 */
    
    /* Temperature dependence */
    double temp_factor = pow(config->Q10_ch4_prod, (temp - 20.0) / 10.0);
    
    /* Base production rate from benthic OC [µmol CH4/m²/day] */
    double prod_areal = config->CH4_prod_rate * TOC_benthic * temp_factor * o2_inhib;
    
    /* Convert to volumetric rate [µmol/L/s] */
    /* prod_areal [µmol/m²/day] / depth [m] / 86400 [s/day] * 1000 [L/m³] */
    double prod_vol = prod_areal / depth / SECONDS_PER_DAY;
    
    /* Add sediment flux if configured */
    prod_vol += config->benthic_CH4_flux / depth / SECONDS_PER_DAY * temp_factor;
    
    return prod_vol;
}

double crive_calc_ch4_oxidation(double CH4_conc, double O2_conc, double temp,
                                 const GHGConfig *config) {
    /* Aerobic CH4 oxidation (methanotrophy)
     * 
     * CH4 + 2O2 → CO2 + 2H2O
     * 
     * Michaelis-Menten kinetics with both CH4 and O2 limitation
     */
    if (!config) return 0.0;
    if (CH4_conc <= 0.0 || O2_conc <= 0.0) return 0.0;
    
    /* Temperature dependence */
    double temp_factor = pow(config->Q10_ch4_ox, (temp - 20.0) / 10.0);
    
    /* Michaelis-Menten for O2 */
    double O2_lim = O2_conc / (O2_conc + config->CH4_ks_o2);
    
    /* First-order in CH4 with O2 limitation */
    double k_ox = config->CH4_ox_rate / SECONDS_PER_DAY;  /* Convert to /s */
    double ox_rate = k_ox * CH4_conc * O2_lim * temp_factor;
    
    return ox_rate;  /* [µmol/L/s] */
}

double crive_calc_ch4_oxidation_anaerobic(double CH4_conc, double SO4_conc, 
                                           double temp, const GHGConfig *config) {
    /* Anaerobic CH4 oxidation (sulfate-mediated)
     * Important in marine/estuarine sediments
     * 
     * CH4 + SO4²⁻ → HCO3⁻ + HS⁻ + H2O
     * 
     * Much slower than aerobic oxidation
     */
    if (!config) return 0.0;
    if (CH4_conc <= 0.0 || SO4_conc <= 0.0) return 0.0;
    
    /* Temperature dependence */
    double temp_factor = pow(config->Q10_ch4_ox, (temp - 20.0) / 10.0);
    
    /* Michaelis-Menten for SO4 */
    double SO4_lim = SO4_conc / (SO4_conc + config->CH4_ks_so4);
    
    /* Anaerobic oxidation ~10x slower than aerobic */
    double k_ox_anaer = config->CH4_ox_rate * 0.1 / SECONDS_PER_DAY;
    double ox_rate = k_ox_anaer * CH4_conc * SO4_lim * temp_factor;
    
    return ox_rate;  /* [µmol/L/s] */
}

double crive_calc_ch4_ebullition(double CH4_conc, double CH4_sat, double depth,
                                  const GHGConfig *config) {
    /* CH4 ebullition (bubble release)
     * 
     * Triggered when CH4 exceeds a threshold relative to saturation
     * From DelSontro et al. (2018)
     * 
     * Ebullition bypasses the water column and releases CH4 directly to atmosphere
     */
    if (!config) return 0.0;
    if (CH4_conc <= config->CH4_ebul_thresh) return 0.0;
    
    /* Ebullition rate proportional to supersaturation */
    double excess = CH4_conc - config->CH4_ebul_thresh;
    double k_ebul = config->CH4_ebul_rate / SECONDS_PER_DAY;
    
    /* Inverse depth dependence (shallow water = more ebullition) */
    double depth_factor = 1.0 / (depth + 0.1);
    
    return k_ebul * excess * depth_factor;  /* [µmol/L/s] */
}

/* ===========================================================================
 * 2-Step Nitrification - From C-RIVE growth_bactn.c
 * ===========================================================================*/

double crive_calc_nitrosation(double NH4_conc, double O2_conc, double temp,
                               const GHGConfig *config) {
    /* Nitrosation: NH4 → NO2 (first step of nitrification)
     * 
     * From C-RIVE growth_nitrosantes() / growth_bactn_1esp()
     * 
     * Rate = k_nit1 * f(T) * mich(NH4) * mich(O2)
     */
    if (!config) return 0.0;
    if (NH4_conc <= 0.0 || O2_conc <= 0.0) return 0.0;
    
    /* Temperature factor (Q10) */
    double temp_factor = pow(config->Q10_nit, (temp - 20.0) / 10.0);
    
    /* Michaelis-Menten limitations */
    double NH4_lim = NH4_conc / (NH4_conc + config->ks_nh4);
    double O2_lim = O2_conc / (O2_conc + config->ks_o2_nit);
    
    /* Rate [/day] converted to [/s] */
    double k1 = config->k_nit1 / SECONDS_PER_DAY;
    
    return k1 * NH4_conc * NH4_lim * O2_lim * temp_factor;  /* [µmol N/L/s] */
}

double crive_calc_nitratation(double NO2_conc, double O2_conc, double temp,
                               const GHGConfig *config) {
    /* Nitratation: NO2 → NO3 (second step of nitrification)
     * 
     * From C-RIVE growth_nitratantes()
     * 
     * Rate = k_nit2 * f(T) * mich(NO2) * mich(O2)
     * 
     * Typically faster than nitrosation, preventing NO2 accumulation
     */
    if (!config) return 0.0;
    if (NO2_conc <= 0.0 || O2_conc <= 0.0) return 0.0;
    
    /* Temperature factor (Q10) */
    double temp_factor = pow(config->Q10_nit, (temp - 20.0) / 10.0);
    
    /* Michaelis-Menten limitations */
    double NO2_lim = NO2_conc / (NO2_conc + config->ks_no2);
    double O2_lim = O2_conc / (O2_conc + config->ks_o2_nit);
    
    /* Rate [/day] converted to [/s] */
    double k2 = config->k_nit2 / SECONDS_PER_DAY;
    
    return k2 * NO2_conc * NO2_lim * O2_lim * temp_factor;  /* [µmol N/L/s] */
}

/* ===========================================================================
 * Integrated GHG System Calculation
 * ===========================================================================*/

int crive_calc_ghg_system(GHGState *state, double temp, double salinity,
                          double depth, double velocity, double O2_conc,
                          double NO3_conc, double NH4_conc, double TOC_conc,
                          double denit_rate, double k600, const GHGConfig *config) {
    if (!state || !config) return -1;
    
    /* Ensure minimum depth */
    if (depth < 0.01) depth = 0.01;
    
    /* =======================================================================
     * N2O Calculations
     * =======================================================================*/
    
    /* N2O saturation and Schmidt number */
    state->N2O_sat = crive_calc_n2o_sat(temp, salinity);
    state->Sc_N2O = crive_calc_n2o_schmidt(temp);
    
    /* 2-step nitrification */
    double nitrosation = crive_calc_nitrosation(NH4_conc, O2_conc, temp, config);
    double nitratation = crive_calc_nitratation(state->NO2, O2_conc, temp, config);
    
    /* N2O production from nitrification (intermediate NO2 pathway) */
    state->N2O_from_nit = crive_calc_n2o_from_nitrification(nitrosation, config->N2O_yield_nit);
    
    /* N2O production from denitrification */
    state->N2O_from_denit = crive_calc_n2o_from_denitrification(denit_rate, config->N2O_yield_denit);
    
    /* N2O air-water flux */
    state->N2O_flux = crive_calc_n2o_flux(state->N2O, state->N2O_sat, depth, 
                                           velocity, state->Sc_N2O);
    
    /* =======================================================================
     * CH4 Calculations
     * =======================================================================*/
    
    /* CH4 saturation and Schmidt number */
    state->CH4_sat = crive_calc_ch4_sat(temp, salinity, config->CH4_atm_ppb);
    state->Sc_CH4 = crive_calc_ch4_schmidt(temp);
    
    /* CH4 production (benthic methanogenesis) */
    state->CH4_prod = crive_calc_ch4_production(TOC_conc * depth, depth, temp, 
                                                 O2_conc, config);
    
    /* CH4 oxidation (aerobic and anaerobic) */
    state->CH4_ox = crive_calc_ch4_oxidation(state->CH4, O2_conc, temp, config);
    
    /* Anaerobic oxidation (assume SO4 from salinity for estuaries) */
    double SO4_est = salinity * 28.0;  /* Rough SO4 estimate from salinity [µmol/L] */
    state->CH4_ox_anaer = crive_calc_ch4_oxidation_anaerobic(state->CH4, SO4_est, 
                                                              temp, config);
    
    /* CH4 ebullition */
    state->CH4_ebul = crive_calc_ch4_ebullition(state->CH4, state->CH4_sat, depth, config);
    
    /* CH4 air-water flux */
    state->CH4_flux = crive_calc_ch4_flux(state->CH4, state->CH4_sat, depth, 
                                           k600, state->Sc_CH4);
    
    return 0;
}
