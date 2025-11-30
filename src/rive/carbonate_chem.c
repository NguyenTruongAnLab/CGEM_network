/**
 * @file carbonate_chem.c
 * @brief C-RIVE Carbonate Chemistry Module Implementation
 * 
 * Full implementation of the C-RIVE unified inorganic carbon module.
 * Directly ported from C-RIVE calc_co2.c and calc_param.c (Unified RIVE v1.0)
 * 
 * CITATION:
 * Wang, S., Flipo, N., Romary, T., 2018. Time-dependent global sensitivity analysis
 * of the C-RIVE biogeochemical model. Water Research 144, 341-355.
 * 
 * COPYRIGHT: Original C-RIVE code (c) 2023 Contributors to the librive library.
 * Eclipse Public License v2.0
 */

#include "carbonate_chem.h"
#include <stdio.h>
#include <stdlib.h>

/* ===========================================================================
 * Local Helper Macros (from C-RIVE functions.h)
 * ===========================================================================*/

/* Polynomial evaluation helpers */
#define CRIVE_POLY2(a, b, c, x) ((a) + (b)*(x) + (c)*(x)*(x))
#define CRIVE_POLY3(a, b, c, d, x) ((a) + (b)*(x) + (c)*(x)*(x) + (d)*(x)*(x)*(x))
#define CRIVE_POLY4(a, b, c, d, e, x) ((a) + (b)*(x) + (c)*(x)*(x) + (d)*(x)*(x)*(x) + (e)*(x)*(x)*(x)*(x))
#define CRIVE_POLY5(a, b, c, d, e, f, x) ((a) + (b)*(x) + (c)*(x)*(x) + (d)*(x)*(x)*(x) + (e)*(x)*(x)*(x)*(x) + (f)*(x)*(x)*(x)*(x)*(x))

/* Time conversion constants */
#define CRIVE_SECONDS_PER_HOUR 3600.0

/* ===========================================================================
 * Initialization
 * ===========================================================================*/

void crive_carbonate_init_config(CarbonateConfig *config) {
    if (!config) return;
    
    config->kb = CRIVE_KB_DEFAULT;
    config->bor0 = CRIVE_BOR0_DEFAULT;
    config->pCO2_atm = CGEM_PCO2_ATM_DEFAULT;
    config->k600_method = K600_ESTUARINE;  /* Use Abril et al. (2009) for estuaries */
    config->k600_user = 1.0e-5;  /* 1e-5 m/s default */
    config->strahler_order = 7;   /* Large estuary default */
    config->wind_speed = 3.0;     /* Default wind speed [m/s] */
    config->fetch = 1000.0;       /* Default fetch [m] */
}

/* ===========================================================================
 * Carbonate Dissociation Constants
 * Direct port from C-RIVE calc_co2.c
 * ===========================================================================*/

double crive_calc_K1(double temp_K) {
    /* Harned & Davis (1943) for freshwater
     * pK1 = -126.34048 + 6320.813/T + 19.568224*ln(T)
     */
    double pK1 = -126.34048 + 6320.813 / temp_K + 19.568224 * log(temp_K);
    return pow(10.0, -pK1);
}

double crive_calc_K2(double temp_K) {
    /* Harned & Scholes (1941) for freshwater
     * pK2 = -90.18333 + 5143.692/T + 14.613358*ln(T)
     */
    double pK2 = -90.18333 + 5143.692 / temp_K + 14.613358 * log(temp_K);
    return pow(10.0, -pK2);
}

/* ===========================================================================
 * pH Calculation - Direct port from C-RIVE calc_pH()
 * Uses analytical Cardano-Vieta solution for cubic [H+] equation
 * ===========================================================================*/

double crive_calc_pH(double DIC_mmol, double TA_mmol, double temp,
                     double kb, double bor0) {
    
    double temp_K = temp + CRIVE_T_KELVIN;
    double K1 = crive_calc_K1(temp_K);
    double K2 = crive_calc_K2(temp_K);
    
    /* Convert DIC from mmol/m³ to µmol/kg for calculation
     * rho_water ~ 1e6 g/m³, so mmol/m³ * 1e6 / rho = µmol/kg
     * For freshwater, this simplifies to just the value */
    double rho = crive_calc_water_density(temp, 0.0, 0.0);  /* g/m³ */
    double DIC_umol_kg = DIC_mmol / rho * 1e6;  /* µmol/kg */
    double TA = TA_mmol;  /* Already in µmol/L ~ µmol/kg for freshwater */
    
    /* Prevent division by zero */
    if (TA <= 0.0) TA = 1.0;
    
    /* Cardano coefficients from C-RIVE calc_pH()
     * X = (1 - bor0/TA)*kb + (1 - DIC/TA)*K1
     * Y = (1 - (bor0 + DIC)/TA)*K1*kb + (1 - 2*DIC/TA)*K1*K2
     * Z = (1 - (bor0 + 2*DIC)/TA)*K1*K2*kb
     */
    double Xdiss = (1.0 - bor0 / TA) * kb + (1.0 - DIC_umol_kg / TA) * K1;
    double Ydiss = (1.0 - (bor0 + DIC_umol_kg) / TA) * K1 * kb + 
                   (1.0 - 2.0 * DIC_umol_kg / TA) * K1 * K2;
    double Zdiss = (1.0 - (bor0 + 2.0 * DIC_umol_kg) / TA) * K1 * K2 * kb;
    
    /* Cardano-Vieta solution for [H+] */
    double aCulb = (Xdiss * Xdiss - 3.0 * Ydiss) / 9.0;
    double bCulb = -(2.0 * pow(Xdiss, 3) - 9.0 * Xdiss * Ydiss + 27.0 * Zdiss) / 54.0;
    
    /* Handle numerical issues */
    if (aCulb <= 0.0) aCulb = 1e-20;
    
    double arg = bCulb / pow(aCulb, 1.5);
    /* Clamp arg to [-1, 1] for acos */
    if (arg > 1.0) arg = 1.0;
    if (arg < -1.0) arg = -1.0;
    
    double phyCulb = acos(arg);
    
    /* Hydrogen ion activity [mol/L] */
    double Hdiss = 2.0 * sqrt(aCulb) * cos(phyCulb / 3.0) - Xdiss / 3.0;
    
    /* Ensure positive [H+] */
    if (Hdiss <= 0.0) Hdiss = 1e-9;
    
    return Hdiss;
}

/* ===========================================================================
 * CO2 Solubility - Direct port from C-RIVE calc_co2_solubility()
 * Weiss (1974) Marine Chemistry
 * ===========================================================================*/

double crive_calc_co2_solubility(double temp_K, double salinity) {
    /* Weiss (1974) coefficients for CO2 solubility in seawater
     * ln(k0) = A1 + A2*(100/T) + A3*ln(T/100) + S*(B1 + B2*(T/100) + B3*(T/100)²)
     */
    const double A1 = -58.0931;
    const double A2 = 90.5069;
    const double A3 = 22.2940;
    
    double T100 = temp_K / 100.0;
    
    /* ln(k0) for freshwater (S=0) */
    double ln_k0 = A1 + A2 * (100.0 / temp_K) + A3 * log(T100);
    
    /* Add salinity correction */
    if (salinity > 0.0) {
        ln_k0 += salinity * CRIVE_POLY2(0.023517, -0.023656, 0.0047036, T100);
    }
    
    return exp(ln_k0);  /* [mol/kg/atm] = [mmol/g/atm] */
}

/* ===========================================================================
 * Dry-Wet Air Conversion - Direct port from C-RIVE calc_co2_dry_wet_conversion()
 * Weiss & Price (1980)
 * ===========================================================================*/

double crive_calc_co2_dry_wet(double temp_air, double pCO2_atm, double salinity) {
    double temp_air_K = temp_air + CRIVE_T_KELVIN;
    
    /* Water vapor pressure (Weiss & Price 1980) */
    double vph2o = exp(24.4543 - 67.4509 * (100.0 / temp_air_K) - 
                       4.8489 * log(temp_air_K / 100.0) - 0.000544 * salinity);
    
    /* Dry air xCO2 from wet air pCO2 (Doe 1994 or Pierrot et al. 2006) */
    double xco2 = (1.0 - vph2o) * pCO2_atm;
    
    return xco2;  /* [µatm] */
}

/* ===========================================================================
 * Water Density - Direct port from C-RIVE calc_water_density()
 * International One Atmosphere Equation (Millero & Poisson 1981)
 * ===========================================================================*/

double crive_calc_water_density(double temp_water, double salinity, double pressure) {
    /* Coefficients for pure water density (rho0) */
    const double a0 = 0.999842594e03, a1 = 6.793952e-02, a2 = -9.095290e-03;
    const double a3 = 1.001685e-04, a4 = -1.120083e-06, a5 = 6.536332e-09;
    
    /* Coefficients for A (salinity correction linear term) */
    const double b0 = 8.24493e-01, b1 = -4.0899e-03, b2 = 7.6438e-05;
    const double b3 = -8.2467e-07, b4 = 5.3875e-09;
    
    /* Coefficients for B (salinity correction sqrt term) */
    const double c0 = -5.72466e-03, c1 = 1.0227e-04, c2 = -1.6546e-06;
    
    /* Coefficient C (salinity correction squared term) */
    const double d0 = 4.8314e-04;
    
    /* Bulk modulus coefficients */
    const double e0 = 1.965221e04, e1 = 1.484206e02, e2 = -2.327105e00;
    const double e3 = 1.360477e-02, e4 = -5.155288e-05;
    const double f0 = 5.46746e01, f1 = -0.603459e00, f2 = 1.09987e-02, f3 = -6.1670e-05;
    const double g0 = 7.944e-02, g1 = 1.6483e-02, g2 = -5.3009e-04;
    const double h0 = 3.239908e00, h1 = 1.43713e-03, h2 = 1.16092e-04, h3 = -5.77905e-07;
    const double i0 = 2.2838e-03, i1 = -1.0981e-05, i2 = -1.6078e-06;
    const double j0 = 1.91075e-04;
    const double k0 = 8.50935e-05, k1 = -6.12293e-06, k2 = 5.2787e-08;
    const double m0 = -9.9348e-07, m1 = 2.0816e-08, m2 = 9.1697e-10;
    
    double T = temp_water;
    double S = salinity;
    double P = pressure;
    
    /* Pure water compressibility factors */
    double kw = CRIVE_POLY4(e0, e1, e2, e3, e4, T);
    double aw = CRIVE_POLY3(h0, h1, h2, h3, T);
    double bw = CRIVE_POLY2(k0, k1, k2, T);
    
    /* Salinity corrections for compressibility */
    double b = bw + CRIVE_POLY2(m0, m1, m2, T) * S;
    double a = aw + CRIVE_POLY2(i0, i1, i2, T) * S + j0 * pow(S, 1.5);
    
    /* Secant bulk modulus k(S,T,0) */
    double knul = kw + CRIVE_POLY3(f0, f1, f2, f3, T) * S;
    knul += CRIVE_POLY2(g0, g1, g2, T) * pow(S, 1.5);
    
    /* Secant bulk modulus k(S,T,P) */
    double kp = CRIVE_POLY2(knul, a, b, P);
    
    /* Pure water density */
    double row = CRIVE_POLY5(a0, a1, a2, a3, a4, a5, T);
    
    /* Seawater density at surface (P=0) */
    double ronul = row + CRIVE_POLY4(b0, b1, b2, b3, b4, T) * S;
    ronul += CRIVE_POLY2(c0, c1, c2, T) * pow(S, 1.5) + d0 * S * S;
    
    /* Final density with pressure correction */
    double density = ronul / (1.0 - P / kp) * 1000.0;  /* [g/m³] */
    
    return density;
}

/* ===========================================================================
 * CO2 Saturation - Direct port from C-RIVE calc_co2_sat()
 * ===========================================================================*/

double crive_calc_co2_sat(double temp_water, double temp_air, 
                          double pCO2_atm, double salinity) {
    /* Get CO2 solubility [mol/kg/atm] */
    double k0 = crive_calc_co2_solubility(temp_water + CRIVE_T_KELVIN, salinity);
    
    /* Convert wet-air pCO2 to dry-air xCO2 */
    double xCO2 = crive_calc_co2_dry_wet(temp_air, pCO2_atm, salinity);
    
    /* Get water density [g/m³] */
    double rho = crive_calc_water_density(temp_water, salinity, 0.0);
    
    /* CO2 saturation [mmol/m³]
     * = k0 [mol/kg/atm] * xCO2 [µatm] * 1e-6 [atm/µatm] * rho [g/m³] * 1e-3 [kg/g]
     * = k0 * xCO2 * 1e-6 * rho * 1e-3
     * = k0 * xCO2 * rho * 1e-9 * 1e6 [mmol/mol]
     * Actually: mmol/m³ = mol/m³ * 1e3 = (mol/kg * kg/m³) * 1e3
     *                   = k0 * xCO2 * 1e-6 * (rho/1000) * 1e3
     *                   = k0 * xCO2 * rho * 1e-6
     */
    double Csat = k0 * xCO2 * 1e-6 * rho;  /* [mmol/m³] */
    
    return Csat;
}

/* ===========================================================================
 * Schmidt Number - From Wanninkhof (1992) / C-RIVE
 * ===========================================================================*/

double crive_calc_co2_schmidt(double temp) {
    /* Schmidt number for CO2 in freshwater
     * Sc = 1911 - 118.11*T + 3.4527*T² - 0.04132*T³
     * From Wanninkhof (1992), adjusted for freshwater
     */
    return CRIVE_POLY3(1911.0, -118.11, 3.4527, -0.04132, temp);
}

/* ===========================================================================
 * k600 Calculation - Direct port from C-RIVE calc_k600_co2()
 * Raymond et al. (2012) Limnology and Oceanography
 * ===========================================================================*/

double crive_calc_k600(double velocity, double depth, K600Method method,
                       int strahler, double k600_user) {
    double k600 = 0.0;
    
    /* Ensure minimum depth */
    if (depth < 0.01) depth = 0.01;
    
    switch (method) {
        case K600_RESERVOIRS:
            /* Alin et al. (2011) - reservoirs formula */
            k600 = (13.82 + 0.35 * velocity * 100.0) / 100.0;  /* [m/h] */
            k600 /= CRIVE_SECONDS_PER_HOUR;  /* [m/s] */
            break;
            
        case K600_STRAHLER:
            /* River formula depending on Strahler order */
            switch (strahler) {
                case 1:
                case 2:
                case 3:
                case 4:
                case 5:
                    /* Small streams (Strahler 1-5): Alin et al. (2011) */
                    k600 = (13.82 + 0.35 * velocity * 100.0) / 100.0;  /* [m/h] */
                    k600 /= CRIVE_SECONDS_PER_HOUR;  /* [m/s] */
                    break;
                    
                case 6:
                    /* Medium rivers (Strahler 6): O'Connor & Dobbins (1958) */
                    k600 = 1.5 * sqrt(velocity / depth) / 100.0;  /* [m/h] */
                    k600 /= CRIVE_SECONDS_PER_HOUR;  /* [m/s] */
                    break;
                    
                default:
                    /* Large rivers (Strahler > 6): Ho et al. (2016) */
                    k600 = 0.55 * sqrt(velocity / depth) / 100.0;  /* [m/h] */
                    k600 /= CRIVE_SECONDS_PER_HOUR;  /* [m/s] */
                    break;
            }
            break;
            
        case K600_ESTUARINE:
            /* Abril et al. (2009) for tidal estuaries
             * Based on measurements in the Scheldt, Gironde, and other European estuaries
             * k600 = 1.0 + 1.719 * sqrt(v) [cm/h]
             * where v = current velocity [m/s]
             * 
             * This formulation accounts for:
             * - Tidal currents generating turbulence at the water surface
             * - Typical estuarine mixing conditions
             * - Valid for large tidal estuaries (Mekong, Saigon, Scheldt)
             * 
             * Reference: Abril, G., Commarieu, M.-V., Sottolichio, A., Bretel, P., Guérin, F.
             * (2009) Turbidity limits gas exchange in a large macrotidal estuary.
             * Estuarine, Coastal and Shelf Science, 83, 342-348.
             */
            {
                double v_abs = fabs(velocity);
                k600 = (1.0 + 1.719 * sqrt(v_abs)) / 100.0;  /* [m/h] */
                k600 /= CRIVE_SECONDS_PER_HOUR;  /* [m/s] */
            }
            break;
            
        case K600_BORGES:
            /* Borges et al. (2004) wind-driven estuarine formula
             * k600 = 1.0 + 1.0 * U10 + 0.25 * U10² [cm/h]
             * where U10 = wind speed at 10m [m/s]
             * 
             * Accounts for wind-driven gas exchange in open estuarine waters
             * Better for wider sections with significant wind fetch
             * 
             * Reference: Borges, A.V., Vanderborght, J.-P., Schiettecatte, L.-S.,
             * Gazeau, F., Ferrón-Smith, S., Delille, B., Frankignoulle, M. (2004)
             * Variability of the gas transfer velocity of CO2 in a macrotidal estuary
             * Estuaries, 27, 593-603.
             */
            {
                /* Use wind speed if provided, otherwise estimate from velocity */
                double U10 = k600_user > 0 ? k600_user : 3.0;  /* Default 3 m/s */
                k600 = (1.0 + 1.0 * U10 + 0.25 * U10 * U10) / 100.0;  /* [m/h] */
                k600 /= CRIVE_SECONDS_PER_HOUR;  /* [m/s] */
            }
            break;
            
        case K600_USER:
        default:
            /* User-defined constant */
            k600 = k600_user;
            break;
    }
    
    return k600;
}

/* ===========================================================================
 * CO2 Air-Water Flux - Direct port from C-RIVE rea_degassing_CO2()
 * ===========================================================================*/

double crive_calc_co2_flux(double CO2_conc, double CO2_sat, double depth,
                           double k600, double Sc) {
    /* Ensure minimum depth */
    if (depth < 0.01) depth = 0.01;
    
    /* Gas transfer velocity corrected for Schmidt number */
    double kg = k600 * sqrt(600.0 / Sc);  /* [m/s] */
    
    /* CO2 flux [mmol/m³/s] (positive = into water, i.e., undersaturation) */
    double Fwa = kg * (CO2_sat - CO2_conc) / depth;
    
    return Fwa;
}

/* ===========================================================================
 * Full Carbonate System Calculation
 * Combines all components from C-RIVE calc_CO2()
 * ===========================================================================*/

int crive_calc_carbonate_system(CarbonateState *state, double temp, 
                                 double salinity, const CarbonateConfig *config) {
    if (!state || !config) return -1;
    
    double temp_K = temp + CRIVE_T_KELVIN;
    
    /* Get dissociation constants */
    double K1 = crive_calc_K1(temp_K);
    double K2 = crive_calc_K2(temp_K);
    
    /* Calculate [H+] from DIC and TA */
    double Hdiss = crive_calc_pH(state->DIC, state->TA, temp, 
                                  config->kb, config->bor0);
    
    /* Store pH */
    state->pH = -log10(Hdiss);
    
    /* Calculate CO2 speciation from DIC (Marescaux et al. 2019)
     * CO2 = (H/K1) * DIC / (1 + K2/H + H/K1)
     * HCO3 = DIC / (1 + H/K1 + K2/H)
     * CO3 = (K2/H) * DIC / (H/K1 + K2/H + 1)
     */
    double denom = 1.0 + K2 / Hdiss + Hdiss / K1;
    if (denom <= 0.0) denom = 1e-10;
    
    state->CO2 = (Hdiss / K1) * state->DIC / denom;
    state->HCO3 = state->DIC / (1.0 + Hdiss / K1 + K2 / Hdiss);
    state->CO3 = (K2 / Hdiss) * state->DIC / (Hdiss / K1 + K2 / Hdiss + 1.0);
    
    /* Calculate CO2 solubility and saturation */
    state->k0 = crive_calc_co2_solubility(temp_K, salinity);
    state->Csat = crive_calc_co2_sat(temp, temp, config->pCO2_atm, salinity);
    
    /* Calculate pCO2 [µatm]
     * pCO2 = CO2 / k0 (Henry's law)
     * Need to convert units properly */
    double rho = crive_calc_water_density(temp, salinity, 0.0);
    state->pCO2 = state->CO2 / (state->k0 * rho * 1e-6);
    if (state->pCO2 < 0.0) state->pCO2 = 0.0;
    
    /* Calculate Schmidt number */
    state->Sc = crive_calc_co2_schmidt(temp);
    
    return 0;
}
