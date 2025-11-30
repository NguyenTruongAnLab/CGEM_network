/**
 * @file rive_common.h
 * @brief Common definitions and utilities for RIVE biogeochemistry
 * 
 * This header provides shared constants, temperature functions, and
 * utility macros used across all RIVE submodules.
 * 
 * Part of C-RIVE (Unified RIVE v1.0)
 * 
 * References:
 * - Wang et al. (2018) Water Research 144, 341-355
 * - Hasanyar et al. (2022) Biogeosciences
 */

#ifndef RIVE_COMMON_H
#define RIVE_COMMON_H

#include <math.h>

/* ===========================================================================
 * Physical Constants
 * ===========================================================================*/

#define RIVE_T_KELVIN       273.15      /**< Kelvin offset [K] */
#define RIVE_GAS_CONSTANT   8.314       /**< Universal gas constant [J/mol/K] */
#define RIVE_GRAVITY        9.81        /**< Gravitational acceleration [m/s²] */
#define RIVE_RHO_WATER      1000.0      /**< Water density [kg/m³] */

/* ===========================================================================
 * Time Conversion Constants
 * ===========================================================================*/

#define RIVE_SECONDS_PER_DAY    86400.0
#define RIVE_SECONDS_PER_HOUR   3600.0
#define RIVE_HOURS_PER_DAY      24.0

/* ===========================================================================
 * Numerical Constants
 * ===========================================================================*/

#define RIVE_EPS            1e-10       /**< Small value to prevent division by zero */
#define RIVE_MIN_DEPTH      0.1         /**< Minimum water depth [m] */
#define RIVE_MIN_CONC       0.0         /**< Minimum concentration (non-negative) */

/* ===========================================================================
 * Stoichiometric Ratios (Redfield-based)
 * ===========================================================================*/

#define RIVE_C_N_REDFIELD   6.625       /**< C:N ratio (106/16) */
#define RIVE_C_P_REDFIELD   106.0       /**< C:P ratio */
#define RIVE_N_P_REDFIELD   16.0        /**< N:P ratio */
#define RIVE_O2_C_RATIO     1.0         /**< O2:C for aerobic respiration */
#define RIVE_O2_NH4_NIT     1.5         /**< O2:NH4 for step 1 nitrification */
#define RIVE_O2_NO2_NIT     0.5         /**< O2:NO2 for step 2 nitrification */

/* ===========================================================================
 * Molecular Weights [g/mol]
 * ===========================================================================*/

#define RIVE_MW_C           12.011      /**< Carbon */
#define RIVE_MW_N           14.007      /**< Nitrogen */
#define RIVE_MW_P           30.974      /**< Phosphorus */
#define RIVE_MW_O2          31.998      /**< Oxygen */
#define RIVE_MW_CO2         44.009      /**< Carbon dioxide */
#define RIVE_MW_CH4         16.043      /**< Methane */
#define RIVE_MW_N2O         44.013      /**< Nitrous oxide */

/* ===========================================================================
 * Utility Macros
 * ===========================================================================*/

#define RIVE_MAX(a, b)      (((a) > (b)) ? (a) : (b))
#define RIVE_MIN(a, b)      (((a) < (b)) ? (a) : (b))
#define RIVE_CLAMP(x, lo, hi) (RIVE_MIN(RIVE_MAX((x), (lo)), (hi)))
#define RIVE_SQUARE(x)      ((x) * (x))

/** Polynomial evaluation helpers */
#define RIVE_POLY2(a, b, c, x) \
    ((a) + (b)*(x) + (c)*(x)*(x))
#define RIVE_POLY3(a, b, c, d, x) \
    ((a) + (b)*(x) + (c)*(x)*(x) + (d)*(x)*(x)*(x))
#define RIVE_POLY4(a, b, c, d, e, x) \
    ((a) + (b)*(x) + (c)*(x)*(x) + (d)*(x)*(x)*(x) + (e)*(x)*(x)*(x)*(x))
#define RIVE_POLY5(a, b, c, d, e, f, x) \
    ((a) + (b)*(x) + (c)*(x)*(x) + (d)*(x)*(x)*(x) + (e)*(x)*(x)*(x)*(x) + (f)*(x)*(x)*(x)*(x)*(x))

/* ===========================================================================
 * Temperature Functions
 * ===========================================================================*/

/**
 * Arrhenius temperature correction factor
 * 
 * @param temp Temperature [°C]
 * @param T_ref Reference temperature [°C]
 * @param theta Temperature coefficient (typically 1.04-1.08)
 * @return Temperature correction factor [-]
 */
static inline double rive_temp_arrhenius(double temp, double T_ref, double theta) {
    return pow(theta, temp - T_ref);
}

/**
 * Q10 temperature correction factor
 * 
 * @param temp Temperature [°C]
 * @param T_ref Reference temperature [°C]
 * @param Q10 Q10 coefficient (typically 2.0)
 * @return Temperature correction factor [-]
 */
static inline double rive_temp_Q10(double temp, double T_ref, double Q10) {
    return pow(Q10, (temp - T_ref) / 10.0);
}

/**
 * Gaussian temperature response (optimum temperature model)
 * From C-RIVE param_RIVE.h
 * 
 * @param temp Temperature [°C]
 * @param T_opt Optimum temperature [°C]
 * @param dT_spread Temperature spread [°C]
 * @return Temperature factor [0-1]
 */
static inline double rive_temp_gaussian(double temp, double T_opt, double dT_spread) {
    double diff = temp - T_opt;
    return exp(-RIVE_SQUARE(diff) / RIVE_SQUARE(dT_spread));
}

/* ===========================================================================
 * Kinetic Functions
 * ===========================================================================*/

/**
 * Monod (Michaelis-Menten) kinetics
 * 
 * @param S Substrate concentration
 * @param Ks Half-saturation constant (same units as S)
 * @return Limitation factor [0-1]
 */
static inline double rive_monod(double S, double Ks) {
    if (S <= 0.0) return 0.0;
    return S / (S + Ks);
}

/**
 * Inhibition function (inverse Monod)
 * 
 * @param I Inhibitor concentration
 * @param Ki Half-inhibition constant (same units as I)
 * @return Inhibition factor [0-1]
 */
static inline double rive_inhibition(double I, double Ki) {
    return Ki / (I + Ki + RIVE_EPS);
}

/**
 * Liebig's law of the minimum (multi-limitation)
 * 
 * @param f1 First limitation factor
 * @param f2 Second limitation factor
 * @return Minimum of the two factors
 */
static inline double rive_liebig(double f1, double f2) {
    return RIVE_MIN(f1, f2);
}

/* ===========================================================================
 * Unit Conversions
 * ===========================================================================*/

/**
 * Convert rate from /day to /second
 */
static inline double rive_day_to_sec(double rate_per_day) {
    return rate_per_day / RIVE_SECONDS_PER_DAY;
}

/**
 * Convert rate from /second to /day
 */
static inline double rive_sec_to_day(double rate_per_sec) {
    return rate_per_sec * RIVE_SECONDS_PER_DAY;
}

/**
 * Convert rate from /hour to /second
 */
static inline double rive_hour_to_sec(double rate_per_hour) {
    return rate_per_hour / RIVE_SECONDS_PER_HOUR;
}

#endif /* RIVE_COMMON_H */
