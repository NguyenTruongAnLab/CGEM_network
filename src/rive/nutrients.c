/**
 * @file nutrients.c
 * @brief RIVE Nutrient Cycling Module Implementation
 * 
 * Reference: Garnier et al. (2002), C-RIVE growth_bactn.c
 */

#include "nutrients.h"
#include "phytoplankton.h"
#include "rive_params.h"
#include "../define.h"
#include <math.h>

/**
 * Calculate nitrification rate
 */
double rive_calc_nitrification(double nh4, double o2, double temp,
                               double knit, double ko2_nit, double knh4) {
    if (nh4 <= 0.0 || o2 <= 0.0) return 0.0;
    
    /* Convert rate from /day to /s */
    double knit_s = knit / RIVE_SECONDS_PER_DAY;
    
    /* Temperature correction */
    double rate = rive_temp_kin(TEMP_KIN_KNIT, knit_s, temp);
    
    /* Monod kinetics for O2 and NH4 */
    double o2_lim = o2 / (o2 + ko2_nit);
    double nh4_lim = nh4 / (nh4 + knh4);
    
    return rate * o2_lim * nh4_lim;
}

/**
 * Calculate denitrification rate
 */
double rive_calc_denitrification(double no3, double o2, double toc, double temp,
                                 double kdenit, double kino2, double kno3, double ktox) {
    if (no3 <= 0.0 || toc <= 0.0) return 0.0;
    
    /* Convert rate from /day to /s */
    double kdenit_s = kdenit / RIVE_SECONDS_PER_DAY;
    
    /* Temperature correction */
    double rate = rive_temp_kin(TEMP_KIN_KDENIT, kdenit_s, temp);
    
    /* Substrate limitation */
    double toc_lim = toc / (toc + ktox);
    double no3_lim = no3 / (no3 + kno3);
    
    /* O2 inhibition (only occurs at low O2) */
    double o2_inhib = kino2 / (o2 + kino2);
    
    return rate * toc_lim * o2_inhib * no3_lim;
}

/**
 * Calculate aerobic organic matter degradation
 */
double rive_calc_aerobic_degradation(double toc, double o2, double temp,
                                     double kox, double ktox, double ko2) {
    if (toc <= 0.0 || o2 <= 0.0) return 0.0;
    
    /* Convert rate from /day to /s */
    double kox_s = kox / RIVE_SECONDS_PER_DAY;
    
    /* Temperature correction (Q10 = 2.5) */
    double rate = rive_temp_kin(TEMP_KIN_KHET, kox_s, temp);
    
    /* Monod kinetics */
    double toc_lim = toc / (toc + ktox);
    double o2_lim = o2 / (o2 + ko2);
    
    return rate * toc_lim * o2_lim;
}

/**
 * Calculate silica consumption by diatoms
 */
double rive_calc_silica_consumption(double npp_phy1, double redsi) {
    /* Silica consumed proportional to diatom NPP */
    return fmin(0.0, -redsi * npp_phy1);
}
