/**
 * @file phytoplankton.h
 * @brief RIVE Phytoplankton Module
 * 
 * Phytoplankton growth, light limitation, nutrient uptake, mortality.
 * Implements PHY1 (diatoms) and PHY2 (green algae).
 * 
 * Reference: Billen et al. (1994), Garnier et al. (2002)
 */

#ifndef RIVE_PHYTOPLANKTON_H
#define RIVE_PHYTOPLANKTON_H

#include "../network.h"

/**
 * Temperature-dependent rate kinetics ID
 */
typedef enum {
    TEMP_KIN_PBMAX = 0,   /**< Max photosynthesis rate */
    TEMP_KIN_KMAINT,      /**< Maintenance respiration */
    TEMP_KIN_KMORT,       /**< Mortality rate */
    TEMP_KIN_KHET,        /**< Heterotrophic decomposition */
    TEMP_KIN_KDENIT,      /**< Denitrification */
    TEMP_KIN_KNIT         /**< Nitrification */
} TempKinID;

/**
 * Temperature-dependent rate function
 * @param id Kinetic type
 * @param base_rate Base rate at 20°C
 * @param temp Temperature [°C]
 * @return Temperature-corrected rate
 */
double rive_temp_kin(TempKinID id, double base_rate, double temp);

/**
 * Light limitation function (depth-integrated)
 * Based on Eilers-Peeters P-I curve
 * 
 * @param alpha Initial slope of P-I curve
 * @param pbmax Maximum photosynthesis rate
 * @param I0 Surface irradiance [W/m²]
 * @param kd Light attenuation coefficient [1/m]
 * @param depth Water depth [m]
 * @return Depth-integrated light limitation factor
 */
double rive_light_limitation(double alpha, double pbmax, double I0, double kd, double depth);

/**
 * Calculate phytoplankton growth rates
 * 
 * @param branch Branch structure with species concentrations
 * @param i Grid index
 * @param temp Temperature [°C]
 * @param depth Water depth [m]
 * @param gpp1 Output: GPP for diatoms
 * @param gpp2 Output: GPP for green algae
 * @param npp_no3_1 Output: NPP from NO3 for diatoms
 * @param npp_no3_2 Output: NPP from NO3 for greens
 * @param npp_nh4_1 Output: NPP from NH4 for diatoms
 * @param npp_nh4_2 Output: NPP from NH4 for greens
 * @param mort1 Output: Mortality for diatoms
 * @param mort2 Output: Mortality for greens
 */
void rive_calc_phytoplankton(Branch *branch, int i, double temp, double depth,
                             double *gpp1, double *gpp2,
                             double *npp_no3_1, double *npp_no3_2,
                             double *npp_nh4_1, double *npp_nh4_2,
                             double *mort1, double *mort2);

#endif /* RIVE_PHYTOPLANKTON_H */
