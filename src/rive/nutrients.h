/**
 * @file nutrients.h
 * @brief RIVE Nutrient Cycling Module
 * 
 * Nitrogen transformations (nitrification, denitrification),
 * phosphorus dynamics, silica cycling.
 * 
 * Reference: Garnier et al. (2002), C-RIVE growth_bactn.c
 */

#ifndef RIVE_NUTRIENTS_H
#define RIVE_NUTRIENTS_H

#include "../network.h"

/**
 * Calculate nitrification rate (NH4 → NO3, or NH4 → NO2 for 2-step)
 * 
 * @param nh4 Ammonium concentration [µmol/L]
 * @param o2 Oxygen concentration [µmol/L]
 * @param temp Temperature [°C]
 * @param knit Base nitrification rate [1/day]
 * @param ko2_nit O2 half-saturation for nitrification [µmol/L]
 * @param knh4 NH4 half-saturation [µmol/L]
 * @return Nitrification rate [µmol/L/s]
 */
double rive_calc_nitrification(double nh4, double o2, double temp,
                               double knit, double ko2_nit, double knh4);

/**
 * Calculate denitrification rate (NO3 → N2)
 * 
 * @param no3 Nitrate concentration [µmol/L]
 * @param o2 Oxygen concentration [µmol/L]
 * @param toc TOC concentration [µmol C/L]
 * @param temp Temperature [°C]
 * @param kdenit Base denitrification rate [1/day]
 * @param kino2 O2 inhibition constant [µmol/L]
 * @param kno3 NO3 half-saturation [µmol/L]
 * @param ktox TOC half-saturation [µmol C/L]
 * @return Denitrification rate [µmol/L/s]
 */
double rive_calc_denitrification(double no3, double o2, double toc, double temp,
                                 double kdenit, double kino2, double kno3, double ktox);

/**
 * Calculate aerobic organic matter degradation
 * 
 * @param toc TOC concentration [µmol C/L]
 * @param o2 Oxygen concentration [µmol/L]
 * @param temp Temperature [°C]
 * @param kox Base oxidation rate [1/day]
 * @param ktox TOC half-saturation [µmol C/L]
 * @param ko2 O2 half-saturation [µmol/L]
 * @return Aerobic degradation rate [µmol C/L/s]
 */
double rive_calc_aerobic_degradation(double toc, double o2, double temp,
                                     double kox, double ktox, double ko2);

/**
 * Calculate silica consumption by diatoms
 * 
 * @param npp_phy1 Net primary production of diatoms [µmol C/L/s]
 * @param redsi Si:C Redfield ratio
 * @return Si consumption rate [µmol Si/L/s]
 */
double rive_calc_silica_consumption(double npp_phy1, double redsi);

#endif /* RIVE_NUTRIENTS_H */
