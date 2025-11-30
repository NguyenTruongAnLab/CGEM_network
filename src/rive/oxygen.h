/**
 * @file oxygen.h
 * @brief RIVE Oxygen Dynamics Module
 * 
 * Dissolved oxygen production, consumption, and air-water exchange.
 * 
 * Reference: Weiss (1970), O'Connor & Dobbins (1958)
 */

#ifndef RIVE_OXYGEN_H
#define RIVE_OXYGEN_H

#include "../network.h"

/**
 * Calculate oxygen saturation concentration
 * 
 * @param temp Temperature [°C]
 * @param salinity Salinity [PSU]
 * @return O2 saturation concentration [µmol/L]
 */
double rive_oxygen_saturation(double temp, double salinity);

/**
 * Calculate piston velocity for gas exchange
 * Combines wind-driven and current-driven components.
 * 
 * @param velocity Water velocity [m/s]
 * @param depth Water depth [m]
 * @param temp Temperature [°C]
 * @return Piston velocity [m/s]
 */
double rive_piston_velocity(double velocity, double depth, double temp);

/**
 * Calculate O2 air-water exchange rate
 * 
 * @param o2 Current O2 concentration [µmol/L]
 * @param o2_sat Saturation concentration [µmol/L]
 * @param piston_vel Piston velocity [m/s]
 * @param depth Water depth [m]
 * @return O2 exchange rate [µmol/L/s] (positive = influx)
 */
double rive_calc_o2_exchange(double o2, double o2_sat, double piston_vel, double depth);

#endif /* RIVE_OXYGEN_H */
