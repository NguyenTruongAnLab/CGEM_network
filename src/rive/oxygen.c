/**
 * @file oxygen.c
 * @brief RIVE Oxygen Dynamics Module Implementation
 * 
 * Reference: Weiss (1970), O'Connor & Dobbins (1958)
 */

#include "oxygen.h"
#include "rive_params.h"
#include "../define.h"
#include <math.h>

/**
 * Calculate oxygen saturation concentration
 * Simplified Weiss (1970) equation
 */
double rive_oxygen_saturation(double temp, double salinity) {
    double A1 = -173.4292, A2 = 249.6339, A3 = 143.3483, A4 = -21.8492;
    double B1 = -0.033096, B2 = 0.014259, B3 = -0.0017;
    
    double Ts = log((298.15 - temp) / (273.15 + temp));
    double O2_sat = exp(A1 + A2*Ts + A3*Ts*Ts + A4*Ts*Ts*Ts + 
                       salinity*(B1 + B2*Ts + B3*Ts*Ts));
    return O2_sat; /* [µmol/kg] */
}

/**
 * Calculate piston velocity for gas exchange
 */
double rive_piston_velocity(double velocity, double depth, double temp) {
    BiogeoParams *p = rive_get_params();
    
    if (depth < CGEM_MIN_DEPTH) depth = CGEM_MIN_DEPTH;
    
    /* Schmidt number for CO2 (proxy for O2) */
    double Sc = 2116.8 - 136.25*temp + 4.7353*temp*temp - 0.092307*temp*temp*temp 
                + 0.0007555*temp*temp*temp*temp;
    if (Sc < 100.0) Sc = 100.0;
    
    /* Wind-driven component: k = a * U10² * (Sc/660)^n */
    double U10 = p->wind_speed;
    double a = p->wind_coeff;
    double n = p->schmidt_exp;
    
    double k_wind = a * U10 * U10 * pow(Sc / 660.0, n);  /* [cm/hr] */
    
    /* Current-driven component */
    double abs_vel = fabs(velocity);
    double k_current = p->current_k_factor * sqrt(abs_vel / depth) * 100.0 * 3600.0;  /* [cm/hr] */
    
    /* Combined (quadratic) */
    double k_total = sqrt(k_wind * k_wind + k_current * k_current);
    
    /* Convert cm/hr to m/s */
    return k_total / (100.0 * 3600.0);
}

/**
 * Calculate O2 air-water exchange rate
 */
double rive_calc_o2_exchange(double o2, double o2_sat, double piston_vel, double depth) {
    if (depth < CGEM_MIN_DEPTH) depth = CGEM_MIN_DEPTH;
    return (piston_vel / depth) * (o2_sat - o2);
}
