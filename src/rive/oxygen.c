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
 * 
 * Garcia & Gordon (1992) combined fit, more accurate than Weiss (1970)
 * Returns µmol/L (assuming seawater density ~1025 kg/m³)
 * 
 * Reference: Garcia, H.E. & Gordon, L.I. (1992) Limnol. Oceanogr. 37(6):1307-1312
 */
double rive_oxygen_saturation(double temp, double salinity) {
    /* Garcia & Gordon (1992) coefficients for µmol/kg */
    double A0 = 5.80871, A1 = 3.20291, A2 = 4.17887, A3 = 5.10006;
    double A4 = -0.0986643, A5 = 3.80369;
    double B0 = -0.00701577, B1 = -0.00770028, B2 = -0.0113864, B3 = -0.00951519;
    double C0 = -2.75915e-7;
    
    /* Scaled temperature */
    double Ts = log((298.15 - temp) / (273.15 + temp));
    
    /* O2 saturation in µmol/kg */
    double O2_sat = exp(A0 + A1*Ts + A2*Ts*Ts + A3*Ts*Ts*Ts + A4*Ts*Ts*Ts*Ts + A5*Ts*Ts*Ts*Ts*Ts
                       + salinity*(B0 + B1*Ts + B2*Ts*Ts + B3*Ts*Ts*Ts)
                       + C0*salinity*salinity);
    
    /* Convert to µmol/L (multiply by density ~1.025 for seawater at 20 PSU) */
    /* For simplicity, use approximate density factor */
    double density = 1.0 + 0.0008 * salinity;  /* Approximate seawater density correction */
    
    return O2_sat * density;  /* [µmol/L] */
}

/**
 * Calculate piston velocity for gas exchange
 * 
 * Uses Wanninkhof (1992) wind parameterization plus O'Connor-Dobbins current term.
 * 
 * CRITICAL FIX (December 2025): The current-driven term was incorrectly scaled.
 * The O'Connor-Dobbins formula gives k in m/day:
 *   k = 3.93 * sqrt(v/H)  [m/day, with v in m/s and H in m]
 * 
 * The current_k_factor parameter scales this (default ~0.3 for moderate turbulence).
 * We must NOT multiply by 100*3600 as that would give cm/hr from an already-converted value.
 */
double rive_piston_velocity(double velocity, double depth, double temp) {
    BiogeoParams *p = rive_get_params();
    
    if (depth < CGEM_MIN_DEPTH) depth = CGEM_MIN_DEPTH;
    
    /* Schmidt number for O2 at given temperature */
    /* Reference: Wanninkhof (1992) for O2 */
    double Sc = 1800.6 - 120.10*temp + 3.7818*temp*temp - 0.047608*temp*temp*temp;
    if (Sc < 100.0) Sc = 100.0;
    
    /* Wind-driven component: k600 = a * U10² [cm/hr] 
     * Reference: Wanninkhof (1992): k600 = 0.31 * U10² for steady winds
     * Convert to actual k using Schmidt number scaling */
    double U10 = p->wind_speed;
    double a = p->wind_coeff;  /* Default ~0.31 */
    double n = p->schmidt_exp; /* Default -0.5 */
    
    double k_wind = a * U10 * U10 * pow(Sc / 660.0, n);  /* [cm/hr] */
    
    /* Current-driven component: O'Connor-Dobbins formula
     * k = 3.93 * sqrt(v/H) [m/day] where v is velocity [m/s], H is depth [m]
     * 
     * Convert to cm/hr: multiply by 100/24 = 4.167
     * current_k_factor is a calibration parameter (default ~0.3-1.0)
     */
    double abs_vel = fabs(velocity);
    double k_current_mday = 3.93 * sqrt(abs_vel / depth);  /* [m/day] */
    double k_current = p->current_k_factor * k_current_mday * (100.0 / 24.0);  /* [cm/hr] */
    
    /* Combined using quadratic addition (independent processes) */
    double k_total = sqrt(k_wind * k_wind + k_current * k_current);
    
    /* Convert cm/hr to m/s for internal calculations */
    return k_total / (100.0 * 3600.0);  /* [m/s] */
}

/**
 * Calculate O2 air-water exchange rate
 * 
 * Positive rate = O2 entering water (undersaturated)
 * Negative rate = O2 leaving water (supersaturated)
 */
double rive_calc_o2_exchange(double o2, double o2_sat, double piston_vel, double depth) {
    if (depth < CGEM_MIN_DEPTH) depth = CGEM_MIN_DEPTH;
    return (piston_vel / depth) * (o2_sat - o2);
}
