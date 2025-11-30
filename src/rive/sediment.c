/**
 * @file sediment.c
 * @brief C-GEM Network Sediment Transport Module
 * 
 * Implements SPM erosion and deposition based on shear stress.
 * Matches Fortran CGEM_Sediment.f90
 * 
 * TROPICAL ESTUARY ENHANCEMENT (Mekong Delta):
 * - Salinity-dependent settling velocity (flocculation)
 * - Flocculation drastically increases settling at 0-5 PSU (salt wedge)
 * - Critical for reproducing Estuarine Turbidity Maximum (ETM)
 * 
 * Reference: Winterwerp (2002), Manning et al. (2010)
 */

#include "../network.h"
#include "../define.h"
#include "rive_params.h"
#include <math.h>
#include <stdio.h>

/* -------------------------------------------------------------------------- */

/**
 * Calculate salinity-dependent settling velocity (flocculation model)
 * 
 * In estuaries, clay particles aggregate (flocculate) when freshwater
 * meets saltwater. Maximum flocculation occurs at 2-5 PSU.
 * 
 * ws_floc = ws_primary * (1 + factor_max * sal / (sal + sal_scale))
 * 
 * For sal = 0:   ws_floc = ws (primary settling)
 * For sal = 2:   ws_floc = ws * (1 + 10 * 2/(2+2)) = ws * 6 (6x faster)
 * For sal = 30:  ws_floc = ws * (1 + 10 * 30/(30+2)) = ws * 10.4 (plateau)
 * 
 * Reference: Winterwerp (2002), Manning et al. (2010), Wolanski et al. (2004)
 * 
 * @param ws_primary Primary settling velocity (non-flocculated) [m/s]
 * @param salinity Local salinity [PSU]
 * @param floc_sal_scale Salinity scale parameter (default 2 PSU)
 * @param floc_factor_max Maximum flocculation enhancement (default 10)
 * @return Effective settling velocity with flocculation [m/s]
 */
static double calc_flocculated_settling(double ws_primary, double salinity,
                                         double floc_sal_scale, double floc_factor_max) {
    if (salinity <= 0.0 || ws_primary <= 0.0) {
        return ws_primary;
    }
    
    /* Hyperbolic flocculation model */
    double floc_factor = 1.0 + floc_factor_max * salinity / (salinity + floc_sal_scale);
    
    return ws_primary * floc_factor;
}

/**
 * Calculate bed shear stress from velocity
 * @param velocity Water velocity [m/s]
 * @param depth Water depth [m]
 * @param chezy Chezy coefficient [m^0.5/s]
 * @return Bed shear stress [Pa]
 */
static double bed_shear_stress(double velocity, double depth, double chezy) {
    double rho = CGEM_RHO_WATER;
    double g = CGEM_GRAVITY;
    
    if (depth < CGEM_MIN_DEPTH) depth = CGEM_MIN_DEPTH;
    if (chezy < 10.0) chezy = 50.0;
    
    /* tau_b = rho * g * U^2 / C^2 */
    return rho * g * velocity * velocity / (chezy * chezy);
}

/**
 * Initialize sediment transport parameters for a branch
 * Sets up default values for erosion/deposition parameters
 */
void InitializeSedimentParameters(Branch *branch) {
    if (!branch) return;
    
    /* Initialize sediment parameters with spatial variation */
    for (int i = 0; i <= branch->M + 1; ++i) {
        /* Linear interpolation from downstream to upstream */
        double tau_ero_down = 0.5;  /* Erosion threshold [Pa] */
        double tau_ero_up = 0.2;
        double tau_dep_down = 0.1;  /* Deposition threshold [Pa] */
        double tau_dep_up = 0.05;
        double mero_down = 1e-6;    /* Erosion coefficient [kg/m²/s] */
        double mero_up = 5e-7;
        
        double s = (double)i / (double)branch->M;
        branch->tau_ero[i] = tau_ero_down + (tau_ero_up - tau_ero_down) * s;
        branch->tau_dep[i] = tau_dep_down + (tau_dep_up - tau_dep_down) * s;
        branch->mero[i] = mero_down + (mero_up - mero_down) * s;
    }
}
int Sediment_Branch(Branch *branch, double dt) {
    if (!branch || branch->M <= 0) return -1;
    
    int M = branch->M;
    
    /* Get global parameters for flocculation */
    BiogeoParams *p = rive_get_params();
    double floc_sal_scale = (p && p->floc_sal_scale > 0.0) ? p->floc_sal_scale : 2.0;
    double floc_factor_max = (p && p->floc_factor_max > 0.0) ? p->floc_factor_max : 10.0;
    int enable_flocculation = (p) ? p->enable_flocculation : 1;  /* Default ON */
    
    /* Get SPM and salinity concentration arrays */
    double *spm = NULL;
    double *salinity = NULL;
    if (branch->conc && branch->num_species > CGEM_SPECIES_SPM) {
        spm = branch->conc[CGEM_SPECIES_SPM];
    }
    if (branch->conc && branch->num_species > CGEM_SPECIES_SALINITY) {
        salinity = branch->conc[CGEM_SPECIES_SALINITY];
    }
    if (!spm) return 0; /* No SPM to transport */
    
    /* Calculate erosion and deposition rates */
    for (int i = 1; i <= M; ++i) {
        /* Calculate bed shear stress */
        double tau_b = bed_shear_stress(branch->velocity[i], branch->depth[i], branch->chezyArray[i]);
        
        /* Erosion rate [kg/m²/s] */
        double erosion_rate = 0.0;
        if (branch->tau_ero[i] > 0.0 && tau_b >= branch->tau_ero[i]) {
            erosion_rate = branch->mero[i] * (tau_b / branch->tau_ero[i] - 1.0);
        }
        
        /* =======================================================================
         * SALINITY-DEPENDENT SETTLING VELOCITY (Flocculation Model)
         * 
         * In the Mekong Delta, sediment settling velocity increases dramatically
         * at the "salt wedge" (0-5 PSU) due to flocculation (clay aggregation).
         * Without this, the model cannot reproduce the Estuarine Turbidity 
         * Maximum (ETM), which determines where light limitation is strongest.
         * 
         * Reference: Winterwerp (2002), Wolanski et al. (2004)
         * =======================================================================*/
        double ws_effective = branch->ws;
        
        if (enable_flocculation && salinity) {
            double sal = (salinity[i] > 0.0) ? salinity[i] : 0.0;
            ws_effective = calc_flocculated_settling(branch->ws, sal, 
                                                      floc_sal_scale, floc_factor_max);
        }
        
        /* Deposition rate [kg/m²/s] - using flocculated settling velocity */
        double deposition_rate = 0.0;
        if (branch->tau_dep[i] > 0.0 && tau_b <= branch->tau_dep[i]) {
            deposition_rate = ws_effective * (1.0 - tau_b / branch->tau_dep[i]) * (spm[i] * 1e-3); /* Convert g/L to kg/m³ */
        }
        
        /* Store volumetric rates [kg/m³/s] for consistency with biogeochemistry */
        double depth = branch->depth[i];
        if (depth < CGEM_MIN_DEPTH) depth = CGEM_MIN_DEPTH;
        
        branch->reaction_rates[CGEM_REACTION_EROSION_S][i] = erosion_rate;
        branch->reaction_rates[CGEM_REACTION_EROSION_V][i] = erosion_rate / depth;
        branch->reaction_rates[CGEM_REACTION_DEPOSITION_S][i] = deposition_rate;
        branch->reaction_rates[CGEM_REACTION_DEPOSITION_V][i] = deposition_rate / depth;
        
        /* Update SPM concentration */
        double net_rate = branch->reaction_rates[CGEM_REACTION_EROSION_V][i] - 
                         branch->reaction_rates[CGEM_REACTION_DEPOSITION_V][i];
        spm[i] += net_rate * dt;
        
        /* Ensure non-negative SPM */
        if (spm[i] < 0.0) spm[i] = 0.0;
    }
    
    return 0;
}