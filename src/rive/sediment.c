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
 * SEDIMENT-BIOGEOCHEMISTRY COUPLING (Based on C-RIVE):
 * - POC/SPM ratio: Particulate organic carbon settles/erodes with mud
 * - P-Adsorption: PO4 adsorbs onto SPM (Langmuir isotherm)
 * - These explain Oxygen Sag in ETM and P buffering
 * 
 * Reference: 
 * - Winterwerp (2002), Manning et al. (2010) - Flocculation
 * - Wang et al. (2018), Hasanyar et al. (2022) - C-RIVE sediment exchanges
 * - House & Warwick (1999) - P adsorption on river sediments
 */

#include "../network.h"
#include "../define.h"
#include "rive_params.h"
#include <math.h>
#include <stdio.h>

/* ============================================================================
 * SEDIMENT-BIOGEOCHEMISTRY COUPLING PARAMETERS
 * ============================================================================*/

/* POC/SPM ratio (fraction of POC bound to SPM)
 * Typical values:
 * - Mekong: 1-3% (Liu et al., 2007)
 * - Amazon: 1-2% (Moreira-Turcq et al., 2003)
 * - Temperate rivers: 2-5%
 * This ratio determines how much TOC settles/erodes with mud
 */
#define DEFAULT_POC_SPM_RATIO       0.02    /* 2% = 20 mg C / g SPM */

/* P-Adsorption Langmuir parameters (House & Warwick, 1999)
 * PIP_eq = Kp * PAC * SPM * PO4 / (1 + Kp * PO4)
 * where:
 *   Kp = equilibrium constant [L/µmol]
 *   PAC = P adsorption capacity [µmol P / mg SPM]
 */
#define DEFAULT_KP_ADS              0.05    /* Langmuir Kp [L/µmol] */
#define DEFAULT_PAC_SPM             0.005   /* P capacity [µmol P / mg SPM] */

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
    
    /* =======================================================================
     * SEDIMENT-TOC COUPLING (From C-RIVE: exchanges.c, sedimentation_erosion.c)
     * 
     * A fraction of TOC is Particulate Organic Carbon (POC) bound to SPM.
     * When SPM settles, POC settles. When SPM erodes, POC is resuspended.
     * 
     * This is CRITICAL for the "Oxygen Sag" at the Estuarine Turbidity
     * Maximum (ETM). Without this coupling:
     * - Model predicts high O2 in turbid zones (wrong)
     * - Reality: Resuspended POC consumes O2 immediately
     * 
     * Implementation:
     *   dTOC = f_POC * (Erosion_SPM - Deposition_SPM) / rho_SPM
     * where f_POC = POC/SPM mass ratio (typically 1-3% for Mekong)
     * 
     * Reference: Wang et al. (2018) Water Research 144, 341-355
     * =======================================================================*/
    double *toc = NULL;
    if (branch->conc && branch->num_species > CGEM_SPECIES_TOC) {
        toc = branch->conc[CGEM_SPECIES_TOC];
    }
    
    /* Get POC/SPM ratio from parameters or use default */
    double f_poc = (p && p->poc_spm_ratio > 0.0) ? p->poc_spm_ratio : DEFAULT_POC_SPM_RATIO;
    
    if (toc && spm) {
        for (int i = 1; i <= M; ++i) {
            /* Net SPM flux [kg/m³/s] */
            double net_spm_flux = branch->reaction_rates[CGEM_REACTION_EROSION_V][i] - 
                                  branch->reaction_rates[CGEM_REACTION_DEPOSITION_V][i];
            
            /* POC flux = f_POC × SPM_flux
             * Units: [kg SPM/m³/s] × [kg C/kg SPM] × [1000 g/kg] × [1000 mg/g] × [µmol/12 mg]
             *      = [kg/m³/s] × 0.02 × 1e6 / 12 = [µmol C/m³/s] × 1e-3 = [µmol/L/s]
             * For mg/L output: [kg/m³/s] × f_poc × 1e6 mg/kg / 1000 L/m³ = f_poc × 1e3 mg/L/s
             */
            double toc_flux = net_spm_flux * f_poc * 1e3;  /* Convert to mg C/L/s */
            
            /* Convert mg C/L to µmol C/L (×83.3 = ×1000/12) for internal units */
            toc_flux *= 83.3;  /* µmol C/L/s */
            
            /* Apply: Erosion adds TOC, Deposition removes TOC */
            toc[i] += toc_flux * dt;
            
            /* Ensure non-negative */
            if (toc[i] < 0.0) toc[i] = 0.0;
        }
    }
    
    /* =======================================================================
     * PHOSPHATE ADSORPTION / DESORPTION (From C-RIVE: ads_desorption.c)
     * 
     * PO4 "sticks" to SPM particles. The suspended sediment acts as a
     * BUFFER for dissolved phosphate:
     * - High PO4 (city discharge): SPM absorbs excess P
     * - Low PO4 (ocean): SPM releases adsorbed P
     * 
     * Without this:
     * - Model over-predicts eutrophication downstream of pollution sources
     * - Reality: Much P is "hidden" on particles and unavailable for algae
     * 
     * Implementation: Langmuir Isotherm (equilibrium approximation)
     *   PIP_eq = Kp × PAC × SPM × PO4 / (1 + Kp × PO4)
     * where:
     *   PIP = Particulate Inorganic Phosphorus [µmol P/L]
     *   Kp = Langmuir equilibrium constant [L/µmol]
     *   PAC = P adsorption capacity [µmol P / mg SPM]
     *   SPM = Suspended sediment [mg/L]
     *   PO4 = Dissolved phosphate [µmol P/L]
     * 
     * We use fast equilibrium (instantaneous partitioning) since adsorption
     * kinetics are fast relative to our timestep (300s).
     * 
     * Reference: 
     * - House & Warwick (1999) Sci. Total Environ. 228, 63-74
     * - Wang et al. (2018) - C-RIVE phosphorus module
     * =======================================================================*/
    double *po4 = NULL;
    double *pip = NULL;
    if (branch->conc && branch->num_species > CGEM_SPECIES_PO4) {
        po4 = branch->conc[CGEM_SPECIES_PO4];
    }
    if (branch->conc && branch->num_species > CGEM_SPECIES_PIP) {
        pip = branch->conc[CGEM_SPECIES_PIP];
    }
    
    /* Get P-adsorption parameters from config or use defaults */
    double kp_ads = (p && p->kpads > 0.0) ? p->kpads : DEFAULT_KP_ADS;
    double pac_spm = (p && p->pac > 0.0) ? p->pac : DEFAULT_PAC_SPM;
    int enable_p_ads = (p) ? !p->skip_p_adsorption : 1;  /* Default ON if not skipped */
    
    if (enable_p_ads && po4 && pip && spm) {
        for (int i = 1; i <= M; ++i) {
            if (spm[i] < 1.0) continue;  /* Skip if very low SPM */
            
            /* Total phosphorus = dissolved + particulate */
            double total_P = po4[i] + pip[i];
            if (total_P < 0.01) continue;  /* Skip if no P */
            
            /* Langmuir equilibrium PIP 
             * PIP_eq = Kp × PAC × SPM × PO4 / (1 + Kp × PO4)
             * 
             * Rearrange for PO4_eq given total P:
             * total_P = PO4 + PIP = PO4 + Kp × PAC × SPM × PO4 / (1 + Kp × PO4)
             * 
             * For simplicity, use iterative approach or quadratic solution.
             * Here we use Newton-Raphson (1 iteration is usually sufficient).
             */
            double spm_mg = spm[i];  /* Already in mg/L */
            double max_capacity = pac_spm * spm_mg;  /* Max adsorbable P [µmol/L] */
            
            /* Initial guess: partition based on current ratio */
            double po4_guess = po4[i];
            if (po4_guess < 0.01) po4_guess = 0.5 * total_P;
            
            /* Newton iteration for equilibrium */
            for (int iter = 0; iter < 3; ++iter) {
                double denom = 1.0 + kp_ads * po4_guess;
                double pip_eq = kp_ads * max_capacity * po4_guess / denom;
                
                /* Clamp to max capacity */
                if (pip_eq > max_capacity) pip_eq = max_capacity;
                
                /* Residual: should sum to total_P */
                double residual = po4_guess + pip_eq - total_P;
                
                /* Derivative: d(PO4 + PIP)/d(PO4) */
                double deriv = 1.0 + kp_ads * max_capacity / (denom * denom);
                
                /* Update */
                if (fabs(deriv) > 1e-10) {
                    po4_guess -= residual / deriv;
                }
                
                /* Clamp */
                if (po4_guess < 0.0) po4_guess = 0.01 * total_P;
                if (po4_guess > total_P) po4_guess = 0.99 * total_P;
            }
            
            /* Final equilibrium */
            double denom = 1.0 + kp_ads * po4_guess;
            double pip_eq = kp_ads * max_capacity * po4_guess / denom;
            if (pip_eq > max_capacity) pip_eq = max_capacity;
            if (pip_eq > total_P) pip_eq = total_P - 0.01;
            
            double po4_eq = total_P - pip_eq;
            if (po4_eq < 0.0) po4_eq = 0.0;
            
            /* Update concentrations */
            po4[i] = po4_eq;
            pip[i] = pip_eq;
        }
    }
    
    return 0;
}