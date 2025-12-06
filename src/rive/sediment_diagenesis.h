/**
 * @file sediment_diagenesis.h
 * @brief Active Sediment Layer (Simplified Diagenesis) Module
 * 
 * CRITICAL MISSING FEATURE IMPLEMENTATION (December 2025 Audit)
 * 
 * PURPOSE:
 * ========
 * Without an active sediment layer, the model has fixed benthic fluxes that
 * don't respond to changes in pollution. If you reduce urban pollution by 50%,
 * the sediment still releases the same fluxes forever - missing the "legacy 
 * load" effect where historical pollution accumulates in sediments.
 * 
 * THE DELFT3D WAY (Complex):
 * - Multi-layer S1/S2 soil chemistry with full diagenesis
 * - Requires extensive parameterization (10+ layers, 20+ parameters)
 * 
 * THE CGEM "SMART" WAY (Parsimonious):
 * - Single "Sediment Organic Carbon" (SOC) state variable per grid cell
 * - INPUT: Deposition of POC + Dead Algae from water column
 * - PROCESS: SOC decays with slow first-order rate (k ~ 0.001-0.005 /day)
 * - OUTPUT: Dynamic benthic fluxes (SOD, NH4, PO4, DIC, CH4, N2O)
 * 
 * This physically links upstream loads to downstream benthic fluxes WITHOUT
 * needing vertical sediment profiles.
 * 
 * STOICHIOMETRY:
 * ==============
 * SOC mineralization produces:
 *   O2 demand:  SOC + O2 → CO2 + H2O           (1 mol O2 / mol C)
 *   NH4 release: SOC → NH4                      (redn mol N / mol C)
 *   PO4 release: SOC → PO4                      (redp mol P / mol C)
 *   DIC release: SOC → DIC                      (1 mol C / mol C × RQ)
 *   CH4 production: SOC → CH4 (anaerobic only)  (fraction of C under anoxia)
 * 
 * Reference:
 * - Chapra (2008) Surface Water-Quality Modeling, Chapter 24
 * - DiToro (2001) Sediment Flux Modeling
 * - Soetaert et al. (1996) Estuarine, Coastal and Shelf Science
 * 
 * @author CGEM Development Team
 * @date December 2025
 */

#ifndef SEDIMENT_DIAGENESIS_H
#define SEDIMENT_DIAGENESIS_H

#include "../network.h"

/* ============================================================================
 * SEDIMENT ORGANIC CARBON (SOC) PARAMETERS
 * ============================================================================
 * These parameters control the active sediment layer dynamics.
 * All are configurable via biogeo_params.txt
 * ============================================================================*/

/**
 * SOC Pool Configuration
 * Contains all parameters for the simplified diagenesis model
 */
typedef struct {
    /* SOC decay kinetics */
    double k_soc_20C;           /* SOC mineralization rate at 20°C [1/day] (default 0.003) */
    double soc_Q10;             /* Temperature coefficient for SOC decay (default 2.0) */
    
    /* Initial conditions */
    double soc_init;            /* Initial SOC pool [g C/m²] (default 500) */
    double soc_max;             /* Maximum SOC capacity [g C/m²] (default 5000) */
    
    /* Burial (permanent removal from active layer) */
    double k_burial;            /* Burial rate coefficient [1/day] (default 0.0001) */
    
    /* Stoichiometry for flux calculations */
    double soc_redn;            /* N:C ratio in sediment organic matter (default 0.151) */
    double soc_redp;            /* P:C ratio in sediment organic matter (default 0.0094) */
    double soc_RQ;              /* Respiratory quotient CO2:O2 (default 1.2) */
    
    /* Anaerobic fraction (for CH4 production) */
    double f_anaerobic;         /* Fraction of SOC decay that's anaerobic (default 0.3) */
    
    /* CH4 production parameters */
    double ch4_yield;           /* CH4 yield from anaerobic OC decay [mol CH4/mol C] (default 0.5) */
    double ch4_o2_inhibit;      /* O2 threshold for CH4 production [µmol/L] (default 50) */
    
    /* N2O production parameters */
    double n2o_yield_sediment;  /* N2O yield from sediment N cycling [mol N2O/mol N] (default 0.02) */
    
    /* Control flags */
    int enable_soc;             /* 1 = Enable dynamic SOC, 0 = Use fixed benthic fluxes */
    int loaded;                 /* Flag: parameters loaded */
} SOCConfig;

/**
 * SOC State for a single grid cell
 */
typedef struct {
    double soc_pool;            /* Current SOC pool [g C/m²] */
    double deposition_flux;     /* POC deposition rate [g C/m²/day] */
    double mineralization_rate; /* SOC mineralization rate [g C/m²/day] */
    double burial_rate;         /* SOC burial rate [g C/m²/day] */
    
    /* Output fluxes (positive = release from sediment to water) */
    double sod_flux;            /* Sediment O2 demand [mmol O2/m²/day] */
    double nh4_flux;            /* NH4 release [mmol N/m²/day] */
    double po4_flux;            /* PO4 release [mmol P/m²/day] */
    double dic_flux;            /* DIC release [mmol C/m²/day] */
    double ch4_flux;            /* CH4 release [µmol CH4/m²/day] */
    double n2o_flux;            /* N2O release [nmol N2O/m²/day] */
} SOCState;

/* ============================================================================
 * FUNCTION DECLARATIONS
 * ============================================================================*/

/**
 * Initialize SOC configuration with defaults
 * @param config Pointer to SOCConfig structure to initialize
 */
void soc_init_config(SOCConfig *config);

/**
 * Get global SOC configuration (loaded from biogeo_params.txt)
 * @return Pointer to global SOCConfig
 */
SOCConfig* soc_get_config(void);

/**
 * Initialize SOC pools for a branch
 * Sets initial SOC values based on configuration
 * 
 * @param branch Branch to initialize
 * @param config SOC configuration
 */
void soc_init_branch(Branch *branch, SOCConfig *config);

/**
 * Calculate POC deposition from water column to sediment
 * 
 * Sources of POC deposition:
 * 1. Dead phytoplankton (phy_death * f_poc_phy)
 * 2. Settling POC (w_s * POC_conc)
 * 3. SPM-bound organic carbon (w_s * SPM * poc_spm_ratio)
 * 
 * @param branch Branch structure
 * @param i Grid cell index
 * @param phy_death Total phytoplankton mortality [µmol C/L/day]
 * @param toc TOC concentration [µmol C/L]
 * @param spm SPM concentration [mg/L]
 * @param ws Settling velocity [m/s]
 * @param depth Water depth [m]
 * @param config SOC configuration
 * @return POC deposition flux [g C/m²/day]
 */
double soc_calc_deposition(Branch *branch, int i, double phy_death, 
                           double toc, double spm, double ws, double depth,
                           SOCConfig *config);

/**
 * Calculate benthic fluxes from SOC mineralization
 * 
 * This is the KEY function that makes benthic fluxes dynamic!
 * Instead of fixed fluxes, they now depend on the accumulated SOC pool.
 * 
 * @param state Output: SOC state and fluxes for this grid cell
 * @param soc_current Current SOC pool [g C/m²]
 * @param deposition POC deposition rate [g C/m²/day]
 * @param temp Water temperature [°C]
 * @param o2 Bottom water O2 concentration [µmol/L]
 * @param config SOC configuration
 */
void soc_calc_fluxes(SOCState *state, double soc_current, double deposition,
                     double temp, double o2, SOCConfig *config);

/**
 * Update SOC pool and apply benthic fluxes for entire branch
 * 
 * This function:
 * 1. Calculates POC deposition at each grid cell
 * 2. Updates SOC pool (+ deposition - mineralization - burial)
 * 3. Calculates dynamic benthic fluxes based on SOC pool
 * 4. Applies fluxes to water column concentrations
 * 
 * @param branch Branch to process
 * @param dt Time step [s]
 * @return 0 on success, -1 on error
 */
int soc_update_branch(Branch *branch, double dt);

/**
 * Get total SOC inventory for a branch
 * Useful for mass balance checking and reporting
 * 
 * @param branch Branch structure
 * @return Total SOC inventory [kg C]
 */
double soc_get_inventory(Branch *branch);

/* ============================================================================
 * REACTION INDICES FOR SOC-RELATED FLUXES
 * These should be added to define.h if not present
 * ============================================================================*/
#ifndef CGEM_REACTION_SOC_DEPOSITION
#define CGEM_REACTION_SOC_DEPOSITION  55  /* POC deposition [g C/m²/day] */
#define CGEM_REACTION_SOC_MINERAL     56  /* SOC mineralization [g C/m²/day] */
#define CGEM_REACTION_SOC_BURIAL      57  /* SOC burial [g C/m²/day] */
#define CGEM_REACTION_SOC_SOD         58  /* Dynamic SOD from SOC [mmol O2/m²/day] */
#define CGEM_REACTION_SOC_NH4         59  /* Dynamic NH4 flux [mmol N/m²/day] */
#define CGEM_REACTION_SOC_PO4         60  /* Dynamic PO4 flux [mmol P/m²/day] */
#define CGEM_REACTION_SOC_DIC         61  /* Dynamic DIC flux [mmol C/m²/day] */
#define CGEM_REACTION_SOC_CH4         62  /* Dynamic CH4 flux [µmol/m²/day] */
#define CGEM_REACTION_SOC_N2O         63  /* Dynamic N2O flux [nmol/m²/day] */
#endif

#endif /* SEDIMENT_DIAGENESIS_H */
