/**
 * @file rive.h
 * @brief Main RIVE Biogeochemistry API
 * 
 * This is the main header for the C-RIVE biogeochemistry module.
 * Include this header to access all RIVE functionality.
 * 
 * C-RIVE is a comprehensive process-based model for carbon, nutrient,
 * and oxygen cycling in aquatic systems, based on Unified RIVE v1.0.
 * 
 * Features:
 * - Complete carbonate chemistry (DIC, TA, pH, pCO2)
 * - Greenhouse gas dynamics (CO2, CH4, N2O)
 * - 2-step nitrification (NH4 → NO2 → NO3)
 * - RK4 adaptive solver
 * 
 * Authors:
 * - Original RIVE: G. Billen, J. Garnier, N. Flipo (Mines Paris)
 * - C-RIVE: S. Wang, L. Vilmin, M. Hasanyar, N. Flipo
 * 
 * Citation:
 * Wang, S., Flipo, N., Romary, T. (2018). Time-dependent global sensitivity
 * analysis of the C-RIVE biogeochemical model. Water Research 144, 341-355.
 */

#ifndef RIVE_H
#define RIVE_H

/* Common definitions and utilities */
#include "rive_common.h"

/* Module headers */
#include "carbonate_chem.h"
#include "ghg_module.h"
#include "rk4_solver.h"

/* ===========================================================================
 * RIVE Module Summary
 * ===========================================================================
 * 
 * biogeo.c          - Main biogeochemistry driver (Biogeo_Branch function)
 *                     Implements: phytoplankton growth, nutrient uptake,
 *                     O2 dynamics, TOC mineralization, carbonate system
 * 
 * carbonate_chem.h/c - Carbonate chemistry
 *                      DIC, TA, pH (Cardano solver), pCO2, CO2 flux
 *                      Henry's law, Schmidt number, k600
 * 
 * ghg_module.h/c    - Greenhouse gas dynamics
 *                     N2O from nitrification/denitrification
 *                     CH4 methanogenesis, oxidation, ebullition
 *                     Air-water exchange for all gases
 * 
 * rk4_solver.h/c    - Numerical solver
 *                     4th-order Runge-Kutta with adaptive time-stepping
 * 
 * ===========================================================================*/

#endif /* RIVE_H */
