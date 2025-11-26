#ifndef CGEM_DEFINE_H
#define CGEM_DEFINE_H

#define CGEM_MAX_SPECIES 20
#define CGEM_MAX_BRANCH_NAME 64
#define CGEM_MAX_PATH 4096

#define CGEM_DEFAULT_BRANCH_CELLS 80
#define CGEM_MAX_BRANCH_CELLS 1000
#define CGEM_MAX_BRANCHES 20
#define CGEM_MAX_TS_STEPS 10000

/* Physical constants (match legacy single-branch solver) */
#define CGEM_GRAVITY 9.81
#define CGEM_PI 3.14159265358979323846
#define CGEM_RHO_WATER 1000.0

/* Simulation defaults (overridden by global/case parameters) */
#define CGEM_DEFAULT_DT_SECONDS 300.0
#define CGEM_DEFAULT_DX_METERS  2000.0

/* Numerical thresholds */
#define CGEM_MIN_DEPTH 0.1
#define CGEM_MIN_WIDTH 1.0
#define CGEM_TOL 1e-5
#define CGEM_MAX_ITER 1000

/* Storage width ratio (RS in Fortran) */
#define CGEM_RS 1.0

/* Tidal parameters */
#define CGEM_M2_FREQ 0.080536912751677847  /* M2 tidal frequency [cycle/hr] */

/* Species indices (matching Fortran CGEMids) */
#define CGEM_SPECIES_SALINITY 0
#define CGEM_SPECIES_PHY1     1
#define CGEM_SPECIES_PHY2     2
#define CGEM_SPECIES_DSI      3
#define CGEM_SPECIES_NO3      4
#define CGEM_SPECIES_NH4      5
#define CGEM_SPECIES_PO4      6
#define CGEM_SPECIES_O2       7
#define CGEM_SPECIES_TOC      8
#define CGEM_SPECIES_SPM      9
#define CGEM_SPECIES_DIC      10
#define CGEM_SPECIES_AT       11
#define CGEM_SPECIES_PCO2     12
#define CGEM_SPECIES_CO2      13
#define CGEM_SPECIES_PH       14
#define CGEM_SPECIES_HS       15
#define CGEM_SPECIES_ALKC     16

/* Reaction indices */
#define CGEM_REACTION_NPP_NO3     0
#define CGEM_REACTION_NPP_NO3_1   1
#define CGEM_REACTION_NPP_NO3_2   2
#define CGEM_REACTION_NPP_NH4     3
#define CGEM_REACTION_NPP_NH4_1   4
#define CGEM_REACTION_NPP_NH4_2   5
#define CGEM_REACTION_GPP_1       6
#define CGEM_REACTION_GPP_2       7
#define CGEM_REACTION_NPP         8
#define CGEM_REACTION_PHY_DEATH   9
#define CGEM_REACTION_PHY_DEATH_1 10
#define CGEM_REACTION_PHY_DEATH_2 11
#define CGEM_REACTION_SI_CONS     12
#define CGEM_REACTION_AER_DEG     13
#define CGEM_REACTION_DENIT       14
#define CGEM_REACTION_NIT         15
#define CGEM_REACTION_O2_EX       16
#define CGEM_REACTION_O2_EX_S     17
#define CGEM_REACTION_CO2_EX      18
#define CGEM_REACTION_CO2_EX_S    19
#define CGEM_REACTION_DIC_REACT   20
#define CGEM_REACTION_TA_REACT    21
#define CGEM_REACTION_HS_REACT    22
#define CGEM_REACTION_EROSION_S   23
#define CGEM_REACTION_EROSION_V   24
#define CGEM_REACTION_DEPOSITION_S 25
#define CGEM_REACTION_DEPOSITION_V 26

#define CGEM_NUM_SPECIES 17
#define CGEM_NUM_REACTIONS 27

/* Utility macros */
#define CGEM_CLAMP(x, lo, hi) (((x) < (lo)) ? (lo) : (((x) > (hi)) ? (hi) : (x)))
#define CGEM_MAX(a, b) (((a) > (b)) ? (a) : (b))
#define CGEM_MIN(a, b) (((a) < (b)) ? (a) : (b))

#endif /* CGEM_DEFINE_H */
