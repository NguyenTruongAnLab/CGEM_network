/**
 * @file biogeo.c
 * @brief C-GEM Network Biogeochemical Module
 * 
 * Implements water quality reactions including phytoplankton growth,
 * nutrient cycling, oxygen dynamics, and carbon chemistry.
 * Matches Fortran CGEM_Biogeochem.f90 (simplified version)
 */

#include "network.h"
#include "define.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* -------------------------------------------------------------------------- */

/**
 * Temperature-dependent rate function
 * Matches Fortran: Temp_kin
 */
static double temp_kin(double base_rate, double temp, double ref_temp) {
    double Q10 = 2.0; /* Temperature coefficient */
    return base_rate * pow(Q10, (temp - ref_temp) / 10.0);
}

/**
 * Light limitation function
 * Simplified version of Fortran LightLim
 */
static double light_limitation(double alpha, double pbmax, double I0, double kd, double depth) {
    if (depth < CGEM_MIN_DEPTH) depth = CGEM_MIN_DEPTH;
    /* Simplified light limitation - kd parameter available for future enhancement */
    (void)kd; /* Suppress unused parameter warning */
    double light_term = 1.0 - exp(-alpha * I0 / pbmax);
    return light_term;
}

/**
 * Oxygen saturation concentration [µM]
 */
static double oxygen_saturation(double temp, double salinity) {
    /* Simplified Weiss (1970) equation */
    double A1 = -173.4292, A2 = 249.6339, A3 = 143.3483, A4 = -21.8492;
    double B1 = -0.033096, B2 = 0.014259, B3 = -0.0017;
    
    double Ts = log((298.15 - temp) / (273.15 + temp));
    double O2_sat = exp(A1 + A2*Ts + A3*Ts*Ts + A4*Ts*Ts*Ts + 
                       salinity*(B1 + B2*Ts + B3*Ts*Ts));
    return O2_sat; /* [µmol/kg] */
}

/**
 * Piston velocity for gas exchange [m/s]
 * Simplified version
 */
static double piston_velocity(double velocity, double depth, double salinity, double temp) {
    if (depth < CGEM_MIN_DEPTH) depth = CGEM_MIN_DEPTH;
    
    /* Wanninkhof (1992) parameterization - salinity not used in simplified version */
    (void)salinity; /* Suppress unused parameter warning */
    
    double Sc = 1923.6 - 125.06*temp + 3.997*temp*temp - 0.050*temp*temp*temp; /* Schmidt number */
    double U10 = velocity * 3600.0 / 1000.0; /* Convert to cm/hr */
    double k = 0.31 * U10 * U10 / sqrt(Sc) * 0.00001; /* Convert to m/s */
    
    return k;
}

/**
 * Initialize biogeochemical parameters for a branch
 * Sets up default values for all water quality parameters
 */
void InitializeBiogeoParameters(Branch *branch) {
    if (!branch) return;
    
    /* Water and environmental parameters */
    branch->water_temp = 25.0;    /* Default temperature [°C] */
    branch->ws = 1e-3;            /* Settling velocity [m/s] */
    
    /* Light parameters */
    branch->I0 = 200.0;           /* Surface irradiance [W/m²] */
    branch->kd1 = 0.1;            /* Base attenuation [1/m] */
    branch->kd2_spm = 0.01;       /* SPM attenuation coefficient */
    branch->kd2_phy1 = 0.005;     /* Phytoplankton attenuation */
    branch->kd2_phy2 = 0.005;
    
    /* Phytoplankton parameters */
    branch->alpha1 = 0.02;        /* Initial slope for photosynthesis */
    branch->alpha2 = 0.015;
    branch->pbmax1 = 2.0;         /* Max photosynthetic rate [1/day] */
    branch->pbmax2 = 1.5;
    branch->kexc1 = 0.1;          /* Excretion fraction */
    branch->kexc2 = 0.1;
    branch->kgrowth1 = 0.1;       /* Growth respiration */
    branch->kgrowth2 = 0.1;
    branch->kmaint1 = 0.02;       /* Maintenance respiration [1/day] */
    branch->kmaint2 = 0.015;
    branch->kmort1 = 0.05;        /* Mortality [1/day] */
    branch->kmort2 = 0.04;
    
    /* Nutrient limitation parameters */
    branch->kdsi1 = 5.0;          /* Half-saturation for Si [µM] */
    branch->kn1 = 2.0;            /* Half-saturation for N [µM] */
    branch->kpo41 = 0.5;          /* Half-saturation for P [µM] */
    branch->kn2 = 3.0;
    branch->kpo42 = 0.7;
    
    /* Decomposition parameters */
    branch->kox = 0.1;            /* Aerobic decay [1/day] */
    branch->kdenit = 0.05;        /* Denitrification [1/day] */
    branch->knit = 0.1;           /* Nitrification [1/day] */
    branch->ktox = 50.0;          /* TOC half-saturation [µM] */
    branch->ko2 = 10.0;           /* O2 half-saturation [µM] */
    branch->ko2_nit = 5.0;        /* O2 half-saturation for nitrification */
    branch->kno3 = 5.0;           /* NO3 half-saturation [µM] */
    branch->knh4 = 2.0;           /* NH4 half-saturation [µM] */
    branch->kino2 = 5.0;          /* Inverse O2 half-saturation */
    
    /* Stoichiometric ratios */
    branch->redn = 0.151;         /* N:C ratio */
    branch->redp = 0.01;          /* P:C ratio */
    branch->redsi = 0.1;          /* Si:C ratio */
    
    /* Gas exchange */
    branch->pco2_atm = 400.0;     /* Atmospheric pCO2 [µatm] */
}
int Biogeo_Branch(Branch *branch, double dt) {
    if (!branch || branch->M <= 0 || branch->num_species < CGEM_NUM_SPECIES) {
        return 0;
    }

    int M = branch->M;
    double temp = branch->water_temp;

    /* Get concentration arrays */
    double *salinity = branch->conc[CGEM_SPECIES_SALINITY];
    double *phy1 = branch->conc[CGEM_SPECIES_PHY1];
    double *phy2 = branch->conc[CGEM_SPECIES_PHY2];
    double *dsi = branch->conc[CGEM_SPECIES_DSI];
    double *no3 = branch->conc[CGEM_SPECIES_NO3];
    double *nh4 = branch->conc[CGEM_SPECIES_NH4];
    double *po4 = branch->conc[CGEM_SPECIES_PO4];
    double *o2 = branch->conc[CGEM_SPECIES_O2];
    double *toc = branch->conc[CGEM_SPECIES_TOC];
    double *dic = branch->conc[CGEM_SPECIES_DIC];
    double *at = branch->conc[CGEM_SPECIES_AT];

    if (!salinity || !phy1 || !phy2 || !dsi || !no3 || !nh4 || !po4 || !o2 || !toc || !dic || !at) {
        return 0;
    }

    /* Calculate piston velocity for gas exchange */
    double *pis_vel = (double *)malloc((M + 2) * sizeof(double));
    for (int i = 1; i <= M; ++i) {
        pis_vel[i] = piston_velocity(branch->velocity[i], branch->depth[i],
                                   salinity[i], temp);
    }

    /* Process each grid cell */
    for (int i = 1; i <= M - 1; i += 2) { /* Odd indices (cell centers) */
        double depth = branch->depth[i];
        if (depth < CGEM_MIN_DEPTH) depth = CGEM_MIN_DEPTH;

        /* Light attenuation */
        double kd = branch->kd1 +
                   branch->kd2_spm * (branch->conc[CGEM_SPECIES_SPM] ? branch->conc[CGEM_SPECIES_SPM][i] * 1e-3 : 0.0) +
                   branch->kd2_phy1 * phy1[i] +
                   branch->kd2_phy2 * phy2[i];

        /* Light limitation */
        double lightlim1 = light_limitation(branch->alpha1, branch->pbmax1, branch->I0, kd, depth);
        double lightlim2 = light_limitation(branch->alpha2, branch->pbmax2, branch->I0, kd, depth);

        /* Nutrient limitation */
        double nlim1 = (dsi[i] / (dsi[i] + branch->kdsi1)) *
                      ((no3[i] + nh4[i]) / ((no3[i] + nh4[i]) + branch->kn1)) *
                      (po4[i] / (po4[i] + branch->kpo41));

        double nlim2 = ((no3[i] + nh4[i]) / ((no3[i] + nh4[i]) + branch->kn2)) *
                      (po4[i] / (po4[i] + branch->kpo42));

        /* Ammonium preference */
        double nswitch = nh4[i] / (5.0 + nh4[i]);

        /* Primary production rates [µM/day] */
        double gpp1 = temp_kin(branch->pbmax1, temp, 20.0) * phy1[i] * nlim1 * lightlim1;
        double gpp2 = temp_kin(branch->pbmax2, temp, 20.0) * phy2[i] * nlim2 * lightlim2;

        double npp_no31 = (1.0 - nswitch) * (gpp1 / depth) * (1.0 - branch->kexc1) * (1.0 - branch->kgrowth1) -
                         temp_kin(branch->kmaint1, temp, 20.0) * phy1[i];
        double npp_no32 = (1.0 - nswitch) * (gpp2 / depth) * (1.0 - branch->kexc2) * (1.0 - branch->kgrowth2) -
                         temp_kin(branch->kmaint2, temp, 20.0) * phy2[i];

        double npp_nh41 = nswitch * (gpp1 / depth) * (1.0 - branch->kexc1) * (1.0 - branch->kgrowth1) -
                         temp_kin(branch->kmaint1, temp, 20.0) * phy1[i];
        double npp_nh42 = nswitch * (gpp2 / depth) * (1.0 - branch->kexc2) * (1.0 - branch->kgrowth2) -
                         temp_kin(branch->kmaint2, temp, 20.0) * phy2[i];

        /* Mortality */
        double mort1 = temp_kin(branch->kmort1, temp, 20.0) * phy1[i];
        double mort2 = temp_kin(branch->kmort2, temp, 20.0) * phy2[i];

        /* Silica consumption */
        double si_cons = -branch->redsi * (npp_nh41 + npp_nh42);

        /* Decomposition rates */
        double aer_deg = temp_kin(branch->kox, temp, 20.0) *
                        (toc[i] / (toc[i] + branch->ktox)) *
                        (o2[i] / (o2[i] + branch->ko2));

        double denit = temp_kin(branch->kdenit, temp, 20.0) *
                      (toc[i] / (toc[i] + branch->ktox)) *
                      (branch->kino2 / (o2[i] + branch->kino2)) *
                      (no3[i] / (no3[i] + branch->kno3));

        double nitrif = temp_kin(branch->knit, temp, 20.0) *
                       (o2[i] / (o2[i] + branch->ko2_nit)) *
                       (nh4[i] / (nh4[i] + branch->knh4));

        /* Oxygen exchange */
        double o2_sat = oxygen_saturation(temp, salinity[i]);
        double o2_ex = (pis_vel[i] / depth) * (o2_sat - o2[i]);

        /* Store reaction rates */
        branch->reaction_rates[CGEM_REACTION_GPP_1][i] = gpp1;
        branch->reaction_rates[CGEM_REACTION_GPP_2][i] = gpp2;
        branch->reaction_rates[CGEM_REACTION_NPP_NO3_1][i] = npp_no31;
        branch->reaction_rates[CGEM_REACTION_NPP_NO3_2][i] = npp_no32;
        branch->reaction_rates[CGEM_REACTION_NPP_NH4_1][i] = npp_nh41;
        branch->reaction_rates[CGEM_REACTION_NPP_NH4_2][i] = npp_nh42;
        branch->reaction_rates[CGEM_REACTION_PHY_DEATH_1][i] = mort1;
        branch->reaction_rates[CGEM_REACTION_PHY_DEATH_2][i] = mort2;
        branch->reaction_rates[CGEM_REACTION_SI_CONS][i] = si_cons;
        branch->reaction_rates[CGEM_REACTION_AER_DEG][i] = aer_deg;
        branch->reaction_rates[CGEM_REACTION_DENIT][i] = denit;
        branch->reaction_rates[CGEM_REACTION_NIT][i] = nitrif;
        branch->reaction_rates[CGEM_REACTION_O2_EX][i] = o2_ex;
        branch->reaction_rates[CGEM_REACTION_O2_EX_S][i] = o2_ex * depth;

        /* Update concentrations */
        phy1[i] += (npp_no31 + npp_nh41 - mort1) * dt;
        phy2[i] += (npp_no32 + npp_nh42 - mort2) * dt;
        dsi[i] += si_cons * dt;
        po4[i] += branch->redp * (aer_deg + denit - (npp_no31 + npp_no32 + npp_nh41 + npp_nh42)) * dt;
        o2[i] += (-aer_deg + (npp_nh41 + npp_nh42) * 138.0/106.0 - 2.0 * nitrif + o2_ex) * dt;
        toc[i] += (-aer_deg - denit + (mort1 + mort2)) * dt;
        nh4[i] += branch->redn * (aer_deg - npp_nh41 - npp_nh42) - nitrif * dt;
        no3[i] += -branch->redn * (npp_no31 + npp_no32) + nitrif - denit * 94.4/106.0 * dt;

        /* Ensure non-negative concentrations */
        if (phy1[i] < 0.0) phy1[i] = 0.0;
        if (phy2[i] < 0.0) phy2[i] = 0.0;
        if (dsi[i] < 0.0) dsi[i] = 0.0;
        if (no3[i] < 0.0) no3[i] = 0.0;
        if (nh4[i] < 0.0) nh4[i] = 0.0;
        if (po4[i] < 0.0) po4[i] = 0.0;
        if (o2[i] < 0.0) o2[i] = 0.0;
        if (toc[i] < 0.0) toc[i] = 0.0;
    }

    free(pis_vel);
    return 0;
}
