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

typedef enum {
    TEMP_KIN_PBMAX = 0,
    TEMP_KIN_KMAINT,
    TEMP_KIN_KMORT,
    TEMP_KIN_KHET,
    TEMP_KIN_KDENIT,
    TEMP_KIN_KNIT
} TempKinID;

/**
 * Temperature-dependent rate function (matches Fortran Temp_kin)
 */
static double temp_kin(TempKinID id, double base_rate, double temp) {
    double delta = temp - 20.0;
    double factor = 1.0;

    switch (id) {
        case TEMP_KIN_KNIT:
            factor = pow(1.08, delta);
            break;
        case TEMP_KIN_KDENIT:
            factor = pow(1.07, delta);
            break;
        case TEMP_KIN_KHET:
            factor = pow(2.5, delta / 10.0);
            break;
        case TEMP_KIN_KMORT:
            factor = exp(0.07 * (delta + 20.0));
            break;
        case TEMP_KIN_KMAINT:
            factor = exp(0.0322 * delta);
            break;
        case TEMP_KIN_PBMAX:
        default:
            factor = pow(1.067, delta);
            break;
    }

    return base_rate * factor;
}

/**
 * Light limitation function
 * Port of the Fortran LightLim implementation
 */
static double light_limitation(double alpha, double pbmax, double I0, double kd, double depth) {
    const double EULER_GAMMA = 0.5772156649015328606;
    if (depth < CGEM_MIN_DEPTH) depth = CGEM_MIN_DEPTH;
    if (alpha <= 0.0 || pbmax <= 0.0 || kd <= 0.0 || I0 <= 0.0) {
        return 0.0;
    }

    double Ebottom = I0 * exp(-kd * depth);
    if (Ebottom < 1e-30) {
        return 0.0;
    }

    double psurf = (I0 * alpha) / pbmax;
    double pbot = (Ebottom * alpha) / pbmax;
    if (psurf <= 0.0 || pbot <= 0.0) {
        return 0.0;
    }

    double ps2 = psurf * psurf;
    double ps3 = ps2 * psurf;
    double ps4 = ps3 * psurf;
    double ps5 = ps4 * psurf;

    double pb2 = pbot * pbot;
    double pb3 = pb2 * pbot;
    double pb4 = pb3 * pbot;
    double pb5 = pb4 * pbot;

    double appGAMMAsurf;
    double appGAMMAbot;

    if (psurf <= 1.0) {
        appGAMMAsurf = -(log(psurf) + EULER_GAMMA - psurf + (ps2 / 4.0) - (ps3 / 18.0) + (ps4 / 96.0) - (ps5 / 600.0));
    } else {
        double frac = psurf + 9.0;
        frac = psurf + 7.0 - (16.0 / frac);
        frac = psurf + 5.0 - (9.0 / frac);
        frac = psurf + 3.0 - (4.0 / frac);
        frac = psurf + 1.0 - (1.0 / frac);
        appGAMMAsurf = exp(-psurf) * (1.0 / frac);
    }

    if (pbot <= 1.0) {
        appGAMMAbot = -(log(pbot) + EULER_GAMMA + pbot - (pb2 / 4.0) + (pb3 / 18.0) - (pb4 / 96.0) + (pb5 / 600.0));
    } else {
        double frac = pbot + 9.0;
        frac = pbot + 7.0 - (16.0 / frac);
        frac = pbot + 5.0 - (9.0 / frac);
        frac = pbot + 3.0 - (4.0 / frac);
        frac = pbot + 1.0 - (1.0 / frac);
        appGAMMAbot = exp(-pbot) * (1.0 / frac);
    }

    double integral = (appGAMMAsurf - appGAMMAbot + log(I0 / Ebottom)) / kd;
    if (!isfinite(integral) || integral < 0.0) {
        return 0.0;
    }
    return integral;
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
 * Henry's constant for CO2 (µmol m⁻³ µatm⁻¹)
 */
static double henry_co2(double temp, double salinity) {
    double Tabs = temp + 273.15;
    if (Tabs <= 0.0) Tabs = 273.15;
    double lnK0 = -574.70126 + (21541.52 / Tabs) - 0.000147759 * Tabs * Tabs + 89.892 * log(Tabs);
    double f = 0.029941 - 0.00027455 * Tabs + 0.00000053407 * Tabs * Tabs;
    return exp(lnK0 + f * salinity) * 1e6;
}

static void carbonate_constants(double temp, double salinity,
                                double *K1, double *K2, double *Kw,
                                double *KH2S, double *KBOH4) {
    double Tabs = temp + 273.15;
    if (Tabs <= 0.0) Tabs = 273.15;
    double sqrtS = sqrt(CGEM_MAX(salinity, 0.0));
    double sal = CGEM_MAX(salinity, 0.0);

    double pK = -14.8425 + (3404.71 / Tabs) + (0.032786 * Tabs);
    double f1 = -0.0230848 - (14.3456 / Tabs);
    double f2 = 0.000691881 + (0.429955 / Tabs);
    *K1 = pow(10.0, -(pK + (f1 * sqrtS) + (f2 * sal)));

    pK = -6.4980 + (2902.39 / Tabs) + (0.02379 * Tabs);
    f1 = -0.458898 + (41.24048 / Tabs);
    f2 = 0.0284743 - (2.55895 / Tabs);
    *K2 = pow(10.0, -(pK + (f1 * sqrtS) + (f2 * sal)));

    pK = (-13847.26 / Tabs) + 148.9652 - (23.6521 * log(Tabs));
    f1 = (118.67 / Tabs) - 5.977 + 1.0495 * log(Tabs);
    f2 = -0.01615;
    *Kw = exp(pK + (f1 * sqrtS) + (f2 * sal));

    pK = (-4276.1 / Tabs) + 141.328 - (23.093 * log(Tabs));
    double ionic = (19.919 * sal) / CGEM_MAX(1000.0 - 1.00198 * sal, 1.0);
    double lnKH2S = (((-13856.0 / Tabs) + 324.57 - 47.986 * log(Tabs)) * sqrt(ionic)) +
                    (((35474.0 / Tabs) - 771.54 + 114.723 * log(Tabs)) * ionic) -
                    ((2698.0 / Tabs) * pow(ionic, 1.5)) +
                    ((1776.0 / Tabs) * ionic * ionic) + pK;
    *KH2S = exp(lnKH2S);

    double term1 = (-8966.90 - 2890.53 * sqrtS - 77.942 * sal + 1.728 * pow(sal, 1.5) - 0.0996 * sal * sal) / Tabs;
    double term2 = (148.0248 + 137.1942 * sqrtS + 1.62142 * sal);
    double term3 = (-24.4344 - 25.085 * sqrtS - 0.2474 * sal) * log(Tabs);
    double term4 = 0.053105 * sqrtS * Tabs;
    *KBOH4 = exp(term1 + term2 + term3 + term4);
}

static double carbonate_residual(double H, double dic, double alk, double sulf, double TotB,
                                 double K1, double K2, double Kw, double KH2S, double KBOH4,
                                 double *xi1, double *xi2, double *gamma1, double *beta1) {
    double Hsafe = CGEM_MAX(H, 1e-12);
    double denom = (Hsafe * Hsafe) + (K1 * Hsafe) + (K1 * K2);
    if (denom <= 0.0) denom = 1e-12;

    *xi1 = (K1 * Hsafe) / denom;
    *xi2 = (*xi1) * K2 / Hsafe;
    *gamma1 = KH2S / (KH2S + Hsafe);
    *beta1 = KBOH4 / (KBOH4 + Hsafe);

    double residual = alk - (((*xi1 + 2.0 * (*xi2)) * dic) + ((*gamma1) * sulf) + ((*beta1) * TotB) +
                              (Kw / Hsafe) - Hsafe);
    return residual;
}

static double calculate_carbonate_system(Branch *branch, int idx, double temp,
                                         double depth, double salinity, double piston_vel) {
    double *dic = branch->conc[CGEM_SPECIES_DIC];
    double *alk = branch->conc[CGEM_SPECIES_AT];
    double *co2 = branch->conc[CGEM_SPECIES_CO2];
    double *pco2 = branch->conc[CGEM_SPECIES_PCO2];
    double *ph = branch->conc[CGEM_SPECIES_PH];
    double *alkc = branch->conc[CGEM_SPECIES_ALKC];
    double *hs = branch->conc[CGEM_SPECIES_HS];

    double henry = henry_co2(temp, salinity);
    if (!dic || !alk || !co2 || !pco2 || !ph || !alkc || !branch->reaction_rates) {
        return henry;
    }

    double sal = CGEM_MAX(salinity, 0.0);
    double Hprot = (ph[idx] > 0.0) ? pow(10.0, -ph[idx]) : 1e-8;
    double sulf = (hs) ? CGEM_MAX(hs[idx], 0.0) : 0.0;
    double TotB = 416.0 * (sal / 35.0);

    double K1, K2, Kw, KH2S, KBOH4;
    carbonate_constants(temp, sal, &K1, &K2, &Kw, &KH2S, &KBOH4);

    double xi1 = 0.0, xi2 = 0.0, gamma1 = 0.0, beta1 = 0.0;
    double guess = carbonate_residual(Hprot, dic[idx], alk[idx], sulf, TotB, K1, K2, Kw, KH2S, KBOH4,
                                      &xi1, &xi2, &gamma1, &beta1);

    double newH = CGEM_MAX(Hprot, 1e-9);
    int iterations = 0;
    while (fabs(guess) > CGEM_TOL && iterations < 50) {
        double carbonAlk = alk[idx] - (gamma1 * sulf) - (beta1 * TotB) - (Kw / CGEM_MAX(newH, 1e-12)) + newH;
        double ratio = dic[idx] / CGEM_MAX(carbonAlk, 1e-12);
        double inner = 1.0 - ratio;
        double discriminant = inner * inner * (K1 * K1) - 4.0 * K1 * K2 * (1.0 - 2.0 * ratio);
        if (discriminant < 0.0) discriminant = 0.0;
        double update = 0.5 * (((ratio) - 1.0) * K1 + sqrt(discriminant));
        if (!isfinite(update) || update <= 0.0) {
            break;
        }
        newH = update;
        guess = carbonate_residual(newH, dic[idx], alk[idx], sulf, TotB, K1, K2, Kw, KH2S, KBOH4,
                                    &xi1, &xi2, &gamma1, &beta1);
        iterations++;
    }

    Hprot = CGEM_MAX(newH, 1e-9);
    ph[idx] = -log10(Hprot);
    alkc[idx] = (double)iterations;

    double denom = 1.0 + (K1 / Hprot) + ((K1 * K2) / (Hprot * Hprot));
    double co2_spec = (denom > 0.0) ? dic[idx] / denom : dic[idx];
    if (!isfinite(co2_spec) || co2_spec < 0.0) {
        co2_spec = CGEM_MAX(dic[idx], 0.0);
    }
    co2[idx] = co2_spec;

    double depth_eff = (depth < CGEM_MIN_DEPTH) ? CGEM_MIN_DEPTH : depth;
    double co2_flux = 0.0;
    if (depth_eff > 0.0) {
        co2_flux = -(piston_vel * 0.913 / depth_eff) * (co2_spec - henry * branch->pco2_atm);
    }
    branch->reaction_rates[CGEM_REACTION_CO2_EX][idx] = co2_flux;
    branch->reaction_rates[CGEM_REACTION_CO2_EX_S][idx] = co2_flux * depth_eff;
    branch->reaction_rates[CGEM_REACTION_HS_REACT][idx] = 0.0;

    double aer_deg = branch->reaction_rates[CGEM_REACTION_AER_DEG][idx];
    double denit = branch->reaction_rates[CGEM_REACTION_DENIT][idx];
    double npp_total = branch->reaction_rates[CGEM_REACTION_NPP][idx];
    double npp_no3 = branch->reaction_rates[CGEM_REACTION_NPP_NO3][idx];
    double npp_nh4 = branch->reaction_rates[CGEM_REACTION_NPP_NH4][idx];
    double nitrif = branch->reaction_rates[CGEM_REACTION_NIT][idx];

    branch->reaction_rates[CGEM_REACTION_DIC_REACT][idx] = co2_flux + aer_deg + denit - npp_total;
    branch->reaction_rates[CGEM_REACTION_TA_REACT][idx] = (15.0 / 106.0) * aer_deg +
                                                          (93.4 / 106.0) * denit -
                                                          2.0 * nitrif +
                                                          (-15.0 / 106.0) * npp_nh4 +
                                                          (17.0 / 106.0) * npp_no3;

    (void)pco2; /* Updated later after flux integration */
    return henry;
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
    if (!branch || branch->M <= 0 || branch->num_species < CGEM_NUM_SPECIES || !branch->reaction_rates) {
        return 0;
    }

    int M = branch->M;
    double temp = branch->water_temp;

    double *salinity = branch->conc[CGEM_SPECIES_SALINITY];
    double *phy1 = branch->conc[CGEM_SPECIES_PHY1];
    double *phy2 = branch->conc[CGEM_SPECIES_PHY2];
    double *dsi = branch->conc[CGEM_SPECIES_DSI];
    double *no3 = branch->conc[CGEM_SPECIES_NO3];
    double *nh4 = branch->conc[CGEM_SPECIES_NH4];
    double *po4 = branch->conc[CGEM_SPECIES_PO4];
    double *o2 = branch->conc[CGEM_SPECIES_O2];
    double *toc = branch->conc[CGEM_SPECIES_TOC];
    double *spm = branch->conc[CGEM_SPECIES_SPM];
    double *dic = branch->conc[CGEM_SPECIES_DIC];
    double *at = branch->conc[CGEM_SPECIES_AT];
    double *co2 = branch->conc[CGEM_SPECIES_CO2];
    double *pco2 = branch->conc[CGEM_SPECIES_PCO2];
    double *ph = branch->conc[CGEM_SPECIES_PH];
    double *alkc = branch->conc[CGEM_SPECIES_ALKC];

    if (!salinity || !phy1 || !phy2 || !dsi || !no3 || !nh4 || !po4 || !o2 || !toc ||
        !dic || !at || !co2 || !pco2 || !ph || !alkc) {
        return 0;
    }

    double *pis_vel = (double *)calloc((size_t)(M + 2), sizeof(double));
    if (!pis_vel) {
        return -1;
    }
    for (int i = 1; i <= M; ++i) {
        pis_vel[i] = piston_velocity(branch->velocity[i], branch->depth[i], salinity[i], temp);
    }

    /* Primary loop: compute reaction rates */
    for (int i = 1; i <= M - 1; i += 2) {
        double depth = CGEM_MAX(branch->depth[i], CGEM_MIN_DEPTH);

        double sal = salinity[i];
        double no3nh4 = no3[i] + nh4[i];

        double kd = branch->kd1 +
                    branch->kd2_spm * (spm ? spm[i] * 1e-3 : 0.0) +
                    branch->kd2_phy1 * phy1[i] +
                    branch->kd2_phy2 * phy2[i];

        double pbmax1 = temp_kin(TEMP_KIN_PBMAX, branch->pbmax1, temp);
        double pbmax2 = temp_kin(TEMP_KIN_PBMAX, branch->pbmax2, temp);
        double lightlim1 = light_limitation(branch->alpha1, pbmax1, branch->I0, kd, depth);
        double lightlim2 = light_limitation(branch->alpha2, pbmax2, branch->I0, kd, depth);

        double dsi_term = dsi[i] + branch->kdsi1;
        double n_term1 = no3nh4 + branch->kn1;
        double n_term2 = no3nh4 + branch->kn2;
        double po4_term1 = po4[i] + branch->kpo41;
        double po4_term2 = po4[i] + branch->kpo42;

        double nlim1 = 0.0;
        if (dsi_term > 0.0 && n_term1 > 0.0 && po4_term1 > 0.0) {
            nlim1 = (dsi[i] / dsi_term) * ((no3nh4) / n_term1) * (po4[i] / po4_term1);
        }
        double nlim2 = 0.0;
        if (n_term2 > 0.0 && po4_term2 > 0.0) {
            nlim2 = (no3nh4 / n_term2) * (po4[i] / po4_term2);
        }

        double nswitch = nh4[i] / (5.0 + CGEM_MAX(nh4[i], 0.0));

        double gpp1 = pbmax1 * phy1[i] * nlim1 * lightlim1;
        double gpp2 = pbmax2 * phy2[i] * nlim2 * lightlim2;

        double maint1 = temp_kin(TEMP_KIN_KMAINT, branch->kmaint1, temp) * phy1[i];
        double maint2 = temp_kin(TEMP_KIN_KMAINT, branch->kmaint2, temp) * phy2[i];

        double npp_no31 = (1.0 - nswitch) * (gpp1 / depth) * (1.0 - branch->kexc1) * (1.0 - branch->kgrowth1) - maint1;
        double npp_no32 = (1.0 - nswitch) * (gpp2 / depth) * (1.0 - branch->kexc2) * (1.0 - branch->kgrowth2) - maint2;
        double npp_nh41 = nswitch * (gpp1 / depth) * (1.0 - branch->kexc1) * (1.0 - branch->kgrowth1) - maint1;
        double npp_nh42 = nswitch * (gpp2 / depth) * (1.0 - branch->kexc2) * (1.0 - branch->kgrowth2) - maint2;

        double mort1 = temp_kin(TEMP_KIN_KMORT, branch->kmort1, temp) * phy1[i];
        double mort2 = temp_kin(TEMP_KIN_KMORT, branch->kmort2, temp) * phy2[i];

        double si_cons = fmin(0.0, -branch->redsi * (npp_nh41 + npp_nh42));

        double aer_deg = temp_kin(TEMP_KIN_KHET, branch->kox, temp) *
                         (toc[i] / (toc[i] + branch->ktox)) *
                         (o2[i] / (o2[i] + branch->ko2));

        double denit = temp_kin(TEMP_KIN_KDENIT, branch->kdenit, temp) *
                       (toc[i] / (toc[i] + branch->ktox)) *
                       (branch->kino2 / (o2[i] + branch->kino2)) *
                       (no3[i] / (no3[i] + branch->kno3));

        double nitrif = temp_kin(TEMP_KIN_KNIT, branch->knit, temp) *
                        (o2[i] / (o2[i] + branch->ko2_nit)) *
                        (nh4[i] / (nh4[i] + branch->knh4));

        double o2_sat = oxygen_saturation(temp, sal);
        double o2_ex = (pis_vel[i] / depth) * (o2_sat - o2[i]);

        double npp_no3_sum = npp_no31 + npp_no32;
        double npp_nh4_sum = npp_nh41 + npp_nh42;
        double npp_total = npp_no3_sum + npp_nh4_sum;
        double phy_death = mort1 + mort2;

        branch->reaction_rates[CGEM_REACTION_GPP_1][i] = gpp1;
        branch->reaction_rates[CGEM_REACTION_GPP_2][i] = gpp2;
        branch->reaction_rates[CGEM_REACTION_NPP_NO3_1][i] = npp_no31;
        branch->reaction_rates[CGEM_REACTION_NPP_NO3_2][i] = npp_no32;
        branch->reaction_rates[CGEM_REACTION_NPP_NH4_1][i] = npp_nh41;
        branch->reaction_rates[CGEM_REACTION_NPP_NH4_2][i] = npp_nh42;
        branch->reaction_rates[CGEM_REACTION_NPP_NO3][i] = npp_no3_sum;
        branch->reaction_rates[CGEM_REACTION_NPP_NH4][i] = npp_nh4_sum;
        branch->reaction_rates[CGEM_REACTION_NPP][i] = npp_total;
        branch->reaction_rates[CGEM_REACTION_PHY_DEATH_1][i] = mort1;
        branch->reaction_rates[CGEM_REACTION_PHY_DEATH_2][i] = mort2;
        branch->reaction_rates[CGEM_REACTION_PHY_DEATH][i] = phy_death;
        branch->reaction_rates[CGEM_REACTION_SI_CONS][i] = si_cons;
        branch->reaction_rates[CGEM_REACTION_AER_DEG][i] = aer_deg;
        branch->reaction_rates[CGEM_REACTION_DENIT][i] = denit;
        branch->reaction_rates[CGEM_REACTION_NIT][i] = nitrif;
        branch->reaction_rates[CGEM_REACTION_O2_EX][i] = o2_ex;
        branch->reaction_rates[CGEM_REACTION_O2_EX_S][i] = o2_ex * depth;
    }

    /* Second loop: carbonate system and state updates */
    for (int i = 1; i <= M - 1; i += 2) {
        double depth = CGEM_MAX(branch->depth[i], CGEM_MIN_DEPTH);
        double sal = salinity[i];

        double henry = calculate_carbonate_system(branch, i, temp, depth, sal, pis_vel[i]);

        double npp_no31 = branch->reaction_rates[CGEM_REACTION_NPP_NO3_1][i];
        double npp_no32 = branch->reaction_rates[CGEM_REACTION_NPP_NO3_2][i];
        double npp_nh41 = branch->reaction_rates[CGEM_REACTION_NPP_NH4_1][i];
        double npp_nh42 = branch->reaction_rates[CGEM_REACTION_NPP_NH4_2][i];
        double npp_no3 = branch->reaction_rates[CGEM_REACTION_NPP_NO3][i];
        double npp_nh4 = branch->reaction_rates[CGEM_REACTION_NPP_NH4][i];
        double npp_total = branch->reaction_rates[CGEM_REACTION_NPP][i];
        double mort1 = branch->reaction_rates[CGEM_REACTION_PHY_DEATH_1][i];
        double mort2 = branch->reaction_rates[CGEM_REACTION_PHY_DEATH_2][i];
        double phy_death = branch->reaction_rates[CGEM_REACTION_PHY_DEATH][i];
        double si_cons = branch->reaction_rates[CGEM_REACTION_SI_CONS][i];
        double aer_deg = branch->reaction_rates[CGEM_REACTION_AER_DEG][i];
        double denit = branch->reaction_rates[CGEM_REACTION_DENIT][i];
        double nitrif = branch->reaction_rates[CGEM_REACTION_NIT][i];
        double o2_ex = branch->reaction_rates[CGEM_REACTION_O2_EX][i];

        phy1[i] = CGEM_MAX(0.0, phy1[i] + (npp_no31 + npp_nh41 - mort1) * dt);
        phy2[i] = CGEM_MAX(0.0, phy2[i] + (npp_no32 + npp_nh42 - mort2) * dt);
        dsi[i] = CGEM_MAX(0.0, dsi[i] + si_cons * dt);
        po4[i] = CGEM_MAX(0.0, po4[i] + branch->redp * (aer_deg + denit - npp_total) * dt);
        o2[i] = CGEM_MAX(0.0, o2[i] + (-aer_deg + npp_nh4 + (138.0 / 106.0) * npp_no3 - 2.0 * nitrif + o2_ex) * dt);
        toc[i] = CGEM_MAX(0.0, toc[i] + (-aer_deg - denit + phy_death) * dt);
        nh4[i] = CGEM_MAX(0.0, nh4[i] + (branch->redn * (aer_deg - npp_nh4) - nitrif) * dt);
        no3[i] = CGEM_MAX(0.0, no3[i] + (-94.4 / 106.0 * denit + nitrif - branch->redn * npp_no3) * dt);

        dic[i] = CGEM_MAX(0.0, dic[i] + branch->reaction_rates[CGEM_REACTION_DIC_REACT][i] * dt);
        at[i] = CGEM_MAX(0.0, at[i] + branch->reaction_rates[CGEM_REACTION_TA_REACT][i] * dt);
        co2[i] = CGEM_MAX(0.0, co2[i] + (branch->reaction_rates[CGEM_REACTION_CO2_EX][i] + denit - npp_total) * dt);
        if (henry > 0.0) {
            pco2[i] = CGEM_MAX(0.0, co2[i] / henry);
        } else {
            pco2[i] = 0.0;
        }
    }

    free(pis_vel);
    return 0;
}
