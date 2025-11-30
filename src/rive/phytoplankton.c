/**
 * @file phytoplankton.c
 * @brief RIVE Phytoplankton Module Implementation
 * 
 * Reference: Billen et al. (1994), Garnier et al. (2002)
 * 
 * TROPICAL ESTUARY ENHANCEMENT (Mekong Delta):
 * - Salinity toxicity for freshwater phytoplankton (PHY1 diatoms)
 * - Freshwater algae die rapidly when salinity > 5-8 PSU
 * - Prevents unrealistic freshwater blooms in marine zones
 * 
 * Reference: Muylaert et al. (2009), Cloern & Jassby (2010)
 */

#include "phytoplankton.h"
#include "rive_params.h"
#include "../define.h"
#include <math.h>

/**
 * Temperature-dependent rate function (matches Fortran Temp_kin)
 */
double rive_temp_kin(TempKinID id, double base_rate, double temp) {
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
 * Light limitation function (depth-integrated)
 */
double rive_light_limitation(double alpha, double pbmax, double I0, double kd, double depth) {
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
 * Calculate phytoplankton growth rates
 * 
 * TROPICAL ENHANCEMENT: Salinity stress on freshwater phytoplankton
 * - PHY1 (diatoms): Freshwater species, stressed above 5-8 PSU
 * - PHY2 (green algae): More salinity tolerant
 * 
 * Reference: Muylaert et al. (2009), Lionard et al. (2005)
 */
void rive_calc_phytoplankton(Branch *branch, int i, double temp, double depth,
                             double *gpp1, double *gpp2,
                             double *npp_no3_1, double *npp_no3_2,
                             double *npp_nh4_1, double *npp_nh4_2,
                             double *mort1, double *mort2) {
    if (!branch || !branch->conc) return;
    
    /* Get global parameters for salinity stress */
    BiogeoParams *p = rive_get_params();
    
    double *phy1 = branch->conc[CGEM_SPECIES_PHY1];
    double *phy2 = branch->conc[CGEM_SPECIES_PHY2];
    double *dsi = branch->conc[CGEM_SPECIES_DSI];
    double *no3 = branch->conc[CGEM_SPECIES_NO3];
    double *nh4 = branch->conc[CGEM_SPECIES_NH4];
    double *po4 = branch->conc[CGEM_SPECIES_PO4];
    double *spm = branch->conc[CGEM_SPECIES_SPM];
    double *salinity = branch->conc[CGEM_SPECIES_SALINITY];
    
    /* Get local salinity */
    double sal = (salinity) ? salinity[i] : 0.0;
    if (sal < 0.0) sal = 0.0;
    
    /* Light attenuation */
    double kd = branch->kd1 +
                branch->kd2_spm * (spm ? spm[i] * 1e-3 : 0.0) +
                branch->kd2_phy1 * phy1[i] +
                branch->kd2_phy2 * phy2[i];
    
    /* Temperature-corrected max photosynthesis */
    double pbmax1 = rive_temp_kin(TEMP_KIN_PBMAX, branch->pbmax1, temp);
    double pbmax2 = rive_temp_kin(TEMP_KIN_PBMAX, branch->pbmax2, temp);
    
    /* Light limitation */
    double lightlim1 = rive_light_limitation(branch->alpha1, pbmax1, branch->I0, kd, depth);
    double lightlim2 = rive_light_limitation(branch->alpha2, pbmax2, branch->I0, kd, depth);
    
    /* Nutrient concentrations */
    double no3nh4 = no3[i] + nh4[i];
    
    /* Nutrient limitation terms */
    double dsi_term = dsi[i] + branch->kdsi1;
    double n_term1 = no3nh4 + branch->kn1;
    double n_term2 = no3nh4 + branch->kn2;
    double po4_term1 = po4[i] + branch->kpo41;
    double po4_term2 = po4[i] + branch->kpo42;
    
    /* Nutrient limitation (Monod) */
    double nlim1 = 0.0;
    if (dsi_term > 0.0 && n_term1 > 0.0 && po4_term1 > 0.0) {
        nlim1 = (dsi[i] / dsi_term) * ((no3nh4) / n_term1) * (po4[i] / po4_term1);
    }
    double nlim2 = 0.0;
    if (n_term2 > 0.0 && po4_term2 > 0.0) {
        nlim2 = (no3nh4 / n_term2) * (po4[i] / po4_term2);
    }
    
    /* =======================================================================
     * SALINITY STRESS ON PHYTOPLANKTON (Tropical Estuary Enhancement)
     * 
     * Freshwater diatoms (PHY1) are killed by high salinity.
     * The mortality increases linearly above a threshold (typically 5 PSU).
     * 
     * f_sal_phy1 = 1.0 + coef * max(0, sal - threshold)
     * 
     * For sal = 0-5 PSU: f_sal = 1.0 (no stress)
     * For sal = 10 PSU:  f_sal = 1.0 + 0.5 * 5 = 3.5 (3.5x mortality)
     * For sal = 20 PSU:  f_sal = 1.0 + 0.5 * 15 = 8.5 (rapid death)
     * 
     * PHY2 (green algae) are more euryhaline but still stressed at extremes.
     * 
     * Reference: Muylaert et al. (2009), Lionard et al. (2005)
     * =======================================================================*/
    double sal_thresh = (p && p->sal_stress_thresh > 0.0) ? p->sal_stress_thresh : 5.0;
    double sal_coef = (p && p->sal_stress_coef > 0.0) ? p->sal_stress_coef : 0.5;
    
    /* PHY1 (freshwater diatoms): Strong salinity stress */
    double f_sal_phy1 = 1.0;
    if (sal > sal_thresh) {
        f_sal_phy1 = 1.0 + sal_coef * (sal - sal_thresh);
    }
    
    /* PHY2 (green algae): Weaker stress, more tolerant */
    double f_sal_phy2 = 1.0;
    if (sal > sal_thresh * 2.0) {  /* Higher threshold for greens */
        f_sal_phy2 = 1.0 + (sal_coef * 0.3) * (sal - sal_thresh * 2.0);
    }
    
    /* N source preference (NH4 preferred) */
    double nswitch = nh4[i] / (5.0 + CGEM_MAX(nh4[i], 0.0));
    
    /* Gross primary production */
    *gpp1 = pbmax1 * phy1[i] * nlim1 * lightlim1;
    *gpp2 = pbmax2 * phy2[i] * nlim2 * lightlim2;
    
    /* Maintenance respiration */
    double maint1 = rive_temp_kin(TEMP_KIN_KMAINT, branch->kmaint1, temp) * phy1[i];
    double maint2 = rive_temp_kin(TEMP_KIN_KMAINT, branch->kmaint2, temp) * phy2[i];
    
    /* Net primary production by N source */
    *npp_no3_1 = (1.0 - nswitch) * (*gpp1 / depth) * (1.0 - branch->kexc1) * (1.0 - branch->kgrowth1) - maint1;
    *npp_no3_2 = (1.0 - nswitch) * (*gpp2 / depth) * (1.0 - branch->kexc2) * (1.0 - branch->kgrowth2) - maint2;
    *npp_nh4_1 = nswitch * (*gpp1 / depth) * (1.0 - branch->kexc1) * (1.0 - branch->kgrowth1) - maint1;
    *npp_nh4_2 = nswitch * (*gpp2 / depth) * (1.0 - branch->kexc2) * (1.0 - branch->kgrowth2) - maint2;
    
    /* Mortality with salinity stress factor */
    *mort1 = rive_temp_kin(TEMP_KIN_KMORT, branch->kmort1, temp) * phy1[i] * f_sal_phy1;
    *mort2 = rive_temp_kin(TEMP_KIN_KMORT, branch->kmort2, temp) * phy2[i] * f_sal_phy2;
}
