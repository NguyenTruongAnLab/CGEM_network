#ifndef CGEM_NETWORK_H
#define CGEM_NETWORK_H

#include <stddef.h>
#include <stdio.h>
#include "define.h"

#define CGEM_MAX_NODE_BRANCHES 8

typedef enum {
    NODE_JUNCTION = 0,
    NODE_DISCHARGE_BC = 1,
    NODE_LEVEL_BC = 2
} NodeType;

typedef struct {
    int id;
    NodeType type;
    double H;
    double H_new;
    double Q_net;

    int connected_branches[CGEM_MAX_NODE_BRANCHES];
    int connection_dir[CGEM_MAX_NODE_BRANCHES];
    int num_connections;

    /* Time-varying hydrodynamic forcing */
    double *forcing_time;
    double *forcing_value;
    size_t forcing_len;

    /* Junction concentration mixing */
    double *mixed_conc;
    
    /* Time-varying species boundary conditions */
    /* species_forcing_time[sp] = array of times for species sp */
    /* species_forcing_value[sp] = array of values for species sp */
    double **species_forcing_time;
    double **species_forcing_value;
    size_t *species_forcing_len;
    int num_species_forcing;  /* Number of species with time-varying BC */
} Node;

typedef struct {
    int id;
    char name[CGEM_MAX_BRANCH_NAME];
    int node_up;
    int node_down;
    NodeType up_node_type;
    NodeType down_node_type;

    double length_m;
    double width_up_m;
    double width_down_m;
    double depth_up_m;      /* Depth at downstream (ocean side) */
    double depth_down_m;    /* Depth at upstream */
    double depth_m;         /* Reference depth (mean) */
    double storage_ratio;   /* RS: Storage width ratio (1.0 for prismatic, 5-10 for mangroves) */
    int group_id;

    double chezy;
    double lc_convergence;

    int M;
    double dx;

    /* Hydrodynamic arrays [0..M+1] - staggered grid */
    double *totalArea;      /* AA: Total cross-sectional area */
    double *freeArea;       /* AAf: Free cross-section (above reference) */
    double *refArea;        /* AAref: Reference cross-section */
    double *width;          /* B: Width at each point */
    double *depth;          /* H: Water depth */
    double *velocity;       /* U: Current velocity */
    double *waterLevel;     /* Water surface elevation */
    double *dispersion;     /* K: Dispersion coefficient */
    double *chezyArray;     /* Chezy coefficient (spatially varying) */

    /* Previous timestep values (for iteration) */
    double *totalArea_old;
    double *totalArea_old2;
    
    /* =========================================================================
     * RESIDUAL VELOCITY FILTER (Critical for Van den Burgh dispersion)
     * Low-pass filtered velocity to extract tidal-mean residual flow.
     * Reference: Savenije (2005), Audit recommendation for Q_f calculation
     * =========================================================================*/
    double *u_residual;         /* Tidally-averaged velocity [m/s] */
    double residual_alpha;      /* Low-pass filter coefficient (0.001 typical) */
    
    /* Manning's n friction (optional, replaces Chezy) */
    double *manning_n;          /* Manning coefficient at each point [s/m^{1/3}] */
    int use_manning;            /* 1 = Use Manning's n instead of Chezy */
    
    /* Storage width ratio dynamics (floodplain effects) */
    double H_bank;              /* Bank-full water level [m] */
    double RS_channel;          /* Storage ratio for in-channel flow */
    double RS_floodplain;       /* Storage ratio when flooded (5-15 typical) */

    /* Concentration arrays [num_species][0..M+1] */
    double **conc;

    /* Tridiagonal solver work arrays */
    double *tri_lower;
    double *tri_diag;
    double *tri_upper;
    double *tri_rhs;
    double *tri_gam;

    /* Transport parameters */
    double D0;              /* Dispersion at mouth */
    double vdb_coef;        /* Van den Burgh coefficient (tuning) */
    double mixing_alpha;    /* Mixing efficiency for D0 = α * U_tidal * B (0.1-1.0) */

    int num_species;

    /* Boundary concentrations */
    double *conc_down;      /* Ocean-side concentrations [num_species] */
    double *conc_up;        /* River-side concentrations [num_species] */

    int has_biogeo;
    char biogeo_params_path[CGEM_MAX_PATH]; /* Path to branch-specific biogeo params (empty=use global) */

    /* Sediment transport parameters */
    double *tau_ero;           /* Erosion shear stress threshold [Pa] */
    double *tau_dep;           /* Deposition shear stress threshold [Pa] */
    double *mero;              /* Erosion coefficient [kg/m²/s] */
    double ws;                 /* Settling velocity [m/s] */
    
    /* Biogeochemical parameters */
    double water_temp;         /* Water temperature [°C] */
    double salinity;           /* Salinity for calculations */
    
    /* Light parameters */
    double I0;                 /* Surface irradiance [W/m²] */
    double kd1, kd2_spm, kd2_phy1, kd2_phy2; /* Light attenuation coefficients */
    double alpha1, alpha2;     /* Initial slopes for photosynthesis */
    
    /* Phytoplankton parameters */
    double pbmax1, pbmax2;     /* Maximum photosynthetic rates */
    double kexc1, kexc2;       /* Excretion fractions */
    double kgrowth1, kgrowth2; /* Growth respiration fractions */
    double kmaint1, kmaint2;   /* Maintenance respiration rates */
    double kmort1, kmort2;     /* Mortality rates */
    
    /* Nutrient limitation parameters */
    double kdsi1, kn1, kpo41;  /* Half-saturation for Phy1 */
    double kn2, kpo42;         /* Half-saturation for Phy2 */
    
    /* Organic matter parameters */
    double kox, kdenit, knit;  /* Decomposition rates */
    double ktox, ko2, ko2_nit, kno3, knh4, kino2; /* Half-saturations */
    double redn, redp, redsi;  /* Stoichiometric ratios */
    
    /* =======================================================================
     * RIVE-enhanced biogeochemistry parameters (from QualNET)
     * Reference: Billen et al. (1994), Garnier et al. (2002)
     * =======================================================================*/
    
    /* RIVE multi-pool organic matter hydrolysis rates [1/day] */
    double khydr1;             /* Labile HD1/HP1 hydrolysis rate (0.5-1.0) */
    double khydr2;             /* Semi-labile HD2/HP2 hydrolysis rate (0.1-0.3) */
    double khydr3;             /* Refractory HD3/HP3 hydrolysis rate (0.01-0.05) */
    double frac_hd1;           /* Fraction of OC input as HD1 (labile) */
    double frac_hd2;           /* Fraction of OC input as HD2 (semi-labile) */
    double frac_hp1;           /* Fraction of particulate OC as HP1 */
    double frac_hp2;           /* Fraction of particulate OC as HP2 */
    
    /* RIVE bacteria parameters (from bacteria.cpp) */
    double bag_bmax20;         /* BAG max growth rate at 20°C [1/h] (0.6) */
    double bap_bmax20;         /* BAP max growth rate at 20°C [1/h] (0.16) */
    double bag_kdb20;          /* BAG mortality rate at 20°C [1/h] (0.05) */
    double bap_kdb20;          /* BAP mortality rate at 20°C [1/h] (0.02) */
    double bac_ks;             /* Bacteria half-saturation for DSS [mg C/L] (0.1) */
    double bac_yield;          /* Bacteria growth yield [-] (0.25) */
    double bag_topt;           /* BAG optimal temperature [°C] (22) */
    double bap_topt;           /* BAP optimal temperature [°C] (20) */
    double bag_dti;            /* BAG temperature spread [°C] (12) */
    double bap_dti;            /* BAP temperature spread [°C] (17) */
    double bag_vs;             /* BAG settling velocity [m/h] (0.02) */
    
    /* RIVE phosphorus adsorption (from phosphorus.cpp) */
    double kpads;              /* P adsorption equilibrium constant (3.43) */
    double pac;                /* P adsorption capacity on SPM [µmol P/mg SPM] (0.37) */
    
    /* RIVE benthic flux parameters */
    double zf_init;            /* Initial fluid sediment layer depth [m] (0.001) */
    double benthic_porosity;   /* Sediment porosity [-] (0.88) */
    double benthic_density;    /* Sediment density [g/m³] (2300000) */
    
    /* Benthic organic matter pools [mg C/m²] */
    double *benthic_hb1;       /* Benthic labile OC pool */
    double *benthic_hb2;       /* Benthic semi-labile OC pool */
    double *benthic_hb3;       /* Benthic refractory OC pool */
    double *benthic_zf;        /* Fluid sediment layer depth [m] */
    
    /* =======================================================================
     * SEDIMENT BED LAYER (Fluid Mud Tracking)
     * Tracks accumulated sediment and organic carbon for benthic-pelagic coupling
     * Reference: Winterwerp (2002), Wang et al. (2018)
     * =======================================================================*/
    double *bed_mass;          /* Deposited sediment mass [kg/m²] per grid cell */
    double *bed_oc;            /* Organic carbon in bed [kg C/m²] per grid cell */
    
    /* Gas exchange */
    double pco2_atm;           /* Atmospheric pCO2 [µatm] */
    
    /* =======================================================================
     * LATERAL LOADS (Land-Use Coupling)
     * Spatially explicit loads from urban, agriculture, aquaculture, mangroves
     * Reference: Garnier et al. (2005), Alongi (2014)
     * =======================================================================*/
    double *lateral_flow;      /* Lateral inflow [m³/s] per grid cell */
    double **lateral_conc;     /* Lateral concentrations [NSPECIES][M+2] [µmol/L or mg/L] */
    int has_lateral_loads;     /* Flag: 1 if lateral loads are active */
    
    /* Reaction arrays for output */
    double **reaction_rates;   /* [num_reactions][M+2] */
    int num_reactions;

    /* Binary output file handle */
    FILE *bin_fp;

    /* CSV output file handles (hydro + species + reactions) */
    FILE *csv_fps[60];
    int num_csv_fps;
} Branch;

typedef struct {
    Branch **branches;
    size_t num_branches;
    Node *nodes;
    size_t num_nodes;
    int num_species;
    
    /* Simulation parameters */
    double dt;              /* Time step [s] */
    double dx_target;       /* Target grid spacing [m] (from case config DELXI) */
    double total_time;      /* Total simulation time [s] */
    double warmup_time;     /* Current simulation time [s] */
    double current_time;    /* Current simulation time [s] */
    
    /* Tidal parameters */
    double tidal_amplitude; /* [m] */
    double tidal_period;    /* [s] */
    
    /* River discharge */
    double Q_river;         /* Upstream discharge [m³/s] */
} Network;

typedef struct {
    char case_name[128];
    char config_path[CGEM_MAX_PATH];
    char base_dir[CGEM_MAX_PATH];
    char topology_path[CGEM_MAX_PATH];
    char boundary_path[CGEM_MAX_PATH];
    char params_override_path[CGEM_MAX_PATH];
    char global_defaults_path[CGEM_MAX_PATH];
    char output_dir[CGEM_MAX_PATH];
    char start_date[32];
    int duration_days;
    int warmup_days;
    int write_csv;
    int write_netcdf;
    int write_reaction_rates;  /* Output reaction rates to CSV (for debugging pCO2) */
    double dt;              /* Time step [s] */
    double dx_meters;       /* Grid spacing [m] specified in case config (DELXI) */
    double tidal_amplitude; /* [m] */
    double Q_river;         /* [m³/s] */
} CaseConfig;

/* Memory management */
Branch *allocate_branch(int M, int num_species);
void free_branch(Branch *branch, int num_species);
Node *allocate_node(int num_species);
void free_node(Node *node);

/* I/O functions */
int LoadCaseConfig(const char *path, CaseConfig *config);
int LoadTopology(const char *path, Network *net);
int LoadBoundaries(const char *path, Network *net);
int cgem_write_netcdf(Network *net, CaseConfig *config);

/* Geometry initialization */
void InitializeBranchGeometry(Branch *branch, double target_dx);

/* Hydrodynamics */
int Hyd_Branch(Branch *branch, double H_down, double H_up, double Q_up, double dt);
double TidalElevation(double t, double amplitude, double period);
double TotalDischarge(Branch *branch, int loc, double Q_river);

/* Transport */
int Transport_Branch(Branch *branch, double dt);
int Transport_Branch_Network(Branch *branch, double dt, void *network_ptr);
void ComputeDispersionCoefficient(Branch *branch, double Q_total);
double ComputeJunctionDispersion(Branch *branch, void *network_ptr);

/* Biogeochemistry */
int LoadBiogeoParams(const char *path);
void InitializeBiogeoParameters(Branch *branch);
int Biogeo_Branch(Branch *branch, double dt);

/* Lateral Loads (Land-Use Coupling) */
int LoadLateralSources(Network *net, const char *case_dir);
void ApplyLateralLoads(Branch *branch, double dt);

/* C-RIVE Enhanced Biogeochemistry (GHG and RK4 solver) */
int Biogeo_GHG_Branch(Branch *branch, double dt);
int Biogeo_Branch_RK4(Branch *branch, double dt);
int calculate_carbonate_crive(Branch *branch, int idx, double dt);

/* Initialization */
int initializeNetwork(Network *net, CaseConfig *config);

/* Network solver */
int solve_network_step(Network *net, double current_time_seconds);
int network_run_simulation(Network *net, CaseConfig *config);

/* Output */
void write_branch_output(FILE *fp, Branch *branch, double time);
int cgem_init_output(Network *net, CaseConfig *config);
int cgem_write_timestep(Network *net, CaseConfig *config, double time_s);
int cgem_close_output(Network *net);
// Removed: write_branch_matrix_csv is no longer used; matrix CSV outputs are disabled.

#endif /* CGEM_NETWORK_H */
