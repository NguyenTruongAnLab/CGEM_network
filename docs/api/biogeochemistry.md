# Biogeochemistry API

## Main Driver

### Biogeo_Branch

Main biogeochemistry driver for a branch.

```c
int Biogeo_Branch(Branch *b, double dt);
```

Calls all sub-modules:
1. Phytoplankton growth/mortality
2. Nutrient uptake/regeneration
3. Oxygen dynamics
4. Carbonate chemistry
5. GHG production

## Phytoplankton

### rive_calc_phytoplankton

Calculate phytoplankton dynamics.

```c
void rive_calc_phytoplankton(
    Branch *b, int i,
    double *gpp1, double *gpp2,
    double *npp1, double *npp2,
    double *mort1, double *mort2
);
```

### rive_temp_kin

Temperature-dependent rate modifier.

```c
double rive_temp_kin(int type, double k20, double temp);
```

$$k(T) = k_{20} \cdot \theta^{(T-20)}$$

## Nutrients

### rive_calc_nutrients

Calculate nutrient transformations.

```c
void rive_calc_nutrients(
    Branch *b, int i,
    double *dNO3, double *dNH4, double *dPO4, double *dDSi,
    double npp1, double npp2, double mort1, double mort2
);
```

## Oxygen

### rive_calc_oxygen

Calculate oxygen dynamics.

```c
void rive_calc_oxygen(
    Branch *b, int i,
    double *dO2, double *o2_ex,
    double gpp, double resp, double nit
);
```

### calc_o2_saturation

Calculate oxygen saturation.

```c
double calc_o2_saturation(double temp, double sal);
```

## Carbonate Chemistry

### rive_calc_carbonate

Calculate carbonate system equilibrium.

```c
void rive_calc_carbonate(
    Branch *b, int i,
    double *pH, double *pCO2, double *CO2aq
);
```

### calc_co2_flux

Calculate CO2 air-water exchange.

```c
double calc_co2_flux(double pCO2, double pCO2_atm, double kw, double K0);
```

## GHG Module

### rive_calc_n2o

Calculate N2O production and emission.

```c
void rive_calc_n2o(
    Branch *b, int i,
    double nit, double denit,
    double *n2o_nit, double *n2o_denit, double *n2o_ex
);
```

### rive_calc_ch4

Calculate CH4 dynamics.

```c
void rive_calc_ch4(
    Branch *b, int i,
    double *ch4_prod, double *ch4_ox, double *ch4_ex
);
```

## Sediment

### Sediment_Branch

Calculate erosion and deposition.

```c
int Sediment_Branch(Branch *b, double dt);
```

### calc_settling_velocity

Flocculation-enhanced settling.

```c
double calc_settling_velocity(double sal, BiogeoParams *params);
```

$$w_s = w_{s,0} \cdot \left(1 + (f_{max}-1) \cdot \tanh(S/S_{scale})\right)$$

## RK4 Solver

### rk4_step

Runge-Kutta 4th order step.

```c
void rk4_step(
    double *y, int n,
    double (*f)(double*, void*), void *params,
    double dt
);
```
