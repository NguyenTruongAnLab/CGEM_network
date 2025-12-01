# Transport API

## Main Functions

### Transport_Branch

Solve advection-dispersion for a species on one branch.

```c
int Transport_Branch(Branch *b, int species, double dt, double C_up, double C_down);
```

**Parameters:**
- `b`: Branch pointer
- `species`: Species index (from `define.h`)
- `dt`: Time step [s]
- `C_up`: Upstream boundary concentration
- `C_down`: Downstream boundary concentration

**Returns:** 0 on success, -1 on error

### transport_all_species

Transport all species for a branch.

```c
int transport_all_species(Branch *b, double dt, double *C_up, double *C_down, int num_species);
```

## TVD Schemes

### calc_flux_limiter

Calculate flux limiter for TVD schemes.

```c
double calc_flux_limiter(double r, int limiter_type);
```

**Limiter types:**
- 0: Minmod
- 1: Superbee (default)
- 2: Van Leer
- 3: MC

Superbee:
$$\phi(r) = \max(0, \min(2r, 1), \min(r, 2))$$

### calc_advective_flux

Calculate advective flux with TVD limiter.

```c
double calc_advective_flux(Branch *b, int species, int i, int limiter);
```

## Dispersion

### calc_dispersion_coefficient

Calculate Van den Burgh dispersion.

```c
double calc_dispersion_coefficient(Branch *b, int i);
```

$$K = K_0 + K_{VDB} \cdot U \cdot A$$

### update_dispersion_field

Update dispersion for entire branch.

```c
void update_dispersion_field(Branch *b);
```
