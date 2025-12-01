# Hydrodynamics API

## Main Functions

### Hyd_Branch

Solve Saint-Venant equations for a single branch.

```c
int Hyd_Branch(Branch *b, double H_down, double H_up, double Q_up, double dt);
```

**Parameters:**
- `b`: Branch pointer
- `H_down`: Water level at downstream boundary [m]
- `H_up`: Water level at upstream boundary [m]
- `Q_up`: Discharge at upstream boundary [m³/s]
- `dt`: Time step [s]

**Returns:** 0 on success, -1 on error

### solve_network_hydro

Solve hydrodynamics for entire network with junction iteration.

```c
int solve_network_hydro(Network *net, double dt, int max_iter, double tol);
```

**Parameters:**
- `net`: Network pointer
- `dt`: Time step [s]
- `max_iter`: Maximum junction iterations
- `tol`: Convergence tolerance

**Returns:** 0 on success, number of non-converged junctions on failure

## Geometry Functions

### calc_width

Calculate channel width at position x.

```c
double calc_width(Branch *b, double x);
```

Exponential geometry:
$$B(x) = B_0 \exp(-x/L_c)$$

### calc_area

Calculate cross-sectional area.

```c
double calc_area(Branch *b, int i);
```

$$A = B \cdot H$$

### calc_hydraulic_radius

Calculate hydraulic radius.

```c
double calc_hydraulic_radius(Branch *b, int i);
```

$$R = A / B \approx H$$ (for wide channels)

## Friction

### calc_friction_slope

Calculate friction slope using Chezy or Manning.

```c
double calc_friction_slope(Branch *b, int i);
```

Chezy:
$$S_f = \frac{U|U|}{C^2 R}$$

Manning:
$$S_f = \frac{n^2 U|U|}{R^{4/3}}$$

## Solver

### solve_tridiagonal

Thomas algorithm for tridiagonal systems.

```c
void solve_tridiagonal(double *a, double *b, double *c, double *d, double *x, int n);
```

Solves: $a_i x_{i-1} + b_i x_i + c_i x_{i+1} = d_i$
