# SLAM.jl

[![Julia](https://img.shields.io/badge/Julia-1.6+-blue.svg)](https://julialang.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Version](https://img.shields.io/badge/version-0.1.0-green.svg)]()

**Scattering Lagrange Asymptotic Matching** - A Julia package for solving nuclear scattering problems using Lagrange mesh methods with Baye's exact differential matrices.

---

## Table of Contents

- [Overview](#overview)
- [Key Features](#key-features)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Usage Examples](#usage-examples)
- [Graphical User Interface](#graphical-user-interface)
- [Mathematical Background](#mathematical-background)
- [Project Structure](#project-structure)
- [API Reference](#api-reference)
- [Performance and Validation](#performance-and-validation)
- [Testing](#testing)
- [Physical Constants](#physical-constants)
- [References](#references)
- [Contributing](#contributing)
- [FAQ](#faq)
- [License](#license)
- [Author](#author)

---

## Overview

SLAM.jl computes scattering cross-sections and S-matrix elements for nuclear reactions (e.g., proton-nucleus, neutron-nucleus collisions) by solving the radial Schrödinger equation with Coulomb interactions. The package implements **Baye's Lagrange-Legendre mesh method** with exact differential matrices operating on expansion coefficients, providing high accuracy with minimal computational cost.

### Why SLAM.jl?

- **High Accuracy**: Uses analytical D and T matrices that provide exact derivative operations in coefficient space
- **Efficiency**: Requires only N ~ 30-60 mesh points for typical calculations (vs. thousands for finite difference methods)
- **Robustness**: Handles both charged (p, d, α...) and neutral (n) particle scattering
- **Validated**: Results verified against COLOSS Fortran code with ~10⁻⁵ precision match
- **User-Friendly**: Includes interactive GUI for real-time calculations and visualization

---

## Key Features

### Core Numerical Methods

| Feature | Description |
|---------|-------------|
| **Lagrange-Legendre mesh** | Gauss-Legendre quadrature points on (0, R) with optimal weight distribution |
| **Baye's exact D matrix** | First derivative matrix with analytical formulas for coefficient-space operations |
| **Baye's exact T matrix** | Kinetic energy matrix (-d²/dr²) derived analytically |
| **x-regularized basis** | Basis functions f_j(r) = (r/r_j) · L_j(r) / √λ_j satisfying correct r→0 behavior |
| **Coordinate-transformed mesh** | Optional power-law or FRESCO-style transformations for improved point distribution |

### Physics Components

| Feature | Description |
|---------|-------------|
| **Coulomb wave functions** | F_l, G_l, H⁺_l via COUL90 Fortran library or pure Julia fallback |
| **Woods-Saxon potentials** | Volume, surface (derivative form), and spin-orbit terms |
| **Finite-size Coulomb** | Uniformly charged sphere model for r < R_c |
| **Koning-Delaroche (KD02)** | Global optical model parameters for nucleons 1 keV - 200 MeV |
| **Outgoing wave BC** | Coulomb-Hankel function H⁺_l for proper asymptotic matching |

### Analysis Tools

| Feature | Description |
|---------|-------------|
| **S-matrix extraction** | Complex scattering amplitude from boundary matching |
| **Phase shifts** | Real and complex phase shifts in degrees |
| **Cross sections** | Partial wave, total reaction, and differential cross sections |
| **Wave functions** | Full radial wave function on mesh with interpolation |
| **Numerov validation** | Built-in Numerov integrator for method comparison |

### User Interface

| Feature | Description |
|---------|-------------|
| **Interactive GUI** | Web-based interface with parameter sliders |
| **Argand diagrams** | S-matrix visualization in complex plane |
| **Cross section plots** | Angular distributions dσ/dΩ vs θ |
| **Wave function display** | Real and imaginary parts vs r |

---

## Installation

### Requirements

- **Julia 1.6 or later** (1.9+ recommended for best performance)
- Operating System: macOS, Linux, or Windows

### From GitHub

```julia
using Pkg
Pkg.add(url="https://github.com/jinlei/SLAM.jl")
```

### Local Development

```bash
# Clone the repository
git clone https://github.com/jinlei/SLAM.jl
cd SLAM.jl

# Install dependencies
julia --project=. -e 'using Pkg; Pkg.instantiate()'

# (Optional) Run tests
julia --project=. -e 'using Pkg; Pkg.test()'
```

### Building the Coulomb Library (Optional)

For highest accuracy Coulomb functions, compile the COUL90 Fortran library:

```bash
cd deps
gfortran -shared -fPIC -O3 -o libcoul90.dylib coul90.f90  # macOS
# or
gfortran -shared -fPIC -O3 -o libcoul90.so coul90.f90     # Linux
```

If the library is not found, SLAM.jl automatically falls back to a pure Julia implementation using continued fractions.

### Dependencies

| Package | Purpose |
|---------|---------|
| [FastGaussQuadrature.jl](https://github.com/JuliaApproximation/FastGaussQuadrature.jl) | Gauss-Legendre nodes and weights |
| [SpecialFunctions.jl](https://github.com/JuliaMath/SpecialFunctions.jl) | Gamma functions for Coulomb phase shifts |
| [LinearAlgebra](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/) | Matrix operations and linear solvers |
| [Blink.jl](https://github.com/JuliaGizmos/Blink.jl) | Interactive GUI (optional) |
| [JSON3.jl](https://github.com/quinnj/JSON3.jl) | GUI data communication |
| [LaTeXStrings.jl](https://github.com/JuliaStrings/LaTeXStrings.jl) | Plot labels |

---

## Quick Start

```julia
using SLAM

# Define optical potential for p + ⁴⁰Ca at 30 MeV
pot = OpticalPotential(
    V_v=50.0, r_v=1.25, a_v=0.65,    # Real volume Woods-Saxon
    W_v=10.0, r_wv=1.25, a_wv=0.65,  # Imaginary volume (absorption)
    r_c=1.25,                         # Coulomb radius parameter
    Z_proj=1.0, A_proj=1.0,          # Proton
    Z_targ=20.0, A_targ=40.0         # ⁴⁰Ca target
)

# Create scattering problem: E_cm = 30 MeV, l = 0
prob = ScatteringProblem(pot, 30.0, 0)

# Solve for S-matrix and observables
result = solve_scattering(prob)

println("S-matrix:     ", result.S_matrix)
println("Phase shift:  ", result.phase_shift, " degrees")
println("|S|:          ", abs(result.S_matrix))
println("Cross section:", result.cross_section, " mb")
```

**Expected output:**
```
S-matrix:     0.4521 + 0.3142im
Phase shift:  17.35 degrees
|S|:          0.5507
Cross section: 28.55 mb
```

---

## Usage Examples

### Neutron Scattering (No Coulomb)

```julia
# Neutron on ⁴⁰Ca at 20 MeV lab energy
E_lab = 20.0
A_proj, A_targ = 1.0, 40.0
E_cm = E_lab * A_targ / (A_proj + A_targ)  # Convert to CM frame

pot = OpticalPotential(
    V_v=46.553, r_v=1.185, a_v=0.672,
    W_v=1.777, r_wv=1.185, a_wv=0.672,
    W_s=7.182, r_ws=1.288, a_ws=0.538,   # Surface absorption
    V_so=5.343, r_so=0.996, a_so=0.590,  # Spin-orbit
    r_c=0.0,                              # No Coulomb
    Z_proj=0.0, A_proj=A_proj,
    Z_targ=0.0, A_targ=A_targ
)

prob = ScatteringProblem(pot, E_cm, 0; j=0.5, N=60, R=20.0)
result = solve_scattering(prob)
```

### Sum Over Partial Waves

Calculate total reaction cross section by summing all partial waves:

```julia
results = solve_all_partial_waves(pot, 30.0; l_max=20, N=60, R=20.0)

# Total reaction cross section
σ_total = sum(r.cross_section for r in results)
println("Total reaction cross section: $σ_total mb")

# Print partial wave contributions
for r in results
    println("l=$(r.l), j=$(r.j): |S|=$(round(abs(r.S_matrix), digits=4)), σ=$(round(r.cross_section, digits=2)) mb")
end
```

### Differential Cross Section

Compute the elastic angular distribution:

```julia
angles = 0:5:180  # degrees
θ_rad = angles .* (π/180)

dσ_dΩ = [elastic_differential_cross_section(pot, 30.0, θ; l_max=25) for θ in θ_rad]

# Plot (requires Plots.jl)
# using Plots
# plot(angles, dσ_dΩ, yscale=:log10, xlabel="θ (deg)", ylabel="dσ/dΩ (mb/sr)")
```

### Wave Function Extraction

Obtain the full radial wave function:

```julia
wf = solve_wavefunction(prob)

# Access data at mesh points
r_mesh = wf.r           # Radial coordinates (fm)
ψ_sc = wf.ψ_sc          # Scattered wave function
ψ_total = wf.ψ_total    # Total: e^{iσ_l}[F_l + ψ_sc]
F_l = wf.F_l            # Regular Coulomb function
c = wf.c                # Lagrange expansion coefficients

# Interpolate to arbitrary radius
ψ_at_5fm = evaluate_wavefunction(wf, 5.0)
```

### Numerov Method Comparison

Validate Lagrange results against standard Numerov integration:

```julia
# Lagrange solution
result_lag = solve_scattering(prob)

# Numerov solution (built-in)
result_num = solve_numerov(prob; h=0.02)  # Step size 0.02 fm

println("Lagrange S: ", result_lag.S_matrix)
println("Numerov S:  ", result_num.S_matrix)
println("Difference: ", abs(result_lag.S_matrix - result_num.S_matrix))
```

### Koning-Delaroche Global Optical Model

Use the KD02 energy-dependent global parameterization:

```julia
include("src/kd02.jl")
using .KD02

# Get KD02 parameters for p + ¹²C at 30 MeV
# k0=2 for protons, k0=1 for neutrons
params = kd02_potential(30.0, 12, 6, 2)

# Create OpticalPotential from KD02 parameters
pot = OpticalPotential(
    V_v=params.v, r_v=params.rv, a_v=params.av,
    W_v=params.w, r_wv=params.rw, a_wv=params.aw,
    W_s=params.wd, r_ws=params.rwd, a_ws=params.awd,
    V_so=params.vso, r_so=params.rvso, a_so=params.avso,
    W_so=params.wso, r_wso=params.rwso, a_wso=params.awso,
    r_c=params.rc,
    Z_proj=1.0, A_proj=1.0,
    Z_targ=6.0, A_targ=12.0
)
```

### Coordinate-Transformed Mesh

Use transformed coordinates to concentrate mesh points in the interaction region:

```julia
# Standard Legendre mesh
mesh_std = init_legendre_mesh(40, 30.0)

# Power transformation (x^2): more points near r=0
mesh_pow = init_transformed_mesh(40, 30.0, transform=TRANSFORM_POWER, power=2.0)

# FRESCO-style transformation
mesh_fre = init_transformed_mesh(40, 30.0, transform=TRANSFORM_FRESCO)

# Compare point distributions
mesh_point_distribution(mesh_std)
mesh_point_distribution(mesh_pow)
```

---

## Graphical User Interface

SLAM.jl includes an interactive web-based GUI for real-time nuclear scattering calculations.

### Features

- **Parameter Sliders**: Adjust optical potential depths, radii, and diffuseness in real-time
- **Argand Diagrams**: Visualize S-matrix elements in the complex plane
- **Cross Section Plots**: Angular distributions dσ/dΩ vs scattering angle
- **Wave Function Display**: Plot |ψ(r)|², Re(ψ), Im(ψ) vs radial coordinate
- **Method Comparison**: Side-by-side Lagrange vs Numerov results

### Quick Start (Recommended)

Use the provided shell scripts:

```bash
# First time: install Julia and dependencies
./setup.sh

# Launch the GUI
./run_gui.sh
```

The GUI will open in your default web browser.

### Shell Scripts

| Script | Description |
|--------|-------------|
| `setup.sh` | **One-time setup.** Checks for Julia installation; if not found, installs via Homebrew (macOS) or juliaup (Linux). Then installs all required Julia packages. |
| `run_gui.sh` | **GUI launcher.** Verifies Julia is installed (runs `setup.sh` if needed), then launches the interactive GUI in your browser. |

### Manual Launch

```bash
# Command line
julia --project=. scripts/gui.jl
```

Or from within Julia:

```julia
include("scripts/gui.jl")
```

### Platform Support

| Platform | Installation Method |
|----------|-------------------|
| **macOS** | Homebrew or juliaup |
| **Linux** | juliaup (via curl or wget) |
| **Windows** | Manual installation from [julialang.org/downloads](https://julialang.org/downloads/) |

---

## Mathematical Background

### Problem Formulation

SLAM.jl solves the radial Schrödinger equation in the center-of-mass frame:

$$\left[-\frac{\hbar^2}{2\mu}\frac{d^2}{dr^2} + \frac{\hbar^2 l(l+1)}{2\mu r^2} + V(r) - E\right]\psi(r) = 0$$

In dimensionless form with k² = 2μE/ℏ² and U = 2μV/ℏ²:

$$\left[-\frac{d^2}{dr^2} + \frac{l(l+1)}{r^2} + U(r) - k^2\right]\psi(r) = 0$$

### Coulomb Scattering Decomposition

For charged-particle scattering, the total wave function is decomposed as:

$$\psi(r) = F_l(\eta, kr) + \psi_{\text{sc}}(r)$$

where:
- **F_l(η, kr)** is the regular Coulomb function (incident wave)
- **η = Z₁Z₂e²μ/(ℏ²k)** is the Sommerfeld parameter
- **ψ_sc(r)** is the scattered wave satisfying the inhomogeneous equation:

$$\left[-\frac{d^2}{dr^2} + \frac{l(l+1)}{r^2} + U_{\text{full}} - k^2\right]\psi_{\text{sc}} = U_{\text{short}} \cdot F_l$$

The short-range potential is:
$$U_{\text{short}} = U_{\text{nuclear}} + U_{\text{Coulomb}}^{\text{finite}} - U_{\text{Coulomb}}^{\text{point}}$$

### Lagrange-Legendre Mesh

The method uses N Gauss-Legendre quadrature points on (0,1), scaled to (0,R):

| Quantity | Formula |
|----------|---------|
| Mesh points | r_j = R · x_j |
| Quadrature weights | λ_j = w_j / 2 |
| x-regularized basis | f_j(r) = (r/r_j) · L_j(r/R) / √λ_j |

The basis functions satisfy:
- **Lagrange property**: f_j(r_i) = δ_{ij}/√λ_j
- **Boundary condition**: f_j(0) = 0 (built into x-regularization)

The key coefficient relation is:
$$c_j = \psi(r_j) \cdot \sqrt{R \cdot \lambda_j}$$

### Baye's Exact Matrices

The differential matrices act on coefficients c_j, not function values:

**D matrix** (first derivative, scaled by R):

$$D_{ij} = (-1)^{i-j} \sqrt{\frac{x_i(1-x_j)}{x_j(1-x_i)}} \frac{1}{x_i - x_j} \quad (i \neq j)$$

$$D_{ii} = \frac{1}{2x_i(1-x_i)}$$

**T matrix** (kinetic energy -d²/dr², scaled by R²):

$$T_{ii} = \frac{N(N+1)x_i(1-x_i) - 3x_i + 1}{3x_i^2(1-x_i)^2}$$

Off-diagonal elements involve more complex analytical expressions (see Baye 2015, Eq. 203).

### Outgoing Wave Boundary Condition

At r = R, the scattered wave must behave as an outgoing wave:

$$\psi_{\text{sc}}(r) \sim A \cdot H^+_l(\eta, kr)$$

This is enforced via the Robin boundary condition:
$$\psi'_{\text{sc}}(R) - \gamma_s \cdot \psi_{\text{sc}}(R) = 0$$

where γ_s = k · H⁺'_l(kR) / H⁺_l(kR) is the logarithmic derivative of the outgoing Hankel function.

### S-Matrix Extraction

The S-matrix is extracted from the scattered wave at the boundary:

$$S_l = 1 + 2ik \cdot f_l$$

with scattering amplitude:
$$f_l = \frac{\psi_{\text{sc}}(R)}{k \cdot H^+_l(\eta, kR)}$$

### Partial Wave Cross Section

The reaction cross section for partial wave (l, j) with spin-1/2 projectile:

$$\sigma_{l,j} = \frac{\pi}{k^2} \cdot \frac{2j+1}{2s+1} \cdot (1 - |S_l|^2) \times 10 \text{ fm}^2/\text{mb}$$

---

## Project Structure

```
SLAM.jl/
├── src/
│   ├── SLAM.jl            # Module root, exports, and initialization
│   ├── mesh.jl            # Lagrange-Legendre mesh (standard + transformed)
│   ├── baye_matrices.jl   # Exact D and T matrices (analytical + numerical)
│   ├── coulomb.jl         # Coulomb functions F, G, H⁺ (COUL90 + Julia fallback)
│   ├── potentials.jl      # Woods-Saxon optical potentials
│   ├── scattering.jl      # Main scattering solver
│   ├── numerov.jl         # Numerov integration for validation
│   └── kd02.jl            # Koning-Delaroche global optical model
├── test/
│   ├── runtests.jl        # Main test suite
│   ├── compare_coloss.jl  # Validation against COLOSS Fortran
│   ├── compare_proton.jl  # Proton scattering tests
│   └── debug_*.jl         # Diagnostic and debugging scripts
├── scripts/
│   ├── gui.jl             # Interactive GUI application
│   ├── run_gui.jl         # GUI entry point
│   ├── setup.jl           # Julia package installer
│   └── compare_p12C.jl    # p + ¹²C benchmark calculation
├── deps/
│   ├── build.jl           # Build script for Fortran library
│   └── libcoul90.dylib    # Compiled COUL90 library (macOS)
├── figures/               # Generated plots and validation figures
├── setup.sh               # One-time setup script
├── run_gui.sh             # GUI launcher script
├── Project.toml           # Julia package manifest
└── README.md              # This file
```

### Source File Descriptions

| File | Lines | Description |
|------|-------|-------------|
| `SLAM.jl` | ~100 | Module definition, includes, exports |
| `mesh.jl` | ~500 | LagrangeMesh, TransformedLagrangeMesh, Lagrange polynomials |
| `baye_matrices.jl` | ~200 | Analytical D/T matrices, numerical alternatives |
| `coulomb.jl` | ~330 | Coulomb functions via COUL90 or continued fractions |
| `potentials.jl` | ~350 | OpticalPotential struct, Woods-Saxon, Coulomb |
| `scattering.jl` | ~720 | solve_scattering, solve_wavefunction, cross sections |
| `numerov.jl` | ~150 | Numerov integrator for validation |
| `kd02.jl` | ~200 | Koning-Delaroche parameterization |

---

## API Reference

### Core Types

| Type | Description |
|------|-------------|
| `LagrangeMesh` | Standard Lagrange-Legendre mesh: N, R, x, r, λ, α |
| `TransformedLagrangeMesh` | Coordinate-transformed mesh with g(x) mapping |
| `OpticalPotential` | Woods-Saxon parameters: V, r, a for volume/surface/spin-orbit + Coulomb |
| `ScatteringProblem` | Problem specification: potential, E_cm, l, j, N, R |
| `ScatteringResult` | Solution: S_matrix, phase_shift, cross_section, k, η |
| `WaveFunctionResult` | Wave function: r, ψ_sc, ψ_total, F_l, coefficients c |
| `NumerovResult` | Numerov integration result |

### Mesh Functions

| Function | Signature | Description |
|----------|-----------|-------------|
| `init_legendre_mesh` | `(N, R) → LagrangeMesh` | Create N-point mesh on (0,R) |
| `init_transformed_mesh` | `(N, R; transform, power) → TransformedLagrangeMesh` | Create transformed mesh |
| `lagrange_polynomial` | `(mesh, j, r) → Float64` | Evaluate L_j(r) |
| `fhat_at_boundary` | `(mesh, j) → Float64` | Basis function f̂_j(1) |
| `dfhat_dx_at_boundary` | `(mesh, j) → Float64` | Derivative df̂_j/dx at x=1 |
| `mesh_point_distribution` | `(mesh) → nothing` | Print mesh statistics |

### Matrix Functions

| Function | Signature | Description |
|----------|-----------|-------------|
| `baye_D_matrix` | `(mesh) → Matrix` | First derivative matrix (N×N) |
| `baye_T_matrix` | `(mesh) → Matrix` | Kinetic energy matrix (N×N) |
| `numerical_D_matrix` | `(mesh) → Matrix` | Numerical D for comparison |
| `numerical_T_matrix` | `(mesh) → Matrix` | Numerical T for comparison |

### Scattering Functions

| Function | Signature | Description |
|----------|-----------|-------------|
| `solve_scattering` | `(prob) → ScatteringResult` | Solve single partial wave |
| `solve_all_partial_waves` | `(pot, E; l_max, N, R) → Vector{ScatteringResult}` | Sum over l |
| `solve_wavefunction` | `(prob) → WaveFunctionResult` | Get wave function |
| `evaluate_wavefunction` | `(wf, r) → ComplexF64` | Interpolate ψ(r) |
| `elastic_differential_cross_section` | `(pot, E, θ) → Float64` | dσ/dΩ at angle θ |
| `solve_numerov` | `(prob; h) → NumerovResult` | Numerov integration |

### Coulomb Functions

| Function | Signature | Description |
|----------|-----------|-------------|
| `coulomb_k` | `(μ, E) → Float64` | Wave number k = √(2μE)/ℏc |
| `coulomb_eta` | `(Z1, Z2, μ, E) → Float64` | Sommerfeld parameter |
| `coulomb_F` | `(l, η, ρ) → Float64` | Regular Coulomb F_l |
| `coulomb_G` | `(l, η, ρ) → Float64` | Irregular Coulomb G_l |
| `coulomb_FG` | `(l, η, ρ) → (F, G, F', G')` | Both functions + derivatives |
| `coulomb_H_plus` | `(l, η, ρ) → ComplexF64` | Outgoing Hankel H⁺ = G + iF |
| `coulomb_gamma_s` | `(l, η, ρ, k) → ComplexF64` | Logarithmic derivative |
| `coulomb_sigma` | `(l, η) → Float64` | Coulomb phase shift σ_l |

### Potential Functions

| Function | Signature | Description |
|----------|-----------|-------------|
| `reduced_mass` | `(pot) → Float64` | μ = m₁m₂/(m₁+m₂) in MeV/c² |
| `woods_saxon` | `(r, R, a) → Float64` | Form factor 1/(1+exp((r-R)/a)) |
| `woods_saxon_derivative` | `(r, R, a) → Float64` | df/dr = -f(1-f)/a |
| `evaluate_potential` | `(pot, r, l, j) → ComplexF64` | Nuclear optical potential |
| `evaluate_coulomb` | `(pot, r) → Float64` | Coulomb potential (finite-size) |
| `evaluate_short_range` | `(pot, r, l, j) → ComplexF64` | V_nuc + V_coul^finite - V_coul^point |
| `evaluate_total_potential` | `(pot, r, l, j) → ComplexF64` | V_nuc + V_coul |

---

## Performance and Validation

### Convergence Guidelines

The Gauss-Legendre mesh has non-uniform point distribution: ~50% of points are in [R/2, R] where potentials are typically negligible.

**Recommended parameters for |S| accuracy ~10⁻³:**

| R (fm) | N (points) | Notes |
|--------|------------|-------|
| 8-10 | 30-40 | Optimal for short-range potentials |
| 15 | 50-60 | Standard calculations |
| 20 | 60-80 | Extended range |
| 40 | 100-120 | For COLOSS comparison |

**Optimization strategy:**
1. Identify where V(r) becomes negligible (typically ~6-8 fm for nuclei)
2. Set R ≈ (potential range) + 2-4 fm
3. Use coordinate-transformed mesh for additional efficiency

### Validation Results

SLAM.jl has been validated against the COLOSS Fortran code (D. Baye):

| Test Case | |S| Agreement | Re(S) Agreement |
|-----------|--------------|-----------------|
| p + ¹²C @ 30 MeV | ~10⁻⁵ | ~10⁻⁵ |
| n + ⁴⁰Ca @ 20 MeV | ~10⁻⁵ | ~10⁻⁵ |
| p + ⁴⁰Ca @ 30 MeV | ~10⁻⁴ | ~10⁻⁴ |

Numerical tests confirm:
- Wronskian identity: |F·G' - F'·G| = 1 to ~10⁻⁶
- |S| ≤ 1 for absorptive potentials (unitarity)
- Smooth convergence with increasing N

---

## Testing

### Run All Tests

```julia
using Pkg
Pkg.test("SLAM")
```

Or directly:

```julia
include("test/runtests.jl")
```

### Test Coverage

| Test Set | Description |
|----------|-------------|
| **Lagrange Mesh** | N, R, x, r, λ, α initialization; points in (0,1) and (0,R) |
| **Baye Matrices** | D and T matrix dimensions; structure verification |
| **Coulomb Functions** | Wronskian identity F·G' - F'·G = 1; H⁺ = G + iF |
| **Optical Potential** | Woods-Saxon limits; μ > 0; V imaginary part < 0 |
| **Scattering Solver** | Convergence flag; |S| ≤ 1.1; finite phase shifts |
| **Legendre Polynomials** | P_0 = 1, P_1 = x, P_2 = (3x²-1)/2 |

### Validation Scripts

```bash
# Compare with COLOSS
julia --project=. test/compare_coloss.jl

# Neutron scattering validation
julia --project=. test/compare_proton.jl

# Source term method comparison
julia --project=. test/compare_source_methods.jl
```

---

## Physical Constants

| Constant | Value | Description |
|----------|-------|-------------|
| ℏc | 197.3269718 MeV·fm | Reduced Planck constant × c |
| e² | 1.43997 MeV·fm | Coulomb constant (e²/4πε₀) |
| AMU | 931.49432 MeV | Atomic mass unit (c²) |

Values from NIST CODATA 2014.

---

## References

### Primary References

1. **D. Baye**, "The Lagrange-mesh method", *Physics Reports* 565 (2015) 1-107
   [doi:10.1016/j.physrep.2014.11.006](https://doi.org/10.1016/j.physrep.2014.11.006)
   *Comprehensive review of Lagrange mesh methods including exact D and T matrices*

2. **A.R. Barnett**, "COULFG: Coulomb and Bessel functions and their derivatives", *Computer Physics Communications* 21 (1981) 297-314
   [doi:10.1016/0010-4655(81)90030-2](https://doi.org/10.1016/0010-4655(81)90030-2)
   *COUL90 Fortran library for Coulomb functions*

3. **A.J. Koning & J.P. Delaroche**, "Local and global nucleon optical models from 1 keV to 200 MeV", *Nuclear Physics A* 713 (2003) 231-310
   [doi:10.1016/S0375-9474(02)01321-0](https://doi.org/10.1016/S0375-9474(02)01321-0)
   *KD02 global optical model parameterization*

### Additional References

4. **NIST Digital Library of Mathematical Functions**, Chapter 33: Coulomb Functions
   [https://dlmf.nist.gov/33](https://dlmf.nist.gov/33)

5. **I.J. Thompson & F.M. Nunes**, "Nuclear Reactions for Astrophysics", Cambridge University Press (2009)
   *Coordinate transformations and FRESCO methodology*

6. **G.R. Satchler**, "Direct Nuclear Reactions", Oxford University Press (1983)
   *Optical model potentials and scattering theory*

---

## Contributing

Contributions are welcome! Here's how to help:

### Reporting Issues

- Open an issue on GitHub with a clear description
- Include: Julia version, OS, minimal reproducing example
- Attach relevant output or error messages

### Pull Requests

1. Fork the repository
2. Create a feature branch: `git checkout -b feature/my-feature`
3. Write tests for new functionality
4. Ensure all tests pass: `julia --project=. -e 'using Pkg; Pkg.test()'`
5. Submit a pull request with clear description

### Development Guidelines

- Follow Julia style conventions
- Add docstrings for all public functions
- Include type annotations where helpful
- Update README for user-facing changes

### Ideas for Contribution

- Additional optical model parameterizations (Chapel-Hill, etc.)
- Coupled-channels extension
- Relativistic kinematics option
- Performance optimizations
- Additional test cases

---

## FAQ

### Q: Why is |S| > 1 for some calculations?

**A:** This typically indicates numerical issues:
- Increase N (number of mesh points)
- Check that R is large enough for the potential to vanish
- Verify potential parameters are physical (V_v > 0 for attractive)

### Q: How do I choose N and R?

**A:** Start with R = (potential range) + 2 fm and N = 50. Increase N until |S| changes by less than your required precision. For COLOSS comparison, use R = 40 fm, N = 100-120.

### Q: COUL90 library not found?

**A:** The pure Julia fallback is used automatically. For best accuracy, compile the Fortran library:
```bash
cd deps && gfortran -shared -fPIC -O3 -o libcoul90.dylib coul90.f90
```

### Q: How accurate is the pure Julia Coulomb implementation?

**A:** The continued fraction method achieves ~10⁻¹⁵ relative accuracy for most cases. It may have issues for very large η or very small ρ.

### Q: Can I use SLAM.jl for deuteron or alpha scattering?

**A:** Yes! Set `A_proj` and `Z_proj` appropriately. For spin > 1/2 projectiles, the spin-orbit calculation needs modification (currently assumes s = 1/2).

### Q: How do I extract the optical potential at a given radius?

**A:**
```julia
V_total = evaluate_total_potential(pot, r, l, j)  # V_nuc + V_coul
V_short = evaluate_short_range(pot, r, l, j)      # For inhomogeneous equation
```

---

## License

MIT License

Copyright (c) 2025 Jin Lei

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

---

## Author

**Jin Lei**
Email: jinlei.phys@gmail.com

---

*Last updated: December 2025*
