# 2D Ising Model — Monte Carlo Simulation & Analysis

This project implements a full Monte Carlo study of the **2D ferromagnetic Ising model** on a square lattice with periodic boundary conditions, using Julia. It covers simulation, statistical error analysis, and the extraction of critical exponents via finite-size scaling.

The coupling constant is set to **J = 1** throughout, so temperatures are given in units of J/k_B.

---

## Repository Structure

```
├── main_ising_2d.jl      # Monte Carlo simulation engine (data generation)
├── ising_analysis.jl     # Post-processing, plotting, and critical exponent fitting
└── README.md
```

---

## Physics Background

The 2D Ising model describes spins s_i = ±1 on a square lattice interacting via the Hamiltonian:

```
H = -J Σ s_i s_j       (sum over nearest-neighbour pairs)
```

The model has an exact analytical solution (Onsager, 1944) with a second-order phase transition at the critical temperature **T_c ≈ 2.269 J/k_B**. Below T_c the system is ferromagnetically ordered; above it the spins are disordered. Near T_c, thermodynamic observables diverge with universal critical exponents:

```
m     ~ |t|^β        (order parameter,       β = 1/8)
c_v   ~ |t|^(-α)     (heat capacity,         α = 0  — logarithmic divergence)
χ     ~ |t|^(-γ)     (susceptibility,        γ = 7/4)
ξ     ~ |t|^(-ν)     (correlation length,    ν = 1)
```

where `t = (T - T_c) / T_c` is the reduced temperature.

---

## `main_ising_2d.jl` — Simulation Engine

### Overview

Runs Monte Carlo simulations for lattice sizes **L = 4, 8, 16, 32, 64, 128** (or any subset defined in `L_list`) over a temperature sweep from T = 3.00 down to T = 2.00 in steps of 0.01. For each (L, T) pair it takes **1000 measurements** and saves all results to CSV files.

### Algorithm Strategy

The temperature range is split into four regions, each with a different update algorithm and measurement spacing. This is necessary because Metropolis suffers from **critical slowing down** near T_c, making it inefficient there:

| Temperature range | Algorithm   | MCS between measurements |
|---|---|---|
| T ∈ [2.80, 3.00]  | Metropolis  | 40  |
| T ∈ [2.61, 2.79]  | Metropolis  | 40  |
| T ∈ [2.30, 2.60]  | Wolff       | 200 |
| T ∈ [2.00, 2.29]  | Wolff       | 200 |

A thermalization phase of **2000 MCS** is run before measurements begin at each temperature. The lattice is initialized in a random spin configuration at the start.

### Key Functions

**`neighbours_creator(L)`**
Builds four arrays encoding the right, top, left, and bottom neighbour index of every spin under periodic boundary conditions. These are precomputed once per L and stored as a (4 × N) integer matrix `neigh` for fast lookup during updates.

**`metropolis(s_array, neigh, T)`**
Standard single-spin-flip Metropolis update. Performs N attempted flips per call (one full sweep). Acceptance probabilities for the two positive energy-change cases (Δbi = 2 and 4) are precomputed and cached to avoid repeated exponential evaluations.

**`wolff(i, s_array, neigh, pa)`**
Iterative (stack-based) implementation of the Wolff cluster algorithm. Starting from a randomly chosen seed spin, it grows a cluster of aligned spins and flips the entire cluster at once. The bond-addition probability is `p_a = 1 - exp(-2J/T)`. This eliminates critical slowing down near T_c. A `visited` boolean array prevents spins from being added to the cluster twice.

**`hamiltonian_and_M(s_array, neigh, N)`**
Computes the total energy H and the absolute magnetisation |M| of the current configuration in a single pass over all spins. The energy includes the standard factor of 1/2 to avoid double-counting bonds.

**`Correlation_time(Observable)`**
Estimates the integrated autocorrelation time τ from the lag-1 autocorrelation ρ(1) using the approximation:
```
τ ≈ ρ(1) / (1 − ρ(1))
```
This is used to correct statistical error bars, since consecutive measurements are correlated.

**`r1k(L)`**
Computes the Euclidean distance from site (1,1) to every other site in the lattice, respecting periodic boundary conditions. Used to bin the spin-spin correlation function g(r).

**`simulation(L, N, s_array, ...)`**
The main simulation loop. For each temperature it:
1. Thermalizes the system.
2. Takes 1000 measurements of H and M, storing also H², M², and M⁴.
3. Computes the spin-spin correlation function g(r) = ⟨s₁ sⱼ⟩ − ⟨m⟩².
4. Fits g(r) to an exponential decay `exp(-r/ξ)` to extract the correlation length ξ.
5. Computes the autocorrelation time τ for H, M, H², M², and M⁴.

### Output CSV Files

For each lattice size L the simulation writes three CSV files:

| File | Contents |
|---|---|
| `simulation_data_L_ising_2d.csv` | Per-temperature means and std of H, M, H², M², M⁴, plus T |
| `correlation_time_L_ising_2d.csv` | Autocorrelation times τ_H, τ_M, τ_H², τ_M², τ_M⁴ per temperature |
| `correlation_function_L_ising_2d.csv` | Fitted correlation length ξ and its error σ_ξ per temperature |

---

## `ising_analysis.jl` — Post-Processing & Plotting

### Overview

Reads the CSV files produced by the simulation and computes physics quantities, fits critical exponents, and generates all plots. The code is structured around a set of helper functions that eliminate repetition across lattice sizes.

### Statistical Error Treatment

Raw standard deviations from the simulation are **not** valid error bars because consecutive measurements are correlated. The corrected error on a mean is:

```
δ⟨X⟩ = σ_X × sqrt( (2τ + 1) / N_meas )
```

where τ is the autocorrelation time loaded from `correlation_time_L_*.csv` and N_meas = 1000. This is applied to all observables via `error_from_tau`.

### Physics Quantities Computed

From the stored means and variances the following per-spin quantities are derived:

- **Energy per spin** `u = ⟨H⟩ / N`
- **Magnetisation per spin** `m = ⟨|M|⟩ / N`
- **Heat capacity per spin** `c_v = (⟨H²⟩ - ⟨H⟩²) / (N T²)`
- **Magnetic susceptibility per spin** `χ = (⟨M²⟩ - ⟨M⟩²) / (N T)`
- **Binder cumulant** `U = 1 - ⟨M⁴⟩ / (3 ⟨M²⟩²)`

### Analysis Functions

**`plot_τ()`**
Plots the correlation length ξ as a function of temperature for L = 64 and L = 128, with error bars. Saves individual PNG files per L.

**`plot_measures()`**
Plots u, m, c_v, and χ vs T for all lattice sizes simultaneously, with autocorrelation-corrected error bars. Also plots the **Binder cumulant** U vs T zoomed near T_c (T ∈ [2.255, 2.275]). The crossing point of U curves for different L gives a direct estimate of T_c that is independent of critical exponents.

**`crit_exponents()`**
Extracts critical exponents by fitting power laws in log-log space near T_c, using the L = 256 dataset restricted to T ∈ [2.00, 2.267):

- **β** from `log(m)` vs `log|t|`
- **−α** from `log(c_v)` vs `log|t|`
- **−γ** from `log(χ)` vs `log|t|`

Additionally fits **γ/ν** from the finite-size scaling of the susceptibility peak χ_max ~ L^(γ/ν), using all L sizes. From this slope and the known value γ = 7/4, the exponent ν is estimated. All fits use ordinary least squares (`GLM.lm`) on the linearised log-log data, reporting slope, intercept, standard errors, and R².

**`finite_size()`**
Produces a **finite-size scaling collapse** plot. If the scaling hypothesis holds, all χ(T, L) curves should collapse onto a single universal function when plotted as:
```
χ L^(−γ/ν)   vs   t L^(1/ν)
```
using the 2D Ising exact exponents γ/ν = 7/4 and ν = 1. A good collapse confirms the exponent values and the location of T_c.

### Running the Analysis

Uncomment the desired function at the bottom of `ising_analysis.jl`:

```julia
# plot_τ()
# plot_measures()
# crit_exponents()
finite_size()
```

All plots are saved as PNG files in the same directory as the script.

---

## Dependencies

Both scripts require the following Julia packages:

```julia
using Random, Plots, LaTeXStrings, GLM, DataFrames
using CSV, Statistics, StatsBase, LsqFit, LinearAlgebra
```

Install them from the Julia REPL with:

```julia
import Pkg
Pkg.add(["Plots", "LaTeXStrings", "GLM", "DataFrames",
         "CSV", "StatsBase", "LsqFit"])
```

---

## Configuration

Before running the analysis, set the path constants at the top of `ising_analysis.jl` to point to the directory where the simulation CSVs were saved:

```julia
const BASE_PATH = "path/to/your/csv/directory"
```

Lattice sizes and temperature ranges are controlled at the top of `main_ising_2d.jl` via `L_list` and the four `T_array_*` ranges.
