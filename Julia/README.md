# nTOFHist.jl

`nTOFHist.jl` is a Julia module designed for simulating and visualizing Time-of-Flight (ToF) spectra, specifically for experiments involving neutron detection (e.g., beta-delayed neutron emission). It provides tools to model detector response, simulate neutron event counts using Monte Carlo methods, and generate publication-quality plots.

## Core Functionalities

* **Neutron Efficiency Calculation ($\epsilon_n$):** Uses a polynomial calibration to determine detection efficiency as a function of neutron energy ($E_n$).
* **ToF Computation:** Calculates the expected arrival time for neutrons at a detector based on path length, $Q_\beta$ value, separation energy ($S_n$), and excitation energies.
* **Transition Handling:** Provides utilities to select specific nuclear transitions and calculate beta-decay intensities (branching ratios) from Gamow-Teller strengths (BGT).
* **Monte Carlo Simulation (`MCHist`):** Generates synthetic ToF histograms using rejection sampling. It simulates three distinct components:
    1.  **Peak events:** Neutron-induced signals.
    2.  **Prompt events:** The "gamma flash" or prompt background.
    3.  **Background events:** Random/steady-state background levels.
* **Visualization:**
    * `plotSim`: Provides a view of the integrated detector response.
    * `plotMCHist`: Produces detailed plots including the summed response, individual transition contributions, and neutron energy annotations. It supports grouping transitions by spin.

## Dependencies

### Julia Packages
* `Plots`
* `Measures`

### Local Modules/Files
The following files must be present in the same directory or reachable via `include`:
* `MonteCarlo.jl`
* `VandleResponses.jl`

### External Dependencies
* `BetaDecayUtils.jl` (must be available in the Julia environment)

## Key Functions

| Function | Description |
| :--- | :--- |
| `εᵥ(Eₙ, εᵦ=0.8)` | Calculates neutron detection efficiency for a given energy $E_n$. |
| `ToF(path, Qᵦ, Sₙ, Eₓ)` | Computes arrival time at distance `path` for neutrons with excitation energy $E_x$. |
| `MCHist(...)` | Performs Monte Carlo simulation to generate a ToF histogram. |
| `plotMCHist(...)` | Generates a detailed histogram plot with various components and annotations. |

## Usage Example (Conceptual)

```julia
using nTOFHist

# Define parameters
path = 50.0    # cm
Qβ = 7.0       # MeV
Sn = 2.0       # MeV
Ex = [2.5, 3.0, 4.0] # Excitation energies (MeV)
BGT = [0.1, 0.05, 0.01] # Gamow-Teller strengths

# Generate simulation
sample = MCHist(0.1, path, 47, Qβ, Sn, Ex, BGT, 0.01, 1000, 0.0, 100.0)

# Plot results
plotMCHist(0.1, 10, path, 47, Qβ, Sn, Ex, BGT, 0.01, 1000, 0.0, 100.0)
```
