# Fluid Dynamics Simulation

This repository contains a Python script for simulating fluid dynamics, specifically focusing on heat and mass transfer in a specific physical domain. The simulation utilizes a finite difference method to solve a set of coupled partial differential equations.

## Features

* **Tridiagonal System Solver:** Implements a robust solver for tridiagonal matrix systems, a common requirement in finite difference methods.
* **Configurable Parameters:** Easily adjust grid dimensions, physical constants (e.g., Prandtl, Schmidt numbers), and simulation parameters (e.g., time step, maximum simulation time) through constants at the top of the script.
* **Time Integration:** Advances the simulation through time steps, updating velocity, temperature, and concentration fields.
* **Convergence Monitoring:** Tracks the maximum differences between successive iterations of the fields to determine simulation convergence.
* **Detailed Output:** Provides informative print statements regarding simulation progress and convergence status.

## Getting Started

### Prerequisites

* Python 3.x
* NumPy library (`pip install numpy`)

### Running the Simulation

1.  **Save the code:** Save the provided Python code as a `.py` file (e.g., `fluid_simulation.py` or `simulation.py`).
2.  **Run from terminal:** Open your terminal or command prompt, navigate to the directory where you saved the file, and run:

    ```bash
    python fluid_simulation.py
    ```

The simulation will start, printing progress updates to the console.

## Code Overview

### Global Parameters

The simulation's behavior is heavily influenced by a set of global constants defined at the beginning of the script:

* `EPSILON`: Convergence criterion for the simulation.
* `Q_PARAM`, `R_PARAM`, `RT_PARAM`: Physical parameters specific to the problem.
* `M_GRID`, `N_GRID`: Number of grid points in the x and y directions, respectively.
* `X_MAX`, `Y_MAX`: Maximum extent of the computational domain in x and y.
* `DTAU`: Time step size for the simulation.
* `TAU_MAX`: Maximum dimensionless time for the simulation.
* `GRASHOF_NUM`, `PRANDTL_NUM`, `SCHMIDT_NUM`: Dimensionless numbers characterizing the fluid flow, heat transfer, and mass transfer.
* `G_MOD`: Modified gravitational parameter.

### Functions

* `solve_tridiagonal_system(first_idx, last_idx, a, b, c, d, x_out)`: Solves a tridiagonal linear system using a specialized algorithm (likely a variant of the Thomas algorithm or TDMA).
* `initialize_simulation_arrays(m_plus_1, n_plus_1)`: Initializes all the NumPy arrays required for storing velocity, temperature, and concentration fields, as well as coefficients for the tridiagonal solver.
* `run_fluid_simulation()`: The main function that orchestrates the simulation. It initializes parameters, runs the time integration loop, and monitors convergence.

### Simulation Logic

The `run_fluid_simulation` function performs the following steps:

1.  **Initialization:** Sets up the computational grid and initializes all field variables (U, V, T, W, C) to their initial conditions.
2.  **Time Loop:** Iterates until either the `TAU_MAX` is reached or the solution converges.
3.  **Field Updates:** Within each time step, it updates the velocity, temperature, and concentration fields based on the governing equations. (The specific equations are not detailed in the provided snippet but are implicitly handled through the coefficient calculations that would precede the `solve_tridiagonal_system` calls, which are not present in this excerpt).
4.  **Convergence Check:** After each update, it calculates the maximum absolute difference between the current and previous time step's field values. If these differences fall below `EPSILON`, the simulation is considered converged.

## Customization

You can modify the global parameters at the beginning of the `simulation.py` file to experiment with different simulation scenarios:

* **Grid Resolution:** Adjust `M_GRID` and `N_GRID` for finer or coarser grids.
* **Physical Parameters:** Change `GRASHOF_NUM`, `PRANDTL_NUM`, `SCHMIDT_NUM` to simulate different fluids or thermal/mass transfer conditions.
* **Time Stepping:** Modify `DTAU` and `TAU_MAX` to control the simulation's speed and duration.

## Contributing

This is a basic simulation setup. Further development could include:

* Implementation of the actual finite difference equations for U, V, T, W, and C within the time loop.
* Adding visualization capabilities (e.g., using Matplotlib) to plot the field variables.
* More sophisticated boundary conditions.
* Error handling and input validation.
