# Lattice Boltzmann Method for Thermal Flow Around a Cylinder

![1](https://github.com/user-attachments/assets/7766006a-be2f-4993-9287-f07c32ee1abb)


This project simulates two-dimensional thermal fluid flow around a cylindrical obstacle using the Lattice Boltzmann Method (LBM). The simulation models how an inlet flow interacts with a circular cylinder, allowing visualization of velocity fields, flow separation, and wake formation. Adjustable parameters let users explore the effects of changing inlet velocity and fluid viscosity.

## Features

- **2D Lattice Boltzmann simulation** of incompressible flow past a cylinder.
- Circular obstacle placed in a channel for realistic flow scenarios.
- **Interactive real-time visualization** of:
  - Velocity magnitude mapped onto the flow field
  - Obstacle shape and position within the grid
- **Adjustable parameters with sliders:**
  - Inlet velocity
  - Viscosity (via omega parameter)
- **Reset and control** functions for instantaneous experimentation.
- Built-in **Numba acceleration** if available for optimized computation.

## How it works

- The domain is a 2D rectangular channel containing a circular obstacle (cylinder).
- **Lattice Boltzmann Method (D2Q9):**
  - Simulates the Boltzmann equation on a square lattice for fluid densities and velocities.
  - Implements **streaming, collision, and bounce-back** steps to advance the solution.
- **Boundary Conditions:**
  - Inlet: Zou-He velocity inlet for controlled flow.
  - Outlet: Open boundary with zero gradient for smooth outflow.
  - Cylinder: No-slip boundary via bounce-back condition.
- **Visualization:**
  - The simulation displays velocity magnitude and highlights the cylinder's position.
  - A colorbar shows the range of velocities present in the channel.


## Usage

1. **Clone the repository:**
   ```bash
   git clone https://github.com/AakaashfromIndia/Lattice-Botlzmann-Method-for-thermal-flow-around-a-Cylinder.git
   ```

2. **Install dependencies:**
   ```bash
   pip install numpy matplotlib
   ```
   (Optional for performance: `pip install numba`)

3. **Run the simulation:**
   ```bash
   python Fluid-Flow-using-Lattice-Boltzmann-Method.py
   ```

4. **Interact with controls:**
   - Use the **Inlet Velocity** and **Omega (viscosity)** sliders to change fluid parameters.
   - Visualize the flow and velocity field evolving around the cylinder in real time.

## Dependencies

- Python 3.x
- NumPy
- Matplotlib
- (Optional) Numba (for performance improvement)
