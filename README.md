# BimBamBum - Bubble Dynamics Solver

## 🧠 High-Level Overview

**BimBamBum** is a high-performance C++ numerical solver used to simulate axisymmetric bubble dynamics near fluid-fluid interfaces. It relies on the **Boundary Integral Method (BIM)** to solve Laplace's equation for velocity potentials, allowing for the simulation of bubble expansion, contraction, and collapse (including jetting phenomena).

**Key Physics Capabilities:**

* **Bubble Types:** Rayleigh bubbles, Rayleigh-Plesset bubbles.
* **Boundary Types:** Flat interfaces, Spherical droplets.
* **Math & Integration:** Uses 1st/2nd order Runge-Kutta for time-stepping, Longuet-Higgins filtering to smooth interfaces, and custom cubic splines for spatial derivatives.

**Dependencies:** C++14, Armadillo (Linear Algebra), GSL (Numerical Integration/Root-finding), Boost (JSON parsing), OpenMP (Multi-threading).

---

## 📂 Codebase Map (Where to find what)

NotebookLM, use this map to navigate the components of the simulation:

### 1. The Engine (Execution & Setup)

* **`solver/main.cpp`**: The entry point. Parses the JSON config, sets up OpenMP threads, initializes the bubble and boundary, and runs the main time-integration `while` loop.
* **`include/core/Inputs.hpp` & `src/stream/ConfigFileParser.cpp`**: Reads the `.json` files to set physical parameters (like node counts `Nb`, `Ns`, pressure constants, and solver types).
* **`Example_JSON/*.json`**: Configuration files (e.g., `A30_14_e0.json`) that dictate the specific simulation parameters.

### 2. The Physics Data Structures

* **`BubbleData` (`.hpp` / `.cpp`)**: Base class for the bubble. Manages nodes (`r_nodes`, `z_nodes`), velocities, volume computation, node remeshing, and surface filtering.
  * *Derived Classes:* `Case_RayleighBubble`, `Case_RayleighPlessetBubble` (Handles specific initial geometries and velocities).
* **`BoundaryData` (`.hpp` / `.cpp`)**: Base class for the fluid-fluid interface. Manages surface nodes, curvatures, and velocities.
  * *Derived Classes:* `Case_FlatBoundary`, `Case_SphericalBoundary`.

### 3. The Math & Solvers

* **`include/core/BIM_solver.hpp`**: The mathematical core. Constructs the Green's function matrices and solves the dense linear system using Armadillo to find the normal and tangential velocities of the boundaries.
* **`include/core/integrands_BIM.hpp` & `integrands_volume.hpp`**: Contains the complex elliptic integrals (G1, G2, H1, H2) needed for the Boundary Integral Method evaluated via GSL.
* **`include/core/time_stepper.hpp`**: Advances the simulation in time. Implements the temporal loop using RK1 (Euler) or RK2 (Heun's method), updating potentials (`phi`) and coordinates using the Bernoulli equation.
* **`cubic_spline` (`.hpp` / `.cpp`)**: Custom clamped-end cubic spline interpolator used to calculate smooth spatial derivatives (tangential vectors) along the discretized bubble and boundary curves.

---

## ⚙️ Execution Flow (How the code runs)

When `main.cpp` executes, it follows this exact sequence:

1. **Read Input:** Loads the `Inputs` object from the provided JSON file.
2. **Initialize:** Creates the `BIM_solver` object. Instantiates a specific `BubbleData` and `BoundaryData` child class based on the JSON strings (e.g., `Rayleigh_Plesset_Bubble` and `spherical`).
3. **Time Stepping Loop:** A `while(compute)` loop begins.
   * **Write Data:** Dumps current node positions/velocities to a `.txt` file.
   * **Solve BIM & Step Time:** `time_integration()` is called. `BIM_solver` calculates normal/tangential velocities. Runge-Kutta advances the node positions (`r`, `z`) and potentials (`phi`) by `dt`.
   * **Maintain Mesh:** Calls `remesh_bubble()` and `remesh_boundary()` to redistribute nodes evenly and prevent numerical instability.
   * **Check Stop Condition:** If `bubble->intersect() == 1` (the bubble wall hits itself or the boundary, forming a jet or pinching off), the simulation terminates.

---

## 🔑 Key Terminology for AI Analysis

* **`Nb` and `Ns`**: Number of nodes discretizing the Bubble (`Nb`) and the Boundary surface (`Ns`).
* **Green's Functions**: The fundamental solution to Laplace's equation, evaluated over the axisymmetric boundary elements to compute surface velocities.
* **Remeshing**: A geometric progression algorithm used to maintain higher node density at the poles (or areas of high curvature) while keeping the mesh stable as the bubble deforms.
* **Filtering**: Implementation of Longuet-Higgins & Cokelet smoothing (applied every few steps) to remove high-frequency numerical noise (sawtooth instabilities) on the node surfaces.
