# Continuum Approximation of a String via Infinite Coupled Pendulums

This repository contains the theoretical derivations, analytical solutions, and numerical simulations exploring the transition from a discrete system of coupled pendulums to a continuous string medium.

## Project Overview

The project models the mechanics of continuous media by taking the infinitesimal limit of a discrete Lagrangian system. Starting with a chain of coupled pendulums, the physics are scaled to analyze wave propagation, normal modes, and chaotic dynamics under specific boundary conditions.

### Key Physics & Mathematical Concepts Covered:
* **Discrete Lagrangian Mechanics:** Formulating the Kinetic ($K$) and Potential ($U$) energies of $N$ coupled pendulums using the Principle of Least Action.
* **The Continuum Limit:** Transitioning from discrete algebraic equations of motion to partial differential equations (PDEs) representing continuous media.
* **Analytical Solutions:** Solving the resulting wave equation via separation of variables and Fourier series analysis.
* **Nonlinear Dynamics & Chaos:** Implementation of numerical tools to model nonlinear perturbations, including the study of phase space trajectories and caothic behaviour through **Poincaré sections**.

## Repository Structure

* `La_cuerda_de_los_infinitos_pendulos.pdf`: The complete formal academic paper containing full algebraic derivations, analytical frameworks, and structured nomenclature.
* `src/`: Directory containing the numerical solvers and scripts (MATLAB) used to simulate the dynamics, solve the ODE systems (`ode45` framework), and plot the Poincaré sections.

## Documentation & Code

The full mathematical paper with all derivations is available directly in this repository:

👉 **[Read the Full Paper (PDF)](./La%20cuerda%20de%20los%20infinitos%20pendulos.pdf)**

## Author

* **Ernest J. Redolat Guillén** – Engineering Physics Undergraduate Student, Universitat Politècnica de Catalunya (UPC).
* *Co-authored with Cesc Barrera (April 2026).*
