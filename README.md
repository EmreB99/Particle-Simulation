# Particle Simulation (C++)

A 2D particle simulation I'm building from scratch to learn both C++ and the physics behind it:
gravity, electric (Coulomb) forces, collisions, and eventually a black hole bending light.
2D first, then 3D.

## What works so far
- `Vec2`: 2D vector type with operator overloading (`+`, `-`, `*`, `/`, `+=`, `-=`, dot product, length, normalize)
- `Particle`: position, velocity, force, mass, charge, radius
- Newtonian gravity and Coulomb force between every pair of particles (softened to avoid infinities at r → 0)
- Pairwise force loop using Newton's 3rd law (each pair computed once)
- Semi-implicit (symplectic) Euler integrator

## Checks so far
- Total momentum stays at 0 and the center of mass doesn't move
- Two equal charges with q² = m² (G = k = 1): gravity and repulsion cancel exactly, so nothing moves
- Two-body free-fall time matches the analytic result t = (π/2)·√(r₀³ / (2GM))
- Smaller time step → the oscillation returns to its starting amplitude, as energy conservation requires

## Build and run
```bash
g++ -std=c++17 -Wall -Wextra main.cpp -o sim
./sim
```

## To-do list
- Energy diagnostic (kinetic + potential) to measure numerical error
- CSV output + Python animation
- Velocity Verlet integrator, adaptive time step
- Move to 3D
- Black hole + photons
