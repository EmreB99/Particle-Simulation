// =========================
//  Particle Simulation
//  Emre Bilgen
// =========================

#include <iostream>
#include <cmath>
#include <random>
#include <chrono>
#include "particle.h"
#include "physics.h"

int main() {

    // Creating particles
    std::mt19937 rng(42);                                       // the generator (seed 42)
    std::vector<Particle> particles = createParticles(50, rng);

    // Constants of experiment
    const size_t N = particles.size();
    const double dt = 0.001;
    const int numSteps = 60000; //Number of calculations

    const double E0 = totalEnergy(particles);                   // energy at t = 0, to measure drift against
    std::cout << "Initial energy: " << E0 << "\n";

    auto startTime = std::chrono::steady_clock::now();          // start the timer

    // Main force-position loop
    for (int step = 0; step < numSteps; step++) {
        computeForces(particles);
        integrate(particles, dt);
        handleCollisions(particles);
        handleWalls(particles);

        if (step % 500 == 0){
            std::cout << "---" << (step + 1) * dt << "s ---\n";
            for (size_t i = 0; i < N; i++) {
                std::cout << "Particle " << i << " Position: " << particles[i].pos << " Force: " << particles[i].force << "\n";
            }
            double E = totalEnergy(particles);                  // once per block: it's the energy of the whole system
            std::cout << "Total energy: " << E << "   drift: " << 100.0 * (E - E0) / std::fabs(E0) << " %\n";
        }
    }

    auto endTime = std::chrono::steady_clock::now();            // stop the timer
    double seconds = std::chrono::duration<double>(endTime - startTime).count();
    std::cout << "Simulated " << numSteps << " steps of " << N << " particles in " << seconds << " s\n";
    return 0;
}