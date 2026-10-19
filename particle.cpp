#include <iostream>
#include "constants.h"
#include "particle.h"

//=========================
//    Particle Creation
// ========================

// Create a random particle
Particle randomParticle(std::mt19937& rng) {
    Particle p;
    p.radius = 0.3; // Create random radius for each particle later

    std::uniform_real_distribution<double> posDist(-boxHalf + p.radius, boxHalf - p.radius);
    p.pos.x = posDist(rng);
    p.pos.y = posDist(rng);

    std::uniform_real_distribution<double> mass(0.5, 3.0);
    p.mass = mass(rng);
    std::uniform_real_distribution<double> charge(-1.0, +1.0);
    p.charge = charge(rng);
    return p;
}

// Create any amount of random particles
std::vector<Particle> createParticles(size_t N, std::mt19937& rng) {
    std::vector<Particle> result;
    size_t attempts = 0;                // Attempt count
    const size_t maxAttempts = 100 * N; // Max number of attempts, if attempt number is exceeded than break the while loop


    while (result.size() < N && attempts < maxAttempts) {
        attempts++;
        Particle candidate = randomParticle(rng);   // Create particle
        bool overlaps = false;

        for (size_t j = 0; j < result.size(); j++) {
            Vec2 r = result[j].pos - candidate.pos;
            double minDist = candidate.radius +result[j].radius;    // Disadvantage between new particle and any particle

            if (r.lengthSquared() < minDist * minDist) {            // Overlap condition
                overlaps = true;
                break;
            }
        }

        if (!overlaps) {    // If it does not overlap push the particle
            result.push_back(candidate);
        }
    }

    if (result.size() < N) {    // If box is too crowded, only place x amount of particles
        std::cout << "Warning: only placed " << result.size() << " of " << N 
                  << " particles (box too crowded)\n";
    }
    return result;
}