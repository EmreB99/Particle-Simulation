#pragma once
#include <vector>
#include <random>
#include "vec2.h"

//=========================
//   Particle Creation
// ========================

struct Particle {
    double mass = 1.0, charge = 0.0, radius = 1.0, restitution = 1.0;
    Vec2 pos, vel, force;  
};
Particle randomParticle(std::mt19937& rng);
std::vector<Particle> createParticles(size_t N, std::mt19937& rng);