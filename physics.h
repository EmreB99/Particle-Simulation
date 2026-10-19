#pragma once
#include <vector>
#include "particle.h"

Vec2   pairForce(const Particle& a, const Particle& b);
void   computeForces(std::vector<Particle>& particles);
void   integrate(std::vector<Particle>& particles, double dt);
void   handleCollisions(std::vector<Particle>& particles);
void   handleWalls(std::vector<Particle>& particles);
double totalEnergy(const std::vector<Particle>& particles);