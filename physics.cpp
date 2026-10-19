#include <iostream>
#include <cmath>
#include <algorithm>
#include "constants.h"
#include "particle.h"
#include "physics.h"


// =========================
//   Physics Computation
// =========================

// Gravitational force
Vec2 gravityForce(const Particle& a, const Particle& b) {
    Vec2 r = b.pos - a.pos;
    double d2 = r.lengthSquared() + eps*eps;                   // distance squared => d^2 => d2
    Vec2 gForce = ((G*a.mass*b.mass)/(d2 * std::sqrt(d2)))*r;  // Scalar: (G*a.mass*b.mass)/(d2 * std::sqrt(d2)) hidden inside is r_hat, makes gForce a vector
    return gForce;
}

// Electrical force
Vec2 electricForce(const Particle& a, const Particle& b) {
    Vec2 r = b.pos - a.pos;
    double d2 = r.lengthSquared() + eps*eps;
    Vec2 eForce = (-(k*a.charge*b.charge)/(d2 * std::sqrt(d2)))*r;  // minus(-) in the numerator because same charges repel, opposites attract
    return eForce;
}

// Pairing forces gravitiyForce + electricForce
Vec2 pairForce(const Particle& a, const Particle& b) {
    Vec2 r = b.pos - a.pos;
    double d2 = r.lengthSquared() + eps*eps;
    double totald = d2 * std::sqrt(d2);
    double scale = (G*a.mass*b.mass - k*a.charge*b.charge) / totald;   // gravity (+) and Coulomb (-) folded into one number
    return scale * r;
}

// Total force acted upon each particle relative to all particles
void computeForces(std::vector<Particle>& particles) {
    const size_t N = particles.size();
    
    for(size_t i = 0; i < N; i++){
        particles[i].force = {0,0}; // Reset all forces to 0 before calculation
    }

    for (size_t i = 0; i < N; i++) {
        for(size_t j = i + 1; j < N; j++) {
            Vec2 F = pairForce(particles[i], particles[j]);
            particles[i].force += F;
            particles[j].force -= F;
        }
    }
}

void integrate(std::vector<Particle>& particles, double dt) {
    const size_t N = particles.size();

    for (size_t i = 0; i < N; i++) {
        //Calculating acceleration, velocity, and displacement
        Vec2 a = particles[i].force / particles[i].mass;
        Vec2 v = particles[i].vel;
        Vec2 x = particles[i].pos;
        v += a * dt;
        x += v * dt;
        particles[i].vel = v;
        particles[i].pos = x;
    }
}

// Collision physics between 2 particles
void handleCollisions(std::vector<Particle>& particles) {
    const size_t N = particles.size();

    for (size_t i = 0; i < N; i++) {
        for (size_t j = i + 1; j < N; j++) {
            Vec2 r = particles[j].pos - particles[i].pos;   // Distance between 2 particles in vector format
            double d2 = r.lengthSquared();                  // Magnitude^2
            
            if (d2 < (particles[j].radius + particles[i].radius) * (particles[j].radius + particles[i].radius)) {   // Collision check
                double distance = std::sqrt(d2);    // Magnitude
                Vec2 n;                             // Normal vector relative to pos i and j
                if (distance < 1e-12) {
                    n = {1, 0};                     // any direction works; they just need to be pushed apart
                    distance = 0.0;
                } else {
                    n = r / distance;
                }   
                double overlap = (particles[j].radius+particles[i].radius) - distance;
                double w_i = (1/particles[i].mass)/(1/particles[i].mass + 1/particles[j].mass);
                double w_j = (1/particles[j].mass)/(1/particles[i].mass + 1/particles[j].mass);
                // if overlap is more than 0, change the position of particles based on the overlap
                particles[i].pos -= n * overlap * w_i;
                particles[j].pos += n * overlap * w_j;
                // Collision physics
                double vn = dot((particles[j].vel - particles[i].vel), n);  // scalar of (vj-vi) of n
                if (vn < 0) {                                               // are they moving towards each other?
                    double J = -(1 + std::min(particles[i].restitution, particles[j].restitution))*vn/(1/particles[i].mass + 1/particles[j].mass);
                    particles[i].vel -= (J/particles[i].mass) * n;
                    particles[j].vel += (J/particles[j].mass) * n;
                }
            }
        }
    }
}

// Collision physics between the particle and any wall
void handleWalls(std::vector<Particle>& particles) {
    const size_t N = particles.size();

    for (size_t i = 0; i < N; i++) {
        Particle& p = particles[i];
        
        // Wall collision check
        if (p.pos.x > boxHalf - p.radius && p.vel.x > 0) {  // right wall
            p.pos.x = boxHalf - p.radius;                   // clamp: put it back on the wall
            p.vel.x = -p.restitution * p.vel.x;             // bounce
        }
        if (p.pos.x < p.radius - boxHalf && p.vel.x < 0) {  // left wall
            p.pos.x = p.radius - boxHalf;             
            p.vel.x = -p.restitution * p.vel.x;                         
        }
        if (p.pos.y > boxHalf - p.radius && p.vel.y > 0) {  // up wall 
            p.pos.y = boxHalf - p.radius;                   
            p.vel.y = -p.restitution * p.vel.y;                         
        }
        if (p.pos.y < p.radius - boxHalf && p.vel.y < 0) {  // down wall
            p.pos.y = p.radius - boxHalf;                   
            p.vel.y = -p.restitution * p.vel.y;                         
        }
    }
}

// Total energy calculation
double totalEnergy(const std::vector<Particle>& particles) {
    const size_t N = particles.size();
    double KE = 0.0;    // Kinetic energy
    double PE = 0.0;    // Potential energy
    
    for (size_t i = 0; i < N; i++) {
        KE += 0.5 * particles[i].mass * (particles[i].vel.lengthSquared());
        for (size_t j = i + 1; j < N; j++) {
            Vec2 r = particles[j].pos - particles[i].pos;
            double d2 = r.lengthSquared() + eps*eps;
            PE += (-G*particles[i].mass*particles[j].mass + k*particles[i].charge*particles[j].charge)/std::sqrt(d2);
        }
    }
    return KE + PE;
}