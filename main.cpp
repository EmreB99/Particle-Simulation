// =========================
//  Particle Simulation
//  Emre Bilgen
// =========================

#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <random>
#include <chrono>     // timing


// =========================
//        Constants
// =========================

const double G = 1.0, k = 1.0, eps = 0.1;                   // G = Gravitational Constant, k = Coulomb Constant, eps = softening length (keeps forces finite as r -> 0)
const double boxHalf = 10.0;                                // Box spans -boxHalf +boxHalf on both axes (20 x 20)


//=========================
//   Particle Creation
// =========================

// 2D Vector creation
struct Vec2 {
    double x = 0.0;
    double y = 0.0;

    // Vector operations
    Vec2 operator+(const Vec2& other) const {
        return {x + other.x, y + other.y};
    }
    Vec2 operator-(const Vec2& other) const {
        return {x - other.x, y - other.y};
    }
    Vec2 operator*(double scalar) const {
        return {x * scalar, y * scalar};
    }
    Vec2 operator/(double scalar) const {
        return {x / scalar, y / scalar};
    }
    Vec2& operator+=(const Vec2& other) {
        x += other.x;
        y += other.y;
        return *this;
    }
    Vec2& operator-=(const Vec2& other) {
        x -= other.x;
        y -= other.y;
        return *this;
    }
    double lengthSquared() const {
        return x*x+y*y;
    }
    double length() const {
        return std::sqrt(lengthSquared());
    }
    Vec2 normalized() const {
        double len = length();
        if (len < 0.001) return{0, 0};
        double i_hat = x/len, j_hat = y/len;
        return {i_hat, j_hat};
    }
};

// Printing cout in Vec2
std::ostream& operator<<(std::ostream& os, const Vec2& v) {
    return os << "(" << v.x << ", " << v.y << ")";
}

// Dot product
double dot(const Vec2& a, const Vec2& b) {
    return a.x*b.x+a.y*b.y;
}

// Scalar multiplication but with scalar*vector
Vec2 operator*(double s, const Vec2& v){
    return {s*v.x, s*v.y};
}

// Particle
struct Particle {
    double mass = 1.0, charge = 0.0, radius = 1.0, restitution = 1.0;
    Vec2 pos, vel, force;  
};

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
