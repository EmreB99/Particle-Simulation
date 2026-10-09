// =====================================================================
//  Particle Simulation
//  Emre Bilgen
// =====================================================================

#include <iostream>
#include <cmath>
#include <vector>
#include <algorithm>

// Constants
const double G = 1.0, k = 1.0, eps = 0.1;  // G = Gravitational Constant, k = Coulomb Constant, eps = softening length (keeps forces finite as r -> 0)
const double boxHalf = 10.0;                // Box spans -boxHalf +boxHalf on both axes (20 x 20)

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

// Gravitational force
Vec2 gravityForce(const Particle& a, const Particle& b) {
    Vec2 r = b.pos - a.pos;
    double d2 = r.lengthSquared() + eps*eps;                // distance squared => d^2 => d2
    Vec2 gForce = ((G*a.mass*b.mass)/(d2 * std::sqrt(d2)))*r;  // Scalar: (G*a.mass*b.mass)/(d2 * std::sqrt(d2)) hidden inside is r_hat, makes gForce a vector
    return gForce;
}

// Electrical force
Vec2 electricForce(const Particle& a, const Particle& b) {
    Vec2 r = b.pos - a.pos;
    double d2 = r.lengthSquared() + eps*eps;
    Vec2 eForce = (-(k*a.charge*b.charge)/(d2 * std::sqrt(d2)))*r; // minus(-) in the numerator because same charges repel, opposites attract
    return eForce;
}

void computeForces(std::vector<Particle>& particles) {
    const size_t N = particles.size();
    
    for(size_t i = 0; i < N; i++){
        particles[i].force = {0,0}; // Reset all forces to 0 before calculation
    }

    for (size_t i = 0; i < N; i++) {
        for(size_t j = i + 1; j < N; j++) {
            Vec2 F = gravityForce(particles[i], particles[j]) + electricForce(particles[i], particles[j]);
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


int main() {

    // Creating particles
    std::vector<Particle> particles;
    //Testing Electrical Force
    Particle a;
    a.pos = {-5, 0};
    a.mass = 1;
    a.charge = +1;
    a.vel = {0, 0};
    particles.push_back(a);
    Particle b;
    b.pos = {5, 0};
    b.mass = 1;
    b.charge = -1;
    b.vel = {0, 0};
    particles.push_back(b);

    //Testing Gravitational Force
        // Particle a;
        // a.pos = {1, 1};
        // a.mass = 2.0;
        // particles.push_back(a);
        // Particle b;
        // b.pos = {2, 0};
        // b.mass = 3.0;
        // particles.push_back(b);
        // Particle c;
        // c.pos = {-1, 3};
        // c.mass = 1.5;
        // particles.push_back(c);
        // Particle d;
        // d.pos = {-3, 2};
        // d.mass = 0.8;
        // particles.push_back(d);

    // Constants of experiment
    const size_t N = particles.size();
    const double dt = 0.001;
    const int numSteps = 60000; //Number of calculations

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
        }
    }
    return 0;
}
