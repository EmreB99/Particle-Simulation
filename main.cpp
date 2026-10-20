// =========================
//  Particle Simulation
//  Emre Bilgen
// =========================
//
// Build & run (from this folder):
//   clang++ -std=c++17 -O2 -Wall -Wextra *.cpp $(pkg-config --cflags --libs raylib) -o sim
//   ./sim

#include <cmath>
#include <random>
#include <vector>
#include "constants.h"
#include "particle.h"
#include "physics.h"
#include "raylib.h"

// =========================
//   Screen settings
// =========================

const int   screenSize = 800;                                   // window is screenSize x screenSize pixels
const float scale      = screenSize / (2.0f * (float)boxHalf);  // pixels per world unit (800 / 20 = 40)

// World (-boxHalf..+boxHalf, y up)  ->  screen (0..screenSize, y down)
Vector2 worldToScreen(const Vec2& world) {
    Vector2 s;
    s.x = (float)((world.x + boxHalf) * scale);   // shift so -boxHalf -> 0, then scale
    s.y = (float)((boxHalf - world.y) * scale);   // distance BELOW the top (+boxHalf): flips y
    return s;
}

// Color by charge: red = positive, blue = negative, white = (almost) neutral
Color chargeColor(double charge) {
    if (charge >  0.05) return RED;
    if (charge < -0.05) return BLUE;
    return WHITE;
}

int main() {
    // Particles
    std::mt19937 rng(102);                                       // the generator (seed 42)
    std::vector<Particle> particles = createParticles(50, rng);

    // Time stepping
    const double dt = 0.001;                                    // physics step
    const int stepsPerFrame = 17;                               // (1/60 s) / dt ≈ 17 -> real-time speed
    double simTime = 0.0;                                       // simulated seconds so far

    const double E0 = totalEnergy(particles);                   // energy at t = 0, to measure drift against

    // Window
    InitWindow(screenSize, screenSize, "Particle Simulation");
    SetTargetFPS(60);

    while (!WindowShouldClose()) {                              // until Esc or the window is closed
        // 1. CALCULATE PHYSICS
        for (int step = 0; step < stepsPerFrame; step++) {
            computeForces(particles);
            integrate(particles, dt);
            handleCollisions(particles);
            handleWalls(particles);
        }
        simTime += stepsPerFrame * dt;

        // 2. DRAWING
        BeginDrawing();
        ClearBackground(BLACK);

        for (const Particle& p : particles) {                  // "for each particle p in particles"
            DrawCircleV(worldToScreen(p.pos), (float)(p.radius * scale), chargeColor(p.charge));
        }

        // Info text: time and energy drift (once per frame)
        // double E = totalEnergy(particles);
        // DrawText(TextFormat("t = %.1f s   energy drift = %.2f %%", simTime, 100.0 * (E - E0) / std::fabs(E0)),
        //          10, 10, 20, RAYWHITE);
        // DrawFPS(10, 35);

        EndDrawing();
    }

    CloseWindow();
    return 0;
}
