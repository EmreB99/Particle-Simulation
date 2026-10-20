#pragma once

// =========================
//        Constants
// =========================

// Physical constants (simulation units)
constexpr double G   = 0.05;     // gravitational constant
constexpr double k   = 1.0;     // Coulomb constant
constexpr double eps = 0.1;     // softening length (keeps forces finite as r -> 0)

// Box: spans -boxHalf .. +boxHalf on both axes
constexpr double boxHalf = 10.0;