#pragma once
#include <iostream>
#include <cmath>

// =========================
//   2D Vector Operations
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
inline std::ostream& operator<<(std::ostream& os, const Vec2& v) {
    return os << "(" << v.x << ", " << v.y << ")";
}

// Dot product
inline double dot(const Vec2& a, const Vec2& b) { 
    return a.x*b.x+a.y*b.y;
}

// Scalar multiplication but with scalar*vector
inline Vec2 operator*(double s, const Vec2& v) { 
    return {s*v.x, s*v.y};
}