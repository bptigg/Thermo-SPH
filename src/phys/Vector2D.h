#pragma once
#include <cmath>
#include <iostream>

struct Vector2D {
    double x = 0.0;
    double y = 0.0;

    Vector2D() = default;
    Vector2D(double x, double y) : x(x), y(y) {}

    Vector2D operator+(const Vector2D& rhs) const;
    Vector2D operator-(const Vector2D& rhs) const;
    Vector2D operator-() const; // Unary negation

    Vector2D operator*(double scalar) const;
    Vector2D operator/(double scalar) const;

    Vector2D& operator+=(const Vector2D& rhs);
    Vector2D& operator-=(const Vector2D& rhs);
    Vector2D& operator*=(double scalar);
    Vector2D& operator/=(double scalar);

    double dot(const Vector2D& rhs) const;
    double normSq() const;
    double norm() const;
    Vector2D normalized() const;

    double length() const {
        return std::hypot(x, y);
    }

    double lengthSquared() const {
        return x * x + y * y;
    }
};

Vector2D operator*(double scalar, const Vector2D& vec);