#include "Vector2D.h"

Vector2D Vector2D::operator+(const Vector2D& rhs) const {
    return Vector2D(x + rhs.x, y + rhs.y);
}

Vector2D Vector2D::operator-(const Vector2D& rhs) const {
    return Vector2D(x - rhs.x, y - rhs.y);
}

Vector2D Vector2D::operator-() const {
    return Vector2D(-x, -y);
}

Vector2D Vector2D::operator*(double scalar) const {
    return Vector2D(x * scalar, y * scalar);
}

Vector2D Vector2D::operator/(double scalar) const {
    return Vector2D(x / scalar, y / scalar);
}

Vector2D& Vector2D::operator+=(const Vector2D& rhs) {
    x += rhs.x;
    y += rhs.y;
    return *this;
}

Vector2D& Vector2D::operator-=(const Vector2D& rhs) {
    x -= rhs.x;
    y -= rhs.y;
    return *this;
}

Vector2D& Vector2D::operator*=(double scalar) {
    x *= scalar;
    y *= scalar;
    return *this;
}

Vector2D& Vector2D::operator/=(double scalar) {
    x /= scalar;
    y /= scalar;
    return *this;
}

double Vector2D::dot(const Vector2D& rhs) const {
    return x * rhs.x + y * rhs.y;
}

double Vector2D::normSq() const {
    return x * x + y * y;
}

double Vector2D::norm() const {
    return std::sqrt(normSq());
}

Vector2D Vector2D::normalized() const {
    double n = norm();
    if (n == 0.0) return Vector2D(0.0, 0.0);
    return *this / n;
}

Vector2D operator*(double scalar, const Vector2D& vec) {
    return vec * scalar;
}