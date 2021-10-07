#include <srep/Vector3d.h>

#include <cmath>
#include <stdexcept>

namespace srep {

Vector3d::Vector3d()
    : X(0.0)
    , Y(0.0)
    , Z(0.0)
{}

Vector3d::Vector3d(double x, double y, double z)
    : X(x)
    , Y(y)
    , Z(z)
{
    if (std::isnan(x) || std::isnan(y) || std::isnan(z)) {
        throw std::invalid_argument("Point cannot have a nan component");
    }
}

Vector3d::Vector3d(const double p[3])
    : Vector3d(p[0], p[1], p[2])
{}

void Vector3d::SetX(double x) {
    if (std::isnan(x)) {
        throw std::invalid_argument("Point cannot have a nan component");
    }
    this->X = x;
}

void Vector3d::SetY(double y) {
    if (std::isnan(y)) {
        throw std::invalid_argument("Point cannot have a nan component");
    }
    this->Y = y;
}

void Vector3d::SetZ(double z) {
    if (std::isnan(z)) {
        throw std::invalid_argument("Point cannot have a nan component");
    }
    this->Z = z;
}

std::array<double, 3> Vector3d::AsArray() const {
    return std::array<double, 3>{this->GetX(), this->GetY(), this->GetZ()};
}

Vector3d::Vector3d(const Point3d& from, const Point3d& to)
    : Vector3d(to.GetX() - from.GetX(),
               to.GetY() - from.GetY(),
               to.GetZ() - from.GetZ())
{}

double Vector3d::GetLength() const {
    return std::sqrt(
        std::pow(this->GetX(), 2)
        + std::pow(this->GetY(), 2)
        + std::pow(this->GetZ(), 2));
}

void Vector3d::Resize(const double length) {
    const double oldLength = this->GetLength();
    if (oldLength == 0) {
        throw std::runtime_error("Cannot resize a zero length vector. Unsure which direction to go.");
    }

    const double factor = length / oldLength;
    (*this) = (*this) * factor;
}

Vector3d Vector3d::Unit() const {
    const auto length =this->GetLength();
    if (length == 0.0) {
        throw std::runtime_error("Cannot make unit vector when current length is 0");
    }
    return (*this) / this->GetLength();
}

bool operator< (const Vector3d& a, const Vector3d& b) {
    return (a.GetX() != b.GetX()) ? (a.GetX() < b.GetX())
        : ((a.GetY() != b.GetY()) ? (a.GetY() < b.GetY()) : (a.GetZ() < b.GetZ()));
}
bool operator==(const Vector3d& a, const Vector3d& b) {
    return !(a < b) && !(b < a);
}
bool operator!=(const Vector3d& a, const Vector3d& b) {
    return (a < b) || (b < a);
}
bool operator> (const Vector3d& a, const Vector3d& b) {
    return b < a;
}
bool operator<=(const Vector3d& a, const Vector3d& b) {
    return !(b < a);
}
bool operator>=(const Vector3d& a, const Vector3d& b) {
    return !(a < b);
}

Vector3d operator+(const Vector3d& a, const Vector3d& b) {
    Vector3d c(a.GetX() + b.GetX(),
               a.GetY() + b.GetY(),
               a.GetZ() + b.GetZ());
    return c;
}
Vector3d operator-(const Vector3d& a, const Vector3d& b) {
    Vector3d c(a.GetX() - b.GetX(),
               a.GetY() - b.GetY(),
               a.GetZ() - b.GetZ());
    return c;
}
Vector3d operator*(const Vector3d& v, const double multiplier) {
    Vector3d c(v.GetX() * multiplier,
               v.GetY() * multiplier,
               v.GetZ() * multiplier);
    return c;
}
Vector3d operator/(const Vector3d& v, const double divisor) {
    Vector3d c(v.GetX() / divisor,
               v.GetY() / divisor,
               v.GetZ() / divisor);
    return c;
}

Point3d operator+(const Point3d& a, const Vector3d& b) {
    Point3d c(a.GetX() + b.GetX(),
              a.GetY() + b.GetY(),
              a.GetZ() + b.GetZ());
    return c;
}
Point3d operator-(const Point3d& a, const Vector3d& b) {
    Point3d c(a.GetX() - b.GetX(),
              a.GetY() - b.GetY(),
              a.GetZ() - b.GetZ());
    return c;
}

std::ostream& operator<<(std::ostream& os, const Vector3d& point) {
    os << "<" << point.GetX() << ", " << point.GetY() << ", " << point.GetZ() << ">";
    return os;
}

}
