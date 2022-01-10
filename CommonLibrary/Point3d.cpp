#include <srep/Point3d.h>

#include <cmath>
#include <stdexcept>

#include <vtkMath.h>

namespace srep {

Point3d::Point3d()
    : X(0.0)
    , Y(0.0)
    , Z(0.0)
{}

Point3d::Point3d(double x, double y, double z)
    : X(x)
    , Y(y)
    , Z(z)
{
    if (std::isnan(x) || std::isnan(y) || std::isnan(z)) {
        throw std::invalid_argument("Point cannot have a nan component");
    }
}

Point3d::Point3d(const double p[3])
    : Point3d(p[0], p[1], p[2])
{}

Point3d::Point3d(const std::array<double, 3>& p)
    : Point3d(p[0], p[1], p[2])
{}

Point3d::Point3d(const vtkVector3d& p)
    : Point3d(p[0], p[1], p[2])
{}

const double& Point3d::operator[](size_t i) const {
    if (i == 0) {
        return this->X;
    }
    if (i == 1) {
        return this->Y;
    }
    if (i == 2) {
        return this->Z;
    }
    throw std::out_of_range("srep::Vector3d brackets only accept 0, 1, or 2. Found " + std::to_string(i));
}

void Point3d::SetX(double x) {
    if (std::isnan(x)) {
        throw std::invalid_argument("Point cannot have a nan component");
    }
    this->X = x;
}

void Point3d::SetY(double y) {
    if (std::isnan(y)) {
        throw std::invalid_argument("Point cannot have a nan component");
    }
    this->Y = y;
}

void Point3d::SetZ(double z) {
    if (std::isnan(z)) {
        throw std::invalid_argument("Point cannot have a nan component");
    }
    this->Z = z;
}

std::array<double, 3> Point3d::AsArray() const {
    return std::array<double, 3>{this->GetX(), this->GetY(), this->GetZ()};
}

double Point3d::Distance(const Point3d& a, const Point3d& b) {
    const auto distSquared = vtkMath::Distance2BetweenPoints(a.AsArray().data(), b.AsArray().data());
    return pow(distSquared, 0.5);
}

bool operator< (const Point3d& a, const Point3d& b) {
    return (a.GetX() != b.GetX()) ? (a.GetX() < b.GetX())
        : ((a.GetY() != b.GetY()) ? (a.GetY() < b.GetY()) : (a.GetZ() < b.GetZ()));
}
bool operator==(const Point3d& a, const Point3d& b) {
    return !(a < b) && !(b < a);
}
bool operator!=(const Point3d& a, const Point3d& b) {
    return (a < b) || (b < a);
}
bool operator> (const Point3d& a, const Point3d& b) {
    return b < a;
}
bool operator<=(const Point3d& a, const Point3d& b) {
    return !(b < a);
}
bool operator>=(const Point3d& a, const Point3d& b) {
    return !(a < b);
}

std::ostream& operator<<(std::ostream& os, const Point3d& point) {
    os << "(" << point.GetX() << ", " << point.GetY() << ", " << point.GetZ() << ")";
    return os;
}

void PlaceInto(const Point3d& p, vtkVector3d& v) {
    v[0] = p[0];
    v[1] = p[1];
    v[2] = p[2];
}

}
