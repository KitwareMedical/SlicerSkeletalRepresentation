#ifndef __srep_Point3d_h
#define __srep_Point3d_h

#include <array>
#include <iostream>
#include <vtkVector.h>
namespace srep {

class Point3d {
public:
    /// Construct new point at origin.
    Point3d();
    /// @throws std::invalid_argument if any nan
    Point3d(double x, double y, double z);

    //I would rather not have this, but may be convienent for interfacing with VTK
    /// @throws std::invalid_argument if any nan
    explicit Point3d(const double p[3]);

    explicit Point3d(const std::array<double, 3>& p);

    explicit Point3d(const vtkVector3d& p);

    //copy and move defined and valid
    Point3d(const Point3d&) = default;
    Point3d& operator=(const Point3d&) = default;
    Point3d(Point3d&&) = default;
    Point3d& operator=(Point3d&&) = default;
    ~Point3d() = default;

    const double& operator[](size_t i) const;

    /// Gets the X component of the point.
    inline double GetX() const {
        return this->X;
    }
    /// Gets the Y component of the point.
    inline double GetY() const {
        return this->Y;
    }
    /// Gets the Z component of the point.
    inline double GetZ() const {
        return this->Z;
    }

    /// Sets the X component of the point.
    /// @throws std::invalid_argument if x is nan.
    void SetX(double x);
    /// Sets the X component of the point.
    /// @throws std::invalid_argument if y is nan.
    void SetY(double y);
    /// Sets the X component of the point.
    /// @throws std::invalid_argument if z is nan.
    void SetZ(double z);

    /// Returns the point as a length 3 array;
    std::array<double, 3> AsArray() const;

    /// Returns the distance from one point to another.
    static double Distance(const Point3d& a, const Point3d& b);
private:
    double X;
    double Y;
    double Z;
};
bool operator==(const Point3d& a, const Point3d& b);
bool operator!=(const Point3d& a, const Point3d& b);
bool operator< (const Point3d& a, const Point3d& b);
bool operator> (const Point3d& a, const Point3d& b);
bool operator<=(const Point3d& a, const Point3d& b);
bool operator>=(const Point3d& a, const Point3d& b);

std::ostream& operator<<(std::ostream& os, const Point3d& point);

void PlaceInto(const Point3d& p, vtkVector3d& v);


}

#endif
