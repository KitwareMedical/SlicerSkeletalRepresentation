#ifndef __srep_Vector3d_h
#define __srep_Vector3d_h

#include <iostream>
#include <vector>

#include <srep/Point3d.h>

namespace srep {

class Vector3d {
public:
    Vector3d();
    /// Construct a vector
    /// @throws std::invalid_argument if x, y, or z is nan
    Vector3d(double x, double y, double z);

    /// Construct vector going from point to another, both direction and magnitude.
    Vector3d(const Point3d& from, const Point3d& to);

    //I would rather not have this, but may be convienent for interfacing with VTK
    /// Construct a vector
    /// @throws std::invalid_argument if any nan
    explicit Vector3d(const double p[3]);

    explicit Vector3d(const std::array<double, 3>& p);

    explicit Vector3d(const vtkVector3d& p);

    //copy and move defined and valid
    Vector3d(const Vector3d&) = default;
    Vector3d& operator=(const Vector3d&) = default;
    Vector3d(Vector3d&&) = default;
    Vector3d& operator=(Vector3d&&) = default;
    ~Vector3d() = default;

    const double& operator[](size_t i) const;

    /// Gets the X component of the vector.
    inline double GetX() const {
        return this->X;
    }
    /// Gets the Y component of the vector.
    inline double GetY() const {
        return this->Y;
    }
    /// Gets the Z component of the vector.
    inline double GetZ() const {
        return this->Z;
    }

    /// Sets the X component of the vector.
    /// @throws std::invalid_argument if x is nan.
    void SetX(double x);
    /// Sets the Y component of the vector.
    /// @throws std::invalid_argument if y is nan.
    void SetY(double y);
    /// Sets the Z component of the vector.
    /// @throws std::invalid_argument if z is nan.
    void SetZ(double z);

    /// Gets the length of the vector.
    double GetLength() const;

    /// Adjusts this vector to be given length size in same direction.
    /// @throws std::runtime_error if this is currently a zero length vector.
    void Resize(double length);

    /// Adjusts given vector to be given length size in same direction.
    /// @throws std::runtime_error if v is a zero length vector.
    inline static Vector3d Resize(const Vector3d& v, double length) {
        Vector3d v2(v);
        v2.Resize(length);
        return v2;
    }

    /// Returns a unit vector in this direction.
    Vector3d Unit() const;

    /// Returns the point as a length 3 array;
    std::array<double, 3> AsArray() const;
private:
    double X;
    double Y;
    double Z;
};

bool operator==(const Vector3d& a, const Vector3d& b);
bool operator!=(const Vector3d& a, const Vector3d& b);
bool operator< (const Vector3d& a, const Vector3d& b);
bool operator> (const Vector3d& a, const Vector3d& b);
bool operator<=(const Vector3d& a, const Vector3d& b);
bool operator>=(const Vector3d& a, const Vector3d& b);

/// Adds two vectors by adding x, y, and z components.
Vector3d operator+(const Vector3d& a, const Vector3d& b);

/// Subtracts two vectors by subtracting x, y, and z components.
Vector3d operator-(const Vector3d& a, const Vector3d& b);

/// Multiplies vector magnitude by multiplier.
Vector3d operator*(const Vector3d& v, double multiplier);

/// Divides vector magnitude by divisor.
Vector3d operator/(const Vector3d& v, double divisor);

/// Gets point at end of vector if vector, b, started at point, a.
Point3d operator+(const Point3d& a, const Vector3d& b);

/// Gets point at end of vector if vector, -b, started at point, a.
Point3d operator-(const Point3d& a, const Vector3d& b);

std::ostream& operator<<(std::ostream& os, const Vector3d& point);

void PlaceInto(const Vector3d& v1, vtkVector3d& v2);

}

#endif
