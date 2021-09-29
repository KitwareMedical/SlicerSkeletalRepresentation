#ifndef __srep_Spoke_h
#define __srep_Spoke_h

#include <srep/Point3d.h>
#include <srep/Vector3d.h>

namespace srep {

class Spoke {
public:
    Spoke();
    /// Create a spoke
    /// \param skeletalPoint Point on the skeleton
    /// \param direction Direction and length of spoke
    Spoke(const Point3d& skeletalPoint, const Vector3d& direction);

    //copy and move defined and valid
    Spoke(const Spoke&) = default;
    Spoke& operator=(const Spoke&) = default;
    Spoke(Spoke&&) = default;
    Spoke& operator=(Spoke&&) = default;
    ~Spoke() = default;

    /// Sets the radius (length) of the spoke.
    void SetRadius(double radius);
    /// Gets the radius (length) of the spoke.
    double GetRadius() const;

    /// Sets the point that functions as the base of the spoke.
    void SetSkeletalPoint(const Point3d& skeletalPoint);
    /// Gets the point that functions as the base of the spoke.
    Point3d GetSkeletalPoint() const;

    /// Sets the direction and radius
    void SetDirectionAndMagnitude(const Vector3d& direction);
    /// Sets the direction while keeping the current radius
    /// @throws std::runtime_error if length of direction is 0
    /// @note direction does not need to be a unit vector
    void SetDirectionOnly(const Vector3d& direction);
    /// \note This is not a unit direction. For that just call GetDirection().Unit().
    Vector3d GetDirection() const;

    /// Gets the point that is at the tip of the spoke.
    /// @note this is obtained by SkeletalPoint + Direction
    Point3d GetBoundaryPoint() const;
private:
    Point3d SkeletalPoint;
    Vector3d Direction;
};

bool operator==(const Spoke& a, const Spoke& b);
bool operator!=(const Spoke& a, const Spoke& b);
bool operator< (const Spoke& a, const Spoke& b);
bool operator> (const Spoke& a, const Spoke& b);
bool operator<=(const Spoke& a, const Spoke& b);
bool operator>=(const Spoke& a, const Spoke& b);

std::ostream& operator<<(std::ostream& os, const Spoke& spoke);

}

#endif
