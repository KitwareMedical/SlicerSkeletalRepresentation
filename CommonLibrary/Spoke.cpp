#include <srep/Spoke.h>

namespace srep {

Spoke::Spoke()
    : SkeletalPoint()
    , Direction()
{}

Spoke::Spoke(const Point3d& skeletalPoint, const Vector3d& direction)
    : SkeletalPoint(skeletalPoint)
    , Direction(direction)
{}

Spoke::Spoke(const Point3d& skeletalPoint, const Point3d& boundaryPoint)
    : SkeletalPoint(skeletalPoint)
    , Direction(skeletalPoint, boundaryPoint)
{}

void Spoke::SetRadius(const double radius) {
    this->Direction.Resize(radius);
}
double Spoke::GetRadius() const {
    return this->Direction.GetLength();
}

void Spoke::SetSkeletalPoint(const Point3d& skeletalPoint) {
    this->SkeletalPoint = skeletalPoint;
}
Point3d Spoke::GetSkeletalPoint() const {
    return this->SkeletalPoint;
}

void Spoke::SetDirectionAndMagnitude(const Vector3d& direction) {
    this->Direction = direction;
}
void Spoke::SetDirectionOnly(const Vector3d& direction) {
    const auto currentRadius = this->GetRadius();
    this->Direction = direction;
    this->Direction.Resize(currentRadius);
}
/// \note This is not a unit direction
Vector3d Spoke::GetDirection() const {
    return this->Direction;
}

Point3d Spoke::GetBoundaryPoint() const {
    return this->SkeletalPoint + this->Direction;
}

bool operator< (const Spoke& a, const Spoke& b) {
    return (a.GetSkeletalPoint() != b.GetSkeletalPoint())
        ? (a.GetSkeletalPoint() < b.GetSkeletalPoint())
        : (a.GetDirection() < b.GetDirection());

}
bool operator==(const Spoke& a, const Spoke& b) {
    return a.GetSkeletalPoint() == b.GetSkeletalPoint()
        && a.GetDirection() == b.GetDirection();
}
bool operator!=(const Spoke& a, const Spoke& b) {
    return !(a == b);
}
bool operator> (const Spoke& a, const Spoke& b) {
    return b < a;
}
bool operator<=(const Spoke& a, const Spoke& b) {
    return !(b < a);
}
bool operator>=(const Spoke& a, const Spoke& b) {
    return !(a < b);
}

std::ostream& operator<<(std::ostream& os, const Spoke& spoke) {
    os << "Spoke{" << spoke.GetSkeletalPoint() << ", " << spoke.GetDirection() << "}";
    return os;
}

}
