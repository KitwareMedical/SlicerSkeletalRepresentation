#include <vtkSRepSpoke.h>
#include <vtkObjectFactory.h>

//----------------------------------------------------------------------
vtkStandardNewMacro(vtkSRepSpoke);

//----------------------------------------------------------------------
vtkSRepSpoke::vtkSRepSpoke() = default;

//----------------------------------------------------------------------
vtkSRepSpoke::~vtkSRepSpoke() = default;

//----------------------------------------------------------------------
double vtkSRepSpoke::GetRadius() const {
  return this->Direction.GetLength();
}

//----------------------------------------------------------------------
void vtkSRepSpoke::SetRadius(const double radius) {
    this->Direction.Resize(radius);
    this->Modified();
}

//----------------------------------------------------------------------
void vtkSRepSpoke::SetSkeletalPoint(const srep::Point3d& skeletalPoint) {
    this->SkeletalPoint = skeletalPoint;
    this->Modified();
}

//----------------------------------------------------------------------
void vtkSRepSpoke::SetSkeletalPoint(const vtkVector3d& skeletalPoint) {
  this->SetSkeletalPoint(srep::Point3d(skeletalPoint));
}

//----------------------------------------------------------------------
srep::Point3d vtkSRepSpoke::GetSkeletalPoint() const {
    return this->SkeletalPoint;
}

//----------------------------------------------------------------------
void vtkSRepSpoke::GetSkeletalPoint(vtkVector3d& out) const {
  srep::PlaceInto(this->GetSkeletalPoint(), out);
}

//----------------------------------------------------------------------
void vtkSRepSpoke::SetDirectionAndMagnitude(const srep::Vector3d& direction) {
  this->Direction = direction;
  this->Modified();
}

//----------------------------------------------------------------------
void vtkSRepSpoke::SetDirectionAndMagnitude(const vtkVector3d& direction) {
  this->SetDirectionAndMagnitude(srep::Vector3d(direction));
}

//----------------------------------------------------------------------
void vtkSRepSpoke::SetDirectionOnly(const srep::Vector3d& direction) {
  const auto currentRadius = this->GetRadius();
  this->Direction = direction;
  this->Direction.Resize(currentRadius);
  this->Modified();
}

//----------------------------------------------------------------------
void vtkSRepSpoke::SetDirectionOnly(const vtkVector3d& direction) {
  this->SetDirectionOnly(srep::Vector3d(direction));
}

//----------------------------------------------------------------------
srep::Vector3d vtkSRepSpoke::GetDirection() const {
  return this->Direction;
}

//----------------------------------------------------------------------
void vtkSRepSpoke::GetDirection(vtkVector3d& out) const {
  srep::PlaceInto(this->GetDirection(), out);
}

//----------------------------------------------------------------------
srep::Point3d vtkSRepSpoke::GetBoundaryPoint() const {
  return this->SkeletalPoint + this->Direction;
}

//----------------------------------------------------------------------
void vtkSRepSpoke::GetBoundaryPoint(vtkVector3d& out) const {
  srep::PlaceInto(this->GetBoundaryPoint(), out);
}

//----------------------------------------------------------------------
vtkSRepSpoke* vtkSRepSpoke::Create(const srep::Point3d& skeletalPoint, const srep::Vector3d& direction) {
  auto spoke = vtkSRepSpoke::New();
  spoke->SetSkeletalPoint(skeletalPoint);
  spoke->SetDirectionAndMagnitude(direction);
  return spoke;
}

//----------------------------------------------------------------------
vtkSRepSpoke* vtkSRepSpoke::Create(const srep::Point3d& skeletalPoint, const srep::Point3d& boundaryPoint) {
  auto spoke = vtkSRepSpoke::New();
  spoke->SetSkeletalPoint(skeletalPoint);
  spoke->SetDirectionAndMagnitude(srep::Vector3d(skeletalPoint, boundaryPoint));
  return spoke;
}

//----------------------------------------------------------------------
vtkSmartPointer<vtkSRepSpoke> vtkSRepSpoke::SmartCreate(const srep::Point3d& skeletalPoint, const srep::Vector3d& direction) {
  return vtkSmartPointer<vtkSRepSpoke>::Take(Create(skeletalPoint, direction));
}

//----------------------------------------------------------------------
vtkSmartPointer<vtkSRepSpoke> vtkSRepSpoke::SmartCreate(const srep::Point3d& skeletalPoint, const srep::Point3d& boundaryPoint) {
  return vtkSmartPointer<vtkSRepSpoke>::Take(Create(skeletalPoint, boundaryPoint));
}

//----------------------------------------------------------------------
vtkSRepSpoke* vtkSRepSpoke::CreatePointToPoint(const vtkVector3d& skeletalPoint, const vtkVector3d& boundaryPoint) {
  return Create(srep::Point3d(skeletalPoint), srep::Point3d(boundaryPoint));
}

//----------------------------------------------------------------------
vtkSRepSpoke* vtkSRepSpoke::CreatePointAndDirection(const vtkVector3d& skeletalPoint, const vtkVector3d& direction) {
  return Create(srep::Point3d(skeletalPoint), srep::Vector3d(direction));
}

//----------------------------------------------------------------------
vtkSRepSpoke* vtkSRepSpoke::Clone() const {
  return Create(this->SkeletalPoint, this->Direction);
}

//----------------------------------------------------------------------
vtkSmartPointer<vtkSRepSpoke> vtkSRepSpoke::SmartClone() const {
  return vtkSmartPointer<vtkSRepSpoke>::Take(Clone());
}

//----------------------------------------------------------------------
void vtkSRepSpoke::PrintSelf(std::ostream& os, vtkIndent indent) {
  this->Superclass::PrintSelf(os, indent);
  os << indent << *this;
}

//----------------------------------------------------------------------
std::ostream& operator<<(std::ostream& os, const vtkSRepSpoke& spoke) {
  os << "Spoke{" << spoke.GetSkeletalPoint() << ", " << spoke.GetDirection() << "}";
  return os;
}
