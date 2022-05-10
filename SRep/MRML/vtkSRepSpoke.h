#ifndef __vtkSRepSpoke_h
#define __vtkSRepSpoke_h

#include "srepPoint3d.h"
#include "srepVector3d.h"
#include <vtkObject.h>
#include <vtkSmartPointer.h>
#include <vtkVector.h>

#include "vtkSlicerSRepModuleMRMLExport.h"

class VTK_SLICER_SREP_MODULE_MRML_EXPORT vtkSRepSpoke : public vtkObject {
public:
  /// Instantiate this class.
  static vtkSRepSpoke *New();

  ~vtkSRepSpoke();

  /// Standard methods for a VTK class.
  vtkTypeMacro(vtkSRepSpoke, vtkObject);

  void PrintSelf(ostream& os, vtkIndent indent) override;

  /// Gets the radius (length) of the spoke.
  double GetRadius() const;
  /// Sets the radius (length) of the spoke.
  void SetRadius(double radius);

  /// Gets the point that functions as the base of the spoke.
  srep::Point3d GetSkeletalPoint() const;
  /// Gets the point that functions as the base of the spoke.
  /// \param[out] v Vector to copy the point into.
  void GetSkeletalPoint(vtkVector3d& v) const;
  /// Sets the point that functions as the base of the spoke.
  void SetSkeletalPoint(const srep::Point3d& point);
  /// Sets the point that functions as the base of the spoke.
  void SetSkeletalPoint(const vtkVector3d& point);

  /// Gets the direction and magnitude of the spoke.
  ///
  /// \note This is not a unit direction. For that just call spoke.GetDirection().Unit().
  srep::Vector3d GetDirection() const;

  /// Gets the direction and magnitude of the spoke.
  /// \param[out] v Vector to copy the direction into.
  /// \note This is not a unit direction. For that call vtkVector3d::Normalize after calling GetDirection.
  void GetDirection(vtkVector3d& v) const;

  /// @{
  /// Sets the direction while keeping the current radius
  /// \throws std::runtime_error if length of direction is 0
  /// \note direction does not need to be a unit vector
  void SetDirectionOnly(const srep::Vector3d& direction);
  void SetDirectionOnly(const vtkVector3d& direction);
  /// @}

  /// Sets the direction and radius
  void SetDirectionAndMagnitude(const srep::Vector3d& direction);
  /// Sets the direction and radius
  void SetDirectionAndMagnitude(const vtkVector3d& direction);

  /// Gets the point that is at the tip of the spoke.
  /// \note this is obtained by SkeletalPoint + Direction
  srep::Point3d GetBoundaryPoint() const;

  /// Gets the point that is at the tip of the spoke.
  /// \param[out] v Vector to copy the point into.
  /// \note this is obtained by SkeletalPoint + Direction
  void GetBoundaryPoint(vtkVector3d& v) const;

  /// @{
  /// Factory functions
  /// \param skeletalPoint Point that acts as the base of the spoke.
  /// \param direction The direction and magnitude of the spoke.
  /// \param boundaryPoint Point that is the boundary of the spoke.
  /// \note The "Smart" versions of the functions return a smart pointer for better RAII and use with auto.
  static VTK_NEWINSTANCE vtkSRepSpoke* Create(const srep::Point3d& skeletalPoint, const srep::Vector3d& direction);
  static VTK_NEWINSTANCE vtkSRepSpoke* Create(const srep::Point3d& skeletalPoint, const srep::Point3d& boundaryPoint);
  static vtkSmartPointer<vtkSRepSpoke> SmartCreate(const srep::Point3d& skeletalPoint, const srep::Vector3d& direction);
  static vtkSmartPointer<vtkSRepSpoke> SmartCreate(const srep::Point3d& skeletalPoint, const srep::Point3d& boundaryPoint);
  static VTK_NEWINSTANCE vtkSRepSpoke* CreatePointToPoint(const vtkVector3d& skeletalPoint, const vtkVector3d& boundaryPoint);
  static VTK_NEWINSTANCE vtkSRepSpoke* CreatePointAndDirection(const vtkVector3d& skeletalPoint, const vtkVector3d& direction);
  /// @}

  /// @{
  /// Warning: If one derives from vtkSRepSpoke (please don't) these calls will slice.
  /// This only clones spoke data. Any observers or similar are ignored.
  /// This is a deep copy with no connection between "this" object and the output.
  VTK_NEWINSTANCE vtkSRepSpoke* Clone() const;
  vtkSmartPointer<vtkSRepSpoke> SmartClone() const;
  /// @}
protected:
  vtkSRepSpoke();
  vtkSRepSpoke(const vtkSRepSpoke&) = delete;
  vtkSRepSpoke(vtkSRepSpoke&&) = delete;
  vtkSRepSpoke& operator=(const vtkSRepSpoke&) = delete;
  vtkSRepSpoke& operator=(vtkSRepSpoke&&) = delete;
private:
  srep::Point3d SkeletalPoint;
  srep::Vector3d Direction;
};

std::ostream& operator<<(std::ostream& os, const vtkSRepSpoke& spoke);

#endif
