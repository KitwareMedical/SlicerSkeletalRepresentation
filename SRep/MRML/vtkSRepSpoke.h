#ifndef __vtkSRepSpoke_h
#define __vtkSRepSpoke_h

#include <srep/Point3d.h>
#include <srep/Vector3d.h>
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

  double GetRadius() const;
  void SetRadius(double radius);

  srep::Point3d GetSkeletalPoint() const;
  void GetSkeletalPoint(vtkVector3d& v) const;
  void SetSkeletalPoint(const srep::Point3d& point);
  void SetSkeletalPoint(const vtkVector3d& point);

  srep::Vector3d GetDirection() const;
  void GetDirection(vtkVector3d& v) const;
  void SetDirectionOnly(const srep::Vector3d& direction);
  void SetDirectionOnly(const vtkVector3d& direction);
  void SetDirectionAndMagnitude(const srep::Vector3d& direction);
  void SetDirectionAndMagnitude(const vtkVector3d& direction);

  srep::Point3d GetBoundaryPoint() const;
  void GetBoundaryPoint(vtkVector3d& v) const;

  // factory functions that act like constructors
  // The Smart* functions are convenience for use with "auto" keyword
  static VTK_NEWINSTANCE vtkSRepSpoke* Create(const srep::Point3d& skeletalPoint, const srep::Vector3d& direction);
  static VTK_NEWINSTANCE vtkSRepSpoke* Create(const srep::Point3d& skeletalPoint, const srep::Point3d& boundaryPoint);
  static vtkSmartPointer<vtkSRepSpoke> SmartCreate(const srep::Point3d& skeletalPoint, const srep::Vector3d& direction);
  static vtkSmartPointer<vtkSRepSpoke> SmartCreate(const srep::Point3d& skeletalPoint, const srep::Point3d& boundaryPoint);
  static VTK_NEWINSTANCE vtkSRepSpoke* CreatePointToPoint(const vtkVector3d& skeletalPoint, const vtkVector3d& boundaryPoint);
  static VTK_NEWINSTANCE vtkSRepSpoke* CreatePointAndDirection(const vtkVector3d& skeletalPoint, const vtkVector3d& direction);

  /// Warning: If one derives from vtkSRepSpoke (please don't) these calls will slice
  /// This only clones spoke data. Any observers or similar are ignored
  VTK_NEWINSTANCE vtkSRepSpoke* Clone() const;
  vtkSmartPointer<vtkSRepSpoke> SmartClone() const;
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
