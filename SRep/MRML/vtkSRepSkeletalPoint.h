#ifndef __vtkSRepSkeletalPoint_h
#define __vtkSRepSkeletalPoint_h

#include <vtkObject.h>
#include <vtkSmartPointer.h>

#include "vtkSRepSpoke.h"

#include "vtkSlicerSRepModuleMRMLExport.h"

class VTK_SLICER_SREP_MODULE_MRML_EXPORT vtkSRepSkeletalPoint
  : public vtkObject
{
public:
  // only using enum instead of enum class because the VTK Python wrapping
  // for the VTK version in SlicerSALT didn't support enum class when this was written
  enum SpokeOrientation {
    UpOrientation,
    DownOrientation,
    CrestOrientation
  };

  /// Instantiate this class.
  static vtkSRepSkeletalPoint *New();

  ~vtkSRepSkeletalPoint();

  /// Standard methods for a VTK class.
  vtkTypeMacro(vtkSRepSkeletalPoint, vtkObject);

  void PrintSelf(std::ostream& os, vtkIndent indent) override;

  const vtkSRepSpoke* GetUpSpoke() const;
  vtkSRepSpoke* GetUpSpoke();
  const vtkSRepSpoke* GetDownSpoke() const;
  vtkSRepSpoke* GetDownSpoke();
  /// @{
  /// \returns Spoke if IsCrest returns true, otherwise nullptr
  const vtkSRepSpoke* GetCrestSpoke() const;
  vtkSRepSpoke* GetCrestSpoke();
  /// @}
  bool IsCrest() const;

  /// @{
  /// Sets the spoke to the given spoke. Does a shallow copy of spoke.
  /// \sa SetSpoke
  void SetUpSpoke(vtkSRepSpoke* spoke)
    VTK_EXPECTS(nullptr != spoke);
  void SetDownSpoke(vtkSRepSpoke* spoke)
    VTK_EXPECTS(nullptr != spoke);
  void SetCrestSpoke(vtkSRepSpoke* spoke);
  /// @}

  const vtkSRepSpoke* GetSpoke(SpokeOrientation spokeType) const;
  vtkSRepSpoke* GetSpoke(SpokeOrientation spokeType);
  /// \sa SetUpSpoke, SetDownSpoke, SetCrestSpoke
  void SetSpoke(SpokeOrientation spokeType, vtkSRepSpoke* spoke)
    VTK_EXPECTS(spokeType == CrestOrientation || spoke != nullptr);

  /// @{
  /// Creates a new SkeletalPoint. Shallow copies the given spokes.
  /// To essentially do a deep copy during create do
  /// `vtkSRepSkeletalPoint::Create(spoke1.Clone(), spoke2.Clone())`
  static VTK_NEWINSTANCE vtkSRepSkeletalPoint* Create(vtkSRepSpoke* upSpoke, vtkSRepSpoke* downSpoke, vtkSRepSpoke* crestSpoke = nullptr)
    VTK_EXPECTS(upSpoke != nullptr && downSpoke != nullptr);
  static vtkSmartPointer<vtkSRepSkeletalPoint> SmartCreate(vtkSRepSpoke* upSpoke, vtkSRepSpoke* downSpoke, vtkSRepSpoke* crestSpoke = nullptr)
    VTK_EXPECTS(upSpoke != nullptr && downSpoke != nullptr);
  /// @}

  /// Warning: If one derives from vtkSRepSkeletalPoint (please don't) these calls will slice
  /// This only clones skeletal point data. Any observers or similar are ignored
  /// This is a deep copy with no connection between "this" object and the output
  VTK_NEWINSTANCE vtkSRepSkeletalPoint* Clone() const;
  vtkSmartPointer<vtkSRepSkeletalPoint> SmartClone() const;

protected:
  vtkSRepSkeletalPoint();
  vtkSRepSkeletalPoint(const vtkSRepSkeletalPoint&) = delete;
  vtkSRepSkeletalPoint(vtkSRepSkeletalPoint&&) = delete;
  vtkSRepSkeletalPoint& operator=(const vtkSRepSkeletalPoint&) = delete;
  vtkSRepSkeletalPoint& operator=(vtkSRepSkeletalPoint&&) = delete;
private:
  vtkSmartPointer<vtkSRepSpoke> UpSpoke;
  unsigned long UpObservationTag;
  vtkSmartPointer<vtkSRepSpoke> DownSpoke;
  unsigned long DownObservationTag;
  vtkSmartPointer<vtkSRepSpoke> CrestSpoke;
  unsigned long CrestObservationTag;

  void onSpokeModified(vtkObject *caller, unsigned long event, void* callData);
};

std::ostream& operator<<(std::ostream& os, const vtkSRepSkeletalPoint& spokeType);


#endif
