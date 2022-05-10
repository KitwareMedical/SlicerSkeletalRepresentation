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
  /// The different spoke orientations a skeletal point may have.
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

  /// @{
  /// Gets the spoke in the up orientation. Never returns nullptr.
  const vtkSRepSpoke* GetUpSpoke() const;
  vtkSRepSpoke* GetUpSpoke();
  /// @}

  /// @{
  /// Gets the spoke in the down orientation. Never returns nullptr.
  const vtkSRepSpoke* GetDownSpoke() const;
  vtkSRepSpoke* GetDownSpoke();
  /// @}

  /// @{
  /// Gets the spoke in the crest orientation, if it exists.
  /// \returns Spoke if IsCrest returns true, otherwise nullptr
  const vtkSRepSpoke* GetCrestSpoke() const;
  vtkSRepSpoke* GetCrestSpoke();
  /// @}

  /// Returns whether this skeletal point is has a crest spoke.
  bool IsCrest() const;

  /// @{
  /// Sets the spoke to the given spoke. Does a shallow copy of spoke.
  /// \throws std::invalid_argument if spoke is nullptr.
  /// \sa SetSpoke
  void SetUpSpoke(vtkSRepSpoke* spoke)
    VTK_EXPECTS(nullptr != spoke);
  void SetDownSpoke(vtkSRepSpoke* spoke)
    VTK_EXPECTS(nullptr != spoke);
  /// @}

  /// Sets the spoke to the given spoke. Does a shallow copy of spoke.
  ///
  /// Spoke may be nullptr.
  void SetCrestSpoke(vtkSRepSpoke* spoke);

  /// @{
  /// Gets the spoke for the given orientation, if it exists.
  ///
  /// If spokeType is UpOrientation or DownOrientation, the returned pointer
  /// will never be nullptr.
  /// \sa GetUpSpoke, GetDownSpoke, GetCrestSpoke
  const vtkSRepSpoke* GetSpoke(SpokeOrientation spokeType) const;
  vtkSRepSpoke* GetSpoke(SpokeOrientation spokeType);
  /// @}

  /// Sets the spoke in the given orientation.
  /// \throws std::invalid_argument if spokeType is UpOrientation or DownOrientation and spoke is nullptr.
  /// \sa SetUpSpoke, SetDownSpoke, SetCrestSpoke
  void SetSpoke(SpokeOrientation spokeType, vtkSRepSpoke* spoke)
    VTK_EXPECTS(spokeType == CrestOrientation || spoke != nullptr);

  /// @{
  /// Creates a new SkeletalPoint. Shallow copies the given spokes.
  ///
  /// To do a deep copy during create do
  /// `vtkSRepSkeletalPoint::Create(spoke1.Clone(), spoke2.Clone())`
  static VTK_NEWINSTANCE vtkSRepSkeletalPoint* Create(vtkSRepSpoke* upSpoke, vtkSRepSpoke* downSpoke, vtkSRepSpoke* crestSpoke = nullptr)
    VTK_EXPECTS(upSpoke != nullptr && downSpoke != nullptr);
  static vtkSmartPointer<vtkSRepSkeletalPoint> SmartCreate(vtkSRepSpoke* upSpoke, vtkSRepSpoke* downSpoke, vtkSRepSpoke* crestSpoke = nullptr)
    VTK_EXPECTS(upSpoke != nullptr && downSpoke != nullptr);
  /// @}

  /// @{
  /// Warning: If one derives from vtkSRepSkeletalPoint (please don't) these calls will slice.
  /// This only clones skeletal point data. Any observers or similar are ignored.
  /// This is a deep copy with no connection between "this" object and the output.
  VTK_NEWINSTANCE vtkSRepSkeletalPoint* Clone() const;
  vtkSmartPointer<vtkSRepSkeletalPoint> SmartClone() const;
  /// @}

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
