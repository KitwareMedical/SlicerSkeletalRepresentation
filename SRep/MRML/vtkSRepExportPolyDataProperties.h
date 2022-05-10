#ifndef __vtkSRepExportPolyDataProperties_h
#define __vtkSRepExportPolyDataProperties_h

#include <vtkDataArray.h>
#include <vtkObject.h>
#include <vtkSmartPointer.h>

#include "vtkSlicerSRepModuleMRMLExport.h"

class VTK_SLICER_SREP_MODULE_MRML_EXPORT vtkSRepExportPolyDataProperties : public vtkObject {
public:
  static const int UpBoundaryPointType = 0;
  static const int UpSkeletalPointType = 1;
  static const int DownBoundaryPointType = 2;
  static const int DownSkeletalPointType = 3;
  static const int CrestBoundaryPointType = 4;
  static const int CrestSkeletalPointType = 5;
  static const int UpSpokeLineType = 6;
  static const int DownSpokeLineType = 7;
  static const int CrestSpokeLineType = 8;
  static const int CrestCurveLineType = 9;
  static const int SkeletalSheetLineType = 10;
  static const int SkeletonToCrestConnectionLineType = 11;
  static const int SpineLineType = 12;
  static const int NumberOfTypes = 13;

  static vtkSRepExportPolyDataProperties *New();

  ~vtkSRepExportPolyDataProperties() = default;

  /// Standard methods for a VTK class.
  vtkTypeMacro(vtkSRepExportPolyDataProperties, vtkObject);

  void PrintSelf(ostream& os, vtkIndent indent) override;

  void SetSRepDataArray(vtkDataArray* array);
  vtkDataArray* GetSRepDataArray() const;

  void SetPointTypeArrayName(const std::string& name);
  std::string GetPointTypeArrayName() const;

  void SetLineTypeArrayName(const std::string& name);
  std::string GetLineTypeArrayName() const;

  void SetIncludeUpSpokes(bool include);
  bool GetIncludeUpSpokes() const;

  void SetIncludeDownSpokes(bool include);
  bool GetIncludeDownSpokes() const;

  void SetIncludeCrestSpokes(bool include);
  bool GetIncludeCrestSpokes() const;

  void SetIncludeCrestCurve(bool include);
  bool GetIncludeCrestCurve() const;

  void SetIncludeSkeletalSheet(bool include);
  bool GetIncludeSkeletalSheet() const;

  void SetIncludeSkeletonToCrestConnection(bool include);
  bool GetIncludeSkeletonToCrestConnection() const;

  void SetIncludeSpine(bool include);
  bool GetIncludeSpine() const;
protected:
  vtkSRepExportPolyDataProperties() = default;
  vtkSRepExportPolyDataProperties(const vtkSRepExportPolyDataProperties&) = delete;
  vtkSRepExportPolyDataProperties(vtkSRepExportPolyDataProperties&&) = delete;
  vtkSRepExportPolyDataProperties& operator=(const vtkSRepExportPolyDataProperties&) = delete;
  vtkSRepExportPolyDataProperties& operator=(vtkSRepExportPolyDataProperties&&) = delete;
private:
  bool IncludeUpSpokes = true;
  bool IncludeDownSpokes = true;
  bool IncludeCrestSpokes = true;
  bool IncludeCrestCurve = true;
  bool IncludeSkeletalSheet = true;
  bool IncludeSkeletonToCrestConnection = true;
  bool IncludeSpine = true;
  std::string PointTypeArrayName = "SRepPointDataArray";
  std::string LineTypeArrayName = "SRepLineDataArray";
  vtkSmartPointer<vtkDataArray> SRepDataArray;
};

#endif
