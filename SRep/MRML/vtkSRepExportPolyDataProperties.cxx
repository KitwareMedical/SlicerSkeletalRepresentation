#include "vtkSRepExportPolyDataProperties.h"
#include <vtkObjectFactory.h>

//----------------------------------------------------------------------
vtkStandardNewMacro(vtkSRepExportPolyDataProperties);

//----------------------------------------------------------------------
void vtkSRepExportPolyDataProperties::PrintSelf(ostream& os, vtkIndent indent) {
  os << indent << "vtkSRepExportPolyDataProperties {" << std::endl
     << indent << "  IncludeUpSpokes: " << IncludeUpSpokes << std::endl
     << indent << "  IncludeDownSpokes: " << IncludeDownSpokes << std::endl
     << indent << "  IncludeCrestSpokes: " << IncludeCrestSpokes << std::endl
     << indent << "  IncludeCrestCurve: " << IncludeCrestCurve << std::endl
     << indent << "  IncludeSkeletalSheet: " << IncludeSkeletalSheet << std::endl
     << indent << "  IncludeSkeletonToCrestConnection: " << IncludeSkeletonToCrestConnection << std::endl
     << indent << "  IncludeSpine: " << IncludeSpine << std::endl
     << indent << "}";
}

//----------------------------------------------------------------------
void vtkSRepExportPolyDataProperties::SetSRepDataArray(vtkDataArray* name) {
  if (this->SRepDataArray != name) {
    this->SRepDataArray = name;
    this->Modified();
  }
}

//----------------------------------------------------------------------
vtkDataArray* vtkSRepExportPolyDataProperties::GetSRepDataArray() const {
  return this->SRepDataArray;
}

//----------------------------------------------------------------------
void vtkSRepExportPolyDataProperties::SetPointTypeArrayName(const std::string& name) {
  if (this->PointTypeArrayName != name) {
    this->PointTypeArrayName = name;
    this->Modified();
  }
}

//----------------------------------------------------------------------
std::string vtkSRepExportPolyDataProperties::GetPointTypeArrayName() const {
  return this->PointTypeArrayName;
}

//----------------------------------------------------------------------
void vtkSRepExportPolyDataProperties::SetLineTypeArrayName(const std::string& name) {
  if (this->LineTypeArrayName != name) {
    this->LineTypeArrayName = name;
    this->Modified();
  }
}

//----------------------------------------------------------------------
std::string vtkSRepExportPolyDataProperties::GetLineTypeArrayName() const {
  return this->LineTypeArrayName;
}

//----------------------------------------------------------------------
#define SREP_EXPORT_POLY_DATA_PROPERTIES_GET_SET(name) \
  void vtkSRepExportPolyDataProperties::Set##name(bool include) { \
    if (this->name != include) { \
      this->name = include; \
      this->Modified(); \
    } \
  } \
  bool vtkSRepExportPolyDataProperties::Get##name() const { \
    return this->name; \
  }

SREP_EXPORT_POLY_DATA_PROPERTIES_GET_SET(IncludeUpSpokes)
SREP_EXPORT_POLY_DATA_PROPERTIES_GET_SET(IncludeDownSpokes)
SREP_EXPORT_POLY_DATA_PROPERTIES_GET_SET(IncludeCrestSpokes)
SREP_EXPORT_POLY_DATA_PROPERTIES_GET_SET(IncludeCrestCurve)
SREP_EXPORT_POLY_DATA_PROPERTIES_GET_SET(IncludeSkeletalSheet)
SREP_EXPORT_POLY_DATA_PROPERTIES_GET_SET(IncludeSkeletonToCrestConnection)
SREP_EXPORT_POLY_DATA_PROPERTIES_GET_SET(IncludeSpine)
