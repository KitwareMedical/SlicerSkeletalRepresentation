#ifndef vtkSlicerSRepWidget_h
#define vtkSlicerSRepWidget_h

#include "vtkSlicerSRepModuleVTKWidgetsExport.h"
#include "vtkMRMLAbstractWidget.h"

#include "vtkMRMLSRepNode.h"

class VTK_SLICER_SREP_MODULE_VTKWIDGETS_EXPORT vtkSlicerSRepWidget
  : public vtkMRMLAbstractWidget
{
public:
  /// Instantiate this class.
  static vtkSlicerSRepWidget *New();
  /// Standard methods for a VTK class.
  vtkTypeMacro(vtkSlicerSRepWidget, vtkMRMLAbstractWidget);
  void PrintSelf(ostream& os, vtkIndent indent) override;

  void CreateDefaultRepresentation(vtkMRMLSRepDisplayNode* srepDisplayNode, vtkMRMLAbstractViewNode* viewNode, vtkRenderer* renderer);
};

#endif