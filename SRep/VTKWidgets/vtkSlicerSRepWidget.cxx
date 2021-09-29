#include "vtkSlicerSRepWidget.h"
#include "vtkSlicerSRepWidgetRepresentation.h"

vtkStandardNewMacro(vtkSlicerSRepWidget);

//----------------------------------------------------------------------
void vtkSlicerSRepWidget::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os, indent);
}

void vtkSlicerSRepWidget::CreateDefaultRepresentation(vtkMRMLSRepDisplayNode* srepDisplayNode, vtkMRMLAbstractViewNode* viewNode, vtkRenderer* renderer) {
  if (!srepDisplayNode) {
    return;
  }
  auto rep = vtkSmartPointer<vtkSlicerSRepWidgetRepresentation>::New();
  this->SetRenderer(renderer);
  this->SetRepresentation(rep);
  rep->SetViewNode(viewNode);
  rep->SetSRepDisplayNode(srepDisplayNode);
}