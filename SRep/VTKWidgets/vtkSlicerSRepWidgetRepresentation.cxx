#include "vtkSlicerSRepWidgetRepresentation.h"

#include <vtkActor.h>
#include <vtkDoubleArray.h>
#include <vtkGlyph3D.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkSphereSource.h>

#include <vtkMRMLFolderDisplayNode.h>

vtkStandardNewMacro(vtkSlicerSRepWidgetRepresentation);

vtkSlicerSRepWidgetRepresentation::SpokesRep::SpokesRep()
  : GlyphSourceSphere(vtkSmartPointer<vtkSphereSource>::New())
  , Glypher(vtkSmartPointer<vtkGlyph3D>::New())
  , BoundaryPoints(vtkSmartPointer<vtkPoints>::New())
  , BoundaryPointsPolyData(vtkSmartPointer<vtkPolyData>::New())
  , Property(vtkSmartPointer<vtkProperty>::New())
  , Mapper(vtkSmartPointer<vtkPolyDataMapper>::New())
  , Actor(vtkSmartPointer<vtkActor>::New())
{
  this->GlyphSourceSphere->SetRadius(0.5);

  vtkNew<vtkDoubleArray> boundaryPointNormals;
  boundaryPointNormals->SetNumberOfComponents(3);

  this->BoundaryPointsPolyData->SetPoints(this->BoundaryPoints);
  this->BoundaryPointsPolyData->GetPointData()->SetNormals(boundaryPointNormals);

  this->Glypher->SetInputData(this->BoundaryPointsPolyData);
  this->Glypher->OrientOn();
  this->Glypher->ScalingOn();
  this->Glypher->SetScaleModeToDataScalingOff();
  this->Glypher->SetScaleFactor(1.0);
  this->Glypher->SetSourceConnection(this->GlyphSourceSphere->GetOutputPort());

  this->Property = vtkSmartPointer<vtkProperty>::New();
  this->Property->SetRepresentationToSurface();
  this->Property->SetAmbient(0.0);
  this->Property->SetDiffuse(1.0);
  this->Property->SetSpecular(0.0);
  this->Property->SetShading(true);
  this->Property->SetSpecularPower(1.0);
  this->Property->SetPointSize(3.);
  this->Property->SetLineWidth(3.);
  this->Property->SetOpacity(1.);

  this->Mapper->SetInputConnection(this->Glypher->GetOutputPort());
  // This turns on resolve coincident topology for everything
  // as it is a class static on the mapper
  vtkMapper::SetResolveCoincidentTopologyToPolygonOffset();
  this->Mapper->ScalarVisibilityOff();

  this->Actor->SetMapper(this->Mapper);
  this->Actor->SetProperty(this->Property);
}

vtkSlicerSRepWidgetRepresentation::SpokesRep::~SpokesRep() = default;

vtkSlicerSRepWidgetRepresentation::vtkSlicerSRepWidgetRepresentation()
  : UpSpokes()
  , DownSpokes()
{
  this->UpSpokes.Property->SetColor(0.4, 1.0, 1.0);
  this->DownSpokes.Property->SetColor(1.0, 0.4, 1.0);
}

vtkSlicerSRepWidgetRepresentation::~vtkSlicerSRepWidgetRepresentation() = default;

void vtkSlicerSRepWidgetRepresentation::PrintSelf(ostream& os, vtkIndent indent) {
  this->Superclass::PrintSelf(os, indent);
}

void vtkSlicerSRepWidgetRepresentation::GetActors(vtkPropCollection* pc) {
  this->UpSpokes.Actor->GetActors(pc);
  this->DownSpokes.Actor->GetActors(pc);
}
void vtkSlicerSRepWidgetRepresentation::ReleaseGraphicsResources(vtkWindow* window) {
  this->UpSpokes.Actor->ReleaseGraphicsResources(window);
  this->DownSpokes.Actor->ReleaseGraphicsResources(window);
}
int vtkSlicerSRepWidgetRepresentation::RenderOverlay(vtkViewport* viewport) {
  int count = 0;
  if (this->UpSpokes.Actor->GetVisibility()) {
    count += this->UpSpokes.Actor->RenderOverlay(viewport);
  }
  if (this->DownSpokes.Actor->GetVisibility()) {
    count += this->DownSpokes.Actor->RenderOverlay(viewport);
  }
  return count;
}
int vtkSlicerSRepWidgetRepresentation::RenderOpaqueGeometry(vtkViewport* viewport) {
  int count = 0;
  if (this->UpSpokes.Actor->GetVisibility()) {
    count += this->UpSpokes.Actor->RenderOpaqueGeometry(viewport);
  }
  if (this->DownSpokes.Actor->GetVisibility()) {
    count += this->DownSpokes.Actor->RenderOpaqueGeometry(viewport);
  }
  return count;
}
int vtkSlicerSRepWidgetRepresentation::RenderTranslucentPolygonalGeometry(vtkViewport* viewport) {
  int count = 0;
  if (this->UpSpokes.Actor->GetVisibility()) {
    count += this->UpSpokes.Actor->RenderTranslucentPolygonalGeometry(viewport);
  }
  if (this->DownSpokes.Actor->GetVisibility()) {
    count += this->DownSpokes.Actor->RenderTranslucentPolygonalGeometry(viewport);
  }
  return count;
}
vtkTypeBool vtkSlicerSRepWidgetRepresentation::HasTranslucentPolygonalGeometry() {
  if ((this->UpSpokes.Actor->GetVisibility() && this->UpSpokes.Actor->HasTranslucentPolygonalGeometry())
    || (this->DownSpokes.Actor->GetVisibility() && this->DownSpokes.Actor->HasTranslucentPolygonalGeometry()))
  {
    return true;
  }
  return false;
}

void vtkSlicerSRepWidgetRepresentation::UpdateFromMRML(vtkMRMLNode* caller, unsigned long event, void *callData) {
  Superclass::UpdateFromMRML(caller, event, callData);

  //TODO: does this always need to be on?
  this->NeedToRenderOn();
  
  vtkMRMLSRepNode* srepNode = this->GetSRepNode();

  if (!srepNode || !this->IsDisplayable()) {
    this->VisibilityOff();
    return;
  }

  this->VisibilityOn();

  srepNode->PopulateUpBoundaryPoints(this->UpSpokes.BoundaryPoints);
  srepNode->PopulateDownBoundaryPoints(this->DownSpokes.BoundaryPoints);

  this->UpSpokes.BoundaryPoints->ComputeBounds();
  auto upBounds = this->UpSpokes.BoundaryPoints->GetBounds();
  this->DownSpokes.BoundaryPoints->ComputeBounds();
  auto downBounds = this->DownSpokes.BoundaryPoints->GetBounds();
  const double minPoint[] = {
    std::min(upBounds[0], downBounds[0]),
    std::min(upBounds[2], downBounds[2]),
    std::min(upBounds[4], downBounds[4])
  };
  const double maxPoint[] = {
    std::max(upBounds[1], downBounds[1]),
    std::max(upBounds[3], downBounds[3]),
    std::max(upBounds[5], downBounds[5])
  };
  const double distSquared = vtkMath::Distance2BetweenPoints(minPoint, maxPoint);

  this->UpSpokes.GlyphSourceSphere->SetRadius(distSquared / 100.);
  this->DownSpokes.GlyphSourceSphere->SetRadius(distSquared / 100.);
}

void vtkSlicerSRepWidgetRepresentation::SetSRepDisplayNode(vtkMRMLSRepDisplayNode* srepDisplayNode) {
  this->SRepDisplayNode = srepDisplayNode;
}
vtkMRMLSRepDisplayNode* vtkSlicerSRepWidgetRepresentation::GetSRepDisplayNode() const {
  return this->SRepDisplayNode;
}
vtkMRMLSRepNode* vtkSlicerSRepWidgetRepresentation::GetSRepNode() {
  if (this->SRepDisplayNode) {
    return this->SRepDisplayNode->GetSRepNode();
  }
  return nullptr;
}

bool vtkSlicerSRepWidgetRepresentation::IsDisplayable()
{
  if (!this->SRepDisplayNode
    || !this->ViewNode
    || !this->SRepDisplayNode->GetVisibility()
    || !this->SRepDisplayNode->IsDisplayableInView(this->ViewNode->GetID()))
  {
    return false;
  }

  // If parent folder visibility is set to false then the srep is not visible
  if (this->SRepDisplayNode->GetFolderDisplayOverrideAllowed()) {
    vtkMRMLDisplayableNode* displayableNode = this->SRepDisplayNode->GetDisplayableNode();
    // Visibility is applied regardless the fact whether there is override or not.
    // Visibility of items defined by hierarchy is off if any of the ancestors is explicitly hidden.
    // However, this does not apply on display nodes that do not allow overrides (FolderDisplayOverrideAllowed)
    if (!vtkMRMLFolderDisplayNode::GetHierarchyVisibility(displayableNode)) {
      return false;
    }
  }
  return this->SRepDisplayNode->GetVisibility3D();
}
