#include "vtkSlicerSRepWidgetRepresentation.h"
#include "vtkSlicerSRepLogic.h"

#include <vtkActor.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkGlyph3D.h>
#include <vtkNamedColors.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyLine.h>
#include <vtkProperty.h>
#include <vtkSphereSource.h>
#include <vtkTubeFilter.h>
#include <vtkType.h>
#include <vtkUnsignedCharArray.h>

#include <vtkMRMLFolderDisplayNode.h>

#include <algorithm>
#include <functional>
#include <iterator>

vtkStandardNewMacro(vtkSlicerSRepWidgetRepresentation);

vtkSlicerSRepWidgetRepresentation::PointsRep::PointsRep()
  : GlyphSourceSphere(vtkSmartPointer<vtkSphereSource>::New())
  , Glypher(vtkSmartPointer<vtkGlyph3D>::New())
  , PointsPolyData(vtkSmartPointer<vtkPolyData>::New())
  , Property(vtkSmartPointer<vtkProperty>::New())
  , Mapper(vtkSmartPointer<vtkPolyDataMapper>::New())
  , Actor(vtkSmartPointer<vtkActor>::New())
  , TubeFilter(vtkSmartPointer<vtkTubeFilter>::New())
  , TubeMapper(vtkSmartPointer<vtkPolyDataMapper>::New())
  , TubeActor(vtkSmartPointer<vtkActor>::New())
{
  this->GlyphSourceSphere->SetRadius(0.5);

  this->Glypher->SetInputData(this->PointsPolyData);
  this->Glypher->ScalingOff();
  this->Glypher->SetScaleModeToDataScalingOff();
  this->Glypher->SetColorModeToColorByScalar();
  this->Glypher->SetSourceConnection(this->GlyphSourceSphere->GetOutputPort());

  this->Property = vtkSmartPointer<vtkProperty>::New();
  this->Property->SetRepresentationToSurface();
  this->Property->SetAmbient(0.0);
  this->Property->SetDiffuse(1.0);
  this->Property->SetSpecular(0.0);
  this->Property->SetShading(true);
  this->Property->SetSpecularPower(1.0);
  this->Property->SetPointSize(3.0);
  this->Property->SetLineWidth(3.0);
  this->Property->SetOpacity(1.);

  // This turns on resolve coincident topology for everything
  // as it is a class static on the mapper
  vtkMapper::SetResolveCoincidentTopologyToPolygonOffset();
  this->Mapper->SetInputConnection(this->Glypher->GetOutputPort());
  this->Mapper->ScalarVisibilityOn();
  this->Mapper->SetScalarModeToUsePointData();
  this->Mapper->SetColorModeToDirectScalars();

  this->Actor->SetMapper(this->Mapper);
  this->Actor->SetProperty(this->Property);

  this->TubeFilter->SetInputData(this->PointsPolyData);
  this->TubeFilter->SetNumberOfSides(10);
  this->TubeFilter->CappingOff();
  this->TubeFilter->SetVaryRadiusToVaryRadiusOff();

  this->TubeMapper->SetInputConnection(this->TubeFilter->GetOutputPort());
  this->TubeMapper->ScalarVisibilityOn();
  this->TubeMapper->SetScalarModeToUseCellData();
  this->TubeMapper->SetColorModeToDirectScalars();

  this->TubeActor->SetMapper(this->TubeMapper);
  this->TubeActor->SetProperty(this->Property);
}

vtkSlicerSRepWidgetRepresentation::PointsRep::~PointsRep() = default;

void vtkSlicerSRepWidgetRepresentation::PointsRep::SetPolyData(vtkSmartPointer<vtkPolyData> polyData) {
  this->PointsPolyData = polyData;
  this->Glypher->SetInputData(this->PointsPolyData);
  this->TubeFilter->SetInputData(this->PointsPolyData);
}

vtkSlicerSRepWidgetRepresentation::vtkSlicerSRepWidgetRepresentation()
  : Skeleton()
  , SRepDisplayNode(nullptr)
{}

vtkSlicerSRepWidgetRepresentation::~vtkSlicerSRepWidgetRepresentation() = default;

void vtkSlicerSRepWidgetRepresentation::PrintSelf(ostream& os, vtkIndent indent) {
  this->Superclass::PrintSelf(os, indent);
}

void vtkSlicerSRepWidgetRepresentation::GetActors(vtkPropCollection* pc) {
  this->Skeleton.Actor->GetActors(pc);
  this->Skeleton.TubeActor->GetActors(pc);
}
void vtkSlicerSRepWidgetRepresentation::ReleaseGraphicsResources(vtkWindow* window) {
  this->Skeleton.Actor->ReleaseGraphicsResources(window);
  this->Skeleton.TubeActor->ReleaseGraphicsResources(window);
}
int vtkSlicerSRepWidgetRepresentation::RenderOverlay(vtkViewport* viewport) {
  int count = 0;
  if (this->Skeleton.Actor->GetVisibility()) {
    count += this->Skeleton.Actor->RenderOverlay(viewport);
  }
  if (this->Skeleton.TubeActor->GetVisibility()) {
    count += this->Skeleton.TubeActor->RenderOverlay(viewport);
  }
  return count;
}
int vtkSlicerSRepWidgetRepresentation::RenderOpaqueGeometry(vtkViewport* viewport) {
  int count = 0;
  if (this->Skeleton.Actor->GetVisibility()) {
    count += this->Skeleton.Actor->RenderOpaqueGeometry(viewport);
  }
  if (this->Skeleton.TubeActor->GetVisibility()) {
    count += this->Skeleton.TubeActor->RenderOpaqueGeometry(viewport);
  }
  return count;
}
int vtkSlicerSRepWidgetRepresentation::RenderTranslucentPolygonalGeometry(vtkViewport* viewport) {
  int count = 0;

  // The internal actor needs to share property keys.
  // This ensures the mapper state is consistent and allows depth peeling to work as expected.
  this->Skeleton.Actor->SetPropertyKeys(this->GetPropertyKeys());
  this->Skeleton.TubeActor->SetPropertyKeys(this->GetPropertyKeys());

  if (this->Skeleton.Actor->GetVisibility()) {
    count += this->Skeleton.Actor->RenderTranslucentPolygonalGeometry(viewport);
  }
  if (this->Skeleton.TubeActor->GetVisibility()) {
    count += this->Skeleton.TubeActor->RenderTranslucentPolygonalGeometry(viewport);
  }
  return count;
}
vtkTypeBool vtkSlicerSRepWidgetRepresentation::HasTranslucentPolygonalGeometry() {
  if ((this->Skeleton.Actor->GetVisibility() && this->Skeleton.Actor->HasTranslucentPolygonalGeometry())
    || (this->Skeleton.TubeActor->GetVisibility() && this->Skeleton.TubeActor->HasTranslucentPolygonalGeometry())
  )
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
  const auto srep = srepNode->GetSRepWorld();
  if (!srep || srep->IsEmpty()) {
    this->VisibilityOff();
    return;
  }

  auto displayNode = this->GetSRepDisplayNode();
  if (!displayNode || !displayNode->GetVisibility()) {
    this->VisibilityOff();
    return;
  }

  this->VisibilityOn();

  auto logic = vtkSmartPointer<vtkSlicerSRepLogic>::New();
  this->Skeleton.SetPolyData(logic->SmartExportSRepToPolyData(*srep, *displayNode->SmartGetSRepExportPolyDataProperties()));

  // set point size
  double radius = 0;
  if (displayNode->GetUseAbsoluteThickness()) {
    radius = displayNode->GetAbsoluteThickness();
  } else {
    this->Skeleton.PointsPolyData->ComputeBounds();
    double bounds[6];
    this->Skeleton.PointsPolyData->GetBounds(bounds);

    const double minPoint[] = {bounds[0], bounds[2], bounds[4]};
    const double maxPoint[] = {bounds[1], bounds[3], bounds[5]};
    const double dist = sqrt(vtkMath::Distance2BetweenPoints(minPoint, maxPoint));
    radius = dist * displayNode->GetRelativeThickness();
  }

  this->Skeleton.GlyphSourceSphere->SetRadius(radius);
  this->Skeleton.TubeFilter->SetRadius(radius);

  this->Skeleton.Property->SetOpacity(displayNode->GetOpacity());
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
