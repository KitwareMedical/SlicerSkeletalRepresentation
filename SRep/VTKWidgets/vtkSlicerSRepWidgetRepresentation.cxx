#include "vtkSlicerSRepWidgetRepresentation.h"

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

vtkStandardNewMacro(vtkSlicerSRepWidgetRepresentation);

vtkSlicerSRepWidgetRepresentation::PointsRep::PointsRep()
  : GlyphSourceSphere(vtkSmartPointer<vtkSphereSource>::New())
  , Glypher(vtkSmartPointer<vtkGlyph3D>::New())
  , Points(vtkSmartPointer<vtkPoints>::New())
  , PointColors(vtkSmartPointer<vtkUnsignedCharArray>::New())
  , Lines(vtkSmartPointer<vtkCellArray>::New())
  , LineColors(vtkSmartPointer<vtkUnsignedCharArray>::New())
  , PointsPolyData(vtkSmartPointer<vtkPolyData>::New())
  , Property(vtkSmartPointer<vtkProperty>::New())
  , Mapper(vtkSmartPointer<vtkPolyDataMapper>::New())
  , Actor(vtkSmartPointer<vtkActor>::New())
  , TubeFilter(vtkSmartPointer<vtkTubeFilter>::New())
  , TubeMapper(vtkSmartPointer<vtkPolyDataMapper>::New())
  , TubeActor(vtkSmartPointer<vtkActor>::New())
{
  this->GlyphSourceSphere->SetRadius(0.5);

  this->PointsPolyData->SetPoints(this->Points);
  this->PointsPolyData->GetPointData()->SetScalars(this->PointColors);
  this->PointsPolyData->SetLines(this->Lines);
  this->PointsPolyData->GetCellData()->SetScalars(this->LineColors);

  this->PointColors->SetNumberOfComponents(3);
  this->LineColors->SetNumberOfComponents(3);

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
  const auto srep = srepNode->GetSRep();
  if (!srep || srep->IsEmpty()) {
    this->VisibilityOff();
    return;
  }

  this->VisibilityOn();

  vtkNew<vtkNamedColors> namedColors;

  this->Skeleton.Points->Reset();
  this->Skeleton.PointColors->Reset();
  this->Skeleton.Lines->Reset();
  this->Skeleton.LineColors->Reset();

  const auto& grid = srep->GetSkeletalPoints();
  std::vector<std::vector<vtkIdType>> gridToPointId;
  gridToPointId.resize(grid.size());
  for (size_t row = 0; row < grid.size(); ++row) {
    gridToPointId[row].resize(grid[row].size());
  }

  const auto insertNextLine = [this](const vtkIdType start, const vtkIdType end, const vtkColor3ub& color) {
    this->Skeleton.Lines->InsertNextCell(2);
    this->Skeleton.Lines->InsertCellPoint(start);
    this->Skeleton.Lines->InsertCellPoint(end);
    this->Skeleton.LineColors->InsertNextTypedTuple(color.GetData());
  };

  const auto insertNextPoint = [this](const srep::Point3d& point, const vtkColor3ub& color) {
    const auto id = this->Skeleton.Points->InsertNextPoint(point.AsArray().data());
    this->Skeleton.PointColors->InsertNextTypedTuple(color.GetData());
    return id;
  };

  const auto upColor = namedColors->GetColor3ub("Tomato");
  const auto downColor = namedColors->GetColor3ub("Mint");
  const auto skeletonColor = namedColors->GetColor3ub("Cornsilk");
  const auto crestColor = namedColors->GetColor3ub("Gold");
  const auto crestToSkeletonConnectionColor = namedColors->GetColor3ub("Black");

  // visualize all points and up/down spoke lines
  for (size_t row = 0; row < grid.size(); ++row) {
    for (size_t col = 0; col < grid[row].size(); ++col) {
      const auto& skeletalPoint = grid[row][col];
      const auto upSkeletalPointId = insertNextPoint(skeletalPoint.GetUpSpoke().GetSkeletalPoint(), skeletonColor);
      //don't really know where to attach the fold if the up and down aren't at the same skeletal point, so just attach to up?
      gridToPointId[row][col] = upSkeletalPointId;

      const auto upBoundaryId = insertNextPoint(skeletalPoint.GetUpSpoke().GetBoundaryPoint(), upColor);
      insertNextLine(upSkeletalPointId, upBoundaryId, upColor); // up spoke

      const auto downSkeletalPointId = insertNextPoint(skeletalPoint.GetDownSpoke().GetSkeletalPoint(), skeletonColor);
      const auto downBoundaryId = insertNextPoint(skeletalPoint.GetDownSpoke().GetBoundaryPoint(), downColor);
      insertNextLine(downSkeletalPointId, downBoundaryId, downColor); // down spoke
    }
  }

  // visualize crest (aka fold)
  // loop through crest clockwise manner for easier connection of the fold
  std::vector<vtkIdType> crestBaseIds;

  const auto insertCrest = [&](size_t row, size_t col) {
    const auto& crestSpoke = grid[row][col].GetCrestSpoke();

    const auto crestBoundaryId = insertNextPoint(crestSpoke.GetBoundaryPoint(), crestColor);
    const auto crestBaseId = insertNextPoint(crestSpoke.GetSkeletalPoint(), crestColor);
    crestBaseIds.push_back(crestBaseId);

    insertNextLine(crestBoundaryId, crestBaseId, crestColor); // crest spoke
    insertNextLine(gridToPointId[row][col], crestBaseId, crestToSkeletonConnectionColor); // connection skeletal grid to fold
  };
  //do top row left to right
  for (size_t col = 0; col < grid[0].size(); ++col) {
    insertCrest(0, col);
  }
  //do right col top to bottom
  for (size_t row = 0; row < grid.size(); ++row) {
    insertCrest(row, grid[row].size() - 1);
  }
  //do bottom row right to left
  for (size_t col = grid.back().size() - 1; ; --col) {
    insertCrest(grid.size() - 1, col);
    if (col == 0) {
      break;
    }
  }
  //do left col bottom to top
  for (size_t row = grid.size() - 1; ; --row) {
    insertCrest(row, 0);
    if (row == 0) {
      break;
    }
  }

  //connect the fold
  for (size_t cur = 0; cur < crestBaseIds.size(); ++cur) {
    const size_t next = (cur + 1) % crestBaseIds.size();
    insertNextLine(crestBaseIds[cur], crestBaseIds[next], crestColor);
  }

  //connect all skeletal points
  for (size_t row = 0; row < grid.size(); ++row) {
    for (size_t col = 0; col < grid[row].size(); ++col) {
      //left
      if (col != 0) {
        insertNextLine(gridToPointId[row][col], gridToPointId[row][col-1], skeletonColor);
      }

      //right
      if (col != grid[row].size() - 1) {
        insertNextLine(gridToPointId[row][col], gridToPointId[row][col+1], skeletonColor);
      }

      //up
      if (row != 0) {
        insertNextLine(gridToPointId[row][col], gridToPointId[row-1][col], skeletonColor);
      }

      //down
      if (row != grid.size() - 1) {
        insertNextLine(gridToPointId[row][col], gridToPointId[row+1][col], skeletonColor);
      }
    }
  }

  // set point size
  this->Skeleton.Points->ComputeBounds();
  const auto bounds = this->Skeleton.Points->GetBounds();

  const double minPoint[] = {bounds[0], bounds[2], bounds[4]};
  const double maxPoint[] = {bounds[1], bounds[3], bounds[5]};
  const double distSquared = vtkMath::Distance2BetweenPoints(minPoint, maxPoint);

  constexpr double divisor = 500;
  const double radius = distSquared / divisor;
  this->Skeleton.GlyphSourceSphere->SetRadius(radius);
  this->Skeleton.TubeFilter->SetRadius(radius);

  auto displayNode = this->GetSRepDisplayNode();
  if (displayNode) {
    this->Skeleton.Property->SetOpacity(displayNode->GetOpacity());
  }
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
