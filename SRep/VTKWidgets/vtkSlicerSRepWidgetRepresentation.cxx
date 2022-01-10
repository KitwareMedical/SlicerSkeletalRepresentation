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

#include <algorithm>
#include <functional>
#include <iterator>

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

namespace {
struct vtkSpokeIds {
  vtkIdType boundaryId;
  vtkIdType skeletonId;
};
}

void vtkSlicerSRepWidgetRepresentation::ConvertSRepToVisualRepresentation(const vtkMeshSRepInterface& srep, const vtkMRMLSRepDisplayNode& displayNode) {
  //-------------------------------
  const auto insertNextPoint = [this](const srep::Point3d& point, const vtkColor3ub& color) {
    const auto id = this->Skeleton.Points->InsertNextPoint(point.AsArray().data());
    this->Skeleton.PointColors->InsertNextTypedTuple(color.GetData());
    return id;
  };

  //-------------------------------
  const auto insertNextLine = [this](const vtkIdType start, const vtkIdType end, const vtkColor3ub& color) {
    this->Skeleton.Lines->InsertNextCell(2);
    this->Skeleton.Lines->InsertCellPoint(start);
    this->Skeleton.Lines->InsertCellPoint(end);
    this->Skeleton.LineColors->InsertNextTypedTuple(color.GetData());
  };

  //-------------------------------
  const auto addSpokeMesh = [this, insertNextPoint, insertNextLine]
    (const vtkSRepSpokeMesh& mesh,
     bool addSpokes,
     const vtkColor3ub& spokeColor,
     bool addConnections,
     const vtkColor3ub& connectionColor,
     const std::vector<vtkMeshSRepInterface::IndexType> spine,
     const vtkColor3ub& spineColor,
     std::function<bool(long i)> forceAddSkeletalPoint) -> std::vector<vtkSpokeIds>{
    std::vector<vtkSpokeIds> spokesToVTKPointIds;

    const auto isSpine = [&spine](long i) {
      return std::find(spine.begin(), spine.end(), i) != spine.end();
    };

    // add all the points and the spoke lines
    for (long i = 0; i < mesh.GetNumberOfSpokes(); ++i) {
      vtkSpokeIds ids;
      if (addSpokes || addConnections || isSpine(i) || forceAddSkeletalPoint(i)) {
        ids.skeletonId = insertNextPoint(mesh[i]->GetSkeletalPoint(), spokeColor);
      } else {
        ids.skeletonId = -1;
      }
      if (addSpokes) {
        ids.boundaryId = insertNextPoint(mesh[i]->GetBoundaryPoint(), spokeColor);
        insertNextLine(ids.skeletonId, ids.boundaryId, spokeColor);
      } else {
        ids.boundaryId = -1;
      }
      spokesToVTKPointIds.push_back(ids);
    }

    // add the connection lines. It is essentially a bidirectional graph, so only add one line between two points, even if
    // it shows up twice "once in each direction"
    std::set<std::pair<vtkIdType, vtkIdType>> spineConnections;
    for (size_t i = 1; i < spine.size(); ++i) {
      //insert a sorted pair
      vtkIdType point1 = spokesToVTKPointIds[spine[i-1]].skeletonId;
      vtkIdType point2 = spokesToVTKPointIds[spine[i]].skeletonId;
      //sort the points
      if (point1 > point2) {
        std::swap(point1, point2);
      }
      spineConnections.insert(std::make_pair(point1, point2));
    }

    for (const auto& spineConnection : spineConnections) {
      insertNextLine(spineConnection.first, spineConnection.second, spineColor);
    }

    if (addConnections) {
      std::set<std::pair<vtkIdType, vtkIdType>> connections;
      for (long i = 0; i < mesh.GetNumberOfSpokes(); ++i) {
        const auto neighbors = mesh.GetNeighbors(i);
        for (size_t neighbor : neighbors) {
          vtkIdType point1 = spokesToVTKPointIds[i].skeletonId;
          vtkIdType point2 = spokesToVTKPointIds[neighbor].skeletonId;
          //sort the points
          if (point1 > point2) {
            std::swap(point1, point2);
          }
          connections.insert(std::make_pair(point1, point2));
        }
      }

      std::set<std::pair<vtkIdType, vtkIdType>> spinelessConnections;
      std::set_difference(
        connections.begin(), connections.end(),
        spineConnections.begin(), spineConnections.end(),
        std::inserter(spinelessConnections, spinelessConnections.begin()));

      for (const auto& connection : spinelessConnections) {
        insertNextLine(connection.first, connection.second, connectionColor);
      }
    }

    return spokesToVTKPointIds;
  };

  ///////////////////////////////////////
  // Start
  ///////////////////////////////////////

  const auto upSpokeColor = displayNode.GetUpSpokeColor();
  const auto downSpokeColor = displayNode.GetDownSpokeColor();
  const auto upSkeletonColor = displayNode.GetSkeletalSheetColor();
  const auto downSkeletonColor = displayNode.GetSkeletalSheetColor();
  const auto crestCurveColor = displayNode.GetCrestCurveColor();
  const auto crestSpokeColor = displayNode.GetCrestSpokeColor();
  const auto crestToSkeletonConnectionColor = displayNode.GetSkeletonToCrestConnectionColor();
  const auto spineColor = displayNode.GetSpineColor();

  this->Skeleton.Points->Reset();
  this->Skeleton.PointColors->Reset();
  this->Skeleton.Lines->Reset();
  this->Skeleton.LineColors->Reset();

  const auto visibleUpSpineIndexes = displayNode.GetSpineVisibility() ? srep.GetUpSpine() : vtkSRepSpokeMesh::NeighborList{};
  const auto visibleDownSpineIndexes = displayNode.GetSpineVisibility() ? srep.GetDownSpine() : vtkSRepSpokeMesh::NeighborList{};

  const auto upSpokeToPointIds = addSpokeMesh(
    *srep.GetUpSpokes(),
    displayNode.GetUpSpokeVisibility(), upSpokeColor,
    displayNode.GetSkeletalSheetVisibility(), upSkeletonColor,
    visibleUpSpineIndexes, spineColor,
    [&](long i) {
      const auto crestSkeletonConnections = srep.GetCrestToUpSpokeConnections();
      return displayNode.GetSkeletonToCrestConnectionVisibility() && crestSkeletonConnections.end() != std::find(crestSkeletonConnections.begin(), crestSkeletonConnections.end(), i);
      });

  const auto downSpokeToPointIds = addSpokeMesh(
    *srep.GetDownSpokes(),
    displayNode.GetDownSpokeVisibility(), downSpokeColor,
    displayNode.GetSkeletalSheetVisibility(), downSkeletonColor,
    visibleDownSpineIndexes, spineColor,
    [&](long i) {
      const auto crestSkeletonConnections = srep.GetCrestToDownSpokeConnections();
      return displayNode.GetSkeletonToCrestConnectionVisibility() && crestSkeletonConnections.end() != std::find(crestSkeletonConnections.begin(), crestSkeletonConnections.end(), i);
      });

  const auto crestSpokeToPointIds = addSpokeMesh(
    *srep.GetCrestSpokes(),
    displayNode.GetCrestSpokeVisibility(), crestSpokeColor,
    displayNode.GetCrestCurveVisibility(), crestCurveColor,
    {}, vtkColor3ub(), // no spine for crest
    [&](long){ return displayNode.GetSkeletonToCrestConnectionVisibility(); });

  // connect the crest to skeleton
  if (displayNode.GetSkeletonToCrestConnectionVisibility()) {
    for (size_t crestIndex = 0; crestIndex < srep.GetCrestToUpSpokeConnections().size(); ++crestIndex) {
      const auto skeletonIndex = srep.GetCrestToUpSpokeConnections()[crestIndex];
      insertNextLine(crestSpokeToPointIds[crestIndex].skeletonId, upSpokeToPointIds[skeletonIndex].skeletonId,
        crestToSkeletonConnectionColor);
    }
    for (size_t crestIndex = 0; crestIndex < srep.GetCrestToDownSpokeConnections().size(); ++crestIndex) {
      const auto skeletonIndex = srep.GetCrestToDownSpokeConnections()[crestIndex];
      insertNextLine(crestSpokeToPointIds[crestIndex].skeletonId, downSpokeToPointIds[skeletonIndex].skeletonId,
        crestToSkeletonConnectionColor);
    }
  }
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

  // TODO: if performance is an issue, save of a MTime or do something to only do the conversion
  // when the content of srepNode->GetSRepWorld() changes
  this->ConvertSRepToVisualRepresentation(*srep, *displayNode);

  // set point size
  double radius = 0;
  if (displayNode->GetUseAbsoluteThickness()) {
    radius = displayNode->GetAbsoluteThickness();
  } else {
    this->Skeleton.Points->ComputeBounds();
    double bounds[6];
    this->Skeleton.Points->GetBounds(bounds);

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
