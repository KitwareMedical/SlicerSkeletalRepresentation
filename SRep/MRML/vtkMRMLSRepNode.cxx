#include "vtkMRMLSRepNode.h"
#include "vtkMRMLSRepStorageNode.h"
#include <vtkAbstractTransform.h>
#include <vtkGeneralTransform.h>
#include <vtkMRMLScene.h>
#include <vtkMRMLTransformNode.h>
#include <vtkBoundingBox.h>

#include <srep/Spoke.h>
#include <srep/SkeletalPoint.h>
#include <srep/RectangularGridSRep.h>

//----------------------------------------------------------------------------
std::vector<std::vector<srep::SkeletalPoint>>
TransformSkeletalPoints(const std::vector<std::vector<srep::SkeletalPoint>>& grid, vtkAbstractTransform* transform) {
  if (!transform) {
    return grid;
  }

  using namespace srep;
  std::vector<std::vector<srep::SkeletalPoint>> transformedGrid;
  transformedGrid.reserve(grid.size());
  for (const auto& vectorOfPoints : grid) {
    transformedGrid.resize(transformedGrid.size() + 1);
    transformedGrid.back().reserve(vectorOfPoints.size());
    for (const auto& skeletalPoint : vectorOfPoints) {
      const auto& upSpoke = skeletalPoint.GetUpSpoke();
      const auto& downSpoke = skeletalPoint.GetDownSpoke();
      double transformedSkeletal[3];
      double transformedBoundary[3];

      transform->TransformPoint(upSpoke.GetSkeletalPoint().AsArray().data(), transformedSkeletal);
      transform->TransformPoint(upSpoke.GetBoundaryPoint().AsArray().data(), transformedBoundary);
      const Spoke transformedUpSpoke(Point3d{transformedSkeletal}, Vector3d{Point3d{transformedSkeletal}, Point3d{transformedBoundary}});

      transform->TransformPoint(downSpoke.GetSkeletalPoint().AsArray().data(), transformedSkeletal);
      transform->TransformPoint(downSpoke.GetBoundaryPoint().AsArray().data(), transformedBoundary);
      const Spoke transformedDownSpoke(Point3d{transformedSkeletal}, Vector3d{Point3d{transformedSkeletal}, Point3d{transformedBoundary}});
      if (skeletalPoint.IsCrest()) {
        const auto& crestSpoke = skeletalPoint.GetCrestSpoke();
        transform->TransformPoint(crestSpoke.GetSkeletalPoint().AsArray().data(), transformedSkeletal);
        transform->TransformPoint(crestSpoke.GetBoundaryPoint().AsArray().data(), transformedBoundary);
        const Spoke transformedCrestSpoke(Point3d{transformedSkeletal}, Vector3d{Point3d{transformedSkeletal}, Point3d{transformedBoundary}});

        transformedGrid.back().emplace_back(transformedUpSpoke, transformedDownSpoke, transformedCrestSpoke);
      } else {
        transformedGrid.back().emplace_back(transformedUpSpoke, transformedDownSpoke);
      }
    }
  }
  return transformedGrid;
}

//----------------------------------------------------------------------------
vtkMRMLSRepNode::vtkMRMLSRepNode()
  : vtkMRMLDisplayableNode()
  , SRepTransform()
{
  this->SRepTransform->Identity();
}

//----------------------------------------------------------------------------
vtkMRMLSRepNode::~vtkMRMLSRepNode() = default;

//----------------------------------------------------------------------------
bool vtkMRMLSRepNode::HasSRep() const {
  return static_cast<bool>(this->GetSRep());
}

//----------------------------------------------------------------------------
void vtkMRMLSRepNode::PrintSelf(ostream& os, vtkIndent indent) {
  Superclass::PrintSelf(os,indent);

  os << "SRep";
}

//----------------------------------------------------------------------------
void vtkMRMLSRepNode::CreateDefaultDisplayNodes() {
  if (this->GetDisplayNode()
    && vtkMRMLSRepDisplayNode::SafeDownCast(this->GetDisplayNode()))
  {
    //display node already exists
    return;
  }
  if (!this->GetScene()) {
    vtkErrorMacro("vtkMRMLSRepNode::CreateDefaultDisplayNodes failed: scene is invalid");
    return;
  }
  vtkMRMLSRepDisplayNode* dispNode = vtkMRMLSRepDisplayNode::SafeDownCast(
    this->GetScene()->AddNewNodeByClass("vtkMRMLSRepDisplayNode"));

  if (!dispNode) {
    vtkErrorMacro("vtkMRMLSRepNode::CreateDefaultDisplayNodes failed: scene failed to instantiate a vtkMRMLSRepDisplayNode node");
    return;
  }
  this->SetAndObserveDisplayNodeID(dispNode->GetID());
}

//----------------------------------------------------------------------------
vtkMRMLStorageNode* vtkMRMLSRepNode::CreateDefaultStorageNode()
{
  vtkMRMLScene* scene = this->GetScene();
  if (scene == nullptr)
    {
    vtkErrorMacro("CreateDefaultStorageNode failed: scene is invalid");
    return nullptr;
    }
  return vtkMRMLStorageNode::SafeDownCast(scene->CreateNodeByClass("vtkMRMLSRepStorageNode"));
}

//----------------------------------------------------------------------------
vtkMRMLSRepDisplayNode* vtkMRMLSRepNode::GetSRepDisplayNode() {
  auto* dispNode = this->GetDisplayNode();
  if (dispNode) {
    return vtkMRMLSRepDisplayNode::SafeDownCast(dispNode);
  }
  return nullptr;
}

//----------------------------------------------------------------------------
void vtkMRMLSRepNode::GetRASBounds(double bounds[6]) {
  vtkMRMLSRepNode::GetSRepBounds(this->GetSRepWorld(), bounds);
}

//----------------------------------------------------------------------------
void vtkMRMLSRepNode::GetBounds(double bounds[6]) {
  vtkMRMLSRepNode::GetSRepBounds(this->GetSRep(), bounds);
}

//----------------------------------------------------------------------------
void vtkMRMLSRepNode::GetSRepBounds(const srep::MeshSRepInterface* srep, double bounds[6]) {
  vtkBoundingBox box;

  if (!srep) {
    box.GetBounds(bounds);
    return;
  }

  const auto addSpokeMesh = [&](const srep::SpokeMesh& mesh) {
    for (long i = 0 ; i < mesh.GetNumberOfSpokes(); ++i) {
      box.AddPoint(mesh[i].GetBoundaryPoint().AsArray().data());
    }
  };

  addSpokeMesh(srep->GetUpSpokes());
  addSpokeMesh(srep->GetDownSpokes());
  addSpokeMesh(srep->GetCrestSpokes());

  box.GetBounds(bounds);
}

//---------------------------------------------------------------------------
bool vtkMRMLSRepNode::CanApplyNonLinearTransforms() const
{
  return true;
}

//----------------------------------------------------------------------------
void vtkMRMLSRepNode::CopyContent(vtkMRMLNode* anode, bool deepCopy/*=true*/) {
  Superclass::CopyContent(anode, deepCopy);
}

//---------------------------------------------------------------------------
void vtkMRMLSRepNode::OnTransformNodeReferenceChanged(vtkMRMLTransformNode* transformNode) {

  //this next line is a GetTransformToWorld one-liner that works even if this->GetParentTransformNode is nullptr
  vtkMRMLTransformNode::GetTransformBetweenNodes(this->GetParentTransformNode(), nullptr, this->SRepTransform);
  this->UpdateSRepWorld();
  Superclass::OnTransformNodeReferenceChanged(transformNode);
  this->Modified();
}

//---------------------------------------------------------------------------
void vtkMRMLSRepNode::ProcessMRMLEvents (vtkObject* caller,
                                         unsigned long event,
                                         void* callData)
{
  if (caller != nullptr && event == vtkMRMLTransformableNode::TransformModifiedEvent) {
    //this next line is a GetTransformToWorld one-liner that works even if this->GetParentTransformNode is nullptr
    vtkMRMLTransformNode::GetTransformBetweenNodes(this->GetParentTransformNode(), nullptr, this->SRepTransform);
    this->UpdateSRepWorld();
  }

  Superclass::ProcessMRMLEvents(caller, event, callData);
}

//---------------------------------------------------------------------------
void vtkMRMLSRepNode::UpdateSRepWorld() {
  this->DoUpdateSRepWorld(this->SRepTransform);
}
