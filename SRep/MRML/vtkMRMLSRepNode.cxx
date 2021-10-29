#include "vtkMRMLSRepNode.h"
#include <vtkAbstractTransform.h>
#include <vtkGeneralTransform.h>
#include <vtkMRMLScene.h>
#include <vtkMRMLTransformNode.h>
#include <vtkBoundingBox.h>

#include <srep/SRepIO.h>
#include <srep/Spoke.h>
#include <srep/SkeletalPoint.h>
#include <srep/SRep.h>

//----------------------------------------------------------------------------
vtkMRMLNodeNewMacro(vtkMRMLSRepNode);

//----------------------------------------------------------------------------
vtkMRMLSRepNode::vtkMRMLSRepNode()
  : vtkMRMLDisplayableNode()
  , SRep()
  , SRepWorld()
  , SRepTransform()
{
  this->SRepTransform->Identity();
}

//----------------------------------------------------------------------------
vtkMRMLSRepNode::~vtkMRMLSRepNode() = default;

//----------------------------------------------------------------------------
void vtkMRMLSRepNode::LoadSRepFromFile(const std::string& filename) {
  this->SRep = std::make_shared<srep::SRep>(srep::io::ReadSRep(filename));
  this->UpdateSRepWorld(this->SRepTransform);
  this->Modified();
}

//----------------------------------------------------------------------------
bool vtkMRMLSRepNode::WriteSRepToFiles(const std::string& headerFilename,
                                       const std::string& upFilename,
                                       const std::string& downFilename,
                                       const std::string& crestFilename)
{
  if (!this->HasSRep()) {
    return false;
  }

  auto srep = this->GetSRep();
  srep::io::WriteSRep(*srep, headerFilename, upFilename, downFilename, crestFilename);
  return true;
}

//----------------------------------------------------------------------------
bool vtkMRMLSRepNode::HasSRep() const {
  return static_cast<bool>(this->SRep);
}

//----------------------------------------------------------------------------
const srep::SRep* vtkMRMLSRepNode::GetSRep() const {
  return this->SRep.get();
}

//----------------------------------------------------------------------------
const srep::SRep* vtkMRMLSRepNode::GetSRepWorld() const {
  return this->SRepWorld.get();
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
void vtkMRMLSRepNode::GetSRepBounds(const srep::SRep* srep, double bounds[6]) {
  vtkBoundingBox box;

  if (!srep) {
    box.GetBounds(bounds);
    return;
  }

  srep::foreachPoint(*srep, [&box](const srep::SkeletalPoint& point) {
    box.AddPoint(point.GetUpSpoke().GetBoundaryPoint().AsArray().data());
    box.AddPoint(point.GetDownSpoke().GetBoundaryPoint().AsArray().data());
    if (point.IsCrest()) {
      box.AddPoint(point.GetCrestSpoke().GetBoundaryPoint().AsArray().data());
    }
  });
  box.GetBounds(bounds);
}

//----------------------------------------------------------------------------
void vtkMRMLSRepNode::CopyContent(vtkMRMLNode* anode, bool deepCopy/*=true*/) {
  MRMLNodeModifyBlocker blocker(this);
  Superclass::CopyContent(anode, deepCopy);

  vtkMRMLSRepNode* node = vtkMRMLSRepNode::SafeDownCast(anode);
  if (node) {
    if (deepCopy) {
      if (node->SRep) {
        if (!this->SRep) {
          this->SRep = std::make_shared<srep::SRep>();
        }
        *this->SRep = *node->SRep;
      } else {
        this->SRep.reset();
      }
    } else {
      //shallow copy
      this->SRep = node->SRep;
    }
  }
}

//---------------------------------------------------------------------------
bool vtkMRMLSRepNode::CanApplyNonLinearTransforms() const
{
  return true;
}

//---------------------------------------------------------------------------
void vtkMRMLSRepNode::UpdateSRepWorld(vtkAbstractTransform* transform) {
  if (!this->SRep) {
    return;
  }

  if (!transform) {
    // no transform, both are the same. Shallow copy.
    this->SRepWorld = this->SRep;
    return;
  }

  using namespace srep;
  const auto& grid = this->SRep->GetSkeletalPoints();
  SRep::SkeletalGrid transformedGrid(grid.size(), std::vector<SkeletalPoint>(grid[0].size()));

  for (size_t i = 0; i < grid.size(); ++i) {
    for (size_t j  = 0; j < grid[i].size(); ++j) {
      const auto& skeletalPoint = grid[i][j];
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

        transformedGrid[i][j] = SkeletalPoint(transformedUpSpoke, transformedDownSpoke, transformedCrestSpoke);
      } else {
        transformedGrid[i][j] = SkeletalPoint(transformedUpSpoke, transformedDownSpoke);
      }
    }
  }

  this->SRepWorld = std::make_shared<srep::SRep>(transformedGrid);
}

//---------------------------------------------------------------------------
void vtkMRMLSRepNode::ApplyTransform(vtkAbstractTransform* transform)
{
  if (!this->SRep) {
    return;
  }

  this->UpdateSRepWorld(transform);
  if (this->SRep.get() != this->SRepWorld.get()) {
    // deep copy unless transform == nullptr so they are the same
    *this->SRep = *this->SRepWorld;
  }
  this->Modified();
}

//---------------------------------------------------------------------------
void vtkMRMLSRepNode::OnTransformNodeReferenceChanged(vtkMRMLTransformNode* transformNode) {

  //this next line is a GetTransformToWorld one-liner that works even if this->GetParentTransformNode is nullptr
  vtkMRMLTransformNode::GetTransformBetweenNodes(this->GetParentTransformNode(), nullptr, this->SRepTransform);
  this->UpdateSRepWorld(this->SRepTransform);
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
    this->UpdateSRepWorld(this->SRepTransform);
  }

  Superclass::ProcessMRMLEvents(caller, event, callData);
}
