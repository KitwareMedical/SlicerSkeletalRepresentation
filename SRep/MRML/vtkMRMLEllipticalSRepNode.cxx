#include "vtkMRMLEllipticalSRepNode.h"
#include <vtkAbstractTransform.h>
#include <srep/SRepIO.h>

//----------------------------------------------------------------------------
vtkMRMLNodeNewMacro(vtkMRMLEllipticalSRepNode);

//----------------------------------------------------------------------------
vtkMRMLEllipticalSRepNode::vtkMRMLEllipticalSRepNode()
  : vtkMRMLSRepNode()
  , SRep()
  , SRepWorld()
{}

//----------------------------------------------------------------------------
vtkMRMLEllipticalSRepNode::~vtkMRMLEllipticalSRepNode() = default;

//----------------------------------------------------------------------------
const srep::EllipticalSRep* vtkMRMLEllipticalSRepNode::GetEllipticalSRep() const {
  return this->SRep.get();
}

//----------------------------------------------------------------------------
const srep::EllipticalSRep* vtkMRMLEllipticalSRepNode::GetEllipticalSRepWorld() const {
  return this->SRepWorld.get();
}

//----------------------------------------------------------------------------
void vtkMRMLEllipticalSRepNode::SetEllipticalSRep(std::unique_ptr<srep::EllipticalSRep> srep) {
  this->SRep = std::move(srep);
  this->UpdateSRepWorld();
  this->Modified();
}

//----------------------------------------------------------------------------
const srep::MeshSRepInterface* vtkMRMLEllipticalSRepNode::GetSRep() const {
  return this->SRep.get();
}

//----------------------------------------------------------------------------
const srep::MeshSRepInterface* vtkMRMLEllipticalSRepNode::GetSRepWorld() const {
  return this->SRepWorld.get();
}

//---------------------------------------------------------------------------
void vtkMRMLEllipticalSRepNode::DoUpdateSRepWorld(vtkAbstractTransform* transform) {
  if (!this->SRep) {
    return;
  }

  if (!transform) {
    // no transform, both are the same. Shallow copy.
    this->SRepWorld = this->SRep;
    return;
  }
  this->SRepWorld = std::make_shared<srep::EllipticalSRep>(
    TransformSkeletalPoints(this->SRep->GetSkeleton(), transform));
}

//---------------------------------------------------------------------------
void vtkMRMLEllipticalSRepNode::ApplyTransform(vtkAbstractTransform* transform)
{
  if (!this->SRep) {
    return;
  }

  this->DoUpdateSRepWorld(transform);
  if (this->SRep.get() != this->SRepWorld.get()) {
    // deep copy unless transform == nullptr so they are the same
    this->SRep = std::shared_ptr<srep::EllipticalSRep>(this->SRepWorld->Clone());
  }
  this->Modified();
}


//----------------------------------------------------------------------------
void vtkMRMLEllipticalSRepNode::CopyContent(vtkMRMLNode* anode, bool deepCopy/*=true*/) {
  MRMLNodeModifyBlocker blocker(this);
  Superclass::CopyContent(anode, deepCopy);

  auto* node = vtkMRMLEllipticalSRepNode::SafeDownCast(anode);
  if (node) {
    if (deepCopy) {
      if (node->SRep) {
        this->SRep = std::shared_ptr<srep::EllipticalSRep>(node->SRep->Clone());
      } else {
        this->SRep.reset();
      }
    } else {
      //shallow copy
      this->SRep = node->SRep;
    }
    this->UpdateSRepWorld();
  }
}
