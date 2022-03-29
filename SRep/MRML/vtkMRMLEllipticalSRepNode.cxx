#include "vtkMRMLEllipticalSRepNode.h"
#include <vtkAbstractTransform.h>
#include <vtkCommand.h>

//----------------------------------------------------------------------------
vtkEllipticalSRep* TransformSRep(vtkEllipticalSRep* srep, vtkAbstractTransform* transform) {
  auto transformed = SmartTransformSRep(srep, transform);
  if (transformed) {
    transformed->Register(nullptr);
  }
  return transformed;
}

//----------------------------------------------------------------------------
vtkSmartPointer<vtkEllipticalSRep> SmartTransformSRep(vtkEllipticalSRep* srep, vtkAbstractTransform* transform) {
  if (!srep) {
    return nullptr;
  }

  auto transformed = vtkSmartPointer<vtkEllipticalSRep>::Take(srep->Clone());
  if (transform) {
    using IndexType = vtkEllipticalSRep::IndexType;

    for (IndexType l = 0; l < transformed->GetNumberOfLines(); ++l) {
      for (IndexType s = 0; s < transformed->GetNumberOfSteps(); ++s) {
        auto* skeletalPoint = transformed->GetSkeletalPoint(l, s);
        std::vector<vtkSRepSpoke*> spokes;
        spokes.push_back(skeletalPoint->GetUpSpoke());
        spokes.push_back(skeletalPoint->GetDownSpoke());
        if (skeletalPoint->IsCrest()) {
          spokes.push_back(skeletalPoint->GetCrestSpoke());
        }

        for (auto* spoke : spokes) {
          double transformedBoundary[3];
          double transformedSkeletal[3];

          transform->TransformPoint(spoke->GetSkeletalPoint().AsArray().data(), transformedSkeletal);
          transform->TransformPoint(spoke->GetBoundaryPoint().AsArray().data(), transformedBoundary);

          spoke->SetSkeletalPoint(srep::Point3d(transformedSkeletal));
          spoke->SetDirectionAndMagnitude(srep::Vector3d(srep::Point3d(transformedSkeletal), srep::Point3d(transformedBoundary)));
        }
      }
    }
  }

  return transformed;
}

//----------------------------------------------------------------------------
vtkMRMLNodeNewMacro(vtkMRMLEllipticalSRepNode);

//----------------------------------------------------------------------------
vtkMRMLEllipticalSRepNode::vtkMRMLEllipticalSRepNode()
  : vtkMRMLSRepNode()
  , SRep()
  , SRepObservationTag()
  , SRepWorld()
{}

//----------------------------------------------------------------------------
vtkMRMLEllipticalSRepNode::~vtkMRMLEllipticalSRepNode() {
  if (this->SRep) {
    this->SRep->RemoveObserver(this->SRepObservationTag);
  }
};

//----------------------------------------------------------------------------
const vtkEllipticalSRep* vtkMRMLEllipticalSRepNode::GetEllipticalSRep() const {
  return this->SRep;
}

//----------------------------------------------------------------------------
vtkEllipticalSRep* vtkMRMLEllipticalSRepNode::GetEllipticalSRep() {
  return this->SRep;
}

//----------------------------------------------------------------------------
const vtkEllipticalSRep* vtkMRMLEllipticalSRepNode::GetEllipticalSRepWorld() const {
  return this->SRepWorld;
}

//----------------------------------------------------------------------------
void vtkMRMLEllipticalSRepNode::SetEllipticalSRep(vtkEllipticalSRep* srep) {
  if (this->SRep) {
    this->SRep->RemoveObserver(this->SRepObservationTag);
  }
  this->SRep = srep;
  if (this->SRep) {
    this->SRepObservationTag = this->SRep->AddObserver(vtkCommand::ModifiedEvent, this, &vtkMRMLEllipticalSRepNode::onSRepModified);
  }
  this->UpdateSRepWorld();
  this->Modified();
}

//----------------------------------------------------------------------------
const vtkMeshSRepInterface* vtkMRMLEllipticalSRepNode::GetSRep() const {
  return this->SRep;
}

//----------------------------------------------------------------------------
const vtkMeshSRepInterface* vtkMRMLEllipticalSRepNode::GetSRepWorld() const {
  return this->SRepWorld;
}

//---------------------------------------------------------------------------
void vtkMRMLEllipticalSRepNode::DoUpdateSRepWorld(vtkAbstractTransform* transform) {
  if (!this->SRep) {
    return;
  }

  this->SRepWorld = SmartTransformSRep(this->SRep, transform);
}

//---------------------------------------------------------------------------
void vtkMRMLEllipticalSRepNode::ApplyTransform(vtkAbstractTransform* transform)
{
  if (!this->SRep) {
    return;
  }

  this->DoUpdateSRepWorld(transform);
  this->SetEllipticalSRep(vtkSmartPointer<vtkEllipticalSRep>::Take(this->SRepWorld->Clone()));
  // SetEllipticalSRep will call modified
}


//----------------------------------------------------------------------------
void vtkMRMLEllipticalSRepNode::CopyContent(vtkMRMLNode* anode, bool deepCopy/*=true*/) {
  MRMLNodeModifyBlocker blocker(this);
  Superclass::CopyContent(anode, deepCopy);

  auto* node = vtkMRMLEllipticalSRepNode::SafeDownCast(anode);
  if (node) {
    if (deepCopy) {
      if (node->SRep) {
        this->SetEllipticalSRep(vtkSmartPointer<vtkEllipticalSRep>::Take(node->SRep->Clone()));
      } else {
        this->SetEllipticalSRep(nullptr);
      }
    } else {
      //shallow copy
      this->SetEllipticalSRep(node->SRep);
    }
  }
}

//----------------------------------------------------------------------------
void vtkMRMLEllipticalSRepNode::onSRepModified(vtkObject* /*caller*/, unsigned long /*event*/, void* /*callData*/) {
  this->UpdateSRepWorld();
  this->Modified();
}
