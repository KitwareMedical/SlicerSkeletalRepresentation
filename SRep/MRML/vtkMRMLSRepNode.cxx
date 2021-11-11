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
