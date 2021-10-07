#include "vtkMRMLSRepNode.h"
#include <vtkMRMLScene.h>

#include <srep/SRep.h>
#include <srep/SRepIO.h>

//----------------------------------------------------------------------------
vtkMRMLNodeNewMacro(vtkMRMLSRepNode);

vtkMRMLSRepNode::vtkMRMLSRepNode()
    : vtkMRMLDisplayableNode()
    , SRep()
{}

vtkMRMLSRepNode::~vtkMRMLSRepNode() = default;

void vtkMRMLSRepNode::LoadSRepFromFile(const std::string& filename) {
    this->SRep.reset(new srep::SRep(srep::io::ReadSRep(filename)));
    this->Modified();
}

bool vtkMRMLSRepNode::HasSRep() const {
    return static_cast<bool>(this->SRep);
}

const srep::SRep* vtkMRMLSRepNode::GetSRep() const {
  if (this->SRep) {
    return this->SRep.get();
  }
  return nullptr;
}

void vtkMRMLSRepNode::PrintSelf(ostream& os, vtkIndent indent) {
  Superclass::PrintSelf(os,indent);

  os << "SRep";
}

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

vtkMRMLSRepDisplayNode* vtkMRMLSRepNode::GetSRepDisplayNode() {
    auto* dispNode = this->GetDisplayNode();
    if (dispNode) {
        return vtkMRMLSRepDisplayNode::SafeDownCast(dispNode);
    }
    return nullptr;
}
