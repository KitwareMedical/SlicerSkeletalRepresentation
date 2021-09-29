#include "vtkMRMLSRepDisplayNode.h"
#include "vtkMRMLSRepNode.h"

//--------------------------------------------------------------------------
// vtkMRMLSRepDisplayNode::VisibilityHelper
//--------------------------------------------------------------------------
vtkMRMLSRepDisplayNode::VisibilityHelper::VisibilityHelper(const bool allOnOff)
    : topSpokes(allOnOff)
    , bottomSpokes(allOnOff)
    , crestSpokes(allOnOff)
    , crestCurve(allOnOff)
    , medialMesh(allOnOff)
    , connectionToFoldCurve(allOnOff)
{}

vtkMRMLSRepDisplayNode::VisibilityHelper::VisibilityHelper()
    : VisibilityHelper(true)
{}

//--------------------------------------------------------------------------
// vtkMRMLSRepDisplayNode
//--------------------------------------------------------------------------
vtkMRMLNodeNewMacro(vtkMRMLSRepDisplayNode);

vtkMRMLSRepDisplayNode::vtkMRMLSRepDisplayNode()
    : vtkMRMLDisplayNode()
    , ThreeDVisibility()
{}

vtkMRMLSRepDisplayNode::~vtkMRMLSRepDisplayNode() = default;

vtkMRMLSRepNode* vtkMRMLSRepDisplayNode::GetSRepNode() {
    return vtkMRMLSRepNode::SafeDownCast(this->GetDisplayableNode());
}

int vtkMRMLSRepDisplayNode::GetVisibility3D() {
    return static_cast<int>(this->ThreeDVisibility.overall);
}

void vtkMRMLSRepDisplayNode::SetVisibility3D(const int visible) {
    this->ThreeDVisibility.overall = visible;
}

int vtkMRMLSRepDisplayNode::GetVisibility2D() {
    return false;
}

void vtkMRMLSRepDisplayNode::SetVisibility2D(int /*visible*/) {}

void vtkMRMLSRepDisplayNode::PrintSelf(ostream& os, vtkIndent indent)
{
  Superclass::PrintSelf(os,indent);

  os << "SRep";
}
