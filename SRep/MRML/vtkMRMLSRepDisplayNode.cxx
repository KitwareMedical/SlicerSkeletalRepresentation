#include "vtkMRMLSRepDisplayNode.h"
#include "vtkMRMLSRepNode.h"

//--------------------------------------------------------------------------
// vtkMRMLSRepDisplayNode::DisplayHelper
//--------------------------------------------------------------------------
vtkMRMLSRepDisplayNode::DisplayHelper::DisplayHelper()
    : visible(true)
    , color(vtkColor3ub{0,0,0})
{}

vtkMRMLSRepDisplayNode::DisplayHelper::DisplayHelper(bool visible_, const vtkColor3ub& color_)
    : visible(visible_)
    , color(color_)
{}

//--------------------------------------------------------------------------
// vtkMRMLSRepDisplayNode
//--------------------------------------------------------------------------
vtkMRMLNodeNewMacro(vtkMRMLSRepDisplayNode);

vtkMRMLSRepDisplayNode::vtkMRMLSRepDisplayNode()
    : vtkMRMLDisplayNode()
    , OverallVisibility(true)
    , UpSpokes(true, vtkColor3ub{255, 99, 71}) // Tomato
    , DownSpokes(true, vtkColor3ub{189, 252, 201}) // Mint
    , CrestSpokes(true, vtkColor3ub{255, 215, 0}) // Gold
    , CrestCurve(true, vtkColor3ub{255, 215, 0}) // Gold
    , SkeletalSheet(true, vtkColor3ub{255, 248, 220}) // Cornsilk
    , SkeletonToCrestConnection(true, vtkColor3ub{0, 0, 0}) // Black
    , RelativeThickness(0.001)
    , AbsoluteThickness(0.25)
    , UseAbsoluteThickness(false)
{}

vtkMRMLSRepDisplayNode::~vtkMRMLSRepDisplayNode() = default;

vtkMRMLSRepNode* vtkMRMLSRepDisplayNode::GetSRepNode() {
    return vtkMRMLSRepNode::SafeDownCast(this->GetDisplayableNode());
}

int vtkMRMLSRepDisplayNode::GetVisibility3D() {
    return static_cast<int>(this->OverallVisibility);
}

void vtkMRMLSRepDisplayNode::SetVisibility3D(const int visible) {
    if (visible != this->OverallVisibility) {
        this->OverallVisibility = visible;
        this->Modified();
    }
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

void vtkMRMLSRepDisplayNode::SetUpSpokeColor(const vtkColor3ub& color) {
    if (this->UpSpokes.color != color) {
        this->UpSpokes.color = color;
        this->Modified();
    }
}
const vtkColor3ub& vtkMRMLSRepDisplayNode::GetUpSpokeColor() const {
    return this->UpSpokes.color;
}

void vtkMRMLSRepDisplayNode::SetDownSpokeColor(const vtkColor3ub& color) {
    if (this->DownSpokes.color != color) {
        this->DownSpokes.color = color;
        this->Modified();
    }
}
const vtkColor3ub& vtkMRMLSRepDisplayNode::GetDownSpokeColor() const {
    return this->DownSpokes.color;
}

void vtkMRMLSRepDisplayNode::SetCrestSpokeColor(const vtkColor3ub& color) {
    if (this->CrestSpokes.color != color) {
        this->CrestSpokes.color = color;
        this->Modified();
    }
}
const vtkColor3ub& vtkMRMLSRepDisplayNode::GetCrestSpokeColor() const {
    return this->CrestSpokes.color;
}

void vtkMRMLSRepDisplayNode::SetCrestCurveColor(const vtkColor3ub& color) {
    if (this->CrestCurve.color != color) {
        this->CrestCurve.color = color;
        this->Modified();
    }
}
const vtkColor3ub& vtkMRMLSRepDisplayNode::GetCrestCurveColor() const {
    return this->CrestCurve.color;
}

void vtkMRMLSRepDisplayNode::SetSkeletalSheetColor(const vtkColor3ub& color) {
    if (this->SkeletalSheet.color != color) {
        this->SkeletalSheet.color = color;
        this->Modified();
    }
}
const vtkColor3ub& vtkMRMLSRepDisplayNode::GetSkeletalSheetColor() const {
    return this->SkeletalSheet.color;
}

void vtkMRMLSRepDisplayNode::SetSkeletonToCrestConnectionColor(const vtkColor3ub& color) {
    if (this->SkeletonToCrestConnection.color != color) {
        this->SkeletonToCrestConnection.color = color;
        this->Modified();
    }
}
const vtkColor3ub& vtkMRMLSRepDisplayNode::GetSkeletonToCrestConnectionColor() const {
    return this->SkeletonToCrestConnection.color;
}

void vtkMRMLSRepDisplayNode::SetRelativeThickness(double relativeThickness) {
    if (relativeThickness < 0 || 1 < relativeThickness) {
        vtkErrorMacro("Relative thickness must be between 0 and 1. Found " << relativeThickness);
        return;
    }
    if (this->RelativeThickness != relativeThickness) {
        this->UseAbsoluteThickness = false;
        this->RelativeThickness = relativeThickness;
        this->Modified();
    }
}
double vtkMRMLSRepDisplayNode::GetRelativeThickness() const {
    return this->RelativeThickness;
}

void vtkMRMLSRepDisplayNode::SetAbsoluteThickness(double absoluteThickness) {
    if (this->AbsoluteThickness != absoluteThickness) {
        this->UseAbsoluteThickness = true;
        this->AbsoluteThickness = absoluteThickness;
        this->Modified();
    }
}
double vtkMRMLSRepDisplayNode::GetAbsoluteThickness() const {
    return this->AbsoluteThickness;
}

bool vtkMRMLSRepDisplayNode::GetUseAbsoluteThickness() const {
    return this->UseAbsoluteThickness;
}
void vtkMRMLSRepDisplayNode::SetUseAbsoluteThickness(bool use) {
    if (this->UseAbsoluteThickness != use) {
        this->UseAbsoluteThickness = use;
        this->Modified();
    }
}
void vtkMRMLSRepDisplayNode::UseAbsoluteThicknessOn() {
    this->SetUseAbsoluteThickness(true);
}
void vtkMRMLSRepDisplayNode::UseAbsoluteThicknessOff() {
    this->SetUseAbsoluteThickness(false);
}
