#include "vtkMRMLRectangularGridSRepNode.h"
#include <vtkAbstractTransform.h>
#include <srep/SRepIO.h>

//----------------------------------------------------------------------------
vtkMRMLNodeNewMacro(vtkMRMLRectangularGridSRepNode);

//----------------------------------------------------------------------------
vtkMRMLRectangularGridSRepNode::vtkMRMLRectangularGridSRepNode()
  : vtkMRMLSRepNode()
  , SRep()
  , SRepWorld()
{}

//----------------------------------------------------------------------------
vtkMRMLRectangularGridSRepNode::~vtkMRMLRectangularGridSRepNode() = default;

//----------------------------------------------------------------------------
const srep::RectangularGridSRep* vtkMRMLRectangularGridSRepNode::GetRectangularGridSRep() const {
  return this->SRep.get();
}

//----------------------------------------------------------------------------
const srep::RectangularGridSRep* vtkMRMLRectangularGridSRepNode::GetRectangularGridSRepWorld() const {
  return this->SRepWorld.get();
}

//----------------------------------------------------------------------------
void vtkMRMLRectangularGridSRepNode::SetRectangularGridSRep(std::unique_ptr<srep::RectangularGridSRep> srep) {
  this->SRep = std::move(srep);
  this->UpdateSRepWorld();
  this->Modified();
}

//----------------------------------------------------------------------------
const srep::MeshSRepInterface* vtkMRMLRectangularGridSRepNode::GetSRep() const {
  return this->SRep.get();
}

//----------------------------------------------------------------------------
const srep::MeshSRepInterface* vtkMRMLRectangularGridSRepNode::GetSRepWorld() const {
  return this->SRepWorld.get();
}

//---------------------------------------------------------------------------
void vtkMRMLRectangularGridSRepNode::DoUpdateSRepWorld(vtkAbstractTransform* transform) {
  if (!this->SRep) {
    return;
  }

  if (!transform) {
    // no transform, both are the same. Shallow copy.
    this->SRepWorld = this->SRep;
    return;
  }

  this->SRepWorld = std::make_shared<srep::RectangularGridSRep>(TransformSkeletalPoints(this->SRep->GetSkeletalPoints(), transform));
}

//---------------------------------------------------------------------------
void vtkMRMLRectangularGridSRepNode::ApplyTransform(vtkAbstractTransform* transform)
{
  if (!this->SRep) {
    return;
  }

  this->DoUpdateSRepWorld(transform);
  if (this->SRep.get() != this->SRepWorld.get()) {
    // deep copy unless transform == nullptr so they are the same
    this->SRep = std::shared_ptr<srep::RectangularGridSRep>(this->SRepWorld->Clone());
  }
  this->Modified();
}


//----------------------------------------------------------------------------
void vtkMRMLRectangularGridSRepNode::CopyContent(vtkMRMLNode* anode, bool deepCopy/*=true*/) {
  MRMLNodeModifyBlocker blocker(this);
  Superclass::CopyContent(anode, deepCopy);

  auto* node = vtkMRMLRectangularGridSRepNode::SafeDownCast(anode);
  if (node) {
    if (deepCopy) {
      if (node->SRep) {
        this->SRep = std::shared_ptr<srep::RectangularGridSRep>(node->SRep->Clone());
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

//----------------------------------------------------------------------------
void vtkMRMLRectangularGridSRepNode::LoadRectangularGridSRepFromFile(const std::string& filename) {
  this->SRep = srep::io::ReadRectangularGridSRep(filename);
  this->UpdateSRepWorld();
  this->Modified();
}

//----------------------------------------------------------------------------
bool vtkMRMLRectangularGridSRepNode::WriteRectangularGridSRepToFiles(const std::string& headerFilename,
                                       const std::string& upFilename,
                                       const std::string& downFilename,
                                       const std::string& crestFilename)
{
  if (!this->HasSRep()) {
    return false;
  }

  auto srep = this->GetRectangularGridSRep();
  if (srep) {
    srep::io::WriteSRep(*srep, headerFilename, upFilename, downFilename, crestFilename);
  }
  return true;
}
