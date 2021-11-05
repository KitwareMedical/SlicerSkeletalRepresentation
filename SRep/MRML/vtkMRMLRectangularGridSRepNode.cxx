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

  using namespace srep;
  const auto& grid = this->SRep->GetSkeletalPoints();
  RectangularGridSRep::SkeletalGrid transformedGrid(grid.size(), std::vector<SkeletalPoint>(grid[0].size()));

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

  this->SRepWorld = std::make_shared<srep::RectangularGridSRep>(std::move(transformedGrid));
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
        if (!this->SRep) {
          this->SRep = std::make_shared<srep::RectangularGridSRep>();
        }
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
