#include "vtkEllipticalSRep.h"

#include <algorithm>
#include <sstream>

#include <vtkCommand.h>
#include <vtkObjectFactory.h>

//----------------------------------------------------------------------
vtkEllipticalSRep::ModifiedBlocker::ModifiedBlocker(vtkEllipticalSRep* srep)
  : Parent(srep)
{
  if(this->Parent) {
    this->Parent->BlockModify();
  }
}

//----------------------------------------------------------------------
vtkEllipticalSRep::ModifiedBlocker::~ModifiedBlocker()
{
  if(this->Parent) {
    this->Parent->UnblockModify();
  }
}

//----------------------------------------------------------------------
vtkStandardNewMacro(vtkEllipticalSRep);

//----------------------------------------------------------------------
vtkEllipticalSRep::vtkEllipticalSRep() = default;

//----------------------------------------------------------------------
vtkEllipticalSRep::~vtkEllipticalSRep() {
  this->Clear();
}

//----------------------------------------------------------------------
void vtkEllipticalSRep::PrintSelf(ostream& os, vtkIndent indent) {
  this->Superclass::PrintSelf(os, indent);
  os << indent << "Number of Lines: " << this->GetNumberOfLines() << std::endl;
  os << indent << "Number of Lines: " << this->GetNumberOfSteps() << std::endl;
}

//----------------------------------------------------------------------
bool vtkEllipticalSRep::IsEmpty() const {
  return GetNumberOfLines() == 0;
}

//----------------------------------------------------------------------
const vtkSRepSpokeMesh* vtkEllipticalSRep::GetUpSpokes() const {
  return this->SkeletonAsMesh.UpSpokes;
}

//----------------------------------------------------------------------
const vtkSRepSpokeMesh* vtkEllipticalSRep::GetDownSpokes() const {
  return this->SkeletonAsMesh.DownSpokes;
}

//----------------------------------------------------------------------
const vtkSRepSpokeMesh* vtkEllipticalSRep::GetCrestSpokes() const {
  return this->SkeletonAsMesh.CrestSpokes;
}

//----------------------------------------------------------------------
const std::vector<vtkEllipticalSRep::IndexType>&
vtkEllipticalSRep::GetCrestToUpSpokeConnections() const {
    return this->SkeletonAsMesh.CrestToUpSpokeConnections;
}

//----------------------------------------------------------------------
const std::vector<vtkEllipticalSRep::IndexType>&
vtkEllipticalSRep::GetCrestToDownSpokeConnections() const {
    return this->SkeletonAsMesh.CrestToDownSpokeConnections;
}

//----------------------------------------------------------------------
const std::vector<vtkEllipticalSRep::IndexType>&
vtkEllipticalSRep::GetUpSpine() const {
  return this->SkeletonAsMesh.UpSpine;
}

//----------------------------------------------------------------------
const std::vector<vtkEllipticalSRep::IndexType>&
vtkEllipticalSRep::GetDownSpine() const {
  return this->SkeletonAsMesh.DownSpine;
}

//----------------------------------------------------------------------
vtkEllipticalSRep::IndexType vtkEllipticalSRep::GetNumberOfLines() const {
  return this->Skeleton.size();
}

//----------------------------------------------------------------------
vtkEllipticalSRep::IndexType vtkEllipticalSRep::GetNumberOfSteps() const {
  return this->Skeleton.empty() ? 0 : this->Skeleton.front().size();
}
 
//----------------------------------------------------------------------
bool vtkEllipticalSRep::InBounds(IndexType line, IndexType step) const {
  return (0 <= line && line < GetNumberOfLines()) && (0 <= step && step < GetNumberOfSteps());
}

//----------------------------------------------------------------------
void vtkEllipticalSRep::CheckInBounds(IndexType line, IndexType step) const {
  if (!InBounds(line, step)) {
    throw std::out_of_range("Point (" + std::to_string(line) + ", " + std::to_string(step)
      + ") is outside of range (" + std::to_string(GetNumberOfLines()) + ", " + std::to_string(GetNumberOfSteps()) + ")");
  }
}

//----------------------------------------------------------------------
bool vtkEllipticalSRep::IsCrestStep(IndexType step) const {
  return step == GetNumberOfSteps() - 1;
}

//----------------------------------------------------------------------
bool vtkEllipticalSRep::CanSet(IndexType line, IndexType step, vtkSRepSkeletalPoint* skeletalPoint) const {
  return InBounds(line, step) && skeletalPoint && skeletalPoint->IsCrest() == IsCrestStep(step);
}

//----------------------------------------------------------------------
void vtkEllipticalSRep::CheckCanSet(IndexType line, IndexType step, vtkSRepSkeletalPoint* skeletalPoint) const {
  // this is kind of a copy of CanSet, but the extra info the exceptions is worth the duplication
  CheckInBounds(line, step);
  if (!skeletalPoint) {
    throw std::invalid_argument("Cannot set a nullptr skeletal point");
  }
  if (skeletalPoint->IsCrest() != IsCrestStep(step)) {
    std::stringstream ss;
    ss << "Cannot set " << (IsCrestStep(step) ? std::string("") : std::string("non-")) << "crest location "
       << "(" << line << ", " << step << ") "
       << "to a " << (skeletalPoint->IsCrest() ? std::string("") : std::string("non-")) << "crest skeletal point";
    throw std::invalid_argument(ss.str());
  }
  if (!CanSet(line, step, skeletalPoint)) {
    vtkWarningMacro("Unknown reason skeletal point cannot be set. This implies an error in this class.");
    throw std::runtime_error("Unknown reason skeletal point cannot be set");
  }
}

//----------------------------------------------------------------------
const vtkSRepSkeletalPoint* vtkEllipticalSRep::GetSkeletalPoint(IndexType line, IndexType step) const {
  CheckInBounds(line, step);
  return this->Skeleton[line][step];
}

//----------------------------------------------------------------------
vtkSRepSkeletalPoint* vtkEllipticalSRep::GetSkeletalPoint(IndexType line, IndexType step) {
  CheckInBounds(line, step);
  return this->Skeleton[line][step];
}

//----------------------------------------------------------------------
void vtkEllipticalSRep::SetSkeletalPoint(IndexType line, IndexType step, vtkSRepSkeletalPoint* skeletalPoint) {
  this->SetSkeletalPointNoMeshUpdate(line, step, skeletalPoint);

  // as far the mesh representation goes, nothing changes as far as neighbors because the neighbors are
  // indices. Just need to update the actual pointer in the mesh
  this->SkeletonAsMesh.UpSpokes->SetSpoke(this->LineStepToUpDownMeshIndex(line, step), this->Skeleton[line][step]->GetUpSpoke());
  this->SkeletonAsMesh.DownSpokes->SetSpoke(this->LineStepToUpDownMeshIndex(line, step), this->Skeleton[line][step]->GetDownSpoke());
  if (skeletalPoint->IsCrest()) {
    this->SkeletonAsMesh.CrestSpokes->SetSpoke(this->LineStepToUpDownMeshIndex(line, step), this->Skeleton[line][step]->GetCrestSpoke());
  }

  this->Modified();
}

//----------------------------------------------------------------------
void vtkEllipticalSRep::SetSkeletalPointNoMeshUpdate(IndexType line, IndexType step, vtkSRepSkeletalPoint* skeletalPoint) {
  CheckCanSet(line, step, skeletalPoint);

  // nullptr check here is so this works during resize
  if (this->Skeleton[line][step]) {
    this->Skeleton[line][step]->RemoveObserver(this->SkeletonObservationTags[line][step]);
  }
  this->Skeleton[line][step] = skeletalPoint;
  this->SkeletonObservationTags[line][step] =
    this->Skeleton[line][step]->AddObserver(vtkCommand::ModifiedEvent, this, &vtkEllipticalSRep::onSkeletalPointModified);
}

//----------------------------------------------------------------------
void vtkEllipticalSRep::TakeSkeletalPoint(IndexType line, IndexType step, vtkSRepSkeletalPoint* skeletalPoint) {
  this->SetSkeletalPoint(line, step, vtkSmartPointer<vtkSRepSkeletalPoint>::Take(skeletalPoint));
}

//----------------------------------------------------------------------
void vtkEllipticalSRep::Modified() {
  if (this->ModifiedBlocks == 0) {
    this->Superclass::Modified();
  } else {
    this->WasModifiedDuringBlock = true;
  }
}

//----------------------------------------------------------------------
void vtkEllipticalSRep::BlockModify() {
  ++this->ModifiedBlocks;
}

//----------------------------------------------------------------------
void vtkEllipticalSRep::UnblockModify() {
  // do <= 0 in case someone unblocks without blocking
  if (--this->ModifiedBlocks <= 0) {
    this->ModifiedBlocks = 0;
    if (this->WasModifiedDuringBlock) {
      this->Superclass::Modified();
    }
    this->WasModifiedDuringBlock = false;
  }
}

//----------------------------------------------------------------------
void vtkEllipticalSRep::Clear() {
  this->Resize(0, 0);
}

//----------------------------------------------------------------------
void vtkEllipticalSRep::Resize(IndexType lines, IndexType steps) {
  if (lines == GetNumberOfLines() && steps == GetNumberOfSteps()) {
    // nothing to do
    return;
  }

  const auto RemovePoint = [&](IndexType line, IndexType step) {
    if (this->Skeleton[line][step]) {
      this->Skeleton[line][step]->RemoveObserver(this->SkeletonObservationTags[line][step]);
      this->Skeleton[line][step] = nullptr;
    }
  };

  // don't emit a modify for each new/deleted point, just one at the end.
  ModifiedBlocker block(this);

  const auto makeEmptySkeletalPoint = [this](IndexType step) {
    auto ret =vtkSmartPointer<vtkSRepSkeletalPoint>::New();
    if (this->IsCrestStep(step)) {
      ret->SetCrestSpoke(vtkSmartPointer<vtkSRepSpoke>::New());
    }
    return ret;
  };

  // first resize by lines
  if (lines < GetNumberOfLines()) {
    //going down in size
    for (IndexType l = lines; l < GetNumberOfLines(); ++l) {
      for (IndexType s = 0; s < GetNumberOfSteps(); ++s) {
        RemovePoint(l,s);
      }
    }
    this->Skeleton.resize(lines);
    this->SkeletonObservationTags.resize(lines);
  } else if (lines > GetNumberOfLines()) {
    // going up in size
    const auto oldLines = GetNumberOfLines();
    this->Skeleton.resize(lines);
    this->SkeletonObservationTags.resize(lines);

    for (IndexType l = oldLines; l < static_cast<IndexType>(this->Skeleton.size()); ++l) {
      for (IndexType s = 0; s < GetNumberOfSteps(); ++s) {
        this->SetSkeletalPointNoMeshUpdate(l, s, makeEmptySkeletalPoint(s));
      }
    }
  }

  // next resize by steps
  if (steps < GetNumberOfSteps()) {
    //going down in size
    for (IndexType l = 0; l < GetNumberOfLines(); ++l) {
      for (IndexType s = steps; s < GetNumberOfSteps(); ++s) {
        RemovePoint(l,s);
      }
      this->Skeleton[l].resize(steps);
      this->SkeletonObservationTags[l].resize(steps);
    }
  } else if (steps > GetNumberOfSteps()) {
    // going up in size
    const auto oldSteps = GetNumberOfSteps();
    for (IndexType l = 0; l < GetNumberOfLines(); ++l) {
      this->Skeleton[l].resize(steps);
      this->SkeletonObservationTags[l].resize(steps);
      for (IndexType s = oldSteps; s < static_cast<IndexType>(this->Skeleton[l].size()); ++s) {
        this->SetSkeletalPointNoMeshUpdate(l, s, makeEmptySkeletalPoint(s));
      }
    }
  }

  this->CreateMeshRepresentation();
  this->Modified();
}

//----------------------------------------------------------------------
void vtkEllipticalSRep::onSkeletalPointModified(vtkObject */*caller*/, unsigned long /*event*/, void* /*callData*/) {
  this->Modified();
}

//----------------------------------------------------------------------
vtkEllipticalSRep* vtkEllipticalSRep::Clone() const {
  vtkNew<vtkEllipticalSRep> clone;

  clone->Resize(this->GetNumberOfLines(), this->GetNumberOfSteps());

  for (IndexType l = 0; l < this->GetNumberOfLines(); ++l) {
    for (IndexType s = 0; s < this->GetNumberOfSteps(); ++s) {
      // Note use of take because Clone returns an owning pointer
      clone->TakeSkeletalPoint(l, s, this->GetSkeletalPoint(l, s)->Clone());
    }
  }

  // update refcount so it doesn't go away after this function ends
  clone->Register(nullptr);
  return clone;
}

//----------------------------------------------------------------------
// need to implement SmartClone in terms of clone and not vice-versa because Clone is virtual
vtkSmartPointer<vtkEllipticalSRep> vtkEllipticalSRep::SmartClone() const {
  return vtkSmartPointer<vtkEllipticalSRep>::Take(this->Clone());
}

//----------------------------------------------------------------------
vtkEllipticalSRep::IndexType vtkEllipticalSRep::NumberOfSpinePointsWithoutDuplicates() const {
  // +1 because we need the rightmost line
  return this->IsEmpty() ? 0 : (this->Skeleton.size() / 2) + 1;
}

//----------------------------------------------------------------------
vtkSRepSpokeMesh::IndexType vtkEllipticalSRep::LineStepToUpDownMeshIndex(IndexType line, IndexType step) const {
  const auto numberOfSpinePointsWithoutDuplicates = NumberOfSpinePointsWithoutDuplicates();
  if (step == 0) {
    if (line < numberOfSpinePointsWithoutDuplicates) {
      return line;
    } else {
      return GetNumberOfLines() - numberOfSpinePointsWithoutDuplicates;
    }
  } else {
    // all the -1s are because we are stripping off step 0 (the spine)
    return numberOfSpinePointsWithoutDuplicates + (line * (GetNumberOfSteps() - 1)) + (step - 1);
  }
}

//----------------------------------------------------------------------
std::vector<vtkSRepSpokeMesh::IndexType> vtkEllipticalSRep::GetNeighbors(IndexType line, IndexType step) const {
  using LineStep = std::pair<size_t, size_t>;
  std::vector<LineStep> neighborLineSteps;
  {
    const auto prevLine = (GetNumberOfLines() + line - 1) % GetNumberOfLines();
    const auto nextLine = (GetNumberOfLines() + line + 1) % GetNumberOfLines();

    neighborLineSteps.push_back(std::make_pair(prevLine, step));
    neighborLineSteps.push_back(std::make_pair(nextLine, step));
    if (step > 0) {
      neighborLineSteps.push_back(std::make_pair(line, step-1));
    }
    if (step < GetNumberOfSteps() - 1) {
      neighborLineSteps.push_back(std::make_pair(line, step+1));
    }
  }

  std::vector<vtkSRepSpokeMesh::IndexType> neighbors;
  std::transform(neighborLineSteps.begin(), neighborLineSteps.end(), std::back_inserter(neighbors),
    [this](const LineStep& ls) { return this->LineStepToUpDownMeshIndex(ls.first, ls.second); });
  return neighbors;
}

//----------------------------------------------------------------------
void vtkEllipticalSRep::CreateMeshRepresentation() {
  this->SkeletonAsMesh.UpSpokes->Clear();
  this->SkeletonAsMesh.DownSpokes->Clear();
  this->SkeletonAsMesh.CrestSpokes->Clear();
  this->SkeletonAsMesh.CrestToUpSpokeConnections.clear();
  this->SkeletonAsMesh.UpSpine.clear();
  this->SkeletonAsMesh.CrestToDownSpokeConnections.clear();
  this->SkeletonAsMesh.DownSpine.clear();

  if (this->IsEmpty()) {
    return;
  }

  //insert in following order
  // 1) spine - left to right on top side - no duplicate points
  // 2) other points - [line][step] in skeleton order

  const IndexType crestStepIndex = GetNumberOfSteps() - 1;
  const auto numberOfSpinePointsWithoutDuplicates = this->NumberOfSpinePointsWithoutDuplicates();
  const auto numberOfLines = this->GetNumberOfLines();

  // do the spine first - need to ignore duplicate points
  // all duplicate points are after the second pole
  for (IndexType line = 0; line < numberOfSpinePointsWithoutDuplicates; ++line) {
    auto neighbors = this->GetNeighbors(line, 0);

    const auto& skeletalPoint = this->Skeleton[line][0];
    this->SkeletonAsMesh.UpSpokes->AddSpoke(skeletalPoint->GetUpSpoke(), neighbors);
    this->SkeletonAsMesh.DownSpokes->AddSpoke(skeletalPoint->GetDownSpoke(), std::move(neighbors));
    const auto index = LineStepToUpDownMeshIndex(line, 0);
    this->SkeletonAsMesh.UpSpine.push_back(index);
    this->SkeletonAsMesh.DownSpine.push_back(index);
  }

  for (IndexType line = 0; line < static_cast<IndexType>(this->Skeleton.size()); ++line) {
    // no duplicate points because we aren't on the spine
    for (IndexType step = 1; step < static_cast<IndexType>(this->Skeleton[line].size()); ++step) {
      auto neighbors = this->GetNeighbors(line, step);

      const auto& skeletalPoint = this->Skeleton[line][step];
      this->SkeletonAsMesh.UpSpokes->AddSpoke(skeletalPoint->GetUpSpoke(), neighbors);
      this->SkeletonAsMesh.DownSpokes->AddSpoke(skeletalPoint->GetDownSpoke(), std::move(neighbors));
    }
  }

  // crest spokes and connections
  for (IndexType line = 0; line < static_cast<IndexType>(this->Skeleton.size()); ++line) {
    const auto& skeletalPoint = this->Skeleton[line][crestStepIndex];

    //manually get neighbors here because we only want neighboring crests
    std::vector<IndexType> neighbors;
    neighbors.push_back((numberOfLines + line - 1) % numberOfLines);
    neighbors.push_back((numberOfLines + line + 1) % numberOfLines);

    this->SkeletonAsMesh.CrestSpokes->AddSpoke(skeletalPoint->GetCrestSpoke(), neighbors);
    const auto index = LineStepToUpDownMeshIndex(line, crestStepIndex);
    this->SkeletonAsMesh.CrestToUpSpokeConnections.push_back(index);
    this->SkeletonAsMesh.CrestToDownSpokeConnections.push_back(index);
  }
}
