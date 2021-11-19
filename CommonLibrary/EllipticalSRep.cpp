#include <srep/EllipticalSRep.h>
#include <memory>

namespace srep {

EllipticalSRep::EllipticalSRep()
  : Skeleton()
  , SkeletonAsMesh()
{}

EllipticalSRep::EllipticalSRep(UnrolledEllipticalGrid skeleton)
  : Skeleton(std::move(skeleton))
  , SkeletonAsMesh()
{
  Validate(this->Skeleton);
  CreateMeshRepresentation();
}

const EllipticalSRep::UnrolledEllipticalGrid& EllipticalSRep::GetSkeleton() const {
  return this->Skeleton;
}

util::owner<EllipticalSRep*> EllipticalSRep::Clone() const {
    // use a unique ptr right up until the very end in case there are any exceptions
    std::unique_ptr<EllipticalSRep> cloned(new EllipticalSRep);
    cloned->Skeleton = this->Skeleton;
    cloned->SkeletonAsMesh = this->SkeletonAsMesh;
    return cloned.release();
}

bool EllipticalSRep::IsEmpty() const {
    return this->Skeleton.empty();
}

const SpokeMesh& EllipticalSRep::GetUpSpokes() const {
    return this->SkeletonAsMesh.UpSpokes;
}
const SpokeMesh& EllipticalSRep::GetDownSpokes() const {
    return this->SkeletonAsMesh.DownSpokes;
}
const SpokeMesh& EllipticalSRep::GetCrestSpokes() const {
    return this->SkeletonAsMesh.CrestSpokes;
}
const std::vector<EllipticalSRep::UpDownIndices>& EllipticalSRep::GetCrestSkeletalConnections() const {
    return this->SkeletonAsMesh.CrestSkeletalConnections;
}
const std::vector<EllipticalSRep::UpDownIndices>& EllipticalSRep::GetSpine() const {
    return this->SkeletonAsMesh.Spine;
}

void EllipticalSRep::Validate(const UnrolledEllipticalGrid& skeleton) {
  // empty is valid
  if (skeleton.empty()) {
    return;
  }

  // TODO: Should we check that the points that are supposed to be duplicates are actually the same?

  // Even number of "cols" (i.e. points on both sides of the spine)
  if (skeleton.size() % 2 != 0) {
    throw InvalidSkeletalGridException("Must have even number of rows (spine points)");
  }
  // Same number of "cols" per "row" (i.e. same number of steps + 1 (for crest) out of each spine point)
  // and the number of "cols" (steps) >= 2. 2 because 1 step minimum + 1 for crest
  const auto numStepsPlusCrest = skeleton[0].size();
  if (numStepsPlusCrest < 2) {
    throw InvalidSkeletalGridException("Expected at least two columns for a step and a crest");
  }
  for (const auto& line : skeleton) {
    if (line.size() != numStepsPlusCrest) {
      throw InvalidSkeletalGridException("All lines going out from the spine must have the same number of steps");
    }
  }
  // Bottom "row" are all crest points while nothing else is
  for (const auto& line : skeleton) {
    for (size_t s = 0; s < line.size(); ++s) {
      if (s == line.size() - 1 && !line[s].IsCrest()) {
        throw InvalidSkeletalGridException("Expected the last step of each line to be a crest point (" + std::to_string(s) + ")");
      }
      if (s < line.size() - 1 && line[s].IsCrest()) {
        throw InvalidSkeletalGridException("Expected all steps except the last to not be crest points (" + std::to_string(s) + ")");
      }
    }
  }
}

void EllipticalSRep::CreateMeshRepresentation() {
  this->SkeletonAsMesh.UpSpokes.Clear();
  this->SkeletonAsMesh.DownSpokes.Clear();
  this->SkeletonAsMesh.CrestSpokes.Clear();
  this->SkeletonAsMesh.CrestSkeletalConnections.clear();
  this->SkeletonAsMesh.Spine.clear();

  if (this->Skeleton.empty()) {
    return;
  }

  //insert in following order
  // 1) spine - left to right on top side - no duplicate points
  // 2) other points - [line][step] in skeleton order

  // +1 because we need the rightmost line
  const size_t numSpinePointsWithoutDuplicates = (this->Skeleton.size() / 2) + 1;
  const size_t numStepsPlusCrest = this->Skeleton[0].size();
  const size_t crestStepIndex = numStepsPlusCrest - 1;
  const size_t numLines = this->Skeleton.size();

  const auto toUpDownMeshIndex = [numSpinePointsWithoutDuplicates, numStepsPlusCrest, numLines]
    (const size_t line, const size_t step)
  {
    if (step == 0) {
      if (line < numSpinePointsWithoutDuplicates) {
        return line;
      } else {
        return numLines - line;
      }
    } else {
      // all the -1s are because we are stripping off step 0 (the spine)
      return numSpinePointsWithoutDuplicates + (line * (numStepsPlusCrest-1)) + step-1;
    }
  };

  // do the spine first - need to ignore duplicate points
  // all duplicate points are after the second pole
  for (size_t line = 0; line < numSpinePointsWithoutDuplicates; ++line) {
    std::vector<IndexType> neighbors;
    //neighbors along the spine
    if (line > 0) neighbors.push_back(toUpDownMeshIndex(line-1, 0));
    if (line < numSpinePointsWithoutDuplicates-1) neighbors.push_back(toUpDownMeshIndex(line+1, 0));
    //neighbor down my line, if statement should always be true
    if (numStepsPlusCrest > 1) neighbors.push_back(toUpDownMeshIndex(line, 1));
    //neighbor down my opposite line, if it exists
    if (0 < line && line < numSpinePointsWithoutDuplicates-1 && numStepsPlusCrest > 1) {
      //opposite line is numLines - 1 - (line - 1) which equates to numLines - line
      // numLines - 1 --> to go from size to last index
      // (line - 1) --> to account for index 0 being the pole
      neighbors.push_back(toUpDownMeshIndex(numLines-line, 1));
    }

    const auto& skeletalPoint = this->Skeleton[line][0];
    this->SkeletonAsMesh.UpSpokes.AddSpoke(skeletalPoint.GetUpSpoke(), neighbors);
    this->SkeletonAsMesh.DownSpokes.AddSpoke(skeletalPoint.GetDownSpoke(), neighbors);
    const auto index = toUpDownMeshIndex(line, 0);
    this->SkeletonAsMesh.Spine.push_back(std::make_pair(index, index));
  }

  for (size_t line = 0; line < this->Skeleton.size(); ++line) {
    // no duplicate points because we aren't on the spine
    for (size_t step = 1; step < this->Skeleton[line].size(); ++step) {
      std::vector<IndexType> neighbors;
      //neighbors from neighboring lines
      neighbors.push_back(toUpDownMeshIndex((numLines + line + 1) % numLines, step));
      neighbors.push_back(toUpDownMeshIndex((numLines + line - 1) % numLines, step));
      //outer neighbor
      if (step < this->Skeleton[line].size() - 1) {
        neighbors.push_back(toUpDownMeshIndex(line, step + 1));
      }
      //inner neighbor
      neighbors.push_back(toUpDownMeshIndex(line, step - 1));

      const auto& skeletalPoint = this->Skeleton[line][step];
      this->SkeletonAsMesh.UpSpokes.AddSpoke(skeletalPoint.GetUpSpoke(), neighbors);
      this->SkeletonAsMesh.DownSpokes.AddSpoke(skeletalPoint.GetDownSpoke(), neighbors);
    }
  }

  // crest spokes and connections
  for (size_t line = 0; line < this->Skeleton.size(); ++line) {
    const auto& skeletalPoint = this->Skeleton[line][crestStepIndex];
    std::vector<IndexType> neighbors;
    neighbors.push_back((numLines + line - 1) % numLines);
    neighbors.push_back((numLines + line + 1) % numLines);

    this->SkeletonAsMesh.CrestSpokes.AddSpoke(skeletalPoint.GetCrestSpoke(), neighbors);
    const auto index = toUpDownMeshIndex(line, crestStepIndex);
    this->SkeletonAsMesh.CrestSkeletalConnections.push_back(std::make_pair(index, index));
  }
}

} // namespace srep
