#include <srep/EllipticalSRep.h>
#include <algorithm>
#include <memory>
#include <sstream>

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

  const auto debugSkeleton = [](size_t line, size_t step, const UnrolledEllipticalGrid& skeleton) {
    //making a stringstream and returning a string isn't the most efficient, but this shouldn't
    //be called often
    std::stringstream ss;
    ss << "Error at (" << line << ", " << step << ") of skeleton size "
      << skeleton.size() << "x" << (skeleton.empty() ? 0 : skeleton[0].size());
    return ss.str();
  };

  for (size_t l = 0; l < skeleton.size(); ++l) {
    const auto& line = skeleton[l];
    for (size_t s = 0; s < line.size(); ++s) {
      if (s == line.size() - 1 && !line[s].IsCrest()) {
        throw InvalidSkeletalGridException("Expected the last step of each line to be a crest point. " + debugSkeleton(l,s,skeleton));
      }
      if (s < line.size() - 1 && line[s].IsCrest()) {
        throw InvalidSkeletalGridException("Expected all steps except the last to not be crest points. " + debugSkeleton(l,s,skeleton));
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

  const auto getNeighbors = [this, toUpDownMeshIndex](size_t line, size_t step) {
    const auto lineStepNeighbors = GetNeighbors(this->Skeleton, line, step);
    std::vector<IndexType> neighbors;
    std::transform(lineStepNeighbors.begin(), lineStepNeighbors.end(), std::back_inserter(neighbors),
      [toUpDownMeshIndex](const LineStep& ls) { return toUpDownMeshIndex(ls.line, ls.step); });
    return neighbors;
  };

  // do the spine first - need to ignore duplicate points
  // all duplicate points are after the second pole
  for (size_t line = 0; line < numSpinePointsWithoutDuplicates; ++line) {
    auto neighbors = getNeighbors(line, 0);

    const auto& skeletalPoint = this->Skeleton[line][0];
    this->SkeletonAsMesh.UpSpokes.AddSpoke(skeletalPoint.GetUpSpoke(), neighbors);
    this->SkeletonAsMesh.DownSpokes.AddSpoke(skeletalPoint.GetDownSpoke(), std::move(neighbors));
    const auto index = toUpDownMeshIndex(line, 0);
    this->SkeletonAsMesh.Spine.push_back(std::make_pair(index, index));
  }

  for (size_t line = 0; line < this->Skeleton.size(); ++line) {
    // no duplicate points because we aren't on the spine
    for (size_t step = 1; step < this->Skeleton[line].size(); ++step) {
      auto neighbors = getNeighbors(line, step);

      const auto& skeletalPoint = this->Skeleton[line][step];
      this->SkeletonAsMesh.UpSpokes.AddSpoke(skeletalPoint.GetUpSpoke(), neighbors);
      this->SkeletonAsMesh.DownSpokes.AddSpoke(skeletalPoint.GetDownSpoke(), std::move(neighbors));
    }
  }

  // crest spokes and connections
  for (size_t line = 0; line < this->Skeleton.size(); ++line) {
    const auto& skeletalPoint = this->Skeleton[line][crestStepIndex];

    //manually get neighbors here because we only want neighboring crests
    std::vector<IndexType> neighbors;
    neighbors.push_back((numLines + line - 1) % numLines);
    neighbors.push_back((numLines + line + 1) % numLines);

    this->SkeletonAsMesh.CrestSpokes.AddSpoke(skeletalPoint.GetCrestSpoke(), neighbors);
    const auto index = toUpDownMeshIndex(line, crestStepIndex);
    this->SkeletonAsMesh.CrestSkeletalConnections.push_back(std::make_pair(index, index));
  }
}

namespace {
  void verifyInGrid(const EllipticalSRep::UnrolledEllipticalGrid& grid, size_t line, size_t step) {
    if (line >= grid.size()) {
      throw std::invalid_argument("Line not in UnrolledEllipticalGrid");
    }
    if (step >= grid[line].size()) {
      throw std::invalid_argument("Step not in UnrolledEllipticalGrid");
    }
  }
} //namespace {}

std::pair<LineStep, bool> GetLeftNeighbor(const EllipticalSRep::UnrolledEllipticalGrid& grid, size_t line, size_t step) {
  verifyInGrid(grid, line, step);
  if (line == 0) {
    // left pole, left-right are steps
    if (step == grid[line].size() - 1) {
      return std::make_pair(LineStep(), false);
    }
    return std::make_pair(LineStep(line, step+1), true);
  } else if (line == grid.size() / 2) {
    //right pole, left-right are steps
    if (step == 0) {
      return std::make_pair(LineStep(), false);
    }
    return std::make_pair(LineStep(line, step-1), true);
  } else {
    if (line < grid.size() / 2) {
      //top half of the ellipse
      return std::make_pair(LineStep(line - 1, step), true);
    } else {
      //bottom half of the ellipse
      return std::make_pair(LineStep((line + 1) % grid.size(), step), true);
    }
  }
}
std::pair<LineStep, bool> GetRightNeighbor(const EllipticalSRep::UnrolledEllipticalGrid& grid, size_t line, size_t step) {
  verifyInGrid(grid, line, step);
  if (line == 0) {
    // left pole, left-right are steps
    if (step == 0) {
      return std::make_pair(LineStep(), false);
    }
    return std::make_pair(LineStep(line, step-1), true);
  } else if (line == grid.size() / 2) {
    //right pole, left-right are steps
    if (step == grid[line].size() - 1) {
      return std::make_pair(LineStep(), false);
    }
    return std::make_pair(LineStep(line, step+1), true);
  } else {
    if (line < grid.size() / 2) {
      //top half of the ellipse
      return std::make_pair(LineStep(line + 1, step), true);
    } else {
      //bottom half of the ellipse
      return std::make_pair(LineStep((grid.size() + line - 1) % grid.size(), step), true);
    }
  }
}
std::pair<LineStep, bool> GetTopNeighbor(const EllipticalSRep::UnrolledEllipticalGrid& grid, size_t line, size_t step) {
  verifyInGrid(grid, line, step);
  if ((line == 0 || line == grid.size() / 2) && step == 0) {
    // poles have a spine point with no top or bottom neighbors
    return std::make_pair(LineStep(), false);
  }

  if (line < grid.size() / 2) {
    //top half of the ellipse
    if (step == grid[line].size() - 1) {
      //top edge
      return std::make_pair(LineStep(), false);
    }
    return std::make_pair(LineStep(line, step+1), true);
  } else {
    //bottom half of the ellipse
    if (step == 0) {
      //we are crossing the spine
      return std::make_pair(LineStep(grid.size() - line, 1), true);
    }
    return std::make_pair(LineStep(line, step-1), true);
  }
}
std::pair<LineStep, bool> GetBottomNeighbor(const EllipticalSRep::UnrolledEllipticalGrid& grid, size_t line, size_t step) {
  verifyInGrid(grid, line, step);
  if ((line == 0 || line == grid.size() / 2) && step == 0) {
    // poles have a spine point with no top or bottom neighbors
    return std::make_pair(LineStep(), false);
  }
  if (line < grid.size() / 2) {
    //top half of the ellipse
    if (step == 0) {
      //we are crossing the spine
      return std::make_pair(LineStep(grid.size() - line, 1), true);
    }
    return std::make_pair(LineStep(line, step-1), true);
  } else {
    //bottom half of the ellipse
    if (step == grid[line].size() - 1) {
      //bottom edge
      return std::make_pair(LineStep(), false);
    }
    return std::make_pair(LineStep(line, step+1), true);
  }
}

std::vector<LineStep> GetNeighbors(const EllipticalSRep::UnrolledEllipticalGrid& grid, size_t line, size_t step) {
  std::vector<LineStep> neighbors;
  const std::array<GetNeighborFunc, 4> funcs {
    GetTopNeighbor, GetBottomNeighbor, GetLeftNeighbor, GetRightNeighbor
  };

  for (auto f : funcs) {
    const auto neighbor = f(grid, line, step);
    if (neighbor.second) {
      neighbors.push_back(neighbor.first);
    }
  }
  return neighbors;
}

LineStep::LineStep()
  : line(0)
  , step(0)
{}

LineStep::LineStep(size_t line_, size_t step_)
  : line(line_)
  , step(step_)
{}

bool operator<(const LineStep& a, const LineStep& b) {
  return a.line != b.line ? a.line < b.line : a.step < b.step;
}
bool operator>(const LineStep& a, const LineStep& b) {
  return b > a;
}
bool operator<=(const LineStep& a, const LineStep& b) {
  return !(b > a);
}
bool operator>=(const LineStep& a, const LineStep& b) {
  return !(a < b);
}
bool operator==(const LineStep& a, const LineStep& b) {
  return !(a < b) && !(b < a);
}
bool operator!=(const LineStep& a, const LineStep& b) {
  return !(a == b);
}

std::ostream& operator<<(std::ostream& os, const LineStep& ls) {
  os << "(" << ls.line << ", " << ls.step << ")";
  return os;
}

} // namespace srep
