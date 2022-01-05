/*==============================================================================

  Program: 3D Slicer

  Portions (c) Copyright Brigham and Women's Hospital (BWH) All Rights Reserved.

  See COPYRIGHT.txt
  or http://www.slicer.org/copyright/copyright.txt for details.

  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.

==============================================================================*/

#include "SRepInterpolation.h"
#include <algorithm>
#include <cmath>
#include <functional>
#include <sstream> //REMOVE

namespace sreplogic {
namespace detail {

namespace {
//----------------------------------------------------------------------------
// don't have to worry about losing precision on doubles and rounding them
size_t IntegerPower(const size_t x, const size_t y) {
  size_t ret = 1;
  for (size_t i = 0; i < y; ++i) {
    ret *= x;
  }
  return ret;
}

//----------------------------------------------------------------------------
double Clamp(double val, double min, double max) {
  return val < min ? min : (val > max ? max : val);
}

//----------------------------------------------------------------------------
bool IsPowerOfTwo(size_t val) {
  while (val > 1) {
    if (val % 2 == 1) {
      return false;
    }
    val = val / 2;
  }
  return true;
}

template <typename T, size_t N>
T MinValue(const std::array<T,N>& arr) {
  static_assert(N > 0, "Can't be zero");
  T minT = arr[0];
  for (const auto& v : arr) {
    if (v < minT) {
      minT = v;
    }
  }
  return minT;
}

} // namespace {}

//----------------------------------------------------------------------------
const srep::Spoke& SRepInterpolateHelper::GetSpoke(const Grid& grid, const LineStep& loc, SpokeType spokeType) {
  return grid[loc.line][loc.step].GetSpoke(spokeType);
}

//----------------------------------------------------------------------------
const srep::Spoke& SRepInterpolateHelper::GetInterpolatedSpoke(const LineStep& loc, SpokeType spokeType) {
  return this->GetSpoke(this->InterpolatedGrid, loc, spokeType);
}

//----------------------------------------------------------------------------
const srep::SkeletalPoint& SRepInterpolateHelper::GetInterpolatedSkeletalPoint(const LineStep& loc) {
  return this->InterpolatedGrid[loc.line][loc.step];
}

//----------------------------------------------------------------------------
srep::Vector3d SRepInterpolateHelper::ComputeLinewiseDerivative(const Grid& grid, const srep::LineStep& loc, SpokeType spokeType) {
  const auto prevLine = (loc.line + grid.size() - 1) % grid.size();
  const auto nextLine = (loc.line + grid.size() + 1) % grid.size();

  return SRepInterpolateHelper::ComputeDerivative(
    grid,
    loc,
    std::make_pair(LineStep(prevLine, loc.step), true),
    std::make_pair(LineStep(nextLine, loc.step), true),
    spokeType);
}

//----------------------------------------------------------------------------
srep::Vector3d SRepInterpolateHelper::ComputeStepwiseDerivative(const Grid& grid, const srep::LineStep& loc, SpokeType spokeType) {
  const auto prevOptLineStep = [&](){
    if (loc.step == 0) {
      return std::make_pair(LineStep(), false);
    }
    return std::make_pair(LineStep(loc.line, loc.step - 1), true);
  }();
  const auto nextOptLineStep = [&]() {
    if (loc.step == grid[loc.line].size() - 1) {
      return std::make_pair(LineStep(), false);
    }
    return std::make_pair(LineStep(loc.line, loc.step + 1), true);
  }();

  return SRepInterpolateHelper::ComputeDerivative(grid, loc, prevOptLineStep, nextOptLineStep, spokeType);
}

//----------------------------------------------------------------------------
srep::Vector3d SRepInterpolateHelper::ComputeDerivative(
  const Grid& grid,
  const LineStep& loc,
  const OptionalLineStep& lesserNeighbor,
  const OptionalLineStep& greaterNeighbor,
  SpokeType spokeType)
{
  if (!greaterNeighbor.second && !lesserNeighbor.second) {
    throw std::runtime_error("SRep location has neither lesser nor greater neighbors: ("
      + std::to_string(loc.line) + ", " + std::to_string(loc.step) + ")" );
  }

  if (!lesserNeighbor.second) {
    // no greaterNeighbor neighbor
    const auto head = SRepInterpolateHelper::GetSpoke(grid, greaterNeighbor.first, spokeType).GetSkeletalPoint();
    const auto tail = SRepInterpolateHelper::GetSpoke(grid, loc, spokeType).GetSkeletalPoint();
    // pout << "    Derivative from " << lesserNeighbor.first << " to " << loc << " is " << srep::Vector3d(tail, head) << "\n";
    return srep::Vector3d(tail, head);

  } else if (!greaterNeighbor.second) {
    // no lesserNeighbor neighbor
    const auto head = SRepInterpolateHelper::GetSpoke(grid, loc, spokeType).GetSkeletalPoint();
    const auto tail = SRepInterpolateHelper::GetSpoke(grid, lesserNeighbor.first, spokeType).GetSkeletalPoint();
    // pout << "    Derivative from " << loc << " to " << greaterNeighbor.first << " is " << srep::Vector3d(tail, head) << "\n";
    return srep::Vector3d(tail, head);
  } else {
    //neighbors on both sides
    const auto head = SRepInterpolateHelper::GetSpoke(grid, greaterNeighbor.first, spokeType).GetSkeletalPoint();
    const auto tail = SRepInterpolateHelper::GetSpoke(grid, lesserNeighbor.first, spokeType).GetSkeletalPoint();
    // pout << "    Derivative from " << lesserNeighbor.first << " to " << greaterNeighbor.first << " is " << (srep::Vector3d(tail, head)/2) << "\n";
    return srep::Vector3d(tail, head) / 2; // divide by 2 because we are going twice the distance
  }
}

//----------------------------------------------------------------------------
SRepInterpolateHelper::SkeletalPointDerivative
SRepInterpolateHelper::ComputeDerivative(const Grid& grid, const LineStep& loc) {
  SkeletalPointDerivative d;
  d.up.u = SRepInterpolateHelper::ComputeLinewiseDerivative(grid, loc, SpokeType::Up);
  d.down.u = SRepInterpolateHelper::ComputeLinewiseDerivative(grid, loc, SpokeType::Down);
  d.up.v = SRepInterpolateHelper::ComputeStepwiseDerivative(grid, loc, SpokeType::Up);
  d.down.v = SRepInterpolateHelper::ComputeStepwiseDerivative(grid, loc, SpokeType::Down);

  return d;
}

//----------------------------------------------------------------------------
SRepInterpolateHelper::DerivativeGridType
SRepInterpolateHelper::ComputeDerivatives(const Grid& grid) {
  DerivativeGridType d;
  d.reserve(grid.size());
  for (size_t line = 0; line < grid.size(); ++line) {
    d.resize(d.size()+1);
    d.back().reserve(grid[line].size());
    for (size_t step = 0; step< grid[line].size(); ++step) {
      d.back().push_back(ComputeDerivative(grid, LineStep(line, step)));
    }
  }
  return d;
}

//----------------------------------------------------------------------------
srep::Vector3d SRepInterpolateHelper::Slerp(const srep::Vector3d& v1, const srep::Vector3d& v2, const double u) {
  const double v1Tv2 = Clamp(v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2], -1, 1);
  const double phi = acos(v1Tv2);
  const auto theComputation = [&](double val1, double val2) {
    return (sin((1-u)*phi) / sin(phi)) * val1 + (sin(u*phi) / sin(phi)) * val2;
  };
  return srep::Vector3d (
    theComputation(v1[0], v2[0]),
    theComputation(v1[1], v2[1]),
    theComputation(v1[2], v2[2])
  );
}

//----------------------------------------------------------------------------
srep::Vector3d SRepInterpolateHelper::Compute2ndDerivative(
  const srep::Vector3d& startVector,
  const srep::Vector3d& endVector,
  const srep::Vector3d& targetVector,
  const double d)
{
  constexpr double del = 1e-5;
  const auto Upv1 = Slerp(startVector.Unit(), endVector.Unit(), d + 2*del);
  const auto Upv5 = Slerp(startVector.Unit(), endVector.Unit(), d - 2*del);
  const auto unitTargetVector = targetVector.Unit();
  return srep::Vector3d(
    0.25 * (Upv5[0] + Upv1[0] - 2.0 * unitTargetVector[0]),
    0.25 * (Upv5[1] + Upv1[1] - 2.0 * unitTargetVector[1]),
    0.25 * (Upv5[2] + Upv1[2] - 2.0 * unitTargetVector[2])
  );
}

//----------------------------------------------------------------------------
SRepInterpolateHelper::LineStep
SRepInterpolateHelper::OriginalLineStepToInterpolatedLineStep(const LineStep& ols) {
  return LineStep(ols.line * this->Density, ols.step * this->Density);
}

//----------------------------------------------------------------------------
SRepInterpolateHelper::Quad SRepInterpolateHelper::OriginalQuadToInterpolatedQuad(const Quad& oQuad) {
  Quad iQuad;
  std::transform(oQuad.begin(), oQuad.end(), iQuad.begin(),
    std::bind(&SRepInterpolateHelper::OriginalLineStepToInterpolatedLineStep, this, std::placeholders::_1));
  return iQuad;
}

//----------------------------------------------------------------------------
size_t SRepInterpolateHelper::LinewiseDistance(const LineStep& from, const LineStep& to, const Grid& grid) {
  const auto forwardDistance = (to.line + grid.size() - from.line) % grid.size();
  return std::min(forwardDistance, grid.size() - forwardDistance);
}

//----------------------------------------------------------------------------
size_t SRepInterpolateHelper::LinewiseDistance(const Quad& iQuad, const Grid& grid) {
  const auto length = LinewiseDistance(iQuad[0], iQuad[1], grid);
  const auto length2 = LinewiseDistance(iQuad[2], iQuad[3], grid);
  if (length != length2) {
    throw std::invalid_argument("Breaking assumptions in interpolation lines: " + std::to_string(length) + " != " + std::to_string(length2));
  }
  if (length == 0) {
    const auto length3 = LinewiseDistance(iQuad[0], iQuad[2], grid);
    const auto length4 = LinewiseDistance(iQuad[1], iQuad[3], grid);
    if (length3 != length4) {
      throw std::invalid_argument("Breaking assumptions in interpolation lines: " + std::to_string(length3) + " != " + std::to_string(length4));
    }
    return length3;
  }
  return length;
}

//----------------------------------------------------------------------------
size_t SRepInterpolateHelper::StepwiseDistance(const LineStep& from, const LineStep& to, const Grid& grid) {
  const auto forwardDistance = (to.step + grid[0].size() - from.step) % grid[0].size();
  return std::min(forwardDistance, grid[0].size() - forwardDistance);
}

//----------------------------------------------------------------------------
size_t SRepInterpolateHelper::StepwiseDistance(const Quad& iQuad, const Grid& grid) {
  const auto length = StepwiseDistance(iQuad[0], iQuad[2], grid);
  const auto length2 = StepwiseDistance(iQuad[1], iQuad[3], grid);
  if (length != length2) {
    throw std::invalid_argument("Breaking assumptions in interpolation steps: " + std::to_string(length) + " != " + std::to_string(length2));
  }
  if (length == 0) {
    const auto length3 = StepwiseDistance(iQuad[0], iQuad[1], grid);
    const auto length4 = StepwiseDistance(iQuad[2], iQuad[3], grid);
    if (length3 != length4) {
      throw std::invalid_argument("Breaking assumptions in interpolation lines: " + std::to_string(length3) + " != " + std::to_string(length4));
    }
    return length3;
  }
  return length;
}

//----------------------------------------------------------------------------
void SRepInterpolateHelper::InterpolateQuad(const Quad& iQuad, const Quad& originalEnclosingQuad, double lambda) {
  // Connectivity as from GetOrientedQuads
  // 0 - 1
  // |   |
  // 2 - 3

  const auto lineLength = this->LinewiseDistance(iQuad, this->InterpolatedGrid);
  const auto stepLength = this->StepwiseDistance(iQuad, this->InterpolatedGrid);
  if (lineLength != stepLength) {
    std::stringstream msg;
    msg << "Breaking assumptions in interpolation. Should be square. Found "
      << lineLength << "x" << stepLength
      << " for quad " << iQuad[0] << iQuad[1] << iQuad[2] << iQuad[3];
    throw std::invalid_argument(msg.str());
  }

  const auto length = lineLength;
  if (length <= 1) {
    // base case, there are no interior points to interpolate
    return;
  }

  if (!IsPowerOfTwo(length)) {
    throw std::invalid_argument("Breaking assumptions that interpolated density is power of two. Found "
      + std::to_string(length));
  }

  // interpolate centers of connection lines and center of whole quad
  // top,left,right,bottom here are used for convenience of understanding and are not necessarily
  // the actual top,left,right,bottom directions on the SRep

  // t == top
  // l == left
  // r == right
  // b == bottom
  // m == middle
  // therefore tm == top-middle, br = bottom-right, mm = middle-middle aka center of the quad, etc

  const auto tl = iQuad[0];
  const auto tr = iQuad[1];
  const auto bl = iQuad[2];
  const auto br = iQuad[3];

  const auto tm = InterpolateMiddleSkeletalPoint(tl, tr, originalEnclosingQuad, lambda);
  const auto lm = InterpolateMiddleSkeletalPoint(tl, bl, originalEnclosingQuad, lambda);
  const auto rm = InterpolateMiddleSkeletalPoint(tr, br, originalEnclosingQuad, lambda);
  const auto bm = InterpolateMiddleSkeletalPoint(bl, br, originalEnclosingQuad, lambda);

  // for the very center interpolate off of two directions and average
  const auto mmLeftRight = InterpolateMiddleSkeletalPoint(lm, rm, originalEnclosingQuad, lambda);
  const auto ptAtMMLeftRight = this->GetInterpolatedSkeletalPoint(mmLeftRight);
  const auto mmUpDown = InterpolateMiddleSkeletalPoint(tm, bm, originalEnclosingQuad, lambda);
  const auto ptAtMMTopBottom = this->GetInterpolatedSkeletalPoint(mmUpDown);

  if (mmUpDown.line != mmLeftRight.line) {
    throw std::runtime_error("bug in getting middle spot linewise");
  }
  if (mmUpDown.step != mmLeftRight.step) {
    throw std::runtime_error("bug in getting middle spot stepwise");
  }

  const auto& mm = mmLeftRight;
  this->InterpolatedGrid[mm.line][mm.step] = Average(ptAtMMLeftRight, ptAtMMTopBottom);

  if (length == 2) {
    // all the recursors will have length == 1 which will quit. Early return for less verbose
    // debugging statements
    return;
  }

  // recurse and interpolate the sub-quads
  InterpolateQuad({tl, tm, lm, mm}, originalEnclosingQuad, lambda / 2); // top left quad
  InterpolateQuad({tm, tr, mm, rm}, originalEnclosingQuad, lambda / 2); // top right quad
  InterpolateQuad({lm, mm, bl, bm}, originalEnclosingQuad, lambda / 2); // bottom left quad
  InterpolateQuad({mm, rm, bm, br}, originalEnclosingQuad, lambda / 2); // bottom right quad
}

//----------------------------------------------------------------------------
SRepInterpolateHelper::LineStep SRepInterpolateHelper::MiddleLineStep(const LineStep& a, const LineStep& b, const Grid& grid) {
  auto lineDist = LinewiseDistance(a, b, grid);
  auto stepDist = StepwiseDistance(a, b, grid);

  size_t line = 0;
  if ((a.line + lineDist) % grid.size() == b.line) {
    line = (a.line + (lineDist / 2)) % grid.size();
  } else {
    line = (b.line + (lineDist / 2)) % grid.size();
  }

  size_t step = 0;
  //there is no step wrap around
  if (a.step + stepDist == b.step) {
    step = a.step + (stepDist / 2);
  } else {
    step = b.step + (stepDist / 2);
  }
  return LineStep(line, step);
}

//----------------------------------------------------------------------------
SRepInterpolateHelper::UVDerivative SRepInterpolateHelper::Average(const UVDerivative& uv1, const UVDerivative& uv2) {
  UVDerivative avg;
  avg.u = Average(uv1.u, uv2.u);
  avg.v = Average(uv1.v, uv2.v);
  return avg;
}

//----------------------------------------------------------------------------
SRepInterpolateHelper::SkeletalPoint
SRepInterpolateHelper::Average(const SkeletalPoint& pt1, const SkeletalPoint& pt2) {
  if (pt1.IsCrest() != pt2.IsCrest()) {
    throw std::invalid_argument("How does one average two skeletal points when only one is a crest point?");
  }

  const auto avgUp = Average(pt1.GetUpSpoke(), pt2.GetUpSpoke());
  const auto avgDown = Average(pt1.GetDownSpoke(), pt2.GetDownSpoke());
  if (pt1.IsCrest()) {
    return SkeletalPoint(avgUp, avgDown, Average(pt1.GetCrestSpoke(), pt2.GetCrestSpoke()));
  } else {
    return SkeletalPoint(avgUp, avgDown);
  }
}

//----------------------------------------------------------------------------
srep::Spoke SRepInterpolateHelper::Average(const srep::Spoke& s1, const srep::Spoke& s2) {
  return srep::Spoke(
    Average(s1.GetSkeletalPoint(), s2.GetSkeletalPoint()),
    Average(s1.GetDirection(), s2.GetDirection()));
}

//----------------------------------------------------------------------------
srep::Vector3d SRepInterpolateHelper::Average(const srep::Vector3d& v1, const srep::Vector3d& v2) {
  return (v1 + v2) / 2;
}

//----------------------------------------------------------------------------
srep::Point3d SRepInterpolateHelper::Average(const srep::Point3d& pt1, const srep::Point3d& pt2) {
  // just pretend the points are vectors for a second
  return srep::Point3d(Average(srep::Vector3d(pt1.AsArray()), srep::Vector3d(pt2.AsArray())).AsArray());
}

//----------------------------------------------------------------------------
std::vector<SRepInterpolateHelper::Quad> SRepInterpolateHelper::GetOrientedQuads(const Grid& grid) {
  //Quad orientation: (smaller-line, smaller-step)(larger-line, smaller-step)(smaller-line, larger-step)(larger-line, larger-step)
  // i.e. clockwise from the inside out

  using LS = LineStep;
  std::vector<Quad> quads;

  // note, line wraps around, step does not
  for (size_t line = 0; line < grid.size(); ++line) {
    const auto nextLine = (line + 1) % grid.size();
    for (size_t step = 0; step < grid[line].size() - 1; ++step) {
      quads.push_back(Quad{LS(line, step), LS(nextLine, step), LS(line, step+1), LS(nextLine, step+1)});
    }
  }

  return quads;
}

//----------------------------------------------------------------------------
SRepInterpolateHelper::LineStep
SRepInterpolateHelper::InterpolateMiddleSkeletalPoint(
  const LineStep& start,
  const LineStep& end,
  const Quad& originalEnclosingQuad,
  double lambda)
{
  const auto midUpDirection = InterpolateMiddleSpokeDirection(start, end, lambda, SpokeType::Up);
  const auto midUpSkeletonPoint = InterpolateMiddleSkeletalPointSkeletonPoint(start, end, originalEnclosingQuad, SpokeType::Up);
  const srep::Spoke midUpSpoke(midUpSkeletonPoint, midUpDirection);

  const auto midDownDirection = InterpolateMiddleSpokeDirection(start, end, lambda, SpokeType::Down);
  const auto midDownSkeletonPoint = InterpolateMiddleSkeletalPointSkeletonPoint(start, end, originalEnclosingQuad, SpokeType::Down);
  const srep::Spoke midDownSpoke(midDownSkeletonPoint, midDownDirection);

  const auto midLS = MiddleLineStep(start, end, this->InterpolatedGrid);
  if (this->GetInterpolatedSkeletalPoint(start).IsCrest()
    && this->GetInterpolatedSkeletalPoint(end).IsCrest())
  {
    const auto midCrestDirection = InterpolateMiddleSpokeDirection(start, end, lambda, SpokeType::Crest);
    const auto midCrestSkeletonPoint = InterpolateMiddleSkeletalPointSkeletonPoint(start, end, originalEnclosingQuad, SpokeType::Crest);
    const srep::Spoke midCrestSpoke(midCrestSkeletonPoint, midCrestDirection);
    this->InterpolatedGrid[midLS.line][midLS.step] = srep::SkeletalPoint(midUpSpoke, midDownSpoke, midCrestSpoke);
  } else {
    this->InterpolatedGrid[midLS.line][midLS.step] = srep::SkeletalPoint(midUpSpoke, midDownSpoke);
  }
  return midLS;
}

//----------------------------------------------------------------------------
srep::Vector3d SRepInterpolateHelper::InterpolateMiddleSpokeDirection(
  const LineStep& start,
  const LineStep& end,
  double lambda,
  SpokeType spokeType)
{
  return this->InterpolateMiddleSpokeDirection(
    this->GetInterpolatedSpoke(start, spokeType),
    this->GetInterpolatedSpoke(end, spokeType),
    lambda);
}

//----------------------------------------------------------------------------
srep::Vector3d SRepInterpolateHelper::InterpolateMiddleSpokeDirection(
  const srep::Spoke& startSpoke,
  const srep::Spoke& endSpoke,
  const double lambda)
{
  const auto startUnitDirection = startSpoke.GetDirection().Unit();
  const auto endUnitDirection = endSpoke.GetDirection().Unit();
  const auto start2ndDerivative = Compute2ndDerivative(startUnitDirection, endUnitDirection, startUnitDirection, 0);
  const auto end2ndDerivative = Compute2ndDerivative(startUnitDirection, endUnitDirection, endUnitDirection, lambda);
  const auto avgSpokeDirection = Average(startSpoke.GetDirection(), endSpoke.GetDirection());
  const double halfDist = lambda / 2;
  //if startSpoke == endSpoke then interpolated spoke == both
  //I don't think this should ever really happen
  if (startSpoke == endSpoke) {
    return startSpoke.GetDirection();
  }
  const auto middleUnitDirection = Slerp(startUnitDirection, endUnitDirection, halfDist);
  const auto innerProd = [](const srep::Vector3d& a, const srep::Vector3d b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
  };
  const double innerProd1 = innerProd(middleUnitDirection, avgSpokeDirection);
  const double innerProd2 = innerProd(startUnitDirection, start2ndDerivative);
  const double innerProd3 = innerProd(endUnitDirection, end2ndDerivative);
  const double interpolatedRadius = innerProd1 - (halfDist * halfDist * 0.25 * (innerProd2 + innerProd3));
  return middleUnitDirection * interpolatedRadius;
}

//----------------------------------------------------------------------------
SRepInterpolateHelper::UVDerivative
SRepInterpolateHelper::GetUVDerivativeFromOriginalLineStep(const LineStep& ols, SpokeType spokeType) {
  const auto& d = this->DerivativeOriginalGrid[ols.line][ols.step];
  if (spokeType == SpokeType::Up) {
    return d.up;
  } else if (spokeType == SpokeType::Down) {
    return d.down;
  } else if (spokeType == SpokeType::Crest) {
    //average of up and down
    return Average(d.up, d.down);
  }
  throw std::invalid_argument("Unknown spoke type");
}

//----------------------------------------------------------------------------
double SRepInterpolateHelper::h1(double s) { return 2*(s * s * s) - 3*(s * s) + 1; }
double SRepInterpolateHelper::h2(double s) { return -2*(s * s * s) + 3*(s * s); }
double SRepInterpolateHelper::h3(double s) { return (s * s * s) - 2*(s * s) + s; }
double SRepInterpolateHelper::h4(double s) { return (s * s * s) - (s * s); }

//----------------------------------------------------------------------------
srep::Point3d SRepInterpolateHelper::InterpolateMiddleSkeletalPointSkeletonPoint(
  const LineStep& start,
  const LineStep& end,
  const Quad& originalEnclosingQuad,
  SpokeType spokeType)
{
  const auto mid = MiddleLineStep(start, end, this->InterpolatedGrid);
  if (spokeType == SpokeType::Up || spokeType == SpokeType::Down) {
    return this->InterpolateSkeletalPointSkeletonPoint(mid, originalEnclosingQuad, spokeType);
  } else if (spokeType == SpokeType::Crest) {
    // interpolation is a little weird for the crest. Create a fake quad that is really just a line and interpolate
    // with that
    std::vector<LineStep> crestLocs;
    for (const auto& loc : originalEnclosingQuad) {
      if (this->OriginalGrid[loc.line][loc.step].IsCrest()) {
        crestLocs.push_back(loc);
      }
    }
    if (crestLocs.size() != 2) {
      throw std::runtime_error("Cannot interpolate crest when the original enclosing quad does not contain two crest points");
    }
    return this->InterpolateSkeletalPointSkeletonPoint(mid, Quad{crestLocs[0], crestLocs[1], crestLocs[0], crestLocs[1]}, spokeType);
  } else {
    throw std::invalid_argument("Unknown spoke type");
  }
}

//----------------------------------------------------------------------------
srep::Point3d SRepInterpolateHelper::InterpolateSkeletalPointSkeletonPoint(
  const LineStep& loc,
  Quad originalEnclosingQuad,
  SpokeType spokeType)
{
  const auto interpolatedEnclosingQuad = OriginalQuadToInterpolatedQuad(originalEnclosingQuad);

  const auto x11 = this->GetInterpolatedSpoke(interpolatedEnclosingQuad[0], spokeType).GetSkeletalPoint();
  const auto x21 = this->GetInterpolatedSpoke(interpolatedEnclosingQuad[1], spokeType).GetSkeletalPoint();
  const auto x12 = this->GetInterpolatedSpoke(interpolatedEnclosingQuad[2], spokeType).GetSkeletalPoint();
  const auto x22 = this->GetInterpolatedSpoke(interpolatedEnclosingQuad[3], spokeType).GetSkeletalPoint();
  const auto dx11 = GetUVDerivativeFromOriginalLineStep(originalEnclosingQuad[0], spokeType);
  const auto dx21 = GetUVDerivativeFromOriginalLineStep(originalEnclosingQuad[1], spokeType);
  const auto dx12 = GetUVDerivativeFromOriginalLineStep(originalEnclosingQuad[2], spokeType);
  const auto dx22 = GetUVDerivativeFromOriginalLineStep(originalEnclosingQuad[3], spokeType);

  const auto dxdu11 = dx11.u;
  const auto dxdv11 = dx11.v;
  const auto dxdu12 = dx12.u;
  const auto dxdv12 = dx12.v;
  const auto dxdu21 = dx21.u;
  const auto dxdv21 = dx21.v;
  const auto dxdu22 = dx22.u;
  const auto dxdv22 = dx22.v;

  const auto oQuadLineLength = this->LinewiseDistance(interpolatedEnclosingQuad, this->InterpolatedGrid);
  const auto oQuadStepLength = this->StepwiseDistance(interpolatedEnclosingQuad, this->InterpolatedGrid);

  const auto startQuadToLocLineLength = this->LinewiseDistance(interpolatedEnclosingQuad[0], loc, this->InterpolatedGrid);
  const auto startQuadToLocStepLength = this->StepwiseDistance(interpolatedEnclosingQuad[0], loc, this->InterpolatedGrid);

  const double u = oQuadLineLength > 0 ? static_cast<double>(startQuadToLocLineLength) / oQuadLineLength : 0.0;
  const double v = oQuadStepLength > 0 ? static_cast<double>(startQuadToLocStepLength) / oQuadStepLength : 0.0;

  // this was pulled as is from original implementation
  double hx[4][4];
  double hy[4][4];
  double hz[4][4];

  hx[0][0] = x11[0];          hx[0][1] = x12[0];
  hx[1][0] = x21[0];          hx[1][1] = x22[0];
  hx[2][0] = dxdu11[0];       hx[2][1] = dxdu12[0];
  hx[3][0] = dxdu21[0];       hx[3][1] = dxdu22[0];
  hx[0][2] = dxdv11[0];       hx[0][3] = dxdv12[0];
  hx[1][2] = dxdv21[0];       hx[1][3] = dxdv22[0];
  hx[2][2] = 0;               hx[2][3] = 0;
  hx[3][2] = 0;               hx[3][3] = 0;

  hy[0][0] = x11[1];          hy[0][1] = x12[1];
  hy[1][0] = x21[1];          hy[1][1] = x22[1];
  hy[2][0] = dxdu11[1];       hy[2][1] = dxdu12[1];
  hy[3][0] = dxdu21[1];       hy[3][1] = dxdu22[1];
  hy[0][2] = dxdv11[1];       hy[0][3] = dxdv12[1];
  hy[1][2] = dxdv21[1];       hy[1][3] = dxdv22[1];
  hy[2][2] = 0;               hy[2][3] = 0;
  hy[3][2] = 0;               hy[3][3] = 0;

  hz[0][0] = x11[2];       hz[0][1] = x12[2];
  hz[1][0] = x21[2];       hz[1][1] = x22[2];
  hz[2][0] = dxdu11[2];    hz[2][1] = dxdu12[2];
  hz[3][0] = dxdu21[2];    hz[3][1] = dxdu22[2];
  hz[0][2] = dxdv11[2];    hz[0][3] = dxdv12[2];
  hz[1][2] = dxdv21[2];    hz[1][3] = dxdv22[2];
  hz[2][2] = 0;            hz[2][3] = 0;
  hz[3][2] = 0;            hz[3][3] = 0;

  double hu[4], hv[4];
  hu[0] = h1(u);
  hu[1] = h2(u);
  hu[2] = h3(u);
  hu[3] = h4(u);
  hv[0] = h1(v);
  hv[1] = h2(v);
  hv[2] = h3(v);
  hv[3] = h4(v);

  // supposed computation is these
  //    vnl_double_1x1 xn = hu.transpose() * hx * hv;
  //    vnl_double_1x1 yn = hu.transpose() * hy * hv;
  //    vnl_double_1x1 zn = hu.transpose() * hz * hv;
  double huThx[4], huThy[4], huThz[4];
  huThx[0] = hu[0] * hx[0][0] + hu[1] * hx[1][0] + hu[2] * hx[2][0] + hu[3] * hx[3][0];
  huThx[1] = hu[0] * hx[0][1] + hu[1] * hx[1][1] + hu[2] * hx[2][1] + hu[3] * hx[3][1];
  huThx[2] = hu[0] * hx[0][2] + hu[1] * hx[1][2] + hu[2] * hx[2][2] + hu[3] * hx[3][2];
  huThx[3] = hu[0] * hx[0][3] + hu[1] * hx[1][3] + hu[2] * hx[2][3] + hu[3] * hx[3][3];

  huThy[0] = hu[0] * hy[0][0] + hu[1] * hy[1][0] + hu[2] * hy[2][0] + hu[3] * hy[3][0];
  huThy[1] = hu[0] * hy[0][1] + hu[1] * hy[1][1] + hu[2] * hy[2][1] + hu[3] * hy[3][1];
  huThy[2] = hu[0] * hy[0][2] + hu[1] * hy[1][2] + hu[2] * hy[2][2] + hu[3] * hy[3][2];
  huThy[3] = hu[0] * hy[0][3] + hu[1] * hy[1][3] + hu[2] * hy[2][3] + hu[3] * hy[3][3];

  huThz[0] = hu[0] * hz[0][0] + hu[1] * hz[1][0] + hu[2] * hz[2][0] + hu[3] * hz[3][0];
  huThz[1] = hu[0] * hz[0][1] + hu[1] * hz[1][1] + hu[2] * hz[2][1] + hu[3] * hz[3][1];
  huThz[2] = hu[0] * hz[0][2] + hu[1] * hz[1][2] + hu[2] * hz[2][2] + hu[3] * hz[3][2];
  huThz[3] = hu[0] * hz[0][3] + hu[1] * hz[1][3] + hu[2] * hz[2][3] + hu[3] * hz[3][3];

  double output[3];
  output[0] = huThx[0] * hv[0] + huThx[1] * hv[1] + huThx[2] * hv[2] + huThx[3] * hv[3];
  output[1] = huThy[0] * hv[0] + huThy[1] * hv[1] + huThy[2] * hv[2] + huThy[3] * hv[3];
  output[2] = huThz[0] * hv[0] + huThz[1] * hv[1] + huThz[2] * hv[2] + huThz[3] * hv[3];

  return srep::Point3d(output);
}

//----------------------------------------------------------------------------
// Public functions
//----------------------------------------------------------------------------

//----------------------------------------------------------------------------
SRepInterpolateHelper::SRepInterpolateHelper(size_t interpolationLevel, const srep::EllipticalSRep& srep)
  : InterpolationLevel(interpolationLevel)
  , Density(IntegerPower(2, InterpolationLevel))
  , OriginalGrid(srep.GetSkeleton())
  , InterpolatedGrid()
  , DerivativeOriginalGrid(this->ComputeDerivatives(this->OriginalGrid))
{
  if (this->InterpolationLevel < 1) {
    throw std::invalid_argument("Invalid interpolation level");
  }

  if (srep.IsEmpty()) {
    throw std::invalid_argument("Can't interpolate empty srep");
  }
}

//----------------------------------------------------------------------------
std::unique_ptr<srep::EllipticalSRep> SRepInterpolateHelper::interpolate() {
  this->InterpolatedGrid = Grid(); //reset the grid back to a known good state
  this->InterpolatedGrid.resize(OriginalGrid.size() * this->Density);
  for (size_t i = 0; i < InterpolatedGrid.size(); ++i) {
    // Rationale for the -1 and +1: consider this line ("O" == skeletal point, "-" == connection)
    //    O---O---O
    // The line has 3 points on it and 2 connections. If we interpolate this to be twice as dense,
    // we will get the following line
    //    O-O-O-O-O
    // because we interpolate, not extrapolate (the start and end don't change)
    // So this is the math:
    //         numStepsFromSpine          *    Density     + SpinePoint
    // (this->OriginalGrid[0].size() - 1) * this->Density) + 1
    this->InterpolatedGrid[i].resize(((this->OriginalGrid[0].size() - 1) * this->Density) + 1);
  }

  const auto originalQuads = GetOrientedQuads(this->OriginalGrid);

  ///////////////////////////////////////////////
  // Copy original quads
  ///////////////////////////////////////////////
  // copy over the actual skeletal points that don't need to be
  // interpolated
  for (const auto& quad : originalQuads) {
    for (const auto& ols : quad) {
      const auto ils = this->OriginalLineStepToInterpolatedLineStep(ols);
      this->InterpolatedGrid[ils.line][ils.step] = this->OriginalGrid[ols.line][ols.step];
    }
  }

  ///////////////////////////////////////////////
  // Interpolate quads
  ///////////////////////////////////////////////
  for (const auto& oQuad : originalQuads) {
    const auto iQuad = this->OriginalQuadToInterpolatedQuad(oQuad);
    this->InterpolateQuad(iQuad, oQuad);
  }

  return std::unique_ptr<srep::EllipticalSRep>(new srep::EllipticalSRep(std::move(this->InterpolatedGrid)));
}


} //namespace detail
} //namespace sreplogic
