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

#ifndef __vtkSlicerSRepLogic_SRepInterpolator_h
#define __vtkSlicerSRepLogic_SRepInterpolator_h

#include <cstdlib>
#include <memory>
#include <vtkEllipticalSRep.h>

namespace sreplogic {

namespace detail {

struct LineStep {
  size_t line;
  size_t step;

  LineStep();
  LineStep(size_t line_, size_t step_);
  LineStep(const LineStep&) = default;
  LineStep& operator=(const LineStep&) = default;
  LineStep(LineStep&&) = default;
  LineStep& operator=(LineStep&&) = default;
  ~LineStep() = default;
};

bool operator<(const LineStep& a, const LineStep& b);
bool operator>(const LineStep& a, const LineStep& b);
bool operator<=(const LineStep& a, const LineStep& b);
bool operator>=(const LineStep& a, const LineStep& b);
bool operator==(const LineStep& a, const LineStep& b);
bool operator!=(const LineStep& a, const LineStep& b);

std::ostream& operator<<(std::ostream& os, const LineStep& ls);

/// Warning: this must not outlive the srep passed into the constructor.
/// Recommend using InterpolateSRep function to ensure this
class SRepInterpolateHelper {
public:
  SRepInterpolateHelper(size_t interpolationLevel, const vtkEllipticalSRep& srep);

  vtkSmartPointer<vtkEllipticalSRep> interpolate();

private:
  using Grid = std::vector<std::vector<vtkSmartPointer<vtkSRepSkeletalPoint>>>;
  using Quad = std::array<LineStep, 4>;
  using OptionalLineStep = std::pair<LineStep, bool>;
  using SpokeType = vtkSRepSkeletalPoint::SpokeOrientation;

  struct UVDerivative {
    srep::Vector3d u;
    srep::Vector3d v;
  };

  struct SkeletalPointDerivative {
    UVDerivative up;
    UVDerivative down;
  };

  using DerivativeGridType = std::vector<std::vector<SkeletalPointDerivative>>;

  static Grid ToGrid(const vtkEllipticalSRep& srep);
  static vtkSmartPointer<vtkEllipticalSRep> FromGrid(Grid grid);

  static std::vector<Quad> GetOrientedQuads(const Grid& grid);

  static vtkSRepSpoke& GetSpoke(const Grid& grid, const LineStep& loc, SpokeType spokeType);
  vtkSRepSpoke& GetInterpolatedSpoke(const LineStep& loc, SpokeType spokeType);
  vtkSRepSkeletalPoint& GetInterpolatedSkeletalPoint(const LineStep& loc);

  //----------------------------------------------------------------------------
  // Computing grid derivatives
  static srep::Vector3d ComputeLinewiseDerivative(const Grid& grid, const LineStep& loc, SpokeType spokeType);
  static srep::Vector3d ComputeStepwiseDerivative(const Grid& grid, const LineStep& loc, SpokeType spokeType);
  static srep::Vector3d ComputeDerivative(
    const Grid& grid,
    const LineStep& loc,
    const OptionalLineStep& lesserNeighbor,
    const OptionalLineStep& greaterNeighbor,
    SpokeType spokeType);

  static SkeletalPointDerivative ComputeDerivative(const Grid& grid, const LineStep& loc);
  static DerivativeGridType ComputeDerivatives(const Grid& grid);

  //----------------------------------------------------------------------------
  // Interpolation
  /// @param u value between 0 and 1 that is kind of a weighting between the two vectors
  static srep::Vector3d Slerp(const srep::Vector3d& v1, const srep::Vector3d& v2, const double u);
  static srep::Vector3d Compute2ndDerivative(
    const srep::Vector3d& startVector,
    const srep::Vector3d& endVector,
    const srep::Vector3d& targetVector,
    const double d);

  LineStep OriginalLineStepToInterpolatedLineStep(const LineStep& ols);
  Quad OriginalQuadToInterpolatedQuad(const Quad& oQuad);
  static size_t LinewiseDistance(const Quad& quad, const Grid& grid);
  static size_t StepwiseDistance(const Quad& quad, const Grid& grid);
  static size_t LinewiseDistance(const LineStep& from, const LineStep& to, const Grid& grid);
  static size_t StepwiseDistance(const LineStep& from, const LineStep& to, const Grid& grid);

  /// Fill in all interpolated skeletal points within the grid
  /// @note The values in this->InterpolatedGrid at iQuad should be valid
  /// @param iQuad Quad for this->InterpolatedGrid. Needs to be oriented correctly. Only use quads
  ///        returned from GetOrientedQuads
  /// @param lambda The percent of the non-interpolated quad we are. So if the quad is
  ///               the spokes from the original srep, then lambda is 1.0, or 100% of the original
  ///               non-interpolated quad. If the quad is the top-left quadrant of the original
  ///               quad, then lambda is 0.5, or 50% of the original quad in each direction
  void InterpolateQuad(const Quad& iQuad, const Quad& originalEnclosingQuad, double lambda = 1.0);

  /// Interpolates the middle point of the given points. Fills it into
  /// this->InterpolatedGrid
  /// @param start The start location
  /// @param end The end location
  /// @param quad The quad that contains start and end
  /// @param lambda The lambda from InterpolateQuad
  /// @returns The location of the middle point
  LineStep InterpolateMiddleSkeletalPoint(
    const LineStep& start,
    const LineStep& end,
    const Quad& originalEnclosingQuad,
    double lambda);
  srep::Vector3d InterpolateMiddleSpokeDirection(
    const LineStep& start,
    const LineStep& end,
    double lambda,
    SpokeType spokeType);
  static srep::Vector3d InterpolateMiddleSpokeDirection(
    const vtkSRepSpoke& startSpoke,
    const vtkSRepSpoke& endSpoke,
    const double lambda);

  /// a and b cannot be across the spine
  static LineStep MiddleLineStep(const LineStep& a, const LineStep& b, const Grid& grid);
  static UVDerivative Average(const UVDerivative& uv1, const UVDerivative& uv2);
  static vtkSmartPointer<vtkSRepSkeletalPoint> Average(const vtkSRepSkeletalPoint& pt1, const vtkSRepSkeletalPoint& pt2);
  static vtkSmartPointer<vtkSRepSpoke> Average(const vtkSRepSpoke& s1, const vtkSRepSpoke& s2);
  static srep::Point3d Average(const srep::Point3d& pt1, const srep::Point3d& pt2);
  static srep::Vector3d Average(const srep::Vector3d& v1, const srep::Vector3d& v2);

  //----------------------------------------------------------------------------
  // Interpolation of Skeletal point Skeleton point

  // I don't know the purpose or reasoning behind these functions
  static double h1(double s);
  static double h2(double s);
  static double h3(double s);
  static double h4(double s);

  srep::Point3d InterpolateMiddleSkeletalPointSkeletonPoint(
    const LineStep& start,
    const LineStep& end,
    const Quad& originalEnclosingQuad,
    SpokeType spokeType);
  srep::Point3d InterpolateSkeletalPointSkeletonPoint(
    const LineStep& loc,
    Quad originalEnclosingQuad,
    SpokeType spokeType);

  UVDerivative GetUVDerivativeFromOriginalLineStep(const LineStep& ols, SpokeType spokeType);

  //----------------------------------------------------------------------------
  // Members
  const size_t InterpolationLevel;
  const size_t Density;
  const Grid OriginalGrid;
  Grid InterpolatedGrid;
  DerivativeGridType DerivativeOriginalGrid;
};
}

VTK_NEWINSTANCE vtkEllipticalSRep* InterpolateSRep(size_t interpolationLevel, const vtkEllipticalSRep& srep);
vtkSmartPointer<vtkEllipticalSRep> SmartInterpolateSRep(size_t interpolationLevel, const vtkEllipticalSRep& srep);

}

#endif
