#ifndef __srep_EllipticalSRep_h
#define __srep_EllipticalSRep_h

#include <functional>

#include <srep/MeshSRepInterface.h>
#include <srep/SkeletalPoint.h>
#include <srep/Exceptions.h>

namespace srep {

/// An SRep with an elliptical grid
///
/// Grid is somewhat elliptical, but not in polar coordinates.
/// Instead of all lines coming from outer boundary of the ellipse
/// meeting at the center, there is a "spine" in the ellipse that
/// the lines go to.
///
/// Instead of a polar system like this (where all lines meet at the center .)
///   \ | /
/// --- . ---
///   / | \                (characters to prevent multi-line comment)
///
/// It looks more like this (. == spine)
///       \  |  /
///  --- . . . . . ---
///       /  |  \          (characters to prevent multi-line comment)
///
/// where each point on the spine has 2 lines going out of it,
/// with the exception of the poles on the longer axis, which have 1
class EllipticalSRep : public MeshSRepInterface {
public:
  /// Unrolled elliptical grid. first "row" is the spine.
  /// There will be duplicate spine points as there is a "column" for
  /// both sides of the spine. These duplicate points will not show
  /// up as two separate points in the mesh representation.
  /// The last "row" in the grid will be the crest (aka fold).
  /// First line is expected to be one of the poles
  using LineOutFromSpine = std::vector<SkeletalPoint>;
  using UnrolledEllipticalGrid = std::vector<LineOutFromSpine>;
  EllipticalSRep();
  EllipticalSRep(UnrolledEllipticalGrid skeleton);
  ~EllipticalSRep() = default;

  //copy and move deleted
  EllipticalSRep(const EllipticalSRep&) = delete;
  EllipticalSRep& operator=(const EllipticalSRep&) = delete;
  EllipticalSRep(EllipticalSRep&&) = delete;
  EllipticalSRep& operator=(EllipticalSRep&&) = delete;

  const UnrolledEllipticalGrid& GetSkeleton() const;

  ///////////////////////////////////////////////////////////
  //
  // MeshSRepInterface overrides
  //
  ///////////////////////////////////////////////////////////
  util::owner<EllipticalSRep*> Clone() const override;
  bool IsEmpty() const override;
  const SpokeMesh& GetUpSpokes() const override;
  const SpokeMesh& GetDownSpokes() const override;
  const SpokeMesh& GetCrestSpokes() const override;
  const std::vector<UpDownIndices>& GetCrestSkeletalConnections() const override;
  const std::vector<UpDownIndices>& GetSpine() const override;
private:
  static void Validate(const UnrolledEllipticalGrid& skeleton);
  void CreateMeshRepresentation();

  UnrolledEllipticalGrid Skeleton;

  struct MeshRepresentation {
      SpokeMesh UpSpokes;
      SpokeMesh DownSpokes;
      SpokeMesh CrestSpokes;
      std::vector<UpDownIndices> CrestSkeletalConnections;
      std::vector<UpDownIndices> Spine;
  };

  MeshRepresentation SkeletonAsMesh;
};

// The layout of the elliptical SRep skeleton visually to help understand up/down/left/right
// especially when crossing the spine. Note there are duplicate points on the spine
//
// 8 fold points, 2 steps
//
// (line, step) == UnrolledEllipticalGrid indicies
//
// Up              (1,2)           (2,2)         (3,2)
//  ^                  (1,1)       (2,1)      (3,1)
//  |   (0,2) (0,1) (0,0) (1&7,0) (2&6,0) (3&5,0) (4,0) (4,1) (4,2)
//  v                  (7,1)       (6,1)      (5,1)
// Down            (7,2)           (6,2)         (5,2)
//
//      Left                       <--->                      Right
//
// Note: on the top of the ellipse, up is a positive step movement, while on the
//       bottom, up is a negative step movement
// Note: on the top of the ellipse, left is a negative line movement, while on the
//       bottom, left is a positive line movement
// Note: to the get the line opposite you on the ellipse, do (numFoldPoints - lineIndex)
//       if you are not at a pole (then there is not opposite line)

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

/// @{
/// Gets the neighbor in the specified direction.
///
/// @return Pair. First is LineStep of neighbor, second is if first is valid. If second is false, then no neighbor in that direction
std::pair<LineStep, bool> GetLeftNeighbor(const EllipticalSRep::UnrolledEllipticalGrid& grid, size_t line, size_t step);
std::pair<LineStep, bool> GetRightNeighbor(const EllipticalSRep::UnrolledEllipticalGrid& grid, size_t line, size_t step);
std::pair<LineStep, bool> GetTopNeighbor(const EllipticalSRep::UnrolledEllipticalGrid& grid, size_t line, size_t step);
std::pair<LineStep, bool> GetBottomNeighbor(const EllipticalSRep::UnrolledEllipticalGrid& grid, size_t line, size_t step);
/// @}

using GetNeighborFunc = std::function<std::pair<LineStep, bool>(const EllipticalSRep::UnrolledEllipticalGrid& grid, size_t line, size_t step)>;

/// Gets all the neighbors in no particular order
std::vector<LineStep> GetNeighbors(const EllipticalSRep::UnrolledEllipticalGrid& grid, size_t line, size_t step);

}

#endif
